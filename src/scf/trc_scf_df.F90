!
! Density-fitted J and K for a restricted (closed-shell) reference,
! resident on the device. Every array THIS object owns -- `bpk`, the
! packed density/J factors, the K half-transform -- is Cartesian-sized
! (`n = b%nao`), whether the orbital basis `b` is spherical or not: a
! spherical caller folds S/H/D/F through `trc_scf_fock.F90`'s `to_cart`/
! `to_sph` at the boundary, and this object never sees the difference.
! The AUXILIARY basis is a separate question (see `df_jk_setup`'s own
! comment): it follows its own `aux%spherical`, independent of `b`'s.
!
! With the auxiliary basis {P}, its metric M = (P|Q) = L L^T and the
! three-index tensor T(mu nu, P) = (mu nu|P),
!
!     B(mu nu, Q) = sum_P T(mu nu, P) (L^-T)_PQ      so that
!     (mu nu|lam sig) = sum_Q B(mu nu, Q) B(lam sig, Q)
!
! B is symmetric in (mu, nu) -- (mu nu|P) does not see which index came
! first -- so it is stored PACKED, one copy of mu >= nu, which halves the
! device memory a dense (nao, nao, naux) tensor would need (21 GB packed
! against 43 GB dense at w40 cc-pVDZ).
!
! THE BUILD, PER ITERATION
! ------------------------
! J, two GEMVs on the packed storage (bandwidth-bound: no transpose, no
! unpacking, 8 bytes moved per (pair, aux function) touched, twice):
!
!     gamma(Q) = sum_{mu>=nu} w(mu,nu) B(mu nu, Q) D(mu,nu),   w = 2 off the
!                                                               diagonal, 1 on it
!     J(mu,nu) = sum_Q B(mu nu, Q) gamma(Q)                    (packed; mu,nu and
!                                                               nu,mu read the same entry)
!
! K needs a density factor, D = Y Y^T. The SCF's density is not always
! idempotent (damping mixes densities; SAD/SAP are not from orbitals), so
! this always diagonalises it (`lp%syev`, one nao_c x nao_c symmetric
! eigenproblem per iteration -- see df_jk_build's own header for pyscf's
! rule and why the orbital-only shortcut is not taken yet) and keeps every
! eigenvalue with |n_k| > EVAL_REL_TOL max|n|: Y_pos = U sqrt(n_k) for the positive ones,
! Y_neg = U sqrt(|n_k|) for the negative ones (K is linear in D, so
! K[D] = K[D_pos] - K[D_neg] exactly). Then, one auxiliary block at a time:
!
!     unpack B_q to a dense (n, n, qb) scratch (both triangles filled)
!     X_q = B_q Y            one strided GEMM batch over the block, giving
!                             Z(mu, (i,q)), n x (o qb), per factor
!     K += Z Z^T             lower block panels only, one GEMM each (cuBLAS
!                             dsyrk ran at 38% of peak here), accumulated
!                             over blocks and mirrored once at the end
!
! and G = j_scale J[D] - (k_scale/2) K[D] for D the TOTAL density, the SAME convention `trc_eri_t%fock_resident`
! returns (see trc_scf_fock.F90's build_g, which reads this object's `build`
! in its own branch) -- `--df` and the four-centre path must agree up to the
! fitting error, and matching the sign/scale convention exactly is how.
!
! PACKED INDEXING. idx(mu, nu) = mu(mu-1)/2 + nu for mu >= nu (1-based,
! lower-triangular, mu the slower index): `pack_idx`, a module-level pure
! function (`!$acc routine seq`), not a type-bound one, for the same reason
! RI-MP2's `pair_energy_body` is module-level -- a `do concurrent` body may
! not call a procedure from another module, but a local pure function is
! exactly the pattern already in use (trc_rimp2_driver.F90's `ncart`).
!
! SETUP, ONCE PER GEOMETRY: metric, Cholesky, (mu nu|P) by auxiliary block
! (trc_df_3c, same blocking RI-MP2 uses for its own three-index build),
! packed into `bpk` as each block lands, then B = T L^-T, one triangular
! solve from the right over every packed row at once (`trc_linalg%trsm`,
! the same call RI-MP2's step 4 makes on its own packed rows).
!
! OWNERSHIP AND REUSE. `setup` is create-once (the `ready` flag, exactly
! trc_linalg_t%init's own pattern): a caller that owns this object across
! more than one `trc_scf_run` -- trc_capi.F90's context handle is the one
! that does -- gets the metric/Cholesky/(mu nu|P)/trsm done exactly once,
! and a correlated step that follows (RI-MP2) can read `bpk` back out
! (`packed_b`) for the same reason: it is already (mu nu|P) contracted
! against L^-T, over the same auxiliary basis, with nothing more to redo.
!
module trc_scf_df
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use trc_sph, only: sph_metric_cart_to_sph, sph_aux3c_to_sph
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use, intrinsic :: iso_fortran_env, only: output_unit
#ifdef TRC_DEVICE
   use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                      acc_property_free_memory
#endif
   implicit none
   private

   public :: df_jk_t

   integer, parameter :: NSTEP = 11
   integer, parameter :: ST_METRIC = 1, ST_CHOL = 2, ST_3C = 3, ST_PACK = 4, ST_TRSM = 5, &
                         ST_JGAM = 6, ST_JBLD = 7, ST_DFAC = 8, ST_UNPK = 9, ST_KHLF = 10, ST_KSYRK = 11
   character(len=18), parameter :: NAMES(NSTEP) = [character(len=18) :: &
                                                   "metric (P|Q)", "Cholesky", "(mu nu|P)", "pack B", &
                                                   "B = T L^-T", "J gamma", "J build", "density factor", &
                                                   "unpack B", "K half-transform", "K accum (gemm)"]
   real(dp), parameter :: EVAL_REL_TOL = 1.0e-10_dp
   !! RELATIVE to the density factor's own largest |eigenvalue|, not an
   !! absolute 1e-12 -- the absolute cut kept every eigenvalue down near
   !! machine noise regardless of scale, which on a damped (non-
   !! idempotent) density's genuinely continuous spectrum between 0 and
   !! the occupation number kept dozens of factor columns with no
   !! physical weight (measured: ~240 average at w40, nocc=200, before
   !! this change). A dropped eigenvalue n_k contributes a K element of
   !! at most |n_k| * max|(mu nu|lam sig)| -- bounded here by max|n|
   !! itself, since (mu nu|lam sig) <= 1 in the orthonormal-orbital basis
   !! Y is built in -- so every dropped eigenvalue biases K by at most
   !! EVAL_REL_TOL * max|n| ~ 1e-10 in the units D itself is in (electron
   !! count), two orders tighter than the 1e-8 |FDS-SDF| convergence gate
   !! this driver already stops at, and dropping O(10) such eigenvalues
   !! per iteration (the typical noise-rank size cut here) still sums to
   !! under the gate. Verified empirically below: w20/w40 cc-pVDZ and
   !! cc-pVTZ energies are unchanged to 1e-10 Eh with this cut against
   !! the old absolute one.

   type :: df_jk_t
      private
      logical :: ready = .false.
      integer :: n = 0, naux = 0, npack = 0
      integer :: qb = 1       !! auxiliary functions per K block, chosen at setup
      integer :: omax = 0     !! current capacity of the density-factor workspace
      real(dp), allocatable :: bpk(:, :)              !! (npack, naux): B, packed mu >= nu
      real(dp), allocatable :: chol(:, :)              !! (naux, naux): L, the metric's Cholesky
      !! factor (M = (P|Q) = L L^T), RETAINED past setup rather than freed --
      !! unlike the energy build, which never needs L again once B = T L^-T
      !! has absorbed it, the DF GRADIENT (trc_scf_grad_df.F90) needs one
      !! more triangular solve against the SAME L to turn the Cholesky-space
      !! coefficients `build` already computes (gamma(Q), and the occupied
      !! half-transform) into the actual fitting coefficients the analytic
      !! derivative's envelope-theorem formula is stationary in. Keeping L
      !! resident is what lets that solve skip rebuilding the metric or its
      !! factorisation -- the paper's point (Stocks/Palethorpe/Barca 2024,
      !! .claude/notes/papers/): never recompute a 2c/3c integral's VALUE.
      real(dp), allocatable :: dpk(:), jpk(:), gq(:)   !! J's packed density/result, and gamma(Q)
      real(dp), allocatable :: dsc(:, :), evals(:)     !! the density-factor eigenproblem
      real(dp), allocatable :: ypos(:, :), yneg(:, :)  !! (n, omax)
      real(dp), allocatable :: bblk(:, :, :)           !! (n, n, qb): one unpacked aux block
      real(dp), allocatable :: zpos(:, :), zneg(:, :)  !! (n, omax*qb): the half-transform
      real(dp), allocatable :: kpos(:, :), kneg(:, :)  !! (n, n): K's two factors, accumulated over blocks
      type(step_table_t) :: tbl
   contains
      procedure :: setup => df_jk_setup
      procedure :: build => df_jk_build
      procedure :: release => df_jk_release
      procedure :: report => df_jk_report
      procedure :: dims => df_jk_dims
      procedure :: packed_b => df_jk_packed_b
      procedure :: cholesky => df_jk_chol
      procedure :: is_ready => df_jk_is_ready
   end type df_jk_t

contains

   !
   ! `b` the orbital basis, `aux` the fitting basis (def2-universal-jkfit by
   ! default, see the CLI/C API), `pl` the orbital pair list, `lp` the
   ! linear-algebra handles (cuBLAS/cuSOLVER under OpenACC). `nocc_hint`
   ! sizes the density-factor workspace up front (the expected occupied
   ! count; `build` grows it if a density ever needs more factors than
   ! that, so a wrong guess costs a reallocation, not a wrong answer).
   ! `aux_block` is the (mu nu|P) build's own block size, in auxiliary
   ! functions; independent of `qb`, the per-ITERATION K block size, which
   ! is chosen here from the device's free memory (see `choose_qb`).
   !
   ! Idempotent: a second call on an object that already built its `bpk`
   ! returns immediately, exactly `trc_linalg_t%init`'s "created once" rule
   ! -- the object is meant to survive across more than one `trc_scf_run`
   ! on the same geometry (trc_capi.F90's handle), and the caller need not
   ! track whether this is the first call or the fifth.
   !
   subroutine df_jk_setup(this, b, aux, pl, lp, err, nocc_hint, aux_block)
      class(df_jk_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_linalg_t), intent(inout) :: lp
      type(error_t), intent(out) :: err
      integer, intent(in), optional :: nocc_hint
      integer, intent(in), optional :: aux_block

      type(trc_basis_t) :: blk
      real(dp), allocatable :: jmc(:, :), tens(:, :, :), tsph(:, :, :)
      integer :: n, naux, naux_c, npack, ablk, omax0
      integer :: s0, s1, ish, p0, p0s, np, nps, mu, nu, q
      type(timer_type) :: tc

      call err%clear()
      if (this%ready) return
      ! The orbital basis `b` is only ever read here at its CARTESIAN size
      ! (`n = b%nao`, never `b%nao_sph`) -- `d_c`/`gmat_c` in `build`
      ! below are the Cartesian-space scratch trc_scf_fock.F90's `to_cart`/
      ! `to_sph` fold the SCF's own (spherical, if `b%spherical`) matrices
      ! through at the call site -- so nothing here needs `b%spherical`
      ! at all, and a spherical orbital basis is not refused.
      if (maxval(b%sh_l) > DF_LMAX) then
         call err%set(ERROR_VALIDATION, "trc_scf_df: the orbital basis carries angular momentum "// &
                      "above what the density-fitting kernels were generated for (gen_df.py --lmax)")
         return
      end if
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         call err%set(ERROR_VALIDATION, "trc_scf_df: the auxiliary basis carries angular momentum "// &
                      "above what the density-fitting kernels were generated for (gen_df.py --lmax-aux)")
         return
      end if

      n = b%nao
      ! The auxiliary basis follows ITS OWN `spherical` tag, independent
      ! of the orbital basis's -- trc_capi.F90's `basis_json`/
      ! `trc_set_jkfit_json` set it to match the orbital basis by
      ! construction (pyscf's own `auxmol.cart = mol.cart` rule), but
      ! this object only ever reads `aux%spherical` itself, same as every
      ! other aux-basis-aware piece below. `naux_c`, the CARTESIAN count,
      ! sizes the (mu nu|P) build's own kernel calls and block walk
      ! (always Cartesian, like every kernel); `naux`, spherical-folded
      ! when `aux%spherical`, is every OTHER dimension below -- the
      ! metric, its Cholesky factor, `bpk`'s own column count, and (via
      ! `this%naux`) what a chained RI-MP2 step checks against its own.
      naux_c = aux%nao
      naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      npack = n*(n + 1)/2
      this%n = n
      this%naux = naux
      this%npack = npack

      ! `lp`'s own eigensolver workspace: a caller that already `%init`'d
      ! it (trc_capi.F90's own `cx%la%init(n)`, at the SCF's own space --
      ! `b%nao_sph` when `b%spherical`) sized it for a SMALLER n than this
      ! object's own density-factor `syev` (build()'s own `this%dsc`,
      ! Cartesian, `n = b%nao`) ever runs at: every OTHER caller already
      ! passes a Cartesian `n` (RI-MP2's own `trc_rimp2_run`/
      ! `trc_urimp2_run`, or `--df` before a spherical orbital basis was
      ! supported, when `b%nao_sph` and `b%nao` happened to always be
      ! equal), so this went unnoticed until a spherical `--df` run first
      ! made the two sizes differ -- `trc_linalg_t%init` only ever GROWS
      ! (`grow_workspace`, never shrinks), so calling it here with the
      ! size `build` actually needs is always safe and a no-op when a
      ! caller already sized `lp` this large or larger.
      call lp%init(n)

      call dbg_mem("before B")
      allocate (this%bpk(npack, naux))
      !$acc enter data create(this%bpk)
      call dbg_mem("after B")
      call this%tbl%init("dfjk", NAMES, bytes_steps=[ST_JGAM, ST_JBLD, ST_UNPK])

      ! --- the metric and its Cholesky factor -------------------------------
      ! Stored in `this%chol`, not a local metric: the energy build never
      ! needs it again after the trsm below absorbs it into B, but the
      ! gradient (trc_scf_grad_df.F90) does -- see this type's own `chol`
      ! field comment. Same lifetime as `bpk` from here on.
      ! M = (P|Q) is built Cartesian (the kernel's only space) and folded,
      ! M_sph = T^T M_cart T, BEFORE the Cholesky factor -- the same
      ! `sph_cart_to_sph` the SCF driver folds S/H/D/F through, with the
      ! auxiliary basis standing in for the orbital one it is normally
      ! called on (trc_sph.F90's own transforms take any `trc_basis_t`).
      allocate (this%chol(naux, naux))
      !$acc enter data create(this%chol)
      call dbg_mem("after metric alloc")
      call this%tbl%step_on(lp, ST_METRIC, tc)
      if (aux%spherical) then
         allocate (jmc(naux_c, naux_c))
         !$acc enter data create(jmc)
         call trc_df_2c(aux, jmc)
         call sph_metric_cart_to_sph(aux, lp, jmc, this%chol)
         !$acc exit data delete(jmc)
         deallocate (jmc)
      else
         call trc_df_2c(aux, this%chol)
      end if
      call this%tbl%step_off(lp, ST_METRIC, tc)
      call this%tbl%step_on(lp, ST_CHOL, tc)
      call lp%potrf(naux, this%chol, naux)
      call this%tbl%step_off(lp, ST_CHOL, tc)
      call this%tbl%add_flop(ST_CHOL, int(naux, int64)**3/3)
      call dbg_mem("after Cholesky")

      ! The setup block size, queried HERE -- after bpk AND the metric's
      ! Cholesky factor (naux*naux*8 bytes, ~226 MB even at w40's
      ! naux=5320, but Cholesky's own cuSOLVER workspace is not that small
      ! and is still resident at this point too) are both already spoken
      ! for, not back when only bpk was -- because that earlier query (one
      ! round's own bug, caught by the local-memory fix above freeing
      ! enough headroom for the naux/4 default to pass a free-memory check
      ! that did not yet know about the Cholesky step's own footprint) let
      ! `choose_setup_block` approve an ablk whose tensor then OOM'd against
      ! what the metric and Cholesky had since consumed. naux/4 at a time (the old,
      ! size-blind default) allocates a dense (n, n, naux/4) tensor for the
      ! (mu nu|P) build that alone can exceed what is left after bpk+the
      ! metric+Cholesky: w40's n=1000, naux=5320 asks for another 10.6 GB. Free-
      ! memory-aware, same pattern as the per-iteration K block
      ! (`choose_qb`) and RI-MP2's own `k_group`.
      ablk = choose_setup_block(n, naux_c)
      if (present(aux_block)) ablk = max(aux_block, 1)
      call dbg_mem("ablk chosen (ablk="//itoa_dbg(ablk)//")")

      ! --- (mu nu|P), one auxiliary block at a time, packed as it lands ----
      ! `p0`/`np` below are the block's EFFECTIVE (spherical, when
      ! `aux%spherical`) offset/width -- what `bpk`'s own column index
      ! runs over -- kept distinct from the block's CARTESIAN `p0`/`np`
      ! only where the two differ, read out of `blk` after `subset` cuts
      ! it (never across a shell, so `blk%nao_sph` is exactly this
      ! block's own spherical width).
      s0 = 1
      do while (s0 <= aux%nshell)
         s1 = s0
         do ish = s0 + 1, aux%nshell
            if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > ablk) exit
            s1 = ish
         end do
         p0 = aux%sh_ao(s0)
         np = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0
         call aux%subset(s0, s1, blk)
         call blk%to_device()
         allocate (tens(n, n, np))
         !$acc enter data create(tens)
         if (s0 == 1) call dbg_mem("after first (mn|P) block alloc (np="//itoa_dbg(np)//")")
         call this%tbl%step_on(lp, ST_3C, tc)
         call trc_df_3c(b, pl, blk, tens)
         call this%tbl%step_off(lp, ST_3C, tc)
         call this%tbl%step_on(lp, ST_PACK, tc)
         if (aux%spherical) then
            p0s = aux%sh_ao_sph(s0)
            nps = blk%nao_sph
            allocate (tsph(n, n, nps))
            !$acc enter data create(tsph)
            call sph_aux3c_to_sph(blk, n, tens, tsph, lp)
            do concurrent(q=1:nps, mu=1:n, nu=1:n, nu <= mu)
               this%bpk(pack_idx(mu, nu), p0s - 1 + q) = tsph(mu, nu, q)
            end do
            !$acc exit data delete(tsph)
            deallocate (tsph)
         else
            do concurrent(q=1:np, mu=1:n, nu=1:n, nu <= mu)
               this%bpk(pack_idx(mu, nu), p0 - 1 + q) = tens(mu, nu, q)
            end do
         end if
         call this%tbl%step_off(lp, ST_PACK, tc)
         !$acc exit data delete(tens)
         deallocate (tens)
         call blk%release()
         s0 = s1 + 1
      end do

      call dbg_mem("after (mn|P) + pack, before trsm")
      ! --- B = T L^-T: one triangular solve, every packed row at once ------
      ! `this%chol` is READ here, not consumed: trsm solves against it and
      ! leaves it untouched, so -- unlike the old local `jm`, freed right
      ! after this call -- it stays resident for the gradient to reuse.
      call this%tbl%step_on(lp, ST_TRSM, tc)
      call lp%trsm("R", "L", "T", npack, naux, 1.0_dp, this%chol, naux, this%bpk, npack)
      call this%tbl%step_off(lp, ST_TRSM, tc)
      call this%tbl%add_flop(ST_TRSM, int(npack, int64)*int(naux, int64)*int(naux, int64))
      call dbg_mem("after trsm")

      ! --- fixed-size per-iteration scratch, and the density-factor one ----
      allocate (this%dpk(npack), this%jpk(npack), this%gq(naux))
      allocate (this%dsc(n, n), this%evals(n))
      allocate (this%kpos(n, n), this%kneg(n, n))
      omax0 = max(n/4, 1)
      if (present(nocc_hint)) omax0 = max(nocc_hint, 1)
      this%qb = choose_qb(n, naux, omax0)
      call dbg_mem("qb chosen (qb="//itoa_dbg(this%qb)//", omax0="//itoa_dbg(omax0)//")")
      allocate (this%bblk(n, n, this%qb))
      !$acc enter data create(this%dpk, this%jpk, this%gq, this%dsc, this%evals, this%kpos, this%kneg, &
      !$acc                    this%bblk)
      call dbg_mem("after fixed scratch + bblk")
      call grow_yz(this, omax0)
      call dbg_mem("after ypos/yneg/zpos/zneg (setup done)")

      this%ready = .true.
   end subroutine df_jk_setup

   !
   ! G = j_scale J[D] - (k_scale/2) K[D], D the total density: the same convention
   ! `trc_eri_t%fock_resident` returns (trc_scf_fock.F90's build_g branches
   ! on `cfg%df` into this instead of the four-centre call, never both).
   ! `gmat_c` is OVERWRITTEN, not accumulated into -- the caller has
   ! nothing else to add for the restricted case `--df` is scoped to (P1).
   !
   ! K is skipped entirely when `k_scale == 0` (a pure functional): J is
   ! two GEMVs and nothing else runs, which is the point of fitting both
   ! under one flag rather than offering a DF-J-only mode.
   !
   subroutine df_jk_build(this, lp, d_c, gmat_c, j_scale, k_scale)
      class(df_jk_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: d_c(this%n, this%n)
      real(dp), intent(inout) :: gmat_c(this%n, this%n)
      real(dp), intent(in) :: j_scale, k_scale
      integer :: n, naux, npack
      integer :: mu, nu, k, q0, qb_cur, qq, i
      integer(int64) :: t, n8
      integer :: npos, nneg, ktot, ib, r0, r1, nb_panel
      integer, parameter :: NB_KPANEL = 8
      !! K's accumulation (below) is split into this many column panels,
      !! each a GEMM over its own diagonal block down to row n -- only the
      !! lower triangle, block by block -- rather than one GEMM over the
      !! whole dense (redundant) square. More panels means less redundant
      !! work (the executed/full-dense ratio is (NB_KPANEL+1)/(2*NB_KPANEL)
      !! -- 0.5625 at 8) but a smaller N per GEMM; 8 keeps every panel's
      !! N >= 125 at w40 (n=1000) while buying back most of the savings a
      !! much larger panel count would.
      real(dp) :: eval_cut
      type(timer_type) :: tc

      n = this%n
      n8 = int(n, int64)
      naux = this%naux
      npack = this%npack
      nb_panel = (n + NB_KPANEL - 1)/NB_KPANEL

      ! --- J: two GEMVs on the packed storage -------------------------------
      call this%tbl%step_on(lp, ST_JGAM, tc)
      do concurrent(mu=1:n, nu=1:n, nu <= mu)
         this%dpk(pack_idx(mu, nu)) = merge(1.0_dp, 2.0_dp, mu == nu)*d_c(mu, nu)
      end do
      call lp%gemv("T", npack, naux, 1.0_dp, this%bpk, npack, this%dpk, 0.0_dp, this%gq)
      call this%tbl%step_off(lp, ST_JGAM, tc)
      call this%tbl%add_bytes(ST_JGAM, 8_int64*int(npack, int64)*int(naux, int64))

      call this%tbl%step_on(lp, ST_JBLD, tc)
      call lp%gemv("N", npack, naux, 1.0_dp, this%bpk, npack, this%gq, 0.0_dp, this%jpk)
      do concurrent(mu=1:n, nu=1:n)
         gmat_c(mu, nu) = j_scale*this%jpk(pack_idx(max(mu, nu), min(mu, nu)))
      end do
      call this%tbl%step_off(lp, ST_JBLD, tc)
      call this%tbl%add_bytes(ST_JBLD, 8_int64*int(npack, int64)*int(naux, int64))

      if (k_scale == 0.0_dp) return

      ! --- the density factor: D = Y_pos Y_pos^T - Y_neg Y_neg^T ------------
      call this%tbl%step_on(lp, ST_DFAC, tc)
      do concurrent(mu=1:n, nu=1:n)
         this%dsc(mu, nu) = d_c(mu, nu)
      end do
      call lp%syev(this%dsc, this%evals)
      !$acc update self(this%evals)
      eval_cut = EVAL_REL_TOL*maxval(abs(this%evals))
      npos = 0
      nneg = 0
      do k = 1, n
         if (this%evals(k) > eval_cut) npos = npos + 1
         if (this%evals(k) < -eval_cut) nneg = nneg + 1
      end do
      call dbg_rank(npos, nneg)
      if (max(npos, nneg) > this%omax) call grow_yz(this, max(npos, nneg))
      do concurrent(i=1:n, k=1:npos)
         this%ypos(i, k) = this%dsc(i, n - npos + k)*sqrt(this%evals(n - npos + k))
      end do
      do concurrent(i=1:n, k=1:nneg)
         this%yneg(i, k) = this%dsc(i, k)*sqrt(-this%evals(k))
      end do
      call this%tbl%step_off(lp, ST_DFAC, tc)

      do concurrent(mu=1:n, nu=1:n)
         this%kpos(mu, nu) = 0.0_dp
         this%kneg(mu, nu) = 0.0_dp
      end do

      ! --- K, one auxiliary block at a time ---------------------------------
      q0 = 1
      do while (q0 <= naux)
         qb_cur = min(this%qb, naux - q0 + 1)

         ! One flat index over the whole (n, n, qb_cur) block, both
         ! triangles read from the one packed copy. The three-index forms
         ! (every ordering measured: 89-118 GB/s) ran their q loop
         ! sequentially inside each thread (nvfortran -Minfo), and the
         ! staging alternative -- triangular copy into a second (n, n, qb)
         ! buffer, then one cuBLAS geam per qq to symmetrise -- reached
         ! only 134 GB/s at w40, launch-bound at 8 MB per call. The DF
         ! gradient's own unpack (trc_scf_grad_df.F90) is the same loop.
         call this%tbl%step_on(lp, ST_UNPK, tc)
         do concurrent(t=0:n8*n8*int(qb_cur, int64) - 1) local(mu, nu, qq)
            mu = int(mod(t, n8)) + 1
            nu = int(mod(t/n8, n8)) + 1
            qq = int(t/(n8*n8)) + 1
            this%bblk(mu, nu, qq) = this%bpk(pack_idx(max(mu, nu), min(mu, nu)), q0 - 1 + qq)
         end do
         call this%tbl%step_off(lp, ST_UNPK, tc)
         call this%tbl%add_bytes(ST_UNPK, 8_int64*int(n, int64)*int(n, int64)*int(qb_cur, int64))

         call this%tbl%step_on(lp, ST_KHLF, tc)
         if (npos > 0) then
            call lp%gemm_strided("N", "N", n, npos, n, 1.0_dp, this%bblk, n, int(n, int64)*int(n, int64), &
                                 this%ypos, n, 0_int64, 0.0_dp, this%zpos, n, &
                                 int(n, int64)*int(npos, int64), qb_cur)
         end if
         if (nneg > 0) then
            call lp%gemm_strided("N", "N", n, nneg, n, 1.0_dp, this%bblk, n, int(n, int64)*int(n, int64), &
                                 this%yneg, n, 0_int64, 0.0_dp, this%zneg, n, &
                                 int(n, int64)*int(nneg, int64), qb_cur)
         end if
         call this%tbl%step_off(lp, ST_KHLF, tc)
         call this%tbl%add_flop(ST_KHLF, 2_int64*int(n, int64)**2*int(npos + nneg, int64)*int(qb_cur, int64))

         ! K += Z Z^T as a plain GEMM, not cuBLAS's own dsyrk (dsyrk only
         ! ever touches the upper triangle, but at a fraction of dgemm's
         ! rate -- measured at w40 cc-pVDZ: dsyrk 2947 GFLOP/s, ~38% of
         ! this V100's FP64 peak, a GEMM at 90%+) -- and split into
         ! NB_KPANEL column panels, each a GEMM over its own diagonal
         ! block down to row n, so the REDUNDANT half of a full dense
         ! GEMM's work (computing the upper triangle too) is mostly not
         ! done at all rather than done cheaply: executed flops fall to
         ! (NB_KPANEL+1)/(2*NB_KPANEL) of the full dense product (0.5625
         ! at 8), while every panel is still GEMM-shaped (large K, M
         ! shrinking panel to panel, N fixed at the panel width) so the
         ! per-flop rate stays close to the dense GEMM's own. `zpos(r0,1)`/
         ! `kpos(r0,r0)` passed to an assumed-size `a(lda,*)` dummy is
         ! ordinary Fortran sequence association -- the same idiom any
         ! BLAS call into a sub-matrix uses -- and OpenACC's
         ! `host_data use_device` resolves it by offset within the whole
         ! array's own device registration, the same as any other strided
         ! sub-array passed to `lp%gemm` elsewhere in this file.
         ! `add_flop` below is EXECUTED flops (2*M*N*K per panel, standard
         ! GEMM count); the USEFUL flops -- the lower triangle only, what
         ! this was computing before NB_KPANEL>1 existed -- are n*(n+1)*k,
         ! essentially half of a full dense GEMM's 2*n*n*k, so dividing
         ! that by the panelled scheme's own executed total gives the
         ! fraction of this round's savings that NB_KPANEL actually buys.
         call this%tbl%step_on(lp, ST_KSYRK, tc)
         if (npos > 0) then
            ktot = npos*qb_cur
            do ib = 1, NB_KPANEL
               r0 = 1 + (ib - 1)*nb_panel
               r1 = min(ib*nb_panel, n)
               if (r0 > n) exit
               call lp%gemm("N", "T", n - r0 + 1, r1 - r0 + 1, ktot, 1.0_dp, this%zpos(r0, 1), n, &
                            this%zpos(r0, 1), n, 1.0_dp, this%kpos(r0, r0), n)
               call this%tbl%add_flop(ST_KSYRK, &
                                      2_int64*int(n - r0 + 1, int64)*int(r1 - r0 + 1, int64)*int(ktot, int64))
            end do
         end if
         if (nneg > 0) then
            ktot = nneg*qb_cur
            do ib = 1, NB_KPANEL
               r0 = 1 + (ib - 1)*nb_panel
               r1 = min(ib*nb_panel, n)
               if (r0 > n) exit
               call lp%gemm("N", "T", n - r0 + 1, r1 - r0 + 1, ktot, 1.0_dp, this%zneg(r0, 1), n, &
                            this%zneg(r0, 1), n, 1.0_dp, this%kneg(r0, r0), n)
               call this%tbl%add_flop(ST_KSYRK, &
                                      2_int64*int(n - r0 + 1, int64)*int(r1 - r0 + 1, int64)*int(ktot, int64))
            end do
         end if
         call this%tbl%step_off(lp, ST_KSYRK, tc)

         q0 = q0 + qb_cur
      end do

      ! The panelled accumulation above only ever wrote the lower
      ! triangle (row >= column) of kpos/kneg -- mirror it before this
      ! object's own K result (both triangles, since gmat_c is dense) is
      ! read out.
      do concurrent(mu=1:n, nu=1:n, nu > mu)
         this%kpos(mu, nu) = this%kpos(nu, mu)
         this%kneg(mu, nu) = this%kneg(nu, mu)
      end do

      do concurrent(mu=1:n, nu=1:n)
         gmat_c(mu, nu) = gmat_c(mu, nu) - 0.5_dp*k_scale*(this%kpos(mu, nu) - this%kneg(mu, nu))
      end do
   end subroutine df_jk_build

   pure subroutine df_jk_dims(this, n, naux, npack)
      !! This object's dimensions, and whether `setup` has already fitted `bpk`
      !! -- what a chained RI-MP2 step checks before reusing it instead of
      !! redoing the metric/Cholesky/3c/trsm on the same auxiliary basis (see
      !! `packed_b` below, and trc_rimp2_driver.F90's own `dfjk` argument).
      class(df_jk_t), intent(in) :: this
      integer, intent(out) :: n, naux, npack
      n = this%n
      naux = this%naux
      npack = this%npack
   end subroutine df_jk_dims

   pure logical function df_jk_is_ready(this) result(r)
      class(df_jk_t), intent(in) :: this
      r = this%ready
   end function df_jk_is_ready

   function df_jk_packed_b(this) result(p)
      !! B itself, packed (mu >= nu, p = mu(mu-1)/2 + nu) and already L^-T
      !! applied -- (npack, naux), device-resident. A chained RI-MP2 step
      !! reads this directly instead of rebuilding it: L^-T acts only on the
      !! auxiliary index and the occupied/virtual half-transform only on
      !! (mu, nu), so the two commute, and the half-transform run on THIS is
      !! exactly RI-MP2's own B, no further triangular solve needed. Null
      !! unless `is_ready()`.
      class(df_jk_t), intent(in), target :: this
      real(dp), pointer :: p(:, :)
      if (this%ready) then
         p => this%bpk
      else
         p => null()
      end if
   end function df_jk_packed_b

   function df_jk_chol(this) result(p)
      !! L, the metric's Cholesky factor (M = (P|Q) = L L^T), lower-triangular,
      !! (naux, naux), device-resident -- see this type's `chol` field comment.
      !! The DF gradient (trc_scf_grad_df.F90) solves against this directly
      !! (one triangular solve, `lp%trsm`/`lp%trsv` with "L","T") to turn the
      !! Cholesky-space coefficients `build` already computes into the actual
      !! fitting coefficients its envelope-theorem derivative is stationary in,
      !! with no metric or Cholesky recomputation. Null unless `is_ready()`.
      class(df_jk_t), intent(in), target :: this
      real(dp), pointer :: p(:, :)
      if (this%ready) then
         p => this%chol
      else
         p => null()
      end if
   end function df_jk_chol

   subroutine df_jk_release(this)
      class(df_jk_t), intent(inout) :: this
      if (.not. this%ready) return
      !$acc exit data delete(this%bpk, this%chol, this%dpk, this%jpk, this%gq, this%dsc, this%evals, &
      !$acc                  this%kpos, this%kneg, this%bblk)
      if (allocated(this%ypos)) then
         !$acc exit data delete(this%ypos, this%yneg, this%zpos, this%zneg)
      end if
      deallocate (this%bpk, this%chol, this%dpk, this%jpk, this%gq, this%dsc, this%evals, this%kpos, this%kneg, &
                 this%bblk)
      if (allocated(this%ypos)) deallocate (this%ypos, this%yneg, this%zpos, this%zneg)
      this%ready = .false.
      this%omax = 0
   end subroutine df_jk_release

   !! The setup/per-iteration table: time, GFLOP(/s) or GB/s per step, and
   ! one "all dense algebra" row over the GEMM/factorisation-shaped steps
   ! (Cholesky, trsm, the half-transform, the K panels) -- J's two GEMVs and the
   ! density factor are bandwidth- or eigensolver-bound and left out of it,
   ! exactly as RI-MP2 leaves its own screening-dependent integral step out.
   subroutine df_jk_report(this, t_total, rank, nranks)
      class(df_jk_t), intent(in) :: this
      real(dp), intent(in) :: t_total
      integer, intent(in) :: rank, nranks
      logical :: mask(NSTEP)
      mask = .false.
      mask([ST_CHOL, ST_TRSM, ST_KHLF, ST_KSYRK]) = .true.
      call this%tbl%report(t_total, rank, nranks, subtotal_mask=mask, subtotal_label="all dense algebra")
   end subroutine df_jk_report

   !
   ! The density-factor workspace: Y_pos/Y_neg (n, omax) and the half-
   ! transform's Z_pos/Z_neg (n, omax*qb). Grows, never shrinks -- the same
   ! rule `trc_linalg_t%grow_workspace` follows -- so a density that
   ! occasionally needs more factors than the setup hint costs one
   ! reallocation, not a wrong answer.
   !
   subroutine grow_yz(this, want_o)
      type(df_jk_t), intent(inout) :: this
      integer, intent(in) :: want_o
      if (want_o <= this%omax .and. allocated(this%ypos)) return
      if (allocated(this%ypos)) then
         !$acc exit data delete(this%ypos, this%yneg, this%zpos, this%zneg)
         deallocate (this%ypos, this%yneg, this%zpos, this%zneg)
      end if
      this%omax = max(want_o, 1)
      allocate (this%ypos(this%n, this%omax), this%yneg(this%n, this%omax))
      allocate (this%zpos(this%n, this%omax*this%qb), this%zneg(this%n, this%omax*this%qb))
      !$acc enter data create(this%ypos, this%yneg, this%zpos, this%zneg)
   end subroutine grow_yz

   !
   ! Auxiliary functions per K block: enough that k = omax*qb (the K
   ! half-transform/accumulation GEMMs' own contraction length) reaches
   ! `K_TARGET`, since a GEMM only clears ~90% of FP64 peak once its
   ! contraction length is in the thousands -- QB_DEFAULT alone cleared that
   ! for omax around a hundred or more (every case this was tuned against),
   ! but not for a molecule with a genuinely small occupied/virtual count
   ! (`omax` is read off D's eigenvalues at every build, not a size this
   ! function controls). Capped, as before, by whatever the device's free
   ! memory (the unpacked block plus the two half-transform factors, `omax`
   ! columns each) will hold with half a gigabyte kept back for the
   ! libraries' own workspace -- the same pattern as RI-MP2's `k_group`. A
   ! fixed default on the host build, where there is no device memory query
   ! and the arrays are small next to host RAM regardless.
   !
   integer function choose_qb(n, naux, omax) result(qb)
      integer, intent(in) :: n, naux, omax
      integer, parameter :: QB_DEFAULT = 64
      integer, parameter :: K_TARGET = 4096
      integer :: qb_target
#ifdef TRC_DEVICE
      integer(int64) :: free, per_q, fixed
#endif
      qb_target = max(QB_DEFAULT, (K_TARGET + max(omax, 1) - 1)/max(omax, 1))
#ifdef TRC_DEVICE
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_q = 8_int64*(int(n, int64)*int(n, int64) + 2_int64*int(n, int64)*int(max(omax, 1), int64))
      fixed = 512_int64*1024*1024
      qb = int(min(int(qb_target, int64), max(1_int64, (free - fixed)/max(per_q, 1_int64))))
#else
      qb = qb_target
#endif
      qb = max(1, min(qb, naux))
   end function choose_qb

   !
   ! Auxiliary functions per (mu nu|P) SETUP block: as many as naux/4 (the
   ! old, size-blind default) unless the device's free memory -- AFTER
   ! `bpk` (p x naux, the largest single allocation this object ever
   ! makes) is already resident -- would not hold the dense (n, n, ablk)
   ! scratch tensor that block's build needs, in which case as many as
   ! fit with half a gigabyte kept back. Same pattern as `choose_qb`.
   !
   !
   ! TERCO_DF_MEM_DEBUG=1: one line per stage, the device's free memory
   ! right then, flushed immediately -- a device OOM kills the process
   ! before an unflushed `print` ever reaches the pipe, which is what
   ! made an earlier attempt at this (a plain `print` with no flush)
   ! print nothing at all for the stage that actually failed. Debug only;
   ! a no-op allocation-free query otherwise, on the host build too
   ! (prints 0, which is honest: there is no device to ask).
   !
   pure function itoa_dbg(i) result(t)
      !! Plain integer-to-string, for `dbg_mem`'s own labels only.
      integer, intent(in) :: i
      character(len=16) :: t
      write (t, "(i0)") i
   end function itoa_dbg

   subroutine dbg_mem(label)
      character(len=*), intent(in) :: label
      character(len=8) :: v
      integer(int64) :: free
      call get_environment_variable("TERCO_DF_MEM_DEBUG", v)
      if (trim(v) /= "1") return
#ifdef TRC_DEVICE
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
#else
      free = 0_int64
#endif
      block
         character(len=1024) :: logline
         write (logline, "(a,a,a,f10.3,a)") "trc_scf_df: ", label, ": free ", real(free, dp)/1024.0_dp**3, " GiB"
         call global_logger%info(trim(logline))
      end block
      flush (output_unit)
   end subroutine dbg_mem

   !
   ! TERCO_DF_RANK_DEBUG=1: the density factor's rank (npos, nneg) every
   ! `build` call, one line each -- the K half-transform/accumulation's
   ! own cost is linear in npos+nneg, so this is what the relative cut
   ! above (EVAL_REL_TOL) actually buys, iteration by iteration, next to
   ! the old absolute 1e-12 cut. Debug only; a no-op otherwise.
   !
   subroutine dbg_rank(npos, nneg)
      integer, intent(in) :: npos, nneg
      character(len=8) :: v
      call get_environment_variable("TERCO_DF_RANK_DEBUG", v)
      if (trim(v) /= "1") return
      block
         character(len=1024) :: logline
         write (logline, "(a,i0,a,i0)") "trc_scf_df: density factor rank: npos=", npos, " nneg=", nneg
         call global_logger%info(trim(logline))
      end block
      flush (output_unit)
   end subroutine dbg_rank

   integer function choose_setup_block(n, naux) result(ablk)
      integer, intent(in) :: n, naux
#ifdef TRC_DEVICE
      integer(int64) :: free, per_q, fixed
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_q = 8_int64*int(n, int64)*int(n, int64)
      fixed = 512_int64*1024*1024
      ablk = int(min(int(max(naux/4, 1), int64), max(1_int64, (free - fixed)/max(per_q, 1_int64))))
#else
      ablk = max(naux/4, 1)
#endif
      ablk = max(1, min(ablk, naux))
   end function choose_setup_block

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   !! idx(mu, nu) = mu(mu-1)/2 + nu, mu >= nu, 1-based: the packed row a
   !! (mu, nu) pair (either order) lives at. Module-level and pure (not
   !! type-bound) so a `do concurrent` body may call it directly.
   pure integer function pack_idx(mu, nu) result(p)
      !$acc routine seq
      integer, intent(in) :: mu, nu
      p = mu*(mu - 1)/2 + nu
   end function pack_idx

end module trc_scf_df
