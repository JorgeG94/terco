!
! Density-fitted J and K for a restricted (closed-shell) reference,
! Cartesian orbital basis, resident on the device.
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
! eigenvalue with |n_k| > 1e-12: Y_pos = U sqrt(n_k) for the positive ones,
! Y_neg = U sqrt(|n_k|) for the negative ones (K is linear in D, so
! K[D] = K[D_pos] - K[D_neg] exactly). Then, one auxiliary block at a time:
!
!     unpack B_q to a dense (n, n, qb) scratch (both triangles filled)
!     X_q = B_q Y            one strided GEMM batch over the block, giving
!                             Z(mu, (i,q)), n x (o qb), per factor
!     K += Z Z^T             `trc_linalg%syrk`, accumulated over blocks
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
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
#ifdef _OPENACC
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
                                                   "unpack B", "K half-transform", "K syrk"]
   real(dp), parameter :: EVAL_TOL = 1.0e-12_dp

   type :: df_jk_t
      private
      logical :: ready = .false.
      integer :: n = 0, naux = 0, npack = 0
      integer :: qb = 1       !! auxiliary functions per K block, chosen at setup
      integer :: omax = 0     !! current capacity of the density-factor workspace
      real(dp), allocatable :: bpk(:, :)              !! (npack, naux): B, packed mu >= nu
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
      real(dp), allocatable :: jm(:, :), tens(:, :, :)
      integer :: n, naux, npack, ablk, omax0
      integer :: s0, s1, ish, p0, np, mu, nu, q
      real(dp) :: tc

      call err%clear()
      if (this%ready) return
      if (b%spherical) then
         call err%set(ERROR_VALIDATION, "trc_scf_df: the orbital basis is spherical; DF-JK only "// &
                      "supports a Cartesian basis today")
         return
      end if
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
      naux = aux%nao
      npack = n*(n + 1)/2
      this%n = n
      this%naux = naux
      this%npack = npack

      allocate (this%bpk(npack, naux))
      !$acc enter data create(this%bpk)
      ! The setup block size, AFTER bpk's own ~p*naux*8 bytes are already
      ! spoken for (w40 cc-pVDZ: ~21 GB of a 32 GB V100) -- naux/4 at a
      ! time (the old, size-blind default) allocates a SECOND dense (n, n,
      ! naux/4) tensor for the (mu nu|P) build that alone can exceed what
      ! bpk left: w40's n=1000, naux=5320 would ask for another 10.6 GB,
      ! and the setup OOMs. Free-memory-aware, same pattern as the
      ! per-iteration K block (`choose_qb`) and RI-MP2's own `k_group`.
      ablk = choose_setup_block(n, naux)
      if (present(aux_block)) ablk = max(aux_block, 1)
      call this%tbl%init("dfjk", NAMES, bytes_steps=[ST_JGAM, ST_JBLD, ST_UNPK])

      ! --- the metric and its Cholesky factor -------------------------------
      allocate (jm(naux, naux))
      !$acc enter data create(jm)
      call this%tbl%step_on(lp, ST_METRIC, tc)
      call trc_df_2c(aux, jm)
      call this%tbl%step_off(lp, ST_METRIC, tc)
      call this%tbl%step_on(lp, ST_CHOL, tc)
      call lp%potrf(naux, jm, naux)
      call this%tbl%step_off(lp, ST_CHOL, tc)
      call this%tbl%add_flop(ST_CHOL, int(naux, int64)**3/3)

      ! --- (mu nu|P), one auxiliary block at a time, packed as it lands ----
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
         call this%tbl%step_on(lp, ST_3C, tc)
         call trc_df_3c(b, pl, blk, tens)
         call this%tbl%step_off(lp, ST_3C, tc)
         call this%tbl%step_on(lp, ST_PACK, tc)
         do concurrent(q=1:np, mu=1:n, nu=1:n, nu <= mu)
            this%bpk(pack_idx(mu, nu), p0 - 1 + q) = tens(mu, nu, q)
         end do
         call this%tbl%step_off(lp, ST_PACK, tc)
         !$acc exit data delete(tens)
         deallocate (tens)
         call blk%release()
         s0 = s1 + 1
      end do

      ! --- B = T L^-T: one triangular solve, every packed row at once ------
      call this%tbl%step_on(lp, ST_TRSM, tc)
      call lp%trsm("R", "L", "T", npack, naux, 1.0_dp, jm, naux, this%bpk, npack)
      call this%tbl%step_off(lp, ST_TRSM, tc)
      call this%tbl%add_flop(ST_TRSM, int(npack, int64)*int(naux, int64)*int(naux, int64))
      !$acc exit data delete(jm)
      deallocate (jm)

      ! --- fixed-size per-iteration scratch, and the density-factor one ----
      allocate (this%dpk(npack), this%jpk(npack), this%gq(naux))
      allocate (this%dsc(n, n), this%evals(n))
      allocate (this%kpos(n, n), this%kneg(n, n))
      omax0 = max(n/4, 1)
      if (present(nocc_hint)) omax0 = max(nocc_hint, 1)
      this%qb = choose_qb(n, naux, omax0)
      allocate (this%bblk(n, n, this%qb))
      !$acc enter data create(this%dpk, this%jpk, this%gq, this%dsc, this%evals, this%kpos, this%kneg, this%bblk)
      call grow_yz(this, omax0)

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
      integer :: mu, nu, k, q0, qb_cur, qq, muu, nuu, i
      integer :: npos, nneg
      real(dp) :: tc

      n = this%n
      naux = this%naux
      npack = this%npack

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
      npos = 0
      nneg = 0
      do k = 1, n
         if (this%evals(k) > EVAL_TOL) npos = npos + 1
         if (this%evals(k) < -EVAL_TOL) nneg = nneg + 1
      end do
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

         call this%tbl%step_on(lp, ST_UNPK, tc)
         do concurrent(qq=1:qb_cur, muu=1:n, nuu=1:n)
            this%bblk(muu, nuu, qq) = this%bpk(pack_idx(max(muu, nuu), min(muu, nuu)), q0 - 1 + qq)
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

         call this%tbl%step_on(lp, ST_KSYRK, tc)
         if (npos > 0) call lp%syrk("U", "N", n, npos*qb_cur, 1.0_dp, this%zpos, n, 1.0_dp, this%kpos, n)
         if (nneg > 0) call lp%syrk("U", "N", n, nneg*qb_cur, 1.0_dp, this%zneg, n, 1.0_dp, this%kneg, n)
         call this%tbl%step_off(lp, ST_KSYRK, tc)
         call this%tbl%add_flop(ST_KSYRK, int(n, int64)**2*int(npos + nneg, int64)*int(qb_cur, int64))

         q0 = q0 + qb_cur
      end do

      do concurrent(mu=1:n, nu=1:n)
         gmat_c(mu, nu) = gmat_c(mu, nu) - 0.5_dp*k_scale*(this%kpos(mu, nu) - this%kneg(mu, nu))
      end do
   end subroutine df_jk_build

   !! This object's dimensions, and whether `setup` has already fitted `bpk`
   !! -- what a chained RI-MP2 step checks before reusing it instead of
   !! redoing the metric/Cholesky/3c/trsm on the same auxiliary basis (see
   !! `packed_b` below, and trc_rimp2_driver.F90's own `dfjk` argument).
   pure subroutine df_jk_dims(this, n, naux, npack)
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

   !! B itself, packed (mu >= nu, p = mu(mu-1)/2 + nu) and already L^-T
   !! applied -- (npack, naux), device-resident. A chained RI-MP2 step
   !! reads this directly instead of rebuilding it: L^-T acts only on the
   !! auxiliary index and the occupied/virtual half-transform only on
   !! (mu, nu), so the two commute, and the half-transform run on THIS is
   !! exactly RI-MP2's own B, no further triangular solve needed. Null
   !! unless `is_ready()`.
   function df_jk_packed_b(this) result(p)
      class(df_jk_t), intent(in), target :: this
      real(dp), pointer :: p(:, :)
      if (this%ready) then
         p => this%bpk
      else
         p => null()
      end if
   end function df_jk_packed_b

   subroutine df_jk_release(this)
      class(df_jk_t), intent(inout) :: this
      if (.not. this%ready) return
      !$acc exit data delete(this%bpk, this%dpk, this%jpk, this%gq, this%dsc, this%evals, &
      !$acc                  this%kpos, this%kneg, this%bblk)
      if (allocated(this%ypos)) then
         !$acc exit data delete(this%ypos, this%yneg, this%zpos, this%zneg)
      end if
      deallocate (this%bpk, this%dpk, this%jpk, this%gq, this%dsc, this%evals, this%kpos, this%kneg, this%bblk)
      if (allocated(this%ypos)) deallocate (this%ypos, this%yneg, this%zpos, this%zneg)
      this%ready = .false.
      this%omax = 0
   end subroutine df_jk_release

   !! The setup/per-iteration table: time, GFLOP(/s) or GB/s per step, and
   ! one "all dense algebra" row over the GEMM/factorisation-shaped steps
   ! (Cholesky, trsm, the half-transform, syrk) -- J's two GEMVs and the
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
   ! Auxiliary functions per K block: `QB_DEFAULT`, or fewer if the device's
   ! free memory (the unpacked block plus the two half-transform factors,
   ! `omax` columns each) would not hold it with half a gigabyte kept back
   ! for the libraries' own workspace -- the same pattern as RI-MP2's
   ! `k_group`. A fixed default on the host build, where there is no device
   ! memory query and the arrays are small next to host RAM regardless.
   !
   integer function choose_qb(n, naux, omax) result(qb)
      integer, intent(in) :: n, naux, omax
      integer, parameter :: QB_DEFAULT = 64
#ifdef _OPENACC
      integer(int64) :: free, per_q, fixed
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_q = 8_int64*(int(n, int64)*int(n, int64) + 2_int64*int(n, int64)*int(max(omax, 1), int64))
      fixed = 512_int64*1024*1024
      qb = int(min(int(QB_DEFAULT, int64), max(1_int64, (free - fixed)/max(per_q, 1_int64))))
#else
      qb = QB_DEFAULT
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
   integer function choose_setup_block(n, naux) result(ablk)
      integer, intent(in) :: n, naux
#ifdef _OPENACC
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
