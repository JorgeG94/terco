!
! RI-MP2 for a closed-shell reference, on the device.
!
! With the auxiliary basis {P} and its metric J_PQ = (P|Q) = L L^T,
!
!     B_ia^Q = sum_P (ia|P) (L^-T)_PQ        so that  (ia|jb) = sum_Q B_ia^Q B_jb^Q
!
! and the correlation energy splits into its opposite- and same-spin parts,
!
!     E_os = sum_ijab (ia|jb)^2 / D_ijab
!     E_ss = sum_ijab [(ia|jb)^2 - (ia|jb)(ib|ja)] / D_ijab,   D = e_i + e_j - e_a - e_b
!
! E_MP2 = E_os + E_ss; SCS and SOS are weights on the two and cost nothing.
!
! THE STEPS
! ---------
!   1. (P|Q) by trc_df_2c, Cholesky on the device (trc_linalg%potrf).
!   2. (mu nu|P) by trc_df_3c for this rank's auxiliary shells, a block of
!      them at a time so the tensor never exceeds `aux_block` functions
!      deep; each block is an auxiliary basis of its own, built from a
!      shell range, and the kernels see nothing unusual.
!   3. (ia|P) = C_occ^T (mu nu|P) C_vir for every correlated i and this
!      rank's P, a fastest: one GEMM and one strided batch per auxiliary
!      block and destination rank. Then one exchange, after which each rank
!      holds X(ia, P) for its own occupied block and every P.
!   4. B = X L^-T: one triangular solve from the right, on this rank's rows.
!   5. K(ia, jb) = sum_Q B(ia, Q) B(jb, Q) for the pairs this rank owns, a
!      group of occupied i per GEMM; then a kernel over (j, a) summing b
!      with the denominators. Both (ia|jb) and (ib|ja) are in K. The pair
!      energy is symmetric -- swapping i with j and a with b leaves
!      (ia|jb)^2/D and (ia|jb)(ib|ja)/D as they were -- so each unordered
!      pair is computed once and an i /= j pair counted twice, which takes
!      the GEMMs, nearly all of RI-MP2, to about half the work.
!
! RANKS: B IS DISTRIBUTED
! -----------------------
! The correlated occupied orbitals are split into one contiguous block per
! rank, O_r, and rank r holds only B_r, the rows ia with i in O_r: nocc nvir
! naux / nranks doubles instead of all of them, which is what bounds the
! system on one device.
!
! The integrals are split the other way, by auxiliary shell: rank r builds
! (mu nu|P) only for P in P_r, a contiguous shell range holding about
! naux/nranks Cartesian functions, and transforms it for every occupied
! orbital, which is the same transform work as before. Built on every
! rank, they had been half of the run: ionic05 cc-pVTZ/cc-pVTZ-RIFIT on
! 16 A100s, 9.50 s of 19.18. That leaves rank r with X(ia, P in P_r) for
! every i. The block of rows for O_q goes to rank q, and the block for
! O_r is written straight into this rank's own X, so a one-rank run sends
! nothing and needs no extra memory. The exchange is nranks - 1 pairwise
! shifts. At step k rank r sends to r + k and receives from r - k, and
! what it receives is X(O_r, P_{r-k}), a contiguous column block of X.
! The send buffer is about one rank's share of B and is freed before the
! ring's two visiting buffers are allocated, so the peak memory does not
! change.
!
! The pairs then come from a ring. At step s = 0 .. nranks/2 rank r holds
! the block of rank r + s (its own at s = 0, then each step passes the held
! block one rank to the left) and computes the pairs (i in O_r, j in
! O_{r+s}):
!   s = 0             j >= i inside its own block;
!   0 < s < n/2       every pair of the two blocks;
!   s = n/2, n even   the blocks r and r + n/2 meet twice, once on each
!                     rank, so the lower rank takes the first half of its
!                     own rows and the higher rank the second half of the
!                     visiting block's: each pair once.
! Every unordered pair (i, j) is then computed by exactly one rank, and
! E_os and E_ss are summed over the ranks at the end. On one rank only
! s = 0 runs, which is the single-device algorithm.
!
! A CUDA-aware MPI is handed the device buffers directly, which between
! GPUs joined by NVLink moves a block in a fraction of the time; any other
! MPI gets them through the host (trc_mpi_device decides, at run time). This
! is the distribution of Stocks, Palethorpe and Barca, JCTC 2024, 20, 7503.
!
! Steps 3-5 are cuBLAS from OpenACC through trc_linalg, pic-blas and LAPACK
! on the host; the two kernels are `do concurrent`.
!
! Module trc_rimp2_driver, not trc_rimp2: the C entry is bound to the
! symbol trc_rimp2, and a binding label may not equal a module name
! (F2018 19.2); gfortran silently resolves the call to the wrong thing.
module trc_rimp2_driver
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use pic_mpi_lib, only: comm_t, allreduce, bcast, isend, recv, request_t, waitall, MPI_SUM, &
                          comm_isend_real_dp_array_n, comm_recv_real_dp_array_n
   use trc_mpi_device, only: trc_mpi_device_direct
   use pic_flop_recorder, only: flop_recorder_type
   use trc_step_table, only: step_table_t
   use trc_scf_df, only: df_jk_t
   use trc_sph, only: sph_metric_cart_to_sph, sph_mo_to_cart, sph_aux3c_to_sph
   implicit none
   private

   public :: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, rimp2_driver_mo_to_cart

   real(dp), parameter :: PAD_EPS = 1.0e6_dp
   !! `rimp2_driver_mo_to_cart`'s orbital energy for a padded, zero-coefficient
   !! "virtual" -- see that routine's own header for why finite (never
   !! literally infinite) is what keeps every denominator it enters finite
   !! rather than NaN, and why the padded terms are exactly zero anyway.

   !
   ! THE STEPS, AS TIMED, through the shared trc_step_table (src/core/
   ! trc_step_table.F90): each is device time to completion -- the clock is
   ! read after a device sync, since cuBLAS returns before its kernels end --
   ! and the floating-point operations it executes, so verbose output can say
   ! how close each step runs to the device's peak. The integral steps are
   ! time only: their operation counts are not GEMM-shaped and depend on the
   ! screening. Each step is also an NVTX range of the same name (a no-op
   ! unless built with TERCO_NVTX), opened and closed at the same device
   ! syncs as the timer, so on an Nsight Systems timeline a range covers its
   ! kernels and not only their launches.
   !
   integer, parameter :: NSTEP = 9
    !! number of algorithmic steps in the precodure
   integer, parameter :: ST_METRIC = 1, ST_CHOL = 2, ST_3C = 3, ST_TRANS = 4, ST_XCHG = 5, ST_TRSM = 6, &
                         ST_KGEMM = 7, ST_PAIR = 8, ST_RING = 9
    !! their map to an int
   character(len=18), parameter :: STEP_NAME(NSTEP) = [character(len=18) :: &
                                                       "metric (P|Q)", "Cholesky", "(mu nu|P)", &
                                                       "transform (ia|P)", "(ia|P) exchange", "B = X L^-T", &
                                                       "K_i GEMMs", "pair energies", "B exchange"]
    !! the names of the stages for FLOP counting

   !
   ! The largest message the ring sends in one call: MPI counts are 32-bit,
   ! and a block of B passes 2^31 doubles near w100 cc-pVDZ.
   !
   integer(int64), parameter :: MSG_CHUNK = 2_int64**30

   type :: trc_rimp2_result_t
    !! holds the rimp2 result
      real(dp) :: e_os = 0.0_dp     !! opposite-spin pair energy
      real(dp) :: e_ss = 0.0_dp     !! same-spin pair energy
      real(dp) :: e_corr = 0.0_dp   !! E_os + E_ss
      integer :: naux = 0, nocc = 0, nvir = 0  !! number of auxiliary, occupied, and virtual orbitals
      character(len=200) :: message = ""
   end type trc_rimp2_result_t

contains

   !
   ! `cmo` (nao, nao) and `eps` (nao) are the RHF orbitals, occupied first;
   ! `nocc` doubly occupied, `nfrozen` of them left out of the correlation.
   ! `pl` is the orbital basis's pair list. `la` optional, made here if
   ! absent; `comm` optional, one rank without it.
   !
   ! `dfjk`, optional: a density-fitted J/K object (trc_scf_df.F90) already
   ! fitted on an auxiliary basis the CALLER has checked is the same one
   ! as `aux` (same JSON path/contents -- trc_capi.F90's trc_run_rimp2 is
   ! the one caller that does; this routine only checks the dimensions
   ! agree, not that the basis itself is really the same one). When
   ! present, steps 1 (metric, Cholesky), 2 (the (mu nu|P) build) and 4
   ! (B = X L^-T) are skipped entirely: `dfjk%packed_b()` is already B
   ! (packed mu >= nu, L^-T applied), and the occupied/virtual half-
   ! transform of step 3 run on it directly IS RI-MP2's own B -- L^-T acts
   ! only on the auxiliary index, the half-transform only on (mu, nu), and
   ! the two commute. `aux` is still used for its shell structure (the
   ! block loop below only needs auxiliary function COUNTS and offsets,
   ! identical either way).
   !
   subroutine trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, nfrozen, aux_block, la, comm, verbose, dfjk)
      !! run a restricted ri-mp2 calculation
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      type(trc_rimp2_result_t), intent(out) :: res
      integer, intent(in), optional :: nfrozen, aux_block
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm
      logical, intent(in), optional :: verbose   !! one line per step with its wall time, rank 0
      type(df_jk_t), intent(in), optional, target :: dfjk

      type(timer_type) :: tw0, tw1
      logical :: talk, chained
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      type(trc_basis_t) :: aux_subset
      real(dp), pointer :: bchain(:, :)
      integer :: nao, naux, naux_c, nfz, no, nv, nov, nblk, s0, s1, p0, np, p0c, npc, ish
      integer :: i, rank, nranks, ng, ldk, i0, i1
      integer :: lo, nb, nmax, novr, step, nstep_ring, qv, cur, ncolmax, ldv
      integer :: ilo, ihi, jlo, jhi, jbase, jstart, ncol, roff, coff, jfirst, ncol_i
      integer :: n_chk, naux_chk, npack_chk, mu, nu, q
      integer, allocatable :: blo(:), bn(:), shlo(:), shhi(:), plo(:)
      integer(int64), allocatable :: soff(:)
      logical :: diag, direct
      real(dp), allocatable :: jm(:, :), jmc(:, :), cocc(:, :), cvir(:, :)
      real(dp), allocatable :: tens(:, :, :), tcart3(:, :, :), half(:, :, :)
      real(dp), allocatable :: xr(:), vis(:, :), ysend(:)
      real(dp), allocatable :: ki(:, :), kt(:, :), part(:), eo(:), ev(:)
      integer(int64) :: novr8, nsend8, nq8
      integer :: npr
      real(dp) :: e_os, e_ss
      real(dp) :: es(2)
      type(timer_type) :: tc
      type(step_table_t) :: stab
      logical :: dense_mask(NSTEP)
      integer(int64) :: nv8, nao8, naux8

      nao = b%nao
      ! `cmo`/`eps` above are declared (b%nao, b%nao)/(b%nao): ALWAYS the
      ! Cartesian AO space, whether `b` itself is spherical or not. A
      ! caller whose SCF ran spherical (`b%spherical`) does not hand this
      ! its raw (b%nao_sph, b%nao_sph) orbitals -- that would read past
      ! the end of a smaller array under sequence association -- but
      ! `rimp2_driver_mo_to_cart`'s own fold-and-pad result instead: every real
      ! MO transformed into the Cartesian AO space (C_cart = T C_sph, MO
      ! space itself does not care which AO representation it came from),
      ! padded out to `b%nao` columns with zero-coefficient, PAD_EPS
      ! "virtuals" that contribute exactly zero to every pair energy
      ! below. This routine never reads `b%spherical` itself.
      !
      ! The auxiliary basis `aux` is a separate question: it follows its
      ! own `aux%spherical` (trc_capi.F90's `basis_json` sets it to match
      ! the orbital basis, pyscf's own `auxmol.cart = mol.cart` rule, but
      ! this routine only ever reads `aux%spherical`). `naux_c`, the
      ! CARTESIAN count, is every kernel call and block-walk below
      ! (always Cartesian); `naux`, folded when `aux%spherical`, is every
      ! other dimension -- the metric, its Cholesky factor, `xr`'s own
      ! width, and (via `dfjk%dims`) what a chained call's own aux count
      ! is checked against.
      naux_c = aux%nao
      naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      nfz = 0
      if (present(nfrozen)) nfz = max(nfrozen, 0)
      no = nocc - nfz
      ! For a spherical caller, `nao` (Cartesian) is bigger than the
      ! true spherical MO count `b%nao_sph` -- rimp2_driver_mo_to_cart's own
      ! header explains why this routine is never told that -- so `nv`
      ! below runs over `b%nao - b%nao_sph` zero-coefficient pad
      ! "virtuals" that cost real GEMM flops and contribute nothing.
      ! Measured from shell counts alone (no kernel change, not
      ! profiled): water/cc-pVDZ(sph) is nv=20 Cartesian vs. 19 true
      ! (1 pad column, one O d-shell's 6-vs-5), formamide/cc-pVDZ(sph)
      ! is nv=48 vs. 45 (3 pad columns, one d-shell per heavy atom) --
      ! a 5-7% larger nv, which an O(nv^2) GEMM (the dominant scaling
      ! of the half-transform/pair-energy steps below) turns into an
      ! 11-14% larger flop count; any O(nv^3) term costs proportionally
      ! more still. Don't change the kernels for this -- it is the
      ! price of keeping `trc_rimp2_run` Cartesian-only and reusing it
      ! unchanged for a spherical caller (rimp2_driver_mo_to_cart's own
      ! header).
      nv = nao - nocc
      nov = no*nv
      res%naux = naux
      res%nocc = no
      res%nvir = nv
      if (no < 1 .or. nv < 1) then
         res%message = "trc_rimp2: nothing to correlate"
         return
      end if
      if (maxval(b%sh_l) > DF_LMAX) then
         res%message = "trc_rimp2: the orbital basis carries angular momentum above what the "// &
                       "density-fitting kernels were generated for (see gen_df.py --lmax)"
         return
      end if
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         res%message = "trc_rimp2: the auxiliary basis carries angular momentum above what the "// &
                       "density-fitting kernels were generated for (see gen_df.py --lmax-aux)"
         return
      end if
      chained = .false.
      bchain => null()
      if (present(dfjk)) then
         if (dfjk%is_ready()) then
            call dfjk%dims(n_chk, naux_chk, npack_chk)
            if (n_chk /= nao .or. naux_chk /= naux) then
               res%message = "trc_rimp2: the pre-fitted B does not match this orbital/auxiliary basis"
               return
            end if
            bchain => dfjk%packed_b()
            chained = .true.
         end if
      end if
      nblk = naux_c
      if (present(aux_block)) nblk = max(aux_block, 1)
      rank = 0
      nranks = 1
      if (present(comm)) then
         rank = comm%rank()
         nranks = comm%size()
      end if
      if (present(la)) then
         lp => la
      else
         lp => la_own
      end if
      ! This rank's occupied block, O_rank = lo .. lo + nb - 1, and every
      ! rank's, which the ring needs to size what it receives.
      allocate (blo(0:nranks - 1), bn(0:nranks - 1))
      do i = 0, nranks - 1
         bn(i) = no/nranks + merge(1, 0, i < mod(no, nranks))
         blo(i) = i*(no/nranks) + min(i, mod(no, nranks)) + 1
      end do
      lo = blo(rank)
      nb = bn(rank)
      nmax = maxval(bn)
      novr = nb*nv
      novr8 = novr
      ! Every rank's auxiliary shells, P_q = shells shlo(q) .. shhi(q) (empty
      ! when shhi < shlo): a shell goes to rank q while fewer than
      ! naux_c (q + 1) / nranks Cartesian functions come before it, so the
      ! integral work, which is Cartesian, is split evenly. plo(q) is where
      ! P_q starts among X's columns (spherical when `aux` is) and
      ! plo(q + 1) - plo(q) is its width.
      allocate (shlo(0:nranks - 1), shhi(0:nranks - 1), plo(0:nranks))
      ish = 1
      do q = 0, nranks - 1
         shlo(q) = ish
         do while (ish <= aux%nshell)
            if (int(aux%sh_ao(ish) - 1, int64) >= (int(naux_c, int64)*(q + 1))/nranks) exit
            ish = ish + 1
         end do
         shhi(q) = ish - 1
         if (shlo(q) <= aux%nshell) then
            plo(q) = merge(aux%sh_ao_sph(shlo(q)), aux%sh_ao(shlo(q)), aux%spherical)
         else
            plo(q) = naux + 1
         end if
      end do
      plo(nranks) = naux + 1
      ! Create-once: a caller's handles are kept, the workspace only grows.
      call lp%init(naux)

      talk = .false.
      if (present(verbose)) talk = verbose .and. rank == 0
      nv8 = nv
      nao8 = nao
      naux8 = naux
      ! The metric and its factor run in full on every rank; the integrals are
      ! split by auxiliary shell.
      call stab%init("rimp2", STEP_NAME, replicated_steps=[ST_METRIC, ST_CHOL])
      call lp%sync()
      call tw0%start()
      ! --- 1. the metric and its Cholesky factor (skipped when chained) ----
      if (.not. chained) then
         allocate (jm(naux, naux))
         !$acc enter data create(jm)
         call stab%step_on(lp, ST_METRIC, tc)
         if (aux%spherical) then
            allocate (jmc(naux_c, naux_c))
            !$acc enter data create(jmc)
            call trc_df_2c(aux, jmc)
            call sph_metric_cart_to_sph(aux, lp, jmc, jm)
            !$acc exit data delete(jmc)
            deallocate (jmc)
         else
            call trc_df_2c(aux, jm)
         end if
         call stab%step_off(lp, ST_METRIC, tc)
         call stab%step_on(lp, ST_CHOL, tc)
         call lp%potrf(naux, jm, naux)
         call stab%step_off(lp, ST_CHOL, tc)
         call stab%add_flop(ST_CHOL, naux8**3/3)
         if (talk) then
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,f8.2,a)") "  rimp2: metric (", naux, ") and Cholesky ", since(lp, tw0), " s"
               call global_logger%performance(trim(logline))
            end block
         end if
      else if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a)") "  rimp2: metric, Cholesky and (mu nu|P) skipped -- chained off --df's fit"
            call global_logger%performance(trim(logline))
         end block
      end if

      ! --- 3. (ia|P) for every correlated i and this rank's P, an auxiliary
      ! block at a time. The rows for O_q, q /= rank, go to `ysend`, one
      ! contiguous (bn(q) nv) x |P_rank| segment per q starting at soff(q);
      ! the rows for O_rank go straight into `xr`, in the columns they keep.
      direct = trc_mpi_device_direct()
      npr = plo(rank + 1) - plo(rank)
      allocate (soff(0:nranks - 1))
      nsend8 = 0
      do q = 0, nranks - 1
         soff(q) = nsend8
         if (q /= rank) nsend8 = nsend8 + int(bn(q), int64)*nv8*int(npr, int64)
      end do
      allocate (cocc(nao, no), cvir(nao, nv), xr(max(novr8*naux8, 1_int64)), ysend(max(nsend8, 1_int64)))
      cocc = cmo(:, nfz + 1:nocc)
      cvir = cmo(:, nocc + 1:nao)
      !$acc enter data copyin(cocc, cvir) create(xr, ysend)
      s0 = shlo(rank)
      do while (s0 <= shhi(rank))
         ! Shells s0..s1 hold at most nblk functions (and at least one shell).
         s1 = s0
         do ish = s0 + 1, shhi(rank)
            if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > nblk) exit
            s1 = ish
         end do
         p0c = aux%sh_ao(s0)
         npc = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0c
         call aux%subset(s0, s1, aux_subset)
         call aux_subset%to_device()
         ! Effective (spherical, when `aux%spherical`) offset/width this
         ! block's columns land at in `xr` below -- `aux_subset%nao_sph` is
         ! exactly this block's own spherical width, since `subset` never
         ! cuts a shell in two.
         if (aux%spherical) then
            p0 = aux%sh_ao_sph(s0)
            np = aux_subset%nao_sph
         else
            p0 = p0c
            np = npc
         end if
         allocate (tens(nao, nao, np), half(no, nao, np))
         !$acc enter data create(tens, half)
         call stab%step_on(lp, ST_3C, tc)
         if (chained) then
            ! B, already L^-T applied (and, when `aux` is spherical, the
            ! P-fold too -- `dfjk%setup` does both before packing),
            ! packed mu >= nu: unpack this block's columns straight into
            ! the dense tensor the half-transform below reads either way
            ! -- no (mu nu|P) build, and no P-fold, at all.
            do concurrent(q=1:np, mu=1:nao, nu=1:nao)
               tens(mu, nu, q) = bchain(rimp2_pack_index(max(mu, nu), min(mu, nu)), p0 - 1 + q)
            end do
         else if (aux%spherical) then
            ! (mu nu|Pc), built Cartesian, folded into `tens` on the P
            ! index alone -- same `sph_aux3c_to_sph` `df_jk_setup` folds
            ! its own block through, before the half-transform below
            ! reads it (which never needs to know the difference).
            allocate (tcart3(nao, nao, npc))
            !$acc enter data create(tcart3)
            call trc_df_3c(b, pl, aux_subset, tcart3)
            call sph_aux3c_to_sph(aux_subset, nao, tcart3, tens, lp)
            !$acc exit data delete(tcart3)
            deallocate (tcart3)
         else
            call trc_df_3c(b, pl, aux_subset, tens)
         end if
         call stab%step_off(lp, ST_3C, tc)
         call stab%step_on(lp, ST_TRANS, tc)
         ! A few library calls per block, not two per auxiliary function:
         ! the occupied side of the whole block is one GEMM,
         !   half(i, nu q) = sum_mu C_occ(mu, i) T(mu, nu q)   (no x nao np),
         ! and the virtual side one strided batch over q with C_vir shared
         ! per destination rank, on that rank's rows of `half`,
         !   X(:, P) = C_vir^T half_q^T                        (nv x bn(q), a fastest).
         ! 2 x 3840 calls at w40 cc-pVDZ ran at 41% of the device's peak.
         call lp%gemm("T", "N", no, nao*np, nao, 1.0_dp, cocc, nao, tens, nao, 0.0_dp, half, no)
         do q = 0, nranks - 1
            if (bn(q) == 0) cycle
            if (q == rank) then
               call lp%gemm_strided("T", "T", nv, nb, nao, 1.0_dp, cvir, nao, 0_int64, half(lo, 1, 1), no, &
                                    int(no, int64)*nao8, 0.0_dp, xr(1 + int(p0 - 1, int64)*novr8), nv, novr8, np)
            else
               nq8 = int(bn(q), int64)*nv8
               call lp%gemm_strided("T", "T", nv, bn(q), nao, 1.0_dp, cvir, nao, 0_int64, half(blo(q), 1, 1), no, &
                                    int(no, int64)*nao8, 0.0_dp, ysend(1 + soff(q) + int(p0 - plo(rank), int64)*nq8), &
                                    nv, nq8, np)
            end if
         end do
         call stab%step_off(lp, ST_TRANS, tc)
         call stab%add_flop(ST_TRANS, int(np, int64)*(2*int(no, int64)*nao8*nao8 + 2*nv8*int(no, int64)*nao8))
         !$acc exit data delete(tens, half)
         deallocate (tens, half)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      if (nranks > 1) then
         call stab%step_on(lp, ST_XCHG, tc)
         call exchange_x(comm, rank, nranks, nv8, bn, plo, soff, ysend, xr, direct)
         call stab%step_off(lp, ST_XCHG, tc)
      end if
      !$acc exit data delete(ysend)
      deallocate (ysend)
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,f8.2,a,f8.2,a,f8.2,a)") "  rimp2: (mu nu|P) ", stab%t(ST_3C), " s, transform to (ia|P) ", &
               stab%t(ST_TRANS), " s, exchange ", stab%t(ST_XCHG), " s"
            call global_logger%performance(trim(logline))
         end block
      end if

      ! --- 4. B = X L^-T: solve B L^T = X from the right, this rank's rows ---
      ! (skipped when chained: `xr` is already built from B, which has
      ! L^-T applied already -- see this routine's own `dfjk` comment)
      if (.not. chained) then
         call stab%step_on(lp, ST_TRSM, tc)
         if (nb > 0) call lp%trsm("R", "L", "T", novr, naux, 1.0_dp, jm, naux, xr, novr)
         call stab%step_off(lp, ST_TRSM, tc)
         call stab%add_flop(ST_TRSM, novr8*naux8*naux8)
         if (talk) then
            block
               character(len=1024) :: logline
               write (logline, "(a,f8.2,a)") "  rimp2: B = X L^-T ", stab%t(ST_TRSM), " s"
               call global_logger%performance(trim(logline))
            end block
         end if
      end if
      call lp%sync()
      call tw1%start()

      ! --- 5. the pair energies, over the ring ------------------------------
      ! The visiting blocks need two buffers, the one in use and the one
      ! arriving; one rank needs neither.
      nstep_ring = nranks/2
      if (talk .and. nranks > 1) then
         block
            character(len=1024) :: logline
            write (logline, "(a,a)") "  rimp2: B exchange ", &
               trim(merge("device to device (CUDA-aware MPI)", "through the host                 ", direct))
            call global_logger%performance(trim(logline))
         end block
      end if
      if (nranks > 1) then
         allocate (vis(int(nmax, int64)*nv8*naux8, 2))
         !$acc enter data create(vis)
      end if
      ncolmax = nmax*nv
      ng = k_group(nmax, nv, ncolmax)
      ldk = ng*nv
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a,i0,a)") "  rimp2: K GEMMs over ", ng, " occupied orbitals at a time, ", &
               nranks, " rank(s)"
            call global_logger%performance(trim(logline))
         end block
      end if
      allocate (ki(ldk, ncolmax), kt(ncolmax, nv), part(2*ncolmax), eo(no), ev(nv))
      eo = eps(nfz + 1:nocc)
      ev = eps(nocc + 1:nao)
      !$acc enter data create(ki, kt, part) copyin(eo, ev)
      e_os = 0.0_dp
      e_ss = 0.0_dp
      cur = 0
      do step = 0, nstep_ring
         qv = mod(rank + step, nranks)
         if (step > 0) then
            call stab%step_on(lp, ST_RING, tc)
            call ring_shift(comm, rank, nranks, xr, vis, cur, &
                            int(bn(mod(rank + step - 1, nranks)), int64)*nv8*naux8, &
                            int(bn(qv), int64)*nv8*naux8, step == 1, direct)
            call stab%step_off(lp, ST_RING, tc)
         end if
         ! Rows i in ilo..ihi of this rank's block against columns j in
         ! jlo..jhi of block qv; `diag` keeps j >= i.
         diag = step == 0
         ilo = lo
         ihi = lo + nb - 1
         jbase = blo(qv)
         jlo = jbase
         jhi = jbase + bn(qv) - 1
         if (2*step == nranks) then
            if (rank < qv) then
               ihi = lo + nb/2 - 1
            else
               jlo = jbase + bn(qv)/2
            end if
         end if
         if (step == 0) then
            ldv = novr
         else
            ldv = bn(qv)*nv
         end if
         ! Groups of `ng` occupied orbitals: one GEMM per group, K(i a, jb)
         ! for i in the group and j from the group's first column on, so
         ! cuBLAS gets (ng nv) rows instead of nv. On the diagonal block it
         ! also computes the j < i corner inside the group, which the pair
         ! kernel skips.
         do i0 = ilo, ihi, ng
            i1 = min(ihi, i0 + ng - 1)
            jstart = jlo
            if (diag) jstart = i0
            ncol = (jhi - jstart + 1)*nv
            if (ncol <= 0) cycle
            call stab%step_on(lp, ST_KGEMM, tc)
            if (step == 0) then
               call lp%gemm("N", "T", (i1 - i0 + 1)*nv, ncol, naux, 1.0_dp, xr(1 + (i0 - lo)*nv), novr, &
                            xr(1 + (jstart - jbase)*nv), ldv, 0.0_dp, ki, ldk)
            else
               call lp%gemm("N", "T", (i1 - i0 + 1)*nv, ncol, naux, 1.0_dp, xr(1 + (i0 - lo)*nv), novr, &
                            vis(1 + (jstart - jbase)*nv, cur), ldv, 0.0_dp, ki, ldk)
            end if
            call stab%step_off(lp, ST_KGEMM, tc)
            call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nv8*int(ncol, int64)*naux8)
            do i = i0, i1
               ! K_i(a, jb) is the block of rows roff + a and, on the
               ! diagonal, columns from j = i; off it, every column.
               roff = (i - i0)*nv
               if (diag) then
                  coff = (i - i0)*nv
                  jfirst = i
               else
                  coff = 0
                  jfirst = jstart
               end if
               ncol_i = (jhi - jfirst + 1)*nv
               call stab%step_on(lp, ST_PAIR, tc)
               ! K_i^T, so the kernel reads (ia|jb) and (ib|ja) both coalesced
               call lp%transpose(nv, ncol_i, ki(roff + 1, coff + 1), ldk, kt, ncolmax)
               call pair_energy(i, jfirst, no, nv, ldk, roff, coff, ncolmax, ncol_i, eo, ev, ki, kt, part)
               !$acc update self(part(1:2*ncol_i))
               call stab%step_off(lp, ST_PAIR, tc)
               call stab%add_flop(ST_PAIR, 10*nv8*int(ncol_i, int64))
               e_os = e_os + sum(part(1:2*ncol_i:2))
               e_ss = e_ss + sum(part(2:2*ncol_i:2))
            end do
         end do
      end do
      !$acc exit data delete(ki, kt, part, eo, ev, xr, cocc, cvir)
      if (.not. chained) then
         !$acc exit data delete(jm)
      end if
      if (nranks > 1) then
         !$acc exit data delete(vis)
      end if
      if (nranks > 1) then
         es = [e_os, e_ss]
         call allreduce(comm, es, op=MPI_SUM)
         call bcast(comm, es, 2, 0)
         e_os = es(1)
         e_ss = es(2)
      end if
      ! The job's table, not rank 0's share: the slowest rank's time and the
      ! work of every rank. The energy reduction above has already lined the
      ! ranks up, so rank 0's wall time below is the job's.
      if (nranks > 1) call stab%reduce(comm)
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,f8.2,a,f8.2,a)") "  rimp2: pair energies ", since(lp, tw1), " s;  total ", since(lp, tw0), " s"
            call global_logger%performance(trim(logline))
         end block
         dense_mask = .false.
         dense_mask(ST_CHOL) = .true.
         dense_mask(ST_TRANS) = .true.
         dense_mask(ST_TRSM) = .true.
         dense_mask(ST_KGEMM) = .true.
         call stab%report(since(lp, tw0), rank, nranks, subtotal_mask=dense_mask, subtotal_label="all dense algebra")
      end if
      res%e_os = e_os
      res%e_ss = e_ss
      res%e_corr = e_os + e_ss
      if (.not. present(la)) call la_own%release()
   end subroutine trc_rimp2_run

   !
   ! `trc_rimp2_run`/`trc_urimp2_run`'s own fold for a spherical orbital
   ! basis: `csph`/`epssph` are an SCF's own (b%nao_sph, b%nao_sph)/
   ! (b%nao_sph) orbitals and energies (`trc_scf_result_t%cmo`/`eps`,
   ! already at the SCF's own `scf_space_nao` size); `ccart`/`epscart`
   ! come back (b%nao, b%nao)/(b%nao), the Cartesian AO space RI-MP2's own
   ! three-centre tensor is built in.
   !
   ! Every real MO is `sph_mo_to_cart`'s T C_sph, trc_sph.F90's one-sided
   ! fold (an MO's AO row crosses the spherical/Cartesian boundary, its
   ! column -- the MO label -- does not: MO-space quantities do not care
   ! which AO representation produced them). The Cartesian AO count is
   ! always at least the spherical MO count (a Cartesian shell never has
   ! fewer functions than its spherical one), so `b%nao - b%nao_sph`
   ! columns are left over; they are filled with a ZERO coefficient
   ! vector and `PAD_EPS`, a finite energy well outside any real
   ! spectrum. `trc_rimp2_run`'s own `nv = nao - nocc` (nao = b%nao here)
   ! then runs over these padded columns exactly as it would a real
   ! virtual, but every (ia|jb) with a or b among them is a GEMM
   ! contraction against a zero `C_vir` column and comes out exactly
   ! zero -- PAD_EPS only has to keep `1/(e_i + e_j - e_a - e_b)` FINITE
   ! (0 * finite = 0), not physically meaningful, since it is always
   ! multiplied by that zero. This is what lets `trc_rimp2_run` itself
   ! stay unchanged, Cartesian-interface, for a spherical caller.
   !
   subroutine rimp2_driver_mo_to_cart(b, csph, epssph, ccart, epscart)
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: csph(b%nao_sph, b%nao_sph), epssph(b%nao_sph)
      real(dp), intent(out) :: ccart(b%nao, b%nao), epscart(b%nao)
      integer :: ns, nc

      ns = b%nao_sph
      nc = b%nao
      !$acc enter data copyin(csph) create(ccart)
      call sph_mo_to_cart(b, ns, csph, ccart(:, 1:ns))
      !$acc update self(ccart(:, 1:ns))
      !$acc exit data delete(csph, ccart)
      ! These `nc - ns` zero columns are exactly the pad virtuals
      ! trc_rimp2_run's own `nv = nao - nocc` then costs real GEMM
      ! flops on -- see that routine's own comment for the measured
      ! fraction.
      if (nc > ns) ccart(:, ns + 1:nc) = 0.0_dp
      epscart(1:ns) = epssph
      if (nc > ns) epscart(ns + 1:nc) = PAD_EPS
   end subroutine rimp2_driver_mo_to_cart

   !
   ! RI-MP2 for an UNRESTRICTED (UHF or UKS) reference, on the device, one
   ! rank only.
   !
   ! The metric, its Cholesky factor and (mu nu|P) do not know about spin --
   ! they are built exactly once, as trc_rimp2_run's own steps 1-2 -- and
   ! transformed TWICE, into B_a (from C_occ_a, C_vir_a) and B_b (from
   ! C_occ_b, C_vir_b), by the same GEMM-plus-strided-batch half-transform
   ! and the same triangular solve as the restricted routine's steps 3-4,
   ! called once per spin.
   !
   ! Three pair blocks follow, all from B_a and B_b:
   !
   !   E_aa = 1/2 sum_{ij in a, ab in a} (ia|jb)[(ia|jb) - (ib|ja)] / D
   !   E_bb = 1/2 sum_{ij in b, ab in b} (ia|jb)[(ia|jb) - (ib|ja)] / D
   !   E_ab =     sum_{i in a, J in b, a in a, B in b} (ia|JB)^2 / D
   !
   ! E_aa and E_bb reuse trc_rimp2_run's own K GEMM and `pair_energy` kernel
   ! verbatim -- same occupied-pair triangle (j >= i), same weight (2 off
   ! the diagonal, 1 on it) -- and only the same-spin half of what that
   ! kernel returns (`pss`) is kept; the opposite-spin half it also computes
   ! (`pos`) has no meaning for two electrons of the SAME spin and is
   ! discarded. The overall 1/2 is the one difference from the restricted
   ! energy: expanding the standard spin-orbital result
   !
   !     E_aa = 1/4 sum_{ijab} [(ia|jb) - (ib|ja)]^2 / D     (i, j, a, b all
   !                                                           independent,
   !                                                           no triangle)
   !
   ! and using (ia|jb)^2 = (ib|ja)^2 over that full range gives exactly
   ! 1/2 sum_{ijab} (ia|jb)[(ia|jb)-(ib|ja)]/D -- restricted's own full
   ! double sum, weighted down to the triangle the same way restricted's
   ! E_ss is -- hence the 1/2 that sits outside the triangle sum below. See
   ! metalquicha's run_czt_uri_mp2 (mqc_czt_mp2.f90) for the same algebra
   ! checked against a conventional (non-fitted) transform.
   !
   ! E_ab has no such symmetry -- i (alpha) and J (beta) are already
   ! distinguishable electrons, so there is no exchange term and no pair to
   ! fold in half: the full rectangle of i against J, one GEMM
   ! K(iA, JB) = sum_Q B_a(iA, Q) B_b(JB, Q) per group of alpha occupied
   ! orbitals, and a kernel that only needs K itself (no transpose, no
   ! second read) since there is nothing to subtract.
   !
   ! MULTI-RANK: refused. The restricted ring distributes one occupied
   ! block of B and lets it meet itself (aa) through a ring of nranks/2
   ! steps; an open-shell run would need that ring run independently for aa
   ! (over alpha's block) and bb (over beta's, a different size whenever
   ! nalpha /= nbeta), plus a THIRD distribution for ab that is not a ring
   ! at all -- alpha's block never needs to meet another alpha block, only
   ! every beta block, so the "pass one step left" structure the ring
   ! exists for buys nothing there. Three different distributions behind
   ! one call is not a clean generalisation of the one the restricted path
   ! has, so this round runs one rank only and refuses anything else below,
   ! with a message that says so.
   !
   subroutine trc_urimp2_run(b, aux, pl, nalpha, nbeta, cmo_a, cmo_b, eps_a, eps_b, res, &
                             nfrozen, aux_block, la, comm, verbose)
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nalpha, nbeta
      real(dp), intent(in) :: cmo_a(b%nao, b%nao), cmo_b(b%nao, b%nao), eps_a(b%nao), eps_b(b%nao)
      type(trc_rimp2_result_t), intent(out) :: res
      integer, intent(in), optional :: nfrozen, aux_block
      type(trc_linalg_t), intent(inout), optional, target :: la
      type(comm_t), intent(in), optional :: comm
      logical, intent(in), optional :: verbose

      type(timer_type) :: tw0, tw1
      logical :: talk
      type(trc_linalg_t), target :: la_own
      type(trc_linalg_t), pointer :: lp
      type(trc_basis_t) :: aux_subset
      integer :: nao, naux, naux_c, nfz, noa, nob, nva, nvb, nblk, s0, s1, p0, np, p0c, npc, ish
      integer :: rank, nranks
      integer :: novr_a, novr_b, ng_aa, ng_bb, ng_ab, ldk_aa, ldk_bb, ldk_ab
      integer :: i0, i1, i, jstart, jfirst, ncol, ncol_i, roff, coff, ncolmax_aa, ncolmax_bb
      real(dp), allocatable :: jm(:, :), jmc(:, :)
      real(dp), allocatable :: cocc_a(:, :), cvir_a(:, :), cocc_b(:, :), cvir_b(:, :)
      real(dp), allocatable :: tens(:, :, :), tcart3(:, :, :), half_a(:, :, :), half_b(:, :, :)
      real(dp), allocatable :: xr_a(:), xr_b(:)
      real(dp), allocatable :: ki_aa(:, :), kt_aa(:, :), part_aa(:)
      real(dp), allocatable :: ki_bb(:, :), kt_bb(:, :), part_bb(:)
      real(dp), allocatable :: ki_ab(:, :), part_ab(:)
      real(dp), allocatable :: eo_a(:), ev_a(:), eo_b(:), ev_b(:)
      real(dp) :: e_aa, e_bb, e_ab
      type(timer_type) :: tc
      type(step_table_t) :: stab
      logical :: dense_mask(NSTEP)
      integer(int64) :: nao8, naux8, nva8, nvb8, noa8, nob8, novr_a8, novr_b8

      nao = b%nao
      ! `cmo_a`/`cmo_b`/`eps_a`/`eps_b` are the Cartesian-space, fold-and-
      ! pad result of `rimp2_driver_mo_to_cart` when `b` is spherical -- see
      ! `trc_rimp2_run`'s own header on that routine, which this one reads
      ! exactly as it does. `naux_c`/`naux` split the auxiliary basis the
      ! same way too (CARTESIAN count for the kernel/block-walk below,
      ! folded when `aux%spherical` for every other dimension).
      naux_c = aux%nao
      naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      rank = 0
      nranks = 1
      if (present(comm)) then
         rank = comm%rank()
         nranks = comm%size()
      end if
      if (nranks > 1) then
         res%message = "trc_urimp2: multi-rank UHF RI-MP2 is not implemented this round -- "// &
                       "the restricted ring distributes one block of B and lets it meet "// &
                       "itself; an open-shell run needs that for aa AND bb (different sizes "// &
                       "whenever nalpha /= nbeta) plus a third, non-ring distribution for the "// &
                       "mixed ab block. Run on one rank (no communicator, or a single-rank one)"
         return
      end if
      nfz = 0
      if (present(nfrozen)) nfz = max(nfrozen, 0)
      noa = nalpha - nfz
      nob = nbeta - nfz
      nva = nao - nalpha
      nvb = nao - nbeta
      res%naux = naux
      res%nocc = noa + nob    ! spin orbitals, not spatial -- see this
      res%nvir = nva + nvb    ! routine's own header and fill_ump2_result's
      ! equivalent in metalquicha
      if (noa < 1 .or. nob < 1 .or. nva < 1 .or. nvb < 1) then
         res%message = "trc_urimp2: nothing to correlate in one spin channel"
         return
      end if
      if (maxval(b%sh_l) > DF_LMAX) then
         res%message = "trc_urimp2: the orbital basis carries angular momentum above what the "// &
                       "density-fitting kernels were generated for (see gen_df.py --lmax)"
         return
      end if
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         res%message = "trc_urimp2: the auxiliary basis carries angular momentum above what the "// &
                       "density-fitting kernels were generated for (see gen_df.py --lmax-aux)"
         return
      end if
      nblk = naux_c
      if (present(aux_block)) nblk = max(aux_block, 1)
      if (present(la)) then
         lp => la
      else
         lp => la_own
      end if
      novr_a = noa*nva
      novr_b = nob*nvb
      novr_a8 = novr_a
      novr_b8 = novr_b
      nao8 = nao
      naux8 = naux
      nva8 = nva
      nvb8 = nvb
      noa8 = noa
      nob8 = nob
      call lp%init(naux)

      talk = .false.
      if (present(verbose)) talk = verbose .and. rank == 0
      call stab%init("rimp2u", STEP_NAME)
      call lp%sync()
      call tw0%start()

      ! --- 1. the metric and its Cholesky factor (spin-independent) --------
      allocate (jm(naux, naux))
      !$acc enter data create(jm)
      call stab%step_on(lp, ST_METRIC, tc)
      if (aux%spherical) then
         allocate (jmc(naux_c, naux_c))
         !$acc enter data create(jmc)
         call trc_df_2c(aux, jmc)
         call sph_metric_cart_to_sph(aux, lp, jmc, jm)
         !$acc exit data delete(jmc)
         deallocate (jmc)
      else
         call trc_df_2c(aux, jm)
      end if
      call stab%step_off(lp, ST_METRIC, tc)
      call stab%step_on(lp, ST_CHOL, tc)
      call lp%potrf(naux, jm, naux)
      call stab%step_off(lp, ST_CHOL, tc)
      call stab%add_flop(ST_CHOL, naux8**3/3)
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a,f8.2,a)") "  rimp2u: metric (", naux, ") and Cholesky ", since(lp, tw0), " s"
            call global_logger%performance(trim(logline))
         end block
      end if

      ! --- 2-3. (mu nu|P) once per aux block, (ia|P) for BOTH spins from it -
      allocate (cocc_a(nao, noa), cvir_a(nao, nva), xr_a(max(novr_a8*naux8, 1_int64)))
      allocate (cocc_b(nao, nob), cvir_b(nao, nvb), xr_b(max(novr_b8*naux8, 1_int64)))
      cocc_a = cmo_a(:, nfz + 1:nalpha)
      cvir_a = cmo_a(:, nalpha + 1:nao)
      cocc_b = cmo_b(:, nfz + 1:nbeta)
      cvir_b = cmo_b(:, nbeta + 1:nao)
      !$acc enter data copyin(cocc_a, cvir_a, cocc_b, cvir_b) create(xr_a, xr_b)
      s0 = 1
      do while (s0 <= aux%nshell)
         s1 = s0
         do ish = s0 + 1, aux%nshell
            if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > nblk) exit
            s1 = ish
         end do
         p0c = aux%sh_ao(s0)
         npc = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0c
         call aux%subset(s0, s1, aux_subset)
         call aux_subset%to_device()
         if (aux%spherical) then
            p0 = aux%sh_ao_sph(s0)
            np = aux_subset%nao_sph
         else
            p0 = p0c
            np = npc
         end if
         allocate (tens(nao, nao, np))
         !$acc enter data create(tens)
         call stab%step_on(lp, ST_3C, tc)
         if (aux%spherical) then
            allocate (tcart3(nao, nao, npc))
            !$acc enter data create(tcart3)
            call trc_df_3c(b, pl, aux_subset, tcart3)
            call sph_aux3c_to_sph(aux_subset, nao, tcart3, tens, lp)
            !$acc exit data delete(tcart3)
            deallocate (tcart3)
         else
            call trc_df_3c(b, pl, aux_subset, tens)
         end if
         call stab%step_off(lp, ST_3C, tc)

         allocate (half_a(noa, nao, np))
         !$acc enter data create(half_a)
         call stab%step_on(lp, ST_TRANS, tc)
         call lp%gemm("T", "N", noa, nao*np, nao, 1.0_dp, cocc_a, nao, tens, nao, 0.0_dp, half_a, noa)
         call lp%gemm_strided("T", "T", nva, noa, nao, 1.0_dp, cvir_a, nao, 0_int64, half_a, noa, &
                              int(noa, int64)*nao8, 0.0_dp, xr_a(1 + int(p0 - 1, int64)*novr_a8), nva, &
                              novr_a8, np)
         call stab%step_off(lp, ST_TRANS, tc)
         call stab%add_flop(ST_TRANS, int(np, int64)*(2*noa8*nao8*nao8 + 2*nva8*noa8*nao8))
         !$acc exit data delete(half_a)
         deallocate (half_a)

         allocate (half_b(nob, nao, np))
         !$acc enter data create(half_b)
         call stab%step_on(lp, ST_TRANS, tc)
         call lp%gemm("T", "N", nob, nao*np, nao, 1.0_dp, cocc_b, nao, tens, nao, 0.0_dp, half_b, nob)
         call lp%gemm_strided("T", "T", nvb, nob, nao, 1.0_dp, cvir_b, nao, 0_int64, half_b, nob, &
                              int(nob, int64)*nao8, 0.0_dp, xr_b(1 + int(p0 - 1, int64)*novr_b8), nvb, &
                              novr_b8, np)
         call stab%step_off(lp, ST_TRANS, tc)
         call stab%add_flop(ST_TRANS, int(np, int64)*(2*nob8*nao8*nao8 + 2*nvb8*nob8*nao8))
         !$acc exit data delete(half_b)
         deallocate (half_b)

         !$acc exit data delete(tens)
         deallocate (tens)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,f8.2,a,f8.2,a)") "  rimp2u: (mu nu|P) ", stab%t(ST_3C), " s, transform to (ia|P) ", &
               stab%t(ST_TRANS), " s"
            call global_logger%performance(trim(logline))
         end block
      end if

      ! --- 4. B = X L^-T, one triangular solve per spin ---------------------
      call stab%step_on(lp, ST_TRSM, tc)
      call lp%trsm("R", "L", "T", novr_a, naux, 1.0_dp, jm, naux, xr_a, novr_a)
      call stab%step_off(lp, ST_TRSM, tc)
      call stab%add_flop(ST_TRSM, novr_a8*naux8*naux8)
      call stab%step_on(lp, ST_TRSM, tc)
      call lp%trsm("R", "L", "T", novr_b, naux, 1.0_dp, jm, naux, xr_b, novr_b)
      call stab%step_off(lp, ST_TRSM, tc)
      call stab%add_flop(ST_TRSM, novr_b8*naux8*naux8)
      !$acc exit data delete(jm)
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,f8.2,a)") "  rimp2u: B = X L^-T ", stab%t(ST_TRSM), " s"
            call global_logger%performance(trim(logline))
         end block
      end if
      call lp%sync()
      call tw1%start()

      ! --- 5a/5b. the like-spin pair energies, restricted's own triangle ---
      allocate (eo_a(noa), ev_a(nva), eo_b(nob), ev_b(nvb))
      eo_a = eps_a(nfz + 1:nalpha)
      ev_a = eps_a(nalpha + 1:nao)
      eo_b = eps_b(nfz + 1:nbeta)
      ev_b = eps_b(nbeta + 1:nao)
      !$acc enter data copyin(eo_a, ev_a, eo_b, ev_b)

      ncolmax_aa = noa*nva
      ng_aa = k_group(noa, nva, ncolmax_aa)
      ldk_aa = ng_aa*nva
      allocate (ki_aa(ldk_aa, ncolmax_aa), kt_aa(ncolmax_aa, nva), part_aa(2*ncolmax_aa))
      !$acc enter data create(ki_aa, kt_aa, part_aa)
      e_aa = 0.0_dp
      do i0 = 1, noa, ng_aa
         i1 = min(noa, i0 + ng_aa - 1)
         jstart = i0
         ncol = (noa - jstart + 1)*nva
         call stab%step_on(lp, ST_KGEMM, tc)
         call lp%gemm("N", "T", (i1 - i0 + 1)*nva, ncol, naux, 1.0_dp, xr_a(1 + (i0 - 1)*nva), novr_a, &
                      xr_a(1 + (jstart - 1)*nva), novr_a, 0.0_dp, ki_aa, ldk_aa)
         call stab%step_off(lp, ST_KGEMM, tc)
         call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nva8*int(ncol, int64)*naux8)
         do i = i0, i1
            roff = (i - i0)*nva
            coff = (i - i0)*nva
            jfirst = i
            ncol_i = (noa - jfirst + 1)*nva
            call stab%step_on(lp, ST_PAIR, tc)
            call lp%transpose(nva, ncol_i, ki_aa(roff + 1, coff + 1), ldk_aa, kt_aa, ncolmax_aa)
            call pair_energy(i, jfirst, noa, nva, ldk_aa, roff, coff, ncolmax_aa, ncol_i, eo_a, ev_a, &
                             ki_aa, kt_aa, part_aa)
            !$acc update self(part_aa(1:2*ncol_i))
            call stab%step_off(lp, ST_PAIR, tc)
            call stab%add_flop(ST_PAIR, 10*nva8*int(ncol_i, int64))
            e_aa = e_aa + sum(part_aa(2:2*ncol_i:2))   ! the same-spin half only; see this
         end do                                        ! routine's own header for the discarded "pos" half
      end do
      !$acc exit data delete(ki_aa, kt_aa, part_aa)
      deallocate (ki_aa, kt_aa, part_aa)
      e_aa = 0.5_dp*e_aa

      ncolmax_bb = nob*nvb
      ng_bb = k_group(nob, nvb, ncolmax_bb)
      ldk_bb = ng_bb*nvb
      allocate (ki_bb(ldk_bb, ncolmax_bb), kt_bb(ncolmax_bb, nvb), part_bb(2*ncolmax_bb))
      !$acc enter data create(ki_bb, kt_bb, part_bb)
      e_bb = 0.0_dp
      do i0 = 1, nob, ng_bb
         i1 = min(nob, i0 + ng_bb - 1)
         jstart = i0
         ncol = (nob - jstart + 1)*nvb
         call stab%step_on(lp, ST_KGEMM, tc)
         call lp%gemm("N", "T", (i1 - i0 + 1)*nvb, ncol, naux, 1.0_dp, xr_b(1 + (i0 - 1)*nvb), novr_b, &
                      xr_b(1 + (jstart - 1)*nvb), novr_b, 0.0_dp, ki_bb, ldk_bb)
         call stab%step_off(lp, ST_KGEMM, tc)
         call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nvb8*int(ncol, int64)*naux8)
         do i = i0, i1
            roff = (i - i0)*nvb
            coff = (i - i0)*nvb
            jfirst = i
            ncol_i = (nob - jfirst + 1)*nvb
            call stab%step_on(lp, ST_PAIR, tc)
            call lp%transpose(nvb, ncol_i, ki_bb(roff + 1, coff + 1), ldk_bb, kt_bb, ncolmax_bb)
            call pair_energy(i, jfirst, nob, nvb, ldk_bb, roff, coff, ncolmax_bb, ncol_i, eo_b, ev_b, &
                             ki_bb, kt_bb, part_bb)
            !$acc update self(part_bb(1:2*ncol_i))
            call stab%step_off(lp, ST_PAIR, tc)
            call stab%add_flop(ST_PAIR, 10*nvb8*int(ncol_i, int64))
            e_bb = e_bb + sum(part_bb(2:2*ncol_i:2))
         end do
      end do
      !$acc exit data delete(ki_bb, kt_bb, part_bb)
      deallocate (ki_bb, kt_bb, part_bb)
      e_bb = 0.5_dp*e_bb

      ! --- 5c. the mixed-spin block: the full rectangle, no exchange term --
      ng_ab = k_group(noa, nva, nob*nvb)
      ldk_ab = ng_ab*nva
      allocate (ki_ab(ldk_ab, nob*nvb), part_ab(nva*nob))
      !$acc enter data create(ki_ab, part_ab)
      e_ab = 0.0_dp
      do i0 = 1, noa, ng_ab
         i1 = min(noa, i0 + ng_ab - 1)
         ncol = nob*nvb
         call stab%step_on(lp, ST_KGEMM, tc)
         call lp%gemm("N", "T", (i1 - i0 + 1)*nva, ncol, naux, 1.0_dp, xr_a(1 + (i0 - 1)*nva), novr_a, &
                      xr_b, novr_b, 0.0_dp, ki_ab, ldk_ab)
         call stab%step_off(lp, ST_KGEMM, tc)
         call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nva8*int(ncol, int64)*naux8)
         do i = i0, i1
            roff = (i - i0)*nva
            call stab%step_on(lp, ST_PAIR, tc)
            call mixed_energy(i, nva, nob, nvb, ldk_ab, roff, eo_a, ev_a, eo_b, ev_b, ki_ab, part_ab)
            !$acc update self(part_ab(1:nva*nob))
            call stab%step_off(lp, ST_PAIR, tc)
            call stab%add_flop(ST_PAIR, 4_int64*nva8*nob8*nvb8)
            e_ab = e_ab + sum(part_ab(1:nva*nob))
         end do
      end do
      !$acc exit data delete(ki_ab, part_ab)
      deallocate (ki_ab, part_ab)

      !$acc exit data delete(eo_a, ev_a, eo_b, ev_b, xr_a, xr_b, cocc_a, cvir_a, cocc_b, cvir_b)
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,f8.2,a,f8.2,a)") "  rimp2u: pair energies ", since(lp, tw1), " s;  total ", since(lp, tw0), " s"
            call global_logger%performance(trim(logline))
         end block
         dense_mask = .false.
         dense_mask(ST_CHOL) = .true.
         dense_mask(ST_TRANS) = .true.
         dense_mask(ST_TRSM) = .true.
         dense_mask(ST_KGEMM) = .true.
         call stab%report(since(lp, tw0), rank, nranks, subtotal_mask=dense_mask, subtotal_label="all dense algebra")
      end if
      res%e_os = e_ab
      res%e_ss = e_aa + e_bb
      res%e_corr = res%e_os + res%e_ss
      if (.not. present(la)) call la_own%release()
   end subroutine trc_urimp2_run

   !
   ! The time `tm` has run, read after everything queued on the device has
   ! finished.
   !
   function since(lp, tm) result(t)
      type(trc_linalg_t), intent(in) :: lp
      type(timer_type), intent(in) :: tm
      real(dp) :: t
      call lp%sync()
      t = tm%get_elapsed_time()
   end function since

   !
   ! Occupied orbitals per K GEMM. One i gives cuBLAS nvir rows, which it
   ! tiles small; four give it a matrix tall enough for its larger tiles.
   ! Measured at w40 cc-pVDZ on a V100, K GEMM seconds for 1, 2, 4, 8, 16
   ! per group: 9.37, 9.09, 8.78, 9.24, 9.38 -- four is the best, and
   ! larger groups lose to the j < i corner each one computes and throws
   ! away. K is that many times larger, so on the device the group shrinks
   ! to what the free memory holds next to K^T and the energies, keeping
   ! 512 MB back for the libraries' own workspace.
   !
   !
   ! One step of the ring: pass the block this rank holds to the rank on its
   ! left and take the one on its right, through the host. At the first
   ! step the held block is this rank's own (`xr`), after that the visiting
   ! buffer `cur`, and what arrives goes into the other buffer, which
   ! becomes `cur`. Blocking calls in an order that cannot deadlock: even
   ! ranks send first and odd ranks receive first, so every send meets a
   ! receive that is already posted or about to be. In pieces of at most
   ! MSG_CHUNK doubles, since MPI counts are 32-bit.
   !
   subroutine ring_shift(comm, rank, nranks, xr, vis, cur, nsend, nrecv, first, direct)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: rank, nranks
      real(dp), intent(inout) :: xr(:), vis(:, :)
      integer, intent(inout) :: cur
      integer(int64), intent(in) :: nsend, nrecv
      logical, intent(in) :: first, direct
      integer :: left, right, nxt, nreq
      type(request_t), allocatable :: reqs(:)
      left = mod(rank - 1 + nranks, nranks)
      right = mod(rank + 1, nranks)
      nxt = merge(1, 3 - cur, first)
      nreq = int((nsend + MSG_CHUNK - 1)/MSG_CHUNK)
      allocate (reqs(max(nreq, 1)))
      if (direct) then
         ! The device addresses straight to MPI: inside host_data, through
         ! pic-mpi's explicit-count procedures, which pass the buffer on
         ! without a Fortran temporary (a temporary would be a host copy of
         ! a device address).
         !$acc host_data use_device(xr, vis)
         call send_block_dev(comm, xr, vis, cur, nsend, left, first, reqs)
         call recv_block_dev(comm, vis, nxt, nrecv, right)
         !$acc end host_data
         if (nreq > 0) call waitall(reqs(1:nreq))
         cur = nxt
         return
      end if
      if (first) then
         !$acc update self(xr(1:nsend))
      else
         !$acc update self(vis(1:nsend, cur))
      end if
      ! The send is posted first and non-blocking, so both directions of
      ! every link move at once; the receive then blocks, and the send is
      ! waited for before its buffer can be reused.
      call send_block(comm, xr, vis, cur, nsend, left, first, reqs)
      call recv_block(comm, vis, nxt, nrecv, right)
      if (nreq > 0) call waitall(reqs(1:nreq))
      !$acc update device(vis(1:nrecv, nxt))
      cur = nxt
   end subroutine ring_shift

   !
   ! The held block to `dest`, from the host copy, in MSG_CHUNK pieces.
   ! Module procedures rather than internal ones: see app/main.F90's
   ! trc_run for what a sibling internal procedure did under NVHPC's MPI.
   !
   subroutine send_block(comm, xr, vis, cur, n, dest, first, reqs)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in), target :: xr(:), vis(:, :)
      integer, intent(in) :: cur, dest
      integer(int64), intent(in) :: n
      logical, intent(in) :: first
      type(request_t), intent(inout) :: reqs(:)
      integer(int64) :: k, m
      integer :: r
      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         if (first) then
            call isend(comm, xr(k:k + m - 1), dest, 7, reqs(r))
         else
            call isend(comm, vis(k:k + m - 1, cur), dest, 7, reqs(r))
         end if
         k = k + m
      end do
   end subroutine send_block

   !
   ! The next block from `source` into buffer `c` of `vis`, in MSG_CHUNK pieces.
   !
   subroutine recv_block(comm, vis, c, n, source)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: vis(:, :)
      integer, intent(in) :: c, source
      integer(int64), intent(in) :: n
      integer(int64) :: k, m
      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call recv(comm, vis(k:k + m - 1, c), int(m), source, 7)
         k = k + m
      end do
   end subroutine recv_block

   !
   ! send_block and recv_block for device addresses, called inside
   ! `host_data use_device`: pic-mpi's explicit-count procedures by their
   ! specific names, since its generic interfaces misbehave there under
   ! nvfortran.
   !
   subroutine send_block_dev(comm, xr, vis, cur, n, dest, first, reqs)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in), target :: xr(:), vis(:, :)
      integer, intent(in) :: cur, dest
      integer(int64), intent(in) :: n
      logical, intent(in) :: first
      type(request_t), intent(inout) :: reqs(:)
      integer(int64) :: k, m
      integer :: r
      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         if (first) then
            call comm_isend_real_dp_array_n(comm, xr(k:k + m - 1), int(m), dest, 7, reqs(r))
         else
            call comm_isend_real_dp_array_n(comm, vis(k:k + m - 1, cur), int(m), dest, 7, reqs(r))
         end if
         k = k + m
      end do
   end subroutine send_block_dev

   subroutine recv_block_dev(comm, vis, c, n, source)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout), target :: vis(:, :)
      integer, intent(in) :: c, source
      integer(int64), intent(in) :: n
      integer(int64) :: k, m
      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call comm_recv_real_dp_array_n(comm, vis(k:k + m - 1, c), int(m), source, 7)
         k = k + m
      end do
   end subroutine recv_block_dev

   !
   ! Step 3's exchange. Rank r holds X(ia, P in P_r) for every i, with the
   ! rows for O_q (q /= r) packed in `ysend` from soff(q). It needs
   ! X(ia in O_r, P) for every P. The exchange runs nranks - 1 shifts: at
   ! step k rank r sends its segment to rank r + k and receives from rank
   ! r - k the columns P_{r-k} of its own X, which are contiguous in `xr`.
   ! As in ring_shift, each send is posted non-blocking before the blocking
   ! receive, so neither side can wait on the other. Messages go in
   ! MSG_CHUNK pieces. Sender and receiver compute the same count,
   ! bn(dest) nv |P_r|, so an empty one is skipped on both sides.
   !
   subroutine exchange_x(comm, rank, nranks, nv8, bn, plo, soff, ysend, xr, direct)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: rank, nranks
      integer, intent(in) :: bn(0:), plo(0:)
      integer(int64), intent(in) :: nv8
      integer(int64), intent(in) :: soff(0:)
      real(dp), intent(inout) :: ysend(:), xr(:)
      logical, intent(in) :: direct
      integer :: k, dest, src, nreq
      integer(int64) :: nsend, nrecv, roff, novr8
      type(request_t), allocatable :: reqs(:)
      novr8 = int(bn(rank), int64)*nv8
      if (.not. direct) then
         !$acc update self(ysend)
      end if
      do k = 1, nranks - 1
         dest = mod(rank + k, nranks)
         src = mod(rank - k + nranks, nranks)
         nsend = int(bn(dest), int64)*nv8*int(plo(rank + 1) - plo(rank), int64)
         nrecv = novr8*int(plo(src + 1) - plo(src), int64)
         roff = int(plo(src) - 1, int64)*novr8
         nreq = int((nsend + MSG_CHUNK - 1)/MSG_CHUNK)
         allocate (reqs(max(nreq, 1)))
         if (direct) then
            !$acc host_data use_device(ysend, xr)
            call send_seg_dev(comm, ysend, soff(dest), nsend, dest, reqs)
            call recv_seg_dev(comm, xr, roff, nrecv, src)
            !$acc end host_data
         else
            call send_seg(comm, ysend, soff(dest), nsend, dest, reqs)
            call recv_seg(comm, xr, roff, nrecv, src)
         end if
         if (nreq > 0) call waitall(reqs(1:nreq))
         if (.not. direct .and. nrecv > 0) then
            !$acc update device(xr(roff + 1:roff + nrecv))
         end if
         deallocate (reqs)
      end do
   end subroutine exchange_x

   !
   ! buf(off + 1 : off + n) to `dest` and from `source`, in MSG_CHUNK pieces,
   ! host addresses (send_seg, recv_seg) or device ones inside `host_data`
   ! (the _dev pair, by pic-mpi's specific names as send_block_dev does).
   !
   subroutine send_seg(comm, buffer, off, n, dest, reqs)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in), target :: buffer(:)
      integer(int64), intent(in) :: off, n
      integer, intent(in) :: dest
      type(request_t), intent(inout) :: reqs(:)
      integer(int64) :: k, m
      integer :: r
      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         call isend(comm, buffer(off + k:off + k + m - 1), dest, 8, reqs(r))
         k = k + m
      end do
   end subroutine send_seg

   subroutine recv_seg(comm, buffer, off, n, source)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: buffer(:)
      integer(int64), intent(in) :: off, n
      integer, intent(in) :: source
      integer(int64) :: k, m
      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call recv(comm, buffer(off + k:off + k + m - 1), int(m), source, 8)
         k = k + m
      end do
   end subroutine recv_seg

   subroutine send_seg_dev(comm, buffer, off, n, dest, reqs)
      type(comm_t), intent(in) :: comm
      real(dp), intent(in), target :: buffer(:)
      integer(int64), intent(in) :: off, n
      integer, intent(in) :: dest
      type(request_t), intent(inout) :: reqs(:)
      integer(int64) :: k, m
      integer :: r
      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         call comm_isend_real_dp_array_n(comm, buffer(off + k:off + k + m - 1), int(m), dest, 8, reqs(r))
         k = k + m
      end do
   end subroutine send_seg_dev

   subroutine recv_seg_dev(comm, buffer, off, n, source)
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout), target :: buffer(:)
      integer(int64), intent(in) :: off, n
      integer, intent(in) :: source
      integer(int64) :: k, m
      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call comm_recv_real_dp_array_n(comm, buffer(off + k:off + k + m - 1), int(m), source, 8)
         k = k + m
      end do
   end subroutine recv_seg_dev

   integer function k_group(no, nv, ncol) result(ng)
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, &
                         acc_property_free_memory
#endif
      integer, intent(in) :: no, nv, ncol
      integer, parameter :: NG_BEST = 4
#ifdef TRC_DEVICE
      integer(int64) :: free, per_i, fixed
      free = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                  acc_property_free_memory), int64)
      per_i = 8_int64*int(nv, int64)*int(ncol, int64)
      fixed = per_i + 16_int64*int(ncol, int64) + 512_int64*1024*1024
      ng = int(min(int(NG_BEST, int64), max(1_int64, (free - fixed)/per_i)))
#else
      ng = NG_BEST
#endif
      ng = max(1, min(ng, no))
   end function k_group

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   !! idx(mu, nu) = mu(mu-1)/2 + nu, mu >= nu, 1-based: trc_scf_df.F90's own
   !! `pack_idx`, duplicated here rather than used across the module
   !! boundary -- a `do concurrent` body may not call a procedure from
   !! another module (the same reason `ncart` above is its own copy and
   !! not trc_scf_df.F90's).
   pure integer function rimp2_pack_index(mu, nu) result(p)
      !$acc routine seq
      integer, intent(in) :: mu, nu
      p = mu*(mu - 1)/2 + nu
   end function rimp2_pack_index

   !
   ! For occupied i and the pairs j = jfirst.., one thread per (j, a), with a
   ! fastest and b summed inside: part(2 ja - 1) = w sum_b (ia|jb)^2 / D and
   ! part(2 ja) = w sum_b [(ia|jb)^2 - (ia|jb)(ib|ja)] / D, w = 1 for j = i
   ! and 2 otherwise -- the (j, i) pair is the same numbers again and is
   ! never computed. No reduction crosses threads; the host adds 2 nji
   ! numbers. Diagonal terms are in the sum: (ia|ia)^2 - (ia|ia)^2 = 0 for
   ! the same-spin part, as it must.
   !
   ! (ia|jb) is K(roff + a, coff + (j - jfirst) nv + b) and (ib|ja) is
   ! K(roff + b, coff + (j - jfirst) nv + a) = K_i^T(ja, b): with a fastest
   ! across threads, both are unit-stride reads. From K alone one of the
   ! two is always strided by nv across a warp, which held this kernel to a
   ! fifth of the device's bandwidth.
   !
   subroutine pair_energy(i, jfirst, no, nv, ldk, roff, coff, ldt, nji, eo, ev, ki, kt, part)
      integer, intent(in) :: i, jfirst, no, nv, ldk, roff, coff, ldt, nji
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: part(*)
      integer :: ja
      do concurrent(ja=1:nji)
         call pair_energy_body(i, jfirst, ja, no, nv, ldk, roff, coff, ldt, eo, ev, ki, kt, &
                               part(2*ja - 1), part(2*ja))
      end do
   end subroutine pair_energy

   pure subroutine pair_energy_body(i, jfirst, ja, no, nv, ldk, roff, coff, ldt, eo, ev, ki, kt, pos, pss)
      !$acc routine seq
      integer, intent(in) :: i, jfirst, ja, no, nv, ldk, roff, coff, ldt
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: pos, pss
      integer :: jj, j, a, b
      real(dp) :: t, u, r, eija, w
      jj = (ja - 1)/nv               ! j - jfirst
      j = jfirst + jj
      a = ja - jj*nv
      eija = eo(i) + eo(j) - ev(a)
      w = merge(1.0_dp, 2.0_dp, j == i)
      pos = 0.0_dp
      pss = 0.0_dp
      do b = 1, nv
         t = ki(roff + a, coff + jj*nv + b)  ! (ia|jb)
         u = kt(ja, b)                       ! (ib|ja)
         r = 1.0_dp/(eija - ev(b))
         pos = pos + t*t*r
         pss = pss + (t*t - t*u)*r
      end do
      pos = w*pos
      pss = w*pss
   end subroutine pair_energy_body

   !
   ! The mixed-spin (ab) block for occupied alpha orbital `i`: one thread
   ! per (a, J), B summed inside. `aJ` decomposes the same way `pair_energy`'s
   ! `ja` does (the row selector `a` fastest, the column group `J` next),
   ! here against a rectangle that is never folded by any symmetry, so
   ! there is no weight, no second read of a transpose, and nothing to
   ! subtract -- only (ia|JB)^2, never (ia|JB)(iB|Ja), which does not exist
   ! for two distinguishable electrons.
   !
   subroutine mixed_energy(i, nva, nob, nvb, ldk, roff, eo_a, ev_a, eo_b, ev_b, ki, part)
      integer, intent(in) :: i, nva, nob, nvb, ldk, roff
      real(dp), intent(in) :: eo_a(*), ev_a(nva), eo_b(nob), ev_b(nvb), ki(ldk, *)
      real(dp), intent(out) :: part(*)
      integer :: aj
      do concurrent(aj=1:nva*nob)
         call mixed_energy_body(i, aj, nva, nob, nvb, roff, eo_a, ev_a, eo_b, ev_b, ldk, ki, part(aj))
      end do
   end subroutine mixed_energy

   pure subroutine mixed_energy_body(i, aj, nva, nob, nvb, roff, eo_a, ev_a, eo_b, ev_b, ldk, ki, pab)
      !$acc routine seq
      integer, intent(in) :: i, aj, nva, nob, nvb, roff, ldk
      real(dp), intent(in) :: eo_a(*), ev_a(nva), eo_b(nob), ev_b(nvb), ki(ldk, *)
      real(dp), intent(out) :: pab
      integer :: jj, a, bb
      real(dp) :: t, eiaj
      jj = (aj - 1)/nva       ! J - 1 (0-based)
      a = aj - jj*nva
      eiaj = eo_a(i) + eo_b(jj + 1) - ev_a(a)
      pab = 0.0_dp
      do bb = 1, nvb
         t = ki(roff + a, jj*nvb + bb)   ! (iA|JB)
         pab = pab + t*t/(eiaj - ev_b(bb))
      end do
   end subroutine mixed_energy_body

end module trc_rimp2_driver
