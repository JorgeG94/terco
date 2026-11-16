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
! THE CODE
! --------
! trc_rimp2_run and trc_urimp2_run are sequencers: each reads as the list
! in THE STEPS and calls one module procedure per step. rimp2_stage_*
! are the restricted run's steps (the open-shell run reuses
! rimp2_stage_metric, rimp2_stage_aux_block, rimp2_stage_block_3c and
! rimp2_stage_solve_b, and the log lines rimp2_log_metric,
! rimp2_log_transform, rimp2_log_solve and rimp2_log_report);
! urimp2_stage_* are the open-shell run's own. What the stages hand each
! other apart from arrays is a few small host types: rimp2_layout_t
! (urimp2_layout_t), the sizes and this rank's share, fixed once the setup
! stages have run, and rimp2_ring_t, the ring step and group in hand with
! the pair energies. Memory: an array that outlives a stage is allocated
! and mapped (`enter data`/`exit data`) in the driver and passed down; a
! stage's own scratch is allocated, mapped and released inside that stage.
! One pair crosses a routine boundary: rimp2_stage_aux_block puts the
! auxiliary block's basis on the device (`to_device`) and the driver
! releases it at the end of the block. In the restricted run the
! half-transformed block and the K workspace are the driver's; in the
! open-shell run the half-transformed block and the same-spin K workspace
! are the stages' own (urimp2_stage_spin_transform, urimp2_stage_same_spin)
! and the mixed-spin K workspace is the driver's.
!
! Module trc_rimp2_driver, not trc_rimp2: the C entry is bound to the
! symbol trc_rimp2, and a binding label may not equal a module name
! (F2018 19.2); gfortran silently resolves the call to the wrong thing.
#include "inc/trc_dc_index.inc"

module trc_rimp2_driver
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_cart, only: ncart => cart_ncart_of
   use trc_integrals, only: trc_df_2c, trc_df_3c
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use pic_mpi_lib, only: comm_t, allreduce, bcast, isend, recv, request_t, waitall, MPI_SUM, &
                          comm_isend_real_dp_array_n, comm_recv_real_dp_array_n
   use trc_mpi_device, only: trc_mpi_device_direct
   use trc_step_table, only: step_table_t
   use trc_sph, only: sph_metric_cart_to_sph, sph_mo_to_cart, sph_aux3c_to_sph
   implicit none
   private

   public :: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, rimp2_driver_mo_to_cart
   public :: trc_rimp2_options_t

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
    !! number of algorithmic steps in the procedure
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
      real(dp) :: t_total = 0.0_dp   !! wall time of the whole call, read after the device queue drained (since()), seconds
      character(len=24), allocatable :: step_name(:)   !! the step table's names, as printed
      real(dp), allocatable :: step_t(:)               !! and each step's seconds (the job's slowest rank)
      character(len=200) :: message = ""
   end type trc_rimp2_result_t

   type :: trc_rimp2_options_t
    !! The knobs `trc_rimp2_run`/`trc_urimp2_run` used to take as separate
    !! optional arguments, bundled the way `trc_scf_options_t`
    !! (src/scf/trc_scf_settings.F90) bundles the SCF's own knobs. Plain
    !! values only; `la` (a caller-owned device resource, mutated and
    !! possibly reused across calls) stays a separate argument.
      integer :: nfrozen = 0
      !! Correlated occupied orbitals start at nfrozen + 1. Absent/default
      !! (0) correlates every occupied orbital.
      integer :: aux_block = 0
      !! Auxiliary functions per (mu nu|P) block. 0 (the default) means
      !! "the whole auxiliary basis at once" (`naux_c`); no real caller
      !! wants a block of exactly zero, so 0 doubles as the "not set"
      !! sentinel the old `present(aux_block)` check used.
      logical :: verbose = .false.
      !! One line per step with its wall time, rank 0 only.
      type(comm_t) :: comm
      !! Ranks, a copy of the caller's communicator handle; read only when
      !! `has_comm`. Without one, one rank.
      logical :: has_comm = .false.
   end type trc_rimp2_options_t

   type :: rimp2_layout_t
      !! The sizes of a restricted run and this rank's share of it. Written
      !! by `rimp2_stage_dimensions` (nao .. nblock, chained, rank, nranks)
      !! and `rimp2_stage_distribution` (lo .. nsend8 and the arrays), and
      !! read-only after that. Host values only. Never built with a structure
      !! constructor.
      integer :: nao = 0      !! Cartesian AOs (`b%nao`), whether `b` is spherical or not
      integer :: naux = 0     !! auxiliary functions: spherical when `aux` is
      integer :: nfz = 0      !! frozen orbitals
      integer :: no = 0       !! correlated occupied orbitals
      integer :: nv = 0       !! virtual orbitals, `nao - nocc`, pad columns included
      integer :: nblock = 0   !! auxiliary functions per (mu nu|P) block
      logical :: chained = .false.   !! B arrives pre-fitted (`packed_b`)
      integer :: rank = 0, nranks = 1
      integer :: lo = 0       !! first occupied orbital of this rank's block O_rank
      integer :: nb = 0       !! its size
      integer :: nmax = 0     !! the largest block of any rank
      integer :: novr = 0     !! `nb nv`, the rows of this rank's B
      integer :: npr = 0      !! columns of X this rank builds, |P_rank|
      integer(int64) :: nsend8 = 0   !! doubles in `ysend`
      integer, allocatable :: blo(:), bn(:)       !! every rank's occupied block: first orbital, size (0:nranks-1)
      integer, allocatable :: shlo(:), shhi(:)    !! every rank's auxiliary shells P_q (0:nranks-1)
      integer, allocatable :: plo(:)              !! where P_q starts among X's columns (0:nranks)
      integer(int64), allocatable :: soff(:)      !! where the rows for O_q start in `ysend` (0:nranks-1)
   end type rimp2_layout_t

   type :: rimp2_ring_t
      !! The pair loop of a restricted run: its shape, the ring step and the
      !! group of occupied orbitals in hand, and the pair energies summed so
      !! far. `rimp2_stage_pair_layout` sets nring .. ncolmax once;
      !! `rimp2_ring_step_range` sets qv, diag, ilo .. ldv at each step,
      !! `rimp2_group_range` sets i0 .. ncol at each group; `ring_shift`
      !! (through the driver) updates `cur`, which the driver starts at 0;
      !! the driver zeroes e_os and e_ss, `rimp2_stage_pair_group` adds this
      !! rank's pairs to them and `rimp2_stage_reduce_ranks` replaces them by
      !! the sum over ranks. The energies live here, not in the driver, so
      !! that `rimp2_stage_pair_group` stays within eight arguments.
      integer :: nring = 0    !! ring steps after the first
      integer :: ng = 0       !! occupied orbitals per K GEMM
      integer :: ldk = 0      !! rows of K, `ng nv`
      integer :: ncolmax = 0  !! columns of K and the width of `kt`, `nmax nv`
      integer :: qv = 0       !! the rank whose block of B is visiting at this ring step
      integer :: cur = 0      !! which of the two visiting buffers holds it
      logical :: diag = .false.   !! the first step: both blocks are this rank's own, and j >= i
      integer :: ilo = 0, ihi = 0, jlo = 0, jhi = 0   !! rows i and columns j of this ring step
      integer :: jbase = 0    !! the first orbital of block `qv`
      integer :: ldv = 0      !! leading dimension of the block the columns come from
      integer :: i0 = 0, i1 = 0   !! the group of occupied orbitals of this K GEMM
      integer :: jstart = 0   !! its first column's orbital
      integer :: ncol = 0     !! its columns
      real(dp) :: e_os = 0.0_dp, e_ss = 0.0_dp   !! this rank's pair energies so far; the job's after the reduction
   end type rimp2_ring_t

   type :: urimp2_layout_t
      !! The sizes of an unrestricted run (one rank), written by
      !! `urimp2_stage_dimensions` and read-only after that.
      integer :: nao = 0, naux = 0, nfz = 0
      integer :: noa = 0, nob = 0   !! correlated occupied orbitals per spin
      integer :: nva = 0, nvb = 0   !! virtual orbitals per spin
      integer :: nblock = 0         !! auxiliary functions per (mu nu|P) block
   end type urimp2_layout_t

contains

   recursive subroutine trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, opts, la, packed_b)
      !!
      !! `cmo` (nao, nao) and `eps` (nao) are the RHF orbitals, occupied first;
      !! `nocc` doubly occupied, `opts%nfrozen` of them left out of the
      !! correlation (see `trc_rimp2_options_t`'s own header for `opts`, optional
      !! throughout, as a whole; its `comm` splits the work over ranks). `pl` is
      !! the orbital basis's pair list. `la` optional, made here if absent.
      !!
      !! `packed_b`: the fitted B of a density-fitted J/K object
      !! (trc_scf_df.F90's `df_jk_t%bpk`, allocated once that object is set
      !! up), fitted on an auxiliary basis the CALLER has checked is the same
      !! one as `aux` (same JSON path/contents -- trc_capi.F90's trc_run_rimp2
      !! is the one library caller that does; `trc_rimp2_run` only checks the
      !! dimensions agree, not that the basis itself is really the same one).
      !! Passed as the array rather than the object so that post_hf does not
      !! depend on scf. When present, steps 1 (metric, Cholesky), 2 (the
      !! (mu nu|P) build) and 4 (B = X L^-T) are skipped entirely: `packed_b`
      !! is already B, (npack, naux) device-resident (packed mu >= nu, L^-T
      !! applied), and the occupied/virtual half-transform of step 3 run on it
      !! directly IS RI-MP2's own B -- L^-T acts only on the auxiliary index,
      !! the half-transform only on (mu, nu), and the two commute. `aux` is
      !! still used for its shell structure (the block loop below only needs
      !! auxiliary function COUNTS and offsets, identical either way).
      !!
      !! run a restricted ri-mp2 calculation
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      type(trc_rimp2_result_t), intent(out) :: res
      type(trc_rimp2_options_t), intent(in), optional :: opts
      type(trc_linalg_t), intent(inout), optional :: la
      real(dp), intent(in), optional, contiguous :: packed_b(:, :)
      !! (npack, naux), device-resident: see the header above.

      type(rimp2_layout_t) :: lay
      type(rimp2_ring_t) :: ring
      type(trc_rimp2_options_t) :: op
      type(trc_basis_t) :: aux_subset
      type(step_table_t) :: stab
      type(timer_type) :: tc, tw0, tw1
      logical :: ok, talk, direct
      integer :: s0, s1, p0, np, i0, step
      real(dp), allocatable :: jm(:, :), cocc(:, :), cvir(:, :)
      real(dp), allocatable :: tens(:, :, :), half(:, :, :)
      real(dp), allocatable :: xr(:), vis(:, :), ysend(:)
      real(dp), allocatable :: ki(:, :), kt(:, :), part(:), eo(:), ev(:)
      integer(int64) :: novr8, nv8, nao8, naux8

      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call trc_rimp2_run(b, aux, pl, nocc, cmo, eps, res, opts, la_local, packed_b)
            call la_local%release()
         end block
         return
      end if
      if (present(opts)) op = opts
      ! `cmo`/`eps` are declared (b%nao, b%nao)/(b%nao): ALWAYS the
      ! Cartesian AO space, see rimp2_stage_dimensions.
      call rimp2_stage_dimensions(b, aux, nocc, op, packed_b, res, lay, ok)
      if (.not. ok) return
      call rimp2_stage_distribution(aux, lay)
      ! Create-once: a caller's handles are kept, the workspace only grows.
      call la%init(lay%naux)

      talk = op%verbose .and. lay%rank == 0
      if (talk) call global_logger%performance("  rimp2: starting, no "//to_string(lay%no)//" nv "//to_string(lay%nv)// &
                                               " naux "//to_string(lay%naux)//", ranks "//to_string(lay%nranks))
      novr8 = lay%novr
      nv8 = lay%nv
      nao8 = lay%nao
      naux8 = lay%naux
      ! The metric and its factor run in full on every rank; the integrals are
      ! split by auxiliary shell.
      call stab%init("rimp2", STEP_NAME, replicated_steps=[ST_METRIC, ST_CHOL])
      call la%sync()
      call tw0%start()
      ! --- 1. the metric and its Cholesky factor (skipped when chained) ----
      ! rimp2_stage_metric opens and closes ST_METRIC and ST_CHOL itself.
      if (.not. lay%chained) then
         allocate (jm(lay%naux, lay%naux))
         !$acc enter data create(jm)
         call rimp2_stage_metric(aux, la, stab, lay%naux, jm)
         if (talk) call rimp2_log_metric("rimp2", lay%naux, since(la, tw0))
      else if (talk) then
         call global_logger%performance("  rimp2: metric, Cholesky and (mu nu|P) skipped -- chained off --df's fit")
      end if

      ! --- 3. (ia|P) for every correlated i and this rank's P, an auxiliary
      ! block at a time. The rows for O_q, q /= rank, go to `ysend`, one
      ! contiguous (bn(q) nv) x |P_rank| segment per q starting at soff(q);
      ! the rows for O_rank go straight into `xr`, in the columns they keep.
      direct = trc_mpi_device_direct()
      allocate (cocc(lay%nao, lay%no), cvir(lay%nao, lay%nv), xr(max(novr8*naux8, 1_int64)), ysend(max(lay%nsend8, 1_int64)))
      cocc = cmo(:, lay%nfz + 1:nocc)
      cvir = cmo(:, nocc + 1:lay%nao)
      !$acc enter data copyin(cocc, cvir) create(xr, ysend)
      s0 = lay%shlo(lay%rank)
      do while (s0 <= lay%shhi(lay%rank))
         call rimp2_stage_aux_block(aux, s0, lay%shhi(lay%rank), lay%nblock, aux_subset, s1, p0, np)
         allocate (tens(lay%nao, lay%nao, np), half(lay%no, lay%nao, np))
         !$acc enter data create(tens, half)
         call stab%step_on(la, ST_3C, tc)
         if (lay%chained) then
            ! B, already L^-T applied (and, when `aux` is spherical, the
            ! P-fold too -- `dfjk%setup` does both before packing),
            ! packed mu >= nu: unpack this block's columns straight into
            ! the dense tensor the half-transform below reads either way
            ! -- no (mu nu|P) build, and no P-fold, at all.
            call rimp2_unpack_block(lay%nao, size(packed_b, 1), lay%naux, np, p0, packed_b, tens)
         else
            call rimp2_stage_block_3c(b, pl, aux_subset, la, tens)
         end if
         call stab%step_off(la, ST_3C, tc)
         call stab%step_on(la, ST_TRANS, tc)
         call rimp2_stage_half_occupied(la, lay, np, cocc, tens, half)
         call rimp2_stage_half_virtual(la, lay, p0, np, cvir, half, xr, ysend)
         call stab%step_off(la, ST_TRANS, tc)
         call stab%add_flop(ST_TRANS, int(np, int64)*(2*int(lay%no, int64)*nao8*nao8 + 2*nv8*int(lay%no, int64)*nao8))
         !$acc exit data delete(tens, half)
         deallocate (tens, half)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      if (lay%nranks > 1) then
         call stab%step_on(la, ST_XCHG, tc)
         call exchange_x(op%comm, nv8, lay%bn, lay%plo, lay%soff, ysend, xr, direct)
         call stab%step_off(la, ST_XCHG, tc)
      end if
      !$acc exit data delete(ysend)
      deallocate (ysend)
      if (talk) call rimp2_log_transform("rimp2", stab, .true.)

      ! --- 4. B = X L^-T: solve B L^T = X from the right, this rank's rows ---
      ! (skipped when chained: `xr` is already built from B, which has
      ! L^-T applied already -- see this routine's own `packed_b` comment)
      ! rimp2_stage_solve_b opens and closes ST_TRSM itself.
      if (.not. lay%chained) then
         call rimp2_stage_solve_b(la, stab, lay%naux, lay%novr, jm, xr)
         if (talk) call rimp2_log_solve("rimp2", stab)
      end if
      call la%sync()
      call tw1%start()

      ! --- 5. the pair energies, over the ring ------------------------------
      ! The visiting blocks need two buffers, the one in use and the one
      ! arriving; one rank needs neither.
      if (talk .and. lay%nranks > 1) call rimp2_log_exchange_mode(direct)
      if (lay%nranks > 1) then
         allocate (vis(int(lay%nmax, int64)*nv8*naux8, 2))
         !$acc enter data create(vis)
      end if
      call rimp2_stage_pair_layout(lay, ring)
      if (talk) call rimp2_log_k_groups(ring%ng, lay%nranks)
      allocate (ki(ring%ldk, ring%ncolmax), kt(ring%ncolmax, lay%nv), part(2*ring%ncolmax), eo(lay%no), ev(lay%nv))
      eo = eps(lay%nfz + 1:nocc)
      ev = eps(nocc + 1:lay%nao)
      !$acc enter data create(ki, kt, part) copyin(eo, ev)
      ring%e_os = 0.0_dp
      ring%e_ss = 0.0_dp
      ring%cur = 0
      do step = 0, ring%nring
         call rimp2_ring_step_range(lay, ring, step)
         if (step > 0) then
            call stab%step_on(la, ST_RING, tc)
            call ring_shift(op%comm, xr, vis, ring%cur, &
                            int(lay%bn(mod(lay%rank + step - 1, lay%nranks)), int64)*nv8*naux8, &
                            int(lay%bn(ring%qv), int64)*nv8*naux8, step == 1, direct)
            call stab%step_off(la, ST_RING, tc)
         end if
         ! Groups of `ng` occupied orbitals: one GEMM per group, K(i a, jb)
         ! for i in the group and j from the group's first column on, so
         ! cuBLAS gets (ng nv) rows instead of nv. On the diagonal block it
         ! also computes the j < i corner inside the group, which the pair
         ! kernel skips.
         do i0 = ring%ilo, ring%ihi, ring%ng
            call rimp2_group_range(lay, ring, i0)
            if (ring%ncol <= 0) cycle
            if (ring%diag) then
               call rimp2_stage_k_gemm_own(la, stab, lay, ring, xr, ki)
            else
               call rimp2_stage_k_gemm_visiting(la, stab, lay, ring, xr, vis, ki)
            end if
            call rimp2_stage_pair_group(la, stab, ring, ki, kt, part, eo, ev)
         end do
      end do
      !$acc exit data delete(ki, kt, part, eo, ev, xr, cocc, cvir)
      if (.not. lay%chained) then
         !$acc exit data delete(jm)
      end if
      if (lay%nranks > 1) then
         !$acc exit data delete(vis)
      end if
      call rimp2_stage_reduce_ranks(op%comm, stab, lay, ring)
      if (talk) call rimp2_log_report("rimp2", la, stab, tw0, tw1, lay%rank, lay%nranks)
      call rimp2_store_result(res, stab, since(la, tw0), ring%e_os, ring%e_ss)
   end subroutine trc_rimp2_run

   subroutine rimp2_unpack_block(nao, npack, naux, np, p0, packed_b, tens)
      !! tens(:, :, q) = B's column p0 - 1 + q, packed mu >= nu, unpacked to
      !! the dense (nao, nao) block the half-transform reads; on the device.
      integer, intent(in) :: nao, npack, naux, np, p0
      real(dp), intent(in) :: packed_b(npack, naux)
      real(dp), intent(out) :: tens(nao, nao, np)
      integer :: q, mu, nu
      do concurrent(q=1:np, mu=1:nao, nu=1:nao)
         tens(mu, nu, q) = packed_b(rimp2_pack_index(max(mu, nu), min(mu, nu)), p0 - 1 + q)
      end do
   end subroutine rimp2_unpack_block

   subroutine rimp2_driver_mo_to_cart(b, csph, epssph, ccart, epscart)
      !!
      !! `trc_rimp2_run`/`trc_urimp2_run`'s own fold for a spherical orbital
      !! basis: `csph`/`epssph` are an SCF's own (b%nao_sph, b%nao_sph)/
      !! (b%nao_sph) orbitals and energies (`trc_scf_result_t%cmo`/`eps`,
      !! already at the SCF's own `scf_space_nao` size); `ccart`/`epscart`
      !! come back (b%nao, b%nao)/(b%nao), the Cartesian AO space RI-MP2's own
      !! three-centre tensor is built in.
      !!
      !! Every real MO is `sph_mo_to_cart`'s T C_sph, trc_sph.F90's one-sided
      !! fold (an MO's AO row crosses the spherical/Cartesian boundary, its
      !! column -- the MO label -- does not: MO-space quantities do not care
      !! which AO representation produced them). The Cartesian AO count is
      !! always at least the spherical MO count (a Cartesian shell never has
      !! fewer functions than its spherical one), so `b%nao - b%nao_sph`
      !! columns are left over; they are filled with a ZERO coefficient
      !! vector and `PAD_EPS`, a finite energy well outside any real
      !! spectrum. `rimp2_stage_dimensions`' `nv = nao - nocc` (nao = b%nao here)
      !! then runs over these padded columns exactly as it would a real
      !! virtual, but every (ia|jb) with a or b among them is a GEMM
      !! contraction against a zero `C_vir` column and comes out exactly
      !! zero -- PAD_EPS only has to keep `1/(e_i + e_j - e_a - e_b)` FINITE
      !! (0 * finite = 0), not physically meaningful, since it is always
      !! multiplied by that zero. This is what lets `trc_rimp2_run` itself
      !! stay unchanged, Cartesian-interface, for a spherical caller.
      !!
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
      ! rimp2_stage_dimensions' `nv = nao - nocc` then costs real GEMM
      ! flops on -- see the comment in `rimp2_stage_dimensions` for the measured
      ! fraction.
      if (nc > ns) ccart(:, ns + 1:nc) = 0.0_dp
      epscart(1:ns) = epssph
      if (nc > ns) epscart(ns + 1:nc) = PAD_EPS
   end subroutine rimp2_driver_mo_to_cart

   recursive subroutine trc_urimp2_run(b, aux, pl, nalpha, nbeta, cmo, eps, res, opts, la)
      !!
      !! RI-MP2 for an UNRESTRICTED (UHF or UKS) reference, on the device, one
      !! rank only.
      !!
      !! The metric, its Cholesky factor and (mu nu|P) do not know about spin --
      !! they are built exactly once, as trc_rimp2_run's own steps 1-2 -- and
      !! transformed TWICE, into B_a (from C_occ_a, C_vir_a) and B_b (from
      !! C_occ_b, C_vir_b), by the same GEMM-plus-strided-batch half-transform
      !! and the same triangular solve as the restricted routine's steps 3-4,
      !! called once per spin.
      !!
      !! Three pair blocks follow, all from B_a and B_b:
      !!
      !!   E_aa = 1/2 sum_{ij in a, ab in a} (ia|jb)[(ia|jb) - (ib|ja)] / D
      !!   E_bb = 1/2 sum_{ij in b, ab in b} (ia|jb)[(ia|jb) - (ib|ja)] / D
      !!   E_ab =     sum_{i in a, J in b, a in a, B in b} (ia|JB)^2 / D
      !!
      !! E_aa and E_bb reuse trc_rimp2_run's own K GEMM and `pair_energy` kernel
      !! verbatim -- same occupied-pair triangle (j >= i), same weight (2 off
      !! the diagonal, 1 on it) -- and only the same-spin half of what that
      !! kernel returns (`pss`) is kept; the opposite-spin half it also computes
      !! (`pos`) has no meaning for two electrons of the SAME spin and is
      !! discarded. The overall 1/2 is the one difference from the restricted
      !! energy: expanding the standard spin-orbital result
      !!
      !!     E_aa = 1/4 sum_{ijab} [(ia|jb) - (ib|ja)]^2 / D     (i, j, a, b all
      !!                                                           independent,
      !!                                                           no triangle)
      !!
      !! and using (ia|jb)^2 = (ib|ja)^2 over that full range gives exactly
      !! 1/2 sum_{ijab} (ia|jb)[(ia|jb)-(ib|ja)]/D -- restricted's own full
      !! double sum, weighted down to the triangle the same way restricted's
      !! E_ss is -- hence the 1/2 that sits outside the triangle sum in
      !! `urimp2_stage_same_spin`. See
      !! metalquicha's run_czt_uri_mp2 (mqc_czt_mp2.f90) for the same algebra
      !! checked against a conventional (non-fitted) transform.
      !!
      !! E_ab has no such symmetry -- i (alpha) and J (beta) are already
      !! distinguishable electrons, so there is no exchange term and no pair to
      !! fold in half: the full rectangle of i against J, one GEMM
      !! K(iA, JB) = sum_Q B_a(iA, Q) B_b(JB, Q) per group of alpha occupied
      !! orbitals, and a kernel that only needs K itself (no transpose, no
      !! second read) since there is nothing to subtract.
      !!
      !! MULTI-RANK: refused. The restricted ring distributes one occupied
      !! block of B and lets it meet itself (aa) through a ring of nranks/2
      !! steps; an open-shell run would need that ring run independently for aa
      !! (over alpha's block) and bb (over beta's, a different size whenever
      !! nalpha /= nbeta), plus a THIRD distribution for ab that is not a ring
      !! at all -- alpha's block never needs to meet another alpha block, only
      !! every beta block, so the "pass one step left" structure the ring
      !! exists for buys nothing there. Three different distributions behind
      !! one call is not a clean generalisation of the one the restricted path
      !! has, so this round runs one rank only and `urimp2_stage_dimensions` refuses anything else,
      !! with a message that says so.
      !!
      !! `cmo`/`eps`: alpha in column/row 1, beta in 2 (trc_scf_result_t's layout). `res%cmo`/`res%eps`
      !! pass straight through for a Cartesian, nspin = 2 result only; fold a spherical one first.
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      integer, intent(in) :: nalpha, nbeta
      real(dp), intent(in) :: cmo(b%nao, b%nao, 2), eps(b%nao, 2)
      type(trc_rimp2_result_t), intent(out) :: res
      type(trc_rimp2_options_t), intent(in), optional :: opts
      type(trc_linalg_t), intent(inout), optional :: la

      type(urimp2_layout_t) :: lay
      type(trc_rimp2_options_t) :: op
      type(timer_type) :: tw0, tw1, tc
      logical :: talk, ok
      type(trc_basis_t) :: aux_subset
      type(step_table_t) :: stab
      integer :: s0, s1, p0, np, i0, i1, i, roff, ng_ab, ldk_ab
      real(dp), allocatable :: jm(:, :)
      real(dp), allocatable :: cocc_a(:, :), cvir_a(:, :), cocc_b(:, :), cvir_b(:, :)
      real(dp), allocatable :: tens(:, :, :)
      real(dp), allocatable :: xr_a(:), xr_b(:)
      real(dp), allocatable :: ki_ab(:, :), part_ab(:)
      real(dp), allocatable :: eo_a(:), ev_a(:), eo_b(:), ev_b(:)
      real(dp) :: e_aa, e_bb, e_ab
      integer(int64) :: novr_a8, novr_b8, nva8, nob8, nvb8

      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call trc_urimp2_run(b, aux, pl, nalpha, nbeta, cmo, eps, res, opts, la_local)
            call la_local%release()
         end block
         return
      end if
      if (present(opts)) op = opts
      ! `cmo`/`eps` are the Cartesian-space, fold-and-pad result of
      ! `rimp2_driver_mo_to_cart` when `b` is spherical -- see `trc_rimp2_run`'s
      ! own header.
      call urimp2_stage_dimensions(b, aux, nalpha, nbeta, op, res, lay, ok)
      if (.not. ok) return
      novr_a8 = int(lay%noa, int64)*lay%nva
      novr_b8 = int(lay%nob, int64)*lay%nvb
      nva8 = lay%nva
      nvb8 = lay%nvb
      nob8 = lay%nob
      call la%init(lay%naux)

      talk = op%verbose   ! one rank, rank 0
      call stab%init("rimp2u", STEP_NAME)
      call la%sync()
      call tw0%start()

      ! --- 1. the metric and its Cholesky factor (spin-independent) --------
      ! rimp2_stage_metric opens and closes ST_METRIC and ST_CHOL itself.
      allocate (jm(lay%naux, lay%naux))
      !$acc enter data create(jm)
      call rimp2_stage_metric(aux, la, stab, lay%naux, jm)
      if (talk) call rimp2_log_metric("rimp2u", lay%naux, since(la, tw0))

      ! --- 2-3. (mu nu|P) once per aux block, (ia|P) for BOTH spins from it -
      allocate (cocc_a(lay%nao, lay%noa), cvir_a(lay%nao, lay%nva), xr_a(max(novr_a8*lay%naux, 1_int64)))
      allocate (cocc_b(lay%nao, lay%nob), cvir_b(lay%nao, lay%nvb), xr_b(max(novr_b8*lay%naux, 1_int64)))
      cocc_a = cmo(:, lay%nfz + 1:nalpha, 1)
      cvir_a = cmo(:, nalpha + 1:lay%nao, 1)
      cocc_b = cmo(:, lay%nfz + 1:nbeta, 2)
      cvir_b = cmo(:, nbeta + 1:lay%nao, 2)
      !$acc enter data copyin(cocc_a, cvir_a, cocc_b, cvir_b) create(xr_a, xr_b)
      s0 = 1
      do while (s0 <= aux%nshell)
         call rimp2_stage_aux_block(aux, s0, aux%nshell, lay%nblock, aux_subset, s1, p0, np)
         allocate (tens(lay%nao, lay%nao, np))
         !$acc enter data create(tens)
         call stab%step_on(la, ST_3C, tc)
         call rimp2_stage_block_3c(b, pl, aux_subset, la, tens)
         call stab%step_off(la, ST_3C, tc)
         ! urimp2_stage_spin_transform opens and closes ST_TRANS itself.
         call urimp2_stage_spin_transform(la, stab, p0, cocc_a, cvir_a, tens, xr_a)
         call urimp2_stage_spin_transform(la, stab, p0, cocc_b, cvir_b, tens, xr_b)
         !$acc exit data delete(tens)
         deallocate (tens)
         call aux_subset%release()
         s0 = s1 + 1
      end do
      if (talk) call rimp2_log_transform("rimp2u", stab, .false.)

      ! --- 4. B = X L^-T, one triangular solve per spin ---------------------
      ! rimp2_stage_solve_b opens and closes ST_TRSM itself.
      call rimp2_stage_solve_b(la, stab, lay%naux, lay%noa*lay%nva, jm, xr_a)
      call rimp2_stage_solve_b(la, stab, lay%naux, lay%nob*lay%nvb, jm, xr_b)
      !$acc exit data delete(jm)
      if (talk) call rimp2_log_solve("rimp2u", stab)
      call la%sync()
      call tw1%start()

      ! --- 5a/5b. the like-spin pair energies, restricted's own triangle ---
      allocate (eo_a(lay%noa), ev_a(lay%nva), eo_b(lay%nob), ev_b(lay%nvb))
      eo_a = eps(lay%nfz + 1:nalpha, 1)
      ev_a = eps(nalpha + 1:lay%nao, 1)
      eo_b = eps(lay%nfz + 1:nbeta, 2)
      ev_b = eps(nbeta + 1:lay%nao, 2)
      !$acc enter data copyin(eo_a, ev_a, eo_b, ev_b)
      call urimp2_stage_same_spin(la, stab, lay%naux, xr_a, eo_a, ev_a, e_aa)
      call urimp2_stage_same_spin(la, stab, lay%naux, xr_b, eo_b, ev_b, e_bb)

      ! --- 5c. the mixed-spin block: the full rectangle, no exchange term --
      ng_ab = k_group(lay%noa, lay%nva, lay%nob*lay%nvb)
      ldk_ab = ng_ab*lay%nva
      allocate (ki_ab(ldk_ab, lay%nob*lay%nvb), part_ab(lay%nva*lay%nob))
      !$acc enter data create(ki_ab, part_ab)
      e_ab = 0.0_dp
      do i0 = 1, lay%noa, ng_ab
         i1 = min(lay%noa, i0 + ng_ab - 1)
         call urimp2_stage_mixed_gemm(la, stab, lay, i0, i1, xr_a, xr_b, ki_ab)
         do i = i0, i1
            roff = (i - i0)*lay%nva
            call stab%step_on(la, ST_PAIR, tc)
            call mixed_energy(i, lay%nva, lay%nob, lay%nvb, ldk_ab, roff, eo_a, ev_a, eo_b, ev_b, ki_ab, part_ab)
            !$acc update self(part_ab(1:lay%nva*lay%nob))
            call stab%step_off(la, ST_PAIR, tc)
            call stab%add_flop(ST_PAIR, 4_int64*nva8*nob8*nvb8)
            e_ab = e_ab + sum(part_ab(1:lay%nva*lay%nob))
         end do
      end do
      !$acc exit data delete(ki_ab, part_ab)
      deallocate (ki_ab, part_ab)

      !$acc exit data delete(eo_a, ev_a, eo_b, ev_b, xr_a, xr_b, cocc_a, cvir_a, cocc_b, cvir_b)
      if (talk) call rimp2_log_report("rimp2u", la, stab, tw0, tw1, 0, 1)
      call rimp2_store_result(res, stab, since(la, tw0), e_ab, e_aa + e_bb)
   end subroutine trc_urimp2_run

   subroutine rimp2_stage_dimensions(b, aux, nocc, op, packed_b, res, lay, ok)
      !! The sizes of a restricted run, from the two bases, `nocc` and `op`,
      !! and the checks that decide whether it can run at all. `ok` is false,
      !! with `res%message` set, when there is nothing to correlate, a basis
      !! carries angular momentum above what the density-fitting kernels were
      !! generated for, or the pre-fitted `packed_b` does not match the bases.
      type(trc_basis_t), intent(in) :: b, aux
      integer, intent(in) :: nocc
      type(trc_rimp2_options_t), intent(in) :: op
      real(dp), intent(in), optional, contiguous :: packed_b(:, :)
      type(trc_rimp2_result_t), intent(inout) :: res
      type(rimp2_layout_t), intent(inout) :: lay
      logical, intent(out) :: ok
      integer :: npack_chk

      ok = .false.
      lay%nao = b%nao
      ! `cmo`/`eps` of `trc_rimp2_run` are declared (b%nao, b%nao)/(b%nao): ALWAYS the
      ! Cartesian AO space, whether `b` itself is spherical or not. A
      ! caller whose SCF ran spherical (`b%spherical`) does not hand this
      ! its raw (b%nao_sph, b%nao_sph) orbitals -- that would read past
      ! the end of a smaller array under sequence association -- but
      ! `rimp2_driver_mo_to_cart`'s own fold-and-pad result instead: every real
      ! MO transformed into the Cartesian AO space (C_cart = T C_sph, MO
      ! space itself does not care which AO representation it came from),
      ! padded out to `b%nao` columns with zero-coefficient, PAD_EPS
      ! "virtuals" that contribute exactly zero to every pair energy.
      ! `trc_rimp2_run` never reads `b%spherical` itself.
      !
      ! The auxiliary basis `aux` is a separate question: it follows its
      ! own `aux%spherical` (trc_capi.F90's `basis_json` sets it to match
      ! the orbital basis, pyscf's own `auxmol.cart = mol.cart` rule, but
      ! `trc_rimp2_run` only ever reads `aux%spherical`). `aux%nao`, the
      ! CARTESIAN count, is every kernel call and block-walk
      ! (always Cartesian); `naux`, folded when `aux%spherical`, is every
      ! other dimension -- the metric, its Cholesky factor, `xr`'s own
      ! width, and (via `packed_b`'s column count) what a chained call's own aux count
      ! is checked against.
      lay%naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      lay%nfz = max(op%nfrozen, 0)
      lay%no = nocc - lay%nfz
      ! For a spherical caller, `nao` (Cartesian) is bigger than the
      ! true spherical MO count `b%nao_sph` -- rimp2_driver_mo_to_cart's own
      ! header explains why `trc_rimp2_run` is never told that -- so `nv`
      ! runs over `b%nao - b%nao_sph` zero-coefficient pad
      ! "virtuals" that cost real GEMM flops and contribute nothing.
      ! Measured from shell counts alone (no kernel change, not
      ! profiled): water/cc-pVDZ(sph) is nv=20 Cartesian vs. 19 true
      ! (1 pad column, one O d-shell's 6-vs-5), formamide/cc-pVDZ(sph)
      ! is nv=48 vs. 45 (3 pad columns, one d-shell per heavy atom) --
      ! a 5-7% larger nv, which an O(nv^2) GEMM (the dominant scaling
      ! of the half-transform/pair-energy steps) turns into an
      ! 11-14% larger flop count; any O(nv^3) term costs proportionally
      ! more still. Don't change the kernels for this -- it is the
      ! price of keeping `trc_rimp2_run` Cartesian-only and reusing it
      ! unchanged for a spherical caller (rimp2_driver_mo_to_cart's own
      ! header).
      lay%nv = lay%nao - nocc
      res%naux = lay%naux
      res%nocc = lay%no
      res%nvir = lay%nv
      if (lay%no < 1 .or. lay%nv < 1) then
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
      lay%chained = present(packed_b)
      if (lay%chained) then
         npack_chk = lay%nao*(lay%nao + 1)/2
         if (size(packed_b, 1) /= npack_chk .or. size(packed_b, 2) /= lay%naux) then
            res%message = "trc_rimp2: the pre-fitted B does not match this orbital/auxiliary basis"
            return
         end if
      end if
      lay%nblock = aux%nao
      if (op%aux_block > 0) lay%nblock = op%aux_block
      lay%rank = 0
      lay%nranks = 1
      if (op%has_comm) then
         lay%rank = op%comm%rank()
         lay%nranks = op%comm%size()
      end if
      ok = .true.
   end subroutine rimp2_stage_dimensions

   subroutine rimp2_stage_distribution(aux, lay)
      !! What each rank owns: its block of occupied orbitals, its auxiliary
      !! shells and where their columns sit in X, and where the rows for each
      !! other rank start in `ysend`. Host integers only.
      type(trc_basis_t), intent(in) :: aux
      type(rimp2_layout_t), intent(inout) :: lay
      integer :: i, q, ish, naux_c
      integer(int64) :: nsend8

      ! This rank's occupied block, O_rank = lo .. lo + nb - 1, and every
      ! rank's, which the ring needs to size what it receives.
      allocate (lay%blo(0:lay%nranks - 1), lay%bn(0:lay%nranks - 1))
      do i = 0, lay%nranks - 1
         lay%bn(i) = lay%no/lay%nranks + merge(1, 0, i < mod(lay%no, lay%nranks))
         lay%blo(i) = i*(lay%no/lay%nranks) + min(i, mod(lay%no, lay%nranks)) + 1
      end do
      lay%lo = lay%blo(lay%rank)
      lay%nb = lay%bn(lay%rank)
      lay%nmax = maxval(lay%bn)
      lay%novr = lay%nb*lay%nv
      ! Every rank's auxiliary shells, P_q = shells shlo(q) .. shhi(q) (empty
      ! when shhi < shlo): a shell goes to rank q while fewer than
      ! naux_c (q + 1) / nranks Cartesian functions come before it, so the
      ! integral work, which is Cartesian, is split evenly. plo(q) is where
      ! P_q starts among X's columns (spherical when `aux` is) and
      ! plo(q + 1) - plo(q) is its width.
      naux_c = aux%nao
      allocate (lay%shlo(0:lay%nranks - 1), lay%shhi(0:lay%nranks - 1), lay%plo(0:lay%nranks))
      ish = 1
      do q = 0, lay%nranks - 1
         lay%shlo(q) = ish
         do while (ish <= aux%nshell)
            if (int(aux%sh_ao(ish) - 1, int64) >= (int(naux_c, int64)*(q + 1))/lay%nranks) exit
            ish = ish + 1
         end do
         lay%shhi(q) = ish - 1
         if (lay%shlo(q) <= aux%nshell) then
            lay%plo(q) = merge(aux%sh_ao_sph(lay%shlo(q)), aux%sh_ao(lay%shlo(q)), aux%spherical)
         else
            lay%plo(q) = lay%naux + 1
         end if
      end do
      lay%plo(lay%nranks) = lay%naux + 1
      lay%npr = lay%plo(lay%rank + 1) - lay%plo(lay%rank)
      allocate (lay%soff(0:lay%nranks - 1))
      nsend8 = 0
      do q = 0, lay%nranks - 1
         lay%soff(q) = nsend8
         if (q /= lay%rank) nsend8 = nsend8 + int(lay%bn(q), int64)*int(lay%nv, int64)*int(lay%npr, int64)
      end do
      lay%nsend8 = nsend8
   end subroutine rimp2_stage_distribution

   subroutine rimp2_stage_metric(aux, la, stab, naux, jm)
      !! The metric (P|Q) of the auxiliary basis, folded to spherical when
      !! `aux` is, and its Cholesky factor, in place, on the device. Opens
      !! and closes ST_METRIC and then ST_CHOL, and counts the factor's flops.
      type(trc_basis_t), intent(in) :: aux
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      integer, intent(in) :: naux
      real(dp), intent(inout) :: jm(naux, naux)
      real(dp), allocatable :: jmc(:, :)
      type(timer_type) :: tc
      integer(int64) :: naux8

      naux8 = naux
      call stab%step_on(la, ST_METRIC, tc)
      if (aux%spherical) then
         allocate (jmc(aux%nao, aux%nao))
         !$acc enter data create(jmc)
         call trc_df_2c(aux, jmc)
         call sph_metric_cart_to_sph(aux, la, jmc, jm)
         !$acc exit data delete(jmc)
         deallocate (jmc)
      else
         call trc_df_2c(aux, jm)
      end if
      call stab%step_off(la, ST_METRIC, tc)
      call stab%step_on(la, ST_CHOL, tc)
      call la%potrf(naux, jm, naux)
      call stab%step_off(la, ST_CHOL, tc)
      call stab%add_flop(ST_CHOL, naux8**3/3)
   end subroutine rimp2_stage_metric

   subroutine rimp2_stage_aux_block(aux, s0, slast, nblock, aux_subset, s1, p0, np)
      !! The next auxiliary block: shells `s0` to `s1`, the most that hold
      !! at most `nblock` functions (and at least one shell), up to the last
      !! shell `slast` the caller owns. `aux_subset` is that range as a basis
      !! of its own, on the device; the caller releases it. `p0` and `np` are the
      !! column of X the block starts at and its width, spherical when `aux` is.
      type(trc_basis_t), intent(in) :: aux
      integer, intent(in) :: s0, slast, nblock
      type(trc_basis_t), intent(inout) :: aux_subset
      integer, intent(out) :: s1, p0, np
      integer :: ish, p0c, npc

      ! Shells s0..s1 hold at most nblock functions (and at least one shell).
      s1 = s0
      do ish = s0 + 1, slast
         if (aux%sh_ao(ish) + ncart(aux%sh_l(ish)) - aux%sh_ao(s0) > nblock) exit
         s1 = ish
      end do
      p0c = aux%sh_ao(s0)
      npc = aux%sh_ao(s1) + ncart(aux%sh_l(s1)) - p0c
      call aux%subset(s0, s1, aux_subset)
      call aux_subset%to_device()
      ! Effective (spherical, when `aux%spherical`) offset/width this
      ! block's columns land at in X -- `aux_subset%nao_sph` is
      ! exactly this block's own spherical width, since `subset` never
      ! cuts a shell in two.
      if (aux%spherical) then
         p0 = aux%sh_ao_sph(s0)
         np = aux_subset%nao_sph
      else
         p0 = p0c
         np = npc
      end if
   end subroutine rimp2_stage_aux_block

   subroutine rimp2_stage_block_3c(b, pl, aux_subset, la, tens)
      !! (mu nu|P) of one auxiliary block, on the device: `tens` is
      !! (nao, nao, np), np the block's width, spherical when the basis is.
      type(trc_basis_t), intent(in) :: b, aux_subset
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(inout) :: tens(b%nao, b%nao, merge(aux_subset%nao_sph, aux_subset%nao, aux_subset%spherical))
      real(dp), allocatable :: tcart3(:, :, :)

      if (aux_subset%spherical) then
         ! (mu nu|Pc), built Cartesian, folded into `tens` on the P
         ! index alone -- same `sph_aux3c_to_sph` `df_jk_setup` folds
         ! its own block through, before the half-transform
         ! reads it (which never needs to know the difference).
         allocate (tcart3(b%nao, b%nao, aux_subset%nao))
         !$acc enter data create(tcart3)
         call trc_df_3c(b, pl, aux_subset, tcart3)
         call sph_aux3c_to_sph(aux_subset, b%nao, tcart3, tens, la)
         !$acc exit data delete(tcart3)
         deallocate (tcart3)
      else
         call trc_df_3c(b, pl, aux_subset, tens)
      end if
   end subroutine rimp2_stage_block_3c

   subroutine rimp2_stage_half_occupied(la, lay, np, cocc, tens, half)
      !! The occupied side of the half transform of one auxiliary block of
      !! width `np`: half(i, nu q) = sum_mu C_occ(mu, i) T(mu, nu q), one GEMM
      !! for every correlated i and every q of the block. `half` is the
      !! block's work space.
      type(trc_linalg_t), intent(inout) :: la
      type(rimp2_layout_t), intent(in) :: lay
      integer, intent(in) :: np
      real(dp), intent(in) :: cocc(lay%nao, lay%no), tens(lay%nao, lay%nao, np)
      real(dp), intent(inout) :: half(lay%no, lay%nao, np)

      ! A few library calls per block, not two per auxiliary function:
      ! the occupied side of the whole block is one GEMM,
      !   half(i, nu q) = sum_mu C_occ(mu, i) T(mu, nu q)   (no x nao np),
      ! and the virtual side (rimp2_stage_half_virtual) one strided batch over q with C_vir shared
      ! per destination rank, on that rank's rows of `half`,
      !   X(:, P) = C_vir^T half_q^T                        (nv x bn(q), a fastest).
      ! 2 x 3840 calls at w40 cc-pVDZ ran at 41% of the device's peak.
      call la%gemm("T", "N", lay%no, lay%nao*np, lay%nao, 1.0_dp, cocc, lay%nao, tens, lay%nao, 0.0_dp, half, lay%no)
   end subroutine rimp2_stage_half_occupied

   subroutine rimp2_stage_half_virtual(la, lay, p0, np, cvir, half, xr, ysend)
      !! (ia|P) for the auxiliary block at column `p0`, width `np`: the rows
      !! for this rank's O_rank go into `xr`, in the columns the block keeps,
      !! the rows for every other rank's block into its segment of `ysend`.
      type(trc_linalg_t), intent(inout) :: la
      type(rimp2_layout_t), intent(in) :: lay
      integer, intent(in) :: p0, np
      real(dp), intent(in) :: cvir(lay%nao, lay%nv), half(lay%no, lay%nao, np)
      real(dp), intent(inout) :: xr(*), ysend(*)
      type(trc_gemm_batch_t) :: gemm_shape
      integer :: q
      integer(int64) :: novr8, nq8, nv8

      novr8 = lay%novr
      nv8 = lay%nv
      do q = 0, lay%nranks - 1
         if (lay%bn(q) == 0) cycle
         if (q == lay%rank) then
            gemm_shape = trc_gemm_batch(lay%nv, lay%nb, lay%nao, np, ldb=lay%no, sa=0_int64, sc=novr8, transa="T", transb="T")
            call la%gemm_batched(gemm_shape, cvir, half(lay%lo, 1, 1), xr(1 + int(p0 - 1, int64)*novr8))
         else
            nq8 = int(lay%bn(q), int64)*nv8
            gemm_shape = trc_gemm_batch(lay%nv, lay%bn(q), lay%nao, np, ldb=lay%no, sa=0_int64, sc=nq8, transa="T", transb="T")
            call la%gemm_batched(gemm_shape, cvir, half(lay%blo(q), 1, 1), &
                                 ysend(1 + lay%soff(q) + int(p0 - lay%plo(lay%rank), int64)*nq8))
         end if
      end do
   end subroutine rimp2_stage_half_virtual

   subroutine rimp2_stage_solve_b(la, stab, naux, novr, jm, xr)
      !! B = X L^-T: solve B L^T = X from the right, in place, for the
      !! `novr` rows of `xr` this rank holds. Opens and closes ST_TRSM and
      !! counts its flops.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      integer, intent(in) :: naux, novr
      real(dp), intent(in) :: jm(naux, naux)
      real(dp), intent(inout) :: xr(*)
      type(timer_type) :: tc

      call stab%step_on(la, ST_TRSM, tc)
      if (novr > 0) call la%trsm("R", "L", "T", novr, naux, 1.0_dp, jm, naux, xr, novr)
      call stab%step_off(la, ST_TRSM, tc)
      call stab%add_flop(ST_TRSM, int(novr, int64)*int(naux, int64)*int(naux, int64))
   end subroutine rimp2_stage_solve_b

   subroutine rimp2_stage_pair_layout(lay, ring)
      !! How the pair energies are batched: the ring's length, the number of
      !! occupied orbitals per K GEMM (the device's free memory decides, so
      !! this runs after every array except K's own workspace is allocated)
      !! and the shape of K.
      type(rimp2_layout_t), intent(in) :: lay
      type(rimp2_ring_t), intent(inout) :: ring

      ring%nring = lay%nranks/2
      ring%ncolmax = lay%nmax*lay%nv
      ring%ng = k_group(lay%nmax, lay%nv, ring%ncolmax)
      ring%ldk = ring%ng*lay%nv
   end subroutine rimp2_stage_pair_layout

   subroutine rimp2_ring_step_range(lay, ring, step)
      !! The pairs of ring step `step`: the visiting rank `qv`, the rows i and
      !! columns j this rank computes against its visitor, whether it is the
      !! diagonal block, and the leading dimension of the block the columns
      !! are read from.
      type(rimp2_layout_t), intent(in) :: lay
      type(rimp2_ring_t), intent(inout) :: ring
      integer, intent(in) :: step

      ring%qv = mod(lay%rank + step, lay%nranks)
      ! Rows i in ilo..ihi of this rank's block against columns j in
      ! jlo..jhi of block qv; `diag` keeps j >= i.
      ring%diag = step == 0
      ring%ilo = lay%lo
      ring%ihi = lay%lo + lay%nb - 1
      ring%jbase = lay%blo(ring%qv)
      ring%jlo = ring%jbase
      ring%jhi = ring%jbase + lay%bn(ring%qv) - 1
      if (2*step == lay%nranks) then
         if (lay%rank < ring%qv) then
            ring%ihi = lay%lo + lay%nb/2 - 1
         else
            ring%jlo = ring%jbase + lay%bn(ring%qv)/2
         end if
      end if
      if (step == 0) then
         ring%ldv = lay%novr
      else
         ring%ldv = lay%bn(ring%qv)*lay%nv
      end if
   end subroutine rimp2_ring_step_range

   subroutine rimp2_group_range(lay, ring, i0)
      !! The group of occupied orbitals starting at `i0`: its last orbital,
      !! the first column's orbital and the columns of its K GEMM (zero or
      !! fewer when the group has no pairs at this ring step).
      type(rimp2_layout_t), intent(in) :: lay
      type(rimp2_ring_t), intent(inout) :: ring
      integer, intent(in) :: i0

      ring%i0 = i0
      ring%i1 = min(ring%ihi, i0 + ring%ng - 1)
      ring%jstart = ring%jlo
      if (ring%diag) ring%jstart = i0
      ring%ncol = (ring%jhi - ring%jstart + 1)*lay%nv
   end subroutine rimp2_group_range

   subroutine rimp2_stage_k_gemm_own(la, stab, lay, ring, xr, ki)
      !! K(i a, j b) = sum_Q B(i a, Q) B(j b, Q) for the group in hand, both
      !! blocks this rank's own (`xr`). Opens and closes ST_KGEMM and counts
      !! its flops.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      type(rimp2_layout_t), intent(in) :: lay
      type(rimp2_ring_t), intent(in) :: ring
      real(dp), intent(in) :: xr(*)
      real(dp), intent(inout) :: ki(ring%ldk, *)
      type(timer_type) :: tc

      call stab%step_on(la, ST_KGEMM, tc)
      call la%gemm("N", "T", (ring%i1 - ring%i0 + 1)*lay%nv, ring%ncol, lay%naux, 1.0_dp, xr(1 + (ring%i0 - lay%lo)*lay%nv), &
                   lay%novr, xr(1 + (ring%jstart - ring%jbase)*lay%nv), ring%ldv, 0.0_dp, ki, ring%ldk)
      call stab%step_off(la, ST_KGEMM, tc)
      call stab%add_flop(ST_KGEMM, 2*int(ring%i1 - ring%i0 + 1, int64)*int(lay%nv, int64)*int(ring%ncol, int64)* &
                         int(lay%naux, int64))
   end subroutine rimp2_stage_k_gemm_own

   subroutine rimp2_stage_k_gemm_visiting(la, stab, lay, ring, xr, vis, ki)
      !! The same K GEMM against the visiting block, in buffer `cur` of `vis`.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      type(rimp2_layout_t), intent(in) :: lay
      type(rimp2_ring_t), intent(in) :: ring
      real(dp), intent(in) :: xr(*), vis(int(lay%nmax, int64)*lay%nv*lay%naux, *)
      real(dp), intent(inout) :: ki(ring%ldk, *)
      type(timer_type) :: tc

      call stab%step_on(la, ST_KGEMM, tc)
      call la%gemm("N", "T", (ring%i1 - ring%i0 + 1)*lay%nv, ring%ncol, lay%naux, 1.0_dp, xr(1 + (ring%i0 - lay%lo)*lay%nv), &
                   lay%novr, vis(1 + (ring%jstart - ring%jbase)*lay%nv, ring%cur), ring%ldv, 0.0_dp, ki, ring%ldk)
      call stab%step_off(la, ST_KGEMM, tc)
      call stab%add_flop(ST_KGEMM, 2*int(ring%i1 - ring%i0 + 1, int64)*int(lay%nv, int64)*int(ring%ncol, int64)* &
                         int(lay%naux, int64))
   end subroutine rimp2_stage_k_gemm_visiting

   subroutine rimp2_stage_pair_group(la, stab, ring, ki, kt, part, eps_occ, eps_vir)
      !! The pair energies of the group in hand, one occupied i at a time:
      !! K_i^T, the kernel, the `part` it wrote back to the host, and its sum
      !! into `ring%e_os` and `ring%e_ss`. Opens and closes ST_PAIR per i.
      !! The numbers of occupied and virtual orbitals are the sizes of `eps_occ`
      !! and `eps_vir`.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      type(rimp2_ring_t), intent(inout) :: ring
      real(dp), intent(in) :: ki(ring%ldk, *)
      real(dp), intent(in), contiguous :: eps_occ(:)
      real(dp), intent(in), contiguous :: eps_vir(:)
      real(dp), intent(inout) :: kt(ring%ncolmax, size(eps_vir)), part(2*ring%ncolmax)
      type(timer_type) :: tc
      integer :: no, nv
      integer :: i, roff, coff, jfirst, ncol_i

      no = size(eps_occ)
      nv = size(eps_vir)
      do i = ring%i0, ring%i1
         ! K_i(a, jb) is the block of rows roff + a and, on the
         ! diagonal, columns from j = i; off it, every column.
         roff = (i - ring%i0)*nv
         if (ring%diag) then
            coff = (i - ring%i0)*nv
            jfirst = i
         else
            coff = 0
            jfirst = ring%jstart
         end if
         ncol_i = (ring%jhi - jfirst + 1)*nv
         call stab%step_on(la, ST_PAIR, tc)
         ! K_i^T, so the kernel reads (ia|jb) and (ib|ja) both coalesced
         call la%transpose(nv, ncol_i, ki(roff + 1, coff + 1), ring%ldk, kt, ring%ncolmax)
         call pair_energy(i, jfirst, no, nv, ring%ldk, roff, coff, ring%ncolmax, ncol_i, eps_occ, eps_vir, ki, kt, part)
         !$acc update self(part(1:2*ncol_i))
         call stab%step_off(la, ST_PAIR, tc)
         call stab%add_flop(ST_PAIR, 10*int(nv, int64)*int(ncol_i, int64))
         ring%e_os = ring%e_os + sum(part(1:2*ncol_i:2))
         ring%e_ss = ring%e_ss + sum(part(2:2*ncol_i:2))
      end do
   end subroutine rimp2_stage_pair_group

   subroutine rimp2_stage_reduce_ranks(comm, stab, lay, ring)
      !! The pair energies and the step table over the ranks. A one-rank run
      !! does nothing.
      type(comm_t), intent(in) :: comm
      type(step_table_t), intent(inout) :: stab
      type(rimp2_layout_t), intent(in) :: lay
      type(rimp2_ring_t), intent(inout) :: ring
      real(dp) :: es(2)

      if (lay%nranks > 1) then
         es = [ring%e_os, ring%e_ss]
         call allreduce(comm, es, op=MPI_SUM)
         call bcast(comm, es, 2, 0)
         ring%e_os = es(1)
         ring%e_ss = es(2)
      end if
      ! The job's table, not rank 0's share: the slowest rank's time and the
      ! work of every rank. The energy reduction has already lined the
      ! ranks up, so rank 0's wall time in `rimp2_log_report` is the job's.
      if (lay%nranks > 1) call stab%reduce(comm)
   end subroutine rimp2_stage_reduce_ranks

   subroutine urimp2_stage_dimensions(b, aux, nalpha, nbeta, op, res, lay, ok)
      !! The sizes of an unrestricted run and the checks that decide whether
      !! it can run at all. `ok` is false, with `res%message` set, when more
      !! than one rank is asked for, a spin channel has nothing to correlate
      !! or a basis carries angular momentum above what the density-fitting
      !! kernels were generated for.
      type(trc_basis_t), intent(in) :: b, aux
      integer, intent(in) :: nalpha, nbeta
      type(trc_rimp2_options_t), intent(in) :: op
      type(trc_rimp2_result_t), intent(inout) :: res
      type(urimp2_layout_t), intent(inout) :: lay
      logical, intent(out) :: ok
      integer :: nranks

      ok = .false.
      lay%nao = b%nao
      ! `naux` is the folded count when `aux%spherical`; the kernels and the
      ! block walk use the CARTESIAN count, `aux%nao`.
      lay%naux = merge(aux%nao_sph, aux%nao, aux%spherical)
      nranks = 1
      if (op%has_comm) nranks = op%comm%size()
      if (nranks > 1) then
         res%message = "trc_urimp2: multi-rank UHF RI-MP2 is not implemented this round -- "// &
                       "the restricted ring distributes one block of B and lets it meet "// &
                       "itself; an open-shell run needs that for aa AND bb (different sizes "// &
                       "whenever nalpha /= nbeta) plus a third, non-ring distribution for the "// &
                       "mixed ab block. Run on one rank (no communicator, or a single-rank one)"
         return
      end if
      lay%nfz = max(op%nfrozen, 0)
      lay%noa = nalpha - lay%nfz
      lay%nob = nbeta - lay%nfz
      lay%nva = lay%nao - nalpha
      lay%nvb = lay%nao - nbeta
      res%naux = lay%naux
      res%nocc = lay%noa + lay%nob    ! spin orbitals, not spatial -- see
      res%nvir = lay%nva + lay%nvb    ! `trc_urimp2_run`'s own header and fill_ump2_result's
      ! equivalent in metalquicha
      if (lay%noa < 1 .or. lay%nob < 1 .or. lay%nva < 1 .or. lay%nvb < 1) then
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
      lay%nblock = aux%nao
      if (op%aux_block > 0) lay%nblock = op%aux_block
      ok = .true.
   end subroutine urimp2_stage_dimensions

   subroutine urimp2_stage_spin_transform(la, stab, p0, cocc, cvir, tens, xr)
      !! (ia|P) of one spin for the auxiliary block in `tens`: `xr` gets the
      !! block's columns, from the column `p0`. Opens and closes ST_TRANS and
      !! counts its flops; the half-transformed `half` lives here.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      integer, intent(in) :: p0
      real(dp), intent(in), contiguous :: cocc(:, :), cvir(:, :), tens(:, :, :)
      real(dp), intent(inout) :: xr(*)
      real(dp), allocatable :: half(:, :, :)
      type(trc_gemm_batch_t) :: gemm_shape
      type(timer_type) :: tc
      integer :: nao, no, nv, np
      integer(int64) :: nao8, no8, nv8, novr8

      nao = size(cocc, 1)
      no = size(cocc, 2)
      nv = size(cvir, 2)
      np = size(tens, 3)
      nao8 = nao
      no8 = no
      nv8 = nv
      novr8 = no8*nv8
      allocate (half(no, nao, np))
      !$acc enter data create(half)
      call stab%step_on(la, ST_TRANS, tc)
      call la%gemm("T", "N", no, nao*np, nao, 1.0_dp, cocc, nao, tens, nao, 0.0_dp, half, no)
      gemm_shape = trc_gemm_batch(nv, no, nao, np, sa=0_int64, sc=novr8, transa="T", transb="T")
      call la%gemm_batched(gemm_shape, cvir, half, xr(1 + int(p0 - 1, int64)*novr8))
      call stab%step_off(la, ST_TRANS, tc)
      call stab%add_flop(ST_TRANS, int(np, int64)*(2*no8*nao8*nao8 + 2*nv8*no8*nao8))
      !$acc exit data delete(half)
      deallocate (half)
   end subroutine urimp2_stage_spin_transform

   subroutine urimp2_stage_same_spin(la, stab, naux, xr, eps_occ, eps_vir, e)
      !! The like-spin energy of one spin from its B (`xr`, `size(eps_occ)` x
      !! `size(eps_vir)` rows): restricted's own K GEMM and `pair_energy` kernel
      !! over the occupied-pair triangle (j >= i, weight 2 off the diagonal, 1 on
      !! it), keeping the same-spin half the kernel returns and the 1/2
      !! `trc_urimp2_run`'s header derives; the opposite-spin half it also
      !! computes has no meaning for two electrons of the same spin and is
      !! discarded. Opens and closes ST_KGEMM and ST_PAIR; `eps_occ` and `eps_vir` are
      !! on the device already.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      integer, intent(in) :: naux
      real(dp), intent(in) :: xr(*)
      real(dp), intent(in), contiguous :: eps_occ(:)
      real(dp), intent(in), contiguous :: eps_vir(:)
      real(dp), intent(out) :: e
      real(dp), allocatable :: kgemm(:, :), kt(:, :), part(:)
      type(timer_type) :: tc
      integer :: no, nv, novr, ncolmax, ng, ldk, i0, i1, i, jstart, jfirst, ncol, ncol_i, roff, coff
      integer(int64) :: nv8, naux8

      no = size(eps_occ)
      nv = size(eps_vir)
      nv8 = nv
      naux8 = naux
      novr = no*nv
      ncolmax = no*nv
      ng = k_group(no, nv, ncolmax)
      ldk = ng*nv
      allocate (kgemm(ldk, ncolmax), kt(ncolmax, nv), part(2*ncolmax))
      !$acc enter data create(kgemm, kt, part)
      e = 0.0_dp
      do i0 = 1, no, ng
         i1 = min(no, i0 + ng - 1)
         jstart = i0
         ncol = (no - jstart + 1)*nv
         call stab%step_on(la, ST_KGEMM, tc)
         call la%gemm("N", "T", (i1 - i0 + 1)*nv, ncol, naux, 1.0_dp, xr(1 + (i0 - 1)*nv), novr, &
                      xr(1 + (jstart - 1)*nv), novr, 0.0_dp, kgemm, ldk)
         call stab%step_off(la, ST_KGEMM, tc)
         call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nv8*int(ncol, int64)*naux8)
         do i = i0, i1
            roff = (i - i0)*nv
            coff = (i - i0)*nv
            jfirst = i
            ncol_i = (no - jfirst + 1)*nv
            call stab%step_on(la, ST_PAIR, tc)
            call la%transpose(nv, ncol_i, kgemm(roff + 1, coff + 1), ldk, kt, ncolmax)
            call pair_energy(i, jfirst, no, nv, ldk, roff, coff, ncolmax, ncol_i, eps_occ, eps_vir, kgemm, kt, part)
            !$acc update self(part(1:2*ncol_i))
            call stab%step_off(la, ST_PAIR, tc)
            call stab%add_flop(ST_PAIR, 10*nv8*int(ncol_i, int64))
            e = e + sum(part(2:2*ncol_i:2))   ! the same-spin half only
         end do
      end do
      !$acc exit data delete(kgemm, kt, part)
      deallocate (kgemm, kt, part)
      e = 0.5_dp*e
   end subroutine urimp2_stage_same_spin

   subroutine urimp2_stage_mixed_gemm(la, stab, lay, i0, i1, xr_a, xr_b, ki)
      !! K(i A, J B) = sum_Q B_a(i A, Q) B_b(J B, Q) for the alpha occupied
      !! orbitals `i0` to `i1` against every beta orbital J, into `ki`, whose
      !! leading dimension is the group size times nva. Opens and closes
      !! ST_KGEMM and counts its flops.
      type(trc_linalg_t), intent(inout) :: la
      type(step_table_t), intent(inout) :: stab
      type(urimp2_layout_t), intent(in) :: lay
      integer, intent(in) :: i0, i1
      real(dp), intent(in) :: xr_a(*), xr_b(*)
      real(dp), intent(inout), contiguous :: ki(:, :)
      type(timer_type) :: tc
      integer :: ncol
      integer(int64) :: nva8, naux8

      nva8 = lay%nva
      naux8 = lay%naux
      ncol = lay%nob*lay%nvb
      call stab%step_on(la, ST_KGEMM, tc)
      call la%gemm("N", "T", (i1 - i0 + 1)*lay%nva, ncol, lay%naux, 1.0_dp, xr_a(1 + (i0 - 1)*lay%nva), lay%noa*lay%nva, &
                   xr_b, lay%nob*lay%nvb, 0.0_dp, ki, size(ki, 1))
      call stab%step_off(la, ST_KGEMM, tc)
      call stab%add_flop(ST_KGEMM, 2*int(i1 - i0 + 1, int64)*nva8*int(ncol, int64)*naux8)
   end subroutine urimp2_stage_mixed_gemm

   subroutine rimp2_log_report(tag, la, stab, tw0, tw1, rank, nranks)
      !! The closing lines: the pair-energy time and the run's total, then the
      !! step table with the dense-algebra subtotal. Rank 0 only, so the
      !! caller asks for it on `talk`. `tag` is the log prefix ("rimp2" or "rimp2u").
      character(len=*), intent(in) :: tag
      type(trc_linalg_t), intent(in) :: la
      type(step_table_t), intent(inout) :: stab
      type(timer_type), intent(in) :: tw0, tw1
      integer, intent(in) :: rank, nranks
      logical :: dense_mask(NSTEP)
      character(len=1024) :: logline

      write (logline, "(a,a,a,f8.2,a,f8.2,a)") "  ", tag, ": pair energies ", since(la, tw1), " s;  total ", since(la, tw0), " s"
      call global_logger%performance(trim(logline))
      dense_mask = .false.
      dense_mask(ST_CHOL) = .true.
      dense_mask(ST_TRANS) = .true.
      dense_mask(ST_TRSM) = .true.
      dense_mask(ST_KGEMM) = .true.
      call stab%report(since(la, tw0), rank, nranks, subtotal_mask=dense_mask, subtotal_label="all dense algebra")
   end subroutine rimp2_log_report

   subroutine rimp2_log_metric(tag, naux, t)
      !! "metric (naux) and Cholesky t s", prefixed by `tag`.
      character(len=*), intent(in) :: tag
      integer, intent(in) :: naux
      real(dp), intent(in) :: t
      character(len=1024) :: logline

      write (logline, "(a,a,a,i0,a,f8.2,a)") "  ", tag, ": metric (", naux, ") and Cholesky ", t, " s"
      call global_logger%performance(trim(logline))
   end subroutine rimp2_log_metric

   subroutine rimp2_log_transform(tag, stab, exchange)
      !! The integral and transform times of the step table so far, with the
      !! exchange's when `exchange` (a run over ranks has one).
      character(len=*), intent(in) :: tag
      type(step_table_t), intent(in) :: stab
      logical, intent(in) :: exchange
      character(len=1024) :: logline

      if (exchange) then
         write (logline, "(a,a,a,f8.2,a,f8.2,a,f8.2,a)") "  ", tag, ": (mu nu|P) ", stab%t(ST_3C), " s, transform to (ia|P) ", &
            stab%t(ST_TRANS), " s, exchange ", stab%t(ST_XCHG), " s"
      else
         write (logline, "(a,a,a,f8.2,a,f8.2,a)") "  ", tag, ": (mu nu|P) ", stab%t(ST_3C), " s, transform to (ia|P) ", &
            stab%t(ST_TRANS), " s"
      end if
      call global_logger%performance(trim(logline))
   end subroutine rimp2_log_transform

   subroutine rimp2_log_solve(tag, stab)
      !! The triangular solves' time of the step table.
      character(len=*), intent(in) :: tag
      type(step_table_t), intent(in) :: stab
      character(len=1024) :: logline

      write (logline, "(a,a,a,f8.2,a)") "  ", tag, ": B = X L^-T ", stab%t(ST_TRSM), " s"
      call global_logger%performance(trim(logline))
   end subroutine rimp2_log_solve

   subroutine rimp2_log_exchange_mode(direct)
      !! How the ring moves B: device to device, or through the host. Only
      !! the restricted run calls it, so the tag is fixed.
      logical, intent(in) :: direct
      character(len=1024) :: logline

      write (logline, "(a,a)") "  rimp2: B exchange ", &
         trim(merge("device to device (CUDA-aware MPI)", "through the host                 ", direct))
      call global_logger%performance(trim(logline))
   end subroutine rimp2_log_exchange_mode

   subroutine rimp2_log_k_groups(ng, nranks)
      !! The occupied orbitals per K GEMM and the ranks of the run. Only the
      !! restricted run calls it, so the tag is fixed.
      integer, intent(in) :: ng, nranks
      character(len=1024) :: logline

      write (logline, "(a,i0,a,i0,a)") "  rimp2: K GEMMs over ", ng, " occupied orbitals at a time, ", &
         nranks, " rank(s)"
      call global_logger%performance(trim(logline))
   end subroutine rimp2_log_k_groups

   subroutine rimp2_store_result(res, stab, t_total, e_os, e_ss)
      !! The pair energies and, for `--json` and the regression bench, the
      !! step table's times into the result (each step waits for the device
      !! before its clock starts and stops).
      type(trc_rimp2_result_t), intent(inout) :: res
      type(step_table_t), intent(in) :: stab
      real(dp), intent(in) :: t_total, e_os, e_ss
      res%e_os = e_os
      res%e_ss = e_ss
      res%e_corr = e_os + e_ss
      res%t_total = t_total
      res%step_name = stab%names
      res%step_t = stab%t
   end subroutine rimp2_store_result

   function since(lp, tm) result(t)
      !!
      !! The time `tm` has run, read after everything queued on the device has
      !! finished.
      !!
      type(trc_linalg_t), intent(in) :: lp
      type(timer_type), intent(in) :: tm
      real(dp) :: t
      call lp%sync()
      t = tm%get_elapsed_time()
   end function since

   subroutine ring_shift(comm, xr, vis, cur, nsend, nrecv, first, direct)
      !!
      !! One step of the ring: pass the block this rank holds to the rank on its
      !! left and take the one on its right, from device buffers when
      !! `direct` and through the host otherwise. At the first step the held
      !! block is this rank's own (`xr`), after that the visiting buffer
      !! `cur`, and what arrives goes into the other buffer, which becomes
      !! `cur`. Every rank posts its send non-blocking, then receives, then
      !! waits for the send, so no order of ranks can deadlock. In pieces of
      !! at most MSG_CHUNK doubles, since MPI counts are 32-bit.
      !!
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: xr(:), vis(:, :)
      integer, intent(inout) :: cur
      integer(int64), intent(in) :: nsend, nrecv
      logical, intent(in) :: first, direct
      integer :: rank, nranks, left, right, nxt, nreq
      type(request_t), allocatable :: reqs(:)
      rank = comm%rank()
      nranks = comm%size()
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

   subroutine send_block(comm, xr, vis, cur, n, dest, first, reqs)
      !!
      !! The held block to `dest`, from the host copy, in MSG_CHUNK pieces.
      !! Module procedures rather than internal ones: see app/main.F90's
      !! trc_run for what a sibling internal procedure did under NVHPC's MPI.
      !!
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

   subroutine recv_block(comm, vis, c, n, source)
      !!
      !! The next block from `source` into buffer `c` of `vis`, in MSG_CHUNK pieces.
      !!
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

   subroutine send_block_dev(comm, xr, vis, cur, n, dest, first, reqs)
      !!
      !! send_block and recv_block for device addresses, called inside
      !! `host_data use_device`: pic-mpi's explicit-count procedures by their
      !! specific names, since its generic interfaces misbehave there under
      !! nvfortran.
      !!
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

   subroutine exchange_x(comm, nv8, bn, plo, soff, ysend, xr, direct)
      !!
      !! Step 3's exchange. Rank r holds X(ia, P in P_r) for every i, with the
      !! rows for O_q (q /= r) packed in `ysend` from soff(q). It needs
      !! X(ia in O_r, P) for every P. The exchange runs nranks - 1 shifts: at
      !! step k rank r sends its segment to rank r + k and receives from rank
      !! r - k the columns P_{r-k} of its own X, which are contiguous in `xr`.
      !! As in ring_shift, each send is posted non-blocking before the blocking
      !! receive, so neither side can wait on the other. Messages go in
      !! MSG_CHUNK pieces. Sender and receiver compute the same count,
      !! bn(dest) nv |P_r|, so an empty one is skipped on both sides.
      !!
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: bn(0:), plo(0:)
      integer(int64), intent(in) :: nv8
      integer(int64), intent(in) :: soff(0:)
      real(dp), intent(inout) :: ysend(:), xr(:)
      logical, intent(in) :: direct
      integer :: rank, nranks, k, dest, src, nreq
      integer(int64) :: nsend, nrecv, roff, novr8
      type(request_t), allocatable :: reqs(:)
      rank = comm%rank()
      nranks = comm%size()
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

   subroutine send_seg(comm, buffer, off, n, dest, reqs)
      !!
      !! buf(off + 1 : off + n) to `dest` and from `source`, in MSG_CHUNK pieces,
      !! host addresses (send_seg, recv_seg) or device ones inside `host_data`
      !! (the _dev pair, by pic-mpi's specific names as send_block_dev does).
      !!
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
      !! Occupied orbitals per K GEMM. One i gives cuBLAS nvir rows, which it
      !! tiles small; four give it a matrix tall enough for its larger tiles.
      !! Measured at w40 cc-pVDZ on a V100, K GEMM seconds for 1, 2, 4, 8, 16
      !! per group: 9.37, 9.09, 8.78, 9.24, 9.38 -- four is the best, and
      !! larger groups lose to the j < i corner each one computes and throws
      !! away. K is that many times larger, so on the device the group shrinks
      !! to what the free memory holds next to K^T and the energies, keeping
      !! 512 MB back for the libraries' own workspace.
      !!
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

   pure integer function rimp2_pack_index(mu, nu) result(p)
      !! idx(mu, nu) = mu(mu-1)/2 + nu, mu >= nu, 1-based: trc_scf_df.F90's own
      !! `pack_idx`, duplicated here rather than used across the module
      !! boundary -- a `do concurrent` body may not call a procedure from
      !! another module.
      !$acc routine seq
      integer, intent(in) :: mu, nu
      p = mu*(mu - 1)/2 + nu
   end function rimp2_pack_index

   subroutine pair_energy(i, jfirst, no, nv, ldk, roff, coff, ldt, nji, eo, ev, ki, kt, part)
      !!
      !! For occupied i and the pairs j = jfirst.., one thread per (j, a), with a
      !! fastest and b summed inside: part(2 ja - 1) = w sum_b (ia|jb)^2 / D and
      !! part(2 ja) = w sum_b [(ia|jb)^2 - (ia|jb)(ib|ja)] / D, w = 1 for j = i
      !! and 2 otherwise -- the (j, i) pair is the same numbers again and is
      !! never computed. No reduction crosses threads; the host adds 2 nji
      !! numbers. Diagonal terms are in the sum: (ia|ia)^2 - (ia|ia)^2 = 0 for
      !! the same-spin part, as it must.
      !!
      !! (ia|jb) is K(roff + a, coff + (j - jfirst) nv + b) and (ib|ja) is
      !! K(roff + b, coff + (j - jfirst) nv + a) = K_i^T(ja, b): with a fastest
      !! across threads, both are unit-stride reads. From K alone one of the
      !! two is always strided by nv across a warp, which held this kernel to a
      !! fifth of the device's bandwidth.
      !!
      integer, intent(in) :: i, jfirst, no, nv, ldk, roff, coff, ldt, nji
      real(dp), intent(in) :: eo(no), ev(nv), ki(ldk, *), kt(ldt, nv)
      real(dp), intent(out) :: part(*)
      integer :: ja
      do concurrent(ja=1:nji)
         call pair_energy_body(i, jfirst, TRC_DC_INDEX(ja), no, nv, ldk, roff, coff, ldt, eo, ev, ki, kt, &
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

   subroutine mixed_energy(i, nva, nob, nvb, ldk, roff, eo_a, ev_a, eo_b, ev_b, ki, part)
      !!
      !! The mixed-spin (ab) block for occupied alpha orbital `i`: one thread
      !! per (a, J), B summed inside. `aJ` decomposes the same way `pair_energy`'s
      !! `ja` does (the row selector `a` fastest, the column group `J` next),
      !! here against a rectangle that is never folded by any symmetry, so
      !! there is no weight, no second read of a transpose, and nothing to
      !! subtract -- only (ia|JB)^2, never (ia|JB)(iB|Ja), which does not exist
      !! for two distinguishable electrons.
      !!
      integer, intent(in) :: i, nva, nob, nvb, ldk, roff
      real(dp), intent(in) :: eo_a(*), ev_a(nva), eo_b(nob), ev_b(nvb), ki(ldk, *)
      real(dp), intent(out) :: part(*)
      integer :: aj
      do concurrent(aj=1:nva*nob)
         call mixed_energy_body(i, TRC_DC_INDEX(aj), nva, nob, nvb, roff, eo_a, ev_a, eo_b, ev_b, ldk, ki, part(aj))
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
