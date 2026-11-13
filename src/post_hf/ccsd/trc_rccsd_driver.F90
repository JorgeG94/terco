module trc_rccsd_driver
   !!
   !! The closed-shell RI-CCSD driver: integrals, MP2 start, amplitude update with DIIS,
   !! convergence, and the converged amplitudes kept for the (T) step.
   !!
   !! ALGORITHM, per iteration (mqc's `run_czt_rccsd`, with PySCF's convergence test):
   !!
   !!     (t1', t2') = update(t1, t2)              trc_rccsd_update, from the SAME input amplitudes
   !!     e = t' - t                               the step; its size is measured on the device
   !!     (t1, t2) <- DIIS(history of t', e)       or t' itself while the history is short (trc_rccsd_diis)
   !!     E = E_corr(t1, t2)                       of the amplitudes the NEXT iteration starts from
   !!
   !! the start is t1 = 0, t2 = MP2 (E_MP2 is its energy), or the amplitudes the caller left in
   !! `state%amps` when `state%use_start` is set. DIIS starts pushing at iteration `diis_start` (default 1, as mqc)
   !! with the history length `diis_size` (default 8; 0 or 1 switches it off).
   !!
   !! CONVERGENCE. The run stops as CONVERGED when BOTH
   !!
   !!     |E_k - E_{k-1}| < conv_energy      (default 1e-8; E_0 is the energy of the starting amplitudes, as PySCF's e_prev)
   !!     max |t' - t| < conv_amplitude      (default 1e-6; the largest element of the step over t1 and t2)
   !!
   !! hold. mqc tests the energy alone, and only from the second iteration; the amplitude criterion is a deliberate
   !! addition. PySCF tests both, with conv_tol 1e-7 and conv_tol_normt 1e-5. Its normt is
   !! ||vector(t1', t2') - vector(t1, t2)||_2 where the vector is t1 followed by the PACKED lower triangle of t2 as an
   !! (o v x o v) matrix: each off-diagonal pair (i,a),(j,b) once, no sqrt(2), N' = o v + o v (o v + 1)/2 elements. The
   !! `result%step_norm2` here is the norm over the FULL t2 (a pair-symmetric step counts every off-diagonal pair
   !! twice), so normt <= step_norm2 <= sqrt(2) normt: it is NOT PySCF's number, and comparing it to 1e-5 is not PySCF's
   !! test. Nor are the two tests nested: PySCF's normt < 1e-5 forces max|e| < 1e-5 only, while the default here
   !! (max|e| < 1e-6) allows normt up to sqrt(N') 1e-6, which exceeds 1e-5 once N' > 100 (every case of interest). It is
   !! the stricter test per element and can be the looser one in the 2-norm, by up to sqrt(N')/10; the energy test (1e-8
   !! against PySCF's 1e-7) is the tighter of the two.
   !! A step or an energy that is NaN or infinite stops the run as NOT converged with a message, whatever
   !! the comparisons would say; the test is written so a NaN can never satisfy it.
   !!
   !! THE T1 DIAGNOSTIC is PySCF's `get_t1_diagnostic`: sqrt(||t1||_F^2 / (2 no)), no the number of correlated
   !! occupied orbitals (2 no correlated electrons), of the converged (or last) t1.
   !!
   !! THE HANDLE. `trc_rccsd_solve` leaves in `state` what the (T) step needs without recomputing: `state%ints` (the
   !! fitted blocks, on the host and the device, stamp unchanged) and `state%amps` (t1, t2 on both, the last iterate,
   !! the converged one when `state%result%converged`). The caller owns them and calls `state%release()`. A refused
   !! run leaves `state%ints` empty (its `status` says why) and, unless `use_start` supplied them, `state%amps` too; a run that
   !! was not converged leaves both in place. The update
   !! object, the DIIS history and the second set of amplitudes are private to the run and always released before it
   !! returns, on every path. `trc_rccsd_run` is the same run for a caller that only wants the numbers: it releases
   !! everything.
   !!
   !! TIMING. Each stage is a step of a `step_table_t` the run owns (trc_rccsd_steps names them), timed on the device
   !! (`la%sync()` before the clock is read, since cuBLAS returns before its kernels end) with an NVTX range, and
   !! carrying the flops of its GEMMs where it has any. The times are summed over the iterations into `result%seconds`
   !! (the ladder and its dressing, and the T1 and T2 equations with the rings, each in one entry as before); `verbose`
   !! prints one line per iteration, then the step table of RI-MP2's shape and the mean time of an iteration, through
   !! global_logger.
   !!
   !! The arguments are those of `trc_rimp2_run` except that the doubly occupied count is `opts%nocc` (more than eight
   !! dummies is a rule of this tree) and the options carry what trc_rccsd_ints refuses on: an ECP, an open shell
   !! (`nbeta`). Frozen core: the caller passes the number (`opts%nfrozen`); the count-by-element helper lives in the API
   !! layer, which this one may not use.
   !!
   !! RANKS. With `opts%has_comm` and more than one rank in `opts%comm` (one rank per device, none shared), the run is
   !! collective: every rank calls it with the same arguments, and no MPI call is made otherwise.
   !!   - ORBITALS. In the CLI and the C API the SCF is collective (rank 0 diagonalises and broadcasts `cmo` and `eps`), so
   !!     every rank already holds the same orbitals. The run still broadcasts rank 0's `cmo` and `eps` first (one nmo^2
   !!     message) as a GUARD for a caller that hands ranks different ones: SCFs run separately on each rank agree to the
   !!     threshold but can differ in the sign of a column, and rank 0's amplitudes in another rank's integrals gave a (T)
   !!     share wrong by 5e-5 Eh on water.
   !!   - INTEGRALS are built on EVERY rank (replicated); the (T) step of trc_rccsd_triples distributes the (T) work.
   !!   - THE ITERATIONS. Rank 0 runs the iteration (the intermediates, W, the T1 and T2 equations with the rings, DIIS, the
   !!     energy and the convergence test); the particle-particle ladder, the largest term, is cut over ALL the ranks by
   !!     the virtual index d (trc_rccsd_ladder, SHARES; the cut is trc_rccsd_share's, smaller for rank 0 because of the
   !!     other work it has). Per iteration, each a collective call in a statement of its own on every rank:
   !!       1. the control word (trc_rccsd_mpi, rccsd_mpi_control): 0 to stop, else 2 it + (1 when tau has the pair symmetry
   !!          the half ladder's one-contraction path needs). Rank 0 sends 0 whenever it leaves the loop for a reason of its
   !!          own (converged, the limit, a NaN, a refusal after the ladder's sum), so no follower is left waiting;
   !!       2. the amplitudes, t1 and t2, device to device (trc_rccsd_xfer): point to point from rank 0 to each other rank,
   !!          the unique half of t2 when tau has the pair symmetry (every production run), through MPI device buffers
   !!          when trc_mpi_device_direct() says the MPI takes them and through the host otherwise (step RCCSD_ST_BCAST).
   !!          Not sent with the stop word;
   !!       3. every rank builds tau and its dressed factor bt and runs the passes of its range d_lo..d_hi on its device.
   !!          Rank 0 does its own range and then EVERYTHING that does not read the ladder (the intermediates and W first,
   !!          then the T1 equation, every T2 term but the ladder, the rings and the t1 t1 correction of its own share):
   !!          the ladder is added to t2 last, after the rings and before the denominators, so rank 0's work overlaps the others'
   !!          passes (trc_rccsd_update, `update_step`); an empty range is a rank that runs no pass;
   !!       4. ONE agreement on failure (rccsd_mpi_agree_text; step RCCSD_ST_WAIT, which is also the barrier that makes the
   !!          idle seconds of a rank that finishes early visible): a rank that refused anywhere above gives its reason to
   !!          all, every rank leaves the iteration and the loop, and no control word is due;
   !!       5. each follower with a nonempty range sends its partial ladder once to rank 0 (the unique half when pair
   !!          symmetric), which receives it into a staging array and adds it on its device with a kernel; nothing comes
   !!          back (step RCCSD_ST_REDUCE). The symmetrising kernel at the end of each rank's passes is applied BEFORE the
   !!          transfer: it is linear, so the sum is the same, and a one-rank run (no transfer) is bitwise as it was;
   !!       6. rank 0 adds the ladder to t2, applies the denominators and finishes the iteration. DIIS and its history
   !!          live on rank 0.
   !!     Transfers per iteration: the amplitudes down to each follower and each follower's ladder up, o^2 v^2 doubles (half
   !!     of it when pair symmetric) each. Only the half form can be cut: AUTO becomes the half form over several ranks, and
   !!     an explicit assembled or direct form is refused with a message (all ranks), not run on rank 0 alone. Energies
   !!     agree with a one-rank run to round-off, not bitwise (the ladder is summed in another order, and the ladder is
   !!     added to t2 after the rings in every run, one rank included); the iteration count is the same.
   !!   - FAILURE is agreed, never inferred. A rank that could not build its blocks, or set up its part of the update (or
   !!     got an explicit form that cannot be cut), sets a flag; one agreement (rccsd_mpi_agree_text) after the integrals and
   !!     one after the setup make every rank refused with the message of the lowest-numbered refusing rank, blocks
   !!     released, and nothing further is collective, before any rank enters a loop that waits for another. A failure
   !!     inside an iteration is agreed at step 4 above.
   !!   - OUTCOME. After the loop one broadcast carries rank 0's whole `state%result` (scalars, timings and message) as a
   !!     packet of doubles, then t1 and t2 through the host (trc_rccsd_mpi) into every rank's host arrays, which each rank
   !!     copies to its device as `amps%init` mapped them (the followers' copies were one iteration stale): afterwards
   !!     every rank holds the same result and the same amplitudes, converged or not. Each rank's ladder seconds are
   !!     gathered then (`state%ladder_seconds`). `verbose` output comes from rank 0 only: the step table (it labels itself
   !!     "rank 0 of N", with the three rows of the distributed ladder) and one line per rank with its range and ladder
   !!     seconds.
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_step_table, only: step_table_t
   use trc_rccsd_steps, only: RCCSD_ST_INTEGRALS, RCCSD_ST_METRIC, RCCSD_ST_BLOCKS, RCCSD_ST_INTER, RCCSD_ST_W, &
                              RCCSD_ST_LADDER, RCCSD_ST_DRESSING, RCCSD_ST_T1, RCCSD_ST_T2, RCCSD_ST_RING, RCCSD_ST_DIIS, &
                              RCCSD_ST_ENERGY, RCCSD_STEP_NAME, RCCSD_MULTI_STEP_NAME, RCCSD_ST_BCAST, RCCSD_ST_WAIT, &
                              RCCSD_ST_REDUCE
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use pic_mpi_lib, only: comm_t, bcast
   use trc_rccsd_mpi, only: rccsd_mpi_agree_text, rccsd_mpi_bcast_amps, rccsd_mpi_bcast_orbitals, rccsd_mpi_control, &
                            rccsd_mpi_gather_values
   use trc_rccsd_dist, only: rccsd_dist_start, rccsd_dist_follow, rccsd_dist_word
   use trc_rccsd_share, only: rccsd_share_ranges, rccsd_share_ladder_flops
   use trc_linalg, only: trc_linalg_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, rccsd_inter_energy, rccsd_inter_mp2_guess
   use trc_rccsd_ladder, only: trc_rccsd_ladder_options_t
   use trc_rccsd_update, only: trc_rccsd_update_t
   use trc_rccsd_ladder, only: RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_HALF, RCCSD_LADDER_AUTO
   use trc_rccsd_diis, only: trc_rccsd_diis_t, rccsd_diis_history_bytes, rccsd_diis_free_device_bytes
   implicit none
   private

   integer, parameter :: MESSAGE_LEN = 300  !! length of `trc_rccsd_result_t%message`, and of the text a rank broadcasts
   real(dp), parameter :: DIIS_FREE_FRACTION = 0.7_dp
   !! the DIIS history stays on the device when it takes at most this share of the free device memory at its start;
   !! the rest is for what the iterations allocate after it. GH200, (H2O)9 cc-pVTZ: a history of 18.9 GB with
   !! 41.0 GB free, and the run peaked 3.3 GB above the two (17% of the history, 8% of the free memory); (H2O)10, a
   !! history of 28.8 GB with 14.0 GB free, goes to the host at 9.02 s an iteration and fits
   integer, parameter, public :: RCCSD_TIME_INTEGRALS = 1  !! `result%seconds`: the fitted integral blocks
   integer, parameter, public :: RCCSD_TIME_INTER = 2      !! tau and the one-particle intermediates
   integer, parameter, public :: RCCSD_TIME_W = 3          !! Woooo, Wvoov, Wvovo
   integer, parameter, public :: RCCSD_TIME_LADDER = 4     !! the ladder and its dressing
   integer, parameter, public :: RCCSD_TIME_AMPS = 5       !! the rest of the T1 and T2 equations
   integer, parameter, public :: RCCSD_TIME_DIIS = 6       !! the step measure, the history and the extrapolation
   integer, parameter, public :: RCCSD_TIME_TOTAL = 7      !! the whole run
   integer, parameter, public :: RCCSD_NTIME = 7

   public :: trc_rccsd_options_t, trc_rccsd_result_t, trc_rccsd_state_t
   public :: trc_rccsd_run, trc_rccsd_solve

   type :: trc_rccsd_options_t
      !! Plain values only
      integer :: nocc = 0                 !! doubly occupied orbitals (alpha electrons)
      integer :: nbeta = -1                !! beta electrons; -1 means nocc; anything else is an open shell and is refused
      integer :: nfrozen = 0              !! the lowest nfrozen occupied orbitals are left out of the correlation
      logical :: has_ecp = .false.        !! the calculation carries an effective core potential: refused
      real(dp) :: conv_energy = 1.0e-8_dp     !! |dE| below this, and the next, converge the run
      real(dp) :: conv_amplitude = 1.0e-6_dp  !! max |t' - t| below this
      integer :: max_iterations = 100
      integer :: diis_size = 8            !! history length; 0 or 1 is no DIIS
      integer :: diis_start = 1           !! first iteration whose step enters the history
      logical :: diis_host = .false.      !! force the DIIS history onto the host; else it goes there only when it does
      !! not fit in DIIS_FREE_FRACTION of the free device memory (trc_rccsd_diis, HOST HISTORY)
      integer :: aux_block = 0            !! auxiliary functions per (mu nu|P) block; 0 means all at once
      type(trc_rccsd_ladder_options_t) :: ladder  !! form and batch sizes of the particle ladder, passed through
      logical :: verbose = .false.        !! one line per iteration and the timing table (rank 0 only)
      type(comm_t) :: comm
      !! The caller's communicator, a copy of its handle; read only when `has_comm`. Without one, one rank.
      logical :: has_comm = .false.       !! the run is collective over `comm` (module header, RANKS)
   end type trc_rccsd_options_t

   type :: trc_rccsd_result_t
      !! What a run reports. Energies are correlation energies in Hartree.
      logical :: converged = .false.
      logical :: refused = .false.       !! the run did not start: `message` says why
      logical :: finite = .true.         !! false when a NaN or an infinity stopped the run
      integer :: iterations = 0
      integer :: no = 0                  !! active occupied orbitals
      integer :: nv = 0                  !! virtual orbitals
      integer :: naux = 0                !! auxiliary functions
      integer :: nfrozen = 0
      integer :: diis_restarts = 0       !! DIIS breakdowns (singular or non-finite B)
      real(dp) :: e_mp2 = 0.0_dp         !! the MP2 correlation energy of the reference
      real(dp) :: e_ccsd = 0.0_dp
      !! e_singles + e_doubles at the last iterate (the last finite one if a NaN
      !! stopped the run; 0 if the start had none)
      real(dp) :: e_singles = 0.0_dp     !! the t1 t1 part of the energy functional
      real(dp) :: e_doubles = 0.0_dp     !! the t2 part
      real(dp) :: de = 0.0_dp            !! E_k - E_{k-1} at the last iteration
      real(dp) :: step_max = 0.0_dp      !! max |t' - t| at the last iteration
      real(dp) :: step_norm2 = 0.0_dp    !! ||t' - t||_2 over the full t2 at the last iteration; not PySCF's normt (header)
      real(dp) :: t1_diagnostic = 0.0_dp  !! sqrt(||t1||^2 / (2 no)), PySCF's get_t1_diagnostic
      real(dp) :: seconds(RCCSD_NTIME) = 0.0_dp  !! device seconds per stage, RCCSD_TIME_*
      logical :: has_triples = .false.   !! (T) was run: `e_t`, `e_ccsd_t` and `seconds_triples` are valid only when true
      real(dp) :: e_t = 0.0_dp            !! the (T) correction; set by `rccsd_triples_finish` (trc_rccsd_triples), else 0
      real(dp) :: e_ccsd_t = 0.0_dp       !! e_ccsd + e_t, the CCSD(T) correlation energy; set with e_t
      real(dp) :: seconds_triples = 0.0_dp  !! host wall seconds of the (T), set with e_t
      character(len=MESSAGE_LEN) :: message = ""
   end type trc_rccsd_result_t

   type :: trc_rccsd_state_t
      !! What a run leaves for the caller (module header): the result, the integral blocks and the amplitudes
      type(trc_rccsd_result_t) :: result
      type(trc_rccsd_ints_t) :: ints   !! empty after a refusal
      type(trc_rccsd_amps_t) :: amps   !! the last iterate, host and device
      type(step_table_t) :: steps     !! time and counted flops per step of the last run (trc_rccsd_steps names the rows)
      real(dp), allocatable :: rank_seconds(:, :)
      !! over several ranks, rank_seconds(:, r + 1) of rank r, on every rank after the run (not allocated for a one-rank
      !! run): the seconds of the ladder (its passes, bt, the correction), of everything computed in the iterations (the
      !! rows tau to energy), of the two transfers, and of the wait at the ladder's agreement
      logical :: use_start = .false.   !! start from `amps`, sized and filled on the host by the caller; the run clears it
   contains
      procedure :: release => state_release
   end type trc_rccsd_state_t

   integer, parameter :: PK_NREAL = 8, PK_NFLAG = 9
   integer, parameter :: PK_SECONDS = PK_NREAL                      !! offset of the RCCSD_NTIME seconds in the result packet
   integer, parameter :: PK_FLAGS = PK_SECONDS + RCCSD_NTIME        !! offset of the PK_NFLAG counts and flags
   integer, parameter :: PK_TEXT = PK_FLAGS + PK_NFLAG              !! offset of the message, one code per character
   integer, parameter :: PACKET_LEN = PK_TEXT + MESSAGE_LEN
   !! the result broadcast: energies and measures, seconds, counts and flags, the message; every offset follows from
   !! the one before, and a slice assignment of the wrong size does not compile

   type :: run_work_t
      !! The parts of a run that never outlive it
      type(trc_rccsd_update_t) :: upd
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: newamps
   contains
      procedure :: release => work_release
   end type run_work_t

contains

   subroutine state_release(this)
      !! Frees the integrals and the amplitudes, host and device; safe on an empty type
      class(trc_rccsd_state_t), intent(inout) :: this

      call this%ints%release()
      call this%amps%release()
      this%use_start = .false.
   end subroutine state_release

   subroutine work_release(this)
      !! Frees everything the run allocated
      class(run_work_t), intent(inout) :: this

      call this%upd%release()
      call this%diis%release()
      call this%newamps%release()
   end subroutine work_release

   subroutine trc_rccsd_run(b, aux, pl, cmo, eps, res, opts, la)
      !! RI-CCSD correlation energy of an RHF reference; nothing is kept. `cmo` (nmo, nmo) and `eps` (nmo) are the
      !! RHF orbitals, occupied first, in the space of the basis (the SPHERICAL count for a spherical basis).
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_rccsd_result_t), intent(out) :: res
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la

      type(trc_rccsd_state_t) :: state

      call trc_rccsd_solve(b, aux, pl, cmo, eps, opts, la, state)
      res = state%result
      call state%release()
   end subroutine trc_rccsd_run

   subroutine trc_rccsd_solve(b, aux, pl, cmo, eps, opts, la, state)
      !! As `trc_rccsd_run`, but the result is `state%result` and `state%ints` and `state%amps` stay for the caller.
      !! With `state%use_start` the iteration starts from `state%amps` instead of MP2. The caller sizes them with
      !! `state%amps%init(no, nv)` (so they are device-mapped), fills t1 and t2 on the host and sets `use_start`;
      !! the run copies them to the device. Amplitudes that are unallocated or of another size are refused.
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state

      type(run_work_t) :: work
      type(timer_type) :: total, tc
      type(trc_rccsd_options_t) :: local
      real(dp), allocatable :: shared_cmo(:, :), shared_eps(:)
      logical :: use_start, multi, agreed
      integer :: rank, nranks

      rank = 0
      nranks = 1
      multi = .false.
      agreed = .false.
      if (opts%has_comm) then
         rank = opts%comm%rank()
         nranks = opts%comm%size()
         multi = nranks > 1
      end if
      local = opts
      if (rank /= 0) local%verbose = .false.
      use_start = state%use_start
      state%result = trc_rccsd_result_t()
      call state%ints%release()
      if (.not. use_start) call state%amps%release()
      call total%start()
      if (multi) then
         call state%steps%init("rccsd", RCCSD_MULTI_STEP_NAME, bytes_steps=[RCCSD_ST_BCAST, RCCSD_ST_REDUCE])
      else
         call state%steps%init("rccsd", RCCSD_STEP_NAME)
      end if
      if (allocated(state%rank_seconds)) deallocate (state%rank_seconds)
      ! the DF integral build is the first silent stage of a large run (the header below waits for its sizes)
      if (local%verbose) call global_logger%info("  rccsd: building the integrals")
      call state%steps%step_on(la, RCCSD_ST_INTEGRALS, tc)
      if (multi) then
         ! A guard (module header, RANKS): a collective SCF already gives every rank these orbitals; this protects
         ! a caller that hands ranks different ones (sign flips of a column)
         shared_cmo = cmo
         shared_eps = eps
         call rccsd_mpi_bcast_orbitals(opts%comm, shared_cmo, shared_eps, 0)
         ! the form of the ladder is checked where the options are, before any integral is built (the same agreement
         ! below makes a refusal of one rank the refusal of all); AUTO becomes the half form for the rest of the run
         call ladder_for_ranks(local%ladder, state%result)
         if (.not. state%result%refused) call build_integrals(b, aux, pl, shared_cmo, shared_eps, local, la, state)
         call agree_refusal(opts%comm, state%result)
         agreed = state%result%refused
      else
         call build_integrals(b, aux, pl, cmo, eps, local, la, state)
      end if
      call state%steps%step_off(la, RCCSD_ST_INTEGRALS, tc)
      if (.not. state%result%refused) call add_integral_pieces(state%ints, state%steps)
      if (.not. state%result%refused) then
         if (multi) then
            call run_over_ranks(local, la, state, work, use_start)
            agreed = state%result%refused
         else
            call run_stages(local, la, state, work, use_start)
         end if
      end if
      if (state%result%refused) then
         ! a refusal from inside run_stages comes after the blocks were built: give them back, as the header promises
         call state%ints%release()
         if (.not. use_start) call state%amps%release()
      end if
      call work%release()
      call la%sync()
      state%result%seconds(RCCSD_TIME_TOTAL) = total%get_elapsed_time()
      call fill_seconds(state%steps, state%result)
      ! an agreed refusal already left every rank with the same result and nothing to share
      if (multi .and. .not. agreed) call share_outcome(opts%comm, rank, state, use_start)
      state%use_start = .false.
      if (local%verbose) call log_summary(state%result, state%steps, nranks)
      if (local%verbose .and. multi .and. allocated(state%rank_seconds)) call log_ranks(state)
   end subroutine trc_rccsd_solve

   subroutine agree_refusal(comm, res)
      !! Collective. If any rank is refused, every rank becomes refused with the message of the lowest-numbered refusing
      !! rank, so none of them goes on to a stage or a collective call the refused rank will not make and the reason
      !! is not lost
      type(comm_t), intent(in) :: comm
      type(trc_rccsd_result_t), intent(inout) :: res

      logical :: mine, anyone
      character(len=MESSAGE_LEN) :: text

      mine = res%refused
      text = res%message
      if (mine .and. len_trim(text) == 0) text = "trc_rccsd: the run was refused"
      ! the collective call in a statement of its own: an `.and.` is free to skip an operand, and a rank that
      ! skipped the call would leave the others waiting
      anyone = rccsd_mpi_agree_text(comm, text, mine)
      if (anyone) then
         res%refused = .true.
         res%converged = .false.
         res%message = text
      end if
   end subroutine agree_refusal

   subroutine follower_amps(state, use_start)
      !! A rank that does not iterate: device-mapped amplitudes of the integrals' size for the broadcast to fill. The
      !! caller's own, when `use_start` and they have the size, are kept (the broadcast overwrites them).
      type(trc_rccsd_state_t), intent(inout) :: state
      logical, intent(in) :: use_start

      logical :: keep

      keep = .false.
      if (use_start) then
         keep = allocated(state%amps%t1) .and. allocated(state%amps%t2)
         if (keep) keep = state%amps%no == state%ints%no .and. state%amps%nv == state%ints%nv
      end if
      if (.not. keep) call state%amps%init(state%ints%no, state%ints%nv)
   end subroutine follower_amps

   subroutine share_outcome(comm, rank, state, use_start)
      !! Collective. Rank 0's `state%result` and amplitudes to every rank (module header, RANKS), and every rank's
      !! ladder seconds to every rank. A refused result ends here, with the other ranks' blocks and amplitudes given
      !! back as rank 0's were.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: rank
      type(trc_rccsd_state_t), intent(inout) :: state
      logical, intent(in) :: use_start

      real(dp) :: packet(PACKET_LEN)

      allocate (state%rank_seconds(4, comm%size()))
      state%rank_seconds = rccsd_mpi_gather_values(comm, [state%result%seconds(RCCSD_TIME_LADDER), &
                                                          sum(state%steps%t(RCCSD_ST_INTER:RCCSD_ST_ENERGY)), &
                                                          state%steps%t(RCCSD_ST_BCAST) + state%steps%t(RCCSD_ST_REDUCE), &
                                                          state%steps%t(RCCSD_ST_WAIT)])
      if (rank == 0) call result_to_packet(state%result, packet)
      call bcast(comm, packet, int(PACKET_LEN), 0)
      if (rank /= 0) call packet_to_result(packet, state%result)
      if (state%result%refused) then
         call state%ints%release()
         if (.not. use_start) call state%amps%release()
         return
      end if
      call rccsd_mpi_bcast_amps(comm, state%amps%t1, state%amps%t2, 0)
      if (rank /= 0) call state%amps%push()
   end subroutine share_outcome

   subroutine result_to_packet(res, packet)
      !! `res` as doubles, in the order `packet_to_result` reads them back (counts and flags are exact in a double)
      type(trc_rccsd_result_t), intent(in) :: res
      real(dp), intent(out) :: packet(PACKET_LEN)

      integer :: i

      packet = 0.0_dp
      packet(1:PK_NREAL) = [res%e_mp2, res%e_ccsd, res%e_singles, res%e_doubles, res%de, res%step_max, &
                            res%step_norm2, res%t1_diagnostic]
      packet(PK_SECONDS + 1:PK_SECONDS + RCCSD_NTIME) = res%seconds
      packet(PK_FLAGS + 1:PK_FLAGS + PK_NFLAG) = [merge(1.0_dp, 0.0_dp, res%converged), &
                                                  merge(1.0_dp, 0.0_dp, res%refused), &
                                                  merge(1.0_dp, 0.0_dp, res%finite), real(res%iterations, dp), &
                                                  real(res%no, dp), real(res%nv, dp), real(res%naux, dp), &
                                                  real(res%nfrozen, dp), real(res%diis_restarts, dp)]
      do i = 1, MESSAGE_LEN
         packet(PK_TEXT + i) = real(iachar(res%message(i:i)), dp)
      end do
   end subroutine result_to_packet

   subroutine packet_to_result(packet, res)
      !! The inverse of `result_to_packet`; the fields (T) fills are left at their defaults
      real(dp), intent(in) :: packet(PACKET_LEN)
      type(trc_rccsd_result_t), intent(out) :: res

      integer :: i
      real(dp) :: f(PK_NFLAG)

      res%e_mp2 = packet(1)
      res%e_ccsd = packet(2)
      res%e_singles = packet(3)
      res%e_doubles = packet(4)
      res%de = packet(5)
      res%step_max = packet(6)
      res%step_norm2 = packet(7)
      res%t1_diagnostic = packet(8)
      res%seconds = packet(PK_SECONDS + 1:PK_SECONDS + RCCSD_NTIME)
      f = packet(PK_FLAGS + 1:PK_FLAGS + PK_NFLAG)
      res%converged = f(1) > 0.5_dp
      res%refused = f(2) > 0.5_dp
      res%finite = f(3) > 0.5_dp
      res%iterations = nint(f(4))
      res%no = nint(f(5))
      res%nv = nint(f(6))
      res%naux = nint(f(7))
      res%nfrozen = nint(f(8))
      res%diis_restarts = nint(f(9))
      do i = 1, MESSAGE_LEN
         res%message(i:i) = achar(nint(packet(PK_TEXT + i)))
      end do
   end subroutine packet_to_result

   subroutine add_integral_pieces(ints, stab)
      !! The pieces of the integral build (time and flops each, trc_rccsd_ints) as rows under its own
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(step_table_t), intent(inout) :: stab

      integer :: k

      do k = RCCSD_ST_METRIC, RCCSD_ST_BLOCKS
         call stab%add_time(k, ints%seconds(k))
         call stab%add_flop(k, ints%flops(k))
      end do
   end subroutine add_integral_pieces

   subroutine fill_seconds(stab, res)
      !! `res%seconds` from the table's totals, in the stages `RCCSD_TIME_*` have always meant
      type(step_table_t), intent(in) :: stab
      type(trc_rccsd_result_t), intent(inout) :: res

      res%seconds(RCCSD_TIME_INTEGRALS) = stab%t(RCCSD_ST_INTEGRALS)
      res%seconds(RCCSD_TIME_INTER) = stab%t(RCCSD_ST_INTER)
      res%seconds(RCCSD_TIME_W) = stab%t(RCCSD_ST_W)
      res%seconds(RCCSD_TIME_LADDER) = stab%t(RCCSD_ST_LADDER) + stab%t(RCCSD_ST_DRESSING)
      res%seconds(RCCSD_TIME_AMPS) = stab%t(RCCSD_ST_T1) + stab%t(RCCSD_ST_T2) + stab%t(RCCSD_ST_RING)
      res%seconds(RCCSD_TIME_DIIS) = stab%t(RCCSD_ST_DIIS)
   end subroutine fill_seconds

   function start_is_pair_symmetric(t1, t2, no, nv) result(sym)
      !! True when tau = t2 + t1 t1 of these host amplitudes has tau(j,i,b,a) = tau(i,j,a,b) to round-off
      !! (t1 t1 always has, so only t2 is looked at). A NaN anywhere answers false: the general path assumes nothing.
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv), t2(no, no, nv, nv)
      logical :: sym

      real(dp), parameter :: REL = 1.0e-12_dp  !! far above the round-off of an update, far below a real asymmetry
      integer :: i, j, a, b
      real(dp) :: gap, big

      gap = 0.0_dp
      big = 0.0_dp
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  gap = max(gap, abs(t2(i, j, a, b) - t2(j, i, b, a)))
                  big = max(big, abs(t2(i, j, a, b)))
               end do
            end do
         end do
      end do
      sym = gap <= REL*big .and. is_finite(sum(t2)) .and. is_finite(sum(t1))
   end function start_is_pair_symmetric

   subroutine refuse(res, text)
      !! The run did not start
      type(trc_rccsd_result_t), intent(inout) :: res
      character(len=*), intent(in) :: text

      res%refused = .true.
      res%converged = .false.
      ! The integral layer's messages already name their module (trc_rccsd_ints: ...).
      if (len(text) >= 4) then
         if (text(1:4) == "trc_") then
            res%message = text
            return
         end if
      end if
      res%message = "trc_rccsd: "//text
   end subroutine refuse

   function options_ok(opts, res) result(ok)
      !! True when the options are usable; else `res` says which is not
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_rccsd_result_t), intent(inout) :: res
      logical :: ok

      ok = .false.
      if (opts%nocc < 1) then
         call refuse(res, "opts%nocc must be the number of doubly occupied orbitals (>= 1)")
      else if (opts%nfrozen < 0) then
         call refuse(res, "opts%nfrozen must be given as a count >= 0 (the count-by-element helper is in the API layer)")
      else if (opts%max_iterations < 1 .or. opts%diis_size < 0 .or. opts%diis_start < 1) then
         call refuse(res, "max_iterations and diis_start must be >= 1 and diis_size >= 0")
      else if (.not. (opts%conv_energy > 0.0_dp .and. opts%conv_amplitude > 0.0_dp)) then
         call refuse(res, "the convergence thresholds must be positive")
      else
         ok = .true.
      end if
   end function options_ok

   subroutine build_integrals(b, aux, pl, cmo, eps, opts, la, state)
      !! The options, then the fitted blocks into `state%ints`; a refusal is in `state%result`
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cmo(:, :), eps(:)
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state

      type(trc_rccsd_ints_options_t) :: bopts

      if (.not. options_ok(opts, state%result)) return
      bopts%nalpha = opts%nocc
      bopts%nbeta = merge(opts%nocc, opts%nbeta, opts%nbeta < 0)
      bopts%nfrozen = opts%nfrozen
      bopts%aux_block = opts%aux_block
      bopts%has_ecp = opts%has_ecp
      call state%ints%build(b, aux, pl, cmo, eps, la, bopts)
      if (state%ints%status%has_error()) then
         call refuse(state%result, state%ints%status%get_message())
         return
      end if
      state%result%no = state%ints%no
      state%result%nv = state%ints%nv
      state%result%naux = state%ints%naux
      state%result%nfrozen = state%ints%nfrozen
      if (opts%verbose) call log_header(state%result, opts)
   end subroutine build_integrals

   subroutine run_stages(opts, la, state, work, use_start)
      !! From built integrals to the end of the iteration of a one-rank run; may return anywhere, `work` is released
      !! by the caller
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work
      logical, intent(in) :: use_start

      logical :: go
      real(dp) :: e_start

      call stages_setup(opts, la, state, work, use_start, go)
      if (state%result%refused) return
      e_start = state%result%e_ccsd
      if (go) call run_iterations(opts, la, state, work, e_start)
      call stages_finish(state, work)
   end subroutine run_stages

   subroutine stages_finish(state, work)
      !! The converged (or last) amplitudes to the host, and what is measured on them
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(in) :: work

      call state%amps%pull()
      state%result%t1_diagnostic = sqrt(sum(state%amps%t1**2)/(2.0_dp*real(state%ints%no, dp)))
      state%result%diis_restarts = work%diis%restarts
   end subroutine stages_finish

   subroutine stages_setup(opts, la, state, work, use_start, go)
      !! The first rank's part before the iterations: the update, the starting amplitudes, the DIIS history and the
      !! energy of the start. A refusal is in `state%result`. `go` is false when no iteration is to run (a refusal, or
      !! a start that has no finite energy: then `state%result` says so)
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work
      logical, intent(in) :: use_start
      logical, intent(out) :: go

      integer :: no, nv
      real(dp) :: e_start
      logical :: diis_host
      integer(int64) :: need, free
      type(trc_rccsd_ladder_options_t) :: ladder

      go = .false.
      no = state%ints%no
      nv = state%ints%nv
      ladder = opts%ladder
      if (is_multi(opts)) call ladder_for_ranks(ladder, state%result)
      if (state%result%refused) return
      call work%upd%init(state%ints, ladder)
      if (work%upd%status%has_error()) then
         call refuse(state%result, work%upd%status%get_message())
         return
      end if
      if (is_multi(opts)) call share_set(opts, state, work)
      if (state%result%refused) return
      ! the step is named for the form that ran, which the options may have left to the selector
      state%steps%names(RCCSD_ST_LADDER) = merge("ladder, assembled", "ladder, direct   ", &
                                               work%upd%ls%form == RCCSD_LADDER_ASSEMBLED)
      if (work%upd%ls%form == RCCSD_LADDER_HALF) state%steps%names(RCCSD_ST_LADDER) = "ladder, half"
      if (use_start) then
         if (.not. allocated(state%amps%t1) .or. .not. allocated(state%amps%t2) .or. state%amps%no /= no .or. &
             state%amps%nv /= nv) then
            call refuse(state%result, "use_start: state%amps is not allocated with the (no, nv) of this calculation")
            return
         end if
         ! The half ladder halves its contraction with tau(j,i,d,c) = tau(i,j,c,d). Amplitudes of the caller's
         ! need not have it, and the terms of the t2 equation carry an asymmetry on, so such a run takes the
         ! general path throughout (one more contraction a pass).
         work%upd%tau_symmetric = start_is_pair_symmetric(state%amps%t1, state%amps%t2, no, nv)
         if (.not. work%upd%tau_symmetric .and. opts%verbose) then
            call global_logger%info("  rccsd: the starting t2 is not symmetric under (i,a) <-> (j,b); the ladder "// &
                                    "takes its general path")
         end if
         call state%amps%push()
      else
         call state%amps%init(no, nv)
      end if
      call work%newamps%init(no, nv)
      ! the history keeps half of t2 exactly when the amplitudes have the pair symmetry the ladder was told about
      ! the history on the host when it does not fit in DIIS_FREE_FRACTION of the free device memory, or when asked
      need = rccsd_diis_history_bytes(no, nv, opts%diis_size, work%upd%tau_symmetric)
      free = rccsd_diis_free_device_bytes()
      diis_host = opts%diis_host
      if (.not. diis_host .and. free > 0 .and. need > 0) then
         if (real(need, dp) > DIIS_FREE_FRACTION*real(free, dp)) then
            diis_host = .true.
            call global_logger%warning("trc_rccsd: the DIIS history ("//to_string(real(need, dp)/1.0e9_dp, "f0.1")// &
                                       " GB) does not fit in "//to_string(nint(100*DIIS_FREE_FRACTION))// &
                                       "% of the free device memory ("//to_string(real(free, dp)/1.0e9_dp, "f0.1")// &
                                       " GB): kept on the host, about 2% slower an iteration")
         end if
      end if
      call work%diis%init(no, nv, opts%diis_size, packed=work%upd%tau_symmetric, host=diis_host)
      if (work%diis%status%has_error()) then
         call refuse(state%result, work%diis%status%get_message())
         return
      end if

      state%result%e_mp2 = state%ints%mp2_energy()
      if (use_start) then
         call rccsd_inter_energy(state%amps, work%upd%scr, state%result%e_singles, state%result%e_doubles)
      else
         call rccsd_inter_mp2_guess(state%ints, state%amps, work%upd%scr, state%result%e_doubles)
         state%result%e_singles = 0.0_dp
      end if
      e_start = state%result%e_singles + state%result%e_doubles
      if (.not. is_finite(e_start)) then
         ! e_ccsd stays 0: no finite energy exists to report
         state%result%e_singles = 0.0_dp
         state%result%e_doubles = 0.0_dp
         state%result%finite = .false.
         state%result%message = "trc_rccsd: the starting amplitudes give a NaN or infinite energy; no finite energy was "// &
                                "obtained, stopped, not converged"
         call global_logger%warning(trim(state%result%message))
      else
         state%result%e_ccsd = e_start
         go = .true.
      end if
   end subroutine stages_setup

   function is_multi(opts) result(multi)
      !! The run is collective: a communicator of more than one rank
      type(trc_rccsd_options_t), intent(in) :: opts
      logical :: multi

      multi = .false.
      if (opts%has_comm) multi = opts%comm%size() > 1
   end function is_multi

   subroutine ladder_for_ranks(ladder, res)
      !! The ladder options of a run over several ranks (called before the integrals are built, and again by the setups as a
      !! backstop, where it changes nothing): only the half form can be cut by its virtual index d, so
      !! the automatic choice becomes the half form (the assembled and direct forms are what AUTO picks only when the
      !! counts favour them: the direct one for v > about 2 o^2, and the assembled one where one d of the half form
      !! does not fit its budget, which `plan_batch` never refuses) and an explicit assembled or direct one is
      !! refused, with the reason, rather than run on the first rank alone with the others idle
      type(trc_rccsd_ladder_options_t), intent(inout) :: ladder
      type(trc_rccsd_result_t), intent(inout) :: res

      if (ladder%form == RCCSD_LADDER_AUTO) ladder%form = RCCSD_LADDER_HALF
      if (ladder%form /= RCCSD_LADDER_HALF) then
         call refuse(res, "the ladder of a run over several ranks is cut by its virtual index and needs the half form; "// &
                     "use --cc-ladder auto (the half form), or one rank")
      end if
   end subroutine ladder_for_ranks

   subroutine share_set(opts, state, work)
      !! This rank's d range of the ladder (trc_rccsd_share) into the update; a refusal is in `state%result`
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work

      integer, allocatable :: lo(:), hi(:)

      allocate (lo(opts%comm%size()), hi(opts%comm%size()))
      call rccsd_share_ranges(state%ints%no, state%ints%nv, state%ints%naux, opts%comm%size(), lo, hi)
      call work%upd%set_share(opts%comm, opts%comm%rank() == 0, lo, hi)
      if (work%upd%status%has_error()) call refuse(state%result, work%upd%status%get_message())
   end subroutine share_set

   subroutine run_over_ranks(opts, la, state, work, use_start)
      !! The stages and iterations over the ranks of `opts%comm` (module header, RANKS). Every rank calls it. Each
      !! rank sets up its part; ONE agreement makes a refusal of any rank the refusal of all (before any rank enters a
      !! loop that waits for another); then rank 0 iterates and the others follow until its stop word.
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work
      logical, intent(in) :: use_start

      logical :: go, lead
      integer :: word
      real(dp) :: e_start

      lead = opts%comm%rank() == 0
      go = .false.
      if (lead) then
         call stages_setup(opts, la, state, work, use_start, go)
      else
         call follower_setup(opts, state, work, use_start)
      end if
      call agree_refusal(opts%comm, state%result)
      if (state%result%refused) return
      if (.not. lead) then
         call rccsd_dist_follow(opts%comm, work%upd, state%ints, state%amps, la, state%steps)
         return
      end if
      e_start = state%result%e_ccsd
      if (go) call run_iterations(opts, la, state, work, e_start)
      ! the followers wait for a word at the top of an iteration: tell them the loop is over, unless the ranks stopped
      ! together on a failure they agreed on
      if (.not. work%upd%agreed_stop) word = rccsd_mpi_control(opts%comm, 0)
      call stages_finish(state, work)
   end subroutine run_over_ranks

   subroutine follower_setup(opts, state, work, use_start)
      !! A rank other than the first before the iterations: its part of the update (tau and the ladder, not W), its
      !! share, and device-mapped amplitudes for the broadcast to fill
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work
      logical, intent(in) :: use_start

      type(trc_rccsd_ladder_options_t) :: ladder

      ladder = opts%ladder
      call ladder_for_ranks(ladder, state%result)
      if (state%result%refused) return
      call work%upd%init(state%ints, ladder, follower=.true.)
      if (work%upd%status%has_error()) then
         call refuse(state%result, work%upd%status%get_message())
         return
      end if
      call share_set(opts, state, work)
      if (state%result%refused) return
      call follower_amps(state, use_start)
   end subroutine follower_setup

   pure function is_finite(x) result(finite)
      !! False for a NaN and for an infinity: a NaN fails every comparison, so it is not <= huge
      real(dp), intent(in) :: x
      logical :: finite

      finite = abs(x) <= huge(1.0_dp)
   end function is_finite

   subroutine run_iterations(opts, la, state, work, e_start)
      !! The loop of the module header, filling `state%result`; `state%amps` is the iterate on entry and on exit
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work
      real(dp), intent(in) :: e_start

      type(timer_type) :: tc, iteration_clock
      real(dp) :: e_old, e_new, es, ed
      integer :: it, restarts_seen, word
      logical :: extrapolated, converged

      e_old = e_start
      restarts_seen = 0
      converged = .false.
      if (opts%verbose) call global_logger%info("  iter        E_corr(CCSD)            dE   max|step|    ||step||_2  diis     time")
      do it = 1, opts%max_iterations
         state%result%iterations = it
         call iteration_clock%start()
         ! over ranks: the followers get the control word and the amplitudes (a statement of its own, on every rank)
         if (work%upd%distributed) then
            word = rccsd_dist_start(opts%comm, .true., rccsd_dist_word(it, work%upd%tau_symmetric), work%upd, &
                                    state%amps, la, state%steps)
         end if
         call work%upd%step(state%ints, state%amps, work%newamps, la, state%steps)
         if (work%upd%status%has_error()) then
            state%result%message = "trc_rccsd: the amplitude update refused: "//work%upd%status%get_message()
            exit
         end if
         call state%steps%step_on(la, RCCSD_ST_DIIS, tc)
         call work%diis%push(la, state%amps, work%newamps, it >= opts%diis_start)
         state%result%step_max = work%diis%step_max
         state%result%step_norm2 = work%diis%step_norm2
         if (.not. work%diis%finite) then
            state%result%finite = .false.
            state%result%message = "trc_rccsd: the amplitude step is NaN or infinite at iteration "//to_string(it)// &
                                   "; stopped, not converged"
            call state%steps%step_off(la, RCCSD_ST_DIIS, tc)
            exit
         end if
         call work%diis%next(la, work%newamps, state%amps, extrapolated)
         call state%steps%step_off(la, RCCSD_ST_DIIS, tc)
         if (work%diis%restarts > restarts_seen) then
            restarts_seen = work%diis%restarts
            call global_logger%warning("trc_rccsd: iteration "//to_string(it)//": "//work%diis%status%get_message())
         end if
         call state%steps%step_on(la, RCCSD_ST_ENERGY, tc)
         call rccsd_inter_energy(state%amps, work%upd%scr, es, ed)
         call state%steps%step_off(la, RCCSD_ST_ENERGY, tc)
         e_new = es + ed
         state%result%de = e_new - e_old
         if (.not. (is_finite(e_new) .and. is_finite(state%result%de))) then
            state%result%finite = .false.
            state%result%message = "trc_rccsd: the energy is NaN or infinite at iteration "//to_string(it)// &
                                   "; stopped, not converged"
            exit
         end if
         state%result%e_singles = es
         state%result%e_doubles = ed
         state%result%e_ccsd = e_new
         call la%sync()
         if (opts%verbose) call log_iteration(it, state%result, extrapolated, work%diis%nstored, iteration_clock)
         ! both criteria; a NaN was excluded above and cannot satisfy either
         converged = abs(state%result%de) < opts%conv_energy .and. state%result%step_max < opts%conv_amplitude
         if (converged) exit
         e_old = e_new
      end do
      state%result%converged = converged
      if (converged) then
         state%result%message = "trc_rccsd: converged in "//to_string(state%result%iterations)//" iterations"
      else if (state%result%finite .and. len_trim(state%result%message) == 0) then
         state%result%message = "trc_rccsd: not converged in "//to_string(state%result%iterations)// &
                                " iterations (dE, max step, ||step||_2: "//to_string(state%result%de, "es10.2")//" "// &
                                to_string(state%result%step_max, "es10.2")//" "//to_string(state%result%step_norm2, "es10.2")//")"
      end if
      if (.not. converged) call global_logger%warning(trim(state%result%message))
   end subroutine run_iterations

   subroutine log_header(res, opts)
      !! The size of the problem, before the first iteration
      type(trc_rccsd_result_t), intent(in) :: res
      type(trc_rccsd_options_t), intent(in) :: opts

      character(len=200) :: line

      write (line, "(a,i0,a,i0,a,i0,a,i0,a)") "  rccsd: ", res%no, " active occupied (", res%nfrozen, " frozen), ", &
         res%nv, " virtual, ", res%naux, " auxiliary functions"
      call global_logger%info(trim(line))
      write (line, "(a,es8.1,a,es8.1,a,i0,a,i0)") "  rccsd: conv_energy ", opts%conv_energy, " conv_amplitude ", &
         opts%conv_amplitude, " max_iterations ", opts%max_iterations, " diis_size ", opts%diis_size
      call global_logger%info(trim(line))
   end subroutine log_header

   subroutine log_iteration(it, res, extrapolated, nstored, clock)
      !! One line of the iteration table
      integer, intent(in) :: it, nstored
      type(trc_rccsd_result_t), intent(in) :: res
      logical, intent(in) :: extrapolated
      type(timer_type), intent(in) :: clock

      character(len=200) :: line

      write (line, "(i6,f24.12,es14.3,es12.3,es14.3,i6,a,f9.2,a)") it, res%e_ccsd, res%de, res%step_max, res%step_norm2, &
         nstored, merge("*", " ", extrapolated), clock%get_elapsed_time(), " s"
      call global_logger%info(trim(line))
   end subroutine log_iteration

   subroutine log_summary(res, stab, nranks)
      !! The result, then the step table and the mean time of an iteration. Over several ranks the table is rank 0's
      !! (it labels itself "rank 0 of N"); no other rank prints one.
      type(trc_rccsd_result_t), intent(in) :: res
      type(step_table_t), intent(in) :: stab
      integer, intent(in) :: nranks

      character(len=200) :: line
      real(dp) :: seconds

      if (res%refused) return
      write (line, "(a,f20.12,a,f20.12)") "  rccsd: E_MP2 ", res%e_mp2, "   E_CCSD(corr) ", res%e_ccsd
      call global_logger%info(trim(line))
      write (line, "(a,f10.6,a,i0)") "  rccsd: T1 diagnostic ", res%t1_diagnostic, "   DIIS restarts ", res%diis_restarts
      call global_logger%info(trim(line))
      call stab%report(res%seconds(RCCSD_TIME_TOTAL), 0, nranks)
      if (res%iterations > 0) then
         ! the stages of the iteration only: the integrals are built once
         seconds = sum(stab%t(RCCSD_ST_INTER:RCCSD_ST_ENERGY))
         if (stab%n >= RCCSD_ST_REDUCE) seconds = seconds + stab%t(RCCSD_ST_BCAST) + stab%t(RCCSD_ST_WAIT) + stab%t(RCCSD_ST_REDUCE)
         write (line, "(a,f10.3,a,i0,a)") "  rccsd: mean time of an iteration ", seconds/real(res%iterations, dp), " s over ", &
            res%iterations, " iterations"
         call global_logger%performance(trim(line))
      end if
   end subroutine log_summary

   subroutine log_ranks(state)
      !! Over several ranks, one line per rank: its d range of the ladder and its share of the ladder's flops, the
      !! seconds of its ladder in all, and per iteration the seconds it computed, transferred and waited (the rows
      !! tau to energy; `RCCSD_ST_BCAST` and `_REDUCE`; `_WAIT`). The ranges are those of trc_rccsd_share, recomputed.
      type(trc_rccsd_state_t), intent(in) :: state

      integer, allocatable :: lo(:), hi(:)
      integer :: r, n
      real(dp) :: all_flops, its
      character(len=220) :: line

      n = size(state%rank_seconds, 2)
      its = real(max(state%result%iterations, 1), dp)
      allocate (lo(n), hi(n))
      call rccsd_share_ranges(state%result%no, state%result%nv, state%result%naux, n, lo, hi)
      all_flops = rccsd_share_ladder_flops(state%result%no, state%result%nv, state%result%naux, 1, state%result%nv)
      do r = 1, n
         write (line, "(a,i0,a,i0,a,i0,a,f6.1,a,f8.3,a,3f8.3,a)") "  rccsd: rank ", r - 1, ": ladder d = ", lo(r), "..", hi(r), &
            " (", 100.0_dp*rccsd_share_ladder_flops(state%result%no, state%result%nv, state%result%naux, lo(r), hi(r))/ &
            max(all_flops, 1.0_dp), " % of its flops) ", state%rank_seconds(1, r), " s in all; per iteration compute, "// &
            "transfer, wait:", state%rank_seconds(2:4, r)/its, " s"
         call global_logger%info(trim(line))
      end do
   end subroutine log_ranks

end module trc_rccsd_driver
