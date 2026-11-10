!
! The closed-shell RI-CCSD driver: integrals, MP2 start, amplitude update with DIIS,
! convergence, and the converged amplitudes kept for the (T) step.
!
! ALGORITHM, per iteration (mqc's `run_czt_rccsd`, with PySCF's convergence test):
!
!     (t1', t2') = update(t1, t2)              trc_rccsd_update, from the SAME input amplitudes
!     e = t' - t                               the step; its size is measured on the device
!     (t1, t2) <- DIIS(history of t', e)       or t' itself while the history is short (trc_rccsd_diis)
!     E = E_corr(t1, t2)                       of the amplitudes the NEXT iteration starts from
!
! the start is t1 = 0, t2 = MP2 (E_MP2 is its energy), or the amplitudes the caller left in
! `state%amps` when `state%use_start` is set. DIIS starts pushing at iteration `diis_start` (default 1, as mqc)
! with the history length `diis_size` (default 8; 0 or 1 switches it off).
!
! CONVERGENCE. The run stops as CONVERGED when BOTH
!
!     |E_k - E_{k-1}| < conv_energy      (default 1e-8; E_0 is the energy of the starting amplitudes, as PySCF's e_prev)
!     max |t' - t| < conv_amplitude      (default 1e-6; the largest element of the step over t1 and t2)
!
! hold. mqc tests the energy alone, and only from the second iteration; the amplitude criterion is a deliberate
! addition. PySCF tests both, with conv_tol 1e-7 and conv_tol_normt 1e-5, the latter on the Euclidean norm
! ||t' - t||_2 of the same concatenated step vector, before the extrapolation: `result%step_norm2` is that
! number, and max|e| <= ||e||_2, so the default here (1e-6 on the max) is not looser than PySCF's 1e-5 on the
! 2-norm. A step or an energy that is NaN or infinite stops the run as NOT converged with a message, whatever
! the comparisons would say; the test is written so a NaN can never satisfy it.
!
! THE T1 DIAGNOSTIC is PySCF's `get_t1_diagnostic`: sqrt(||t1||_F^2 / (2 no)), no the number of correlated
! occupied orbitals (2 no correlated electrons), of the converged (or last) t1.
!
! THE HANDLE. `trc_rccsd_solve` leaves in `state` what the (T) step needs without recomputing: `state%ints` (the
! fitted blocks, on the host and the device, stamp unchanged) and `state%amps` (t1, t2 on both, the last iterate,
! the converged one when `state%result%converged`). The caller owns them and calls `state%release()`. A refused
! run leaves `state%ints` empty (its `status` says why) and, unless `use_start` supplied them, `state%amps` too; a run that
! was not converged leaves both in place. The update
! object, the DIIS history and the second set of amplitudes are private to the run and always released before it
! returns, on every path. `trc_rccsd_run` is the same run for a caller that only wants the numbers: it releases
! everything.
!
! TIMING. Each stage is timed on the device (`la%sync()` before the clock is read, since cuBLAS returns before its
! kernels end), summed over the iterations, in `result%seconds`; `verbose` prints one line per iteration and the
! table at the end through global_logger.
!
! The arguments are those of `trc_rimp2_run` except that the doubly occupied count is `opts%nocc` (more than eight
! dummies is a rule of this tree) and the options carry what trc_rccsd_ints refuses on: more than one rank, an ECP,
! an open shell (`nbeta`). Frozen core: the caller passes the number (`opts%nfrozen`); the count-by-element helper
! lives in the API layer, which this one may not use. Single rank only.
module trc_rccsd_driver
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_linalg, only: trc_linalg_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, rccsd_inter_energy, rccsd_inter_mp2_guess
   use trc_rccsd_ladder, only: trc_rccsd_ladder_options_t
   use trc_rccsd_update, only: trc_rccsd_update_t, RCCSD_STAGE_INTER, RCCSD_STAGE_W, RCCSD_STAGE_LADDER, RCCSD_STAGE_AMPS
   use trc_rccsd_diis, only: trc_rccsd_diis_t
   implicit none
   private

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
      integer :: nranks = 1               !! ranks of the caller's communicator; more than one is refused
      logical :: has_ecp = .false.        !! the calculation carries an effective core potential: refused
      real(dp) :: conv_energy = 1.0e-8_dp     !! |dE| below this, and the next, converge the run
      real(dp) :: conv_amplitude = 1.0e-6_dp  !! max |t' - t| below this
      integer :: max_iterations = 100
      integer :: diis_size = 8            !! history length; 0 or 1 is no DIIS
      integer :: diis_start = 1           !! first iteration whose step enters the history
      integer :: aux_block = 0            !! auxiliary functions per (mu nu|P) block; 0 means all at once
      type(trc_rccsd_ladder_options_t) :: ladder  !! form and batch sizes of the particle ladder, passed through
      logical :: verbose = .false.        !! one line per iteration and the timing table
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
      real(dp) :: step_norm2 = 0.0_dp    !! ||t' - t||_2 at the last iteration: PySCF's normt
      real(dp) :: t1_diagnostic = 0.0_dp  !! sqrt(||t1||^2 / (2 no)), PySCF's get_t1_diagnostic
      real(dp) :: seconds(RCCSD_NTIME) = 0.0_dp  !! device seconds per stage, RCCSD_TIME_*
      character(len=300) :: message = ""
   end type trc_rccsd_result_t

   type :: trc_rccsd_state_t
      !! What a run leaves for the caller (module header): the result, the integral blocks and the amplitudes
      type(trc_rccsd_result_t) :: result
      type(trc_rccsd_ints_t) :: ints   !! empty after a refusal
      type(trc_rccsd_amps_t) :: amps   !! the last iterate, host and device
      logical :: use_start = .false.   !! start from `amps`, sized and filled on the host by the caller; the run clears it
   contains
      procedure :: release => state_release
   end type trc_rccsd_state_t

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
      type(timer_type) :: total
      logical :: use_start

      use_start = state%use_start
      state%result = trc_rccsd_result_t()
      call state%ints%release()
      if (.not. use_start) call state%amps%release()
      call total%start()
      call build_integrals(b, aux, pl, cmo, eps, opts, la, state)
      if (.not. state%result%refused) call run_stages(opts, la, state, work, use_start)
      if (state%result%refused) then
         ! a refusal from inside run_stages comes after the blocks were built: give them back, as the header promises
         call state%ints%release()
         if (.not. use_start) call state%amps%release()
      end if
      call work%release()
      call la%sync()
      state%result%seconds(RCCSD_TIME_TOTAL) = total%get_elapsed_time()
      state%use_start = .false.
      if (opts%verbose) call log_summary(state%result)
   end subroutine trc_rccsd_solve

   subroutine refuse(res, text)
      !! The run did not start
      type(trc_rccsd_result_t), intent(inout) :: res
      character(len=*), intent(in) :: text

      res%refused = .true.
      res%converged = .false.
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
      type(timer_type) :: clock

      if (.not. options_ok(opts, state%result)) return
      bopts%nalpha = opts%nocc
      bopts%nbeta = merge(opts%nocc, opts%nbeta, opts%nbeta < 0)
      bopts%nfrozen = opts%nfrozen
      bopts%aux_block = opts%aux_block
      bopts%nranks = opts%nranks
      bopts%has_ecp = opts%has_ecp
      call clock%start()
      call state%ints%build(b, aux, pl, cmo, eps, la, bopts)
      call la%sync()
      state%result%seconds(RCCSD_TIME_INTEGRALS) = clock%get_elapsed_time()
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
      !! From built integrals to the end of the iteration; may return anywhere, `work` is released by the caller
      type(trc_rccsd_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_state_t), intent(inout) :: state
      type(run_work_t), intent(inout) :: work
      logical, intent(in) :: use_start

      integer :: no, nv
      real(dp) :: e_start

      no = state%ints%no
      nv = state%ints%nv
      call work%upd%init(state%ints, opts%ladder)
      if (work%upd%status%has_error()) then
         call refuse(state%result, work%upd%status%get_message())
         return
      end if
      if (use_start) then
         if (.not. allocated(state%amps%t1) .or. .not. allocated(state%amps%t2) .or. state%amps%no /= no .or. &
             state%amps%nv /= nv) then
            call refuse(state%result, "use_start: state%amps is not allocated with the (no, nv) of this calculation")
            return
         end if
         call state%amps%push()
      else
         call state%amps%init(no, nv)
      end if
      call work%newamps%init(no, nv)
      call work%diis%init(no, nv, opts%diis_size)
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
         call run_iterations(opts, la, state, work, e_start)
      end if
      call state%amps%pull()
      state%result%t1_diagnostic = sqrt(sum(state%amps%t1**2)/(2.0_dp*real(no, dp)))
      state%result%diis_restarts = work%diis%restarts
      state%result%seconds(RCCSD_TIME_INTER) = work%upd%seconds(RCCSD_STAGE_INTER)
      state%result%seconds(RCCSD_TIME_W) = work%upd%seconds(RCCSD_STAGE_W)
      state%result%seconds(RCCSD_TIME_LADDER) = work%upd%seconds(RCCSD_STAGE_LADDER)
      state%result%seconds(RCCSD_TIME_AMPS) = work%upd%seconds(RCCSD_STAGE_AMPS)
   end subroutine run_stages

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

      type(timer_type) :: clock, iteration_clock
      real(dp) :: e_old, e_new, es, ed
      integer :: it, restarts_seen
      logical :: extrapolated, converged

      e_old = e_start
      restarts_seen = 0
      converged = .false.
      if (opts%verbose) call global_logger%info("  iter        E_corr(CCSD)            dE   max|step|    ||step||_2  diis     time")
      do it = 1, opts%max_iterations
         state%result%iterations = it
         call iteration_clock%start()
         call work%upd%step(state%ints, state%amps, work%newamps, la)
         if (work%upd%status%has_error()) then
            state%result%message = "trc_rccsd: the amplitude update refused: "//work%upd%status%get_message()
            exit
         end if
         call la%sync()
         call clock%start()
         call work%diis%push(la, state%amps, work%newamps, it >= opts%diis_start)
         state%result%step_max = work%diis%step_max
         state%result%step_norm2 = work%diis%step_norm2
         if (.not. work%diis%finite) then
            state%result%finite = .false.
            state%result%message = "trc_rccsd: the amplitude step is NaN or infinite at iteration "//to_string(it)// &
                                   "; stopped, not converged"
            call la%sync()
            state%result%seconds(RCCSD_TIME_DIIS) = state%result%seconds(RCCSD_TIME_DIIS) + clock%get_elapsed_time()
            exit
         end if
         call work%diis%next(la, work%newamps, state%amps, extrapolated)
         call la%sync()
         state%result%seconds(RCCSD_TIME_DIIS) = state%result%seconds(RCCSD_TIME_DIIS) + clock%get_elapsed_time()
         if (work%diis%restarts > restarts_seen) then
            restarts_seen = work%diis%restarts
            call global_logger%warning("trc_rccsd: iteration "//to_string(it)//": "//work%diis%status%get_message())
         end if
         call rccsd_inter_energy(state%amps, work%upd%scr, es, ed)
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

   subroutine log_summary(res)
      !! The result, then the stage times
      type(trc_rccsd_result_t), intent(in) :: res

      character(len=200) :: line
      character(len=16), parameter :: NAMES(RCCSD_NTIME) = [character(len=16) :: "integrals", "intermediates", "W", &
                                                           "ladder", "T1/T2 equations", "DIIS", "total"]
      integer :: k

      if (res%refused) return
      write (line, "(a,f20.12,a,f20.12)") "  rccsd: E_MP2 ", res%e_mp2, "   E_CCSD(corr) ", res%e_ccsd
      call global_logger%info(trim(line))
      write (line, "(a,f10.6,a,i0)") "  rccsd: T1 diagnostic ", res%t1_diagnostic, "   DIIS restarts ", res%diis_restarts
      call global_logger%info(trim(line))
      do k = 1, RCCSD_NTIME
         write (line, "(a,a16,f12.3,a)") "  rccsd time ", NAMES(k), res%seconds(k), " s"
         call global_logger%performance(trim(line))
      end do
   end subroutine log_summary

end module trc_rccsd_driver
