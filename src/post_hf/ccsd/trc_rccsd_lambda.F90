!
! The closed-shell, density-fitted CCSD Lambda solver: from a converged RI-CCSD run (trc_rccsd_driver's state) to
! (l1, l2), in pyscf's convention (trc_rccsd_lambda_update, CONVENTION). The densities, the orbital relaxation,
! the gradient and (T) are not here.
!
! ALGORITHM, as pyscf/cc/ccsd_lambda.py `kernel`:
!     build the intermediates once                    trc_rccsd_lambda_inter
!     l <- (t1, t2)
!     l' = update(l)                                  trc_rccsd_lambda_update
!     l <- DIIS(history of l', l' - l)                trc_rccsd_diis, packed: l2 is pair-symmetric
! until BOTH
!     max |l' - l| < conv_amplitude    (default 1e-6, over l1 and l2)
!     ||l' - l||_2 < conv_norm         (default 1e-5, over l1 and the FULL l2)
! of the step before the extrapolation. PySCF tests the norm alone, over l1 and the packed half of l2, with
! tolerance 1e-5; the norm here counts every off-diagonal pair twice, so it is the stricter of the two. A step
! that is not finite stops the run as not converged.
!
! MEMORY. Everything is on the device except, when it does not fit, the DIIS history: after the update's arrays
! exist, a history above DIIS_FREE_FRACTION of the free device memory goes to the host (trc_rccsd_diis, HOST
! HISTORY). Whatever then does not fit is refused with the sizes, before the iterations. The CCSD's integrals and
! amplitudes stay resident (the solver reads them). Besides those, a solve of (no, nv, naux) = (25, 130, 630)
! held 5.0 GB (sampled), 2.6 GB of it the half ladder's work arrays at their default budget
! (`options%ladder%budget_bytes`). From the array sizes, integrals and amplitudes included: about 10 GB at
! (20, 265, 705) and 37 GB at (30, 400, 1100), of which the DIIS history is 9 GB; neither was run.
!
! MEASURED on one V100, adenine cc-pVDZ / cc-pVDZ-RIFIT (25, 130, 630), nvfortran 26.9: 0.140 s an iteration
! against 0.183 s for a CCSD iteration of the same run, 11 iterations; the intermediates 0.097 s. With the DIIS
! history forced onto the host an iteration took 0.216 s.
!
! REFUSALS (a message in `result%message`, nothing left allocated): a CCSD that was refused or did not converge, or
! released; one run over several ranks (this solver is one rank and does not cut its ladder); amplitudes whose t2
! is not pair-symmetric; options out of range; too little memory. An ECP, an open shell and an empty virtual
! space are refused by the CCSD itself. No MPI call.
module trc_rccsd_lambda
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use trc_step_table, only: step_table_t
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_device_memory, only: device_memory_set_cap, device_memory_free_bytes
   use trc_contract, only: contract_permute
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_rccsd_diis, only: trc_rccsd_diis_t, rccsd_diis_history_bytes
   use trc_rccsd_driver, only: trc_rccsd_state_t
   use trc_rccsd_ladder, only: trc_rccsd_ladder_options_t
   use trc_rccsd_lambda_inter, only: trc_rccsd_lambda_options_t, RCCSD_LAMBDA_STEP_NAME
   use trc_rccsd_lambda_update, only: trc_rccsd_lambda_update_t, trc_rccsd_lambda_update_options_t, &
                                      RCCSD_LUPD_NSTEP, RCCSD_LUPD_STEP_NAME
   implicit none
   private

   public :: trc_rccsd_lambda_solve_options_t, trc_rccsd_lambda_result_t, trc_rccsd_lambda_state_t
   public :: trc_rccsd_lambda_solve

   integer, parameter :: MESSAGE_LEN = 300  !! length of `trc_rccsd_lambda_result_t%message`
   real(dp), parameter :: DIIS_FREE_FRACTION = 0.7_dp  !! as trc_rccsd_driver's
   real(dp), parameter :: SYMMETRY_REL = 1.0e-10_dp
   !! t2(ijab) - t2(jiba) allowed, relative to max |t2|: a converged run's asymmetry is round-off
   integer, parameter, public :: RCCSD_LUPD_ST_DIIS = RCCSD_LUPD_NSTEP + 1  !! `state%steps`: the DIIS row

   type :: trc_rccsd_lambda_solve_options_t
      !! Plain values only
      real(dp) :: conv_amplitude = 1.0e-6_dp  !! max |l' - l| below this
      real(dp) :: conv_norm = 1.0e-5_dp       !! and ||l' - l||_2 below this
      integer :: max_iterations = 100
      integer :: diis_size = 8       !! history length; 0 or 1 is no DIIS
      integer :: diis_start = 1      !! first iteration whose step enters the history
      logical :: diis_host = .false.  !! force the DIIS history onto the host
      integer(int64) :: device_free_cap = 0_int64
      !! TEST HOOK, as trc_rccsd_options_t's: above 0, the free device memory every size check of the solve sees
      type(trc_rccsd_ladder_options_t) :: ladder  !! form and batch sizes of the ladder on l2
      type(trc_rccsd_lambda_options_t) :: inter   !! passed to the build of the intermediates
      logical :: verbose = .false.  !! one line per iteration and the step tables
   end type trc_rccsd_lambda_solve_options_t

   type :: trc_rccsd_lambda_result_t
      !! What a solve reports
      logical :: converged = .false.
      logical :: refused = .false.  !! the solve did not start: `message` says why
      logical :: finite = .true.    !! false when a NaN or an infinity stopped it
      logical :: diis_host = .false.  !! the DIIS history was on the host
      integer :: iterations = 0
      integer :: diis_restarts = 0
      real(dp) :: step_max = 0.0_dp     !! max |l' - l| at the last iteration
      real(dp) :: step_norm2 = 0.0_dp   !! ||l' - l||_2 at the last iteration
      real(dp) :: seconds_setup = 0.0_dp       !! the intermediates and the allocations
      real(dp) :: seconds_iterations = 0.0_dp  !! all iterations
      character(len=MESSAGE_LEN) :: message = ""
   end type trc_rccsd_lambda_result_t

   type :: trc_rccsd_lambda_state_t
      !! What a solve leaves for the caller
      type(trc_rccsd_lambda_result_t) :: result
      type(trc_rccsd_amps_t) :: l
      !! l1(i,a) in `l%t1`, l2(i,j,a,b) in `l%t2`, host and device; empty after a refusal, the last iterate otherwise
      type(step_table_t) :: inter_steps  !! the intermediates' build, rows RCCSD_LAMBDA_STEP_NAME
      type(step_table_t) :: steps        !! the iterations, rows RCCSD_LUPD_STEP_NAME and DIIS
   contains
      procedure :: release => state_release
   end type trc_rccsd_lambda_state_t

   type :: work_t
      !! The parts of a solve that never outlive it
      type(trc_rccsd_lambda_update_t) :: upd
      type(trc_rccsd_diis_t) :: diis
      type(trc_rccsd_amps_t) :: lnew
   end type work_t

contains

   subroutine state_release(this)
      !! Frees (l1, l2), host and device; safe on an empty type
      class(trc_rccsd_lambda_state_t), intent(inout) :: this

      call this%l%release()
   end subroutine state_release

   subroutine trc_rccsd_lambda_solve(cc, opts, la, state)
      !! Solves the Lambda equations of the converged RI-CCSD `cc` (module header). `cc%ints` and `cc%amps` are
      !! read, host and device copies of the amplitudes equal, as `trc_rccsd_solve` leaves them.
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_lambda_state_t), intent(inout) :: state

      type(work_t) :: work
      type(timer_type) :: clock

      call state%l%release()
      state%result = trc_rccsd_lambda_result_t()
      call state%inter_steps%init("rccsd lambda", RCCSD_LAMBDA_STEP_NAME)
      call state%steps%init("rccsd lambda", [character(len=len(RCCSD_LUPD_STEP_NAME)) :: RCCSD_LUPD_STEP_NAME, "DIIS"])
      if (.not. options_ok(opts, state%result)) return
      if (.not. ccsd_ok(cc, state%result)) return
      call device_memory_set_cap(opts%device_free_cap)
      call clock%start()
      call setup(cc, opts, la, state, work)
      call la%sync()
      state%result%seconds_setup = clock%get_elapsed_time()
      if (.not. state%result%refused) then
         call clock%start()
         call iterate(cc, opts, la, state, work)
         call la%sync()
         state%result%seconds_iterations = clock%get_elapsed_time()
         state%result%diis_restarts = work%diis%restarts
         call state%l%pull()
      else
         call state%l%release()
      end if
      call work%upd%release()
      call work%diis%release()
      call work%lnew%release()
      call device_memory_set_cap(0_int64)
      if (opts%verbose) call log_summary(state)
   end subroutine trc_rccsd_lambda_solve

   subroutine refuse(res, text)
      !! The solve did not start
      type(trc_rccsd_lambda_result_t), intent(inout) :: res
      character(len=*), intent(in) :: text

      res%refused = .true.
      res%converged = .false.
      res%message = "trc_rccsd_lambda: "//text
      ! the messages of the modules underneath already name them
      if (len(text) >= 4) then
         if (text(1:4) == "trc_") res%message = text
      end if
      call global_logger%warning(trim(res%message))
   end subroutine refuse

   function options_ok(opts, res) result(ok)
      !! True when the options are usable; else `res` says which is not
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      type(trc_rccsd_lambda_result_t), intent(inout) :: res
      logical :: ok

      ok = .false.
      if (opts%max_iterations < 1 .or. opts%diis_size < 0 .or. opts%diis_start < 1) then
         call refuse(res, "max_iterations and diis_start must be >= 1 and diis_size >= 0")
      else if (.not. (opts%conv_amplitude > 0.0_dp .and. opts%conv_norm > 0.0_dp)) then
         call refuse(res, "the convergence thresholds must be positive")
      else
         ok = .true.
      end if
   end function options_ok

   function ccsd_ok(cc, res) result(ok)
      !! True when `cc` is a converged one-rank CCSD whose blocks and amplitudes are there
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_result_t), intent(inout) :: res
      logical :: ok

      ok = .false.
      if (cc%result%refused) then
         call refuse(res, "the CCSD was refused: "//trim(cc%result%message))
         return
      end if
      if (allocated(cc%rank_seconds)) then
         call refuse(res, "the CCSD ran over several ranks: the Lambda solver is one rank only")
         return
      end if
      if (.not. cc%result%converged) then
         call refuse(res, "the CCSD did not converge: the Lambda equations need converged amplitudes")
         return
      end if
      if (.not. allocated(cc%ints%b_ov) .or. .not. allocated(cc%amps%t1) .or. .not. allocated(cc%amps%t2)) then
         call refuse(res, "the CCSD state holds no integrals or amplitudes (released?)")
         return
      end if
      if (cc%amps%no /= cc%ints%no .or. cc%amps%nv /= cc%ints%nv) then
         call refuse(res, "(no, nv) differ between the CCSD's integrals and amplitudes")
         return
      end if
      if (.not. pair_symmetric(cc%amps%t2, cc%amps%no, cc%amps%nv)) then
         call refuse(res, "t2 is not symmetric under (i,a) <-> (j,b), which the Lambda equations assume")
         return
      end if
      ok = .true.
   end function ccsd_ok

   function pair_symmetric(t2, no, nv) result(sym)
      !! t2(i,j,a,b) = t2(j,i,b,a) to round-off, on the host copy
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no, no, nv, nv)
      logical :: sym

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
      sym = gap <= SYMMETRY_REL*big
   end function pair_symmetric

   subroutine setup(cc, opts, la, state, work)
      !! The amplitude sets, the update (intermediates and work arrays), the start l = t and the DIIS history. A
      !! refusal is in `state%result`
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_lambda_state_t), intent(inout) :: state
      type(work_t), intent(inout) :: work

      type(trc_rccsd_lambda_update_options_t) :: uopts
      type(error_t) :: err
      integer :: no, nv

      no = cc%ints%no
      nv = cc%ints%nv
      call state%l%init(no, nv)
      if (state%l%status%has_error()) then
         call refuse(state%result, state%l%status%get_message())
         return
      end if
      call work%lnew%init(no, nv)
      if (work%lnew%status%has_error()) then
         call refuse(state%result, work%lnew%status%get_message())
         return
      end if
      uopts%ladder = opts%ladder
      uopts%inter = opts%inter
      call work%upd%init(cc%ints, cc%amps, la, uopts, state%inter_steps)
      if (work%upd%status%has_error()) then
         call refuse(state%result, work%upd%status%get_message())
         return
      end if
      call setup_diis(opts, no, nv, state%result, work)
      if (state%result%refused) return
      ! pyscf's start
      call contract_permute("ia->ia", [no, nv], cc%amps%t1, state%l%t1, err)
      call contract_permute("ijab->ijab", [no, no, nv, nv], cc%amps%t2, state%l%t2, err)
      if (err%has_error()) call refuse(state%result, err%get_message())
   end subroutine setup

   subroutine setup_diis(opts, no, nv, res, work)
      !! The history, on the host when asked or when it does not fit in DIIS_FREE_FRACTION of what is free now
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      integer, intent(in) :: no, nv
      type(trc_rccsd_lambda_result_t), intent(inout) :: res
      type(work_t), intent(inout) :: work

      integer(int64) :: need, free
      logical :: host

      need = rccsd_diis_history_bytes(no, nv, opts%diis_size, .true.)
      free = device_memory_free_bytes()
      host = opts%diis_host
      if (.not. host .and. free > 0 .and. need > 0) then
         if (real(need, dp) > DIIS_FREE_FRACTION*real(free, dp)) then
            host = .true.
            call global_logger%warning("trc_rccsd_lambda: the DIIS history ("//to_string(real(need, dp)/1.0e9_dp, "f0.3")// &
                                       " GB) does not fit in "//to_string(nint(100*DIIS_FREE_FRACTION))// &
                                       "% of the free device memory ("//to_string(real(free, dp)/1.0e9_dp, "f0.3")// &
                                       " GB): kept on the host")
         end if
      end if
      call work%diis%init(no, nv, opts%diis_size, packed=.true., host=host)
      if (work%diis%status%has_error()) then
         call refuse(res, work%diis%status%get_message())
         return
      end if
      res%diis_host = work%diis%host
   end subroutine setup_diis

   subroutine iterate(cc, opts, la, state, work)
      !! The loop of the module header; `state%l` is the iterate on entry and on exit
      type(trc_rccsd_state_t), intent(in) :: cc
      type(trc_rccsd_lambda_solve_options_t), intent(in) :: opts
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_lambda_state_t), intent(inout) :: state
      type(work_t), intent(inout) :: work

      type(timer_type) :: tc, one
      integer :: it, restarts_seen
      logical :: extrapolated, converged
      character(len=200) :: line

      restarts_seen = 0
      converged = .false.
      if (opts%verbose) call global_logger%info("  iter   max|step|    ||step||_2  diis      time")
      do it = 1, opts%max_iterations
         state%result%iterations = it
         call one%start()
         call work%upd%step(cc%ints, cc%amps, state%l, work%lnew, la, state%steps)
         if (work%upd%status%has_error()) then
            state%result%message = "trc_rccsd_lambda: the update refused: "//work%upd%status%get_message()
            exit
         end if
         call state%steps%step_on(la, RCCSD_LUPD_ST_DIIS, tc)
         call work%diis%push(la, state%l, work%lnew, it >= opts%diis_start)
         state%result%step_max = work%diis%step_max
         state%result%step_norm2 = work%diis%step_norm2
         if (.not. work%diis%finite) then
            state%result%finite = .false.
            state%result%message = "trc_rccsd_lambda: the step is NaN or infinite at iteration "//to_string(it)// &
                                   "; stopped, not converged"
            call state%steps%step_off(la, RCCSD_LUPD_ST_DIIS, tc)
            exit
         end if
         call work%diis%next(la, work%lnew, state%l, extrapolated)
         call state%steps%step_off(la, RCCSD_LUPD_ST_DIIS, tc)
         if (work%diis%restarts > restarts_seen) then
            restarts_seen = work%diis%restarts
            call global_logger%warning("trc_rccsd_lambda: iteration "//to_string(it)//": "//work%diis%status%get_message())
         end if
         if (opts%verbose) then
            write (line, "(i6,es12.3,es14.3,i6,a,f9.3,a)") it, state%result%step_max, state%result%step_norm2, &
               work%diis%nstored, merge("*", " ", extrapolated), one%get_elapsed_time(), " s"
            call global_logger%info(trim(line))
         end if
         ! both criteria; a NaN was excluded above and cannot satisfy either
         converged = state%result%step_max < opts%conv_amplitude .and. state%result%step_norm2 < opts%conv_norm
         if (converged) exit
      end do
      state%result%converged = converged
      if (converged) then
         state%result%message = "trc_rccsd_lambda: converged in "//to_string(state%result%iterations)//" iterations"
      else
         if (len_trim(state%result%message) == 0) then
            state%result%message = "trc_rccsd_lambda: not converged in "//to_string(state%result%iterations)// &
                                   " iterations (max step, ||step||_2: "//to_string(state%result%step_max, "es10.2")// &
                                   " "//to_string(state%result%step_norm2, "es10.2")//")"
         end if
         call global_logger%warning(trim(state%result%message))
      end if
   end subroutine iterate

   subroutine log_summary(state)
      !! The step tables and the mean time of an iteration
      type(trc_rccsd_lambda_state_t), intent(inout) :: state

      character(len=200) :: line

      if (state%result%refused) return
      ! the intermediates' table totals its own rows: the setup also maps every array of the solve
      call state%inter_steps%report(sum(state%inter_steps%t), 0, 1)
      call state%steps%report(state%result%seconds_iterations, 0, 1)
      write (line, "(a,f10.3,a)") "  rccsd lambda: setup (intermediates, allocations) ", state%result%seconds_setup, " s"
      call global_logger%performance(trim(line))
      if (state%result%iterations > 0) then
         write (line, "(a,f10.3,a,i0,a)") "  rccsd lambda: mean time of an iteration ", &
            state%result%seconds_iterations/real(state%result%iterations, dp), " s over ", state%result%iterations, &
            " iterations"
         call global_logger%performance(trim(line))
      end if
   end subroutine log_summary

end module trc_rccsd_lambda
