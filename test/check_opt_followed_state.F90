module check_opt_followed_state_suite
   !!
   !! The state watch of an excited-state optimisation (trc_dlfind_bridge's trc_opt_judge_followed_state), on synthetic
   !! roots, so that a crossing does not have to be found in a molecule: S1 and S2 that swap after a step is a change of
   !! the followed root's dominant orbital pair.
   !!  1. the first evaluation has no reference pair and is not compared;
   !!  2. the same pair, or only a degenerate partner's phase, is silent;
   !!  3. a swapped pair is counted and, by default, only warned about (the evaluation goes on);
   !!  4. in strict mode the same swap is returned as a failure that names the pair change, and is counted;
   !!  5. another root inside the gap window is a counted warning, never a failure, strict or not;
   !!  6. a gap window of zero switches the gap warning off, with two roots at the same energy;
   !!  7. a followed root that is not among the solved ones is an error.
   !! A test-drive suite; no SCF, no DL-FIND call.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_dlfind_bridge, only: trc_opt_watch_t, trc_opt_judge_followed_state
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_check_opt_followed_state

contains

   subroutine collect_check_opt_followed_state(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("the first evaluation is the reference and is not compared", test_first), &
                  new_unittest("the same pair is silent", test_same), &
                  new_unittest("a swapped pair is counted and only warned about", test_swap_warns), &
                  new_unittest("a swapped pair fails the evaluation in strict mode", test_swap_strict), &
                  new_unittest("a close root is a counted warning, never a failure", test_gap), &
                  new_unittest("a zero gap window switches the gap warning off", test_gap_off), &
                  new_unittest("a followed root that was not solved is an error", test_not_solved) &
                  ]
   end subroutine collect_check_opt_followed_state

   subroutine two_roots(w, followed, hole, particle)
      !! S1 (0.20 Eh) and S2 (0.30 Eh), each with its own dominant pair (hole, particle), weight 0.9.
      type(trc_opt_watch_t), intent(out) :: w
      integer, intent(in) :: followed, hole(2), particle(2)
      w%followed = followed
      w%omega = [0.20_dp, 0.30_dp]
      w%weight = [0.9_dp, 0.9_dp]
      w%hole = hole
      w%particle = particle
   end subroutine two_roots

   subroutine test_first(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      call two_roots(w, 1, [1, 2], [1, 1])
      w%strict = .true.
      why = trc_opt_judge_followed_state(w)
      call check(error, len_trim(why) == 0, "no reference pair yet, nothing to compare: "//trim(why))
      if (allocated(error)) return
      call check(error, w%n_state_changes == 0 .and. w%n_gap_warnings == 0, "nothing is counted at the first evaluation")
   end subroutine test_first

   subroutine test_same(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      call two_roots(w, 2, [1, 2], [1, 1])
      w%ref_hole = 2
      w%ref_particle = 1
      w%strict = .true.
      why = trc_opt_judge_followed_state(w)
      call check(error, len_trim(why) == 0, "the same pair as the reference: "//trim(why))
      if (allocated(error)) return
      call check(error, w%n_state_changes == 0, "nothing is counted")
   end subroutine test_same

   subroutine test_swap_warns(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      ! S1 followed; the reference was S1's pair (1,1), the evaluation has S2's pair (2,1) as root 1's.
      call two_roots(w, 1, [2, 1], [1, 1])
      w%ref_hole = 1
      w%ref_particle = 1
      why = trc_opt_judge_followed_state(w)
      call check(error, len_trim(why) == 0, "not strict: the evaluation goes on: "//trim(why))
      if (allocated(error)) return
      call check(error, w%n_state_changes == 1, "one state change counted, got "//to_string(w%n_state_changes))
      if (allocated(error)) return
      why = trc_opt_judge_followed_state(w)
      call check(error, w%n_state_changes == 2, "the counter accumulates over evaluations")
   end subroutine test_swap_warns

   subroutine test_swap_strict(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      call two_roots(w, 1, [2, 1], [1, 1])
      w%ref_hole = 1
      w%ref_particle = 1
      w%strict = .true.
      why = trc_opt_judge_followed_state(w)
      call check(error, index(why, "dominant orbital pair of root 1 changed") > 0, "strict: the pair change fails: '"// &
                 trim(why)//"'")
      if (allocated(error)) return
      call check(error, index(why, "--tddft-strict-root") > 0, "the message names the option")
      if (allocated(error)) return
      call check(error, w%n_state_changes == 1, "the change is counted too")
      if (allocated(error)) return
      ! The particle changing alone is a change as well.
      call two_roots(w, 1, [1, 1], [2, 1])
      w%ref_hole = 1
      w%ref_particle = 1
      w%strict = .true.
      why = trc_opt_judge_followed_state(w)
      call check(error, len_trim(why) > 0, "a changed particle fails too")
   end subroutine test_swap_strict

   subroutine test_gap(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      call two_roots(w, 1, [1, 2], [1, 1])
      w%omega = [0.200_dp, 0.203_dp]
      w%gap_warn = 0.005_dp
      w%strict = .true.
      w%ref_hole = 1
      w%ref_particle = 1
      why = trc_opt_judge_followed_state(w)
      call check(error, len_trim(why) == 0, "a gap is never an error: "//trim(why))
      if (allocated(error)) return
      call check(error, w%n_gap_warnings == 1 .and. w%n_state_changes == 0, "one gap warning, no state change")
      if (allocated(error)) return
      w%omega = [0.200_dp, 0.250_dp]
      why = trc_opt_judge_followed_state(w)
      call check(error, w%n_gap_warnings == 1, "a root 0.05 Eh away is outside a 0.005 window")
   end subroutine test_gap

   subroutine test_gap_off(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      call two_roots(w, 1, [1, 2], [1, 1])
      w%omega = [0.2_dp, 0.2_dp]
      w%gap_warn = 0.0_dp
      why = trc_opt_judge_followed_state(w)
      call check(error, len_trim(why) == 0 .and. w%n_gap_warnings == 0, "gap window 0: no warning, even for two roots at "// &
                 "the same energy")
   end subroutine test_gap_off

   subroutine test_not_solved(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_opt_watch_t) :: w
      character(len=256) :: why
      call two_roots(w, 3, [1, 2], [1, 1])
      why = trc_opt_judge_followed_state(w)
      call check(error, index(why, "not among the 2 solved roots") > 0, "root 3 of 2: '"//trim(why)//"'")
   end subroutine test_not_solved

end module check_opt_followed_state_suite

program check_opt_followed_state
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_followed_state_suite, only: collect_check_opt_followed_state
   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(collect_check_opt_followed_state, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_followed_state
