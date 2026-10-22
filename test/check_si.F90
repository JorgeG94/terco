module check_si_suite
   !! trc_si's prefix choice, as a test-drive suite: each count is printed
   !! with the prefix that keeps its number in [1, 1000), unprefixed below
   !! 1, exa-scaled past it, never rounded up to "1000.00", and
   !! right-justified when a width is given. The cases are the RI-MP2
   !! table's own (2453056.2 GFLOP and friends).
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_si, only: si_value_t, si_string
   implicit none
   private
   public :: collect_si

contains

   subroutine collect_si(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("giga to peta", test_giga_to_peta), &
                  new_unittest("rates keep their unit", test_rates), &
                  new_unittest("under a thousand stays put", test_stays_put), &
                  new_unittest("prefix boundaries", test_boundaries), &
                  new_unittest("below one and past exa", test_extremes), &
                  new_unittest("width pads on the left", test_width) &
                  ]
   end subroutine collect_si

   subroutine test_giga_to_peta(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, si_string(2453056.2e9_dp, "FLOP"), "2.45 PFLOP")
   end subroutine test_giga_to_peta

   subroutine test_rates(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, si_string(53270.3e9_dp, "FLOP/s"), "53.27 TFLOP/s")
      if (allocated(error)) return
      call check(error, si_string(62.0e9_dp, "B/s"), "62.00 GB/s")
   end subroutine test_rates

   subroutine test_stays_put(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, si_string(866.2e9_dp, "FLOP/s"), "866.20 GFLOP/s")
      if (allocated(error)) return
      call check(error, si_string(124.7e9_dp, "FLOP"), "124.70 GFLOP")
   end subroutine test_stays_put

   subroutine test_boundaries(error)
      !! 1000 steps up a prefix; so does 999.996, which two decimals would
      !! otherwise print as "1000.00".
      type(error_type), allocatable, intent(out) :: error
      call check(error, si_string(1000.0_dp, "FLOP"), "1.00 kFLOP")
      if (allocated(error)) return
      call check(error, si_string(999.996e12_dp, "FLOP"), "1.00 PFLOP")
   end subroutine test_boundaries

   subroutine test_extremes(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, si_string(0.5_dp, "FLOP"), "0.50 FLOP")
      if (allocated(error)) return
      call check(error, si_string(3.0e21_dp, "FLOP"), "3000.00 EFLOP")
   end subroutine test_extremes

   subroutine test_width(error)
      !! The type and the one-call form agree, padded to the table's column.
      type(error_type), allocatable, intent(out) :: error
      type(si_value_t) :: q
      q = si_value_t(2453056.2e9_dp, "FLOP")
      call check(error, q%to_string(15), "     2.45 PFLOP")
      if (allocated(error)) return
      call check(error, si_string(2453056.2e9_dp, "FLOP", 15), q%to_string(15))
   end subroutine test_width

end module check_si_suite

program check_si
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_si_suite, only: collect_si
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_si, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_si
