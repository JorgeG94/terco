module check_rccsd_xfer_suite
   !!
   !! The staging array of the RI-CCSD transfers (trc_rccsd_xfer), as a test-drive suite: `init` accepts a size whose
   !! full amplitudes o^2 v^2 and packed half fit a default integer and sizes the half as N (N + 1)/2, N = o v; it
   !! refuses a size whose packed half fits but whose full amplitudes do not (o = 50000, v = 1: N (N + 1)/2 = 1.25e9 <
   !! 2^31 - 1 < N^2 = 2.5e9, which the kernels index with default integers), and one so large that N^2 itself overflows
   !! 64 bits (o = v = 2^31 - 1), without allocating. WHAT THIS PINS: both refusals, and that a
   !! refused init leaves nothing allocated. It does not prove the old guard wrong: on code that checked only the packed
   !! half the o = 50000 case would try to allocate 10 GB and pass or fail depending on whether that succeeds, and for
   !! o = v = 2^31 - 1 a guard that forms N^2 overflows 64 bits (a wrapped value that happens to compare small would
   !! slip through, which is undefined, not measured). With the guard as written (N > huge/N) both are refused in
   !! constant time whatever the allocator does.
   !!
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_rccsd_xfer, only: trc_rccsd_xfer_t
   implicit none
   private
   public :: collect_rccsd_xfer

contains

   subroutine collect_rccsd_xfer(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("init sizes the packed half", test_init_small), &
                  new_unittest("init refuses o^2 v^2 above the default integer before allocating", test_init_oversize), &
                  new_unittest("init refuses a size whose o^2 v^2 overflows 64 bits", test_init_overflow) &
                  ]
   end subroutine collect_rccsd_xfer

   subroutine test_init_small(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_xfer_t) :: xf
      logical :: ok

      call xf%init(3, 4, ok)
      call check(error, ok, "a 3 x 4 staging array must fit")
      if (allocated(error)) return
      call check(error, xf%nhalf == 12*13/2, "nhalf is N (N + 1)/2 with N = 12")
      if (allocated(error)) return
      call check(error, allocated(xf%stage) .and. size(xf%stage) == xf%nhalf, "the staging array holds the half")
      call xf%release()
   end subroutine test_init_small

   subroutine test_init_oversize(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_xfer_t) :: xf
      logical :: ok

      call xf%init(50000, 1, ok)
      call check(error, .not. ok, "o^2 v^2 = 2.5e9 must be refused")
      if (allocated(error)) then
         call xf%release()
         return
      end if
      call check(error, .not. allocated(xf%stage) .and. xf%nhalf == 0, "a refused init leaves nothing allocated")
      call xf%release()
   end subroutine test_init_oversize

   subroutine test_init_overflow(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_rccsd_xfer_t) :: xf
      logical :: ok

      call xf%init(huge(0), huge(0), ok)
      call check(error, .not. ok, "o = v = 2^31 - 1 must be refused")
      if (allocated(error)) then
         call xf%release()
         return
      end if
      call check(error, .not. allocated(xf%stage) .and. xf%nhalf == 0, "a refused init leaves nothing allocated")
      call xf%release()
   end subroutine test_init_overflow

end module check_rccsd_xfer_suite

program check_rccsd_xfer
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_xfer_suite, only: collect_rccsd_xfer
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_rccsd_xfer, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_xfer
