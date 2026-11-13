module check_cli_stabilise_suite
   !!
   !! What the terco CLI says about the two SCF convergence aids, through the
   !! installed binary:
   !!
   !! 1. A run that does not converge and ends on a level shifted iteration
   !!    (--level-shift-until far below the commutator, three iterations) prints
   !!    the SCF's own message, which says its virtual orbital energies include
   !!    the shift. Only the library's trc_message carried it before.
   !!
   !! 2. O3+ UHF/6-31G* with --level-shift-until 1e-3 (never converges unless the
   !!    stall guard drops the shift) converges and says, once, that it dropped
   !!    the shift.
   !!
   !! A test-drive suite.
   !!
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   implicit none
   private
   public :: collect_cli_stabilise, cli_exe, cli_out

   character(len=:), allocatable :: cli_exe
   !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: cli_out
   !! Where the CLI output goes, from the command line.

contains

   subroutine collect_cli_stabilise(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("a non-converged run on a shifted iteration prints the SCF message", test_noconv_message), &
                  new_unittest("O3+ --level-shift-until 1e-3 converges and reports the dropped shift", test_dropped_line) &
                  ]
   end subroutine collect_cli_stabilise

   function run_cli(args) result(text)
      !! The merged stdout and stderr of `terco <args>`.
      character(len=*), intent(in) :: args
      character(len=:), allocatable :: text
      character(len=1024) :: line
      integer :: u, ios, rc
      call execute_command_line(trim(cli_exe)//" "//args//" > "//cli_out//" 2>&1", exitstat=rc, wait=.true.)
      text = ""
      open (newunit=u, file=cli_out, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         text = text//trim(line)//new_line("a")
      end do
      close (u)
   end function run_cli

   subroutine test_noconv_message(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      text = run_cli("water.xyz -b basis_sets/sto-3g.json --level-shift-until 1e-12 --iterations 3")
      call check(error, index(text, "NOT converged") > 0, "the run should not converge in 3 iterations:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text, "level shifted") > 0, &
                 "the SCF message about the shifted virtuals was not printed:"//new_line("a")//text)
   end subroutine test_noconv_message

   subroutine test_dropped_line(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: text
      integer :: first
      text = run_cli("o3.xyz -b basis_sets/6-31gs.json -c 1 -u 1 --level-shift-until 1e-3")
      call check(error, index(text, "NOT converged") == 0 .and. index(text, "iterations, converged") > 0, &
                 "O3+ with --level-shift-until 1e-3 did not converge:"//new_line("a")//text)
      if (allocated(error)) return
      first = index(text, "level_shift_until dropped the level shift")
      call check(error, first > 0, "the dropped shift was not reported:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text(first + 1:), "level_shift_until dropped the level shift") == 0, &
                 "the dropped shift was reported more than once")
   end subroutine test_dropped_line

end module check_cli_stabilise_suite

program check_cli_stabilise
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_stabilise_suite, only: collect_cli_stabilise, cli_exe, cli_out
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); cli_out = trim(arg)
   if (len(cli_exe) == 0 .or. len(cli_out) == 0) then
      write (error_unit, "(a)") "check_cli_stabilise: usage: check_cli_stabilise <terco CLI> <scratch output file>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_stabilise, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_stabilise
