!
! step_table_t%report's layout: a table where no step counted operations or
! bytes prints the time column alone, so a run's XC- and RI-MP2-gradient
! tables do not show a whole column of dashes that reads as a failed count;
! a table where some steps count keeps the work and rate columns, with a
! dash on the steps that do not. The report goes through pic_logger, so the
! logger is pointed at a file and the file is read back.
!
! A test-drive suite.
!
module check_step_table_report_suite
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger, performance_level
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_step_table, only: step_table_t
   implicit none
   private
   public :: collect_step_table_report

   character(len=*), parameter :: LOG_FILE = "check_step_table_report.log"

contains

   subroutine collect_step_table_report(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("nothing counted: time column only", test_uncounted), &
                   new_unittest("something counted: work and rate kept, dashes where missing", test_counted)]
   end subroutine collect_step_table_report

   subroutine test_uncounted(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      character(len=:), allocatable :: text
      call tbl%init("uncounted", [character(len=8) :: "collect", "scatter"])
      call tbl%add_time(1, 0.010_dp)
      call tbl%add_time(2, 0.015_dp)
      text = report_text(tbl)
      call check(error, index(text, "uncounted: step") > 0, "no header line:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text, "work") == 0 .and. index(text, "rate") == 0, &
                 "the header still names work/rate columns:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text, " -") == 0, "a row still has a dash:"//new_line("a")//text)
   end subroutine test_uncounted

   subroutine test_counted(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      character(len=:), allocatable :: text
      call tbl%init("counted", [character(len=8) :: "gemm", "pack"])
      call tbl%add_time(1, 0.5_dp)
      call tbl%add_flop(1, 1000000000_int64)
      call tbl%add_time(2, 0.1_dp)
      text = report_text(tbl)
      call check(error, index(text, "work") > 0 .and. index(text, "rate") > 0, &
                 "the header lost its work/rate columns:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text, "FLOP/s") > 0, "the counted step has no rate:"//new_line("a")//text)
      if (allocated(error)) return
      call check(error, index(text, " -") > 0, "the uncounted step has no dash:"//new_line("a")//text)
   end subroutine test_counted

   function report_text(tbl) result(text)
      !! What `tbl%report` writes, read back from the log file.
      type(step_table_t), intent(in) :: tbl
      character(len=:), allocatable :: text
      character(len=512) :: line
      integer :: u, ios
      open (newunit=u, file=LOG_FILE, status="replace")
      close (u)
      call global_logger%configure_file_output(LOG_FILE, level=performance_level)
      call tbl%report(1.0_dp, 0, 1)
      call global_logger%close_log_file()
      text = ""
      open (newunit=u, file=LOG_FILE, status="old", action="read")
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         text = text//trim(line)//new_line("a")
      end do
      close (u)
   end function report_text

end module check_step_table_report_suite

program check_step_table_report
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_step_table_report_suite, only: collect_step_table_report
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_step_table_report, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_step_table_report
