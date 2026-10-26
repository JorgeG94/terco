!
! --logger's effect on --opt and --opt-freq: this branch (feat/cli-help-
! topics-logger-opt-freq) added --logger but, until feat/opt-info-logging,
! it never touched the SCF driver's own per-iteration table
! (src/scf/trc_scf_driver.F90's "RMS(D)"/"|FDS-SDF|" header and its basis/
! timing lines, all gated by scf_config%talk -- app/main.F90's own
! trc_set_verbose call forced that on unconditionally unless --silent).
! These checks pin the fix:
!
!   - a PLAIN single point (no --opt/--hessian) at the default --logger
!     (info) still prints the SCF driver's own per-iteration table, EXACTLY
!     as before this branch -- --logger must never make an ordinary run
!     quieter. (This is the one case that regressed between this file's
!     first version and the fix: toggling trc_set_verbose globally at
!     startup, rather than only around --opt/--opt-freq's own loops,
!     silenced a plain run too.)
!   - --opt at the default --logger reads like a single SCF instead: a
!     header, then one row per step (step, E(total), dE, max/rms |g|, SCF
!     iteration count, wall time -- app/trc_dlfind_bridge.F90's own
!     consider_step/TRC_OPT_STEP_HEADER) and a closing rule, never the SCF
!     driver's own per-iteration table for any evaluation DL-FIND makes.
!   - at --logger verbose, today's full per-evaluation SCF output is back
!     for those same evaluations, on top of the step table (nothing is
!     REMOVED at verbose).
!   - --opt-freq at info: the opt step table, then trc_run_hessian's own
!     one-line-per-displacement progress (src/api/trc_capi.F90 -- a
!     library call, never gated on "talk" at all, so this was already
!     right), then the frequency table -- never a per-displacement SCF
!     table, and never one for the Hessian's own reference SCF (the
!     reconvergence at the optimised geometry, also run quiet).
!
! water/STO-3G, lbfgs (the CLI's default --opt-algorithm): small enough
! that both --opt and --opt-freq finish in well under a second.
!
! A test-drive suite.
!
module check_cli_opt_logging_suite
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_cli_opt_logging, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the --opt/--opt-freq
      !! output (<geometry>_opt.xyz/_traj.xyz land in the CLI's CWD).
   character(len=*), parameter :: SCF_TABLE_MARKER = 'RMS(D)'
      !! src/scf/trc_scf_driver.F90's own per-iteration header, nowhere
      !! else in the CLI's output -- the one unambiguous sign that a
      !! "talk"-gated SCF ran noisily.
   character(len=*), parameter :: STEP_TABLE_MARKER = 'SCF it'
      !! app/trc_dlfind_bridge.F90's own TRC_OPT_STEP_HEADER, equally
      !! unambiguous for "the opt step table printed".

contains

   subroutine collect_cli_opt_logging(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("a plain single point at the default --logger still prints the SCF "// &
                      "driver's own per-iteration table (--logger must not quiet an ordinary run)", &
                      test_plain_single_point_stays_noisy), &
         new_unittest("--opt at the default --logger (info): the step table, never the SCF "// &
                      "driver's own per-iteration table", test_opt_info_terse), &
         new_unittest("--opt --logger verbose: the full per-iteration SCF table is back, "// &
                      "alongside the step table", test_opt_verbose_full), &
         new_unittest("--opt-freq at info: the step table, the FD-Hessian displacement "// &
                      "progress and the frequency table, never a per-displacement (or reference) "// &
                      "SCF table", test_opt_freq_info) &
         ]
   end subroutine collect_cli_opt_logging

   subroutine setup_work_dir()
      integer :: rc
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
   end subroutine setup_work_dir

   subroutine test_plain_single_point_stays_noisy(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      logical :: has_scf_table

      call setup_work_dir()
      out_path = work_dir//'/plain.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json '// &
            '> plain.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call grep_contains(out_path, SCF_TABLE_MARKER, has_scf_table)
      if (.not. has_scf_table) then
         call test_failed(error, "a plain single point lost its SCF per-iteration table at the "// &
                          "default --logger level -- --logger must only gate --opt/--opt-freq's "// &
                          "own loop, never an ordinary run"); return
      end if
   end subroutine test_plain_single_point_stays_noisy

   subroutine test_opt_info_terse(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      logical :: has_scf_table, has_step_table, has_summary

      call setup_work_dir()
      out_path = work_dir//'/opt_info.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-maxiter 50 --opt-gtol 1e-5 > opt_info.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call grep_contains(out_path, SCF_TABLE_MARKER, has_scf_table)
      if (has_scf_table) then
         call test_failed(error, "the SCF driver's own per-iteration table ("//SCF_TABLE_MARKER// &
                          ") appeared at the default --logger level; it must only show at "// &
                          "verbose/debug for an evaluation inside --opt"); return
      end if
      call grep_contains(out_path, STEP_TABLE_MARKER, has_step_table)
      if (.not. has_step_table) then
         call test_failed(error, "no opt step table at the default --logger level"); return
      end if
      call grep_contains(out_path, 'geometry optimisation (--opt)', has_summary)
      if (.not. has_summary) then
         call test_failed(error, "no final optimisation summary at the default --logger level"); return
      end if
   end subroutine test_opt_info_terse

   subroutine test_opt_verbose_full(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      logical :: has_scf_table, has_step_table

      call setup_work_dir()
      out_path = work_dir//'/opt_verbose.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-maxiter 50 --opt-gtol 1e-5 --logger verbose > opt_verbose.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call grep_contains(out_path, SCF_TABLE_MARKER, has_scf_table)
      if (.not. has_scf_table) then
         call test_failed(error, "--logger verbose did not bring back the SCF driver's own "// &
                          "per-iteration table for --opt's own evaluations"); return
      end if
      call grep_contains(out_path, STEP_TABLE_MARKER, has_step_table)
      if (.not. has_step_table) then
         call test_failed(error, "--logger verbose dropped the opt step table -- it must stay, "// &
                          "on top of the full SCF output"); return
      end if
   end subroutine test_opt_verbose_full

   subroutine test_opt_freq_info(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      logical :: has_scf_table, has_step_table, has_displacement, has_freq_table

      call setup_work_dir()
      out_path = work_dir//'/opt_freq_info.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt-freq '// &
            '--opt-maxiter 50 --opt-gtol 1e-5 > opt_freq_info.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call grep_contains(out_path, SCF_TABLE_MARKER, has_scf_table)
      if (has_scf_table) then
         call test_failed(error, "the SCF driver's own per-iteration table appeared somewhere in "// &
                          "--opt-freq's output at the default --logger level (--opt's own loop, its "// &
                          "reconvergence SCF, or an FD-Hessian displacement)"); return
      end if
      call grep_contains(out_path, STEP_TABLE_MARKER, has_step_table)
      if (.not. has_step_table) then
         call test_failed(error, "no opt step table before the FD Hessian ran"); return
      end if
      call grep_contains(out_path, 'displacement', has_displacement)
      if (.not. has_displacement) then
         call test_failed(error, "no FD-Hessian displacement progress line"); return
      end if
      call grep_contains(out_path, 'vibrational analysis', has_freq_table)
      if (.not. has_freq_table) then
         call test_failed(error, "no final vibrational-analysis/thermochemistry block"); return
      end if
   end subroutine test_opt_freq_info

   subroutine grep_contains(path, needle, found_out)
      character(len=*), intent(in) :: path, needle
      logical, intent(out) :: found_out
      integer :: u, ios
      character(len=4096) :: line
      found_out = .false.
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, needle) > 0) then
            found_out = .true.
            exit
         end if
      end do
      close (u)
   end subroutine grep_contains

end module check_cli_opt_logging_suite

program check_cli_opt_logging
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_opt_logging_suite, only: collect_cli_opt_logging, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_opt_logging: usage: check_cli_opt_logging <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_opt_logging, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_opt_logging
