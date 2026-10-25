!
! Invalid --opt-* combinations are refused by app/main.F90's own CLI
! validation BEFORE trc_dlfind_optimize/api_dl_find ever runs -- the
! point of validating up front at all (DL-FIND itself answers an invalid
! setting with `error stop`, which is unacceptable to let through; see
! app/trc_dlfind_bridge.F90's own header). Every case here must exit
! non-zero with the SPECIFIC message app/main.F90's own `die()` call
! raises for it, and NEVER reach an "SCF :" line (the marker that
! validation passed and a run started) -- asserting only "non-zero exit"
! would let a missing input file or an unrelated crash pass as "refused"
! just as easily as the intended validation.
!
! A test-drive suite.
!
module check_opt_refusals_suite
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_opt_refusals, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the --silent output.

contains

   subroutine collect_opt_refusals(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-saddle-method dimer needs --opt-target saddle", test_dimer_needs_saddle), &
         new_unittest("--opt-connect needs --opt-target saddle", test_connect_needs_saddle), &
         new_unittest("--opt-connect with --opt-saddle-method prfo is refused (libdlfind crash)", &
                      test_connect_prfo_refused), &
         new_unittest("an unknown --opt-algorithm is refused", test_bad_algorithm), &
         new_unittest("an unknown --opt-coords is refused", test_bad_coords), &
         new_unittest("an unknown --opt-target is refused", test_bad_target), &
         new_unittest("an unknown --opt-saddle-method is refused", test_bad_saddle_method), &
         new_unittest("--opt-algorithm together with --opt-target saddle is refused", test_algorithm_with_saddle), &
         new_unittest("--opt-zero-modes above 3*natm is refused", test_zero_modes_too_large), &
         new_unittest("a negative --opt-zero-modes is refused", test_zero_modes_negative), &
         new_unittest("an unknown --opt-init-hessian is refused", test_bad_init_hessian), &
         new_unittest("an unknown --opt-hessian-update is refused", test_bad_hessian_update), &
         new_unittest("an unknown --opt-linesearch is refused", test_bad_linesearch), &
         new_unittest("a --opt-print-level DL-FIND does not accept is refused", test_bad_print_level), &
         new_unittest("--opt-init-hessian with an algorithm that ignores it is refused", &
                      test_init_hessian_wrong_algorithm), &
         new_unittest("--opt-hessian-update with an algorithm that ignores it is refused", &
                      test_hessian_update_wrong_algorithm), &
         new_unittest("--opt-max-updates with an algorithm that ignores it is refused", &
                      test_max_updates_wrong_algorithm), &
         new_unittest("--opt-zero-modes with an algorithm that ignores it is refused", &
                      test_zero_modes_wrong_algorithm), &
         new_unittest("--opt-soft-modes with an algorithm that ignores it is refused", &
                      test_soft_modes_wrong_algorithm), &
         new_unittest("--opt-algorithm without --opt is refused, not silently ignored", &
                      test_opt_option_without_opt), &
         new_unittest("--opt-algorithm prfo on a single atom is refused (no internal DOF)", &
                      test_single_atom_prfo), &
         new_unittest("--opt-max-updates 0 is accepted (DL-FIND's own 'always recompute')", &
                      test_max_updates_zero_accepted), &
         new_unittest("none of the refusal cases above started an SCF or left xyz output in test/ itself", &
                      test_no_stray) &
         ]
   end subroutine collect_opt_refusals

   subroutine test_dimer_needs_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-saddle-method dimer', 'needs --opt-target saddle', error)
   end subroutine test_dimer_needs_saddle

   subroutine test_connect_needs_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-connect', 'needs --opt-target saddle', error)
   end subroutine test_connect_needs_saddle

   subroutine test_connect_prfo_refused(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method prfo --opt-connect', &
                  'is refused: it crashes inside', error)
   end subroutine test_connect_prfo_refused

   subroutine test_bad_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm bogus', '--opt-algorithm must be', error)
   end subroutine test_bad_algorithm

   subroutine test_bad_coords(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-coords bogus', '--opt-coords must be', error)
   end subroutine test_bad_coords

   subroutine test_bad_target(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target bogus', '--opt-target must be', error)
   end subroutine test_bad_target

   subroutine test_bad_saddle_method(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method bogus', &
                  '--opt-saddle-method must be', error)
   end subroutine test_bad_saddle_method

   subroutine test_algorithm_with_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-algorithm lbfgs', &
                  'is not used with --opt-target', error)
   end subroutine test_algorithm_with_saddle

   subroutine test_zero_modes_too_large(error)
      type(error_type), allocatable, intent(out) :: error
      ! nzero out of range: 3*natm = 9 for water.
      call refuse('water.xyz', '--opt --opt-zero-modes 10', 'must be between 0 and 3*natm', error)
   end subroutine test_zero_modes_too_large

   subroutine test_zero_modes_negative(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-zero-modes -1', 'must be between 0 and 3*natm', error)
   end subroutine test_zero_modes_negative

   subroutine test_bad_init_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method prfo --opt-init-hessian bogus', &
                  '--opt-init-hessian must be', error)
   end subroutine test_bad_init_hessian

   subroutine test_bad_hessian_update(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method prfo --opt-hessian-update bogus', &
                  '--opt-hessian-update must be', error)
   end subroutine test_bad_hessian_update

   subroutine test_bad_linesearch(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-linesearch bogus', '--opt-linesearch must be', error)
   end subroutine test_bad_linesearch

   subroutine test_bad_print_level(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-print-level 3', '--opt-print-level must be', error)
   end subroutine test_bad_print_level

   !
   ! The Hessian-related options are meaningful only for P-RFO/NR (or
   ! --opt-target saddle, which forces one of P-RFO/dimer) -- refused
   ! for the CLI's own default algorithm (lbfgs) rather than silently
   ! ignored (app/main.F90's own opt_hessian_relevant/opt_zero_soft_relevant).
   !
   subroutine test_init_hessian_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-init-hessian terco', '--opt-init-hessian needs --opt-algorithm', error)
   end subroutine test_init_hessian_wrong_algorithm

   subroutine test_hessian_update_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-hessian-update none', '--opt-hessian-update needs --opt-algorithm', &
                  error)
   end subroutine test_hessian_update_wrong_algorithm

   subroutine test_max_updates_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-max-updates 10', '--opt-max-updates needs --opt-algorithm', error)
   end subroutine test_max_updates_wrong_algorithm

   subroutine test_zero_modes_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      ! In range (2 <= 3*natm=9) so the bound check does not fire first:
      ! this must be refused by the algorithm-relevance check instead.
      call refuse('water.xyz', '--opt --opt-zero-modes 2', '--opt-zero-modes needs --opt-algorithm', error)
   end subroutine test_zero_modes_wrong_algorithm

   subroutine test_soft_modes_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-soft-modes 0.01', '--opt-soft-modes needs --opt-algorithm', error)
   end subroutine test_soft_modes_wrong_algorithm

   !
   ! Every --opt-* option is parsed only inside `if (r%opt) then` --
   ! given without --opt itself, it must be refused, not silently
   ! ignored (app/main.F90's own check just before that block).
   !
   subroutine test_opt_option_without_opt(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt-algorithm lbfgs', '--opt-algorithm needs --opt', error)
   end subroutine test_opt_option_without_opt

   !
   ! A single atom has no internal degrees of freedom for P-RFO/NR/the
   ! dimer method's own eigenvector-following to search over.
   !
   subroutine test_single_atom_prfo(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('o_atom.xyz', '--opt --opt-algorithm prfo', 'need more than one atom', error)
   end subroutine test_single_atom_prfo

   !
   ! DL-FIND itself accepts maxupd=0 (dlf_formstep.f90's own
   ! iupd>=maxupd test, iupd starting at 0: it recomputes the Hessian
   ! fresh every step) -- NOT refused, unlike the cases above.
   !
   subroutine test_max_updates_zero_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      call execute_command_line('mkdir -p '//work_dir, exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -rf basis_sets '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      out_path = 'opt_refusal_zero.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-target saddle --opt-saddle-method prfo --opt-max-updates 0 --opt-maxiter 1 --silent > '// &
            out_path//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      ! --opt-maxiter 1 will not converge (exit non-zero), but NOT with
      ! the "--opt-max-updates must be at least 0" validation message --
      ! this is a normal non-convergence exit, never a CLI refusal.
      if (scan_for(work_dir//'/'//out_path, '--opt-max-updates must be')) then
         call test_failed(error, "--opt-max-updates 0 was refused at the CLI, but DL-FIND itself accepts it")
      end if
   end subroutine test_max_updates_zero_accepted

   !
   ! None of the refusal cases above ever reaches --opt's own
   ! write-the-xyz-outputs code (every one is refused before DL-FIND
   ! runs at all) -- checked anyway, for the same reason every other opt
   ! check checks it: a future refusal that stops catching an invalid
   ! combination should fail LOUDLY here, not just silently start
   ! optimising and littering test/ again.
   !
   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: exists
      inquire (file='water_opt.xyz', exist=exists)
      if (exists) then
         call test_failed(error, "water_opt.xyz was left in test/ itself"); return
      end if
      inquire (file='water_traj.xyz', exist=exists)
      if (exists) call test_failed(error, "water_traj.xyz was left in test/ itself")
   end subroutine test_no_stray

   !
   ! Runs `xyz --opt <flags> --silent` (flags must include --opt unless
   ! the case is testing an --opt-* option given WITHOUT --opt) and
   ! requires: a non-zero exit, no "SCF :" line (validation refused
   ! before any run started), and `expect_substr` present somewhere in
   ! the redirected stdout/stderr -- the specific message, not just "it
   ! failed somehow".
   !
   subroutine refuse(xyz, flags, expect_substr, error)
      character(len=*), intent(in) :: xyz, flags, expect_substr
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc, u, ios
      character(len=1024) :: line
      logical :: saw_scf_line, saw_message

      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f '//xyz//' '//work_dir//'/'//xyz, exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      out_path = work_dir//'/opt_refusal.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' '//xyz//' -b basis_sets/sto-3g.json '//flags// &
            ' --silent > opt_refusal.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)

      saw_scf_line = .false.
      saw_message = .false.
      open (newunit=u, file=out_path, status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'SCF              :') > 0) saw_scf_line = .true.
            if (index(line, expect_substr) > 0) saw_message = .true.
         end do
         close (u)
      end if

      if (rc == 0) then
         call test_failed(error, "'"//flags//"' exited 0 -- should have been refused"); return
      end if
      if (saw_scf_line) then
         call test_failed(error, "'"//flags//"' ran an SCF before being refused"); return
      end if
      if (.not. saw_message) then
         call test_failed(error, "'"//flags//"' was refused, but the output did not contain '"// &
                          expect_substr//"' -- it failed for a different reason than the one this case tests")
      end if
   end subroutine refuse

   !
   ! Whether `substr` appears on any line of `path` -- used where a case
   ! must assert a message's ABSENCE (test_max_updates_zero_accepted).
   !
   function scan_for(path, substr) result(found)
      character(len=*), intent(in) :: path, substr
      logical :: found
      integer :: u, ios
      character(len=1024) :: line
      found = .false.
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, substr) > 0) found = .true.
      end do
      close (u)
   end function scan_for

end module check_opt_refusals_suite

program check_opt_refusals
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_refusals_suite, only: collect_opt_refusals, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_refusals: usage: check_opt_refusals <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_refusals, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_refusals
