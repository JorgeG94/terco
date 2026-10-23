!
! Invalid --opt-* combinations are refused by app/main.F90's own CLI
! validation BEFORE trc_dlfind_optimize/api_dl_find ever runs -- the
! point of validating up front at all (DL-FIND itself answers an invalid
! setting with `error stop`, which is unacceptable to let through; see
! app/trc_dlfind_bridge.F90's own header). Every case here must exit
! non-zero with a clear message and NEVER reach an "SCF :" line (the
! marker that validation passed and a run started) -- same contract
! check_cli_opt's sibling refusal tests use.
!
program check_opt_refusals
   implicit none

   character(len=:), allocatable :: cli_exe, work_dir
   character(len=1024) :: arg
   integer :: rc, nfail

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_opt_refusals: usage: check_opt_refusals <terco CLI> <scratch dir>'
      stop 1
   end if
   ! Run with WORKING_DIRECTORY set to work_dir (not test/ itself, this
   ! program's own WORKING_DIRECTORY, ctest's convention): none of these
   ! cases ever reach --opt's own write-the-xyz-outputs code (every one is
   ! refused before DL-FIND runs at all), but the redirected --silent
   ! stdout still lands in the CLI's CWD, so it gets the same treatment as
   ! every other opt check -- check_opt_connect.F90's own pattern.
   call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)

   ! dimer with a minimum target (--opt-target defaults to minimum)
   call refuse(cli_exe, work_dir, '--opt-saddle-method dimer', nfail)
   ! --opt-connect with a minimum target
   call refuse(cli_exe, work_dir, '--opt-connect', nfail)
   ! --opt-connect with P-RFO: crashes inside libdlfind itself (confirmed
   ! with a debug build, see app/main.F90's own refusal message) --
   ! refused before DL-FIND ever runs, not left to crash the process.
   call refuse(cli_exe, work_dir, '--opt-target saddle --opt-saddle-method prfo --opt-connect', nfail)
   ! an unknown algorithm
   call refuse(cli_exe, work_dir, '--opt-algorithm bogus', nfail)
   ! an unknown coordinate system
   call refuse(cli_exe, work_dir, '--opt-coords bogus', nfail)
   ! an unknown target
   call refuse(cli_exe, work_dir, '--opt-target bogus', nfail)
   ! an unknown saddle method
   call refuse(cli_exe, work_dir, '--opt-target saddle --opt-saddle-method bogus', nfail)
   ! --opt-algorithm together with a saddle target
   call refuse(cli_exe, work_dir, '--opt-target saddle --opt-algorithm lbfgs', nfail)
   ! nzero out of range (3*natm = 9 for water)
   call refuse(cli_exe, work_dir, '--opt-zero-modes 10', nfail)
   call refuse(cli_exe, work_dir, '--opt-zero-modes -1', nfail)
   ! an unknown Hessian source/update/line search
   call refuse(cli_exe, work_dir, '--opt-init-hessian bogus', nfail)
   call refuse(cli_exe, work_dir, '--opt-hessian-update bogus', nfail)
   call refuse(cli_exe, work_dir, '--opt-linesearch bogus', nfail)
   ! a print level DL-FIND does not accept
   call refuse(cli_exe, work_dir, '--opt-print-level 3', nfail)

   ! None of the above ever runs DL-FIND, so none of them could have
   ! written water_opt.xyz/water_traj.xyz -- checked anyway, for the same
   ! reason every other opt check checks it: a future refusal that stops
   ! catching an invalid combination should fail LOUDLY here, not just
   ! silently start optimising and littering test/ again.
   block
      logical :: exists
      inquire (file='water_opt.xyz', exist=exists)
      call report(nfail, .not. exists, 'check_opt_refusals: water_opt.xyz was not left in test/ itself')
      inquire (file='water_traj.xyz', exist=exists)
      call report(nfail, .not. exists, 'check_opt_refusals: water_traj.xyz was not left in test/ itself')
   end block

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_refusals: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_refusals: RESULT: PASS'

contains

   subroutine refuse(exe, wdir, flags, nf)
      character(len=*), intent(in) :: exe, wdir, flags
      integer, intent(inout) :: nf
      character(len=:), allocatable :: cmd, out_path
      integer :: rc, u, ios
      character(len=1024) :: line
      logical :: saw_scf_line
      out_path = wdir//'/opt_refusal.out'
      cmd = 'cd '//wdir//' && '//trim(exe)//' water.xyz -b basis_sets/sto-3g.json --opt '//flags//' --silent > '// &
            'opt_refusal.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      saw_scf_line = .false.
      open (newunit=u, file=out_path, status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'SCF              :') > 0) saw_scf_line = .true.
         end do
         close (u)
      end if
      if (rc == 0) then
         print '(a,a)', '  FAIL  (exit 0, should have been refused) ', flags
         nf = nf + 1
      else
         print '(a,a)', '  PASS  (refused) ', flags
      end if
      if (saw_scf_line) then
         print '(a,a)', '  FAIL  (an SCF ran before the refusal) ', flags
         nf = nf + 1
      end if
   end subroutine refuse

   subroutine report(nf, cond, what)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (cond) then
         print '(a,a)', '  PASS  ', what
      else
         print '(a,a)', '  FAIL  ', what
         nf = nf + 1
      end if
   end subroutine report

end program check_opt_refusals
