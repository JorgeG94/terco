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
   call execute_command_line('cp -f hcn.xyz '//work_dir//'/hcn.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)

   ! water.xyz with its two H atoms swapped: refuse()'s own start geometry
   ! is always water.xyz, so linear interpolation toward THIS as the NEB
   ! product endpoint puts both H atoms on the SAME point (the midpoint
   ! of their swapped positions, water's own mirror plane) at the image
   ! halfway along the path -- a cheap, deterministic way to make
   ! check_neb_initial_path's "two atoms closer than 0.5 Bohr" pre-check
   ! actually fire (app/trc_dlfind_bridge.F90), not just reason about it.
   block
      integer :: u
      open (newunit=u, file=work_dir//'/water_hswap.xyz', status='replace', action='write')
      write (u, '(a)') '3'
      write (u, '(a)') 'water.xyz with its two H atoms swapped'
      write (u, '(a)') 'O   0.0000000   0.0000000   0.1173000'
      write (u, '(a)') 'H   0.0000000  -0.7572000  -0.4692000'
      write (u, '(a)') 'H   0.0000000   0.7572000  -0.4692000'
      close (u)
   end block

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

   ! --opt-algorithm stochastic/genetic: refused unconditionally --
   ! libdlfind's own dlf_parallel_opt.f90 segfaults (a missing `kiter`
   ! argument at every one of its three dlf_get_gradient calls, confirmed
   ! by reading the pinned source -- app/main.F90's own refusal message
   ! and app/trc_dlfind_bridge.F90's header).
   call refuse(cli_exe, work_dir, '--opt-algorithm stochastic', nfail)
   call refuse(cli_exe, work_dir, '--opt-algorithm genetic', nfail)

   ! --opt-neb: mutually exclusive with a saddle target/algorithm this
   ! bridge does not offer for a chain.
   call refuse(cli_exe, work_dir, '--opt-neb water.xyz --opt-target saddle', nfail)
   call refuse(cli_exe, work_dir, '--opt-neb water.xyz --opt-algorithm prfo', nfail)
   ! a NEB sub-option without --opt-neb itself
   call refuse(cli_exe, work_dir, '--opt-neb-images 10', nfail)
   call refuse(cli_exe, work_dir, '--opt-neb-climb', nfail)
   ! below DL-FIND's own minimum (3 -- dlf_neb.f90: "More than 2 images
   ! have to be used in NEB"); water.xyz as its own endpoint is fine here,
   ! this never reaches DL-FIND either way.
   call refuse(cli_exe, work_dir, '--opt-neb water.xyz --opt-neb-images 2', nfail)
   ! an unknown --opt-neb-endpoints
   call refuse(cli_exe, work_dir, '--opt-neb water.xyz --opt-neb-endpoints bogus', nfail)
   ! a missing endpoint file
   call refuse(cli_exe, work_dir, '--opt-neb does_not_exist.xyz', nfail)
   ! mismatched atoms/order: hcn.xyz (H,C,N) as the endpoint of a water
   ! (O,H,H) start -- different atom count AND different elements.
   call refuse(cli_exe, work_dir, '--opt-neb hcn.xyz', nfail)

   ! --opt-neb-via: an orphan sub-option, a missing file, a via frame
   ! with mismatched atoms (hcn.xyz again, same reason as above), and the
   ! "two atoms closer than 0.5 Bohr" pre-check actually firing (the
   ! water_hswap.xyz product above, no via point needed to trigger it).
   call refuse(cli_exe, work_dir, '--opt-neb-via hcn.xyz', nfail)
   call refuse(cli_exe, work_dir, '--opt-neb water.xyz --opt-neb-via does_not_exist.xyz', nfail)
   call refuse(cli_exe, work_dir, '--opt-neb water.xyz --opt-neb-via hcn.xyz', nfail)
   call refuse(cli_exe, work_dir, '--opt-neb water_hswap.xyz --opt-neb-images 10', nfail)

   ! --opt-frozen: refused together with pure DLC/DLC-TC (DL-FIND itself
   ! rejects a frozen atom there -- dlf_coords.f90's own check), and a bad
   ! index.
   call refuse(cli_exe, work_dir, '--opt-frozen 1 --opt-coords dlc', nfail)
   call refuse(cli_exe, work_dir, '--opt-frozen 1 --opt-coords dlc-tc', nfail)
   call refuse(cli_exe, work_dir, '--opt-frozen 0', nfail)
   call refuse(cli_exe, work_dir, '--opt-frozen 99', nfail)
   call refuse(cli_exe, work_dir, '--opt-frozen bogus', nfail)

   ! --opt-constrain: refused with Cartesian coordinates outright, and
   ! 'cartesian i' refused with DLC/DLC-TC specifically (dlf_hdlc_
   ! interface.f90's own hdlc%internal=(icoord>=3) check); plus bad types
   ! and bad indices.
   call refuse(cli_exe, work_dir, '--opt-constrain "angle 2 1 3"', nfail)   !! default --opt-coords cartesian
   call refuse(cli_exe, work_dir, '--opt-constrain "cartesian 1" --opt-coords dlc', nfail)
   call refuse(cli_exe, work_dir, '--opt-constrain "cartesian 1" --opt-coords dlc-tc', nfail)
   call refuse(cli_exe, work_dir, '--opt-constrain "bogus 1 2" --opt-coords hdlc', nfail)
   call refuse(cli_exe, work_dir, '--opt-constrain "bond 1 99" --opt-coords hdlc', nfail)
   call refuse(cli_exe, work_dir, '--opt-constrain "bond 1" --opt-coords hdlc', nfail)   !! too few indices

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
