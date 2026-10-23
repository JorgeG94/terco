!
! `terco water.xyz -b basis_sets/sto-3g.json --opt --json <probe>`, end to end: the
! JSON "optimization" block, and that <geometry>_opt.xyz/_traj.xyz were
! actually written with the right shape, in the style of
! check_cli_grad/check_cli_hessian.
! Numerical agreement with an external reference is check_opt_vs_pyscf.F90
! (tight) and check_opt_vs_mqc.F90 (loose sanity); this one is the CLI
! plumbing: does --opt --json write what it says it writes, in a shape a
! caller could actually parse.
!
program check_cli_opt
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   implicit none

   character(len=:), allocatable :: cli_exe, probe_json, probe_out, cmd, work_dir
   character(len=1024) :: arg
   integer :: rc, nfail, natm_json, steps, evals, u, ios
   logical :: converged, found, exists_opt, exists_traj
   type(json_file) :: json
   real(dp) :: final_energy, gmax, grms
   real(dp), allocatable :: row(:), traj_e(:)
   character(len=512) :: line
   character(len=:), allocatable :: coord_str

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   call get_command_argument(3, arg); work_dir = trim(arg)   !! a scratch dir this check owns, for the xyz outputs
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_cli_opt: usage: check_cli_opt <terco CLI> <scratch .json> <scratch dir>'
      stop 1
   end if
   call execute_command_line('mkdir -p '//work_dir, exitstat=rc, wait=.true.)
   call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
   probe_out = probe_json//'.out'

   cmd = trim(cli_exe)//' '//work_dir//'/water.xyz -b basis_sets/sto-3g.json --opt --opt-maxiter 50 --silent --json '// &
         probe_json//' > '//probe_out//' 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   call expect(rc == 0, 'check_cli_opt: the CLI run itself failed')

   call json%initialize()
   call json%load(filename=probe_json)
   call expect(.not. json%failed(), 'check_cli_opt: cannot read '//probe_json)
   call json%get('optimization.converged', converged, found)
   call expect(found, 'check_cli_opt: no optimization.converged')
   call report(nfail, converged, 'check_cli_opt: optimization.converged is true', 0.0_dp)
   call json%get('optimization.steps', steps, found)
   call expect(found, 'check_cli_opt: no optimization.steps')
   call report(nfail, steps > 0, 'check_cli_opt: optimization.steps', real(steps, dp))
   call json%get('optimization.energy_evaluations', evals, found)
   call expect(found, 'check_cli_opt: no optimization.energy_evaluations')
   call report(nfail, evals >= steps, 'check_cli_opt: energy_evaluations - steps', real(evals - steps, dp))
   call json%get('optimization.final_energy', final_energy, found)
   call expect(found, 'check_cli_opt: no optimization.final_energy')
   call report(nfail, final_energy < -70.0_dp .and. final_energy > -80.0_dp, &
              'check_cli_opt: optimization.final_energy is a plausible water RHF/STO-3G number, Eh', final_energy)
   call json%get('optimization.final_gradient_max', gmax, found)
   call expect(found, 'check_cli_opt: no optimization.final_gradient_max')
   call json%get('optimization.final_gradient_rms', grms, found)
   call expect(found, 'check_cli_opt: no optimization.final_gradient_rms')
   call report(nfail, gmax >= grms, 'check_cli_opt: final_gradient_max - final_gradient_rms, Eh/Bohr', gmax - grms)
   call json%get('optimization.coordinates', coord_str, found)
   call expect(found, 'check_cli_opt: no optimization.coordinates')
   call report(nfail, coord_str == 'cartesian', 'check_cli_opt: optimization.coordinates is the default (cartesian)', &
              0.0_dp)
   call json%get('optimization.trajectory_energies', traj_e, found)
   call expect(found, 'check_cli_opt: no optimization.trajectory_energies')
   call report(nfail, size(traj_e) == steps, 'check_cli_opt: size(trajectory_energies) - steps', &
              real(size(traj_e) - steps, dp))
   call json%get('optimization.final_geometry_angstrom(1)', row, found)
   call expect(found .and. size(row) == 3, 'check_cli_opt: no final_geometry_angstrom(1)')
   natm_json = 0
   do
      call json%get('optimization.final_geometry_angstrom('//to_string(natm_json + 1)//')', row, found)
      if (.not. found) exit
      natm_json = natm_json + 1
   end do
   call report(nfail, natm_json == 3, 'check_cli_opt: final_geometry_angstrom rows (water)', real(natm_json, dp))
   call json%destroy()

   ! --- the xyz files themselves -------------------------------------------
   inquire (file=work_dir//'/water_opt.xyz', exist=exists_opt)
   call report(nfail, exists_opt, 'check_cli_opt: water_opt.xyz was written', 0.0_dp)
   inquire (file=work_dir//'/water_traj.xyz', exist=exists_traj)
   call report(nfail, exists_traj, 'check_cli_opt: water_traj.xyz was written', 0.0_dp)
   if (exists_opt) then
      open (newunit=u, file=work_dir//'/water_opt.xyz', status='old', action='read', iostat=ios)
      read (u, '(a)', iostat=ios) line
      call report(nfail, ios == 0 .and. trim(line) == '3', 'check_cli_opt: water_opt.xyz starts with atom count 3', &
                 0.0_dp)
      read (u, '(a)', iostat=ios) line
      call report(nfail, ios == 0 .and. index(line, 'converged') == 1, &
                 'check_cli_opt: water_opt.xyz''s comment line starts with converged/NOT converged', 0.0_dp)
      close (u)
   end if
   if (exists_traj) then
      open (newunit=u, file=work_dir//'/water_traj.xyz', status='old', action='read', iostat=ios)
      block
         integer :: nlines, frames
         nlines = 0
         do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            nlines = nlines + 1
         end do
         frames = nlines/5   !! each frame is natm(1) + comment(1) + 3 atom lines
         call report(nfail, frames == steps, 'check_cli_opt: traj.xyz frames - steps', real(frames - steps, dp))
      end block
      close (u)
   end if

   if (nfail > 0) then
      print '(a,i0,a)', 'check_cli_opt: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_cli_opt: RESULT: PASS'

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   subroutine report(nf, cond, what, value)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: value
      if (cond) then
         print '(a,a,1x,es10.2)', '  PASS  ', what, value
      else
         print '(a,a,1x,es10.2)', '  FAIL  ', what, value
         nf = nf + 1
      end if
   end subroutine report

end program check_cli_opt
