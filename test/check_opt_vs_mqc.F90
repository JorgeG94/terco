!
! A LOOSE sanity check against metalquicha's own reference optimisations
! (water and HCN, RHF/STO-3G, Cartesian, DL-FIND's own engine default gtol
! 4.5e-4 Eh/Bohr -- terco's own CLI default, so no --opt-gtol override is
! passed here, unlike check_opt_vs_pyscf's tight comparison).
!
! mqc's own STO-3G basis file differs from pyscf's (and so from terco's
! test/basis_sets/sto-3g.json, generated from pyscf's table) by ~2e-8 Eh
! at this geometry, measured when mqc's reference runs were made. This
! check's tolerances carry that difference plus the two optimisers'
! independent convergence paths
! at a loose gtol; check_opt_vs_pyscf.F90, against a pyscf+scipy
! reference generated through terco's OWN basis file at a tight gtol, is
! the real anchor. This one exists only to show terco's --opt lands on
! the same answer mqc's own run (and so, by that run's own validation,
! DL-FIND itself) does.
!
program check_opt_vs_mqc
   use pic_types, only: dp
   use json_module, only: json_file
   implicit none

   character(len=:), allocatable :: cli_exe, probe_json, probe_out, cmd, work_dir
   character(len=1024) :: arg
   integer :: rc, nfail
   real(dp) :: e_json
   logical :: converged_json, found
   type(json_file) :: json
   real(dp), allocatable :: row(:)
   real(dp) :: oh
   real(dp), parameter :: ENERGY_TOL = 2.0e-6_dp
      !! mqc's own ~2e-8 Eh basis difference, plus generous headroom for
      !! two different optimisers both sitting somewhere under a loose
      !! 4.5e-4 Eh/Bohr gradient tolerance rather than at a shared point.
   real(dp), parameter :: BOND_TOL = 2.0e-4_dp   !! Angstrom

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   call get_command_argument(3, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_opt_vs_mqc: usage: check_opt_vs_mqc <terco CLI binary> <scratch .json path> <scratch dir>'
      stop 1
   end if
   ! Run with WORKING_DIRECTORY set to work_dir (not test/ itself, this
   ! program's own WORKING_DIRECTORY, ctest's convention) -- --opt writes
   ! <geometry>_opt.xyz/_traj.xyz into the CLI's CWD (app/main.F90);
   ! check_opt_connect.F90's own `cd work_dir &&` pattern.
   call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f hcn.xyz '//work_dir//'/hcn.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)
   probe_out = probe_json//'.water.out'

   ! --- water -------------------------------------------------------------
   cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt --silent --json '// &
         probe_json//' > '//probe_out//' 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   call expect(rc == 0, 'check_opt_vs_mqc: the water CLI run itself failed')
   call json%initialize(); call json%load(filename=probe_json)
   call expect(.not. json%failed(), 'check_opt_vs_mqc: cannot read '//probe_json)
   call json%get('optimization.converged', converged_json, found)
   call expect(found .and. converged_json, 'check_opt_vs_mqc: water run did not converge')
   call json%get('optimization.final_energy', e_json, found)
   call expect(found, 'check_opt_vs_mqc: no optimization.final_energy')
   call report(nfail, abs(e_json - (-74.96590121729567_dp)) < ENERGY_TOL, &
              'check_opt_vs_mqc: water final energy vs mqc, Eh', abs(e_json - (-74.96590121729567_dp)))
   call json%get('optimization.final_geometry_angstrom(1)', row, found)
   call expect(found .and. size(row) == 3, 'check_opt_vs_mqc: no water final_geometry_angstrom(1) (O)')
   block
      real(dp) :: o(3), h1(3)
      o = row
      call json%get('optimization.final_geometry_angstrom(2)', row, found)
      call expect(found .and. size(row) == 3, 'check_opt_vs_mqc: no water final_geometry_angstrom(2) (H)')
      h1 = row
      oh = norm2(o - h1)
   end block
   call json%destroy()
   call report(nfail, abs(oh - 0.98941_dp) < BOND_TOL, 'check_opt_vs_mqc: water O-H vs mqc, Angstrom', &
              abs(oh - 0.98941_dp))

   ! --- hcn -----------------------------------------------------------------
   probe_out = probe_json//'.hcn.out'
   cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn.xyz -b basis_sets/sto-3g.json --opt --silent --json '// &
         probe_json//' > '//probe_out//' 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   call expect(rc == 0, 'check_opt_vs_mqc: the hcn CLI run itself failed')
   call json%initialize(); call json%load(filename=probe_json)
   call expect(.not. json%failed(), 'check_opt_vs_mqc: cannot read '//probe_json)
   call json%get('optimization.converged', converged_json, found)
   call expect(found .and. converged_json, 'check_opt_vs_mqc: hcn run did not converge')
   call json%get('optimization.final_energy', e_json, found)
   call expect(found, 'check_opt_vs_mqc: no optimization.final_energy (hcn)')
   call json%destroy()
   call report(nfail, abs(e_json - (-91.67520899058249_dp)) < ENERGY_TOL, &
              'check_opt_vs_mqc: hcn final energy vs mqc, Eh', abs(e_json - (-91.67520899058249_dp)))

   ! water_opt.xyz/water_traj.xyz/hcn_opt.xyz/hcn_traj.xyz must not have
   ! landed in test/ itself (this program's own WORKING_DIRECTORY, ctest's
   ! convention) -- the whole point of the `cd work_dir &&` above.
   block
      character(len=*), parameter :: stray(4) = &
         [character(len=14) :: 'water_opt.xyz', 'water_traj.xyz', 'hcn_opt.xyz', 'hcn_traj.xyz']
      integer :: k
      logical :: exists
      do k = 1, size(stray)
         inquire (file=trim(stray(k)), exist=exists)
         call report(nfail, .not. exists, 'check_opt_vs_mqc: '//trim(stray(k))//' was not left in test/ itself', &
                    merge(1.0_dp, 0.0_dp, exists))
      end do
   end block

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_vs_mqc: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_vs_mqc: RESULT: PASS'

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   subroutine report(nf, cond, what, diff)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: diff
      if (cond) then
         print '(a,a,es10.2)', '  PASS  ', what, diff
      else
         print '(a,a,es10.2)', '  FAIL  ', what, diff
         nf = nf + 1
      end if
   end subroutine report

end program check_opt_vs_mqc
