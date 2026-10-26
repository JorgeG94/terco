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
! A test-drive suite.
!
module check_opt_vs_mqc_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use pyscf_ref_io, only: ref_report
   implicit none
   private
   public :: collect_opt_vs_mqc, cli_exe, probe_json, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: probe_json
      !! Where the CLI writes its --json, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz outputs (--opt
      !! writes <geometry>_opt.xyz/_traj.xyz into the CLI's own CWD,
      !! app/main.F90, so each case runs with its own `cd work_dir &&`).
   real(dp), parameter :: ENERGY_TOL = 2.0e-6_dp
      !! mqc's own ~2e-8 Eh basis difference, plus generous headroom for
      !! two different optimisers both sitting somewhere under a loose
      !! 4.5e-4 Eh/Bohr gradient tolerance rather than at a shared point.
   real(dp), parameter :: BOND_TOL = 2.0e-4_dp
      !! Angstrom.

contains

   subroutine collect_opt_vs_mqc(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt vs mqc's own reference, water", test_water), &
         new_unittest("--opt vs mqc's own reference, hcn", test_hcn) &
         ]
   end subroutine collect_opt_vs_mqc

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out
      integer :: rc, nfail
      real(dp) :: e_json, oh
      logical :: converged_json, found
      type(json_file) :: json
      real(dp), allocatable :: row(:)

      nfail = 0
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      probe_out = probe_json//'.water.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt --silent --json '// &
            probe_json//' > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the water CLI run itself failed, exit code "//to_string(rc)); return
      end if
      call json%initialize(); call json%load(filename=probe_json)
      if (json%failed()) then
         call test_failed(error, "cannot read "//probe_json); return
      end if
      call json%get('optimization.converged', converged_json, found)
      if (.not. (found .and. converged_json)) then
         call json%destroy(); call test_failed(error, "the water run did not converge"); return
      end if
      call json%get('optimization.final_energy', e_json, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.final_energy"); return
      end if
      call ref_report(nfail, abs(e_json - (-74.96590121729567_dp)) < ENERGY_TOL, &
                      'check_opt_vs_mqc: water final energy vs mqc, Eh', abs(e_json - (-74.96590121729567_dp)))
      call json%get('optimization.final_geometry_angstrom(1)', row, found)
      if (.not. (found .and. size(row) == 3)) then
         call json%destroy(); call test_failed(error, "no water final_geometry_angstrom(1) (O)"); return
      end if
      block
         real(dp) :: o(3), h1(3)
         o = row
         call json%get('optimization.final_geometry_angstrom(2)', row, found)
         if (.not. (found .and. size(row) == 3)) then
            call json%destroy(); call test_failed(error, "no water final_geometry_angstrom(2) (H)"); return
         end if
         h1 = row
         oh = norm2(o - h1)
      end block
      call json%destroy()
      call ref_report(nfail, abs(oh - 0.98941_dp) < BOND_TOL, 'check_opt_vs_mqc: water O-H vs mqc, Angstrom', &
                      abs(oh - 0.98941_dp))
      call check_stray(nfail, 'water')
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_water

   subroutine test_hcn(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, probe_out
      integer :: rc, nfail
      real(dp) :: e_json
      logical :: converged_json, found
      type(json_file) :: json

      nfail = 0
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn.xyz '//work_dir//'/hcn.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      probe_out = probe_json//'.hcn.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn.xyz -b basis_sets/sto-3g.json --opt --silent --json '// &
            probe_json//' > '//probe_out//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the hcn CLI run itself failed, exit code "//to_string(rc)); return
      end if
      call json%initialize(); call json%load(filename=probe_json)
      if (json%failed()) then
         call test_failed(error, "cannot read "//probe_json); return
      end if
      call json%get('optimization.converged', converged_json, found)
      if (.not. (found .and. converged_json)) then
         call json%destroy(); call test_failed(error, "the hcn run did not converge"); return
      end if
      call json%get('optimization.final_energy', e_json, found)
      call json%destroy()
      if (.not. found) then
         call test_failed(error, "no optimization.final_energy (hcn)"); return
      end if
      call ref_report(nfail, abs(e_json - (-91.67520899058249_dp)) < ENERGY_TOL, &
                      'check_opt_vs_mqc: hcn final energy vs mqc, Eh', abs(e_json - (-91.67520899058249_dp)))
      call check_stray(nfail, 'hcn')
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_hcn

   !
   ! <base>_opt.xyz/<base>_traj.xyz must not have landed in test/ itself
   ! (this program's own WORKING_DIRECTORY, ctest's convention) -- the
   ! whole point of each case's own `cd work_dir &&` above.
   !
   subroutine check_stray(nfail, base)
      integer, intent(inout) :: nfail
      character(len=*), intent(in) :: base
      logical :: exists
      inquire (file=base//'_opt.xyz', exist=exists)
      call ref_report(nfail, .not. exists, 'check_opt_vs_mqc: '//base//'_opt.xyz was not left in test/ itself', &
                      merge(1.0_dp, 0.0_dp, exists))
      inquire (file=base//'_traj.xyz', exist=exists)
      call ref_report(nfail, .not. exists, 'check_opt_vs_mqc: '//base//'_traj.xyz was not left in test/ itself', &
                      merge(1.0_dp, 0.0_dp, exists))
   end subroutine check_stray

end module check_opt_vs_mqc_suite

program check_opt_vs_mqc
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_vs_mqc_suite, only: collect_opt_vs_mqc, cli_exe, probe_json, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); probe_json = trim(arg)
   call get_command_argument(3, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_vs_mqc: usage: check_opt_vs_mqc <terco CLI binary> <scratch .json path> "// &
         "<scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_vs_mqc, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_vs_mqc
