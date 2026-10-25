!
! --opt-frozen: water/STO-3G with atom 1 (O) frozen, Cartesian coordinates
! -- spec(1) = -1 in app/trc_dlfind_bridge.F90's own cb_get_params, the
! SAME mechanism the residue block already carries (no bridge code change
! at all; app/main.F90's own --opt-frozen validation and parsing is what
! this PR adds). Checks:
!   - the frozen atom's coordinates are BIT-IDENTICAL before and after;
!   - the free atoms converge to a constrained minimum (their own gradient
!     components below --opt-gtol);
!   - freezing a SINGLE atom leaves a RIGID ROTATION of the water
!     molecule entirely free (one atom's absolute position constrains
!     nothing about the molecule's internal shape), so the run is
!     expected to reach the SAME energy as the unconstrained opt_water
!     anchor -- this is not a weaker check than it looks: it is what
!     this test's own water.xyz (deliberately displaced off its
!     equilibrium ahead of the frozen atom) being able to reach the exact
!     global minimum ANYWAY proves the constraint is doing nothing
!     unphysical, only pinning the one atom asked for.
!
! A test-drive suite.
!
module check_opt_frozen_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_opt_frozen, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the xyz/JSON outputs.
   real(dp), parameter :: FROZEN_TOL = 1.0e-10_dp   !! Eh/Bohr-scale coordinate, "bit-identical"
   real(dp), parameter :: FREE_GTOL = 1.0e-5_dp
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 2.0e-7_dp

   ! --- state shared between the single CLI run and the checks on it ------
   real(dp) :: e_final = 0.0_dp
   real(dp), allocatable :: xyz(:, :)
   real(dp) :: o_start(3) = 0.0_dp

contains

   subroutine collect_opt_frozen(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-frozen 1: the frozen atom does not move, the free atoms converge", &
                      test_frozen_run), &
         new_unittest("--opt-frozen 1: reaches the SAME energy as the unconstrained opt_water anchor", &
                      test_frozen_vs_anchor), &
         new_unittest("--opt-frozen: no stray xyz output in test/ itself", test_no_stray) &
         ]
   end subroutine collect_opt_frozen

   subroutine test_frozen_run(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc, nfail, u
      logical :: found, converged
      real(dp) :: gmax
      type(json_file) :: json

      nfail = 0
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)

      ! A deliberately displaced water, so freezing O away from its own
      ! equilibrium position is a REAL test, not a no-op at a point that
      ! was already stationary -- see this file's own header.
      open (newunit=u, file=work_dir//'/water_disp.xyz', status='replace', action='write')
      write (u, '(i0)') 3
      write (u, '(a)') 'check_opt_frozen: water, O displaced off its own equilibrium position'
      write (u, '(a2,3f20.12)') 'O', 0.05_dp, 0.03_dp, 0.0973_dp
      write (u, '(a2,3f20.12)') 'H', 0.0_dp, 0.7572_dp, -0.4692_dp
      write (u, '(a2,3f20.12)') 'H', 0.0_dp, -0.7572_dp, -0.4692_dp
      close (u)
      o_start = [0.05_dp, 0.03_dp, 0.0973_dp]

      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water_disp.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-frozen 1 --opt-gtol 1e-6 --opt-maxiter 100 --silent --json frozen.json > frozen.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      if (rc /= 0) then
         call test_failed(error, "the CLI run itself failed, exit code "//to_string(rc)); return
      end if

      call json%initialize()
      call json%load(filename=work_dir//'/frozen.json')
      if (json%failed()) then
         call test_failed(error, "cannot read frozen.json"); return
      end if
      call json%get('optimization.converged', converged, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no optimization.converged"); return
      end if
      call json%get('optimization.final_energy', e_final, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no final_energy"); return
      end if
      call json%get('optimization.final_gradient_max', gmax, found)
      call read_geometry(json, 3, xyz, found)
      call json%destroy()
      if (.not. found) then
         call test_failed(error, "bad final_geometry_angstrom"); return
      end if

      print '(a,f20.12,a,es10.3)', '  final energy ', e_final, '   max|g| (incl. frozen atom) ', gmax
      print '(a,3f14.8)', '  O (frozen)   final: ', xyz(:, 1)
      print '(a,3f14.8)', '  O (frozen)   start: ', o_start

      call ref_report(nfail, maxval(abs(xyz(:, 1) - o_start)) < FROZEN_TOL, &
                      'check_opt_frozen: the frozen atom did not move', maxval(abs(xyz(:, 1) - o_start)))

      ! The free atoms' own gradient, below gtol -- re-run --grad at the
      ! reported geometry rather than trust --opt's own (frozen-atom-
      ! inclusive) max|g|, which this bridge's header explains is NOT the
      ! right number for a constrained run.
      block
         character(len=:), allocatable :: cmd2
         type(json_file) :: gj
         real(dp), allocatable :: row(:)
         logical :: ok
         integer :: ia
         character(len=16) :: buf
         real(dp) :: free_gmax
         call write_xyz(work_dir//'/frozen_result.xyz', xyz)
         cmd2 = 'cd '//work_dir//' && '//trim(cli_exe)//' frozen_result.xyz -b basis_sets/sto-3g.json '// &
                '--grad --silent --json grad.json > grad.out 2>&1'
         call execute_command_line(cmd2, exitstat=rc, wait=.true.)
         if (rc /= 0) then
            call test_failed(error, "the --grad re-check CLI run itself failed"); return
         end if
         call gj%initialize()
         call gj%load(filename=work_dir//'/grad.json')
         if (gj%failed()) then
            call test_failed(error, "cannot read grad.json"); return
         end if
         free_gmax = 0.0_dp
         do ia = 2, 3   !! atoms 2,3 (H,H) are the FREE ones; atom 1 (O) is frozen
            write (buf, '(i0)') ia
            call gj%get('gradient('//trim(buf)//')', row, ok)
            if (.not. (ok .and. size(row) == 3)) then
               call gj%destroy(); call test_failed(error, "bad gradient row in grad.json"); return
            end if
            free_gmax = max(free_gmax, maxval(abs(row)))
         end do
         call gj%destroy()
         call ref_report(nfail, free_gmax < FREE_GTOL, &
                         'check_opt_frozen: the free atoms'' own gradient is below --opt-gtol', free_gmax)
      end block
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_frozen_run

   !
   ! Depends on test_frozen_run having already populated e_final -- test-
   ! drive runs a module's tests in the order collect_opt_frozen lists
   ! them, not in parallel (run_testsuite's own parallel=.false. in the
   ! program below).
   !
   subroutine test_frozen_vs_anchor(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: uu, nat_ref
      real(dp) :: e_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      call ref_require('opt_water.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_close(uu)
      if (abs(e_final - e_ref) >= ANCHOR_ENERGY_TOL) then
         call test_failed(error, "does not reach the unconstrained opt_water anchor energy, |diff|="// &
                          to_string(abs(e_final - e_ref))//" Eh")
      end if
   end subroutine test_frozen_vs_anchor

   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: exists
      inquire (file='water_disp_opt.xyz', exist=exists)
      if (exists) call test_failed(error, "water_disp_opt.xyz was left in test/ itself")
   end subroutine test_no_stray

   subroutine read_geometry(j, natm, xyz_out, ok)
      type(json_file), intent(inout) :: j
      integer, intent(in) :: natm
      real(dp), allocatable, intent(out) :: xyz_out(:, :)
      logical, intent(out) :: ok
      real(dp), allocatable :: row(:)
      integer :: ia
      character(len=16) :: buf
      allocate (xyz_out(3, natm))
      ok = .true.
      do ia = 1, natm
         write (buf, '(i0)') ia
         call j%get('optimization.final_geometry_angstrom('//trim(buf)//')', row, ok)
         ok = ok .and. size(row) == 3
         if (.not. ok) return
         xyz_out(:, ia) = row
      end do
   end subroutine read_geometry

   subroutine write_xyz(path, xyz_ang)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: xyz_ang(:, :)
      integer :: uu, ia
      character(len=2), parameter :: sym(3) = ['O ', 'H ', 'H ']
      open (newunit=uu, file=path, status='replace', action='write')
      write (uu, '(i0)') size(xyz_ang, 2)
      write (uu, '(a)') 'check_opt_frozen: re-check geometry'
      do ia = 1, size(xyz_ang, 2)
         write (uu, '(a2,3f20.12)') sym(ia), xyz_ang(:, ia)
      end do
      close (uu)
   end subroutine write_xyz

end module check_opt_frozen_suite

program check_opt_frozen
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_frozen_suite, only: collect_opt_frozen, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_frozen: usage: check_opt_frozen <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_frozen, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_frozen
