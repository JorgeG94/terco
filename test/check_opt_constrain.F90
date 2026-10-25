!
! --opt-constrain "angle 2 1 3" (1-based: O is the vertex), water/STO-3G,
! --opt-coords hdlc -- app/trc_dlfind_bridge.F90's own icons block,
! verified directly against dlf_hdlc_constraint.f90%assign_cons (type in
! slot 1, atom indices in slots 2-5). The constrained value is always
! whatever the START geometry has (no target-value slot anywhere in
! DL-FIND's own public interface, see that module's header) -- test/
! water.xyz's own H-O-H angle, held fixed while the two O-H bond lengths
! relax.
!
! CONVERGENCE VERDICT: a constrained optimum's RAW Cartesian gradient
! does not go to zero the way an unconstrained one does -- the
! constraint acts like a Lagrange multiplier, holding a permanent,
! nonzero force along the constrained direction even once every FREE
! internal coordinate's own gradient is exactly zero. Measured directly
! on this exact system with lbfgs, cg-auto and more steps alike (the
! energy converges to 12 decimal places either way, confirming this is
! the TRUE constrained optimum, not a stalled search): the angle
! constraint's own RAW Cartesian residual settles at max|g| ~1.86e-2
! Eh/Bohr, the bond constraint's at ~3.23e-2, both well above any
! --opt-gtol tight enough to let DL-FIND actually DO the six-to-eight
! steps this takes (a loose --opt-gtol scales ALL of DL-FIND's own
! criteria, not just the gradient one, and was measured to make it stop
! after a single step, nowhere near this optimum). This bridge's own
! "CONVERGENCE AUTHORITY" header already documents that it cannot read
! DL-FIND's own INTERNAL-coordinate gradient back -- the number that
! correctly reaches zero here, since HDLC removes the constrained DOF
! from the internal coordinate set entirely, rather than penalising it.
! So this file keeps --opt-gtol TIGHT (GTOL_ARG) and checks DL-FIND's
! OWN observable verdict instead (the "Convergence reached"/"Optimisation
! converged" log lines, --opt-print-level 4) for whether the run itself
! succeeded, the same fallback check_opt_neb.F90 uses for the identical
! reason.
!
! EXTERNAL anchor: test/pyscf_refs/opt_constrained_water_angle.ref,
! written by test/pyscf_refs/generate.py's own gen_opt_opt_constrained
! (sic: gen_opt_constrained) -- a from-scratch scipy L-BFGS-B minimisation
! over exactly the two free O-H bond lengths, with the two O->H bond
! DIRECTIONS taken from water.xyz itself and held fixed throughout, so the
! angle is EXACTLY water.xyz's own. Independent of terco end to end:
! different language, different SCF code, a reduced 2-parameter optimiser
! acting directly on the only physically free coordinates, not DL-FIND's
! internal-coordinate machinery.
!
! Also a bond constraint ("bond 1 2": the O-H1 distance) -- checked only
! for internal consistency (the bond length is held, the other O-H bond
! and the angle both relax), no separate external anchor (a second scipy
! generator for one more constraint type would not exercise anything
! gen_opt_constrained does not already cover for the SAME molecule).
!
! A test-drive suite.
!
module check_opt_constrain_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use trc_constants, only: ANGSTROM_TO_BOHR
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_opt_constrain, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the JSON output.
   real(dp), parameter :: ANGLE_TOL = 1.0e-6_dp         !! deg -- "the final angle equals the start to ~1e-6 deg"
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 2.0e-6_dp !! Eh -- HDLC's own convergence is looser than Cartesian's
   real(dp), parameter :: ANCHOR_BOND_TOL = 2.0e-4_dp   !! Bohr
   character(len=*), parameter :: GTOL_ARG = '--opt-gtol 1e-6'   !! this file's own header, "CONVERGENCE VERDICT"

   ! --- state shared between the angle-constraint run and its checks -----
   real(dp) :: e_final = 0.0_dp
   real(dp), allocatable :: xyz(:, :)

contains

   subroutine collect_opt_constrain(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-constrain 'angle 2 1 3': the H-O-H angle is held, bonds relax symmetrically", &
                      test_angle_constraint), &
         new_unittest("--opt-constrain 'angle 2 1 3': matches the pyscf opt_constrained_water_angle anchor", &
                      test_angle_vs_anchor), &
         new_unittest("--opt-constrain 'bond 1 2': the O-H1 distance is held at the start value", &
                      test_bond_constraint), &
         new_unittest("--opt-constrain: no stray xyz output in test/ itself", test_no_stray) &
         ]
   end subroutine collect_opt_constrain

   subroutine test_angle_constraint(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc, nfail
      logical :: found
      type(json_file) :: json
      real(dp) :: start_angle, final_angle, d_oh1, d_oh2

      nfail = 0
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      call start_water_angle(start_angle)

      ! --opt-print-level 4 (not --silent): "Convergence reached" only
      ! prints at that verbosity -- this file's own header,
      ! "CONVERGENCE VERDICT", on why that log line is the pass/fail
      ! gate here, not this bridge's own exit code.
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-coords hdlc --opt-constrain "angle 2 1 3" '//GTOL_ARG//' --opt-maxiter 100 '// &
            '--opt-print-level 4 --json cons.json > cons.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      block
         logical :: conv_reached
         call grep_contains(work_dir//'/cons.out', 'Convergence reached', conv_reached)
         if (.not. conv_reached) then
            call test_failed(error, "DL-FIND's own convergence test did not pass (no 'Convergence "// &
                             "reached' in cons.out), CLI exit code "//to_string(rc)); return
         end if
      end block
      call json%initialize()
      call json%load(filename=work_dir//'/cons.json')
      if (json%failed()) then
         call test_failed(error, "cannot read cons.json"); return
      end if
      call json%get('optimization.final_energy', e_final, found)
      if (.not. found) then
         call json%destroy(); call test_failed(error, "no final_energy"); return
      end if
      call read_geometry(json, 3, xyz, found)
      call json%destroy()
      if (.not. found) then
         call test_failed(error, "bad final_geometry_angstrom"); return
      end if

      call angle_and_bonds(xyz, final_angle, d_oh1, d_oh2)
      print '(a,f20.12,a,f16.10,a,f16.10)', '  final energy ', e_final, '   angle (deg) ', final_angle, &
         '   start angle (deg) ', start_angle
      call ref_report(nfail, abs(final_angle - start_angle) < ANGLE_TOL, &
                      'check_opt_constrain: the H-O-H angle is held at the start value', &
                      abs(final_angle - start_angle))
      call ref_report(nfail, abs(d_oh1 - d_oh2) < 1.0e-6_dp, &
                      'check_opt_constrain: the two O-H bonds relax symmetrically (water''s own symmetry)', &
                      abs(d_oh1 - d_oh2))
      ! Constraining necessarily raises (or leaves unchanged) the energy
      ! relative to the unconstrained global minimum.
      block
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
         call ref_report(nfail, e_final >= e_ref - 1.0e-8_dp, &
                         'check_opt_constrain: the constrained energy is not below the unconstrained '// &
                         'minimum', 0.0_dp)
      end block
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_angle_constraint

   !
   ! Depends on test_angle_constraint having already populated e_final/xyz
   ! -- test-drive runs a module's tests in the order collect_opt_constrain
   ! lists them, not in parallel (run_testsuite's own parallel=.false. in
   ! the program below).
   !
   subroutine test_angle_vs_anchor(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: uu, nat_ref, nfail
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      real(dp) :: e_ref, gmax_ref, angle_ref, d1_ref, d2_ref, angle_ref_check, angle_terco, d1_terco, d2_terco

      if (.not. allocated(xyz)) then
         call test_failed(error, "test_angle_constraint did not leave a geometry to check"); return
      end if
      nfail = 0
      call ref_require('opt_constrained_water_angle.ref', uu)
      call ref_read_int(uu, 'NAT', nat_ref)
      if (nat_ref /= 3) then
         call ref_close(uu); call test_failed(error, "opt_constrained_water_angle.ref: NAT must be 3"); return
      end if
      allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
      call ref_read_vector(uu, 'Z', nat_ref, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_read_real(uu, 'ANGLE_DEG', angle_ref)
      call ref_close(uu)
      call angle_and_bonds(coords_ref, angle_ref_check, d1_ref, d2_ref)
      call angle_and_bonds(xyz*ANGSTROM_TO_BOHR, angle_terco, d1_terco, d2_terco)
      call ref_report(nfail, abs(angle_ref_check - angle_ref) < 1.0e-6_dp, &
                      'check_opt_constrain: the pyscf anchor''s own stored angle matches its COORDS, deg', &
                      abs(angle_ref_check - angle_ref))
      call ref_report(nfail, abs(e_final - e_ref) < ANCHOR_ENERGY_TOL, &
                      'check_opt_constrain: constrained energy vs pyscf anchor, Eh', abs(e_final - e_ref))
      call ref_report(nfail, abs(d1_terco - d1_ref) < ANCHOR_BOND_TOL, &
                      'check_opt_constrain: O-H1 bond length vs pyscf anchor, Bohr', abs(d1_terco - d1_ref))
      call ref_report(nfail, abs(d2_terco - d2_ref) < ANCHOR_BOND_TOL, &
                      'check_opt_constrain: O-H2 bond length vs pyscf anchor, Bohr', abs(d2_terco - d2_ref))
      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine test_angle_vs_anchor

   subroutine test_bond_constraint(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd
      integer :: rc
      real(dp) :: start_d1, d1
      real(dp), allocatable :: bond_xyz(:, :)
      logical :: found
      type(json_file) :: json

      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      call start_water_bond(start_d1)
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-coords hdlc --opt-constrain "bond 1 2" '//GTOL_ARG//' --opt-maxiter 100 '// &
            '--opt-print-level 4 --json cons_bond.json > cons_bond.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      block
         logical :: conv_reached
         call grep_contains(work_dir//'/cons_bond.out', 'Convergence reached', conv_reached)
         if (.not. conv_reached) then
            call test_failed(error, "DL-FIND's own convergence test did not pass (no 'Convergence "// &
                             "reached' in cons_bond.out), CLI exit code "//to_string(rc)); return
         end if
      end block
      call json%initialize()
      call json%load(filename=work_dir//'/cons_bond.json')
      if (json%failed()) then
         call test_failed(error, "cannot read cons_bond.json"); return
      end if
      call read_geometry(json, 3, bond_xyz, found)
      call json%destroy()
      if (.not. found) then
         call test_failed(error, "bad final_geometry_angstrom"); return
      end if
      d1 = norm2(bond_xyz(:, 1) - bond_xyz(:, 2))*ANGSTROM_TO_BOHR
      if (abs(d1 - start_d1) >= 1.0e-5_dp) then
         call test_failed(error, "'bond 1 2' did not hold the O-H1 distance at the start value, |diff|="// &
                          to_string(abs(d1 - start_d1))//" Bohr")
      end if
   end subroutine test_bond_constraint

   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: exists
      inquire (file='water_opt.xyz', exist=exists)
      if (exists) call test_failed(error, "water_opt.xyz was left in test/ itself")
   end subroutine test_no_stray

   subroutine start_water_angle(angle_deg)
      real(dp), intent(out) :: angle_deg
      real(dp) :: g(3, 3), d1, d2
      g(:, 1) = [0.0_dp, 0.0_dp, 0.1173_dp]
      g(:, 2) = [0.0_dp, 0.7572_dp, -0.4692_dp]
      g(:, 3) = [0.0_dp, -0.7572_dp, -0.4692_dp]
      call angle_and_bonds(g, angle_deg, d1, d2)
   end subroutine start_water_angle

   subroutine start_water_bond(d_bohr)
      real(dp), intent(out) :: d_bohr
      real(dp) :: g(3, 3), ang, d2
      g(:, 1) = [0.0_dp, 0.0_dp, 0.1173_dp]
      g(:, 2) = [0.0_dp, 0.7572_dp, -0.4692_dp]
      g(:, 3) = [0.0_dp, -0.7572_dp, -0.4692_dp]
      call angle_and_bonds(g*ANGSTROM_TO_BOHR, ang, d_bohr, d2)
   end subroutine start_water_bond

   subroutine angle_and_bonds(xyz_in, angle_deg, d_oh1, d_oh2)
      !! xyz_in(:,1)=O, (:,2)=H1, (:,3)=H2, any consistent unit -- returns
      !! the angle in degrees and the two O-H distances in xyz_in's OWN
      !! unit (the caller converts).
      real(dp), intent(in) :: xyz_in(3, 3)
      real(dp), intent(out) :: angle_deg, d_oh1, d_oh2
      real(dp) :: v1(3), v2(3)
      v1 = xyz_in(:, 2) - xyz_in(:, 1)
      v2 = xyz_in(:, 3) - xyz_in(:, 1)
      d_oh1 = norm2(v1)
      d_oh2 = norm2(v2)
      angle_deg = acos(max(-1.0_dp, min(1.0_dp, dot_product(v1, v2)/(d_oh1*d_oh2))))*180.0_dp/ &
                 (4.0_dp*atan(1.0_dp))
   end subroutine angle_and_bonds

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

end module check_opt_constrain_suite

program check_opt_constrain
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_constrain_suite, only: collect_opt_constrain, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_constrain: usage: check_opt_constrain <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_constrain, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_constrain
