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
! Also a bond constraint ("bond 1 2": the O-H1 distance), cheap to add and
! asked for by this PR's brief -- checked only for internal consistency
! (the bond length is held, the other O-H bond and the angle both relax),
! no separate external anchor (a second scipy generator for one more
! constraint type would not exercise anything gen_opt_constrained does
! not already cover for the SAME molecule).
!
program check_opt_constrain
   use pic_types, only: dp
   use json_module, only: json_file
   use trc_constants, only: ANGSTROM_TO_BOHR
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   implicit none

   character(len=:), allocatable :: cli_exe, work_dir, cmd
   character(len=1024) :: arg
   integer :: rc, nfail
   logical :: found
   real(dp) :: e_final
   real(dp), allocatable :: xyz(:, :)
   type(json_file) :: json
   real(dp), parameter :: ANGLE_TOL = 1.0e-6_dp         !! deg -- "the final angle equals the start to ~1e-6 deg"
   real(dp), parameter :: ANCHOR_ENERGY_TOL = 2.0e-6_dp !! Eh -- HDLC's own convergence is looser than Cartesian's
   real(dp), parameter :: ANCHOR_BOND_TOL = 2.0e-4_dp   !! Bohr

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      print '(a)', 'check_opt_constrain: usage: check_opt_constrain <terco CLI> <scratch dir>'
      stop 1
   end if
   call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f water.xyz '//work_dir//'/water.xyz', exitstat=rc, wait=.true.)
   call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                             exitstat=rc, wait=.true.)

   block
      real(dp) :: start_angle
      call start_water_angle(start_angle)

      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-coords hdlc --opt-constrain "angle 2 1 3" --opt-gtol 1e-7 --opt-maxiter 100 '// &
            '--silent --json cons.json > cons.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call expect(rc == 0, 'check_opt_constrain: the CLI run itself failed')
      call json%initialize()
      call json%load(filename=work_dir//'/cons.json')
      call expect(.not. json%failed(), 'check_opt_constrain: cannot read cons.json')
      call json%get('optimization.final_energy', e_final, found)
      call expect(found, 'check_opt_constrain: no final_energy')
      call read_geometry(json, 3, xyz)
      call json%destroy()

      block
         real(dp) :: final_angle, d_oh1, d_oh2
         call angle_and_bonds(xyz, final_angle, d_oh1, d_oh2)
         print '(a,f20.12,a,f16.10,a,f16.10)', '  final energy ', e_final, '   angle (deg) ', final_angle, &
            '   start angle (deg) ', start_angle
         call report(nfail, abs(final_angle - start_angle) < ANGLE_TOL, &
                    'check_opt_constrain: the H-O-H angle is held at the start value', &
                    abs(final_angle - start_angle))
         call report(nfail, abs(d_oh1 - d_oh2) < 1.0e-6_dp, &
                    'check_opt_constrain: the two O-H bonds relax symmetrically (water''s own symmetry)', &
                    abs(d_oh1 - d_oh2))
         ! Constraining necessarily raises (or leaves unchanged) the
         ! energy relative to the unconstrained global minimum.
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
            call report(nfail, e_final >= e_ref - 1.0e-8_dp, &
                       'check_opt_constrain: the constrained energy is not below the unconstrained minimum')
         end block
      end block

      ! --- the external anchor: test/pyscf_refs/opt_constrained_water_angle.ref
      block
         integer :: uu, nat_ref
         real(dp), allocatable :: z_ref(:), coords_ref(:, :)
         real(dp) :: e_ref, angle_ref, d1_ref, d2_ref, angle_ref_check, angle_terco, d1_terco, d2_terco
         call ref_require('opt_constrained_water_angle.ref', uu)
         call ref_read_int(uu, 'NAT', nat_ref)
         call expect(nat_ref == 3, 'check_opt_constrain: opt_constrained_water_angle.ref: NAT must be 3')
         allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
         call ref_read_vector(uu, 'Z', nat_ref, z_ref)
         call ref_read_matrix(uu, 'COORDS', 3, nat_ref, coords_ref)
         call ref_read_real(uu, 'E_FINAL', e_ref)
         block
            real(dp) :: gmax_ref
            call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
         end block
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
      end block
   end block

   ! --- a bond constraint: internal-consistency only (no anchor) ---------
   block
      real(dp) :: start_d1
      call start_water_bond(start_d1)
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' water.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-coords hdlc --opt-constrain "bond 1 2" --opt-gtol 1e-6 --opt-maxiter 100 '// &
            '--silent --json cons_bond.json > cons_bond.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      call expect(rc == 0, 'check_opt_constrain: the bond-constraint CLI run itself failed')
      call json%initialize()
      call json%load(filename=work_dir//'/cons_bond.json')
      call expect(.not. json%failed(), 'check_opt_constrain: cannot read cons_bond.json')
      call read_geometry(json, 3, xyz)
      call json%destroy()
      block
         real(dp) :: d1
         d1 = norm2(xyz(:, 1) - xyz(:, 2))*ANGSTROM_TO_BOHR
         call report(nfail, abs(d1 - start_d1) < 1.0e-5_dp, &
                    'check_opt_constrain: "bond 1 2" holds the O-H1 distance at the start value', &
                    abs(d1 - start_d1))
      end block
   end block

   block
      logical :: exists
      inquire (file='water_opt.xyz', exist=exists)
      call report(nfail, .not. exists, 'check_opt_constrain: water_opt.xyz was not left in test/ itself')
   end block

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_constrain: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_constrain: RESULT: PASS'

contains

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

   subroutine read_geometry(j, natm, xyz_out)
      type(json_file), intent(inout) :: j
      integer, intent(in) :: natm
      real(dp), allocatable, intent(out) :: xyz_out(:, :)
      real(dp), allocatable :: row(:)
      logical :: ok
      integer :: ia
      character(len=16) :: buf
      allocate (xyz_out(3, natm))
      do ia = 1, natm
         write (buf, '(i0)') ia
         call j%get('optimization.final_geometry_angstrom('//trim(buf)//')', row, ok)
         call expect(ok .and. size(row) == 3, 'check_opt_constrain: bad final_geometry_angstrom row')
         xyz_out(:, ia) = row
      end do
   end subroutine read_geometry

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   subroutine report(nf, cond, what, val)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in), optional :: val
      if (cond) then
         print '(a,a)', '  PASS  ', what
      else
         if (present(val)) then
            print '(a,a,a,es12.4)', '  FAIL  ', what, '   value=', val
         else
            print '(a,a)', '  FAIL  ', what
         end if
         nf = nf + 1
      end if
   end subroutine report

end program check_opt_constrain
