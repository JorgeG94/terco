!
! terco's --opt against an independent optimisation: pyscf's own RHF/STO-3G
! analytic gradient driving scipy's L-BFGS-B (test/pyscf_refs/generate.py's
! gen_opt, committed as opt_water.ref/opt_hcn.ref), from the SAME start
! geometries as /home/jorge/dev/metalquicha/DLFIND_INTEGRATION_NOTE.md
! section 6 (water.xyz, hcn.xyz). Different language, different SCF code,
! different optimiser (scipy's Cartesian L-BFGS-B vs DL-FIND's own, here in
! Cartesian too) -- agreement is independent validation of trc_dlfind_bridge
! end to end, not a restatement of the same algorithm.
!
! Both references optimised to max|g| < 1e-6 Eh/Bohr (actually ~1e-9,
! generate.py's own log); this CLI run uses --opt-gtol 1e-5, two orders
! looser, so the comparison's error is dominated by terco's OWN
! convergence slack, not pyscf's.
!
! TOLERANCES, both derived from --opt-gtol via a harmonic estimate at the
! geometry (deltaE ~ g^2/(2k), deltaR ~ g/k) with k ~ 0.3-0.6 Eh/Bohr^2 for
! a light-atom bond stretch (a generous range: O-H/C-H/C-N/C-triple-N all
! fall in it) and g = 1e-5 Eh/Bohr:
!   ENERGY_TOL = 5e-8 Eh    (g^2/(2*0.3) ~ 1.7e-10, ~300x headroom for the
!                            optimiser's own slack below its own gtol)
!   BOND_TOL   = 3e-4 Bohr  (g/0.4 ~ 2.5e-5, ~12x headroom)
!   ANGLE_TOL  = 0.05 deg   (water's H-O-H; HCN's is ~180 deg, checked the
!                            same way -- a linear molecule relaxing to a
!                            slightly bent one would show up here too)
! Measured water/hcn agreement at the time this was written: energy to
! ~5e-11 Eh, O-H/H-C/C-N to <1e-5 Bohr (test/pyscf_refs/opt_water.ref,
! opt_hcn.ref).
!
program check_opt_vs_pyscf
   use, intrinsic :: iso_c_binding, only: c_double
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use trc_constants, only: ANGSTROM_TO_BOHR
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   implicit none

   character(len=:), allocatable :: cli_exe, xyzfile, ref_base, probe_json, probe_out, cmd
   character(len=1024) :: arg
   integer :: rc, nfail, nat, k
   real(dp), allocatable :: z_ref(:), coords_ref(:, :), coords_json(:, :)
   real(dp) :: e_ref, gmax_ref, e_json
   logical :: converged_json, found
   type(json_file) :: json
   real(dp), allocatable :: row(:)
   real(dp) :: d1_ref, d2_ref, d3_ref, ang_ref
   real(dp) :: d1_json, d2_json, d3_json, ang_json
   real(dp), parameter :: ENERGY_TOL = 5.0e-8_dp
   real(dp), parameter :: BOND_TOL = 3.0e-4_dp
   real(dp), parameter :: ANGLE_TOL = 0.05_dp

   nfail = 0
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); xyzfile = trim(arg)
   call get_command_argument(3, arg); ref_base = trim(arg)
   call get_command_argument(4, arg); probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(xyzfile) == 0 .or. len(ref_base) == 0 .or. len(probe_json) == 0) then
      print '(a)', 'check_opt_vs_pyscf: usage: check_opt_vs_pyscf <terco CLI> <xyz> <ref base> <scratch .json>'
      stop 1
   end if
   probe_out = probe_json//'.out'

   cmd = trim(cli_exe)//' '//xyzfile//' -b basis_sets/sto-3g.json --opt --opt-gtol 1e-5 --silent --json '// &
         probe_json//' > '//probe_out//' 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   if (rc /= 0) then
      print '(a,i0)', 'check_opt_vs_pyscf: the CLI run itself failed, exit code ', rc; stop 1
   end if

   call json%initialize()
   call json%load(filename=probe_json)
   if (json%failed()) then
      print '(a)', 'check_opt_vs_pyscf: cannot read '//probe_json; stop 1
   end if
   call json%get('optimization.converged', converged_json, found)
   call expect(found, 'check_opt_vs_pyscf: no optimization.converged in the CLI''s JSON output')
   call expect(converged_json, 'check_opt_vs_pyscf: the CLI run did not converge at --opt-gtol 1e-5')
   call json%get('optimization.final_energy', e_json, found)
   call expect(found, 'check_opt_vs_pyscf: no optimization.final_energy in the CLI''s JSON output')

   ! Read NAT from the reference FIRST, so the geometry array below is
   ! sized correctly before either set of coordinates is read.
   block
      integer :: uu
      call ref_require(ref_base//'.ref', uu)
      call ref_read_int(uu, 'NAT', nat)
      allocate (z_ref(nat), coords_ref(3, nat))
      call ref_read_vector(uu, 'Z', nat, z_ref)
      call ref_read_matrix(uu, 'COORDS', 3, nat, coords_ref)
      call ref_read_real(uu, 'E_FINAL', e_ref)
      call ref_read_real(uu, 'GMAX_FINAL', gmax_ref)
      call ref_close(uu)
   end block
   call expect(nat == 3, 'check_opt_vs_pyscf: this check only knows water/hcn (nat=3)')

   allocate (coords_json(3, nat))
   do k = 1, nat
      call json%get('optimization.final_geometry_angstrom('//to_string(k)//')', row, found)
      if (.not. found .or. size(row) /= 3) then
         print '(a,i0,a)', 'check_opt_vs_pyscf: no final_geometry_angstrom row ', k, ' in the CLI''s JSON output'
         stop 1
      end if
      coords_json(:, k) = row*ANGSTROM_TO_BOHR   !! back to Bohr, the .ref's own unit
   end do
   call json%destroy()

   call ref_report(nfail, abs(e_json - e_ref) < ENERGY_TOL, &
                   'check_opt_vs_pyscf ('//ref_base//'): final energy vs pyscf+scipy, Eh', abs(e_json - e_ref))

   call two_shortest_and_angle(coords_ref, d1_ref, d2_ref, d3_ref, ang_ref)
   call two_shortest_and_angle(coords_json, d1_json, d2_json, d3_json, ang_json)
   call ref_report(nfail, abs(d1_json - d1_ref) < BOND_TOL, &
                   'check_opt_vs_pyscf ('//ref_base//'): shorter bond vs pyscf+scipy, Bohr', abs(d1_json - d1_ref))
   call ref_report(nfail, abs(d2_json - d2_ref) < BOND_TOL, &
                   'check_opt_vs_pyscf ('//ref_base//'): longer bond vs pyscf+scipy, Bohr', abs(d2_json - d2_ref))
   call ref_report(nfail, abs(d3_json - d3_ref) < BOND_TOL, &
                   'check_opt_vs_pyscf ('//ref_base//'): third distance vs pyscf+scipy, Bohr', abs(d3_json - d3_ref))
   if (ang_ref > 0.0_dp) then
      call ref_report(nfail, abs(ang_json - ang_ref) < ANGLE_TOL, &
                     'check_opt_vs_pyscf ('//ref_base//'): bond angle vs pyscf+scipy, deg', abs(ang_json - ang_ref))
   end if

   if (nfail > 0) then
      print '(a,i0,a)', 'check_opt_vs_pyscf ('//ref_base//'): ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_opt_vs_pyscf ('//ref_base//'): RESULT: PASS'

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   !
   ! For a 3-atom molecule only (water, hcn): the two SHORTEST pairwise
   ! distances (the two bonds, whichever atom is their shared vertex),
   ! the remaining (third) distance, and the angle at the shared vertex
   ! -- 0 if the three atoms have no shared vertex (never happens for 3
   ! atoms: one pair is always the longest and the other two share a
   ! vertex), kept as a dummy guard rather than asserted away.
   !
   subroutine two_shortest_and_angle(xyz, d1, d2, d3, angle_deg)
      real(dp), intent(in) :: xyz(3, 3)
      real(dp), intent(out) :: d1, d2, d3, angle_deg
      real(dp) :: d(3), v1(3), v2(3)
      integer :: apex
      d(1) = norm2(xyz(:, 1) - xyz(:, 2))   !! 1-2
      d(2) = norm2(xyz(:, 1) - xyz(:, 3))   !! 1-3
      d(3) = norm2(xyz(:, 2) - xyz(:, 3))   !! 2-3
      ! The longest distance is the one NOT touching the shared vertex;
      ! the apex is whichever atom is common to the other two pairs.
      if (d(1) >= d(2) .and. d(1) >= d(3)) then
         apex = 3; d3 = d(1); d1 = min(d(2), d(3)); d2 = max(d(2), d(3))
      else if (d(2) >= d(1) .and. d(2) >= d(3)) then
         apex = 2; d3 = d(2); d1 = min(d(1), d(3)); d2 = max(d(1), d(3))
      else
         apex = 1; d3 = d(3); d1 = min(d(1), d(2)); d2 = max(d(1), d(2))
      end if
      v1 = 0.0_dp; v2 = 0.0_dp
      select case (apex)
      case (1)
         v1 = xyz(:, 2) - xyz(:, 1); v2 = xyz(:, 3) - xyz(:, 1)
      case (2)
         v1 = xyz(:, 1) - xyz(:, 2); v2 = xyz(:, 3) - xyz(:, 2)
      case (3)
         v1 = xyz(:, 1) - xyz(:, 3); v2 = xyz(:, 2) - xyz(:, 3)
      end select
      if (norm2(v1) > 0.0_dp .and. norm2(v2) > 0.0_dp) then
         angle_deg = acos(max(-1.0_dp, min(1.0_dp, dot_product(v1, v2)/(norm2(v1)*norm2(v2)))))*180.0_dp/PI_LOCAL()
      else
         angle_deg = 0.0_dp
      end if
   end subroutine two_shortest_and_angle

   pure function PI_LOCAL() result(p)
      real(dp) :: p
      p = 4.0_dp*atan(1.0_dp)
   end function PI_LOCAL

end program check_opt_vs_pyscf
