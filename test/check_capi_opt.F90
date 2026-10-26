!
! The geometry optimiser driven ONLY through the C ABI: terco.h for the
! context (molecule, basis, method -- exactly a C caller like cyclone would
! use) and terco_opt.h for the optimisation itself (trc_opt_options_create/
! trc_opt_set_gtol/trc_opt_run/trc_opt_result_*), never trc_dlfind_bridge's
! Fortran type trc_opt_options_t/trc_opt_result_t directly. This is the one
! check that would not even COMPILE before libterco_opt and
! trc_opt_capi.F90 existed -- the generated bindings it links
! (trc_opt_c_interfaces, include/trc_opt_c_interfaces.f90) did not exist
! either.
!
! Same systems and the same stored references check_opt_vs_pyscf.F90 uses
! for its own water/hcn cases (test/pyscf_refs/opt_water.ref, opt_hcn.ref:
! pyscf's RHF/STO-3G analytic gradient driving scipy's L-BFGS-B, a
! different language, SCF implementation and optimiser entirely), at the
! SAME tolerances -- so this check is independent evidence that the C ABI
! path reaches the exact same answer check_cli_opt/check_opt_vs_pyscf_hcn
! already established for terco's own CLI, not a new numerical claim.
!
! Unlike check_opt_vs_pyscf.F90, there is no terco_cli subprocess and no
! --json/--opt-*.xyz file to inspect: trc_opt_result_geometry/_energy/
! _gradient_max hand the answer straight back in the caller's own buffers,
! which is the whole point of a library entry point over a CLI.
!
! A test-drive suite.
!
module check_capi_opt_suite
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                               TRC_OK
   use trc_opt_c_interfaces, only: trc_opt_options_create, trc_opt_options_destroy, trc_opt_set_gtol, &
                                   trc_opt_result_create, trc_opt_result_destroy, trc_opt_run, &
                                   trc_opt_result_converged, trc_opt_result_energy, trc_opt_result_gradient_max, &
                                   trc_opt_result_geometry, trc_opt_result_message
   use trc_test_basis, only: read_xyz
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_report
   implicit none
   private
   public :: collect_capi_opt

   real(dp), parameter :: ENERGY_TOL = 5.0e-8_dp
   real(dp), parameter :: BOND_TOL = 3.0e-4_dp
   real(dp), parameter :: ANGLE_TOL = 0.05_dp
      !! Identical to check_opt_vs_pyscf.F90's own (see that file's own
      !! header for the harmonic-estimate derivation from --opt-gtol).
   real(dp), parameter :: RUN_GTOL = 1.0e-5_dp
      !! Same --opt-gtol check_opt_vs_pyscf.F90 passes the CLI.

contains

   subroutine collect_capi_opt(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("terco_opt.h vs pyscf+scipy (opt_water.ref), through the C ABI only", test_water), &
         new_unittest("terco_opt.h vs pyscf+scipy (opt_hcn.ref), through the C ABI only", test_hcn) &
         ]
   end subroutine collect_capi_opt

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case("water.xyz", "opt_water", error)
   end subroutine test_water

   subroutine test_hcn(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case("hcn.xyz", "opt_hcn", error)
   end subroutine test_hcn

   subroutine run_case(xyzfile, ref_base, error)
      character(len=*), intent(in) :: xyzfile, ref_base
      type(error_type), allocatable, intent(out) :: error

      integer :: nat, rc, nfail
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      integer(c_int), allocatable :: z_c(:)
      type(c_ptr) :: ctx, opts, res
      integer(c_int) :: converged
      real(c_double) :: energy, gmax
      real(c_double), allocatable :: xyz_out(:)
      character(kind=c_char) :: msg(256)
      real(dp), allocatable :: coords_out(:, :)

      integer :: nat_ref
      real(dp), allocatable :: z_ref(:), coords_ref(:, :)
      real(dp) :: e_ref, gmax_ref
      real(dp) :: d1_ref, d2_ref, d3_ref, ang_ref
      real(dp) :: d1_out, d2_out, d3_out, ang_out

      nfail = 0

      ! --- context setup, exactly what a C caller (cyclone) would do: -----
      ! molecule, basis and method through terco.h, never trc_dlfind_bridge
      ! or trc_api directly.
      call read_xyz(xyzfile, nat, z, r)
      allocate (z_c(nat), source=int(z, c_int))

      rc = trc_create(ctx)
      if (rc /= TRC_OK) then
         call test_failed(error, "check_capi_opt ("//ref_base//"): trc_create failed, status "//to_string(rc))
         return
      end if
      rc = trc_set_molecule(ctx, int(nat, c_int), real(z_c, c_double), &
                            reshape(real(r, c_double), [3*nat]), 0_c_int, 1_c_int)
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_set_molecule", real(rc, dp))
      rc = trc_set_basis_json(ctx, "basis_sets/sto-3g.json"//c_null_char)
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_set_basis_json", real(rc, dp))
      rc = trc_set_method(ctx, ""//c_null_char, 3_c_int)   ! empty functional: Hartree-Fock
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_set_method", real(rc, dp))

      ! --- the optimiser, through terco_opt.h ONLY -------------------------
      rc = trc_opt_options_create(opts)
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_opt_options_create", real(rc, dp))
      rc = trc_opt_set_gtol(opts, real(RUN_GTOL, c_double))
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_opt_set_gtol", real(rc, dp))
      rc = trc_opt_result_create(res)
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_opt_result_create", real(rc, dp))

      rc = trc_opt_run(ctx, int(nat, c_int), z_c, reshape(real(r, c_double), [3*nat]), opts, res)
      if (rc /= TRC_OK) then
         block
            integer :: st2
            st2 = trc_opt_result_message(res, msg, int(size(msg), c_int))
            call test_failed(error, "check_capi_opt ("//ref_base//"): trc_opt_run returned "//to_string(rc)// &
                             ": "//f_message(msg))
         end block
         return
      end if

      rc = trc_opt_result_converged(res, converged)
      call ref_report(nfail, converged /= 0, "check_capi_opt ("//ref_base//"): converged", real(converged, dp))
      rc = trc_opt_result_energy(res, energy)
      rc = trc_opt_result_gradient_max(res, gmax)

      allocate (xyz_out(3*nat))
      rc = trc_opt_result_geometry(res, int(nat, c_int), xyz_out)
      call ref_report(nfail, rc == TRC_OK, "check_capi_opt ("//ref_base//"): trc_opt_result_geometry", real(rc, dp))
      allocate (coords_out(3, nat), source=reshape(real(xyz_out, dp), [3, nat]))

      ! --- compare against the SAME pyscf+scipy anchor check_opt_vs_pyscf
      ! (both water and hcn cases) uses, at the SAME tolerances. ------------
      block
         integer :: uu
         call ref_require(ref_base//".ref", uu)
         call ref_read_int(uu, "NAT", nat_ref)
         allocate (z_ref(nat_ref), coords_ref(3, nat_ref))
         call ref_read_vector(uu, "Z", nat_ref, z_ref)
         call ref_read_matrix(uu, "COORDS", 3, nat_ref, coords_ref)
         call ref_read_real(uu, "E_FINAL", e_ref)
         call ref_read_real(uu, "GMAX_FINAL", gmax_ref)
         call ref_close(uu)
      end block
      if (nat_ref /= nat) then
         call test_failed(error, "check_capi_opt ("//ref_base//"): NAT disagrees with the reference")
         return
      end if

      call ref_report(nfail, abs(real(energy, dp) - e_ref) < ENERGY_TOL, &
                      "check_capi_opt ("//ref_base//"): final energy vs pyscf+scipy, Eh", &
                      abs(real(energy, dp) - e_ref))

      call two_shortest_and_angle(coords_ref, d1_ref, d2_ref, d3_ref, ang_ref)
      call two_shortest_and_angle(coords_out, d1_out, d2_out, d3_out, ang_out)
      call ref_report(nfail, abs(d1_out - d1_ref) < BOND_TOL, &
                      "check_capi_opt ("//ref_base//"): shorter bond vs pyscf+scipy, Bohr", abs(d1_out - d1_ref))
      call ref_report(nfail, abs(d2_out - d2_ref) < BOND_TOL, &
                      "check_capi_opt ("//ref_base//"): longer bond vs pyscf+scipy, Bohr", abs(d2_out - d2_ref))
      call ref_report(nfail, abs(d3_out - d3_ref) < BOND_TOL, &
                      "check_capi_opt ("//ref_base//"): third distance vs pyscf+scipy, Bohr", &
                      abs(d3_out - d3_ref))
      if (ang_ref > 0.0_dp) then
         call ref_report(nfail, abs(ang_out - ang_ref) < ANGLE_TOL, &
                         "check_capi_opt ("//ref_base//"): bond angle vs pyscf+scipy, deg", &
                         abs(ang_out - ang_ref))
      end if

      rc = trc_opt_result_destroy(res)
      rc = trc_opt_options_destroy(opts)
      rc = trc_destroy(ctx)

      if (nfail > 0) call test_failed(error, to_string(nfail)//" check(s) failed (see the FAIL lines above)")
   end subroutine run_case

   pure function f_message(buf) result(s)
      !! The text of a NUL-terminated C buffer -- trc_opt_result_message's
      !! own output, read back here the same way a real C caller would
      !! (strlen then copy), not through any terco-internal helper.
      character(kind=c_char), intent(in) :: buf(:)
      character(len=:), allocatable :: s
      integer :: i, n
      n = size(buf)
      do i = 1, size(buf)
         if (buf(i) == c_null_char) then
            n = i - 1
            exit
         end if
      end do
      allocate (character(len=n) :: s)
      do i = 1, n
         s(i:i) = buf(i)
      end do
   end function f_message

   !
   ! For a 3-atom molecule only (water, hcn) -- the two SHORTEST pairwise
   ! distances (the two bonds), the remaining distance, and the angle at
   ! their shared vertex. A direct copy of check_opt_vs_pyscf.F90's own
   ! helper of the same name: small enough, and stable enough, that
   ! duplicating it here was preferred over exporting it from that file's
   ! own private module scope.
   !
   subroutine two_shortest_and_angle(xyz, d1, d2, d3, angle_deg)
      real(dp), intent(in) :: xyz(3, 3)
      real(dp), intent(out) :: d1, d2, d3, angle_deg
      real(dp), parameter :: PI_LOCAL = 3.14159265358979323846_dp
      real(dp) :: d(3), v1(3), v2(3)
      integer :: apex
      d(1) = norm2(xyz(:, 1) - xyz(:, 2))   ! 1-2
      d(2) = norm2(xyz(:, 1) - xyz(:, 3))   ! 1-3
      d(3) = norm2(xyz(:, 2) - xyz(:, 3))   ! 2-3
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
         angle_deg = acos(max(-1.0_dp, min(1.0_dp, dot_product(v1, v2)/(norm2(v1)*norm2(v2)))))*180.0_dp/PI_LOCAL
      else
         angle_deg = 0.0_dp
      end if
   end subroutine two_shortest_and_angle

end module check_capi_opt_suite

program check_capi_opt
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use, intrinsic :: iso_fortran_env, only: error_unit
   use check_capi_opt_suite, only: collect_capi_opt
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat
   stat = 0
   testsuites = [new_testsuite("check_capi_opt", collect_capi_opt)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat)
   if (stat > 0) then
      write (error_unit, '(i0,1x,a)') stat, "test(s) failed"
      stop 1
   end if
end program check_capi_opt
