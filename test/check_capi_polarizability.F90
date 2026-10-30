module check_capi_polarizability_suite
   !! trc_run_polarizability/trc_polarizability (CPHF for RHF, CPKS for
   !! RKS) and trc_nmo/trc_mo_occupations, through the context API, on
   !! water/6-31G*: the same case test/scf_cphf.F90 and test/scf_cpks.F90
   !! check their own Fortran drivers against (test/pyscf_refs/
   !! cphf_water_6-31gs.ref, cpks_water_6-31gs_svwn.ref). Forced Cartesian,
   !! same reason test/context_scf.F90 is: trc_run_polarizability refuses
   !! a spherical basis outright (trc_capi.F90's own guard, CPHF/CPKS know
   !! nothing of nao_sph).
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, &
                               trc_run_scf, trc_energy, trc_nao, trc_nmo, trc_mo_occupations, &
                               trc_run_polarizability, trc_polarizability, trc_nspin, &
                               TRC_OK
   implicit none
   private
   public :: collect_capi_polarizability_tests

contains

   subroutine collect_capi_polarizability_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("CPHF (RHF) energy and polarizability match pyscf", test_cphf), &
                  new_unittest("CPKS (SVWN) energy and polarizability match pyscf", test_cpks), &
                  new_unittest("trc_nmo/trc_mo_occupations agree with water's electron count", &
                               test_nmo_and_occupations) &
                  ]
   end subroutine collect_capi_polarizability_tests

   subroutine water_geometry(natm, zint, at_r)
      integer, intent(out) :: natm
      integer, allocatable, intent(out) :: zint(:)
      real(dp), allocatable, intent(out) :: at_r(:, :)
      call read_xyz('water.xyz', natm, zint, at_r)
   end subroutine water_geometry

   !
   ! test/pyscf_refs/cphf_water_6-31gs.ref: the SAME reference
   ! test/scf_cphf.F90 checks its own Fortran CPHF driver against, here
   ! read back through the context API's trc_run_polarizability/
   ! trc_polarizability (RHF, CPHF).
   !
   subroutine test_cphf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_scf, alpha(3, 3)
      integer :: u
      real(dp) :: e_scf_ref, alpha_ref(3, 3), worst_alpha

      call water_geometry(natm, zint, at_r)

      h = c_null_ptr
      ! Every C call's status is checked, trc_destroy's included.
      steps: block
         rc = trc_create(h)
         call check(error, rc == TRC_OK, "cphf: create")
         if (allocated(error)) exit steps
         rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
         call check(error, rc == TRC_OK, "cphf: set_molecule")
         if (allocated(error)) exit steps
         rc = trc_set_force_cartesian(h, 1_c_int)
         call check(error, rc == TRC_OK, "cphf: set_force_cartesian")
         if (allocated(error)) exit steps
         rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
         call check(error, rc == TRC_OK, "cphf: set_basis_json")
         if (allocated(error)) exit steps
         rc = trc_set_method(h, ''//c_null_char, 3_c_int)
         call check(error, rc == TRC_OK, "cphf: set_method (HF)")
         if (allocated(error)) exit steps
         rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 200_c_int, 10_c_int)
         call check(error, rc == TRC_OK, "cphf: set_convergence")
         if (allocated(error)) exit steps
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "cphf: run_scf must converge")
         if (allocated(error)) exit steps
         rc = trc_energy(h, e_scf)
         call check(error, rc == TRC_OK, "cphf: energy")
         if (allocated(error)) exit steps

         rc = trc_run_polarizability(h)
         call check(error, rc == TRC_OK, "cphf: run_polarizability status")
         if (allocated(error)) exit steps
         rc = trc_polarizability(h, alpha)
         call check(error, rc == TRC_OK, "cphf: polarizability status")
         if (allocated(error)) exit steps

         call ref_require('cphf_water_6-31gs.ref', u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
         call ref_close(u)

         call check(error, abs(real(e_scf, dp) - e_scf_ref) < 1.0e-9_dp, &
                    "cphf: E_RHF matches pyscf, terco "//to_string(real(e_scf, dp), "f20.12")// &
                    "  pyscf "//to_string(e_scf_ref, "f20.12")//"  diff "// &
                    to_string(abs(real(e_scf, dp) - e_scf_ref), "es9.2"))
         if (allocated(error)) exit steps
         ! all() per element, as the old per-element checks did: maxval
         ! skips a NaN unless every element is one.
         worst_alpha = maxval(abs(real(alpha, dp) - alpha_ref))
         call check(error, all(abs(real(alpha, dp) - alpha_ref) < 1.0e-6_dp), &
                    "cphf: alpha matches pyscf, worst |dalpha| = "//to_string(worst_alpha, "es10.2"))
         if (allocated(error)) exit steps

         rc = trc_destroy(h)
         h = c_null_ptr
         call check(error, rc == TRC_OK, "cphf: destroy")
      end block steps

      if (allocated(error)) rc = trc_destroy(h)
   end subroutine test_cphf

   !
   ! test/pyscf_refs/cpks_water_6-31gs_svwn.ref: the SAME reference
   ! test/scf_cpks.F90 checks its own Fortran CPKS driver against, here
   ! read back through the context API (RKS/SVWN, CPKS).
   !
   subroutine test_cpks(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_scf, alpha(3, 3)
      integer :: u
      real(dp) :: e_scf_ref, alpha_ref(3, 3), worst_alpha

      call water_geometry(natm, zint, at_r)

      h = c_null_ptr
      ! Every C call's status is checked, trc_destroy's included.
      steps: block
         rc = trc_create(h)
         call check(error, rc == TRC_OK, "cpks: create")
         if (allocated(error)) exit steps
         rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
         call check(error, rc == TRC_OK, "cpks: set_molecule")
         if (allocated(error)) exit steps
         rc = trc_set_force_cartesian(h, 1_c_int)
         call check(error, rc == TRC_OK, "cpks: set_force_cartesian")
         if (allocated(error)) exit steps
         rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
         call check(error, rc == TRC_OK, "cpks: set_basis_json")
         if (allocated(error)) exit steps
         rc = trc_set_method(h, 'svwn'//c_null_char, 3_c_int)
         call check(error, rc == TRC_OK, "cpks: set_method (svwn)")
         if (allocated(error)) exit steps
         rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 200_c_int, 10_c_int)
         call check(error, rc == TRC_OK, "cpks: set_convergence")
         if (allocated(error)) exit steps
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "cpks: run_scf must converge")
         if (allocated(error)) exit steps
         rc = trc_energy(h, e_scf)
         call check(error, rc == TRC_OK, "cpks: energy")
         if (allocated(error)) exit steps

         rc = trc_run_polarizability(h)
         call check(error, rc == TRC_OK, "cpks: run_polarizability status")
         if (allocated(error)) exit steps
         rc = trc_polarizability(h, alpha)
         call check(error, rc == TRC_OK, "cpks: polarizability status")
         if (allocated(error)) exit steps

         call ref_require('cpks_water_6-31gs_svwn.ref', u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
         call ref_close(u)

         call check(error, abs(real(e_scf, dp) - e_scf_ref) < 1.0e-8_dp, &
                    "cpks: E_SCF matches pyscf, terco "//to_string(real(e_scf, dp), "f20.12")// &
                    "  pyscf "//to_string(e_scf_ref, "f20.12")//"  diff "// &
                    to_string(abs(real(e_scf, dp) - e_scf_ref), "es9.2"))
         if (allocated(error)) exit steps
         ! all() per element, as the old per-element checks did: maxval
         ! skips a NaN unless every element is one.
         worst_alpha = maxval(abs(real(alpha, dp) - alpha_ref))
         call check(error, all(abs(real(alpha, dp) - alpha_ref) < 1.0e-5_dp), &
                    "cpks: alpha matches pyscf, worst |dalpha| = "//to_string(worst_alpha, "es10.2"))
         if (allocated(error)) exit steps

         rc = trc_destroy(h)
         h = c_null_ptr
         call check(error, rc == TRC_OK, "cpks: destroy")
      end block steps

      if (allocated(error)) rc = trc_destroy(h)
   end subroutine test_cpks

   !
   ! trc_nmo/trc_mo_occupations: water/6-31G* has no linear dependency, so
   ! nmo = nao, and the occupation vector must be exactly 2.0 for the
   ! lowest `nocc` orbitals, 0 above.
   !
   subroutine test_nmo_and_occupations(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm, nocc
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(c_ptr) :: h
      integer(c_int) :: rc, nao_c, nmo_c, nspin_c
      real(c_double), allocatable :: occ(:)
      integer :: nocc_in
      real(dp) :: worst_occ, worst_virt

      call water_geometry(natm, zint, at_r)
      nocc = sum(zint)/2

      h = c_null_ptr
      steps: block
         rc = trc_create(h)
         call check(error, rc == TRC_OK, "nmo: create")
         if (allocated(error)) exit steps
         rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
         rc = trc_set_force_cartesian(h, 1_c_int)
         rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
         rc = trc_nmo(h, nmo_c)
         call check(error, rc == TRC_OK, "nmo: before any SCF, status OK")
         if (allocated(error)) exit steps
         rc = trc_nao(h, nao_c)
         call check(error, nmo_c == nao_c, "nmo: equals nao before any SCF (n_lindep unknown => 0)")
         if (allocated(error)) exit steps

         rc = trc_set_method(h, ''//c_null_char, 3_c_int)
         rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "nmo: run_scf must converge")
         if (allocated(error)) exit steps
         rc = trc_nao(h, nao_c)
         call check(error, rc == TRC_OK, "nmo: nao status")
         if (allocated(error)) exit steps
         rc = trc_nmo(h, nmo_c)
         call check(error, rc == TRC_OK, "nmo: nmo status")
         if (allocated(error)) exit steps
         call check(error, nmo_c == nao_c, "nmo: water/6-31G* has no linear dependency, nmo = nao")
         if (allocated(error)) exit steps
         rc = trc_nspin(h, nspin_c)
         call check(error, nspin_c == 1, "nmo: restricted, nspin 1")
         if (allocated(error)) exit steps

         allocate (occ(nmo_c))
         rc = trc_mo_occupations(h, occ)
         call check(error, rc == TRC_OK, "nmo: mo_occupations status")
         if (allocated(error)) exit steps

         ! One check per logical assertion (occupied vs virtual), each an
         ! all() over the orbitals so a NaN fails it; the worst case is
         ! only reported.
         nocc_in = min(nocc, int(nmo_c))
         worst_occ = maxval(abs(real(occ(1:nocc_in), dp) - 2.0_dp))
         call check(error, all(abs(real(occ(1:nocc_in), dp) - 2.0_dp) < 1.0e-12_dp), &
                    "nmo: every occupied orbital is 2.0, worst |occ-2| = "//to_string(worst_occ, "es10.2"))
         if (allocated(error)) exit steps
         if (nocc_in < int(nmo_c)) then
            worst_virt = maxval(abs(real(occ(nocc_in + 1:), dp)))
         else
            worst_virt = 0.0_dp
         end if
         call check(error, all(abs(real(occ(nocc_in + 1:), dp)) < 1.0e-12_dp), &
                    "nmo: every virtual orbital is 0.0, worst |occ| = "//to_string(worst_virt, "es10.2"))
         if (allocated(error)) exit steps

         rc = trc_destroy(h)
         h = c_null_ptr
         call check(error, rc == TRC_OK, "nmo: destroy")
      end block steps

      if (allocated(error)) rc = trc_destroy(h)
   end subroutine test_nmo_and_occupations

end module check_capi_polarizability_suite

program check_capi_polarizability
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_capi_polarizability_suite, only: collect_capi_polarizability_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   ! Serial: each case builds its own context through the C ABI, which is
   ! not thread-safe, and an -fopenmp build runs test-drive cases at once.
   testsuites = [new_testsuite("check_capi_polarizability", collect_capi_polarizability_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_capi_polarizability
