program check_capi_polarizability
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
   use trc_api, only: trc_bind_device
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix, ref_report
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, &
                               trc_run_scf, trc_energy, trc_nao, trc_nmo, trc_mo_occupations, &
                               trc_run_polarizability, trc_polarizability, trc_nspin, &
                               TRC_OK, TRC_ERR_STATE
   implicit none

   integer :: natm, dev, nocc
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nocc = sum(zint)/2
   ok = .true.

   call check_cphf()
   call check_cpks()
   call check_nmo_and_occupations()

   if (.not. ok) then
      print '(a)', "check_capi_polarizability: FAIL"
      stop 1
   end if
   print '(a)', "check_capi_polarizability: PASS"

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "check_capi_polarizability: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

   subroutine check_cphf()
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_scf, alpha(3, 3)
      integer :: u, x, y, nfail
      real(dp) :: e_scf_ref, alpha_ref(3, 3)

      nfail = 0
      rc = trc_create(h);                     call expect(rc == TRC_OK, "cphf: create")
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc == TRC_OK, "cphf: set_molecule")
      rc = trc_set_force_cartesian(h, 1_c_int); call expect(rc == TRC_OK, "cphf: set_force_cartesian")
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      call expect(rc == TRC_OK, "cphf: set_basis_json")
      rc = trc_set_method(h, ''//c_null_char, 3_c_int); call expect(rc == TRC_OK, "cphf: set_method (HF)")
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 200_c_int, 10_c_int)
      call expect(rc == TRC_OK, "cphf: set_convergence")
      rc = trc_run_scf(h);                    call expect(rc == TRC_OK, "cphf: run_scf")
      rc = trc_energy(h, e_scf);               call expect(rc == TRC_OK, "cphf: energy")

      rc = trc_run_polarizability(h);          call expect(rc == TRC_OK, "cphf: run_polarizability")
      rc = trc_polarizability(h, alpha);       call expect(rc == TRC_OK, "cphf: polarizability")

      call ref_require('cphf_water_6-31gs.ref', u)
      call ref_read_real(u, 'E_SCF', e_scf_ref)
      call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
      call ref_close(u)
      call ref_report(nfail, abs(real(e_scf, dp) - e_scf_ref) < 1.0e-9_dp, 'cphf: E_RHF vs pyscf', &
                      abs(real(e_scf, dp) - e_scf_ref))
      call expect(abs(real(e_scf, dp) - e_scf_ref) < 1.0e-9_dp, "cphf: E_RHF matches pyscf")
      do y = 1, 3
         do x = 1, 3
            call expect(abs(real(alpha(x, y), dp) - alpha_ref(x, y)) < 1.0e-6_dp, "cphf: alpha matches pyscf")
         end do
      end do

      rc = trc_destroy(h);                    call expect(rc == TRC_OK, "cphf: destroy")
   end subroutine check_cphf

   subroutine check_cpks()
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_scf, alpha(3, 3)
      integer :: u, x, y
      real(dp) :: e_scf_ref, alpha_ref(3, 3)

      rc = trc_create(h);                     call expect(rc == TRC_OK, "cpks: create")
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call expect(rc == TRC_OK, "cpks: set_molecule")
      rc = trc_set_force_cartesian(h, 1_c_int); call expect(rc == TRC_OK, "cpks: set_force_cartesian")
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      call expect(rc == TRC_OK, "cpks: set_basis_json")
      rc = trc_set_method(h, 'svwn'//c_null_char, 3_c_int); call expect(rc == TRC_OK, "cpks: set_method (svwn)")
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 200_c_int, 10_c_int)
      call expect(rc == TRC_OK, "cpks: set_convergence")
      rc = trc_run_scf(h);                    call expect(rc == TRC_OK, "cpks: run_scf")
      rc = trc_energy(h, e_scf);               call expect(rc == TRC_OK, "cpks: energy")

      rc = trc_run_polarizability(h);          call expect(rc == TRC_OK, "cpks: run_polarizability")
      rc = trc_polarizability(h, alpha);       call expect(rc == TRC_OK, "cpks: polarizability")

      call ref_require('cpks_water_6-31gs_svwn.ref', u)
      call ref_read_real(u, 'E_SCF', e_scf_ref)
      call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
      call ref_close(u)
      call expect(abs(real(e_scf, dp) - e_scf_ref) < 1.0e-8_dp, "cpks: E_SCF matches pyscf")
      do y = 1, 3
         do x = 1, 3
            call expect(abs(real(alpha(x, y), dp) - alpha_ref(x, y)) < 1.0e-5_dp, "cpks: alpha matches pyscf")
         end do
      end do

      rc = trc_destroy(h);                    call expect(rc == TRC_OK, "cpks: destroy")
   end subroutine check_cpks

   !
   ! trc_nmo/trc_mo_occupations: water/6-31G* has no linear dependency, so
   ! nmo = nao, and the occupation vector must be exactly 2.0 for the
   ! lowest `nocc` orbitals, 0 above.
   !
   subroutine check_nmo_and_occupations()
      type(c_ptr) :: h
      integer(c_int) :: rc, nao_c, nmo_c, nspin_c
      real(c_double), allocatable :: occ(:)
      integer :: i

      rc = trc_create(h);                     call expect(rc == TRC_OK, "nmo: create")
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      rc = trc_nmo(h, nmo_c)
      call expect(rc == TRC_OK, "nmo: before any SCF, status OK")
      rc = trc_nao(h, nao_c)
      call expect(nmo_c == nao_c, "nmo: equals nao before any SCF (n_lindep unknown => 0)")

      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h);                    call expect(rc == TRC_OK, "nmo: run_scf")
      rc = trc_nao(h, nao_c);                 call expect(rc == TRC_OK, "nmo: nao")
      rc = trc_nmo(h, nmo_c);                 call expect(rc == TRC_OK, "nmo: nmo")
      call expect(nmo_c == nao_c, "nmo: water/6-31G* has no linear dependency, nmo = nao")
      rc = trc_nspin(h, nspin_c);              call expect(nspin_c == 1, "nmo: restricted, nspin 1")

      allocate (occ(nmo_c))
      rc = trc_mo_occupations(h, occ);        call expect(rc == TRC_OK, "nmo: mo_occupations")
      do i = 1, int(nmo_c)
         if (i <= nocc) then
            call expect(abs(real(occ(i), dp) - 2.0_dp) < 1.0e-12_dp, "nmo: occupied orbital is 2.0")
         else
            call expect(abs(real(occ(i), dp)) < 1.0e-12_dp, "nmo: virtual orbital is 0.0")
         end if
      end do
      deallocate (occ)

      rc = trc_destroy(h);                    call expect(rc == TRC_OK, "nmo: destroy")
   end subroutine check_nmo_and_occupations

end program check_capi_polarizability
