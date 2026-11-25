module ccsd_grad_support
   !!
   !! What the RI-CCSD(T) gradient-reference checks share (check_ccsd_grad_ref, check_ccsd_grad_fd): the reader of
   !! test/pyscf_refs/ccsd_grad_water_*.ref (generate.py's gen_ccsd_grad) and a driver that runs terco's own RI-CCSD(T)
   !! through the C ABI at an arbitrary geometry, with the tight convergence of the references.
   !!
   !! The molecule is test/ccgrad_water.xyz (water distorted to C1), read in Angstrom and converted to Bohr with CODATA
   !! 2018, as terco and the pyscf reference do. The orbital basis is STO-3G or cc-pVDZ, the auxiliary basis
   !! cc-pVDZ-RIFIT in the orbital basis's convention
   !! (Cartesian with STO-3G, spherical with cc-pVDZ).
   !!
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_matrix
   use rccsd_test_support, only: read_xyz_bohr
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_convergence, trc_set_incremental, trc_set_rccsd, trc_set_rccsd_convergence, &
                               trc_run_scf, trc_run_rccsd, trc_rccsd_energy, trc_energy, TRC_OK
   implicit none
   private
   public :: NCASE, CASE_REF, CASE_BASIS, CASE_NFROZEN, NPROBE, PROBE_NAME, XYZ
   public :: ccsd_grad_ref_t, read_ccsd_grad_ref, terco_energies, geometry, supported

   character(len=*), parameter :: XYZ = "ccgrad_water.xyz"
   integer, parameter :: NCASE = 4
   character(len=32), parameter :: CASE_REF(NCASE) = [character(len=32) :: "ccsd_grad_water_sto3g_fc0.ref", &
      "ccsd_grad_water_sto3g_fc1.ref", "ccsd_grad_water_ccpvdz_fc0.ref", "ccsd_grad_water_ccpvdz_fc1.ref"]
   character(len=8), parameter :: CASE_BASIS(NCASE) = [character(len=8) :: "sto-3g", "sto-3g", "cc-pvdz", "cc-pvdz"]
   integer, parameter :: CASE_NFROZEN(NCASE) = [0, 1, 0, 1]

   integer, parameter :: NPROBE = 10
   character(len=16), parameter :: PROBE_NAME(NPROBE) = [character(len=16) :: "L1_NORM", "L2_NORM", "L1_T1", "L2_T2", &
      "GAMMA_SQ", "GAMMA_TR_OCC", "GAMMA_TR_VIR", "G3_NORM", "GAMMA2_NORM", "E_FROM_DENSITIES"]

   type :: ccsd_grad_ref_t
      !! one reference file; probes(:, 1) are those of CCSD and probes(:, 2) of CCSD(T), in PROBE_NAME order
      integer :: natm = 0, no = 0, nv = 0, naux = 0, nfrozen = 0
      real(dp) :: e_scf = 0.0_dp, e_ccsd_corr = 0.0_dp, e_t = 0.0_dp, e_ccsd_t = 0.0_dp
      real(dp) :: probes(NPROBE, 2) = 0.0_dp
      real(dp), allocatable :: grad(:, :, :)   !! (natm, 3, 2): CCSD, then CCSD(T), Eh/Bohr
   end type ccsd_grad_ref_t

contains

   function supported() result(ok)
      !! cc-pVDZ is d and cc-pVDZ-RIFIT has f on O: both within the build's density-fitting ceilings
      logical :: ok
      ok = DF_LMAX >= 2 .and. DF_LMAX_AUX >= 3
   end function supported

   subroutine read_ccsd_grad_ref(name, ref)
      !! Reads `pyscf_refs/<name>`; a missing file stops the program (run generate.py ccsd_grad)
      character(len=*), intent(in) :: name
      type(ccsd_grad_ref_t), intent(out) :: ref

      integer :: u, k, v
      real(dp) :: value
      character(len=2) :: suffix

      call ref_require(name, u)
      call ref_read_int(u, "NATM", ref%natm)
      call ref_read_int(u, "NO", ref%no)
      call ref_read_int(u, "NV", ref%nv)
      call ref_read_int(u, "NAUX", ref%naux)
      call ref_read_int(u, "NFROZEN", ref%nfrozen)
      call ref_read_real(u, "E_SCF", ref%e_scf)
      call ref_read_real(u, "E_CCSD_CORR", ref%e_ccsd_corr)
      call ref_read_real(u, "E_T", ref%e_t)
      call ref_read_real(u, "E_CCSD_T", ref%e_ccsd_t)
      allocate (ref%grad(ref%natm, 3, 2))
      do v = 1, 2
         suffix = merge("  ", "_T", v == 1)
         do k = 1, NPROBE
            call ref_read_real(u, trim(PROBE_NAME(k))//trim(suffix), value)
            ref%probes(k, v) = value
         end do
         call ref_read_matrix(u, "GRAD"//trim(suffix), ref%natm, 3, ref%grad(:, :, v))
      end do
      call ref_close(u)
   end subroutine read_ccsd_grad_ref

   subroutine geometry(natm, zint, at_r)
      !! test/ccgrad_water.xyz in Bohr
      integer, intent(out) :: natm
      integer, allocatable, intent(out) :: zint(:)
      real(dp), allocatable, intent(out) :: at_r(:, :)
      call read_xyz_bohr(XYZ, natm, zint, at_r)
   end subroutine geometry

   subroutine terco_energies(basis, nfrozen, zint, at_r, e_scf, e_corr, e_t, rc, iterations)
      !! terco's RI-CCSD(T) at the geometry `at_r` (3, natm) in Bohr: exact-ERI RHF to energy 1e-13 and commutator
      !! 1e-11, full Fock builds, then the coupled-cluster iteration to 1e-12 / 1e-10. e_scf is the SCF energy
      !! (trc_energy less the correlation), e_corr the CCSD correlation energy and e_t the (T) correction.
      character(len=*), intent(in) :: basis
      integer, intent(in) :: nfrozen, zint(:)
      real(dp), intent(in) :: at_r(:, :)
      real(dp), intent(out) :: e_scf, e_corr, e_t
      integer, intent(out) :: rc
      integer, intent(out), optional :: iterations

      type(c_ptr) :: h
      real(c_double) :: e_mp2, e_ccsd, e_triples, t1, e_total
      integer(c_int) :: stat, iters, converged, has_triples
      integer :: natm

      natm = size(zint)
      e_scf = 0.0_dp
      e_corr = 0.0_dp
      e_t = 0.0_dp
      if (present(iterations)) iterations = 0
      rc = int(trc_create(h))
      if (rc /= TRC_OK) return
      stat = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      if (stat == TRC_OK) stat = trc_set_basis_json(h, 'basis_sets/'//trim(basis)//'.json'//c_null_char)
      if (stat == TRC_OK) stat = trc_set_aux_json(h, 'basis_sets/cc-pvdz-rifit.json'//c_null_char)
      if (stat == TRC_OK) stat = trc_set_convergence(h, 1.0e-13_c_double, 1.0e-11_c_double, 100_c_int, 10_c_int)
      ! Full Fock builds: an incremental build's screening error sits near these thresholds.
      if (stat == TRC_OK) stat = trc_set_incremental(h, 0_c_int)
      if (stat == TRC_OK) stat = trc_set_rccsd(h, 1_c_int, int(nfrozen, c_int), 0_c_int)
      if (stat == TRC_OK) stat = trc_set_rccsd_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 100_c_int, 0_c_int)
      if (stat == TRC_OK) stat = trc_run_scf(h)
      if (stat == TRC_OK) stat = trc_run_rccsd(h)
      if (stat == TRC_OK) stat = trc_rccsd_energy(h, e_mp2, e_ccsd, e_triples, t1, iters, converged, has_triples)
      if (stat == TRC_OK .and. (converged /= 1 .or. has_triples /= 1)) stat = -1_c_int
      if (stat == TRC_OK) stat = trc_energy(h, e_total)
      if (stat == TRC_OK) then
         e_corr = e_ccsd
         e_t = e_triples
         e_scf = e_total - e_ccsd - e_triples
         if (present(iterations)) iterations = int(iters)
      end if
      rc = int(stat)
      stat = trc_destroy(h)
   end subroutine terco_energies

end module ccsd_grad_support
