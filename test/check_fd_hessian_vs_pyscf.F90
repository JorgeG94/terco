!
! trc_run_hessian/trc_hessian (the C API's finite-difference Cartesian
! Hessian, src/api/trc_capi.F90) against pyscf's own ANALYTIC RHF/STO-3G
! Hessian, the SAME test/pyscf_refs/vib_thermo_{water,hcn,nh3_planar}.ref
! check_vib_thermo_vs_pyscf.F90 already reads (HESS there is pyscf's
! `mf.Hessian().kernel()`, not a finite difference of anything) -- no new
! generate.py case needed, see this program's own header below for why.
!
! BASIS. test/basis_sets/sto-3g.json, which generate.py's gen_vib_thermo
! also loads for pyscf (it is on the .ref's files: line), so there is no
! basis mismatch left to measure beyond the usual SCF/gradient roundoff:
! the electronic energy at the reference geometry agrees with E_ELEC in
! the .ref to << 1e-9 Eh.
!
! TOLERANCES. The hessian_tol of each case is roughly 1.5-2x what this
! build measures (reported on stdout every run). The only claim this
! test makes is "close to the h=0.005 truncation-error floor", a property
! of the finite difference, not of terco. The frequency tolerance is set
! the same way: mode frequencies are much less sensitive to Hessian noise
! than the raw elements (projecting out translation/rotation and
! diagonalising averages over the whole matrix), about 0.06 cm^-1 here.
!
! A test-drive suite: each case is one unit test, which fails when any of
! its checks does; every check still prints its value.
!
module check_fd_hessian_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
   use pic_types, only: dp
   use trc_error, only: error_t
   use trc_constants, only: trc_atomic_mass
   use trc_vib, only: trc_vibrational_analysis
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, &
                               trc_set_screening, trc_set_guess, trc_set_verbose, trc_run_scf, &
                               trc_run_hessian, trc_hessian, trc_set_hessian_step, trc_energy, &
                               trc_iterations, TRC_OK, TRC_ERR_UNSUPPORTED, TRC_ERR_STATE, TRC_GUESS_SAD
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, &
                           ref_read_vector, ref_read_matrix, ref_report
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_fd_hessian

contains

   subroutine collect_fd_hessian(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water FD Hessian vs pyscf", test_water), &
                  new_unittest("hcn FD Hessian vs pyscf", test_hcn), &
                  new_unittest("planar nh3 FD Hessian vs pyscf", test_nh3_planar), &
                  new_unittest("a meta-GGA is refused before any SCF", test_refusal) &
                  ]
   end subroutine collect_fd_hessian

   subroutine verdict(error, nfail, what)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nfail
      character(len=*), intent(in) :: what
      character(len=12) :: buf
      if (nfail == 0) return
      write (buf, "(i0)") nfail
      call test_failed(error, what//": "//trim(buf)//" check(s) failed (see the FAIL lines above)")
   end subroutine verdict

   ! Tolerances: ~1.5-2x what this build measures (printed on every run):
   ! water 1.441e-5, hcn 7.402e-5, nh3_planar 2.331e-5 Eh/Bohr^2 for the
   ! Hessian; 0.061, 0.058, 0.059 cm^-1 for the frequencies (gfortran 15
   ! host build). All are the h=0.005 truncation error of a central
   ! difference of an analytic gradient.
   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_case(nfail, "water", 2.5e-5_dp, 0.12_dp)
      call verdict(error, nfail, "water")
   end subroutine test_water

   subroutine test_hcn(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_case(nfail, "hcn", 1.3e-4_dp, 0.12_dp)
      call verdict(error, nfail, "hcn")
   end subroutine test_hcn

   subroutine test_nh3_planar(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_case(nfail, "nh3_planar", 4.0e-5_dp, 0.12_dp)
      call verdict(error, nfail, "nh3_planar")
   end subroutine test_nh3_planar

   subroutine test_refusal(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_refusal(nfail)
      call verdict(error, nfail, "refusal")
   end subroutine test_refusal

   !
   ! One case: build the context at pyscf's own geometry/charge/mult
   ! (read straight out of the .ref, Bohr, no xyz round trip), run
   ! trc_run_scf/trc_run_hessian/trc_hessian through the C ABI
   ! (trc_c_interfaces, the generated Fortran declaration an external
   ! caller compiles, not trc_capi directly), and compare against pyscf's
   ! analytic HESS and FREQ_CM1 already in the .ref.
   !
   subroutine check_case(nfail, case_name, hessian_tol, freq_tol)
      integer, intent(inout) :: nfail
      character(len=*), intent(in) :: case_name
      real(dp), intent(in) :: hessian_tol, freq_tol
      integer :: u, natm, ndof, n_modes_ref, sigma, mult
      integer, allocatable :: zat(:)
      real(dp), allocatable :: coords(:, :), mass(:), hess_ref(:, :), freq_ref(:)
      real(dp) :: temperature_k, pressure_pa, e_elec_ref
      real(dp), allocatable :: ztmp(:), redmass_ref(:), fc_ref(:)
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(dp), allocatable :: z_c(:), xyz_c(:)
      real(c_double), allocatable :: hessbuf(:)
      real(c_double) :: e_scf
      real(dp), allocatable :: hess(:, :)
      real(dp), allocatable :: freq(:), redmass(:), fc(:)
      real(dp), allocatable :: displacement(:, :, :)
      logical :: is_linear
      integer :: n_tr, k
      type(error_t) :: err
      real(dp) :: worst_hess, worst_asym, worst_freq, de_elec

      call ref_require('vib_thermo_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', natm)
      ndof = 3*natm
      allocate (zat(natm), ztmp(natm))
      call ref_read_vector(u, 'Z', natm, ztmp)
      zat = nint(ztmp)
      allocate (coords(3, natm))
      call ref_read_matrix(u, 'COORDS', 3, natm, coords)
      allocate (mass(natm))
      call ref_read_vector(u, 'MASS', natm, mass)
      call ref_read_int(u, 'SIGMA', sigma)
      call ref_read_int(u, 'MULT', mult)
      call ref_read_real(u, 'TEMPERATURE_K', temperature_k)
      call ref_read_real(u, 'PRESSURE_PA', pressure_pa)
      call ref_read_real(u, 'E_ELEC', e_elec_ref)
      allocate (hess_ref(ndof, ndof))
      call ref_read_matrix(u, 'HESS', ndof, ndof, hess_ref)
      call ref_read_int(u, 'N_MODES', n_modes_ref)
      allocate (freq_ref(n_modes_ref), redmass_ref(n_modes_ref), fc_ref(n_modes_ref))
      call ref_read_vector(u, 'FREQ_CM1', n_modes_ref, freq_ref)
      call ref_close(u)

      allocate (z_c(natm), xyz_c(3*natm))
      z_c = real(zat, dp)
      do k = 1, natm
         xyz_c(3*k - 2:3*k) = coords(:, k)
      end do

      rc = trc_create(h); call expect(nfail, rc == TRC_OK, case_name//': create')
      rc = trc_set_verbose(h, 0_c_int)
      rc = trc_set_molecule(h, int(natm, c_int), z_c, xyz_c, 0_c_int, int(mult, c_int))
      call expect(nfail, rc == TRC_OK, case_name//': set_molecule')
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/sto-3g.json'//c_null_char)
      call expect(nfail, rc == TRC_OK, case_name//': set_basis_json')
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      call expect(nfail, rc == TRC_OK, case_name//': set_method')
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-10_c_double, 150_c_int, 10_c_int)
      call expect(nfail, rc == TRC_OK, case_name//': set_convergence')
      rc = trc_set_screening(h, 1.0e-13_c_double)
      call expect(nfail, rc == TRC_OK, case_name//': set_screening')
      rc = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
      call expect(nfail, rc == TRC_OK, case_name//': set_guess')
      rc = trc_run_scf(h); call expect(nfail, rc == TRC_OK, case_name//': run_scf')
      rc = trc_energy(h, e_scf); call expect(nfail, rc == TRC_OK, case_name//': energy')
      de_elec = abs(real(e_scf, dp) - e_elec_ref)
      call ref_report(nfail, de_elec < 1.0e-7_dp, case_name//': E_elec vs pyscf (basis-table precision), Eh', &
                      de_elec)

      rc = trc_set_hessian_step(h, 0.005_c_double); call expect(nfail, rc == TRC_OK, case_name//': set_hessian_step')
      rc = trc_run_hessian(h); call expect(nfail, rc == TRC_OK, case_name//': run_hessian')
      allocate (hessbuf(ndof*ndof))
      rc = trc_hessian(h, hessbuf); call expect(nfail, rc == TRC_OK, case_name//': trc_hessian')
      allocate (hess(ndof, ndof))
      hess = reshape(real(hessbuf, dp), [ndof, ndof])

      ! The context must come back converged at the reference geometry.
      rc = trc_energy(h, e_scf); call expect(nfail, rc == TRC_OK, case_name//': energy after hessian')
      call ref_report(nfail, abs(real(e_scf, dp) - e_elec_ref) < 1.0e-7_dp, &
                      case_name//': E_elec after trc_run_hessian still matches pyscf (context restored), Eh', &
                      abs(real(e_scf, dp) - e_elec_ref))

      worst_asym = maxval(abs(hess - transpose(hess)))
      call ref_report(nfail, worst_asym == 0.0_dp, case_name//': H is exactly symmetric (H+H^T)/2', worst_asym)

      worst_hess = maxval(abs(hess - hess_ref))
      print '(a,a,a,es10.3)', '  check_fd_hessian_vs_pyscf: ', trim(case_name), &
         ': worst |H_FD - H_analytic| (Eh/Bohr^2) = ', worst_hess
      call ref_report(nfail, worst_hess < hessian_tol, &
                      case_name//': FD Hessian vs pyscf analytic, Eh/Bohr^2', worst_hess)

      rc = trc_destroy(h)

      ! Frequencies: terco's OWN vibrational analysis (trc_vib) on the FD
      ! Hessian just computed, against pyscf's analytic-Hessian
      ! frequencies already in the .ref (the same FREQ_CM1
      ! check_vib_thermo_vs_pyscf.F90 compares its own trc_vib run, on
      ! pyscf's Hessian, against). Any discrepancy here is the FD
      ! truncation error alone, not trc_vib's own eigensolver.
      call trc_vibrational_analysis(mass, coords, hess, freq, redmass, fc, displacement, is_linear, n_tr, err)
      call expect(nfail, .not. err%has_error(), case_name//': trc_vibrational_analysis on the FD Hessian')
      call expect(nfail, size(freq) == n_modes_ref, case_name//': FD Hessian gives the same mode count as pyscf')
      if (size(freq) == n_modes_ref) then
         worst_freq = maxval(abs(freq - freq_ref))
         print '(a,a,a,es10.3,a)', '  check_fd_hessian_vs_pyscf: ', trim(case_name), &
            ': worst |freq_FD - freq_analytic| = ', worst_freq, ' cm^-1'
         call ref_report(nfail, worst_freq < freq_tol, &
                         case_name//': FD Hessian frequencies vs pyscf analytic, cm^-1', worst_freq)
      end if
   end subroutine check_case

   !
   ! --hessian's up-front refusal: trc_run_hessian must refuse a meta-GGA
   ! BEFORE spending any SCF time on a displacement, the same way
   ! trc_run_gradient already does (gradient_capability_check, shared
   ! between the two). Checked by requiring TRC_ERR_UNSUPPORTED, the SCF
   ! iteration count unchanged (any displaced or restoring SCF would
   ! overwrite it, and a warm-started one takes far fewer iterations than
   ! the SAD start), and NO Hessian left to read (trc_hessian then reports
   ! TRC_ERR_STATE, not stale or partial data).
   !
   subroutine check_refusal(nfail)
      integer, intent(inout) :: nfail
      type(c_ptr) :: h
      integer(c_int) :: rc, it_before, it_after
      real(dp) :: z_c(3), xyz_c(9)
      real(c_double) :: hessbuf(81)

      z_c = [8.0_dp, 1.0_dp, 1.0_dp]
      xyz_c = [0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 1.8_dp, &
               1.7_dp, 0.0_dp, -0.5_dp]

      rc = trc_create(h)
      rc = trc_set_verbose(h, 0_c_int)
      rc = trc_set_molecule(h, 3_c_int, z_c, xyz_c, 0_c_int, 1_c_int)
      rc = trc_set_force_cartesian(h, 1_c_int)
      rc = trc_set_basis_json(h, 'basis_sets/sto-3g.json'//c_null_char)
      rc = trc_set_method(h, 'tpss'//c_null_char, 2_c_int)
      call expect(nfail, rc == TRC_OK, 'refusal: set_method (tpss, a runnable meta-GGA)')
      rc = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)
      rc = trc_set_screening(h, 1.0e-12_c_double)
      rc = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
      rc = trc_run_scf(h)
      call expect(nfail, rc == TRC_OK, 'refusal: run_scf (tpss energy runs fine)')

      rc = trc_iterations(h, it_before)
      rc = trc_run_hessian(h)
      call expect(nfail, rc == TRC_ERR_UNSUPPORTED, 'refusal: trc_run_hessian refuses a meta-GGA up front')
      rc = trc_iterations(h, it_after)
      call expect(nfail, it_after == it_before, 'refusal: no SCF ran before the refusal')

      rc = trc_hessian(h, hessbuf)
      call expect(nfail, rc == TRC_ERR_STATE, 'refusal: trc_hessian has nothing to report after the refusal')

      rc = trc_destroy(h)
   end subroutine check_refusal

   subroutine expect(nfail, cond, what)
      integer, intent(inout) :: nfail
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', 'check_fd_hessian_vs_pyscf: FAILED: '//what
         nfail = nfail + 1
      end if
   end subroutine expect

end module check_fd_hessian_suite

program check_fd_hessian_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_fd_hessian_suite, only: collect_fd_hessian
   implicit none
   integer :: stat
   stat = 0
   ! Serial: each case drives its own context through 6N SCFs, and their
   ! log lines would otherwise interleave.
   call run_testsuite(collect_fd_hessian, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_fd_hessian_vs_pyscf
