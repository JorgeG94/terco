!
! trc_vib/trc_thermo against pyscf's own analytic RHF/STO-3G Hessian,
! harmonic_analysis and thermo (test/pyscf_refs/vib_thermo_{water,hcn,
! nh3_planar}.ref, generate.py's gen_vib_thermo -- see src/scf/trc_vib.F90
! and src/scf/trc_thermo.F90's own headers for every formula: sigma and
! mult as inputs, rotational constants in GHz, imaginary modes counted
! among the true vibrations only, one mode set, and every constant from
! trc_constants.F90 (CODATA 2018).
!
! Water (sigma=2) and HCN (linear, sigma=1) are minima at this geometry/
! basis; planar NH3 (sigma=6) has one genuine imaginary umbrella mode.
! None of the three is at a stationary point of terco's own SCF -- the
! same geometry pyscf used, not terco's optimum -- which is fine for a
! like-for-like regression against a fixed, external Hessian; see the
! generator's own docstring.
!
! Also checks two properties no stored number can: a pure rigid-
! translation perturbation added to the raw Hessian leaves every
! frequency unchanged (trc_vib's TR removal is exact, not approximate --
! see trc_vib.F90's own header), HCN's linear path gives exactly
! 3*3-5 = 4 modes, and a slightly bent HCN is called linear by trc_vib
! and trc_thermo alike.
!
! A test-drive suite: each case is one unit test, which fails when any of
! its checks does; every check still prints its value through ref_report.
!
module check_vib_thermo_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use trc_error, only: error_t
   use trc_constants, only: trc_atomic_mass, BOHR_TO_ANGSTROM
   use trc_vib, only: trc_vibrational_analysis
   use trc_thermo, only: trc_thermo_result_t, trc_thermochemistry
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, &
                           ref_read_vector, ref_read_matrix, ref_read_matrix3, ref_report
   implicit none
   private
   public :: collect_vib_thermo

contains

   subroutine collect_vib_thermo(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water vs pyscf", test_water), &
                  new_unittest("hcn vs pyscf", test_hcn), &
                  new_unittest("planar nh3 vs pyscf", test_nh3_planar), &
                  new_unittest("a rigid translation leaves frequencies unchanged", test_tr_invariance), &
                  new_unittest("trc_vib and trc_thermo agree on a near-linear molecule", test_near_linear) &
                  ]
   end subroutine collect_vib_thermo

   subroutine verdict(error, nfail, what)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nfail
      character(len=*), intent(in) :: what
      character(len=12) :: buf
      if (nfail == 0) return
      write (buf, "(i0)") nfail
      call test_failed(error, what//": "//trim(buf)//" check(s) failed (see the FAIL lines above)")
   end subroutine verdict

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_case(nfail, "water")
      call verdict(error, nfail, "water")
   end subroutine test_water

   subroutine test_hcn(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_case(nfail, "hcn")
      call verdict(error, nfail, "hcn")
   end subroutine test_hcn

   subroutine test_nh3_planar(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_case(nfail, "nh3_planar")
      call verdict(error, nfail, "nh3_planar")
   end subroutine test_nh3_planar

   subroutine test_tr_invariance(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail
      nfail = 0
      call check_tr_invariance(nfail)
      call verdict(error, nfail, "tr_invariance")
   end subroutine test_tr_invariance

   !
   ! HCN with the H moved 1e-5 Angstrom off the axis, the noise a geometry
   ! pasted from another code carries. trc_vib used to keep a sixth rigid-
   ! body direction for it (a singular-value test of its own), so returned
   ! 3 modes, while trc_thermo's moment test still called it linear and
   ! applied the linear rotor: the ZPE came out one bend short. Both now
   ! decide from the same moments, so a linear call means 3N-5 modes. The
   ! Hessian is the linear one from the reference; only the geometry
   ! moves, which is all the rotor test reads.
   !
   subroutine test_near_linear(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nfail, u, natm, ndof, sigma, mult
      real(dp), allocatable :: coords(:, :), mass(:), hess(:, :), ztmp(:)
      real(dp), allocatable :: freq(:), redmass(:), fc(:), displacement(:, :, :)
      real(dp) :: temperature_k, pressure_pa, e_elec
      logical :: is_linear
      integer :: n_tr
      type(error_t) :: err
      type(trc_thermo_result_t) :: th
      nfail = 0
      call ref_require('vib_thermo_hcn.ref', u)
      call ref_read_int(u, 'NAT', natm)
      ndof = 3*natm
      allocate (coords(3, natm), mass(natm), hess(ndof, ndof), ztmp(natm))
      call ref_read_vector(u, 'Z', natm, ztmp)
      call ref_read_matrix(u, 'COORDS', 3, natm, coords)
      call ref_read_vector(u, 'MASS', natm, mass)
      call ref_read_int(u, 'SIGMA', sigma)
      call ref_read_int(u, 'MULT', mult)
      call ref_read_real(u, 'TEMPERATURE_K', temperature_k)
      call ref_read_real(u, 'PRESSURE_PA', pressure_pa)
      call ref_read_real(u, 'E_ELEC', e_elec)
      call ref_read_matrix(u, 'HESS', ndof, ndof, hess)
      call ref_close(u)
      coords(1, 1) = coords(1, 1) + 1.0e-5_dp/BOHR_TO_ANGSTROM   ! the H, off the axis
      call trc_vibrational_analysis(mass, coords, hess, freq, redmass, fc, displacement, is_linear, n_tr, err)
      call ref_report(nfail, .not. err%has_error(), 'near_linear: trc_vibrational_analysis: no error', &
                      real(err%get_code(), dp))
      call trc_thermochemistry(mass, coords, freq, e_elec, sigma, mult, th, err)
      call ref_report(nfail, .not. err%has_error(), 'near_linear: trc_thermochemistry: no error', &
                      real(err%get_code(), dp))
      call ref_report(nfail, is_linear .eqv. th%is_linear, 'near_linear: trc_vib and trc_thermo agree on linear', &
                      real(merge(1, 0, is_linear) - merge(1, 0, th%is_linear), dp))
      call ref_report(nfail, size(freq) == 3*natm - 5, 'near_linear: 3N-5 modes, as for the linear rotor', &
                      real(size(freq) - (3*natm - 5), dp))
      call verdict(error, nfail, "near_linear")
   end subroutine test_near_linear

   !
   ! One case's full pipeline: read the reference, run trc_vibrational_
   ! analysis and trc_thermochemistry on pyscf's own Hessian and masses,
   ! and compare every quantity the task asked for.
   !
   subroutine check_case(nfail, case_name)
      integer, intent(inout) :: nfail
      character(len=*), intent(in) :: case_name
      integer :: u, natm, ndof, n_modes_ref, is_linear_ref, n_imag_ref, k
      integer, allocatable :: zat(:)
      real(dp), allocatable :: coords(:, :), mass(:), hess(:, :)
      real(dp), allocatable :: freq_ref(:), redmass_ref(:), fc_ref(:)
      real(dp) :: moments_ref(3), rotconst_ref(3)
      real(dp) :: zpe_ref, de_ref, dh_ref, dg_ref, etot_ref, htot_ref, gtot_ref
      real(dp) :: strans_ref, srot_ref, svib_ref, selec_ref, stot_ref
      real(dp) :: cvtrans_ref, cvrot_ref, cvvib_ref, cvelec_ref, cvtot_ref
      integer :: sigma, mult
      real(dp) :: temperature_k, pressure_pa, e_elec
      real(dp), allocatable :: freq(:), redmass(:), fc(:)
      real(dp), allocatable :: displacement(:, :, :)
      logical :: is_linear
      integer :: n_tr
      type(error_t) :: err
      type(trc_thermo_result_t) :: th
      real(dp) :: worst, tol

      call ref_require('vib_thermo_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', natm)
      ndof = 3*natm
      allocate (zat(natm))
      block
         real(dp) :: ztmp(natm)
         call ref_read_vector(u, 'Z', natm, ztmp)
         zat = nint(ztmp)
      end block
      allocate (coords(3, natm))
      call ref_read_matrix(u, 'COORDS', 3, natm, coords)
      allocate (mass(natm))
      call ref_read_vector(u, 'MASS', natm, mass)
      call ref_read_int(u, 'SIGMA', sigma)
      call ref_read_int(u, 'MULT', mult)
      call ref_read_real(u, 'TEMPERATURE_K', temperature_k)
      call ref_read_real(u, 'PRESSURE_PA', pressure_pa)
      call ref_read_real(u, 'E_ELEC', e_elec)
      allocate (hess(ndof, ndof))
      call ref_read_matrix(u, 'HESS', ndof, ndof, hess)
      call ref_read_int(u, 'N_MODES', n_modes_ref)
      allocate (freq_ref(n_modes_ref), redmass_ref(n_modes_ref), fc_ref(n_modes_ref))
      call ref_read_vector(u, 'FREQ_CM1', n_modes_ref, freq_ref)
      call ref_read_vector(u, 'REDUCED_MASS_AMU', n_modes_ref, redmass_ref)
      call ref_read_vector(u, 'FORCE_CONST_AU', n_modes_ref, fc_ref)
      block
         real(dp) :: norm_mode_ref(3, natm, n_modes_ref)
         call ref_read_matrix3(u, 'NORM_MODE', 3, natm, n_modes_ref, norm_mode_ref)
      end block
      call ref_read_int(u, 'IS_LINEAR', is_linear_ref)
      call ref_read_int(u, 'N_IMAGINARY', n_imag_ref)
      call ref_read_vector(u, 'MOMENTS_AMUA2', 3, moments_ref)
      call ref_read_vector(u, 'ROTCONST_GHZ', 3, rotconst_ref)
      call ref_read_real(u, 'ZPE', zpe_ref)
      call ref_read_real(u, 'DELTA_E', de_ref)
      call ref_read_real(u, 'DELTA_H', dh_ref)
      call ref_read_real(u, 'DELTA_G', dg_ref)
      call ref_read_real(u, 'E_TOT', etot_ref)
      call ref_read_real(u, 'H_TOT', htot_ref)
      call ref_read_real(u, 'G_TOT', gtot_ref)
      call ref_read_real(u, 'S_TRANS_CALMOLK', strans_ref)
      call ref_read_real(u, 'S_ROT_CALMOLK', srot_ref)
      call ref_read_real(u, 'S_VIB_CALMOLK', svib_ref)
      call ref_read_real(u, 'S_ELEC_CALMOLK', selec_ref)
      call ref_read_real(u, 'S_TOT_CALMOLK', stot_ref)
      call ref_read_real(u, 'CV_TRANS_CALMOLK', cvtrans_ref)
      call ref_read_real(u, 'CV_ROT_CALMOLK', cvrot_ref)
      call ref_read_real(u, 'CV_VIB_CALMOLK', cvvib_ref)
      call ref_read_real(u, 'CV_ELEC_CALMOLK', cvelec_ref)
      call ref_read_real(u, 'CV_TOT_CALMOLK', cvtot_ref)
      call ref_close(u)

      ! Masses: terco's own default table (trc_atomic_mass, the
      ! most-abundant isotope) must reproduce exactly what the generator
      ! passed to pyscf -- both are read from the SAME pyscf 2.14 table
      ! (trc_constants.F90's own header), so this should agree to several
      ! more digits than a vibrational frequency needs.
      block
         real(dp) :: mdiff
         integer :: a
         mdiff = 0.0_dp
         do a = 1, natm
            mdiff = max(mdiff, abs(trc_atomic_mass(zat(a), err) - mass(a)))
         end do
         call ref_report(nfail, mdiff < 1.0e-10_dp, case_name//': trc_atomic_mass matches the masses pyscf used', &
                         mdiff)
      end block

      call trc_vibrational_analysis(mass, coords, hess, freq, redmass, fc, displacement, is_linear, n_tr, err)
      call ref_report(nfail, .not. err%has_error(), case_name//': trc_vibrational_analysis: no error', &
                      real(err%get_code(), dp))
      call ref_report(nfail, size(freq) == n_modes_ref, case_name//': n_modes matches pyscf', &
                      real(size(freq) - n_modes_ref, dp))
      call ref_report(nfail, merge(1, 0, is_linear) == is_linear_ref, case_name//': is_linear matches pyscf', &
                      real(merge(1, 0, is_linear) - is_linear_ref, dp))

      ! Frequencies: same Hessian, same masses, both LAPACK dsyev-family
      ! eigensolvers underneath (pyscf's numpy.linalg.eigh, terco's
      ! pic_syev) -- the only differences left are eigensolver roundoff
      ! and the ~1e-9-relative CODATA-2014-vs-2018 difference in the
      ! cm^-1 conversion factor (trc_constants.F90's
      ! HARTREE_PER_AMU_BOHR2_TO_CM1 header): 1e-9 x 4800 cm^-1 = 5e-6,
      ! measured 4.9e-6 at worst.
      tol = 2.0e-5_dp
      worst = maxval(abs(freq - freq_ref))
      call ref_report(nfail, worst < tol, case_name//': frequencies vs pyscf, cm^-1', worst)
      worst = maxval(abs(redmass - redmass_ref))
      call ref_report(nfail, worst < 1.0e-8_dp, case_name//': reduced masses vs pyscf, amu', worst)
      worst = maxval(abs(fc - fc_ref))
      call ref_report(nfail, worst < 1.0e-8_dp, case_name//': force constants vs pyscf, Eh/a0^2', worst)

      call trc_thermochemistry(mass, coords, freq, e_elec, sigma, mult, th, err, &
                                temperature_k=temperature_k, pressure_pa=pressure_pa)
      call ref_report(nfail, .not. err%has_error(), case_name//': trc_thermochemistry: no error', &
                      real(err%get_code(), dp))
      call ref_report(nfail, th%n_imaginary == n_imag_ref, case_name// &
                      ': n_imaginary (among true vibrations only) matches pyscf', &
                      real(th%n_imaginary - n_imag_ref, dp))

      worst = maxval(abs(th%moments_of_inertia_amua2 - moments_ref))
      call ref_report(nfail, worst < 1.0e-6_dp, case_name//': moments of inertia vs pyscf, amu*Ang^2', worst)
      do k = 1, 3
         ! pyscf's own rotation_const divides by the moment, so a linear
         ! molecule's near-zero first moment comes back as +Infinity
         ! there (and, after the generator's %.17e round trip, here too);
         ! trc_thermo.F90 reports 0 for that same direction instead (this
         ! module's own header). Neither is a rotational constant a real
         ! molecule has, so there is nothing to compare for that index.
         if (rotconst_ref(k) > 0.0_dp .and. rotconst_ref(k) < 1.0e8_dp) then
            worst = abs(th%rot_const_ghz(k) - rotconst_ref(k))
            call ref_report(nfail, worst < 1.0e-6_dp*max(1.0_dp, abs(rotconst_ref(k))), &
                            case_name//': rotational constant vs pyscf, GHz', worst)
         end if
      end do

      ! ZPE carries the same 1e-9-relative cm^-1 factor: ~3e-11 Eh at 0.025 Eh
      ! (measured 4e-12).
      call ref_report(nfail, abs(th%zpe_hartree - zpe_ref) < 1.0e-10_dp, case_name//': ZPE vs pyscf, Eh', &
                      abs(th%zpe_hartree - zpe_ref))
      call ref_report(nfail, abs(th%delta_e_hartree - de_ref) < 1.0e-8_dp, case_name//': delta E vs pyscf, Eh', &
                      abs(th%delta_e_hartree - de_ref))
      call ref_report(nfail, abs(th%delta_h_hartree - dh_ref) < 1.0e-8_dp, case_name//': delta H vs pyscf, Eh', &
                      abs(th%delta_h_hartree - dh_ref))
      call ref_report(nfail, abs(th%delta_g_hartree - dg_ref) < 1.0e-8_dp, case_name//': delta G vs pyscf, Eh', &
                      abs(th%delta_g_hartree - dg_ref))
      call ref_report(nfail, abs(th%e_tot_hartree - etot_ref) < 1.0e-8_dp, case_name//': E_tot vs pyscf, Eh', &
                      abs(th%e_tot_hartree - etot_ref))
      call ref_report(nfail, abs(th%h_tot_hartree - htot_ref) < 1.0e-8_dp, case_name//': H_tot vs pyscf, Eh', &
                      abs(th%h_tot_hartree - htot_ref))
      call ref_report(nfail, abs(th%g_tot_hartree - gtot_ref) < 1.0e-8_dp, case_name//': G_tot vs pyscf, Eh', &
                      abs(th%g_tot_hartree - gtot_ref))

      ! S_trans/S_tot/Cv_trans/Cv_rot carry a genuine, QUANTIFIED
      ! CODATA-year difference that no amount of agreement elsewhere
      ! removes: this module's R_CALMOLK/KB_SI/AMU_TO_KG/H_SI are CODATA
      ! 2018 (trc_constants.F90); the generator's pyscf reference used
      ! pyscf's own pyscf.data.nist module, whose BOLTZMANN/PLANCK/
      ! ATOMIC_MASS/AVOGADRO are CODATA-2014-era. Measured directly (both
      ! constant sets, same formula, same geometry): R_CALMOLK itself
      ! differs by 6.59e-7 (3.3e-7 relative), which alone gives
      ! 1.5*that = 9.88e-7 for Cv_trans/Cv_rot (confirmed: every case
      ! above measures exactly 9.88e-7 on both), and S_trans's
      ! ln(q_trans) term, propagating the SAME constants through
      ! (kB*T/P)*(2 pi m kB T/h^2)^1.5, differs by 1.31-1.36e-5 cal/mol/K
      ! (confirmed by recomputing q_trans with both constant sets in
      ! Python).
      ! Tolerances below sit at roughly 1.5x the measured figure, not
      ! loosened to paper over a mismatch.
      call ref_report(nfail, abs(th%s_trans_calmolk - strans_ref) < 2.0e-5_dp, &
                      case_name//': S_trans vs pyscf, cal/mol/K (CODATA-year, see above)', &
                      abs(th%s_trans_calmolk - strans_ref))
      call ref_report(nfail, abs(th%s_rot_calmolk - srot_ref) < 1.0e-5_dp, &
                      case_name//': S_rot vs pyscf, cal/mol/K', abs(th%s_rot_calmolk - srot_ref))
      call ref_report(nfail, abs(th%s_vib_calmolk - svib_ref) < 1.0e-5_dp, &
                      case_name//': S_vib vs pyscf, cal/mol/K', abs(th%s_vib_calmolk - svib_ref))
      call ref_report(nfail, abs(th%s_elec_calmolk - selec_ref) < 1.0e-8_dp, &
                      case_name//': S_elec vs pyscf, cal/mol/K', abs(th%s_elec_calmolk - selec_ref))
      call ref_report(nfail, abs(th%s_tot_calmolk - stot_ref) < 3.0e-5_dp, &
                      case_name//': S_tot vs pyscf, cal/mol/K (CODATA-year, see above)', &
                      abs(th%s_tot_calmolk - stot_ref))
      call ref_report(nfail, abs(th%cv_trans_calmolk - cvtrans_ref) < 2.0e-6_dp, &
                      case_name//': Cv_trans vs pyscf, cal/mol/K (CODATA-year, see above)', &
                      abs(th%cv_trans_calmolk - cvtrans_ref))
      call ref_report(nfail, abs(th%cv_rot_calmolk - cvrot_ref) < 2.0e-6_dp, &
                      case_name//': Cv_rot vs pyscf, cal/mol/K (CODATA-year, see above)', &
                      abs(th%cv_rot_calmolk - cvrot_ref))
      call ref_report(nfail, abs(th%cv_vib_calmolk - cvvib_ref) < 1.0e-5_dp, &
                      case_name//': Cv_vib vs pyscf, cal/mol/K', abs(th%cv_vib_calmolk - cvvib_ref))
      call ref_report(nfail, abs(th%cv_elec_calmolk - cvelec_ref) < 1.0e-8_dp, &
                      case_name//': Cv_elec vs pyscf, cal/mol/K', abs(th%cv_elec_calmolk - cvelec_ref))
      call ref_report(nfail, abs(th%cv_tot_calmolk - cvtot_ref) < 3.0e-6_dp, &
                      case_name//': Cv_tot vs pyscf, cal/mol/K (CODATA-year, see above)', &
                      abs(th%cv_tot_calmolk - cvtot_ref))

      if (case_name == 'hcn') then
         call ref_report(nfail, size(freq) == 3*natm - 5, 'hcn: linear path gives 3N-5 modes', &
                         real(size(freq) - (3*natm - 5), dp))
      end if
      if (case_name == 'nh3_planar') then
         call ref_report(nfail, th%n_imaginary == 1, 'nh3_planar: exactly one imaginary mode', &
                         real(th%n_imaginary - 1, dp))
         call ref_report(nfail, count(freq < 0.0_dp) == 1 .and. freq(1) < -1000.0_dp, &
                         'nh3_planar: the imaginary mode is the umbrella one (about -1118 cm^-1)', freq(1))
      end if
   end subroutine check_case

   !
   ! A pure rigid-translation perturbation, added to the RAW (un-mass-
   ! weighted) Hessian, lies entirely inside the TR null space trc_vib
   ! projects out -- see this program's own header -- so it must leave
   ! every vibrational frequency unchanged, to the eigensolver's own
   ! precision, not merely "close". Built directly from water's own
   ! reference Hessian and masses.
   !
   subroutine check_tr_invariance(nfail)
      integer, intent(inout) :: nfail
      integer :: u, natm, ndof, pa, pb, sigma, mult
      real(dp), allocatable :: coords(:, :), mass(:), hess(:, :), hess_pert(:, :), ztmp(:)
      real(dp), allocatable :: freq(:), freq_pert(:), redmass(:), fc(:)
      real(dp), allocatable :: displacement(:, :, :)
      logical :: is_linear
      integer :: n_tr
      type(error_t) :: err
      real(dp), parameter :: ALPHA = 1.0e-3_dp
      real(dp) :: worst, temperature_k, pressure_pa, e_elec

      ! pyscf_ref_io is a strictly sequential reader (one tag after
      ! another, see that module's own header) -- every tag ahead of the
      ! two this check needs (COORDS, HESS) must still be consumed, in
      ! the order the generator wrote them, even though their values are
      ! unused here.
      call ref_require('vib_thermo_water.ref', u)
      call ref_read_int(u, 'NAT', natm)
      ndof = 3*natm
      allocate (coords(3, natm), mass(natm), hess(ndof, ndof), ztmp(natm))
      call ref_read_vector(u, 'Z', natm, ztmp)
      call ref_read_matrix(u, 'COORDS', 3, natm, coords)
      call ref_read_vector(u, 'MASS', natm, mass)
      call ref_read_int(u, 'SIGMA', sigma)
      call ref_read_int(u, 'MULT', mult)
      call ref_read_real(u, 'TEMPERATURE_K', temperature_k)
      call ref_read_real(u, 'PRESSURE_PA', pressure_pa)
      call ref_read_real(u, 'E_ELEC', e_elec)
      call ref_read_matrix(u, 'HESS', ndof, ndof, hess)
      call ref_close(u)

      allocate (hess_pert(ndof, ndof))
      hess_pert = hess
      do pb = 1, natm
         do pa = 1, natm
            ! x component of atom pa, pb: row/col 3*(a-1)+1
            hess_pert(3*(pa - 1) + 1, 3*(pb - 1) + 1) = hess_pert(3*(pa - 1) + 1, 3*(pb - 1) + 1) + &
                                                         ALPHA*mass(pa)*mass(pb)
         end do
      end do

      call trc_vibrational_analysis(mass, coords, hess, freq, redmass, fc, displacement, is_linear, n_tr, err)
      call ref_report(nfail, .not. err%has_error(), 'tr_invariance: baseline run: no error', &
                      real(err%get_code(), dp))
      call trc_vibrational_analysis(mass, coords, hess_pert, freq_pert, redmass, fc, displacement, is_linear, &
                                     n_tr, err)
      call ref_report(nfail, .not. err%has_error(), 'tr_invariance: perturbed run: no error', &
                      real(err%get_code(), dp))

      worst = maxval(abs(freq - freq_pert))
      call ref_report(nfail, worst < 1.0e-6_dp, &
                      'tr_invariance: a pure translation added to the Hessian leaves frequencies unchanged', worst)
   end subroutine check_tr_invariance

end module check_vib_thermo_suite

program check_vib_thermo_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_vib_thermo_suite, only: collect_vib_thermo
   implicit none
   integer :: stat
   stat = 0
   ! Serial, so the cases' ref_report lines do not interleave.
   call run_testsuite(collect_vib_thermo, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_vib_thermo_vs_pyscf
