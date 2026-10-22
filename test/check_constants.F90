module check_constants_suite
   !! trc_constants' derived relations, checked against the OTHER CODATA
   !! 2018 values on the same file (so a mistyped literal shows up as a
   !! failed identity, not just "looks plausible"), and the two atomic-mass
   !! tables, spot-checked against pyscf 2.14.0's own arrays plus the
   !! out-of-range refusal.
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_constants, only: PI, TWO_PI, SQRT_PI, ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM, &
                            HARTREE_TO_J, HARTREE_TO_KCALMOL, HARTREE_TO_CM1, HARTREE_TO_CALMOL, &
                            KB_SI, KB_HARTREE, H_SI, C_SI, NA_SI, CAL_TO_J, AMU_TO_KG, N_ELEMENTS, &
                            CM1_TO_KELVIN, HARTREE_PER_AMU_BOHR2_TO_CM1, ROTTEMP_AMUA2_TO_K, &
                            ROTCONST_AMUA2_TO_GHZ, trc_atomic_mass, trc_atomic_mass_average
   implicit none
   private
   public :: collect_constants

contains

   subroutine collect_constants(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("angstrom/bohr round trip to one", test_bohr_round_trip), &
                  new_unittest("sqrt(pi) matches pi", test_sqrt_pi), &
                  new_unittest("kB in Hartree matches kB in J / Hartree in J", test_kb_hartree), &
                  new_unittest("Hartree->kcal/mol matches Hartree->J * NA / 4184", test_hartree_kcalmol), &
                  new_unittest("Hartree->cm-1 matches Hartree->J / (h c 100)", test_hartree_cm1), &
                  new_unittest("Hartree->cal/mol matches Hartree->kcal/mol * 1000", test_hartree_calmol), &
                  new_unittest("cm-1 -> K matches h c 100 / kB", test_cm1_to_kelvin), &
                  new_unittest("vib. freq. coefficient matches sqrt(Eh/(amu a0^2))/(2 pi c) 1e-2", &
                               test_freq_coeff), &
                  new_unittest("rotational-temperature coefficient matches h^2/(8 pi^2 amu Ang^2 kB)", &
                               test_rottemp), &
                  new_unittest("rotational-constant coefficient matches h/(8 pi^2 amu Ang^2) in GHz", &
                               test_rotconst_ghz), &
                  new_unittest("isotope masses: H, C, O, P", test_isotope_masses), &
                  new_unittest("standard atomic weights: H, C, O, P", test_average_masses), &
                  new_unittest("an atomic number outside [1, N_ELEMENTS] is refused", test_bad_z) &
                  ]
   end subroutine collect_constants

   subroutine test_bohr_round_trip(error)
      !! ANGSTROM_TO_BOHR is 1/BOHR_TO_ANGSTROM: the product must be 1 to
      !! a couple of ULPs, not merely "close".
      type(error_type), allocatable, intent(out) :: error
      call check(error, ANGSTROM_TO_BOHR*BOHR_TO_ANGSTROM, 1.0_dp, thr=4.0_dp*epsilon(1.0_dp))
   end subroutine test_bohr_round_trip

   subroutine test_sqrt_pi(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, SQRT_PI*SQRT_PI, PI, thr=1.0e-12_dp, rel=.true.)
   end subroutine test_sqrt_pi

   subroutine test_kb_hartree(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, KB_HARTREE, KB_SI/HARTREE_TO_J, thr=1.0e-12_dp, rel=.true.)
   end subroutine test_kb_hartree

   subroutine test_hartree_kcalmol(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, HARTREE_TO_KCALMOL, HARTREE_TO_J*NA_SI/(1000.0_dp*CAL_TO_J), &
                 thr=1.0e-10_dp, rel=.true.)
   end subroutine test_hartree_kcalmol

   subroutine test_hartree_cm1(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, HARTREE_TO_CM1, HARTREE_TO_J/(H_SI*C_SI*100.0_dp), thr=1.0e-10_dp, rel=.true.)
   end subroutine test_hartree_cm1

   subroutine test_hartree_calmol(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, HARTREE_TO_CALMOL, HARTREE_TO_KCALMOL*1000.0_dp, thr=1.0e-12_dp, rel=.true.)
   end subroutine test_hartree_calmol

   subroutine test_cm1_to_kelvin(error)
      !! trc_vib.F90/trc_thermo.F90's own vibrational-temperature factor:
      !! theta_v = CM1_TO_KELVIN * wavenumber.
      type(error_type), allocatable, intent(out) :: error
      call check(error, CM1_TO_KELVIN, H_SI*C_SI*1.0e2_dp/KB_SI, thr=1.0e-10_dp, rel=.true.)
   end subroutine test_cm1_to_kelvin

   subroutine test_freq_coeff(error)
      !! trc_vib.F90's own cm^-1 conversion: freq_cm1 = sign(lambda) *
      !! sqrt(|lambda|) * HARTREE_PER_AMU_BOHR2_TO_CM1, lambda a
      !! mass-weighted-Hessian eigenvalue in Hartree/(amu bohr^2).
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: expect
      expect = sqrt(HARTREE_TO_J/(AMU_TO_KG*(BOHR_TO_ANGSTROM*1.0e-10_dp)**2))/(TWO_PI*C_SI)*1.0e-2_dp
      call check(error, HARTREE_PER_AMU_BOHR2_TO_CM1, expect, thr=1.0e-10_dp, rel=.true.)
      ! And the value itself: CODATA 2018 gives ~5140.487 (mqc's own ad hoc
      ! sqrt(AU_TO_CM1) = 5140.48733 cm^-1, 3.6e-8 relative higher --
      ! HESSIAN_VIB_THERMO_SPEC.md section 3.4).
      if (allocated(error)) return
      call check(error, HARTREE_PER_AMU_BOHR2_TO_CM1, 5140.487_dp, thr=1.0e-3_dp, rel=.true.)
   end subroutine test_freq_coeff

   subroutine test_rottemp(error)
      !! mqc's defect D5: its own ROTTEMP_AMUA2_TO_K literal is 24.2637,
      !! wrong in the 4th digit; this one evaluates to 24.2544.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: expect
      expect = H_SI**2/(8.0_dp*PI**2*AMU_TO_KG*(1.0e-10_dp)**2*KB_SI)
      call check(error, ROTTEMP_AMUA2_TO_K, expect, thr=1.0e-10_dp, rel=.true.)
      if (allocated(error)) return
      call check(error, ROTTEMP_AMUA2_TO_K, 24.2544_dp, thr=1.0e-4_dp, rel=.true.)
   end subroutine test_rottemp

   subroutine test_rotconst_ghz(error)
      !! mqc's defect D3: its own ROTCONST_AMUA2_TO_GHZ literal (505379.07)
      !! is the MHz coefficient mislabelled GHz; this one evaluates to
      !! 505.379, the actual GHz value.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: expect
      expect = H_SI/(8.0_dp*PI**2*AMU_TO_KG*(1.0e-10_dp)**2)*1.0e-9_dp
      call check(error, ROTCONST_AMUA2_TO_GHZ, expect, thr=1.0e-10_dp, rel=.true.)
      if (allocated(error)) return
      call check(error, ROTCONST_AMUA2_TO_GHZ, 505.379_dp, thr=1.0e-4_dp, rel=.true.)
   end subroutine test_rotconst_ghz

   subroutine test_isotope_masses(error)
      !! pyscf 2.14.0's COMMON_ISOTOPE_MASSES[1, 6, 8, 15].
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      call check(error, trc_atomic_mass(1, err), 1.007825_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, trc_atomic_mass(6, err), 12.0_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, trc_atomic_mass(8, err), 15.994915_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, trc_atomic_mass(15, err), 30.973762_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, .not. err%has_error(), "a Z in range must not raise an error")
   end subroutine test_isotope_masses

   subroutine test_average_masses(error)
      !! pyscf 2.14.0's MASSES[1, 6, 8, 15].
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      call check(error, trc_atomic_mass_average(1, err), 1.008_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, trc_atomic_mass_average(6, err), 12.011_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, trc_atomic_mass_average(8, err), 15.999_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, trc_atomic_mass_average(15, err), 30.973761998_dp, thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, .not. err%has_error(), "a Z in range must not raise an error")
   end subroutine test_average_masses

   subroutine test_bad_z(error)
      !! Z = 0 and Z past N_ELEMENTS are both refused, not silently given
      !! whatever garbage sits past the table's end.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: mass

      mass = trc_atomic_mass(0, err)
      call check(error, err%has_error(), "Z = 0 must be refused")
      if (allocated(error)) return

      call err%clear()
      mass = trc_atomic_mass(N_ELEMENTS + 1, err)
      call check(error, err%has_error(), "Z past N_ELEMENTS must be refused")
      if (allocated(error)) return

      call err%clear()
      mass = trc_atomic_mass_average(-5, err)
      call check(error, err%has_error(), "a negative Z must be refused")
      if (allocated(error)) return

      ! Silence "set but never used": the refusal path is what is tested,
      ! not the value it leaves in `mass`.
      call check(error, mass >= 0.0_dp, "the refused call still returns a finite value")
   end subroutine test_bad_z

end module check_constants_suite

program check_constants
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_constants_suite, only: collect_constants
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_constants, error_unit, stat)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_constants
