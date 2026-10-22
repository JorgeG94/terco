module check_constants_suite
   !! trc_constants' derived relations, checked against the OTHER CODATA
   !! 2018 values on the same file (so a mistyped literal shows up as a
   !! failed identity, not just "looks plausible"), and the two atomic-mass
   !! tables, spot-checked against pyscf 2.14.0's own arrays plus the
   !! out-of-range refusal.
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_constants, only: PI, SQRT_PI, ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM, &
                            HARTREE_TO_J, HARTREE_TO_KCALMOL, HARTREE_TO_CM1, &
                            KB_SI, KB_HARTREE, H_SI, C_SI, NA_SI, CAL_TO_J, N_ELEMENTS, &
                            trc_atomic_mass, trc_atomic_mass_average
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
