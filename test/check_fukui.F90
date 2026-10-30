!
! Hirshfeld-condensed Fukui indices (trc_capi's trc_fukui_hirshfeld) against
! metalquicha's own fukui_indices, at its own test_fukui case.
!
! mqc's test/test_mqc_hirshfeld.f90 (test_fukui) runs RHF/cc-pVDZ water at
! the Bohr geometry below, calls fukui_indices("hirshfeld", ...) with no
! neutral orbitals offered to seed the ions, and pins its own F_PLUS/F_MINUS
! there. "No neutral orbitals offered" is mqc's unseeded path, which starts
! each ion from GWH; trc_fukui_hirshfeld instead starts the ions from the
! handle's own guess (SAD here, so SADQ for the ions -- see its header). Both
! are unseeded and converge to the same doublet states, which is what
! F_PLUS/F_MINUS pin, so this is terco's actual default against mqc's.
!
! Also: the sum rule, the handle's state after the call (on success and when
! an ion fails to converge), and the refusal of an open-shell neutral.
!
module check_fukui_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_bind_device
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_method, trc_set_convergence, trc_run_scf, trc_energy, &
                               trc_nspin, trc_fukui_hirshfeld, trc_set_lindep, trc_set_scf_accel, &
                               TRC_OK, TRC_ERR_BADARG, TRC_ERR_NOCONV
   implicit none
   private
   public :: collect_fukui_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   real(dp), parameter :: WATER_R(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])
      !! Bohr, mqc's test_mqc_hirshfeld.f90's WATER (its test_fukui and
      !! test_mqc_water_neutral/test_mqc_water_cation geometry) -- the same
      !! array test/check_hirshfeld.F90's MQC_WATER_R already uses, passed
      !! straight to trc_set_molecule with no Angstrom conversion, exactly
      !! as mqc's own build_czt_molecule call takes it.

   real(dp), parameter :: F_PLUS(3) = [0.304597030610_dp, 0.347701484691_dp, &
                                       0.347701484691_dp]
   real(dp), parameter :: F_MINUS(3) = [0.669189280600_dp, 0.165405359742_dp, &
                                        0.165405359742_dp]
      !! mqc's own test_mqc_hirshfeld.f90's test_fukui, verbatim: RHF/
      !! cc-pVDZ, 10 electrons, closed-shell neutral, max_iteration=100,
      !! energy_tol=1e-11, density_tol=1e-9, unseeded ions.

   real(dp), parameter :: CHARGE_TOL = 1.0e-7_dp
      !! mqc's own CHARGE_TOL (test_mqc_hirshfeld.f90), for the same reason
      !! test/check_hirshfeld.F90's MQC_CHARGE_TOL uses it: two independent
      !! SCF implementations (terco direct, mqc through its CPU backend)
      !! converging the same three states and condensing each through a
      !! Hirshfeld partition that is the same RECIPE but not the same
      !! CODE (terco's cubic-Hermite proatom table + trc_dft_grid vs mqc's
      !! own radial/angular tables and run-time free-atom solve). Measured
      !! (this file, forcing CHARGE_TOL to 1e-12 to print the gap, as that
      !! file's own comment does): worst |f+ - F_PLUS| = 1.568e-9, worst
      !! |f- - F_MINUS| = 4.512e-10 -- both far inside mqc's own 1e-7 and
      !! the same order of magnitude as check_hirshfeld.F90's own
      !! MQC_CHARGE_TOL agreement on the neutral and cation alone, so the
      !! extra anion SCF and the finite difference are not costing much
      !! beyond what the neutral/cation comparison already pays.
   real(dp), parameter :: SUM_TOL = 1.0e-9_dp
      !! f+ and f- each sum to one electron moved; trc_fukui_indices's own
      !! check_sum_rule already refuses a run outside 1e-6, so this is a
      !! tighter, separate check on the numbers trc_fukui_hirshfeld handed
      !! back, round-off above that internal gate.

contains

   subroutine collect_fukui_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water matches mqc's F_PLUS/F_MINUS", test_water_fukui), &
                  new_unittest("f+ and f- each sum to one electron", test_sum_rule), &
                  new_unittest("the handle's charge and multiplicity are restored", test_restore), &
                  new_unittest("an ion that does not converge is reported and the handle restored", &
                               test_ion_noconv), &
                  new_unittest("a non-singlet neutral is refused", test_radical_refused), &
                  new_unittest("trc_set_lindep refuses NaN and out-of-range", test_lindep_refused), &
                  new_unittest("trc_set_scf_accel refuses NaN and out-of-range", test_scf_accel_refused) &
                  ]
   end subroutine collect_fukui_tests

   subroutine run_water(h, charge, multiplicity, error)
      !! A converged RHF (or UHF, for test_radical_refused) handle at
      !! WATER_R/cc-pVDZ, mqc's own convergence for this case.
      type(c_ptr), intent(out) :: h
      integer, intent(in) :: charge, multiplicity
      type(error_type), allocatable, intent(out) :: error

      integer(c_int) :: rc

      rc = trc_create(h)
      rc = trc_set_molecule(h, 3_c_int, real(WATER_Z, c_double), reshape(WATER_R, [9]), &
                            int(charge, c_int), int(multiplicity, c_int))
      call check(error, rc == TRC_OK, "trc_set_molecule failed")
      if (allocated(error)) return
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      call check(error, rc == TRC_OK, "trc_set_basis_json failed")
      if (allocated(error)) return
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      call check(error, rc == TRC_OK, "trc_set_method failed")
      if (allocated(error)) return
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 10_c_int)
      call check(error, rc == TRC_OK, "trc_set_convergence failed")
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "trc_run_scf (neutral) did not converge")
   end subroutine run_water

   subroutine test_water_fukui(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: fplus_c(3), fminus_c(3), fzero_c(3), qn_c(3), qnp1_c(3), qnm1_c(3)
      real(dp) :: worst_plus, worst_minus
      integer :: ia

      call run_water(h, 0, 1, error)
      if (allocated(error)) return

      rc = trc_fukui_hirshfeld(h, fplus_c, fminus_c, fzero_c, qn_c, qnp1_c, qnm1_c)
      call check(error, rc == TRC_OK, "trc_fukui_hirshfeld failed")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      worst_plus = 0.0_dp
      worst_minus = 0.0_dp
      do ia = 1, 3
         worst_plus = max(worst_plus, abs(real(fplus_c(ia), dp) - F_PLUS(ia)))
         worst_minus = max(worst_minus, abs(real(fminus_c(ia), dp) - F_MINUS(ia)))
      end do
      call check(error, worst_plus < CHARGE_TOL, "f+: worst |f_terco - F_PLUS| = "// &
                to_text(worst_plus)//" exceeds mqc's own CHARGE_TOL = "//to_text(CHARGE_TOL))
      if (.not. allocated(error)) then
         call check(error, worst_minus < CHARGE_TOL, "f-: worst |f_terco - F_MINUS| = "// &
                   to_text(worst_minus)//" exceeds mqc's own CHARGE_TOL = "//to_text(CHARGE_TOL))
      end if
      if (.not. allocated(error)) then
         call check(error, all(fzero_c > 0.0_c_double), "f0 must be positive for water (mqc's own observation)")
      end if
      rc = trc_destroy(h)
   end subroutine test_water_fukui

   subroutine test_sum_rule(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: fplus_c(3), fminus_c(3), fzero_c(3), qn_c(3), qnp1_c(3), qnm1_c(3)

      call run_water(h, 0, 1, error)
      if (allocated(error)) return

      rc = trc_fukui_hirshfeld(h, fplus_c, fminus_c, fzero_c, qn_c, qnp1_c, qnm1_c)
      call check(error, rc == TRC_OK, "trc_fukui_hirshfeld failed")
      if (.not. allocated(error)) then
         call check(error, abs(sum(real(fplus_c, dp)) - 1.0_dp), 0.0_dp, thr=SUM_TOL, &
                   message="f+ does not sum to one electron")
      end if
      if (.not. allocated(error)) then
         call check(error, abs(sum(real(fminus_c, dp)) - 1.0_dp), 0.0_dp, thr=SUM_TOL, &
                   message="f- does not sum to one electron")
      end if
      rc = trc_destroy(h)
   end subroutine test_sum_rule

   subroutine test_restore(error)
      !! trc_fukui_hirshfeld must leave the handle converged on the SAME
      !! closed-shell neutral it was called with: nspin back to 1 (a
      !! doublet ion left on the handle would read nspin 2) and the energy
      !! back to what it was before the call, to the re-convergence this
      !! context's own restore does (trc_run_hessian's pattern, copied).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, nspin_before, nspin_after
      real(c_double) :: e_before, e_after
      real(c_double) :: fplus_c(3), fminus_c(3), fzero_c(3), qn_c(3), qnp1_c(3), qnm1_c(3)

      call run_water(h, 0, 1, error)
      if (allocated(error)) return
      rc = trc_nspin(h, nspin_before)
      rc = trc_energy(h, e_before)

      rc = trc_fukui_hirshfeld(h, fplus_c, fminus_c, fzero_c, qn_c, qnp1_c, qnm1_c)
      call check(error, rc == TRC_OK, "trc_fukui_hirshfeld failed")
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if

      rc = trc_nspin(h, nspin_after)
      call check(error, nspin_after == nspin_before, "nspin was not restored to the neutral's")
      if (.not. allocated(error)) then
         rc = trc_energy(h, e_after)
         call check(error, abs(real(e_after, dp) - real(e_before, dp)), 0.0_dp, thr=1.0e-9_dp, &
                   message="the neutral's energy was not restored after the ion SCFs")
      end if
      ! And the handle still runs an ordinary SCF afterwards, at the
      ! charge/multiplicity it was created with.
      if (.not. allocated(error)) then
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "the handle could not run another SCF after trc_fukui_hirshfeld")
      end if
      rc = trc_destroy(h)
   end subroutine test_restore

   subroutine test_ion_noconv(error)
      !! An ion SCF that runs out of iterations must come back as
      !! TRC_ERR_NOCONV, and the handle must still be the closed-shell
      !! neutral afterwards: with the iteration limit lifted, an ordinary
      !! SCF reconverges to the neutral's own energy with nspin 1. Two
      !! iterations cannot converge an ion from an atomic guess at 1e-12.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc, nspin_after
      real(c_double) :: e_before, e_after
      real(c_double) :: fplus_c(3), fminus_c(3), fzero_c(3), qn_c(3), qnp1_c(3), qnm1_c(3)

      call run_water(h, 0, 1, error)
      if (allocated(error)) return
      rc = trc_energy(h, e_before)

      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 2_c_int, 10_c_int)
      rc = trc_fukui_hirshfeld(h, fplus_c, fminus_c, fzero_c, qn_c, qnp1_c, qnm1_c)
      call check(error, rc == TRC_ERR_NOCONV, "an ion out of iterations was not reported as "// &
                 "TRC_ERR_NOCONV")
      if (.not. allocated(error)) then
         rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 10_c_int)
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "the handle could not reconverge the neutral")
      end if
      if (.not. allocated(error)) then
         rc = trc_nspin(h, nspin_after)
         call check(error, nspin_after == 1, "the handle was left open-shell (an ion's state)")
      end if
      if (.not. allocated(error)) then
         rc = trc_energy(h, e_after)
         call check(error, abs(real(e_after, dp) - real(e_before, dp)), 0.0_dp, thr=1.0e-9_dp, &
                    message="the handle was left at another charge: its SCF reconverged elsewhere")
      end if
      rc = trc_destroy(h)
   end subroutine test_ion_noconv

   subroutine test_radical_refused(error)
      !! fukui_indices's own rule, copied verbatim: a neutral that is not a
      !! closed-shell singlet is refused outright, because its ions are
      !! not both doublets and nothing here guesses what they should be.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: fplus_c(3), fminus_c(3), fzero_c(3), qn_c(3), qnp1_c(3), qnm1_c(3)

      call run_water(h, 0, 3, error)
      if (allocated(error)) return

      rc = trc_fukui_hirshfeld(h, fplus_c, fminus_c, fzero_c, qn_c, qnp1_c, qnm1_c)
      call check(error, rc == TRC_ERR_BADARG, "a triplet neutral must be refused, not computed")
      rc = trc_destroy(h)
   end subroutine test_radical_refused

   function to_text(x) result(s)
      real(dp), intent(in) :: x
      character(len=32) :: s
      write (s, "(es12.4)") x
   end function to_text

   subroutine test_lindep_refused(error)
      !! The linear-dependence threshold Fiona's anion needs (--lindep 5e-6)
      !! goes through trc_set_lindep, which must refuse anything outside
      !! [0, 1) -- NaN included, which a plain `< 0 .or. >= 1` test lets by.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: nan

      nan = ieee_value(nan, ieee_quiet_nan)
      rc = trc_create(h)
      rc = trc_set_lindep(h, nan)
      call check(error, rc == TRC_ERR_BADARG, "trc_set_lindep accepted NaN")
      if (.not. allocated(error)) then
         rc = trc_set_lindep(h, -1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_lindep accepted -1")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_lindep(h, 1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_lindep accepted 1")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_lindep(h, 5.0e-6_c_double)
         call check(error, rc == TRC_OK, "trc_set_lindep refused 5e-6")
      end if
      rc = trc_destroy(h)
   end subroutine test_lindep_refused

   subroutine test_scf_accel_refused(error)
      !! trc_set_scf_accel's four refusals (trc_capi.F90): adiis outside
      !! {0,1}, ndamp < -1, damp outside (0,1] and diis_start outside
      !! {-1} union (0, infinity) -- each written as a positive test
      !! ("is it in range") rather than a `<`/`>` exclusion, so a NaN
      !! comparison, always false, falls into the refusal instead of
      !! being let through. Default (adiis=0, ndamp=-1, damp=0.3,
      !! diis_start=-1) and one in-range non-default combination are both
      !! accepted.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: nan

      nan = ieee_value(nan, ieee_quiet_nan)
      rc = trc_create(h)

      ! adiis must be 0 or 1.
      rc = trc_set_scf_accel(h, 2_c_int, -1_c_int, 0.3_c_double, -1.0_c_double)
      call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted adiis = 2")

      ! ndamp must be >= 0, or -1 for the by-spin default.
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -2_c_int, 0.3_c_double, -1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted ndamp = -2")
      end if

      ! damp must be in (0, 1]: 0, above 1, and NaN are all refused.
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, 0.0_c_double, -1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted damp = 0")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, 1.5_c_double, -1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted damp = 1.5")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, nan, -1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted damp = NaN")
      end if

      ! diis_start must be > 0, or exactly -1 for the by-spin default: 0, a
      ! negative value other than -1, and NaN are all refused.
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, 0.3_c_double, 0.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted diis_start = 0")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, 0.3_c_double, -2.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted diis_start = -2")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, 0.3_c_double, nan)
         call check(error, rc == TRC_ERR_BADARG, "trc_set_scf_accel accepted diis_start = NaN")
      end if

      ! The defaults are accepted.
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 0_c_int, -1_c_int, 0.3_c_double, -1.0_c_double)
         call check(error, rc == TRC_OK, "trc_set_scf_accel refused the defaults")
      end if
      ! And a valid non-default combination is accepted.
      if (.not. allocated(error)) then
         rc = trc_set_scf_accel(h, 1_c_int, 2_c_int, 0.5_c_double, 0.1_c_double)
         call check(error, rc == TRC_OK, "trc_set_scf_accel refused a valid non-default combination")
      end if

      rc = trc_destroy(h)
   end subroutine test_scf_accel_refused

end module check_fukui_suite

program check_fukui
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_fukui_suite, only: collect_fukui_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   testsuites = [new_testsuite("fukui", collect_fukui_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_fukui
