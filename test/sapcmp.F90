module sapcmp_suite
   !! Fitted SAP against the grid SAP: two independent constructions of the
   !! same screening potential, which is the only check that settles the
   !! normalisation without trusting either one.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_sap, only: trc_sap_build, trc_sap_grid
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sapcmp, xyz_arg, bjs_arg

   character(len=256) :: xyz_arg = 'water.xyz'
   !! The molecule; optionally overridden by the program from argv(1).
   character(len=256) :: bjs_arg = 'basis_sets/6-31g.json'
   !! The basis JSON; optionally overridden by the program from argv(2).

contains

   subroutine collect_sapcmp(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fitted SAP vs grid SAP normalisation", test_normalisation), &
                  new_unittest("far field is the screened point charge", test_far_field) &
                  ]
   end subroutine collect_sapcmp

   subroutine test_normalisation(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(error_t) :: err
      integer :: natm, dev, i, nocc
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), smat(:, :), tmat(:, :), vmat(:, :)
      real(dp), allocatable :: hgrid(:, :), hfit(:, :), sgrid(:, :)
      real(dp) :: mx, mxg, ratio, sign_dot
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: rg, rf

      dev = trc_bind_device(0)
      print '(a)', "sapcmp: "//trim(xyz_arg)//"  "//trim(bjs_arg)
      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      ! trc_sap_build/trc_sap_grid build directly at `bas%nao` (Cartesian;
      ! SAP's own guess machinery is not spherical-aware, see trc_scf_driver.F90's
      ! header), so this comparison forces Cartesian regardless of what basis
      ! a caller names on the command line.
      call trc_basis_from_json(trim(bjs_arg), natm, zint, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
      allocate (hgrid(bas%nao, bas%nao), hfit(bas%nao, bas%nao), sgrid(bas%nao, bas%nao))

      call trc_sap_grid(bas, tmat, vmat, hgrid, err)
      call check(error, .not. err%has_error(), "grid: "//err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if
      ! strip T + V_ne so only the screening potential is compared
      sgrid = hgrid - tmat - vmat

      call trc_sap_build(bas, pl, tmat, vmat, hfit, err, verbose=.true.)
      call check(error, .not. err%has_error(), "fitted: "//err%get_message())
      if (allocated(error)) then
         call pl%release(); call bas%release()
         return
      end if

      mxg = maxval(abs(sgrid)); mx = maxval(abs(sgrid - (hfit - tmat - vmat)))
      print '(a,es14.6)', "max |Vscr_grid|          ", mxg
      print '(a,es14.6)', "max |Vscr_fit|           ", maxval(abs(hfit - tmat - vmat))
      print '(a,es14.6)', "max |grid - fit|         ", mx
      print '(a,f12.6)',  "relative                 ", mx/mxg
      ratio = sum(sgrid*(hfit - tmat - vmat))/max(sum((hfit - tmat - vmat)**2), 1.0e-300_dp)
      print '(a,f12.6)',  "least-squares grid/fit   ", ratio
      sign_dot = sum(sgrid*(hfit - tmat - vmat))
      print '(a)', "diagonal, first six:"
      do i = 1, min(6, bas%nao)
         print '(i4,2f16.8)', i, sgrid(i, i), hfit(i, i) - tmat(i, i) - vmat(i, i)
      end do

      ! And does an SCF started from each actually behave the same? The
      ! matrices agreeing says nothing about that if the guess is used
      ! differently. Printed, not gated: the normalisation checks below
      ! are the actual test.
      nocc = sum(zint)/2
      opts%conv_energy = 1.0e-8_dp
      opts%max_iteration = 30
      call trc_scf_run(bas, nocc, nocc, opts, rg, hguess=hgrid)
      call trc_scf_run(bas, nocc, nocc, opts, rf, hguess=hfit)
      print '(a,f20.10,a,i0,a,l1)', "grid  E ", rg%energy, " (", rg%iterations, ")  conv ", rg%converged
      print '(a,f20.10,a,i0,a,l1)', "fit   E ", rf%energy, " (", rf%iterations, ")  conv ", rf%converged

      call pl%release(); call bas%release()
      deallocate (smat, tmat, vmat, hgrid, hfit, sgrid)

      !
      ! The two are not the same potential and must not be asserted equal. The
      ! grid one is built from terco's own Hartree-Fock free atom with a
      ! Dirac-Slater exchange potential; the fitted one is Lehtola's, from his
      ! atomic calculation and his fit. They agree to about 1%, which is the
      ! methodological difference between them.
      !
      ! What this test is for is the NORMALISATION, where the failures are not
      ! one percent but factors: 2*sqrt(pi) for the coefficient convention the
      ! builder applies, (2 alpha/pi)^(3/4) against (alpha/pi)^(3/2) for
      ! wavefunction against charge normalisation, and a sign for which way the
      ! electrons point. Every one of those lands far outside this window and
      ! none of them changes the converged energy, so a bound this loose is
      ! still the sharpest instrument available.
      !
      call check(error, mx/mxg <= 0.05_dp, "the two constructions disagree by more than 5%")
      if (allocated(error)) return
      call check(error, abs(ratio - 1.0_dp) <= 0.05_dp, "there is a scale factor between them")
      if (allocated(error)) return
      call check(error, sign_dot > 0.0_dp, "the fitted potential has the opposite sign")
   end subroutine test_normalisation

   subroutine test_far_field(error)
      !! FAR FIELD, WHICH IS THE ONLY ABSOLUTE CHECK HERE.
      !!
      !! Everything in test_normalisation compares two constructions
      !! against each other, so a factor common to both would survive it.
      !! This does not: one oxygen at the origin, and a tight s function
      !! 40 Bohr away carrying no atom of its own. At that distance the
      !! oxygen is a point charge of eight electrons, so the matrix
      !! element over the normalised far function is exactly
      !!
      !!     Vscr = Z / R = 8 / 40 = 0.2
      !!
      !! and nothing about the fit, the normalisation or the sign can be
      !! wrong and still land there.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: fb
      type(trc_pairlist_t) :: fpl
      type(error_t) :: err
      real(dp) :: f_e(1, 2), f_c(1, 2), f_r(3, 2), f_z(1), f_ar(3, 1)
      integer  :: f_l(2), f_np(2)
      real(dp), allocatable :: ft(:, :), fv(:, :), fs(:, :), fh(:, :)
      real(dp) :: want, got

      f_l = 0; f_np = 1
      f_e(1, 1) = 4.0_dp; f_c(1, 1) = 1.0_dp
      f_e(1, 2) = 4.0_dp; f_c(1, 2) = 1.0_dp
      f_r = 0.0_dp; f_r(1, 2) = 40.0_dp
      f_z(1) = 8.0_dp; f_ar = 0.0_dp
      call fb%build(2, f_l, f_np, f_e, f_c, f_r, 1, f_z, f_ar, 1)
      call fb%to_device()
      call fpl%build(fb, 1.0e-14_dp); call fpl%to_device()
      allocate (fs(fb%nao, fb%nao), ft(fb%nao, fb%nao), fv(fb%nao, fb%nao), fh(fb%nao, fb%nao))
      call trc_1e(fb, fpl, fs, ft, fv)
      call trc_sap_build(fb, fpl, ft, fv, fh, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) then
         call fpl%release(); call fb%release()
         return
      end if
      ! strip T and V_ne, and divide by the overlap of the far function with
      ! itself, which the builder has normalised to one but is asked for
      ! rather than assumed
      got = (fh(2, 2) - ft(2, 2) - fv(2, 2))/fs(2, 2)
      want = 8.0_dp/40.0_dp
      print '(a,f14.9,a,f14.9,a,es11.3)', "far field: got ", got, "  want ", want, &
         "  rel ", abs(got - want)/want
      call fpl%release(); call fb%release()
      call check(error, abs(got - want)/want <= 1.0e-4_dp, "the far field is not the charge it should be")
   end subroutine test_far_field

end module sapcmp_suite

program sapcmp
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use sapcmp_suite, only: collect_sapcmp, xyz_arg, bjs_arg
   implicit none
   integer :: stat
   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bjs_arg)
   stat = 0
   call run_testsuite(collect_sapcmp, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program sapcmp
