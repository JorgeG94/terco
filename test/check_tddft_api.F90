module check_tddft_api_suite
   !!
   !! The Fortran entry of TDDFT/TDHF (trc_tddft_run: trc_tddft_run_solve and
   !! trc_tddft_run_gradient), which builds the pair list, the engine and the
   !! XC grid itself, against the pyscf references of the solver and gradient
   !! checks (same water geometry, basis, SCF settings and grid level 3):
   !!  1. TDA and RPA excitation energies, oscillator strengths and ||Y||^2 of
   !!     water PBE/6-31G* against tddft_water_6-31gs_pbe.ref and
   !!     rpa_water_6-31gs_pbe.ref;
   !!  2. (arguments `grad_tda` and `grad_rpa`, one run each) the TDA and RPA gradient of the first and the last of
   !!     three roots of water PBE/6-31G against tddft_grad_water_6-31g_pbe.ref and tddft_rpa_grad_water_6-31g_pbe.ref
   !!     (the root number reaches the gradient through `opts%target_root`: the two ends would show an off-by-one,
   !!     the middle root adds a gradient call for nothing; the gradients of every root against pyscf are
   !!     check_tddft_grad_* and check_tddft_rpa_grad_*).
   !!     Those are Pulay-only (pyscf has no grid response), so this runs with
   !!     `grid_response = .false.`;
   !!  3. (argument `fd`) `grid_response` is ON by default for Kohn-Sham, differs from the
   !!     Pulay-only gradient by the grid's own derivative (measured below), and
   !!     is the exact derivative: it matches a Richardson finite difference of
   !!     E_SCF + omega over one atom's coordinates, the grid and SCF rebuilt at
   !!     every displacement, while the Pulay-only one misses it;
   !!  4. what is refused before any work: triplets, meta-GGA, range-separated,
   !!     VV10 and double-hybrid functionals, a spherical basis, an unknown
   !!     functional, a root outside 1..nroots, --df and a continuum solvent.
   !! Run from test/. A test-drive suite.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger, info_level
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_eri, only: trc_eri_t
   use trc_tddft, only: trc_tddft_op_t, trc_tddft_reference_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft_run, only: trc_tddft_run_options_t, trc_tddft_run_result_t, trc_tddft_run_solve, &
                            trc_tddft_run_gradient, trc_tddft_run_refusal, trc_tddft_run_gap_note
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector, ref_read_matrix3
   implicit none
   private
   public :: collect_check_tddft_api, mode_arg

   character(len=16) :: mode_arg = ""
   !! "grad_tda" and "grad_rpa" run one gradient check each, "fd" only the grid response finite difference (it
   !! rebuilds four SCFs); empty, the energies and the refusals. Separate runs of one program keep each around ten
   !! seconds (or, for "fd", LONG).

   real(dp), parameter :: TOL_OMEGA = 1.0e-6_dp
   real(dp), parameter :: TOL_OSC = 1.0e-5_dp
   real(dp), parameter :: TOL_YY = 1.0e-6_dp
   !! The solver checks' gates against pyscf (check_tddft, check_tddft_rpa).
   real(dp), parameter :: SCREEN_CONV_TOL = 1.0e-7_dp
   !! The Davidson residual tolerance of the density-screen checks.
   real(dp), parameter :: SCREEN_ERI_THRESH = 1.0e-8_dp
   !! The integral threshold of the matvec work check: loose enough that a small transition-density block drops quartets.
   real(dp), parameter :: GATE_GRAD = 5.0e-9_dp
   !! The gradient checks' gate against pyscf (check_tddft_grad).
   real(dp), parameter :: FD_STEP = 2.5e-3_dp
   real(dp), parameter :: FD_GATE = 2.0e-7_dp
   !! Kohn-Sham finite difference with the grid response (check_tddft_grad's FD_ABSTOL_DFT).

   type :: fixture_t
      !! One SCF of water, with what the TDDFT entry needs.
      integer :: nocc = 0, natm = 0
      integer, allocatable :: zint(:)
      type(trc_basis_t) :: bas
      type(trc_linalg_t) :: la
      type(trc_scf_options_t) :: sopts
      type(trc_scf_result_t) :: scf
      logical :: ready = .false.
   end type fixture_t

contains

   subroutine collect_check_tddft_api(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (trim(mode_arg) == "fd") then
         testsuite = [new_unittest("grid response is the default and is the exact derivative", test_grid_response)]
      else if (trim(mode_arg) == "missed_root") then
         testsuite = [ &
                     new_unittest("missed-root budget: TDA of formaldehyde/6-31+G HF", test_missed_root_tda), &
                     new_unittest("missed-root budget: RPA of formaldehyde/6-31+G HF", test_missed_root_rpa) &
                     ]
      else if (trim(mode_arg) == "density_screen") then
         testsuite = [ &
                     new_unittest("density screen: TDA of trans-glyoxal/3-21G HF, on and off", test_density_screen_tda), &
                     new_unittest("density screen: RPA of trans-glyoxal/3-21G HF, on and off", test_density_screen_rpa), &
                     new_unittest("density screen: the TDA matvec drops work with the screen on", test_screen_work_tda), &
                     new_unittest("density screen: the RPA matvec changes at threshold level with the screen on", &
                                  test_screen_work_rpa) &
                     ]
      else if (trim(mode_arg) == "grad_tda") then
         testsuite = [new_unittest("TDA gradient of the first and last root vs pyscf (Pulay only)", test_gradient_tda)]
      else if (trim(mode_arg) == "grad_rpa") then
         testsuite = [new_unittest("RPA gradient of the first and last root vs pyscf (Pulay only)", test_gradient_rpa)]
      else
         testsuite = [ &
                     new_unittest("TDA energies vs pyscf (water PBE/6-31G*)", test_energies_tda), &
                     new_unittest("RPA energies and ||Y||^2 vs pyscf (water PBE/6-31G*)", test_energies_rpa), &
                     new_unittest("refusals before any work", test_refusals) &
                     ]
      end if
   end subroutine collect_check_tddft_api

   subroutine fixture_build(fx, bas_json, functional, error, xyz)
      !! water (or the molecule of `xyz`), the basis and one converged SCF at the settings of the pyscf references.
      type(fixture_t), intent(out) :: fx
      character(len=*), intent(in) :: bas_json, functional
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in), optional :: xyz
      type(error_t) :: err
      real(dp), allocatable :: at_r(:, :)
      integer :: dev

      dev = trc_bind_device(0)
      if (present(xyz)) then
         call read_xyz(xyz, fx%natm, fx%zint, at_r)
      else
         call read_xyz("water.xyz", fx%natm, fx%zint, at_r)
      end if
      call trc_basis_from_json(bas_json, fx%natm, fx%zint, at_r, fx%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%la%init(fx%bas%nao)
      fx%ready = .true.
      fx%nocc = nint(sum(real(fx%zint, dp)))/2
      fx%sopts%conv_energy = 1.0e-12_dp
      fx%sopts%conv_density = 1.0e-10_dp
      fx%sopts%conv_diis = 1.0e-11_dp
      fx%sopts%max_iteration = 200
      fx%sopts%grid_level = 3
      fx%sopts%functional = functional
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, fx%sopts, fx%scf)
      call check(error, fx%scf%converged, "SCF did not converge")
   end subroutine fixture_build

   subroutine fixture_release(fx)
      type(fixture_t), intent(inout) :: fx
      if (.not. fx%ready) return
      call fx%la%release()
      call fx%bas%release()
      fx%ready = .false.
   end subroutine fixture_release

   subroutine compare_energies(error, ref_name, tda, with_y)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: ref_name
      logical, intent(in) :: tda, with_y
      type(fixture_t) :: fx
      type(trc_tddft_run_options_t) :: opts
      type(trc_tddft_run_result_t) :: res
      real(dp), allocatable :: om_ref(:), osc_ref(:), yy_ref(:)
      real(dp) :: dom, dosc, dyy, e_ref
      integer :: u, nr

      nr = 4
      call fixture_build(fx, "basis_sets/6-31gs.json", "pbe", error)
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      opts%nroots = nr
      opts%tda = tda
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call check(error, res%converged, "solve failed: "//res%err%get_message())
      if (.not. allocated(error)) then
         call ref_require(ref_name, u)
         allocate (om_ref(nr), osc_ref(nr), yy_ref(nr))
         call ref_read_real(u, "E_SCF", e_ref)
         call ref_read_vector(u, "OMEGA", nr, om_ref)
         call ref_read_vector(u, "OSC", nr, osc_ref)
         if (with_y) call ref_read_vector(u, "YY", nr, yy_ref)
         call ref_close(u)
         dom = maxval(abs(res%omega - om_ref))
         dosc = maxval(abs(res%osc - osc_ref))
         call global_logger%info("tddft api "//merge("TDA", "RPA", tda)//": max|d omega| "//to_string(dom, "es10.2")// &
                                 ", max|d f| "//to_string(dosc, "es10.2")//", "//to_string(res%iterations)//" iterations")
         call check(error, abs(res%e_scf - e_ref) < 1.0e-9_dp, "E_SCF vs pyscf")
         if (.not. allocated(error)) call check(error, dom < TOL_OMEGA, "omega vs pyscf, "//to_string(dom, "es10.2"))
         if (.not. allocated(error)) call check(error, dosc < TOL_OSC, "oscillator strength vs pyscf, "// &
                                                to_string(dosc, "es10.2"))
         if (.not. allocated(error) .and. with_y) then
            dyy = maxval(abs(res%ynorm - yy_ref))
            call check(error, dyy < TOL_YY, "||Y||^2 vs pyscf, "//to_string(dyy, "es10.2"))
         end if
         if (.not. allocated(error) .and. .not. with_y) then
            call check(error, maxval(res%ynorm) == 0.0_dp, "the TDA has Y = 0")
         end if
         if (.not. allocated(error)) then
            call check(error, all(res%hole >= 1) .and. all(res%particle >= 1) .and. all(res%weight > 0.0_dp), &
                       "every root has a dominant orbital pair")
         end if
      end if
      call fixture_release(fx)
   end subroutine compare_energies

   subroutine test_energies_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_energies(error, "tddft_water_6-31gs_pbe.ref", .true., .false.)
   end subroutine test_energies_tda

   subroutine test_energies_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_energies(error, "rpa_water_6-31gs_pbe.ref", .false., .true.)
   end subroutine test_energies_rpa

   subroutine budget_solve(fx, tda, budget, res)
      !! One solve of the fixture's reference, three roots, with `missed_root_iterations = budget` (a negative
      !! budget leaves the option at its default).
      type(fixture_t), intent(inout) :: fx
      logical, intent(in) :: tda
      integer, intent(in) :: budget
      type(trc_tddft_run_result_t), intent(out) :: res
      type(trc_tddft_run_options_t) :: opts
      opts%nroots = 3
      opts%tda = tda
      opts%conv_tol = 1.0e-7_dp
      if (budget >= 0) opts%missed_root_iterations = budget
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
   end subroutine budget_solve

   subroutine missed_root_budget(error, tda)
      !! The budget of the candidate Ritz pairs after the probe, on formaldehyde/6-31+G Hartree-Fock, where the
      !! third root (0.348172 Eh) sits in a symmetry sector only that chase reaches (trc_tddft's MISSED ROOTS):
      !! the default of 24 is the explicit 24 (same iterations and matvecs, same roots); a budget of 0 uses fewer
      !! matvecs, here loses that root, and says so: `window_open` is set (and the default's is not).
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: tda
      type(fixture_t) :: fx
      type(trc_tddft_run_result_t) :: dflt, b24, b0
      character(len=3) :: what

      what = merge("TDA", "RPA", tda)
      call fixture_build(fx, "../basis_sets/6-31+g.json", "", error, xyz="formaldehyde.xyz")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call budget_solve(fx, tda, -1, dflt)
      call budget_solve(fx, tda, 24, b24)
      call budget_solve(fx, tda, 0, b0)
      call check(error, dflt%converged .and. b24%converged .and. b0%converged, what//": a solve did not converge")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call global_logger%info("tddft api "//what//" missed roots: default "//to_string(dflt%td%iterations)//" iterations "// &
                              to_string(dflt%td%nmatvec)//" matvecs open "//to_string(dflt%td%nwindow_open)//"; budget 0 "// &
                              to_string(b0%td%iterations)//" iterations "//to_string(b0%td%nmatvec)//" matvecs open "// &
                              to_string(b0%td%nwindow_open)//"; roots "//to_string(dflt%omega(3), "f10.6")//" / "// &
                              to_string(b0%omega(3), "f10.6"))
      ! Same solver path twice: the counts agree, the roots to the device's run-to-run noise (atomics, ~1e-12).
      call check(error, dflt%td%iterations == b24%td%iterations .and. dflt%td%nmatvec == b24%td%nmatvec .and. &
                 maxval(abs(dflt%omega - b24%omega)) < 1.0e-9_dp, what//": the default is not a budget of 24")
      if (.not. allocated(error)) call check(error, b0%td%nmatvec < dflt%td%nmatvec, what//": budget 0 used "// &
                                             to_string(b0%td%nmatvec)//" matvecs, the default "//to_string(dflt%td%nmatvec))
      if (.not. allocated(error)) call check(error, .not. dflt%td%window_open, what//": the default left candidates open")
      if (.not. allocated(error)) call check(error, maxval(abs(b0%omega - dflt%omega)) > 1.0e-3_dp, what// &
                                             ": the chase changed no root here, so this case does not exercise the budget")
      if (.not. allocated(error)) call check(error, b0%td%window_open .and. b0%td%nwindow_open > 0, what// &
                                             ": budget 0 lost a root and did not report an open window")
      call fixture_release(fx)
   end subroutine missed_root_budget

   subroutine test_missed_root_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call missed_root_budget(error, .true.)
   end subroutine test_missed_root_tda

   subroutine test_missed_root_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call missed_root_budget(error, .false.)
   end subroutine test_missed_root_rpa


   function logged_text(file) result(text)
      !! A log file read back whole.
      character(len=*), intent(in) :: file
      character(len=:), allocatable :: text
      character(len=512) :: line
      integer :: u, ios
      text = ""
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      do
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) exit
         text = text//trim(line)//new_line("a")
      end do
      close (u)
   end function logged_text

   integer function count_text(text, word) result(n)
      !! Occurrences of `word` in `text`.
      character(len=*), intent(in) :: text, word
      integer :: pos, hit
      n = 0
      pos = 1
      do
         hit = index(text(pos:), word)
         if (hit == 0) exit
         n = n + 1
         pos = pos + hit + len(word) - 1
      end do
   end function count_text

   subroutine screen_solve(fx, tda, screen, res, text)
      !! Three roots of the fixture's reference, with the matvec's density screen `screen`; what the solve logged.
      type(fixture_t), intent(inout) :: fx
      logical, intent(in) :: tda, screen
      type(trc_tddft_run_result_t), intent(out) :: res
      character(len=:), allocatable, intent(out) :: text
      type(trc_tddft_run_options_t) :: opts
      character(len=*), parameter :: LOG_FILE = "check_tddft_api_density_screen.log"
      integer :: u
      opts%nroots = 3
      opts%tda = tda
      opts%conv_tol = SCREEN_CONV_TOL
      opts%density_screen = screen
      open (newunit=u, file=LOG_FILE, status="replace")
      close (u)
      call global_logger%configure_file_output(LOG_FILE, level=info_level)
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call global_logger%close_log_file()
      text = logged_text(LOG_FILE)
   end subroutine screen_solve

   subroutine density_screen_roots(error, tda)
      !! The matvec's density screen on (the default) and off, on trans-glyoxal in 3-21G (Hartree-Fock: the exchange
      !! is the whole two-electron part, and the molecule is extended enough for the screen to drop quartets).
      !! The roots agree to twice the residual tolerance: a Ritz value with residual norm r lies within r of an eigenvalue
      !! of its own operator (r = conv_tol at convergence), and the two operators differ by quartets below the integral
      !! threshold, which are orders of magnitude under conv_tol (the observed difference is logged). The unscreened run
      !! logs ONE warning that says so and the default none, which also shows the option reaches the solver.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: tda
      type(fixture_t) :: fx
      type(trc_tddft_run_result_t) :: on, off, dflt
      character(len=:), allocatable :: log_on, log_off, log_dflt
      character(len=*), parameter :: WORDS = "density screen of the excited-state matvec is OFF"
      type(trc_tddft_run_options_t) :: defaults
      real(dp) :: dom
      character(len=3) :: what

      what = merge("TDA", "RPA", tda)
      call check(error, defaults%density_screen, "the default of density_screen is not true")
      if (allocated(error)) return
      call fixture_build(fx, "../basis_sets/3-21g.json", "", error, xyz="glyoxal_trans.xyz")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call screen_solve(fx, tda, .true., on, log_on)
      call screen_solve(fx, tda, .false., off, log_off)
      call check(error, on%converged .and. off%converged, what//": a solve did not converge")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      dom = maxval(abs(on%omega - off%omega))
      call global_logger%info("tddft api "//what//" density screen: max|d omega| on - off "//to_string(dom, "es10.3")// &
                              ", "//to_string(on%td%nmatvec)//" and "//to_string(off%td%nmatvec)//" matvecs")
      call check(error, dom <= 2.0_dp*SCREEN_CONV_TOL, what//": screened and unscreened roots differ by "// &
                 to_string(dom, "es10.3"))
      ! What the operator DID, not what the options said: the matvec passes that ran screened.
      if (.not. allocated(error)) call check(error, on%td%nscreened_apply > 0 .and. on%td%nscreened_apply <= on%td%iterations, &
                                             what//": the screened run's matvec passes did not run screened: "// &
                                             to_string(on%td%nscreened_apply))
      if (.not. allocated(error)) call check(error, off%td%nscreened_apply == 0, what// &
                                             ": the unscreened run's matvec ran screened "//to_string(off%td%nscreened_apply)// &
                                             " times")
      if (.not. allocated(error)) call check(error, count_text(log_on, WORDS) == 0, what//": the screened run warned")
      if (.not. allocated(error)) call check(error, count_text(log_off, WORDS) == 1, what//": the unscreened run logged "// &
                                             to_string(count_text(log_off, WORDS))//" warnings, not one")
      if (.not. allocated(error)) then
         ! The default is the screened run.
         call screen_default(fx, tda, dflt, log_dflt)
         call check(error, count_text(log_dflt, WORDS) == 0 .and. maxval(abs(dflt%omega - on%omega)) < SCREEN_CONV_TOL, &
                    what//": the default run is not the screened one")
         if (.not. allocated(error)) call check(error, dflt%td%nscreened_apply == on%td%nscreened_apply, &
                                                what//": the default run's matvec is not the screened one")
      end if
      call fixture_release(fx)
   end subroutine density_screen_roots

   subroutine screen_default(fx, tda, res, text)
      !! The same solve with the options at their defaults (density_screen not touched).
      type(fixture_t), intent(inout) :: fx
      logical, intent(in) :: tda
      type(trc_tddft_run_result_t), intent(out) :: res
      character(len=:), allocatable, intent(out) :: text
      type(trc_tddft_run_options_t) :: opts
      character(len=*), parameter :: LOG_FILE = "check_tddft_api_density_screen.log"
      integer :: u
      opts%nroots = 3
      opts%tda = tda
      opts%conv_tol = SCREEN_CONV_TOL
      open (newunit=u, file=LOG_FILE, status="replace")
      close (u)
      call global_logger%configure_file_output(LOG_FILE, level=info_level)
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call global_logger%close_log_file()
      text = logged_text(LOG_FILE)
   end subroutine screen_default

   subroutine test_density_screen_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call density_screen_roots(error, .true.)
   end subroutine test_density_screen_tda

   subroutine test_density_screen_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call density_screen_roots(error, .false.)
   end subroutine test_density_screen_rpa


   subroutine screen_work(error, tda)
      !! The matvec itself, through `trc_tddft_op_t`, with `density_screen` on and off, on trial vectors the screen can
      !! drop work for: X = the lowest occupied orbital (a core orbital, concentrated on the oxygens) into the lowest
      !! two virtuals, so the transition densities are small on most shell pairs of trans-glyoxal. Where a deleted
      !! solver copy (`op%density_screen = opts%density_screen`) or a hard-coded flag in either apply leaves the roots
      !! alone, this reads the work done:
      !!  - TDA: `count_survivors` (the operator's test hook) reads back the quartets that survived both screens,
      !!    which must fall with the screen on, and the result moves by a nonzero amount no larger than the dropped
      !!    quartets allow, 2 x 8 x ncart^2 (|j| + |k|) t each (j = 2, k = 2 c_x for HF; the derivation is in
      !!    check_fock_split's run_screen_threshold);
      !!  - RPA: the enumerated kernel has no survivor counter, so the observable is the result, which must differ
      !!    from the unscreened one (the screen engaged). No bound is asserted: without a survivor count there is no
      !!    number of dropped quartets to derive one from;
      !!  - both: `napply_screened` counts the screened passes.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: tda
      type(fixture_t) :: fx
      type(trc_eri_t) :: eri
      type(trc_tddft_reference_t) :: ref
      type(trc_tddft_op_t) :: op
      type(error_t) :: err
      real(dp), allocatable :: x(:, :, :), ax_on(:, :, :, :), ax_off(:, :, :, :)
      real(dp) :: diff, bound
      integer(int64) :: nk_on, nk_off
      integer :: nvir, k, ncart, lmax
      character(len=3) :: what

      what = merge("TDA", "RPA", tda)
      call fixture_build(fx, "../basis_sets/3-21g.json", "", error, xyz="glyoxal_trans.xyz")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call eri%build(fx%bas, SCREEN_ERI_THRESH)
      call ref%set(fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1))
      call op%setup(fx%bas, fx%la, ref, err)
      call check(error, .not. err%has_error(), "operator setup: "//err%get_message())
      if (allocated(error)) then
         call eri%release()
         call fixture_release(fx)
         return
      end if
      nvir = op%nvir
      allocate (x(fx%nocc, nvir, 2), ax_on(fx%nocc, nvir, 2, 2), ax_off(fx%nocc, nvir, 2, 2))
      x = 0.0_dp
      do k = 1, 2
         x(1, k, k) = 1.0_dp
      end do
      lmax = maxval(fx%bas%sh_l)
      ncart = (lmax + 1)*(lmax + 2)/2
      op%count_survivors = .true.
      ax_on = 0.0_dp
      ax_off = 0.0_dp
      op%density_screen = .true.
      if (tda) then
         call op%apply(fx%bas, eri, fx%la, 2, x, ax_on(:, :, :, 1))
      else
         call op%apply_rpa(fx%bas, eri, fx%la, 2, x, ax_on)
      end if
      nk_on = eri%nkept
      call check(error, op%napply_screened == 1, what//": the screened apply was not counted as screened")
      op%density_screen = .false.
      if (tda) then
         call op%apply(fx%bas, eri, fx%la, 2, x, ax_off(:, :, :, 1))
      else
         call op%apply_rpa(fx%bas, eri, fx%la, 2, x, ax_off)
      end if
      nk_off = eri%nkept
      if (.not. allocated(error)) call check(error, op%napply_screened == 1, what// &
                                             ": the unscreened apply was counted as screened")
      diff = maxval(abs(ax_on - ax_off))
      call global_logger%info("tddft api "//what//" matvec work: survivors "//to_string(nk_on)//" screened, "// &
                              to_string(nk_off)//" unscreened; max|screened - unscreened| "//to_string(diff, "es10.3"))
      if (.not. allocated(error)) call check(error, diff > 0.0_dp, what//": the screen changed nothing")
      if (tda .and. .not. allocated(error)) then
         call check(error, nk_on < nk_off, "the screened TDA matvec dropped no quartet ("//to_string(nk_on)//" of "// &
                    to_string(nk_off)//")")
         ! j = 2, k = 2 c_x = 2 (Hartree-Fock): the folded kernels' bound, per dropped quartet.
         bound = 2.0_dp*8.0_dp*real(ncart*ncart, dp)*(2.0_dp + 2.0_dp)*SCREEN_ERI_THRESH*real(nk_off - nk_on, dp)
         if (.not. allocated(error)) call check(error, diff <= bound, "screened vs unscreened "//to_string(diff, "es10.3")// &
                                                " exceeds the bound "//to_string(bound, "es10.3"))
      end if
      call op%release()
      call eri%release()
      call fixture_release(fx)
   end subroutine screen_work

   subroutine test_screen_work_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call screen_work(error, .true.)
   end subroutine test_screen_work_tda

   subroutine test_screen_work_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call screen_work(error, .false.)
   end subroutine test_screen_work_rpa

   subroutine gradient_vs_ref(error, ref_name, tda)
      !! The first and the last root's Pulay-only gradient against the pyscf reference (3 roots).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: ref_name
      logical, intent(in) :: tda
      type(fixture_t) :: fx
      type(trc_tddft_run_options_t) :: opts
      type(trc_tddft_run_result_t) :: res
      real(dp), allocatable :: gref(:, :, :)
      real(dp) :: worst, d, e_ref
      integer :: u, k, nr

      nr = 3
      call fixture_build(fx, "basis_sets/6-31g.json", "pbe", error)
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      opts%nroots = nr
      opts%tda = tda
      opts%conv_tol = 1.0e-10_dp
      opts%cpks_tol = 1.0e-11_dp
      opts%grid_response = .false.
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call check(error, res%converged, "solve failed: "//res%err%get_message())
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call ref_require(ref_name, u)
      allocate (gref(3, fx%natm, nr))
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_matrix3(u, "GRAD", 3, fx%natm, nr, gref)
      call ref_close(u)
      worst = 0.0_dp
      do k = 1, nr, nr - 1
         opts%target_root = k
         call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
         call check(error, .not. res%err%has_error(), "gradient of root "//to_string(k)//": "//res%err%get_message())
         if (allocated(error)) exit
         call check(error, res%target_root == k, "the result names its root")
         if (allocated(error)) exit
         d = maxval(abs(res%grad - gref(:, :, k)))
         worst = max(worst, d)
         call global_logger%info("tddft api "//merge("TDA", "RPA", tda)//" gradient, root "//to_string(k)// &
                                 ": max|d grad| vs pyscf "//to_string(d, "es10.2")//" (gate "// &
                                 to_string(GATE_GRAD, "es8.1")//")")
      end do
      if (.not. allocated(error)) then
         call check(error, worst < GATE_GRAD, "gradient vs pyscf, worst "//to_string(worst, "es10.2"))
      end if
      if (.not. allocated(error)) then
         call check(error, abs(fx%scf%energy - e_ref) < 1.0e-9_dp, "E_SCF vs pyscf")
      end if
      call fixture_release(fx)
   end subroutine gradient_vs_ref

   subroutine test_gradient_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_vs_ref(error, "tddft_grad_water_6-31g_pbe.ref", .true.)
   end subroutine test_gradient_tda

   subroutine test_gradient_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call gradient_vs_ref(error, "tddft_rpa_grad_water_6-31g_pbe.ref", .false.)
   end subroutine test_gradient_rpa

   subroutine excited_energy(at_r, root, energy, error)
      !! E_SCF + omega of `root` (TDA, PBE/6-31G) with the nuclei at `at_r`, all rebuilt.
      real(dp), intent(in) :: at_r(:, :)
      integer, intent(in) :: root
      real(dp), intent(out) :: energy
      type(error_type), allocatable, intent(out) :: error
      type(fixture_t) :: fx
      type(trc_tddft_run_options_t) :: opts
      type(trc_tddft_run_result_t) :: res
      type(error_t) :: err
      integer :: dev

      energy = 0.0_dp
      dev = trc_bind_device(0)
      allocate (fx%zint(size(at_r, 2)))
      fx%zint = [8, 1, 1]
      fx%natm = size(at_r, 2)
      call trc_basis_from_json("basis_sets/6-31g.json", fx%natm, fx%zint, at_r, fx%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%la%init(fx%bas%nao)
      fx%ready = .true.
      fx%nocc = 5
      fx%sopts%conv_energy = 1.0e-12_dp
      fx%sopts%conv_density = 1.0e-10_dp
      fx%sopts%conv_diis = 1.0e-11_dp
      fx%sopts%max_iteration = 200
      fx%sopts%grid_level = 3
      fx%sopts%functional = "pbe"
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, fx%sopts, fx%scf)
      call check(error, fx%scf%converged, "displaced SCF did not converge")
      if (.not. allocated(error)) then
         opts%nroots = root
         opts%conv_tol = 1.0e-10_dp
         call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
         call check(error, res%converged, "displaced solve: "//res%err%get_message())
         if (.not. allocated(error)) energy = fx%scf%energy + res%omega(root)
      end if
      call fixture_release(fx)
   end subroutine excited_energy

   subroutine test_grid_response(error)
      type(error_type), allocatable, intent(out) :: error
      type(fixture_t) :: fx
      type(trc_tddft_run_options_t) :: opts, pulay
      type(trc_tddft_run_result_t) :: res, res_pulay
      real(dp), allocatable :: rp(:, :)
      real(dp) :: ep, em, coarse, fd, diff, exact_miss, pulay_miss
      integer, parameter :: ROOT = 2, ATOM = 1, DIR = 3

      call fixture_build(fx, "basis_sets/6-31g.json", "pbe", error)
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      opts%nroots = ROOT
      opts%target_root = ROOT
      opts%conv_tol = 1.0e-10_dp
      opts%cpks_tol = 1.0e-11_dp
      call check(error, opts%grid_response, "grid_response is on by default")
      if (allocated(error)) then
         call fixture_release(fx)
         return
      end if
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call check(error, .not. res%err%has_error(), "gradient: "//res%err%get_message())
      if (.not. allocated(error)) then
         pulay = opts
         pulay%grid_response = .false.
         call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, pulay, res_pulay)
         call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, pulay, res_pulay)
         call check(error, .not. res_pulay%err%has_error(), "Pulay-only gradient: "//res_pulay%err%get_message())
      end if
      if (.not. allocated(error)) then
         diff = maxval(abs(res%grad - res_pulay%grad))
         call global_logger%info("tddft api: grid response moves the gradient by "//to_string(diff, "es10.2")// &
                                 " (max over atoms and components)")
         call check(error, diff > 1.0e-6_dp .and. diff < 1.0e-3_dp, "the grid response should move the water PBE "// &
                    "gradient by about 1e-5, got "//to_string(diff, "es10.2"))
      end if
      if (.not. allocated(error)) then
         ! Richardson (4 D(h/2) - D(h))/3 of the central difference of E_SCF + omega, one coordinate.
         rp = fx%bas%at_r
         rp(DIR, ATOM) = fx%bas%at_r(DIR, ATOM) + FD_STEP
         call excited_energy(rp, ROOT, ep, error)
         if (.not. allocated(error)) then
            rp(DIR, ATOM) = fx%bas%at_r(DIR, ATOM) - FD_STEP
            call excited_energy(rp, ROOT, em, error)
         end if
         coarse = (ep - em)/(2.0_dp*FD_STEP)
         if (.not. allocated(error)) then
            rp(DIR, ATOM) = fx%bas%at_r(DIR, ATOM) + 0.5_dp*FD_STEP
            call excited_energy(rp, ROOT, ep, error)
         end if
         if (.not. allocated(error)) then
            rp(DIR, ATOM) = fx%bas%at_r(DIR, ATOM) - 0.5_dp*FD_STEP
            call excited_energy(rp, ROOT, em, error)
         end if
         fd = (4.0_dp*(ep - em)/FD_STEP - coarse)/3.0_dp
      end if
      if (.not. allocated(error)) then
         exact_miss = abs(res%grad(DIR, ATOM) - fd)
         pulay_miss = abs(res_pulay%grad(DIR, ATOM) - fd)
         call global_logger%info("tddft api: d(E_SCF + omega_2)/dz(O) finite difference "//to_string(fd, "f14.10")// &
                                 "; grid response misses by "//to_string(exact_miss, "es10.2")// &
                                 ", Pulay only by "//to_string(pulay_miss, "es10.2"))
         call check(error, exact_miss < FD_GATE, "the default gradient vs finite difference, "// &
                    to_string(exact_miss, "es10.2"))
      end if
      if (.not. allocated(error)) then
         call check(error, pulay_miss > 10.0_dp*FD_GATE, "the Pulay-only gradient should miss the finite difference "// &
                    "by more than ten gates, so the test tells the two apart, got "//to_string(pulay_miss, "es10.2"))
      end if
      call fixture_release(fx)
   end subroutine test_grid_response

   subroutine expect_refusal(error, opts, sopts, word)
      !! The run is refused, and the message contains `word`.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_run_options_t), intent(in) :: opts
      type(trc_scf_options_t), intent(in) :: sopts
      character(len=*), intent(in) :: word
      character(len=256) :: msg
      msg = trc_tddft_run_refusal(sopts, opts)
      call check(error, len_trim(msg) > 0, word//" must be refused")
      if (allocated(error)) return
      call check(error, index(msg, word) > 0, "the message for '"//word//"' says: "//trim(msg))
   end subroutine expect_refusal

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_run_options_t) :: opts
      type(trc_scf_options_t) :: sopts
      type(trc_basis_t) :: bas
      type(error_t) :: err
      real(dp), allocatable :: at_r(:, :)
      integer, allocatable :: zint(:)
      integer :: natm
      character(len=256) :: msg

      sopts%functional = "pbe"
      call check(error, len_trim(trc_tddft_run_refusal(sopts, opts)) == 0, "a plain GGA is accepted")
      if (allocated(error)) return
      opts%missed_root_iterations = 0
      call check(error, len_trim(trc_tddft_run_refusal(sopts, opts)) == 0, "a missed-root budget of zero is accepted")
      if (allocated(error)) return
      opts%missed_root_iterations = -1
      call expect_refusal(error, opts, sopts, "missed_root_iterations")
      if (allocated(error)) return
      opts%missed_root_iterations = 24
      opts%singlet = .false.
      call expect_refusal(error, opts, sopts, "triplet")
      if (allocated(error)) return
      opts%singlet = .true.
      opts%target_root = opts%nroots + 1
      call expect_refusal(error, opts, sopts, "target root")
      if (allocated(error)) return
      opts%target_root = 1
      sopts%functional = "tpss"
      call expect_refusal(error, opts, sopts, "meta-GGA")
      if (allocated(error)) return
      sopts%functional = "wb97x"
      call expect_refusal(error, opts, sopts, "range-separated")
      if (allocated(error)) return
      sopts%functional = "wb97x-v"
      call expect_refusal(error, opts, sopts, "range-separated")
      if (allocated(error)) return
      sopts%functional = "b2plyp"
      call expect_refusal(error, opts, sopts, "double hybrid")
      if (allocated(error)) return
      sopts%functional = "nosuchfunc"
      call expect_refusal(error, opts, sopts, "nosuchfunc")
      if (allocated(error)) return
      sopts%functional = "pbe"
      sopts%df = .true.
      call expect_refusal(error, opts, sopts, "density-fitted")
      if (allocated(error)) return
      sopts%df = .false.
      sopts%pcm_model = "cpcm"
      call expect_refusal(error, opts, sopts, "continuum")
      if (allocated(error)) return
      sopts%pcm_model = ""
      ! A spherical basis: the check sees it on the basis.
      call read_xyz("water.xyz", natm, zint, at_r)
      call trc_basis_from_json("basis_sets/cc-pvtz.json", natm, zint, at_r, bas, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      msg = trc_tddft_run_refusal(sopts, opts, bas)
      call bas%release()
      call check(error, index(msg, "spherical") > 0, "a spherical basis must be refused, got: "//trim(msg))
      if (allocated(error)) return
      ! The gap helper the optimiser's state-following check uses.
      call check(error, trc_tddft_run_gap_note([0.30_dp, 0.31_dp, 0.50_dp], 1, 0.02_dp) == 2, "root 2 is within 0.02 of root 1")
      if (allocated(error)) return
      call check(error, trc_tddft_run_gap_note([0.30_dp, 0.31_dp, 0.50_dp], 3, 0.02_dp) == 0, "root 3 has no near neighbour")
   end subroutine test_refusals

end module check_tddft_api_suite

program check_tddft_api
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_tddft_api_suite, only: collect_check_tddft_api, mode_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, mode_arg)
   stat = 0
   call run_testsuite(collect_check_tddft_api, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft_api
