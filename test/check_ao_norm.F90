!
! GAMESS's per-component AO normalisation, at the C API boundary.
!
! terco (like libcint) normalises a Cartesian shell per ANGULAR MOMENTUM: a
! d shell's six components share one constant, so <XY|XY> is 1/3, not 1.
! GAMESS normalises every component INDIVIDUALLY to 1. trc_set_ao_norm
! switches trc_capi's own AO-indexed entries between the two; see
! src/core/trc_ao_norm.F90's header for the closed form and src/api/trc_capi.F90's
! "AO NORMALISATION" section for exactly which entries convert.
!
! FORCED CARTESIAN. cc-pVDZ's d shells are tagged "gto_spherical" in BSE, so
! since the spherical AO space was added this basis defaults to running
! spherical (trc_api.F90's `b%spherical`). This check's entire subject is
! the per-COMPONENT GAMESS convention, which is defined on CARTESIAN
! components only and refuses to combine with a spherical basis at all
! (trc_set_ao_norm's own refusal, trc_capi.F90) -- so every basis load here
! forces Cartesian (`force_cart=.true.`/`trc_set_force_cartesian`), which is
! what keeps this check testing the same thing it always has rather than
! either failing outright on the refusal or silently changing subject.
!
! The basis is built once, in setup_ao_norm, before run_testsuite: a basis
! that fails to build is an environment problem the original plain-program
! version stopped on outright, not a per-test assertion, so this keeps that
! behaviour rather than inventing a check for it. SIX CHECKS follow, as
! test-drive unittests (parallel=.false., so the shell/unit SCF energies
! each test stores for (b)'s comparison are safe to share through
! module-level state).
!
!   (e) THE FACTOR ITSELF, l = 0..4, against test/reference_ao_norm.dat --
!       sqrt(gamess_norm_sq(lx,ly,lz)) from scripts/cartorder.py, generated
!       once and recorded (the same reason check_binfock reads reference.dat
!       rather than calling libfint from Fortran: no Python at ctest time).
!       Indexed by the EXPONENT TRIPLE, not by within-shell offset, so the
!       one file checks whichever TRC_CART_ORDER this build was made with --
!       order and normalisation are independent on purpose (both headers
!       above say why), and this test would rather prove that than assume
!       it. Through trc_ao_norm's `ao_norm_gamess` directly -- the pure ratio,
!       not `ao_norm_vector`'s full boundary conversion, which ALSO folds
!       in terco's own libcint-inherited `common_fac_sp` scale for s and p
!       (see trc_ao_norm.F90's header; that is a real but SEPARATE effect
!       this check is not the place to re-derive). Nothing here exercises
!       the actual kernels, because ONE_E_LMAX in this generated kernel set
!       is 2 and the formula does not need an integral to be checked
!       against itself -- see (a) below for where the kernels DO get
!       exercised, up to the d this build's one-electron path supports.
!
!   (a) THE OVERLAP DIAGONAL, water/cc-pVDZ (real d functions), through
!       trc_context_basis + trc_compute_1e: exactly 1 for every s, p and
!       "pure power" d component (XX/YY/ZZ), under EITHER mode -- terco's
!       own raw <s|s>/<p_x|p_x> already equals GAMESS's (and pyscf's)
!       unit convention, see trc_ao_norm.F90's header for why. Only the
!       three "mixed" d components (XY/XZ/YZ) differ by mode: 1/3 under
!       TRC_AO_NORM_SHELL (terco shares one normalisation constant across
!       a shell's components), 1 under TRC_AO_NORM_UNIT (GAMESS normalises
!       every component individually) -- found by looking up each AO's own
!       (px,py,pz) via batch_cart_pow, the same way (e) does, so this check does
!       not depend on knowing TRC_CART_ORDER's slot assignment by hand
!       either. Every AO in the basis is walked and the worst mismatch
!       against its own expected value/tolerance is what (a)'s check
!       reports, rather than one check per AO.
!
!       AN EARLIER VERSION OF THIS CHECK asserted 1/(4*pi) and 3/(4*pi) for
!       SHELL-mode s and p. That number was the trc_basis_json bug this
!       basis (loaded through trc_basis_from_json, below) used to have --
!       folding `common_fac_sp` twice into an s or p shell's coefficients --
!       reverse-engineered into the test as if it were the convention,
!       without checking it against the OTHER basis constructor
!       (test/trc_test_basis.F90's `build_631g`) or against pyscf, both of
!       which give 1 directly. It stayed green because
!       trc_ao_norm.F90's own `ao_norm_vector` carried a SECOND
!       `common_fac_sp` divide that canceled it for (b)/(c)/(d)/(f), which
!       test self-consistency and cannot see an overall basis-wide scale
!       error; only this check's hard-coded numbers could, and did, the
!       moment the JSON reader's fold was fixed on its own.
!
!   (b) THE SCF ENERGY, SHELL and UNIT, must agree to ~1e-12 Eh: the whole
!       point of a change of AO NORMALISATION is that it is not a change of
!       physics, only of what a getter hands back.
!
!   (c) Tr(D S) = N_electron, in EACH mode's own D and S -- the physical
!       invariant a contravariant-covariant pair is built to preserve
!       (Tr(N^-1 D N^-1 . N S N) = Tr(D S) exactly, for any diagonal N), so
!       this is a self-consistency check on ao_norm_scale_cov/ao_norm_scale_contra
!       agreeing about which direction is which.
!
!   (d) ROUND TRIP: the density trc_density hands back in a mode is fed
!       straight into trc_set_guess in the SAME mode; the SCF must already
!       be converged (few iterations, same energy). A mismatched
!       conversion (covariant where it should be contravariant, or a sign
!       error in the exponent) would still pass (b) and (c) -- both are
!       symmetric in a way that a transposed mistake cannot break -- but
!       would feed trc_run_scf a density that is not actually the converged
!       one, and it would take real iterations to get back.
!
!   (f) THE FOCK BUILD (`check_fock`), the entry GAMESS's own driver calls
!       for J/K (run_Fockbuild's role in hpcchem_driver.src): ONE physical
!       density, D_shell from a converged SHELL-mode SCF, fed to trc_fock
!       in SHELL mode directly and, EXPLICITLY converted here (D_unit =
!       N^-1 D_shell N^-1, not from a second, independently-converged SCF
!       -- see the subroutine's own comment for why that distinction is
!       the whole point), to trc_fock in UNIT mode. Checks F_unit = N
!       F_shell N to ~1e-12 relative, and that Tr(D F) -- the two-electron
!       energy expression -- agrees between the two mode's own (D, F)
!       pairs, which holds by the SAME congruence identity (c) uses,
!       independent of any SCF.
!
module check_ao_norm_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char, c_loc
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_batch, only: batch_cart_pow
   use trc_ao_norm, only: ao_norm_gamess, ao_norm_vector, ao_norm_scale_cov, ao_norm_scale_contra
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_ao_norm, trc_get_ao_norm, trc_set_force_cartesian, &
                               trc_context_basis, trc_context_eri, &
                               trc_pairs_create, trc_pairs_destroy, trc_compute_1e, trc_fock, &
                               trc_set_method, trc_set_convergence, trc_set_guess, trc_run_scf, &
                               trc_energy, trc_density, trc_nao, trc_converged, trc_iterations, &
                               TRC_OK, TRC_GUESS_GIVEN, &
                               TRC_AO_NORM_SHELL, TRC_AO_NORM_UNIT
   implicit none
   private
   public :: collect_ao_norm_tests, setup_ao_norm, teardown_ao_norm

   type(trc_basis_t) :: b
   integer :: natm
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   real(c_double) :: e_shell, e_unit
      !! e_shell/e_unit are written by test_shell_mode/test_unit_mode and
      !! read by test_energy_agreement -- safe only because run_testsuite
      !! here always runs parallel=.false., in this module's own list
      !! order.

contains

   subroutine setup_ao_norm()
      !! Built once, as the original plain program did at the top level:
      !! a basis that fails to build is an environment problem, so this
      !! stops the program exactly as before rather than reporting through
      !! test-drive.
      type(error_t) :: err
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, b, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', "check_ao_norm: "//err%get_message()
         stop 1
      end if
   end subroutine setup_ao_norm

   subroutine teardown_ao_norm()
      call b%release()
   end subroutine teardown_ao_norm

   subroutine collect_ao_norm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("(e) ao_norm_gamess matches scripts/cartorder.py, l=0..4", test_factors), &
                  new_unittest("shell mode: (a) overlap diagonal, SCF, (c) Tr(D S), (d) round trip", &
                               test_shell_mode), &
                  new_unittest("unit mode: (a) overlap diagonal, SCF, (c) Tr(D S), (d) round trip", &
                               test_unit_mode), &
                  new_unittest("(b) SCF energy agrees between shell and unit", test_energy_agreement), &
                  new_unittest("(f) Fock build: F_unit = N F_shell N, Tr(D F) agrees", test_fock) &
                  ]
   end subroutine collect_ao_norm_tests

   subroutine test_factors(error)
      !! (e) The closed form, l = 0..4, every component, against
      !! test/reference_ao_norm.dat -- indexed by (lx,ly,lz), not by
      !! within-shell offset, so this build's own TRC_CART_ORDER (via
      !! `batch_cart_pow`) decides which offset is checked against which line, and
      !! the one file serves every build.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NREF = 35
      integer :: ref_p(3, NREF)
      real(dp) :: ref_n(NREF)
      integer :: u, i, l, ic, ncart, px, py, pz, j
      real(dp) :: got, worst
      character(len=256) :: line
      logical :: found, all_found

      open (newunit=u, file='reference_ao_norm.dat', status='old', action='read')
      read (u, '(a)') line  ! header comment
      do i = 1, NREF
         read (u, *) ref_p(1, i), ref_p(2, i), ref_p(3, i), ref_n(i)
      end do
      close (u)

      all_found = .true.
      worst = 0.0_dp
      do l = 0, 4
         ncart = (l + 1)*(l + 2)/2
         do ic = 0, ncart - 1
            call batch_cart_pow(l, ic, px, py, pz)
            found = .false.
            do j = 1, NREF
               if (ref_p(1, j) == px .and. ref_p(2, j) == py .and. ref_p(3, j) == pz) then
                  found = .true.; exit
               end if
            end do
            all_found = all_found .and. found
            if (.not. found) cycle
            got = ao_norm_gamess(px, py, pz)
            worst = max(worst, abs(got - ref_n(j)))
         end do
      end do
      call check(error, all_found, "(e) every l=0..4 cartesian component is in reference_ao_norm.dat")
      if (allocated(error)) return
      call check(error, worst < 1.0e-12_dp, &
                 "(e) ao_norm_gamess matches scripts/cartorder.py, worst diff "//to_string(worst, "es10.2"))
   end subroutine test_factors

   subroutine test_shell_mode(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mode(TRC_AO_NORM_SHELL, "shell", e_shell, error)
   end subroutine test_shell_mode

   subroutine test_unit_mode(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mode(TRC_AO_NORM_UNIT, "unit", e_unit, error)
   end subroutine test_unit_mode

   subroutine test_energy_agreement(error)
      !! (b) THE SCF ENERGY, SHELL and UNIT, must agree to ~1e-12 Eh: the
      !! whole point of a change of AO NORMALISATION is that it is not a
      !! change of physics, only of what a getter hands back.
      type(error_type), allocatable, intent(out) :: error
      call check(error, abs(real(e_shell, dp) - real(e_unit, dp)) < 1.0e-12_dp, &
                 "(b) SCF energy agrees between shell ("//to_string(real(e_shell, dp), "f20.12")// &
                 ") and unit ("//to_string(real(e_unit, dp), "f20.12")//")")
   end subroutine test_energy_agreement

   subroutine run_mode(mode, label, e_out, error)
      !! Runs (a), the SCF, (c) and (d) in one AO-normalisation mode; hands
      !! back the SCF energy for (b)'s comparison across modes.
      integer(c_int), intent(in) :: mode
      character(len=*), intent(in) :: label
      real(c_double), intent(out) :: e_out
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h, basis, pairs
      integer(c_int) :: rc, nao, got_mode, flag
      real(c_double), allocatable, target :: smat(:), tmat(:), vmat(:), dmat(:)
      real(dp) :: trds, worst_a, expect_diag, tol_diag, diag_val
      integer :: i, ish, ic, ncart, base, px, py, pz, nbad_a
      integer(c_int) :: niter2

      rc = trc_create(h); call check(error, rc == TRC_OK, label//": create")
      if (allocated(error)) return
      rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call check(error, rc == TRC_OK, label//": set_molecule")
      if (allocated(error)) return
      rc = trc_set_ao_norm(h, mode); call check(error, rc == TRC_OK, label//": set_ao_norm")
      if (allocated(error)) return
      rc = trc_get_ao_norm(h, got_mode)
      call check(error, rc == TRC_OK .and. got_mode == mode, label//": get_ao_norm reads back what was set")
      if (allocated(error)) return
      rc = trc_set_force_cartesian(h, 1_c_int)
      call check(error, rc == TRC_OK, label//": set_force_cartesian")
      if (allocated(error)) return
      rc = trc_set_basis_json(h, 'basis_sets/cc-pvdz.json'//c_null_char)
      call check(error, rc == TRC_OK, label//": set_basis_json")
      if (allocated(error)) return
      rc = trc_nao(h, nao); call check(error, rc == TRC_OK .and. nao == b%nao, label//": nao")
      if (allocated(error)) return

      ! (a) THE OVERLAP DIAGONAL, through trc_context_basis + trc_compute_1e.
      rc = trc_context_basis(h, basis); call check(error, rc == TRC_OK, label//": context_basis")
      if (allocated(error)) return
      rc = trc_pairs_create(basis, 1.0e-14_c_double, pairs)
      call check(error, rc == TRC_OK, label//": pairs_create")
      if (allocated(error)) return
      allocate (smat(nao*nao), tmat(nao*nao), vmat(nao*nao))
      rc = trc_compute_1e(basis, pairs, smat, tmat, vmat)
      call check(error, rc == TRC_OK, label//": compute_1e")
      if (allocated(error)) return

      ! Every AO's own expected diagonal and tolerance (s, p, pure-power d
      ! are 1 under either mode; a mixed d component is 1/3 under SHELL,
      ! 1 under UNIT; anything else -- not reached by cc-pVDZ -- is skipped,
      ! exactly as the original per-AO branch fell through with no check).
      nbad_a = 0; worst_a = 0.0_dp
      do ish = 1, b%nshell
         base = b%sh_ao(ish) - 1
         ncart = (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2
         do ic = 0, ncart - 1
            call batch_cart_pow(b%sh_l(ish), ic, px, py, pz)
            i = base + ic + 1
            if (mode == TRC_AO_NORM_UNIT) then
               expect_diag = 1.0_dp; tol_diag = 1.0e-12_dp
            else if (b%sh_l(ish) == 0) then
               ! s: shell mode's raw diagonal is 1, same as GAMESS and
               ! pyscf -- see trc_ao_norm.F90's header for why.
               expect_diag = 1.0_dp; tol_diag = 1.0e-10_dp
            else if (b%sh_l(ish) == 1) then
               ! p: same reason, same answer.
               expect_diag = 1.0_dp; tol_diag = 1.0e-10_dp
            else if (px == 2 .or. py == 2 .or. pz == 2) then
               ! a "pure power" d/f/g component (XX-type): GAMESS's and
               ! terco's conventions coincide there, diagonal is 1.
               expect_diag = 1.0_dp; tol_diag = 1.0e-10_dp
            else if (b%sh_l(ish) == 2) then
               ! a mixed d component (XY/XZ/YZ): shell normalisation
               ! leaves <chi|chi> = 1/gamess_norm_sq = 1/3.
               expect_diag = 1.0_dp/3.0_dp; tol_diag = 1.0e-10_dp
            else
               cycle
            end if
            diag_val = real(smat((i - 1)*nao + i), dp)
            worst_a = max(worst_a, abs(diag_val - expect_diag))
            if (abs(diag_val - expect_diag) >= tol_diag) nbad_a = nbad_a + 1
         end do
      end do
      call check(error, nbad_a == 0, &
                 label//": (a) overlap diagonal matches the GAMESS/shell convention, worst diff "// &
                 to_string(worst_a, "es10.2"))
      if (allocated(error)) return

      ! --- the SCF, SAD guess, default convergence ---------------------
      rc = trc_set_method(h, ''//c_null_char, 3_c_int); call check(error, rc == TRC_OK, label//": set_method")
      if (allocated(error)) return
      rc = trc_set_convergence(h, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      call check(error, rc == TRC_OK, label//": set_convergence")
      if (allocated(error)) return
      rc = trc_run_scf(h); call check(error, rc == TRC_OK, label//": run_scf")
      if (allocated(error)) return
      rc = trc_converged(h, flag); call check(error, flag == 1, label//": converged")
      if (allocated(error)) return
      rc = trc_energy(h, e_out)

      ! (c) Tr(D S) = N_electron, THIS mode's own D and S.
      allocate (dmat(nao*nao))
      rc = trc_density(h, dmat); call check(error, rc == TRC_OK, label//": density")
      if (allocated(error)) return
      trds = 0.0_dp
      do i = 1, nao*nao
         trds = trds + real(dmat(i), dp)*real(smat(i), dp)
      end do
      ! S above is from BEFORE the SCF touched anything, so it is the same
      ! S under this same mode -- trc_compute_1e depends only on the basis
      ! and its ao_norm_kind, neither of which the SCF changes.
      call check(error, abs(trds - real(sum(zint), dp)) < 1.0e-8_dp, &
                 label//": (c) Tr(D S) = N_electron")
      if (allocated(error)) return

      ! (d) ROUND TRIP: this mode's own converged density, fed back as the
      ! guess in this SAME mode, must already be the answer.
      rc = trc_set_guess(h, TRC_GUESS_GIVEN, c_loc(dmat(1)), 1_c_int)
      call check(error, rc == TRC_OK, label//": (d) set_guess given")
      if (allocated(error)) return
      rc = trc_run_scf(h); call check(error, rc == TRC_OK, label//": (d) second run_scf")
      if (allocated(error)) return
      rc = trc_iterations(h, niter2)
      call check(error, niter2 <= 3, &
                 label//": (d) round trip converges in a few iterations, got "//to_string(int(niter2)))
      if (allocated(error)) return
      block
         real(c_double) :: e2
         rc = trc_energy(h, e2)
         call check(error, abs(real(e2, dp) - real(e_out, dp)) < 1.0e-10_dp, label//": (d) round-trip energy matches")
      end block

      rc = trc_pairs_destroy(pairs)
      rc = trc_destroy(h)
      deallocate (smat, tmat, vmat, dmat)
   end subroutine run_mode

   subroutine test_fock(error)
      !! (f) THE FOCK BUILD, the entry GAMESS's own driver calls for J/K
      !! (run_Fockbuild's role in hpcchem_driver.src) -- see the module
      !! header's own paragraph on this check.
      !!
      !! ONE physical density: D_shell, from a converged SHELL-mode SCF, fed
      !! to trc_fock in SHELL mode directly, and D_unit = N^-1 D_shell N^-1,
      !! built HERE by explicit conversion -- NOT from a second, independently
      !! converged SCF in UNIT mode. Two independent SCF runs agree only to
      !! the SCF's own convergence threshold (~1e-8 to ~1e-10 here), which
      !! would cap this check's precision at that, not the ~1e-12 the
      !! transform itself is exact to; an explicit conversion of the SAME
      !! array sidesteps that entirely. h_unit therefore never runs an SCF at
      !! all -- it exists only to lend trc_context_basis/trc_context_eri
      !! handles for trc_fock's UNIT-mode path to run on.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h_shell, h_unit, basis_shell, basis_unit, eri_shell, eri_unit
      integer(c_int) :: rc, nao
      real(c_double), allocatable, target :: d_shell(:), f_shell(:), d_unit(:), f_unit(:)
      real(dp), allocatable :: nrm(:), tmp(:, :)
      real(dp) :: trdf_shell, trdf_unit, maxdiff, maxf

      rc = trc_create(h_shell); call check(error, rc == TRC_OK, "fock: shell create")
      if (allocated(error)) return
      rc = trc_set_molecule(h_shell, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                            0_c_int, 1_c_int)
      rc = trc_set_ao_norm(h_shell, TRC_AO_NORM_SHELL)
      rc = trc_set_force_cartesian(h_shell, 1_c_int)
      rc = trc_set_basis_json(h_shell, 'basis_sets/cc-pvdz.json'//c_null_char)
      call check(error, rc == TRC_OK, "fock: shell set_basis_json")
      if (allocated(error)) return
      rc = trc_nao(h_shell, nao)
      rc = trc_set_method(h_shell, ''//c_null_char, 3_c_int)
      rc = trc_set_convergence(h_shell, 1.0e-11_c_double, 1.0e-8_c_double, 100_c_int, 10_c_int)
      rc = trc_run_scf(h_shell); call check(error, rc == TRC_OK, "fock: shell run_scf")
      if (allocated(error)) return

      allocate (d_shell(nao*nao), f_shell(nao*nao))
      rc = trc_density(h_shell, d_shell); call check(error, rc == TRC_OK, "fock: shell density")
      if (allocated(error)) return
      rc = trc_context_basis(h_shell, basis_shell); call check(error, rc == TRC_OK, "fock: shell context_basis")
      if (allocated(error)) return
      rc = trc_context_eri(h_shell, eri_shell); call check(error, rc == TRC_OK, "fock: shell context_eri")
      if (allocated(error)) return
      rc = trc_fock(eri_shell, basis_shell, d_shell, f_shell, 1.0_c_double, 1.0_c_double, 1_c_int, 0.0_c_double)
      call check(error, rc == TRC_OK, "fock: shell fock")
      if (allocated(error)) return

      ! D_unit = N^-1 D_shell N^-1: CONTRAVARIANT, the direction trc_density
      ! itself uses (see that entry's own comment) -- a density is built by
      ! CONTRACTING two AO indices, not by pairing two AO functions the way
      ! an integral does, so this is `ao_norm_scale_contra`, not `ao_norm_scale_cov`.
      allocate (nrm(nao), tmp(nao, nao))
      call ao_norm_vector(b%nshell, b%sh_l, b%sh_ao, nao, TRC_AO_NORM_UNIT, nrm)
      tmp = reshape(real(d_shell, dp), [nao, nao])
      call ao_norm_scale_contra(nrm, tmp)
      allocate (d_unit(nao*nao), f_unit(nao*nao))

      rc = trc_create(h_unit); call check(error, rc == TRC_OK, "fock: unit create")
      if (allocated(error)) return
      rc = trc_set_molecule(h_unit, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                            0_c_int, 1_c_int)
      rc = trc_set_ao_norm(h_unit, TRC_AO_NORM_UNIT)
      rc = trc_set_force_cartesian(h_unit, 1_c_int)
      rc = trc_set_basis_json(h_unit, 'basis_sets/cc-pvdz.json'//c_null_char)
      call check(error, rc == TRC_OK, "fock: unit set_basis_json")
      if (allocated(error)) return
      rc = trc_context_basis(h_unit, basis_unit); call check(error, rc == TRC_OK, "fock: unit context_basis")
      if (allocated(error)) return
      rc = trc_context_eri(h_unit, eri_unit); call check(error, rc == TRC_OK, "fock: unit context_eri")
      if (allocated(error)) return
      d_unit = real(reshape(tmp, [nao*nao]), c_double)
      rc = trc_fock(eri_unit, basis_unit, d_unit, f_unit, 1.0_c_double, 1.0_c_double, 1_c_int, 0.0_c_double)
      call check(error, rc == TRC_OK, "fock: unit fock")
      if (allocated(error)) return

      ! F_unit must equal N F_shell N.
      tmp = reshape(real(f_shell, dp), [nao, nao])
      call ao_norm_scale_cov(nrm, tmp)
      maxf = maxval(abs(tmp))
      maxdiff = maxval(abs(tmp - reshape(real(f_unit, dp), [nao, nao])))
      call check(error, maxdiff < 1.0e-12_dp*maxf, &
                 "(f) Fock: max|F_unit - N F_shell N| "//to_string(maxdiff, "es10.2")// &
                 " relative to max|F| "//to_string(maxdiff/maxf, "es10.2")//" exceeds ~1e-12")
      if (allocated(error)) return

      ! Tr(D F), the two-electron energy expression, agrees between the
      ! two modes' own (D, F) pairs by the same congruence identity (c)
      ! uses: Tr(N^-1 D N^-1 . N F N) = Tr(D F) exactly.
      trdf_shell = sum(real(d_shell, dp)*real(f_shell, dp))
      trdf_unit = sum(real(d_unit, dp)*real(f_unit, dp))
      call check(error, abs(trdf_shell - trdf_unit) < 1.0e-10_dp*max(abs(trdf_shell), 1.0_dp), &
                 "(f) Tr(D F) agrees between shell and unit, |diff| "// &
                 to_string(abs(trdf_shell - trdf_unit), "es10.2"))

      rc = trc_destroy(h_shell); rc = trc_destroy(h_unit)
      deallocate (d_shell, f_shell, d_unit, f_unit, nrm, tmp)
   end subroutine test_fock

end module check_ao_norm_suite

program check_ao_norm
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_ao_norm_suite, only: collect_ao_norm_tests, setup_ao_norm, teardown_ao_norm
   implicit none
   integer :: stat, dev
   type(testsuite_type), allocatable :: testsuites(:)

   dev = trc_bind_device(0)
   call setup_ao_norm()

   stat = 0
   testsuites = [new_testsuite("check_ao_norm", collect_ao_norm_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)

   call teardown_ao_norm()

   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ao_norm
