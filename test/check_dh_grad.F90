!
! The analytic double-hybrid nuclear gradient (src/post_hf/rimp2/trc_rimp2_grad.F90
! over a Kohn-Sham reference: E = E_KS + c E_RI-MP2 on the KS orbitals)
! against:
!
! A. Central finite differences of terco's OWN double-hybrid energy
!    (trc_scf_run with the functional, then trc_rimp2_run, c E_corr: the
!    pair test/scf_dh.F90 anchors to pyscf's B2PLYP/B2GP-PLYP/mPW2-PLYP +
!    DF-MP2 at 1e-9..1e-8), Richardson-extrapolated over steps 2h and h:
!    (4 D(h) - D(2h))/3 removes the h^2 truncation a DFT energy's step has
!    to live with. The grid moves with the atoms at every displaced
!    geometry, as the analytic grid response assumes. --df (the priority:
!    Cartesian with one set doing both fits, spherical with
!    def2-universal-jkfit for the SCF) and the exact-ERI reference.
!
! B. pyscf 2.14 (/home/jorge/install/miniconda3/envs/gpu4pyscf): Richardson
!    finite differences (h = 1e-3 Bohr, the same formula) of pyscf's own
!    density-fitted RKS with the double hybrid's semilocal composition
!    (test/pyscf_refs/generate.py's DH_HYB strings, e.g. '0.53*HF +
!    0.47*B88, 0.73*LYP'; grids.level 3, the grid terco builds, conv_tol
!    1e-13, conv_tol_grad 1e-9) plus c DF-MP2 (mp.dfmp2.DFMP2, with_df on
!    the same auxiliary basis, the same frozen core), from the same basis
!    JSON files (generate.py's bse_json_basis) and the same geometry, in
!    Bohr through terco's CODATA 2018 constant (pyscf's own is CODATA
!    2010). PySCF has no double-hybrid gradient of its own. Its energy
!    agrees with terco's to 2e-11 (water) at the reference geometry; its
!    finite-difference noise shows in the components symmetry makes zero
!    (up to 4e-10 for water, 2e-9 for formamide), which sets the gates.
!
! C. The identity that defines the new exchange-correlation term (src/xc/
!    trc_xc_grad_potential.F90): E_xc[D + t P] = E_xc[D] + t Tr(P V_xc[D]) +
!    O(t^2) at every geometry, so d/dR Tr(P V_xc[D]) is the t-derivative of
!    trc_xc_grad's own gradient at D + t P, grid response included.
!
! Plus the --df Z-vector's AO fallback against its MO path, and the refusals
! the KS reference brings (a range-separated functional). Every gradient
! case also checks translational invariance and that the gradient's
! amplitudes give the energy's correlation energy.
!
! A test-drive suite.
!
module check_dh_grad_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char, c_char
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t, trc_rimp2_options_t, rimp2_driver_mo_to_cart
   use trc_rimp2_grad, only: trc_rimp2_gradient, trc_rimp2_grad_t
   use trc_scf_df, only: df_jk_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_grad, only: xc_grad_add, xc_grad_options, xc_grad_options_t
   use trc_xc_grad_potential, only: xc_grad_potential_add
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, trc_set_rimp2, &
                               trc_run_scf, trc_run_rimp2, trc_run_gradient, trc_gradient, trc_set_df, &
                               trc_set_jkfit_json, trc_set_comm, trc_set_lindep, trc_message, TRC_OK, &
                               TRC_ERR_UNSUPPORTED
   implicit none
   private
   public :: collect_dh_grad, collect_dh_grad_long

   real(dp), parameter :: FD_STEP = 1.0e-3_dp
   !! h of the Richardson pair (2h, h), Bohr.
   real(dp), parameter :: FD_TOL = 1.0e-8_dp
   !! Measured worst 2.4e-10 (CPU); the plain central difference at h is
   !! off by up to 2.8e-7, so this gate still needs the extrapolation.
   real(dp), parameter :: PYSCF_TOL = 5.0e-9_dp
   !! Water: measured worst 7.6e-10 against pyscf's finite differences.
   real(dp), parameter :: PYSCF_TOL_FORMAMIDE = 1.0e-8_dp
   !! Formamide: pyscf's own finite-difference noise reaches 2e-9.
   real(dp), parameter :: G_PYSCF_WATER(3, 3) = reshape([ &
                                              -1.894780628694e-11_dp, -1.515824502955e-10_dp, -1.332378636339e-02_dp, &
                                              3.552713678801e-10_dp, -3.896076487801e-03_dp, 6.661893121892e-03_dp, &
                                              1.421085471520e-11_dp, 3.896076588461e-03_dp, 6.661893099391e-03_dp], &
                                              [3, 3])
   !! pyscf: water.xyz, Cartesian cc-pVDZ, b2plyp, cc-pVDZ-RIFIT for the
   !! SCF fit and the DF-MP2, frozen 1. E = -76.355985445232.
   real(dp), parameter :: E_PYSCF_WATER = -76.355985445232_dp
   real(dp), parameter :: G_PYSCF_FORMAMIDE(3, 6) = reshape([ &
                                                  -2.772253125253e-02_dp, -1.823512976766e-02_dp, -1.622405913319e-09_dp, &
                                                  3.639436725204e-02_dp, 3.206997381729e-04_dp, -1.373715955803e-09_dp, &
                                                  2.182021750248e-03_dp, 1.975698219780e-02_dp, 9.544957417044e-10_dp, &
                                                  -9.098323412597e-03_dp, 6.270324886278e-03_dp, 5.542233338929e-10_dp, &
                                                  -6.171912550220e-03_dp, -1.720528596631e-03_dp, -2.998490344908e-09_dp, &
                                                  4.416377469359e-03_dp, -6.392352448851e-03_dp, -3.417710559006e-09_dp], &
                                                  [3, 6])
   !! pyscf: formamide.xyz, spherical cc-pVDZ, b2gp-plyp, def2-universal-
   !! jkfit for the SCF fit and cc-pVDZ-RIFIT for the DF-MP2, frozen 3.
   real(dp), parameter :: E_PYSCF_FORMAMIDE = -169.698488127896_dp
   real(dp), parameter :: ECORR_TOL = 1.0e-11_dp
   real(dp), parameter :: TRANS_TOL = 1.0e-9_dp
   real(dp), parameter :: IDENTITY_STEP = 1.0e-3_dp
   real(dp), parameter :: IDENTITY_TOL = 1.0e-9_dp
   real(dp), parameter :: AO_MO_TOL = 1.0e-10_dp

   type :: case_t
      character(len=64) :: name = ""
      character(len=64) :: xyz = ""
      character(len=64) :: basis = ""
      character(len=64) :: aux = ""
      character(len=64) :: jkfit = ""
      !! The --df reference's fitting set; empty for the exact reference.
      character(len=16) :: functional = ""
      real(dp) :: pt2 = 0.0_dp
      !! The literal PT2 fraction (c_os = c_ss), never read back from terco's table.
      integer :: nfrozen = 0
      logical :: cart = .false.
   end type case_t

   type :: system_t
      type(trc_basis_t) :: bas, aux, jk
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: lp
      type(df_jk_t) :: dfjk
      type(trc_scf_result_t) :: res
      integer :: nocc = 0
      logical :: df = .false.
   end type system_t

contains

   subroutine collect_dh_grad(testsuite)
      !! The default tier: every new term once, water only.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("XC: d/dR Tr(P V_xc[D]) is the t-derivative of the XC gradient at D + tP", &
                               test_xc_potential_identity), &
                  new_unittest("--df: water cc-pVDZ (Cartesian) b2plyp, frozen 1, vs pyscf", test_pyscf_water), &
                  new_unittest("--df: water cc-pVDZ (Cartesian) b2plyp, frozen 1, vs finite differences", &
                               test_fd_df_water_cart), &
                  new_unittest("--df: the AO fallback of the fitted KS orbital Hessian agrees with the MO one", &
                               test_df_response_ao), &
                  new_unittest("a range-separated KS reference is refused with its reason", test_refuses_rsh), &
                  new_unittest("C ABI: trc_run_gradient after b2plyp + RI-MP2 (--df) matches pyscf", test_capi_gradient), &
                  new_unittest("C ABI: the double-hybrid gradient is refused on more than one rank", &
                               test_capi_refuses_multirank), &
                  new_unittest("C ABI: the double-hybrid gradient is refused with functions dropped", &
                               test_capi_refuses_lindep) &
                  ]
   end subroutine collect_dh_grad

   subroutine collect_dh_grad_long(testsuite)
      !! The long tier (ctest label LONG): a molecule without symmetry,
      !! spherical, against pyscf and finite differences, and the exact
      !! reference.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("--df: formamide cc-pVDZ (sph) b2gp-plyp, universal-jkfit, frozen 3, vs pyscf", &
                               test_pyscf_formamide), &
                  new_unittest("--df: formamide cc-pVDZ (sph) b2gp-plyp, universal-jkfit, frozen 3, vs finite "// &
                               "differences", test_fd_df_formamide_sph), &
                  new_unittest("exact: water cc-pVDZ (Cartesian) mpw2plyp, frozen 1, vs finite differences", &
                               test_fd_exact_water_cart), &
                  new_unittest("--df: water cc-pVDZ (Cartesian) PBE + RI-MP2 (c = 1, a_x = 0) vs finite differences", &
                               test_fd_df_water_pbe) &
                  ]
   end subroutine collect_dh_grad_long

   ! ---- the C ABI ------------------------------------------------------------

   subroutine capi_water_dh(error, h, lindep)
      !! A handle through trc_run_rimp2: water, Cartesian cc-pVDZ, b2plyp,
      !! --df with cc-pVDZ-RIFIT doing both fits, frozen 1 -- the pyscf case.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(out) :: h
      real(dp), intent(in), optional :: lindep
      integer :: natm
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer(c_int) :: rc
      call read_xyz("water.xyz", natm, zint, at_r)
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, int(natm, c_int), real(zint, c_double), reshape(at_r, [3*natm]), &
                                              0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_force_cartesian(h, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, "basis_sets/cc-pvdz.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_aux_json(h, "basis_sets/cc-pvdz-rifit.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_jkfit_json(h, "basis_sets/cc-pvdz-rifit.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_df(h, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_method(h, "b2plyp"//c_null_char, 3_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 10_c_int)
      if (rc == TRC_OK .and. present(lindep)) rc = trc_set_lindep(h, real(lindep, c_double))
      if (rc == TRC_OK) rc = trc_set_rimp2(h, 1_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_run_rimp2(h)
      call check(error, rc == TRC_OK, "C ABI setup through trc_run_rimp2 failed: "//message_of(h))
   end subroutine capi_water_dh

   function message_of(h) result(text)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char, len=1024) :: buffer
      integer(c_int) :: rc
      integer :: k
      buffer = ""
      rc = trc_message(h, buffer, 1024_c_int)
      k = index(buffer, c_null_char)
      if (k == 0) k = len_trim(buffer) + 1
      text = buffer(1:k - 1)
   end function message_of

   subroutine test_capi_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(c_double) :: gbuf(9)
      real(dp) :: worst
      integer(c_int) :: rc
      call capi_water_dh(error, h)
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      call check(error, rc == TRC_OK, "trc_run_gradient: "//message_of(h))
      if (allocated(error)) return
      rc = trc_gradient(h, gbuf)
      worst = maxval(abs(reshape(real(gbuf, dp), [3, 3]) - G_PYSCF_WATER))
      print "(a,es10.2)", "  C ABI b2plyp --df: max |g - g_pyscf| ", worst
      rc = trc_destroy(h)
      call check(error, worst < PYSCF_TOL, "the C ABI's double-hybrid gradient differs from pyscf's")
   end subroutine test_capi_gradient

   subroutine test_capi_refuses_multirank(error)
      !! trc_set_comm with any handle but -1 makes the context multi-rank;
      !! the refusal comes before any MPI call, so a serial build tests it.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(len=:), allocatable :: why
      call capi_water_dh(error, h)
      if (allocated(error)) return
      rc = trc_set_comm(h, 0_c_int)
      rc = trc_run_gradient(h)
      why = message_of(h)
      rc = trc_destroy(h)
      call check(error, index(why, "one rank only") > 0, "not refused, or for the wrong reason: '"//why//"'")
   end subroutine test_capi_refuses_multirank

   subroutine test_capi_refuses_lindep(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(len=:), allocatable :: why
      call capi_water_dh(error, h, lindep=0.05_dp)
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      why = message_of(h)
      rc = trc_destroy(h)
      call check(error, index(why, "needs every basis function") > 0, &
                 "not refused, or for the wrong reason: '"//why//"'")
   end subroutine test_capi_refuses_lindep

   ! ---- B. pyscf --------------------------------------------------------------

   subroutine test_pyscf_water(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cart --df b2plyp fc1"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/cc-pvdz-rifit.json"
      c%functional = "b2plyp"
      c%pt2 = 0.27_dp
      c%nfrozen = 1
      c%cart = .true.
      call against_pyscf(error, c, E_PYSCF_WATER, G_PYSCF_WATER, PYSCF_TOL)
   end subroutine test_pyscf_water

   subroutine test_pyscf_formamide(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "formamide sph --df b2gp-plyp fc3"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/def2-universal-jkfit.json"
      c%functional = "b2gp-plyp"
      c%pt2 = 0.36_dp
      c%nfrozen = 3
      call against_pyscf(error, c, E_PYSCF_FORMAMIDE, G_PYSCF_FORMAMIDE, PYSCF_TOL_FORMAMIDE)
   end subroutine test_pyscf_formamide

   subroutine against_pyscf(error, c, e_ref, g_ref, tol)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      real(dp), intent(in) :: e_ref, g_ref(:, :), tol
      real(dp), allocatable :: grad(:, :), r(:, :)
      integer, allocatable :: z(:)
      integer :: nat
      real(dp) :: e_tot, worst

      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, grad, e_tot)
      if (allocated(error)) return
      worst = maxval(abs(grad - g_ref))
      print "(a,a,es10.2,a,es10.2,a,es10.2)", "  ", trim(c%name)//": |E - E_pyscf| ", abs(e_tot - e_ref), &
         "  max |g - g_pyscf| ", worst, "  |sum g| ", maxval(abs(sum(grad, dim=2)))
      call check(error, abs(e_tot - e_ref) < 1.0e-8_dp, trim(c%name)// &
                 ": the double-hybrid energy differs from pyscf's -- not the same energy")
      if (allocated(error)) return
      call check(error, worst < tol, trim(c%name)//": gradient differs from pyscf's finite differences")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(grad, dim=2))) < TRANS_TOL, trim(c%name)//": not translationally invariant")
   end subroutine against_pyscf

   ! ---- C. the exchange-correlation identity --------------------------------

   subroutine test_xc_potential_identity(error)
      !! P = C_o M C_o^T, M symmetric with both signs (indefinite, as a
      !! relaxed density is) and inside the occupied span, so D + t P keeps
      !! the reference's tails: the difference stays where the functional is
      !! smooth.
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      type(system_t) :: sys
      type(trc_xc_functional_t) :: func
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(xc_grad_options_t) :: xopts
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), cocc(:, :), mmat(:, :), dens(:, :), pmat(:, :), dpm(:, :, :)
      real(dp), allocatable :: g_an(:, :), g_p(:, :), g_m(:, :), g_p2(:, :), g_m2(:, :), g_fd(:, :)
      integer, allocatable :: z(:), at_elem(:)
      integer :: nat, n, i, j
      real(dp) :: e_tot, h, worst

      c%name = "water cart b2plyp identity"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%functional = "b2plyp"
      c%cart = .true.
      call read_xyz(trim(c%xyz), nat, z, r)
      call setup_system(error, c, nat, z, r, sys, e_tot, with_pt2=.false.)
      if (allocated(error)) return
      call check(error, xc_functional_by_name("b2plyp", func), "b2plyp is not registered")
      if (allocated(error)) return
      n = sys%bas%nao
      cocc = sys%res%cmo(:, 1:sys%nocc, 1)
      ! Each orbital's phase fixed (largest coefficient positive), so P --
      ! not invariant to them, M being off-diagonal -- is the same P on
      ! every run.
      do j = 1, sys%nocc
         if (cocc(maxloc(abs(cocc(:, j)), dim=1), j) < 0.0_dp) cocc(:, j) = -cocc(:, j)
      end do
      allocate (mmat(sys%nocc, sys%nocc))
      do j = 1, sys%nocc
         do i = 1, sys%nocc
            mmat(i, j) = 0.3_dp*cos(real(i + 2*j, dp)) + 0.3_dp*cos(real(j + 2*i, dp))
         end do
      end do
      dens = 2.0_dp*matmul(cocc, transpose(cocc))
      pmat = matmul(cocc, matmul(mmat, transpose(cocc)))

      ! The analytic term on the level-3 grid trc_xc_grad builds for itself.
      allocate (at_elem(nat), source=sys%bas%at_elem)
      call dft_grid_build(sys%bas%at_r, at_elem, grid, err, level=3)
      call check(error, .not. err%has_error(), "grid: "//err%get_message())
      if (allocated(error)) return
      call xg%build_from_grid(grid, sys%bas, max_pts=512)
      call xg%to_device()
      xopts = xc_grad_options(grid_level=3, grid_max_pts=512, grid_response=.true., partition=grid%partition)
      allocate (g_an(3, nat), g_p(3, nat), g_m(3, nat), g_p2(3, nat), g_m2(3, nat), dpm(n, n, 1))
      g_an = 0.0_dp
      !$acc enter data copyin(g_an)
      call xc_grad_potential_add(sys%bas, sys%lp, func, xg, dens, pmat, g_an, xopts)
      !$acc update self(g_an)
      !$acc exit data delete(g_an)
      call xg%release()
      call grid%destroy()

      h = IDENTITY_STEP
      call xc_gradient_at(sys, func, dens + h*pmat, g_p)
      call xc_gradient_at(sys, func, dens - h*pmat, g_m)
      call xc_gradient_at(sys, func, dens + 2.0_dp*h*pmat, g_p2)
      call xc_gradient_at(sys, func, dens - 2.0_dp*h*pmat, g_m2)
      ! Richardson over t-steps 2h and h.
      g_fd = (4.0_dp*(g_p - g_m)/(2.0_dp*h) - (g_p2 - g_m2)/(4.0_dp*h))/3.0_dp
      worst = maxval(abs(g_an - g_fd))
      print "(a,a,es10.2,a,es10.2)", "  ", trim(c%name)//": max |analytic - d/dt| ", worst, &
         "  max |term| ", maxval(abs(g_an))
      call check(error, maxval(abs(g_an)) > 1.0e-4_dp, trim(c%name)//": the term is too small to test anything")
      if (allocated(error)) return
      call check(error, worst < IDENTITY_TOL, trim(c%name)// &
                 ": d/dR Tr(P V_xc[D]) differs from the t-derivative of the XC gradient")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(g_an, dim=2))) < TRANS_TOL, trim(c%name)//": not translationally invariant")
      call release_system(sys)
   end subroutine test_xc_potential_identity

   subroutine xc_gradient_at(sys, func, dens, grad)
      !! trc_xc_grad's own E_xc gradient (grid response on) at the density
      !! `dens`, on its own level-3 grid.
      type(system_t), intent(inout) :: sys
      type(trc_xc_functional_t), intent(in) :: func
      real(dp), intent(in) :: dens(:, :)
      real(dp), intent(out) :: grad(:, :)
      real(dp), allocatable :: d3(:, :, :)
      allocate (d3(size(dens, 1), size(dens, 2), 1))
      d3(:, :, 1) = dens
      grad = 0.0_dp
      !$acc enter data copyin(d3, grad)
      call xc_grad_add(sys%bas, sys%lp, func, 1, d3, grad, xc_grad_options(grid_level=3, grid_max_pts=512))
      !$acc update self(grad)
      !$acc exit data delete(d3, grad)
   end subroutine xc_gradient_at

   ! ---- A. finite differences of terco's own energy ------------------------

   subroutine test_fd_df_water_cart(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cart --df b2plyp fc1"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/cc-pvdz-rifit.json"
      c%functional = "b2plyp"
      c%pt2 = 0.27_dp
      c%nfrozen = 1
      c%cart = .true.
      ! O z, H y, H z.
      call against_fd(error, c, reshape([3, 1, 2, 2, 3, 3], [2, 3]))
   end subroutine test_fd_df_water_cart

   subroutine test_fd_df_formamide_sph(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "formamide sph --df b2gp-plyp fc3"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/def2-universal-jkfit.json"
      c%functional = "b2gp-plyp"
      c%pt2 = 0.36_dp
      c%nfrozen = 3
      ! C y, N x, an amide H y: no symmetry relates them.
      call against_fd(error, c, reshape([2, 1, 1, 3, 2, 5], [2, 3]))
   end subroutine test_fd_df_formamide_sph

   subroutine test_fd_exact_water_cart(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cart exact mpw2plyp fc1"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%functional = "mpw2plyp"
      c%pt2 = 0.25_dp
      c%nfrozen = 1
      c%cart = .true.
      call against_fd(error, c, reshape([3, 1, 2, 2, 3, 3], [2, 3]))
   end subroutine test_fd_exact_water_cart

   subroutine test_fd_df_water_pbe(error)
      !! Not a double hybrid: a pure GGA's E_KS + E_RI-MP2, which trc_energy
      !! reports after trc_run_rimp2 -- c = 1, and no exact exchange, so
      !! every K scale in the response and the 2e derivative is zero.
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cart --df pbe + RI-MP2 fc1"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/cc-pvdz-rifit.json"
      c%functional = "pbe"
      c%pt2 = 1.0_dp
      c%nfrozen = 1
      c%cart = .true.
      call against_fd(error, c, reshape([3, 1, 2, 2, 3, 3], [2, 3]))
   end subroutine test_fd_df_water_pbe

   subroutine test_df_response_ao(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      real(dp), allocatable :: r(:, :), g_mo(:, :), g_ao(:, :)
      integer, allocatable :: z(:)
      integer :: nat
      real(dp) :: e_tot
      c%name = "water sph --df b2plyp AO fallback"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/def2-universal-jkfit.json"
      c%functional = "b2plyp"
      c%pt2 = 0.27_dp
      c%nfrozen = 1
      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, g_mo, e_tot)
      if (allocated(error)) return
      call analytic_gradient(error, c, nat, z, r, g_ao, e_tot, force_ao=.true.)
      if (allocated(error)) return
      print "(a,a,es10.2)", "  ", trim(c%name)//": max |g(AO) - g(MO)| ", maxval(abs(g_ao - g_mo))
      call check(error, maxval(abs(g_ao - g_mo)) < AO_MO_TOL, trim(c%name)// &
                 ": the two forms of the fitted KS orbital Hessian give different gradients")
   end subroutine test_df_response_ao

   subroutine test_refuses_rsh(error)
      !! wb97x has no CPKS response here (trc_cpks's cpks_refuse_reason): the
      !! gradient must refuse with that reason, not run the wrong physics.
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      type(system_t) :: sys
      type(trc_rimp2_grad_t) :: job
      real(dp), allocatable :: r(:, :)
      integer, allocatable :: z(:)
      integer :: nat
      real(dp) :: e_tot
      c%name = "water wb97x refusal"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%functional = "wb97x"
      c%pt2 = 1.0_dp
      c%cart = .true.
      call read_xyz(trim(c%xyz), nat, z, r)
      call setup_system(error, c, nat, z, r, sys, e_tot, with_pt2=.false.)
      if (allocated(error)) return
      job%opts%nocc = sys%nocc
      call check(error, xc_functional_by_name("wb97x", job%opts%func), "wb97x is not registered")
      if (allocated(error)) return
      call trc_rimp2_gradient(sys%bas, sys%aux, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), job, sys%lp)
      call release_system(sys)
      call check(error, .not. job%res%ok, trim(c%name)//": a range-separated reference was not refused")
      if (allocated(error)) return
      call check(error, index(job%res%message, "range-separated") > 0, trim(c%name)// &
                 ": refused for the wrong reason: "//trim(job%res%message))
   end subroutine test_refuses_rsh

   subroutine against_fd(error, c, comps)
      !! `comps(:, k)` = (component, atom) of the k-th differenced coordinate.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: comps(:, :)
      real(dp), allocatable :: grad(:, :), r(:, :)
      integer, allocatable :: z(:)
      integer :: nat, k, ic, ia
      real(dp) :: e_tot, fd, worst, d1, d2

      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, grad, e_tot)
      if (allocated(error)) return
      worst = 0.0_dp
      do k = 1, size(comps, 2)
         ic = comps(1, k)
         ia = comps(2, k)
         call central_difference(error, c, nat, z, r, ic, ia, FD_STEP, d1)
         if (allocated(error)) return
         call central_difference(error, c, nat, z, r, ic, ia, 2.0_dp*FD_STEP, d2)
         if (allocated(error)) return
         fd = (4.0_dp*d1 - d2)/3.0_dp
         print "(a,a,i0,a,i0,a,f16.10,a,f16.10,a,es9.2,a,es9.2)", "  ", trim(c%name)//": atom ", ia, &
            " comp ", ic, "  analytic ", grad(ic, ia), "  fd ", fd, "  diff ", abs(grad(ic, ia) - fd), &
            "  (plain h: ", abs(grad(ic, ia) - d1)
         worst = max(worst, abs(grad(ic, ia) - fd))
      end do
      call check(error, worst < FD_TOL, trim(c%name)//": analytic gradient differs from finite differences")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(grad, dim=2))) < TRANS_TOL, trim(c%name)//": not translationally invariant")
   end subroutine against_fd

   subroutine central_difference(error, c, nat, z, r, ic, ia, h, d)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat), ic, ia
      real(dp), intent(in) :: r(3, nat), h
      real(dp), intent(out) :: d
      real(dp), allocatable :: rd(:, :)
      real(dp) :: e_plus, e_minus
      d = 0.0_dp
      rd = r
      rd(ic, ia) = r(ic, ia) + h
      call total_energy(error, c, nat, z, rd, e_plus)
      if (allocated(error)) return
      rd(ic, ia) = r(ic, ia) - h
      call total_energy(error, c, nat, z, rd, e_minus)
      if (allocated(error)) return
      d = (e_plus - e_minus)/(2.0_dp*h)
   end subroutine central_difference

   ! ---- shared ------------------------------------------------------------

   subroutine setup_system(error, c, nat, z, r, sys, e_tot, with_pt2)
      !! The bases, the converged KS reference (four-centre, or --df over
      !! `c%jkfit`) and E_KS + c E_RI-MP2 at geometry r, the way
      !! test/scf_dh.F90 computes them; `with_pt2` false stops after the SCF.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      type(system_t), intent(inout) :: sys
      real(dp), intent(out) :: e_tot
      logical, intent(in) :: with_pt2
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      real(dp) :: e_corr

      e_tot = 0.0_dp
      sys%nocc = sum(z)/2
      sys%df = len_trim(c%jkfit) > 0
      call trc_basis_from_json(trim(c%basis), nat, z, r, sys%bas, err, force_cart=c%cart, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (len_trim(c%aux) > 0 .and. .not. err%has_error()) then
         call trc_basis_from_json(trim(c%aux), nat, z, r, sys%aux, err, force_cart=c%cart, &
                                  overrides=trc_basis_json_overrides(prune=.false.))
      end if
      if (sys%df .and. .not. err%has_error()) then
         call trc_basis_from_json(trim(c%jkfit), nat, z, r, sys%jk, err, force_cart=c%cart, &
                                  overrides=trc_basis_json_overrides(prune=.false.))
      end if
      call check(error, .not. err%has_error(), trim(c%name)//": basis: "//err%get_message())
      if (allocated(error)) return
      call sys%bas%to_device()
      if (len_trim(c%aux) > 0) call sys%aux%to_device()
      call sys%pl%build(sys%bas, 1.0e-12_dp)
      call sys%pl%to_device()
      call sys%lp%init(sys%bas%nao)
      opts%functional = trim(c%functional)
      opts%grid_level = 3
      opts%conv_energy = 1.0e-13_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      if (sys%df) then
         call sys%jk%to_device()
         call sys%dfjk%setup(sys%bas, sys%jk, sys%pl, sys%lp, err, nocc_hint=sys%nocc)
         call check(error, .not. err%has_error(), trim(c%name)//": fit: "//err%get_message())
         if (allocated(error)) return
         opts%df = .true.
         call trc_scf_run(sys%bas, sys%nocc, sys%nocc, opts, sys%res, la=sys%lp, scf_dfjk=sys%dfjk)
      else
         call trc_scf_run(sys%bas, sys%nocc, sys%nocc, opts, sys%res, la=sys%lp)
      end if
      call check(error, sys%res%converged, trim(c%name)//": the SCF did not converge: "//trim(sys%res%message))
      if (allocated(error)) return
      call check(error, sys%res%n_lindep == 0, trim(c%name)//": functions dropped; the gradient needs them all")
      if (allocated(error)) return
      e_tot = sys%res%energy
      if (.not. with_pt2) return
      call rimp2_energy(error, c, sys, e_corr)
      if (allocated(error)) return
      e_tot = sys%res%energy + c%pt2*e_corr
   end subroutine setup_system

   subroutine release_system(sys)
      type(system_t), intent(inout) :: sys
      if (sys%df) then
         call sys%dfjk%release()
         call sys%jk%release()
      end if
      call sys%lp%release()
      call sys%pl%release()
      call sys%aux%release()
      call sys%bas%release()
   end subroutine release_system

   subroutine total_energy(error, c, nat, z, r, e_tot)
      !! E_KS + c E_RI-MP2 at geometry r.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      real(dp), intent(out) :: e_tot
      type(system_t) :: sys
      call setup_system(error, c, nat, z, r, sys, e_tot, with_pt2=.true.)
      if (allocated(error)) return
      call release_system(sys)
   end subroutine total_energy

   subroutine rimp2_energy(error, c, sys, e_corr)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      type(system_t), intent(inout) :: sys
      real(dp), intent(out) :: e_corr
      type(trc_rimp2_result_t) :: mp
      type(trc_rimp2_options_t) :: ropts
      real(dp), allocatable :: cmo_c(:, :), eps_c(:)
      allocate (cmo_c(sys%bas%nao, sys%bas%nao), eps_c(sys%bas%nao))
      if (sys%bas%spherical) then
         call rimp2_driver_mo_to_cart(sys%bas, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), cmo_c, eps_c)
      else
         cmo_c = sys%res%cmo(:, :, 1)
         eps_c = sys%res%eps(:, 1)
      end if
      ropts%nfrozen = c%nfrozen
      call trc_rimp2_run(sys%bas, sys%aux, sys%pl, sys%nocc, cmo_c, eps_c, mp, opts=ropts, la=sys%lp)
      e_corr = mp%e_corr
      call check(error, len_trim(mp%message) == 0, trim(c%name)//": RI-MP2: "//trim(mp%message))
   end subroutine rimp2_energy

   subroutine analytic_gradient(error, c, nat, z, r, grad, e_tot, force_ao)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      real(dp), allocatable, intent(out) :: grad(:, :)
      real(dp), intent(out) :: e_tot
      logical, intent(in), optional :: force_ao
      type(system_t) :: sys
      type(trc_rimp2_grad_t) :: job
      real(dp) :: e_corr

      call setup_system(error, c, nat, z, r, sys, e_tot, with_pt2=.true.)
      if (allocated(error)) return
      e_corr = (e_tot - sys%res%energy)/c%pt2
      job%opts%nocc = sys%nocc
      job%opts%nfrozen = c%nfrozen
      job%opts%pt2_scale = c%pt2
      job%opts%grid_level = 3
      call check(error, xc_functional_by_name(trim(c%functional), job%opts%func), &
                 trim(c%name)//": functional not registered")
      if (allocated(error)) return
      if (present(force_ao)) job%opts%df_response_ao = force_ao
      if (sys%df) then
         call trc_rimp2_gradient(sys%bas, sys%aux, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), job, sys%lp, &
                                 dfjk=sys%dfjk, jkaux=sys%jk)
      else
         call trc_rimp2_gradient(sys%bas, sys%aux, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), job, sys%lp)
      end if
      call check(error, job%res%ok, trim(c%name)//": "//trim(job%res%message))
      if (allocated(error)) return
      grad = job%grad
      print "(a,a,f20.12,a,es10.2,a,i0)", "  ", trim(c%name)//": E_DH ", e_tot, &
         "  E_corr (gradient) - E_corr (energy) ", job%res%e_corr - e_corr, "  Z-vector iterations ", &
         job%res%cphf_iterations
      call check(error, abs(job%res%e_corr - e_corr) < ECORR_TOL, trim(c%name)// &
                 ": the gradient's amplitudes give a different correlation energy than trc_rimp2_run")
      call release_system(sys)
   end subroutine analytic_gradient

end module check_dh_grad_suite

program check_dh_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_dh_grad_suite, only: collect_dh_grad, collect_dh_grad_long
   implicit none
   character(len=16) :: tier
   integer :: stat
   ! No argument: the default tier; "long": the LONG-labelled one.
   call get_command_argument(1, tier)
   stat = 0
   if (trim(tier) == "long") then
      call run_testsuite(collect_dh_grad_long, error_unit, stat, parallel=.false.)
   else
      call run_testsuite(collect_dh_grad, error_unit, stat, parallel=.false.)
   end if
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_dh_grad
