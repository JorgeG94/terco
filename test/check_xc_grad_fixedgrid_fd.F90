!
! trc_xc_grad's PULAY TERM ALONE, isolated from everything else a full
! SCF-based finite difference mixes in: the grid points and weights are
! built ONCE (`xg0`, from `dft_grid_build`/`trc_xc_grid_t%build` at R0) and
! reused unchanged at every displaced geometry, so only the basis
! functions move -- exactly what trc_xc_grad.F90 differentiates and
! nothing else (no grid response to omit because there is none to omit
! here: the grid never moves in this comparison, on either side). The
! density is a fixed, synthetic, non-SCF matrix (same deterministic
! construction as check_rys_grad_vs_pyscf.F90's own `synth_density`), so
! there is no SCF to reconverge and no SCF-path noise either.
!
! This is the test that caught the bug check_xc_grad_fd.F90's full,
! SCF-based comparison could not pin down on its own: trc_xc_grad's point
! kernels built the force tiles without the grid WEIGHT (every other
! quantity -- collocation, the density gather, the functional itself --
! was right), so a two-atom case gave a gradient three orders of
! magnitude too large while check_xc_grad_fd's translational-invariance
! check, computed the same way for every atom, stayed numerically
! unaffected by the missing common factor and did not flag it. Fixed by
! multiplying the "P"/"Q" tiles by the point's weight before the force
! GEMM, where the energy path's own `zg`/`zv` do the same (trc_xc.F90's
! xc_point_body); this test is what confirms it stays fixed.
!
! H2 (two well-separated, non-collinear-with-any-axis H atoms, so no
! Cartesian direction is accidentally invariant by symmetry) with LDA
! (SVWN) restricted, GGA (PBE) restricted, and GGA (PBE) unrestricted
! (two independent, asymmetric alpha/beta densities, to exercise the
! cross-spin Q1/Q2 terms trc_xc_grad.F90's header describes). Agreement
! is to ~1e-9 absolute on a gradient of order 1: no fitting, no tolerance
! chosen to let a residual through.
!
! Three independent cases below, each its own test-drive unittest, each
! building its own basis/grid.
!
module check_xc_grad_fixedgrid_fd_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc, only: trc_xc_rks, trc_xc_uks, trc_xc_options, trc_xc_result_t
   use trc_xc_grad, only: xc_grad_add, xc_grad_options
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_xc_grad_fixedgrid_fd

   integer, parameter :: NAT = 2
   real(dp), parameter :: H = 1.0e-4_dp
   real(dp), parameter :: ABSTOL = 1.0e-7_dp

contains

   subroutine collect_xc_grad_fixedgrid_fd(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RKS SVWN H2/6-31G*, fixed grid", test_rks_svwn), &
                  new_unittest("RKS PBE H2/6-31G*, fixed grid", test_rks_pbe), &
                  new_unittest("UKS PBE H2/6-31G*, fixed grid", test_uks_pbe) &
                  ]
   end subroutine collect_xc_grad_fixedgrid_fd

   !
   ! Same deterministic, dense, symmetric construction as
   ! check_rys_grad_vs_pyscf.F90/check_rys_grad_fd.F90.
   !
   subroutine synth_density(n, phase, d)
      integer, intent(in) :: n, phase
      real(dp), intent(out) :: d(n, n)
      integer :: p, q
      do q = 1, n
         do p = 1, n
            d(p, q) = 0.4_dp + 0.1_dp*cos(real(3*p + 5*q + 7*phase, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

   subroutine setup(bas, grid0, xg0, r0, at_z)
      type(trc_basis_t), intent(out) :: bas
      type(dft_grid_t), intent(out) :: grid0
      type(trc_xc_grid_t), intent(out) :: xg0
      real(dp), intent(out) :: r0(3, NAT)
      integer, intent(out) :: at_z(NAT)
      type(error_t) :: err
      integer :: dev
      dev = trc_bind_device(0)
      at_z = [1, 1]
      r0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      r0(:, 2) = [0.3_dp, 0.5_dp, 1.4_dp]   ! off every axis: no direction is invariant by symmetry

      call trc_basis_from_json('basis_sets/6-31gs.json', NAT, at_z, r0, bas, err)
      call bas%to_device()

      call dft_grid_build(bas%at_r, nint(bas%at_z), grid0, err, level=3)
      call xg0%build(grid0%n_points, grid0%coords, grid0%weights, bas, max_pts=512)
      call xg0%to_device()
   end subroutine setup

   subroutine test_rks_svwn(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rks_case(error, 'svwn')
   end subroutine test_rks_svwn

   subroutine test_rks_pbe(error)
      type(error_type), allocatable, intent(out) :: error
      call run_rks_case(error, 'pbe')
   end subroutine test_rks_pbe

   subroutine run_rks_case(error, func_name)
      type(error_type), allocatable, intent(out) :: error
      character(*), intent(in) :: func_name
      type(trc_basis_t) :: bas, basp
      type(dft_grid_t) :: grid0
      type(trc_xc_grid_t) :: xg0
      type(trc_linalg_t) :: lp
      type(trc_xc_functional_t) :: func
      real(dp) :: r0(3, NAT)
      integer :: at_z(NAT)
      real(dp), allocatable :: dmat(:, :), grad(:, :)
      real(dp) :: worst
      logical :: ok

      call setup(bas, grid0, xg0, r0, at_z)
      allocate (dmat(bas%nao, bas%nao), grad(3, NAT))
      call synth_density(bas%nao, 1, dmat)
      call lp%init(bas%nao)

      ok = xc_functional_by_name(func_name, func)
      worst = fixedgrid_worst_diff(bas, basp, xg0, lp, func, 1, dmat, r0, at_z, NAT, H, grad)
      print '(a,es10.2)', '  RKS '//trim(func_name)//' H2/6-31G*, fixed grid: worst|analytic-FD| ', worst

      call lp%release()
      call xg0%release()
      call bas%release()
      deallocate (dmat, grad)

      call check(error, ok, "functional not available")
      if (allocated(error)) return
      call check(error, worst <= ABSTOL, "worst|analytic-FD| = "//to_string(worst, "es10.2")// &
                 " exceeds "//to_string(ABSTOL, "es10.2"))
   end subroutine run_rks_case

   subroutine test_uks_pbe(error)
      !! GGA (PBE), unrestricted: two independent, asymmetric densities,
      !! so the cross-spin Q1/Q2 terms (trc_xc_grad.F90's header) are live.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas, basp
      type(dft_grid_t) :: grid0
      type(trc_xc_grid_t) :: xg0
      type(trc_linalg_t) :: lp
      type(trc_xc_functional_t) :: func
      real(dp) :: r0(3, NAT)
      integer :: at_z(NAT)
      real(dp), allocatable :: dmat2(:, :, :), grad2(:, :)
      real(dp) :: worst
      logical :: ok

      call setup(bas, grid0, xg0, r0, at_z)
      allocate (dmat2(bas%nao, bas%nao, 2), grad2(3, NAT))
      call synth_density(bas%nao, 1, dmat2(:, :, 1))
      call synth_density(bas%nao, 2, dmat2(:, :, 2))
      call lp%init(bas%nao)

      ok = xc_functional_by_name('pbe', func)
      worst = fixedgrid_worst_diff2(bas, basp, xg0, lp, func, dmat2, r0, at_z, NAT, H, grad2)
      print '(a,es10.2)', '  UKS PBE H2/6-31G*, fixed grid: worst|analytic-FD| ', worst

      call lp%release()
      call xg0%release()
      call bas%release()
      deallocate (dmat2, grad2)

      call check(error, ok, "functional not available")
      if (allocated(error)) return
      call check(error, worst <= ABSTOL, "worst|analytic-FD| = "//to_string(worst, "es10.2")// &
                 " exceeds "//to_string(ABSTOL, "es10.2"))
   end subroutine test_uks_pbe

   !
   ! RKS: the analytic fixed-grid XC gradient (into `grad`) and the worst
   ! |analytic - central difference| over every atom/direction, the SAME
   ! `xg` reused unperturbed at every displacement.
   !
   function fixedgrid_worst_diff(bas0, basp0, xg, lp0, func0, nspin, d, r00, z, n, hh, grad0) result(w)
      type(trc_basis_t), intent(in) :: bas0
      type(trc_basis_t), intent(inout) :: basp0
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_linalg_t), intent(inout) :: lp0
      type(trc_xc_functional_t), intent(in) :: func0
      integer, intent(in) :: nspin, n, z(n)
      real(dp), intent(in) :: d(bas0%nao, bas0%nao), r00(3, n), hh
      real(dp), intent(inout) :: grad0(3, n)
      real(dp) :: w
      real(dp), allocatable :: vv(:, :)
      real(dp) :: rt(3, n), ep, em, df
      type(trc_xc_result_t) :: xr
      type(error_t) :: e
      integer :: ia, id

      allocate (vv(bas0%nao, bas0%nao))
      grad0 = 0.0_dp
      !$acc enter data copyin(d, grad0)
      call xc_grad_add(bas0, lp0, func0, nspin, d, grad0, &
                       xc_grad_options(grid_level=3, rho_tol=0.0_dp, grid_response=.false.))
      !$acc update self(grad0)
      !$acc exit data delete(d, grad0)

      w = 0.0_dp
      do ia = 1, n
         do id = 1, 3
            rt = r00; rt(id, ia) = r00(id, ia) + hh
            call trc_basis_from_json('basis_sets/6-31gs.json', n, z, rt, basp0, e)
            call basp0%to_device()
            call trc_xc_rks(basp0, xg, func0, d, vv, xr, trc_xc_options(rho_tol=0.0_dp))
            ep = xr%exc
            call basp0%release()
            rt = r00; rt(id, ia) = r00(id, ia) - hh
            call trc_basis_from_json('basis_sets/6-31gs.json', n, z, rt, basp0, e)
            call basp0%to_device()
            call trc_xc_rks(basp0, xg, func0, d, vv, xr, trc_xc_options(rho_tol=0.0_dp))
            em = xr%exc
            call basp0%release()
            df = (ep - em)/(2.0_dp*hh)
            w = max(w, abs(grad0(id, ia) - df))
         end do
      end do
      deallocate (vv)
   end function fixedgrid_worst_diff

   !
   ! UKS counterpart: `d` is (nao, nao, 2).
   !
   function fixedgrid_worst_diff2(bas0, basp0, xg, lp0, func0, d, r00, z, n, hh, grad0) result(w)
      type(trc_basis_t), intent(in) :: bas0
      type(trc_basis_t), intent(inout) :: basp0
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_linalg_t), intent(inout) :: lp0
      type(trc_xc_functional_t), intent(in) :: func0
      integer, intent(in) :: n, z(n)
      real(dp), intent(in) :: d(bas0%nao, bas0%nao, 2), r00(3, n), hh
      real(dp), intent(inout) :: grad0(3, n)
      real(dp) :: w
      real(dp), allocatable :: vv(:, :, :)
      real(dp) :: rt(3, n), ep, em, df
      type(trc_xc_result_t) :: xr
      type(error_t) :: e
      integer :: ia, id

      allocate (vv(bas0%nao, bas0%nao, 2))
      grad0 = 0.0_dp
      !$acc enter data copyin(d, grad0)
      call xc_grad_add(bas0, lp0, func0, 2, d, grad0, &
                       xc_grad_options(grid_level=3, rho_tol=0.0_dp, grid_response=.false.))
      !$acc update self(grad0)
      !$acc exit data delete(d, grad0)

      w = 0.0_dp
      do ia = 1, n
         do id = 1, 3
            rt = r00; rt(id, ia) = r00(id, ia) + hh
            call trc_basis_from_json('basis_sets/6-31gs.json', n, z, rt, basp0, e)
            call basp0%to_device()
            call trc_xc_uks(basp0, xg, func0, d, vv, xr, trc_xc_options(rho_tol=0.0_dp))
            ep = xr%exc
            call basp0%release()
            rt = r00; rt(id, ia) = r00(id, ia) - hh
            call trc_basis_from_json('basis_sets/6-31gs.json', n, z, rt, basp0, e)
            call basp0%to_device()
            call trc_xc_uks(basp0, xg, func0, d, vv, xr, trc_xc_options(rho_tol=0.0_dp))
            em = xr%exc
            call basp0%release()
            df = (ep - em)/(2.0_dp*hh)
            w = max(w, abs(grad0(id, ia) - df))
         end do
      end do
      deallocate (vv)
   end function fixedgrid_worst_diff2

end module check_xc_grad_fixedgrid_fd_suite

program check_xc_grad_fixedgrid_fd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_xc_grad_fixedgrid_fd_suite, only: collect_xc_grad_fixedgrid_fd
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_xc_grad_fixedgrid_fd, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_xc_grad_fixedgrid_fd
