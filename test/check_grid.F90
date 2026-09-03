!
! The molecular quadrature: closed forms and identities.
!
! The grid modules came from metalquicha with a test-drive suite of fifty
! cases. terco does not carry test-drive, so this is the subset of those
! cases that a wrong grid cannot pass, as one program:
!
!   * a Lebedev order integrates low polynomials over the sphere exactly,
!     its weights sum to one, and every point has its antipode -- the
!     quadrature property, which a point set that merely lies on a sphere
!     would fail;
!   * the Becke cell weights sum to one at every point, which fails for an
!     asymmetric pair loop, a size adjustment with the wrong sign, or an
!     atom dropped from a product;
!   * the radial mesh integrates exp(-r^2) exactly, probing the Jacobian;
!   * the assembled grid integrates a sum of normalised Gaussians and of
!     normalised Slaters over water to the atom count, and refining the
!     level reduces the error;
!   * bad requests are refused rather than half-applied.
!
! Every tolerance is one the metalquicha suite measured, not one chosen to
! pass: level 3 reaches 5.7e-9 on the Gaussian sum and the bound is 1e-8.
!
program check_grid
   use pic_types, only: dp
   use trc_error, only: error_t
   use trc_lebedev, only: lebedev_grid, lebedev_available_orders, lebedev_has_negative_weights
   use trc_dft_radial, only: treutler_ahlrichs_radial, radial_volume_weights
   use trc_dft_partition, only: becke_partition_weights, PARTITION_BECKE, PARTITION_STRATMANN, &
                                ADJUST_TREUTLER, ADJUST_NONE
   use trc_dft_prune, only: PRUNE_NONE
   use trc_dft_grid, only: dft_grid_t, build_dft_grid, grid_level_radial, grid_level_angular
   use trc_test_basis, only: read_xyz
   implicit none

   real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
   integer :: nfail

   nfail = 0
   call lebedev_checks()
   call radial_checks()
   call partition_checks()
   call grid_checks()
   call screening_checks()

   if (nfail == 0) then
      print '(a)', "check_grid: all checks passed"
   else
      print '(a,i0,a)', "check_grid: ", nfail, " check(s) FAILED"
      stop 1
   end if

contains

   subroutine report(ok, what, value)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: what
      real(dp), intent(in), optional :: value
      if (.not. ok) nfail = nfail + 1
      if (present(value)) then
         print '(a,1x,a,t56,es10.2)', merge("ok  ", "FAIL", ok), what, value
      else
         print '(a,1x,a)', merge("ok  ", "FAIL", ok), what
      end if
   end subroutine report

   subroutine water(coords, z)
      real(dp), allocatable, intent(out) :: coords(:, :)
      integer, allocatable, intent(out) :: z(:)
      allocate (coords(3, 3), z(3))
      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
      z = [8, 1, 1]
   end subroutine water

   ! ------------------------------------------------------------------------

   subroutine lebedev_checks()
      integer, allocatable :: orders(:)
      real(dp), allocatable :: p(:, :), w(:)
      type(error_t) :: err
      integer :: io, n, k, j
      real(dp) :: worst_sum, worst_sphere, worst_moment, worst_antipode, d
      logical :: negative_ok

      orders = lebedev_available_orders()
      worst_sum = 0.0_dp; worst_sphere = 0.0_dp; worst_moment = 0.0_dp; worst_antipode = 0.0_dp
      negative_ok = .true.
      do io = 1, size(orders)
         n = orders(io)
         call lebedev_grid(n, p, w, err)
         if (err%has_error()) then
            call report(.false., "lebedev order builds: "//err%get_message())
            cycle
         end if
         worst_sum = max(worst_sum, abs(sum(w) - 1.0_dp))
         do k = 1, n
            worst_sphere = max(worst_sphere, abs(norm2(p(:, k)) - 1.0_dp))
         end do
         ! Averages over the sphere: x, xy, x^2 - y^2, xyz vanish; x^2+y^2+z^2 = 1.
         worst_moment = max(worst_moment, abs(sum(w*p(1, :))), abs(sum(w*p(1, :)*p(2, :))), &
                            abs(sum(w*(p(1, :)**2 - p(2, :)**2))), abs(sum(w*p(1, :)*p(2, :)*p(3, :))), &
                            abs(sum(w*(p(1, :)**2 + p(2, :)**2 + p(3, :)**2)) - 1.0_dp))
         ! Degree 4 too: <x^4> = 1/5 on the unit sphere, exact for every order >= 14.
         if (n >= 14) worst_moment = max(worst_moment, abs(sum(w*p(1, :)**4) - 0.2_dp))
         ! Every point has its antipode with the same weight.
         do k = 1, n
            d = huge(1.0_dp)
            do j = 1, n
               if (abs(w(j) - w(k)) < 1.0e-14_dp) d = min(d, norm2(p(:, j) + p(:, k)))
            end do
            worst_antipode = max(worst_antipode, d)
         end do
         ! Orders 74, 230 and 266 carry negative weights and are declared to.
         negative_ok = negative_ok .and. (any(w < 0.0_dp) .eqv. lebedev_has_negative_weights(n))
      end do
      call report(size(orders) == 32, "lebedev: 32 published orders")
      call report(worst_sum < 1.0e-13_dp, "lebedev: weights sum to one", worst_sum)
      call report(worst_sphere < 1.0e-14_dp, "lebedev: points on the unit sphere", worst_sphere)
      call report(worst_moment < 1.0e-13_dp, "lebedev: low moments integrate exactly", worst_moment)
      call report(worst_antipode < 1.0e-13_dp, "lebedev: closed under inversion", worst_antipode)
      call report(negative_ok, "lebedev: negative-weight orders declared as such")
      call lebedev_grid(100, p, w, err)
      call report(err%has_error(), "lebedev: an order that does not exist is refused")
   end subroutine lebedev_checks

   ! ------------------------------------------------------------------------

   subroutine radial_checks()
      real(dp), allocatable :: r(:), dr(:), w(:)
      type(error_t) :: err
      real(dp) :: total, reach_small, reach_big
      integer :: n

      ! exp(-r^2) over all space is pi^(3/2); at 75 points the mesh contains it.
      call treutler_ahlrichs_radial(75, 8, r, dr, err)
      call report(.not. err%has_error(), "radial: mesh builds")
      call report(all(r > 0.0_dp) .and. all(r(2:) > r(:size(r) - 1)), "radial: nodes ascend and are positive")
      call report(all(dr > 0.0_dp), "radial: weights are positive")
      allocate (w(size(r)))
      call radial_volume_weights(r, dr, w)
      total = sum(w*exp(-r*r))
      call report(abs(total - PI**1.5_dp)/PI**1.5_dp < 1.0e-10_dp, &
                  "radial: integrates exp(-r^2) to pi^(3/2)", abs(total - PI**1.5_dp)/PI**1.5_dp)
      ! 4*pi*r^2*dr and nothing else.
      call report(all(abs(w - 4.0_dp*PI*r*r*dr) <= 1.0e-15_dp*w), "radial: volume weight is 4 pi r^2 dr")
      reach_small = maxval(r)
      call treutler_ahlrichs_radial(200, 8, r, dr, err)
      reach_big = maxval(r)
      call report(reach_big > reach_small, "radial: more points reach further out")
      n = 0
      call treutler_ahlrichs_radial(n, 8, r, dr, err)
      call report(err%has_error(), "radial: an empty mesh is refused")
   end subroutine radial_checks

   ! ------------------------------------------------------------------------

   subroutine partition_checks()
      real(dp), allocatable :: coords(:, :), pts(:, :), w(:)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err
      integer :: k, ia, np
      real(dp) :: worst, mid(3), w_o, w_h, x, w1(1)

      call water(coords, z)
      ! A cloud of points around the molecule, each partitioned once per owner.
      np = 60
      allocate (pts(3, np), owner(np), w(np))
      do k = 1, np
         x = real(k, dp)/real(np, dp)
         pts(:, k) = [3.0_dp*sin(7.0_dp*x), 3.0_dp*cos(5.0_dp*x) - 0.5_dp, 2.5_dp*sin(3.0_dp*x) + 0.4_dp]
      end do
      ! Sum over owners at every point must be exactly one.
      worst = 0.0_dp
      block
         real(dp) :: total(np)
         total = 0.0_dp
         do ia = 1, 3
            owner = ia
            call becke_partition_weights(pts, coords, z, owner, PARTITION_BECKE, ADJUST_TREUTLER, w, err)
            if (err%has_error()) exit
            total = total + w
         end do
         worst = maxval(abs(total - 1.0_dp))
      end block
      call report(.not. err%has_error() .and. worst < 1.0e-13_dp, "partition: cell weights sum to one", worst)

      ! Size adjustment: at the O-H midpoint oxygen takes more than half; without it, exactly half.
      mid = 0.5_dp*(coords(:, 1) + coords(:, 2))
      call becke_partition_weights(reshape(mid, [3, 1]), coords(:, 1:2), z(1:2), [1], PARTITION_BECKE, ADJUST_NONE, w1, err)
      w_o = w1(1)
      call report(abs(w_o - 0.5_dp) < 1.0e-13_dp, "partition: unadjusted midpoint splits evenly", abs(w_o - 0.5_dp))
      call becke_partition_weights(reshape(mid, [3, 1]), coords(:, 1:2), z(1:2), [1], PARTITION_BECKE, ADJUST_TREUTLER, w1, err)
      w_o = w1(1)
      call becke_partition_weights(reshape(mid, [3, 1]), coords(:, 1:2), z(1:2), [2], PARTITION_BECKE, ADJUST_TREUTLER, w1, err)
      w_h = w1(1)
      call report(w_o > 0.5_dp .and. abs(w_o + w_h - 1.0_dp) < 1.0e-13_dp, "partition: size adjustment favours oxygen", w_o)

      ! Stratmann reaches exactly one beside a nucleus.
      call becke_partition_weights(reshape(coords(:, 1) + [0.05_dp, 0.0_dp, 0.0_dp], [3, 1]), coords, z, [1], &
                                   PARTITION_STRATMANN, ADJUST_NONE, w1, err)
      call report(w1(1) == 1.0_dp, "partition: stratmann saturates at exactly one")

      call becke_partition_weights(pts, coords, z(1:2), owner, PARTITION_BECKE, ADJUST_TREUTLER, w, err)
      call report(err%has_error(), "partition: mismatched atom arrays are refused")
   end subroutine partition_checks

   ! ------------------------------------------------------------------------

   subroutine grid_checks()
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      integer :: expected, ia, level
      real(dp) :: g, s, prev, now
      logical :: monotone

      call water(coords, z)

      ! Unpruned, the count is exactly sum over atoms of n_radial * n_angular.
      call build_dft_grid(coords, z, grid, err, level=1, prune=PRUNE_NONE)
      expected = 0
      do ia = 1, size(z)
         expected = expected + grid_level_radial(z(ia), 1)*grid_level_angular(z(ia), 1)
      end do
      call report(.not. err%has_error() .and. grid%n_points == expected .and. &
                  all(grid%atom >= 1 .and. grid%atom <= 3), "grid: unpruned count matches the tables")
      call grid%destroy()

      ! Level 3, the production default, on smooth and on cusped sums.
      call build_dft_grid(coords, z, grid, err, level=3)
      g = abs(integrate_gaussians(grid, coords, 1.3_dp) - 3.0_dp)/3.0_dp
      s = abs(integrate_slaters(grid, coords, 1.7_dp) - 3.0_dp)/3.0_dp
      call report(g < 1.0e-8_dp, "grid: level 3 integrates a Gaussian sum to the atom count", g)
      call report(s < 1.0e-8_dp, "grid: level 3 integrates a Slater sum to the atom count", s)
      call grid%destroy()

      ! Refining must improve it, unpruned.
      prev = huge(1.0_dp)
      monotone = .true.
      do level = 1, 4
         call build_dft_grid(coords, z, grid, err, level=level, prune=PRUNE_NONE)
         now = abs(integrate_gaussians(grid, coords, 1.3_dp) - 3.0_dp)
         monotone = monotone .and. now < prev
         prev = now
         call grid%destroy()
      end do
      call report(monotone, "grid: refining the level reduces the error", prev)

      ! An explicit size applies to every atom; half of one is refused.
      call build_dft_grid(coords, z, grid, err, n_radial=20, n_angular=50, prune=PRUNE_NONE)
      call report(.not. err%has_error() .and. grid%n_points == 3*20*50, "grid: an override applies to every atom")
      call grid%destroy()
      call build_dft_grid(coords, z, grid, err, n_radial=20)
      call report(err%has_error(), "grid: a half-given override is refused")
      call build_dft_grid(coords(:, 1:2), [8], grid, err)
      call report(err%has_error(), "grid: mismatched atom arrays are refused")
   end subroutine grid_checks

   !
   ! The partition's neighbour screening against the unscreened loop, on a
   ! molecule big enough for it to matter: gly10, 73 atoms. Both the weights
   ! themselves and an integral over them. The unscreened build is the
   ! O(n_atoms^2)-per-point loop this screening replaces, and the numbers
   ! here are what the default threshold costs.
   !
   subroutine screening_checks()
      integer :: nat, k, n
      integer, allocatable :: z(:)
      real(dp), allocatable :: coords(:, :), w_full(:), p_full(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp) :: dmax, rmax, t0, t1, t2, g_scr, g_full, wscale
      integer(kind=8) :: cc, rate

      call read_xyz('gly10.xyz', nat, z, coords)
      call system_clock(cc, rate); t0 = real(cc, dp)/real(rate, dp)
      call build_dft_grid(coords, z, grid, err, level=3)
      call system_clock(cc, rate); t1 = real(cc, dp)/real(rate, dp)
      if (err%has_error()) then
         call report(.false., "screening: gly10 grid builds")
         return
      end if
      n = grid%n_points
      allocate (w_full(n), p_full(n))
      call becke_partition_weights(grid%coords, coords, z, grid%atom, grid%scheme, grid%adjust, &
                                   p_full, err, nu_max=2.0_dp)
      call system_clock(cc, rate); t2 = real(cc, dp)/real(rate, dp)
      w_full = grid%quad_weights*p_full
      dmax = maxval(abs(grid%weights - w_full))
      wscale = maxval(abs(w_full))
      rmax = 0.0_dp
      do k = 1, n
         if (abs(w_full(k)) > 1.0e-10_dp*wscale) rmax = max(rmax, abs(grid%weights(k) - w_full(k))/abs(w_full(k)))
      end do
      g_scr = integrate_gaussians(grid, coords, 1.3_dp)
      g_full = 0.0_dp
      block
         type(dft_grid_t) :: g2
         g2 = grid
         g2%weights = w_full
         g_full = integrate_gaussians(g2, coords, 1.3_dp)
      end block
      print '(a,i0,a,i0,a,f7.2,a,f7.2,a)', "screening: gly10, ", nat, " atoms, ", n, " points; grid build ", &
         t1 - t0, " s screened, partition alone ", t2 - t1, " s unscreened"
      call report(dmax < 1.0e-8_dp*wscale, "screening: max |dw| / max w", dmax/wscale)
      call report(rmax < 1.0e-6_dp, "screening: max relative dw on weights above 1e-10", rmax)
      call report(abs(g_scr - g_full) < 1.0e-9_dp*real(nat, dp), "screening: Gaussian sum, screened vs full", &
                  abs(g_scr - g_full))
      call report(abs(g_scr - real(nat, dp))/real(nat, dp) < 1.0e-6_dp, "screening: Gaussian sum vs atom count", &
                  abs(g_scr - real(nat, dp))/real(nat, dp))
      call grid%destroy()
   end subroutine screening_checks

   function integrate_gaussians(grid, coords, alpha) result(total)
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: coords(:, :), alpha
      real(dp) :: total, f, d2, norm
      integer :: k, ia
      norm = (alpha/PI)**1.5_dp
      total = 0.0_dp
      do k = 1, grid%n_points
         f = 0.0_dp
         do ia = 1, size(coords, 2)
            d2 = sum((grid%coords(:, k) - coords(:, ia))**2)
            f = f + norm*exp(-alpha*d2)
         end do
         total = total + grid%weights(k)*f
      end do
   end function integrate_gaussians

   function integrate_slaters(grid, coords, zeta) result(total)
      type(dft_grid_t), intent(in) :: grid
      real(dp), intent(in) :: coords(:, :), zeta
      real(dp) :: total, f, d, norm
      integer :: k, ia
      norm = zeta**3/PI
      total = 0.0_dp
      do k = 1, grid%n_points
         f = 0.0_dp
         do ia = 1, size(coords, 2)
            d = norm2(grid%coords(:, k) - coords(:, ia))
            f = f + norm*exp(-2.0_dp*zeta*d)
         end do
         total = total + grid%weights(k)*f
      end do
   end function integrate_slaters

end program check_grid
