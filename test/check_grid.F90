module check_grid_suite
   !! The molecular quadrature: closed forms and identities.
   !!
   !! The grid modules came from metalquicha with a test-drive suite of fifty
   !! cases. This is the subset of those cases that a wrong grid cannot pass:
   !!
   !!   * a Lebedev order integrates low polynomials over the sphere exactly,
   !!     its weights sum to one, and every point has its antipode -- the
   !!     quadrature property, which a point set that merely lies on a sphere
   !!     would fail;
   !!   * the Becke cell weights sum to one at every point, which fails for an
   !!     asymmetric pair loop, a size adjustment with the wrong sign, or an
   !!     atom dropped from a product;
   !!   * the radial mesh integrates exp(-r^2) exactly, probing the Jacobian;
   !!   * the assembled grid integrates a sum of normalised Gaussians and of
   !!     normalised Slaters over water to the atom count, and refining the
   !!     level reduces the error;
   !!   * bad requests are refused rather than half-applied.
   !!
   !! Every tolerance is one the metalquicha suite measured, not one chosen to
   !! pass: level 3 reaches 5.7e-9 on the Gaussian sum and the bound is 1e-8.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_lebedev, only: lebedev_grid, lebedev_available_orders, lebedev_has_negative_weights
   use trc_dft_radial, only: dft_radial_treutler_ahlrichs, dft_radial_volume_weights
   use trc_dft_partition, only: dft_partition_becke_weights, dft_partition_options_t, PARTITION_STRATMANN, &
                                ADJUST_NONE
   use trc_dft_prune, only: PRUNE_NONE
   use trc_dft_grid, only: dft_grid_t, dft_grid_options_t, dft_grid_build, dft_grid_level_radial, &
                           dft_grid_level_angular
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_grid_tests

   real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp

contains

   subroutine collect_grid_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Lebedev orders: quadrature, sphere, antipode", test_lebedev), &
                  new_unittest("radial mesh integrates exp(-r^2)", test_radial), &
                  new_unittest("Becke/Stratmann cell weights", test_partition), &
                  new_unittest("assembled grid integrates to the atom count", test_grid), &
                  new_unittest("neighbour screening vs the unscreened partition (gly10)", test_screening) &
                  ]
   end subroutine collect_grid_tests

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

   subroutine test_lebedev(error)
      type(error_type), allocatable, intent(out) :: error
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
            call check(error, .false., "lebedev order builds: "//err%get_message())
            return
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

      call check(error, size(orders) == 32, "lebedev: 32 published orders, got "//to_string(size(orders)))
      if (allocated(error)) return
      call check(error, worst_sum < 1.0e-13_dp, "lebedev: weights sum to one, worst = "// &
                 to_string(worst_sum, "es10.2"))
      if (allocated(error)) return
      call check(error, worst_sphere < 1.0e-14_dp, "lebedev: points on the unit sphere, worst = "// &
                 to_string(worst_sphere, "es10.2"))
      if (allocated(error)) return
      call check(error, worst_moment < 1.0e-13_dp, "lebedev: low moments integrate exactly, worst = "// &
                 to_string(worst_moment, "es10.2"))
      if (allocated(error)) return
      call check(error, worst_antipode < 1.0e-13_dp, "lebedev: closed under inversion, worst = "// &
                 to_string(worst_antipode, "es10.2"))
      if (allocated(error)) return
      call check(error, negative_ok, "lebedev: negative-weight orders (74, 230, 266) declared as such")
      if (allocated(error)) return

      call lebedev_grid(100, p, w, err)
      call check(error, err%has_error(), "lebedev: an order that does not exist is refused")
   end subroutine test_lebedev

   ! ------------------------------------------------------------------------

   subroutine test_radial(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:), dr(:), w(:)
      type(error_t) :: err
      real(dp) :: total, reach_small, reach_big, rel
      integer :: n

      ! exp(-r^2) over all space is pi^(3/2); at 75 points the mesh contains it.
      call dft_radial_treutler_ahlrichs(75, 8, r, dr, err)
      call check(error, .not. err%has_error(), "radial: mesh builds")
      if (allocated(error)) return
      call check(error, all(r > 0.0_dp) .and. all(r(2:) > r(:size(r) - 1)), &
                 "radial: nodes ascend and are positive")
      if (allocated(error)) return
      call check(error, all(dr > 0.0_dp), "radial: weights are positive")
      if (allocated(error)) return

      allocate (w(size(r)))
      call dft_radial_volume_weights(r, dr, w)
      total = sum(w*exp(-r*r))
      rel = abs(total - PI**1.5_dp)/PI**1.5_dp
      call check(error, rel < 1.0e-10_dp, "radial: integrates exp(-r^2) to pi^(3/2), rel = "// &
                 to_string(rel, "es10.2"))
      if (allocated(error)) return
      ! 4*pi*r^2*dr and nothing else.
      call check(error, all(abs(w - 4.0_dp*PI*r*r*dr) <= 1.0e-15_dp*w), &
                 "radial: volume weight is 4 pi r^2 dr")
      if (allocated(error)) return

      reach_small = maxval(r)
      call dft_radial_treutler_ahlrichs(200, 8, r, dr, err)
      reach_big = maxval(r)
      call check(error, reach_big > reach_small, "radial: more points reach further out")
      if (allocated(error)) return

      n = 0
      call dft_radial_treutler_ahlrichs(n, 8, r, dr, err)
      call check(error, err%has_error(), "radial: an empty mesh is refused")
   end subroutine test_radial

   ! ------------------------------------------------------------------------

   subroutine test_partition(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :), pts(:, :), w(:)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err
      integer :: k, ia, np
      real(dp) :: worst, mid(3), w_o, w_h, x, w1(1)
      type(dft_partition_options_t) :: becke, becke_unadjusted, stratmann_unadjusted

      becke_unadjusted%adjust = ADJUST_NONE
      stratmann_unadjusted%scheme = PARTITION_STRATMANN
      stratmann_unadjusted%adjust = ADJUST_NONE
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
            call dft_partition_becke_weights(pts, coords, z, owner, becke, w, err)
            if (err%has_error()) exit
            total = total + w
         end do
         worst = maxval(abs(total - 1.0_dp))
      end block
      call check(error, .not. err%has_error() .and. worst < 1.0e-13_dp, &
                 "partition: cell weights sum to one, worst = "//to_string(worst, "es10.2"))
      if (allocated(error)) return

      ! Size adjustment: at the O-H midpoint oxygen takes more than half; without it, exactly half.
      mid = 0.5_dp*(coords(:, 1) + coords(:, 2))
      call dft_partition_becke_weights(reshape(mid, [3, 1]), coords(:, 1:2), z(1:2), [1], becke_unadjusted, w1, err)
      w_o = w1(1)
      call check(error, abs(w_o - 0.5_dp) < 1.0e-13_dp, "partition: unadjusted midpoint splits evenly, |dw| = "// &
                 to_string(abs(w_o - 0.5_dp), "es10.2"))
      if (allocated(error)) return

      call dft_partition_becke_weights(reshape(mid, [3, 1]), coords(:, 1:2), z(1:2), [1], becke, w1, err)
      w_o = w1(1)
      call dft_partition_becke_weights(reshape(mid, [3, 1]), coords(:, 1:2), z(1:2), [2], becke, w1, err)
      w_h = w1(1)
      call check(error, w_o > 0.5_dp .and. abs(w_o + w_h - 1.0_dp) < 1.0e-13_dp, &
                 "partition: size adjustment favours oxygen, w_o = "//to_string(w_o, "f10.6"))
      if (allocated(error)) return

      ! Stratmann reaches exactly one beside a nucleus.
      call dft_partition_becke_weights(reshape(coords(:, 1) + [0.05_dp, 0.0_dp, 0.0_dp], [3, 1]), coords, z, [1], &
                                       stratmann_unadjusted, w1, err)
      call check(error, w1(1) == 1.0_dp, "partition: stratmann saturates at exactly one")
      if (allocated(error)) return

      call dft_partition_becke_weights(pts, coords, z(1:2), owner, becke, w, err)
      call check(error, err%has_error(), "partition: mismatched atom arrays are refused")
   end subroutine test_partition

   ! ------------------------------------------------------------------------

   subroutine test_grid(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: z(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      integer :: expected, ia, level
      real(dp) :: g, s, prev, now
      logical :: monotone
      type(dft_grid_options_t) :: unpruned, override, half_override

      unpruned%prune = PRUNE_NONE
      override%prune = PRUNE_NONE
      override%n_radial = 20
      override%n_angular = 50
      half_override%n_radial = 20
      call water(coords, z)

      ! Unpruned, the count is exactly sum over atoms of n_radial * n_angular.
      call dft_grid_build(coords, z, grid, err, level=1, options=unpruned)
      expected = 0
      do ia = 1, size(z)
         expected = expected + dft_grid_level_radial(z(ia), 1)*dft_grid_level_angular(z(ia), 1)
      end do
      call check(error, .not. err%has_error() .and. grid%n_points == expected .and. &
                 all(grid%atom >= 1 .and. grid%atom <= 3), "grid: unpruned count matches the tables")
      call grid%destroy()
      if (allocated(error)) return

      ! Level 3, the production default, on smooth and on cusped sums.
      call dft_grid_build(coords, z, grid, err, level=3)
      g = abs(integrate_gaussians(grid, coords, 1.3_dp) - 3.0_dp)/3.0_dp
      s = abs(integrate_slaters(grid, coords, 1.7_dp) - 3.0_dp)/3.0_dp
      call check(error, g < 1.0e-8_dp, "grid: level 3 integrates a Gaussian sum to the atom count, rel = "// &
                 to_string(g, "es10.2"))
      if (.not. allocated(error)) then
         call check(error, s < 1.0e-8_dp, "grid: level 3 integrates a Slater sum to the atom count, rel = "// &
                    to_string(s, "es10.2"))
      end if
      call grid%destroy()
      if (allocated(error)) return

      ! Refining must improve it, unpruned.
      prev = huge(1.0_dp)
      monotone = .true.
      do level = 1, 4
         call dft_grid_build(coords, z, grid, err, level=level, options=unpruned)
         now = abs(integrate_gaussians(grid, coords, 1.3_dp) - 3.0_dp)
         monotone = monotone .and. now < prev
         prev = now
         call grid%destroy()
      end do
      call check(error, monotone, "grid: refining the level reduces the error")
      if (allocated(error)) return

      ! An explicit size applies to every atom; half of one is refused.
      call dft_grid_build(coords, z, grid, err, options=override)
      call check(error, .not. err%has_error() .and. grid%n_points == 3*20*50, &
                 "grid: an override applies to every atom")
      call grid%destroy()
      if (allocated(error)) return

      call dft_grid_build(coords, z, grid, err, options=half_override)
      call check(error, err%has_error(), "grid: a half-given override is refused")
      if (allocated(error)) return

      call dft_grid_build(coords(:, 1:2), [8], grid, err)
      call check(error, err%has_error(), "grid: mismatched atom arrays are refused")
   end subroutine test_grid

   !
   ! The partition's neighbour screening against the unscreened loop, on a
   ! molecule big enough for it to matter: gly10, 73 atoms. Both the weights
   ! themselves and an integral over them. The unscreened build is the
   ! O(n_atoms^2)-per-point loop this screening replaces, and the numbers
   ! here are what the default threshold costs.
   !
   subroutine test_screening(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nat, k, n
      integer, allocatable :: z(:)
      real(dp), allocatable :: coords(:, :), w_full(:), p_full(:)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp) :: dmax, rmax, t1, t2, g_scr, g_full, wscale
      type(timer_type) :: tm
      type(dft_partition_options_t) :: unscreened

      call read_xyz('gly10.xyz', nat, z, coords)
      call tm%start()
      call dft_grid_build(coords, z, grid, err, level=3)
      t1 = tm%get_elapsed_time()
      call check(error, .not. err%has_error(), "screening: gly10 grid builds")
      if (allocated(error)) return

      call tm%start()
      n = grid%n_points
      allocate (w_full(n), p_full(n))
      unscreened = grid%partition
      unscreened%nu_max = 2.0_dp
      call dft_partition_becke_weights(grid%coords, coords, z, grid%atom, unscreened, p_full, err)
      t2 = tm%get_elapsed_time()
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
         t1, " s screened, partition alone ", t2, " s unscreened"

      call check(error, dmax < 1.0e-8_dp*wscale, "screening: max |dw| / max w = "// &
                 to_string(dmax/wscale, "es10.2"))
      if (.not. allocated(error)) then
         call check(error, rmax < 1.0e-6_dp, "screening: max relative dw on weights above 1e-10 = "// &
                    to_string(rmax, "es10.2"))
      end if
      if (.not. allocated(error)) then
         call check(error, abs(g_scr - g_full) < 1.0e-9_dp*real(nat, dp), &
                    "screening: Gaussian sum, screened vs full, |diff| = "//to_string(abs(g_scr - g_full), "es10.2"))
      end if
      if (.not. allocated(error)) then
         call check(error, abs(g_scr - real(nat, dp))/real(nat, dp) < 1.0e-6_dp, &
                    "screening: Gaussian sum vs atom count, rel = "// &
                    to_string(abs(g_scr - real(nat, dp))/real(nat, dp), "es10.2"))
      end if
      call grid%destroy()
   end subroutine test_screening

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

end module check_grid_suite

program check_grid
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_grid_suite, only: collect_grid_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat

   stat = 0
   testsuites = [new_testsuite("check_grid", collect_grid_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_grid
