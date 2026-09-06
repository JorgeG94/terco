!! The primitive-shell view against the segmented path: one Fock build each
!! on water/cc-pVDZ, whose s shells carry two general columns over nine
!! shared primitives and whose uncontracted columns come out of `build` as
!! one-primitive shells. The two must agree to the accumulation order.
program check_gc
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   integer :: natm, n, i, j, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), d(:, :), g1(:, :), g2(:, :)
   real(dp) :: worst
   type(trc_basis_t) :: b
   type(trc_eri_t) :: e1, e2
   type(error_t) :: err

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, b, err)
   if (err%has_error()) error stop 'check_gc: basis'
   call b%to_device()
   n = b%nao
   allocate (d(n, n), g1(n, n), g2(n, n))
   do j = 1, n
      do i = 1, n
         d(i, j) = 0.05_dp*exp(-0.3_dp*abs(i - j))
      end do
   end do
   !$acc enter data copyin(d) create(g1, g2)
   call e1%build(b, 1.0e-12_dp)
   call e1%fock_resident(b, d, g1, k_scale=1.0_dp)
   call e2%build(b, 1.0e-12_dp, general=.true.)
   call e2%fock_resident(b, d, g2, k_scale=1.0_dp)
   !$acc update self(g1, g2)
   worst = maxval(abs(g1 - g2))
   print '(a,i0,a,i0,a,i0)', '  nao ', n, '  segmented launches ', e1%nlaunch, '  general launches ', e2%nlaunch
   print '(a,es10.2,a,es10.2)', '  worst |G_seg - G_gen| ', worst, '  scale ', maxval(abs(g1))
   if (worst > 1.0e-10_dp*max(1.0_dp, maxval(abs(g1)))) then
      print '(a)', '  RESULT: FAIL'
      error stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_gc
