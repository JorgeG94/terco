program check_one_e_grad_ceiling
   !!
   !! Negative test: trc_1e_gradient refuses a basis above ONE_E_GRAD_LMAX.
   !!
   !! The gradient twin of check_one_e_ceiling.F90: one_e_grad_dispatch has no
   !! case above ONE_E_GRAD_LMAX either, and trc_1e_gradient had no guard.
   !!
   !! A basis at l = ONE_E_GRAD_LMAX + 1, one above THIS build's ceiling, so the test
   !! holds at any ceiling. Registered with PASS_REGULAR_EXPRESSION on the
   !! error-stop message (test/CMakeLists.txt), not WILL_FAIL: a different crash
   !! must not count as a pass. check_df_3c_ceiling.F90 is the same test for
   !! trc_df_3c.
   !!
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e_gradient
   use trc_1e_grad_kernels, only: ONE_E_GRAD_LMAX
   use trc_test_basis, only: build_aux
   implicit none
   integer, parameter :: nat = 1
   integer  :: at_z(nat), nsh, maxnp
   real(dp) :: at_r(3, nat)
   integer,  allocatable :: sh_l(:), sh_np(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: zat(:), w(:, :), d(:, :), g(:, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl

   at_z = [2]
   at_r = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
   allocate (zat(nat)); zat = real(at_z, dp)
   call build_aux(nat, at_z, at_r, ONE_E_GRAD_LMAX + 1, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call pl%build(bas, 1.0e-14_dp)
   print '(a,i0,a,i0)', '  built a basis with l = ', maxval(sh_l), &
      ', this build has ONE_E_GRAD_LMAX = ', ONE_E_GRAD_LMAX
   allocate (w(bas%nao, bas%nao), d(bas%nao, bas%nao), g(3, nat))
   w = 0.0_dp; d = 0.0_dp; g = 0.0_dp
   call trc_1e_gradient(bas, pl, w, d, g)   ! must error stop before the next line runs
   print '(a)', '  RESULT: FAIL (trc_1e_gradient did not refuse the basis)'
end program check_one_e_grad_ceiling
