!
! Negative test: trc_df_2c refuses an auxiliary basis above DF_LMAX_AUX.
!
! Nothing else in this suite calls trc_df_2c with a basis past this build's
! own ceiling -- every other check builds within it on purpose. The guard
! trc_api.F90 carries exists because the DF dispatchers are generated only
! up to DF_LMAX_AUX (gen_df.py --lmax-aux): a class above that has no
! written-out block, so without the guard trc_df_2c would read one
! silently rather than error stop. This is the control that the guard
! itself still fires.
!
! One shell, l = DF_LMAX_AUX + 1 -- one above THIS build's ceiling, so the
! test is valid whatever DF_LMAX_AUX happens to be -- on a single atom; the
! content of the basis does not matter, only that trc_df_2c sees an l it
! has no kernel for.
!
! Registered with PASS_REGULAR_EXPRESSION matching the error-stop message
! itself (test/CMakeLists.txt), not WILL_FAIL: a different crash, e.g. a
! segfault before the guard runs, must not count as this test passing.
!
program check_df_2c_ceiling
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_df_2c
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_test_basis, only: build_aux
   implicit none

   integer, parameter :: nat = 1
   integer  :: at_z(nat), nsh, maxnp
   real(dp) :: at_r(3, nat)
   integer,  allocatable :: sh_l(:), sh_np(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), zat(:)
   real(dp), allocatable :: jmat(:, :)
   type(trc_basis_t) :: aux

   at_z = [2]
   at_r = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
   allocate (zat(nat)); zat = real(at_z, dp)

   call build_aux(nat, at_z, at_r, DF_LMAX_AUX + 1, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   call aux%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)

   print '(a,i0,a,i0)', '  built an auxiliary basis with l = ', maxval(sh_l), &
      ', this build has DF_LMAX_AUX = ', DF_LMAX_AUX

   allocate (jmat(aux%nao, aux%nao))
   call trc_df_2c(aux, jmat)   ! must error stop before the next line runs
   print '(a)', '  RESULT: FAIL (trc_df_2c did not refuse the basis)'
end program check_df_2c_ceiling
