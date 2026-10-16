!
! Negative test: trc_df_3c refuses an orbital basis above DF_LMAX.
!
! Same reasoning as check_df_2c_ceiling.F90's own header, one level up:
! trc_df_3c's own guard (trc_api.F90) exists because gen_df.py's VRR+HRR
! and Rys bodies are generated only up to DF_LMAX on the orbital side, and
! this is the control that the guard still fires rather than reading an
! unwritten block.
!
! Orbital basis at l = DF_LMAX + 1 -- one above THIS build's ceiling, valid
! at any ceiling -- paired with an auxiliary basis at l = 0, always within
! DF_LMAX_AUX, so only the orbital check is under test here.
!
! Registered with PASS_REGULAR_EXPRESSION matching the error-stop message
! itself (test/CMakeLists.txt), not WILL_FAIL: a different crash must not
! count as this test passing.
!
program check_df_3c_ceiling
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_df_3c
   use trc_df_kernels, only: DF_LMAX
   use trc_test_basis, only: build_aux
   implicit none

   integer, parameter :: nat = 1
   integer  :: at_z(nat), nsh, maxnp, naux, maxnp_aux
   real(dp) :: at_r(3, nat)
   integer,  allocatable :: sh_l(:), sh_np(:), sh_l_aux(:), sh_np_aux(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: sh_e_aux(:, :), sh_c_aux(:, :), sh_r_aux(:, :)
   real(dp), allocatable :: zat(:), tens(:, :, :)
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl

   at_z = [2]
   at_r = reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1])
   allocate (zat(nat)); zat = real(at_z, dp)

   call build_aux(nat, at_z, at_r, DF_LMAX + 1, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call pl%build(bas, 1.0e-14_dp)

   call build_aux(nat, at_z, at_r, 0, naux, sh_l_aux, sh_np_aux, sh_e_aux, sh_c_aux, sh_r_aux, maxnp_aux)
   call aux%build(naux, sh_l_aux, sh_np_aux, sh_e_aux, sh_c_aux, sh_r_aux, nat, zat, at_r, maxnp_aux)

   print '(a,i0,a,i0)', '  built an orbital basis with l = ', maxval(sh_l), &
      ', this build has DF_LMAX = ', DF_LMAX

   allocate (tens(bas%nao, bas%nao, aux%nao))
   call trc_df_3c(bas, pl, aux, tens)   ! must error stop before the next line runs
   print '(a)', '  RESULT: FAIL (trc_df_3c did not refuse the basis)'
end program check_df_3c_ceiling
