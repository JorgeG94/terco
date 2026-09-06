!
! The (P|Q) driver and its per-item routine, as an INCLUDE for the
! module that holds the kernels they dispatch to. They used to live in
! trc_api and call the dispatcher across the module boundary, which is a
! real device call per item -- stdpar only inlines within a translation
! unit. Included after `contains` in the kernel module, they inline; the
! guard keeps a unit from getting two copies. NC_MAX and the NBLK sizes
! come from the including module.
!
#ifndef TRC_DF_2C_DRIVER_INC
#define TRC_DF_2C_DRIVER_INC
   subroutine df_2c_driver(nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                           sh_e, sh_c, sh_r, jmat)
      integer,  intent(in)  :: nshell, nao, maxnp
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(out) :: jmat(nao, nao)

      integer :: ij, i, j

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(j=1:nao, i=1:nao)
         jmat(i, j) = 0.0_dp
      end do

      do concurrent(ij = 0:nshell*nshell - 1)
         call df_2c_item(ij, nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                         sh_e, sh_c, sh_r, jmat)
      end do
   end subroutine df_2c_driver

   pure subroutine df_2c_item(ij, nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                              sh_e, sh_c, sh_r, jmat)
      !$acc routine seq
      integer,  intent(in)    :: ij, nshell, nao, maxnp
      integer,  intent(in)    :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      real(dp), intent(inout) :: jmat(nao, nao)

      integer  :: ii, jj, li, lj, ni, nj, a, c, mu, nu, k
      real(dp) :: blk(NBLK_MAX)

      ii = ij/nshell + 1
      jj = mod(ij, nshell) + 1
      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      call df_2c_dispatch(li, lj, sh_np(ii), sh_e(:, ii), sh_c(:, ii), &
                          sh_r(:, ii), sh_np(jj), sh_e(:, jj), &
                          sh_c(:, jj), sh_r(:, jj), blk)
      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            jmat(mu, nu) = blk(k)
         end do
      end do
   end subroutine df_2c_item
#endif
