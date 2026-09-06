!
! The multipole driver and its per-item routine, as an INCLUDE for the
! module that holds the kernels they dispatch to. They used to live in
! trc_api and call the dispatcher across the module boundary, which is a
! real device call per item -- stdpar only inlines within a translation
! unit. Included after `contains` in the kernel module, they inline; the
! guard keeps a unit from getting two copies. NC_MAX and the NBLK sizes
! come from the including module.
!
#ifndef TRC_MULT_DRIVER_INC
#define TRC_MULT_DRIVER_INC
   subroutine mult_driver(nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                          sh_e, sh_c, sh_r, npair, p_i, p_j, orig, mmat)
      integer,  intent(in)  :: nshell, nao, maxnp, npair
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell), orig(3)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      real(dp), intent(out) :: mmat(nao, nao, TRC_NMULT)

      integer :: ij, i, j, k

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(k=1:TRC_NMULT, j=1:nao, i=1:nao)
         mmat(i, j, k) = 0.0_dp
      end do

      do concurrent(ij = 1:npair)
         call mult_item(p_i(ij), p_j(ij), nshell, nao, maxnp, sh_l, sh_np, &
                        sh_ao, sh_e, sh_c, sh_r, orig, mmat)
      end do
   end subroutine mult_driver

   pure subroutine mult_item(ii, jj, nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                             sh_e, sh_c, sh_r, orig, mmat)
      !$acc routine seq
      integer,  intent(in)    :: ii, jj, nshell, nao, maxnp
      integer,  intent(in)    :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell), orig(3)
      real(dp), intent(inout) :: mmat(nao, nao, TRC_NMULT)

      integer  :: li, lj, ni, nj, a, c, mu, nu, k, ic
      !
      ! ONE DIMENSION, INDEXED BY HAND. The generated kernel declares its
      ! output `mout(ni*nj, 39)` with the ACTUAL block size of its class,
      ! and reaches this routine through a `mout(*)` dummy -- so it writes
      ! component ic starting at (ic-1)*ni*nj. Declaring the buffer here as
      ! `mblk(NBLK_MAX, TRC_NMULT)` reads it back with stride NBLK_MAX
      ! instead, which agrees only for the one class whose block happens to
      ! be NBLK_MAX and silently scrambles every other one. That bug cost a
      ! day: the dipole x component is the first ni*nj elements and so came
      ! out right, which made it look like an ordering problem in the
      ! generator rather than a stride mismatch here.
      !
      real(dp) :: mblk(NBLK_MAX*TRC_NMULT)

      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      call multipole_dispatch(li, lj, sh_np(ii), sh_np(jj), &
                              sh_e(:, ii), sh_c(:, ii), &
                              sh_e(:, jj), sh_c(:, jj), &
                              sh_r(:, ii), sh_r(:, jj), orig, mblk)
      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            do ic = 1, TRC_NMULT
               ! the multipole matrices are symmetric: the operator is
               ! multiplicative and carries no derivative
               mmat(mu, nu, ic) = mblk(k + (ic - 1)*ni*nj)
               mmat(nu, mu, ic) = mblk(k + (ic - 1)*ni*nj)
            end do
         end do
      end do
   end subroutine mult_item
#endif
