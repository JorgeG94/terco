!
! The (mu nu|P) driver and its per-item routine, as an INCLUDE for the
! module that holds the kernels they dispatch to. They used to live in
! trc_api and call the dispatcher across the module boundary, which is a
! real device call per item -- stdpar only inlines within a translation
! unit. Included after `contains` in the kernel module, they inline; the
! guard keeps a unit from getting two copies. NC_MAX and the NBLK sizes
! come from the including module.
!
#ifndef TRC_DF_3C_DRIVER_INC
#define TRC_DF_3C_DRIVER_INC
   subroutine df_3c_driver(nsh, nao, mnp, sh_l, sh_ao, sh_r, &
                           npair, npp, p_i, p_j, pp_off, pp_n, &
                           pp_z, pp_k, pp_p, pp_a, &
                           nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                           tens)
      integer,  intent(in)  :: nsh, nao, mnp, npair, npp, nsx, nax, mnx
      integer,  intent(in)  :: sh_l(nsh), sh_ao(nsh)
      real(dp), intent(in)  :: sh_r(3, nsh)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)  :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)  :: ax_l(nsx), ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)  :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      real(dp), intent(out) :: tens(nao, nao, nax)

      integer(kind=8) :: t
      integer :: i, j, k

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(k=1:nax, j=1:nao, i=1:nao)
         tens(i, j, k) = 0.0_dp
      end do

      ! One thread per (SURVIVING shell pair, auxiliary shell). The pair index
      ! runs fastest so that neighbouring threads share the auxiliary shell,
      ! and therefore its primitive data, across a warp.
      do concurrent(t = 0:int(npair, 8)*nsx - 1)
         call df_3c_item(t, nsh, nao, sh_l, sh_ao, sh_r, npair, npp, &
                         p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                         nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                         tens)
      end do
   end subroutine df_3c_driver

   pure subroutine df_3c_item(t, nsh, nao, sh_l, sh_ao, sh_r, npair, npp, &
                              p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                              nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, &
                              ax_r, tens)
      !$acc routine seq
      integer(kind=8), intent(in) :: t
      integer,  intent(in)    :: nsh, nao, npair, npp, nsx, nax, mnx
      integer,  intent(in)    :: sh_l(nsh), sh_ao(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: p_i(npair), p_j(npair)
      integer,  intent(in)    :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp)
      real(dp), intent(in)    :: pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)    :: ax_l(nsx), ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      real(dp), intent(inout) :: tens(nao, nao, nax)

      integer :: ip, kk, ii, jj, li, lj, lk, ni, nj, nk
      integer :: a, c, e, mu, nu, pq, m, o
      real(dp) :: blk(NBLK3_MAX), abx, aby, abz

      ip = int(mod(t, int(npair, 8))) + 1
      kk = int(t/npair) + 1
      ii = p_i(ip); jj = p_j(ip)
      li = sh_l(ii); lj = sh_l(jj); lk = ax_l(kk)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      nk = (lk + 1)*(lk + 2)/2
      abx = sh_r(1, ii) - sh_r(1, jj)
      aby = sh_r(2, ii) - sh_r(2, jj)
      abz = sh_r(3, ii) - sh_r(3, jj)
      o = pp_off(ip)

      call df_3c_dispatch(li, lj, lk, pp_n(ip), &
                          pp_z(o + 1:o + pp_n(ip)), pp_k(o + 1:o + pp_n(ip)), &
                          pp_p(:, o + 1:o + pp_n(ip)), &
                          pp_a(:, o + 1:o + pp_n(ip)), &
                          abx, aby, abz, &
                          ax_np(kk), ax_e(:, kk), ax_c(:, kk), ax_r(:, kk), blk)

      ! The pair list is lower-triangular in (ii, jj), and (mu nu|P) is
      ! symmetric in mu <-> nu, so each thread fills both. On a diagonal
      ! shell pair the second write repeats the first with the same value.
      m = 0
      do e = 0, nk - 1
         pq = ax_ao(kk) + e
         do c = 0, nj - 1
            nu = sh_ao(jj) + c
            do a = 0, ni - 1
               mu = sh_ao(ii) + a
               m = m + 1
               tens(mu, nu, pq) = blk(m)
               tens(nu, mu, pq) = blk(m)
            end do
         end do
      end do
   end subroutine df_3c_item
#endif
