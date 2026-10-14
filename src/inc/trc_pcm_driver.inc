!
! The PCM surface-integral drivers and their per-item routines, as an
! INCLUDE for the generated module that holds the kernels (trc_pcm_kernels,
! from scripts/gen_1e.py --pcm). Included there for the reason
! trc_1e_driver.inc is: stdpar inlines only within a translation unit.
! NC_MAX and NBLK_MAX come from the including module.
!
! Both contract the (npair x npts) integral tensor as it is made and never
! store it: a cavity has 10^3 to 10^5 points, and the tensor at that size
! is the whole memory budget of the device several times over.
!
! The primitive-pair data is the pair list's own (trc_pairlist_t's pp_*),
! so the products, P and K_ab are formed once per geometry rather than once
! per (pair, point).
!
#ifndef TRC_PCM_DRIVER_INC
#define TRC_PCM_DRIVER_INC
   subroutine pcm_potential_driver(nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                   pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, &
                                   dmat, v)
      !! `v(k) = sum_{mu nu} D_mu,nu <mu| erf(zeta_k|r - C_k|)/|r - C_k| |nu>`
      integer,  intent(in)  :: nshell, nao, npair, npp, npts
      integer,  intent(in)  :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      integer,  intent(in)  :: p_i(npair), p_j(npair), pp_off(npair), pp_n(npair)
      real(dp), intent(in)  :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      real(dp), intent(in)  :: pts(3, npts)    !! point centres
      real(dp), intent(in)  :: zk(npts)      !! zeta_k; the exponent is zeta_k^2
      real(dp), intent(in)  :: dmat(nao, nao)
      real(dp), intent(out) :: v(npts)

      integer :: k

      ! One thread per POINT, walking every surviving shell pair: the point
      ! owns v(k), so nothing is shared and nothing needs an atomic. A point
      ! per thread is also the wide axis -- thousands of points against a
      ! few hundred pairs on a small molecule.
      do concurrent(k = 1:npts)
         call pcm_potential_item(k, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                 pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, &
                                 dmat, v)
      end do
   end subroutine pcm_potential_driver

   pure subroutine pcm_potential_item(k, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                      pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, &
                                      dmat, v)
      !$acc routine seq
      integer,  intent(in)    :: k, nshell, nao, npair, npp, npts
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: p_i(npair), p_j(npair), pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      real(dp), intent(in)    :: pts(3, npts), zk(npts)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: v(npts)

      integer  :: ij, ii, jj, li, lj, ni, nj, a, c, mu, nu, o, n
      real(dp) :: blk(NBLK_MAX), one(1), acc, abx, aby, abz

      one(1) = 1.0_dp
      acc = 0.0_dp
      do ij = 1, npair
         ii = p_i(ij); jj = p_j(ij)
         li = sh_l(ii); lj = sh_l(jj)
         ni = (li + 1)*(li + 2)/2
         nj = (lj + 1)*(lj + 2)/2
         abx = sh_r(1, ii) - sh_r(1, jj)
         aby = sh_r(2, ii) - sh_r(2, jj)
         abz = sh_r(3, ii) - sh_r(3, jj)
         o = pp_off(ij); n = pp_n(ij)
         call pcm_dispatch(li, lj, n, pp_z(o + 1:o + n), pp_k(o + 1:o + n), &
                           pp_p(:, o + 1:o + n), pp_a(:, o + 1:o + n), abx, aby, abz, &
                           1, pts(:, k:k), zk(k:k), one, blk)
         ! The pair list is lower-triangular, and an off-diagonal block
         ! stands for itself and its transpose; D is taken as given rather
         ! than assumed symmetric, so both of its halves are read.
         ! `blk` addressed by the kernel's own affine slot formula rather
         ! than a running counter, as the DF drivers do.
         if (ii == jj) then
            do c = 0, nj - 1
               nu = sh_ao(jj) + c
               do a = 0, ni - 1
                  mu = sh_ao(ii) + a
                  acc = acc + blk(1 + a + ni*c)*dmat(mu, nu)
               end do
            end do
         else
            do c = 0, nj - 1
               nu = sh_ao(jj) + c
               do a = 0, ni - 1
                  mu = sh_ao(ii) + a
                  acc = acc + blk(1 + a + ni*c)*(dmat(mu, nu) + dmat(nu, mu))
               end do
            end do
         end if
      end do
      v(k) = acc
   end subroutine pcm_potential_item

   subroutine pcm_fock_driver(nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                              pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, &
                              q, vmat)
      !! `vmat(mu,nu) += sum_k q_k <mu| erf(zeta_k|r - C_k|)/|r - C_k| |nu>`
      integer,  intent(in)    :: nshell, nao, npair, npp, npts
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: p_i(npair), p_j(npair), pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      real(dp), intent(in)    :: pts(3, npts), zk(npts)
      real(dp), intent(in)    :: q(npts)       !! weight per point
      real(dp), intent(inout) :: vmat(nao, nao)

      integer :: ij

      ! One thread per shell PAIR, walking every point: the pair owns its
      ! block and the transposed one, so the accumulation is race-free with
      ! no atomics, and the kernel sums the points before its single HRR.
      do concurrent(ij = 1:npair)
         call pcm_fock_item(ij, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                            pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, q, vmat)
      end do
   end subroutine pcm_fock_driver

   pure subroutine pcm_fock_item(ij, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                 pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, q, vmat)
      !$acc routine seq
      integer,  intent(in)    :: ij, nshell, nao, npair, npp, npts
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: p_i(npair), p_j(npair), pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      real(dp), intent(in)    :: pts(3, npts), zk(npts), q(npts)
      real(dp), intent(inout) :: vmat(nao, nao)

      integer  :: ii, jj, li, lj, ni, nj, a, c, mu, nu, o, n
      real(dp) :: blk(NBLK_MAX), abx, aby, abz

      ii = p_i(ij); jj = p_j(ij)
      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      abx = sh_r(1, ii) - sh_r(1, jj)
      aby = sh_r(2, ii) - sh_r(2, jj)
      abz = sh_r(3, ii) - sh_r(3, jj)
      o = pp_off(ij); n = pp_n(ij)
      call pcm_dispatch(li, lj, n, pp_z(o + 1:o + n), pp_k(o + 1:o + n), &
                        pp_p(:, o + 1:o + n), pp_a(:, o + 1:o + n), abx, aby, abz, &
                        npts, pts, zk, q, blk)
      ! On a diagonal pair the block already holds both (mu,nu) and
      ! (nu,mu), so it is added once; off the diagonal the transpose is
      ! this pair's to write too. Adding, not storing: the caller may be
      ! accumulating onto something.
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            vmat(mu, nu) = vmat(mu, nu) + blk(1 + a + ni*c)
         end do
      end do
      if (ii /= jj) then
         do c = 0, nj - 1
            nu = sh_ao(jj) + c
            do a = 0, ni - 1
               mu = sh_ao(ii) + a
               vmat(nu, mu) = vmat(nu, mu) + blk(1 + a + ni*c)
            end do
         end do
      end if
   end subroutine pcm_fock_item
#endif
