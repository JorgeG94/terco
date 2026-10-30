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
! THE WORK SPLIT. Each driver is one `do concurrent` over the product of the
! two axes, one of them chunked, with atomics where items meet. Over either
! axis alone the launch is too narrow: on w20/6-31G* (9169 points, 6.9k
! surviving shell pairs, V100) a thread per pair walking every point took
! 666 ms for the Fock term -- 54 blocks, under one per SM -- and a thread
! per point walking every pair 32 ms for the potential, against a 22 ms
! Fock build. Chunked, both run at full width.
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

      integer :: k, c, nch
      real(dp), allocatable :: part(:, :)

      ! One item per (pair chunk, point), the point fastest: a warp shares
      ! its pairs and so walks the same classes in step. Each item leaves
      ! its chunk's sum in part(k, c), and a second pass adds the chunks up
      ! in order -- deterministic, and not an atomic into v(k): that was
      ! tried, and with v(k) indexed by the loop's own index nvfortran gave
      ! wrong sums (w5/6-31G* at E_diel = -4580 Eh on the GPU, right on the
      ! host, and right on both with a single chunk).
      nch = (npair + PCM_PAIR_CHUNK - 1)/PCM_PAIR_CHUNK
      allocate (part(npts, nch))
      !$acc enter data create(part)
      do concurrent(c = 1:nch, k = 1:npts)
         call pcm_potential_item(k, c, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                 pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, &
                                 dmat, part(k, c))
      end do
      do concurrent(k = 1:npts)
         v(k) = sum_chunks(nch, npts, part, k)
      end do
      !$acc exit data delete(part)
   end subroutine pcm_potential_driver

   pure real(dp) function sum_chunks(nch, npts, part, k)
      !$acc routine seq
      integer, intent(in) :: nch, npts, k
      real(dp), intent(in) :: part(npts, nch)
      integer :: c
      sum_chunks = 0.0_dp
      do c = 1, nch
         sum_chunks = sum_chunks + part(k, c)
      end do
   end function sum_chunks

   pure subroutine pcm_potential_item(k, c, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                      pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, &
                                      dmat, vpart)
      !$acc routine seq
      integer,  value         :: k, c, nshell, nao, npair, npp, npts
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: p_i(npair), p_j(npair), pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      real(dp), intent(in)    :: pts(3, npts), zk(npts)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(out)   :: vpart        !! this chunk's share of v(k)

      integer  :: ij, ii, jj, li, lj, ni, nj, a, b, mu, nu, o, n
      real(dp) :: block(NBLK_MAX), one(1), acc, abx, aby, abz

      one(1) = 1.0_dp
      acc = 0.0_dp
      do ij = (c - 1)*PCM_PAIR_CHUNK + 1, min(c*PCM_PAIR_CHUNK, npair)
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
                           1, 1, 1, pts(:, k:k), zk(k:k), one, block)
         ! The pair list is lower-triangular, and an off-diagonal block
         ! stands for itself and its transpose; D is taken as given rather
         ! than assumed symmetric, so both of its halves are read.
         ! `block` addressed by the kernel's own affine slot formula rather
         ! than a running counter, as the DF drivers do.
         if (ii == jj) then
            do b = 0, nj - 1
               nu = sh_ao(jj) + b
               do a = 0, ni - 1
                  mu = sh_ao(ii) + a
                  acc = acc + block(1 + a + ni*b)*dmat(mu, nu)
               end do
            end do
         else
            do b = 0, nj - 1
               nu = sh_ao(jj) + b
               do a = 0, ni - 1
                  mu = sh_ao(ii) + a
                  acc = acc + block(1 + a + ni*b)*(dmat(mu, nu) + dmat(nu, mu))
               end do
            end do
         end if
      end do
      vpart = acc
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

      integer :: ij, c, nch

      ! One item per (shell pair, share of the points), the share fastest so
      ! a warp shares its pair and class. Share c of nch is the points c,
      ! c + nch, c + 2 nch, ...: dealt round-robin, so the warp's threads
      ! read neighbouring points at every step. The kernel sums its share
      ! before its single HRR, and the block goes into vmat by atomics,
      ! since every share of a pair adds into the same block.
      nch = min(npts, (npts + PCM_PT_CHUNK - 1)/PCM_PT_CHUNK)
      do concurrent(ij = 1:npair, c = 1:nch)
         call pcm_fock_item(ij, c, nch, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                            pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, q, vmat)
      end do
   end subroutine pcm_fock_driver

   pure subroutine pcm_fock_item(ij, c, nch, nshell, nao, sh_l, sh_ao, sh_r, npair, npp, p_i, p_j, &
                                 pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, npts, pts, zk, q, vmat)
      !$acc routine seq
      integer,  value         :: ij, c, nch, nshell, nao, npair, npp, npts
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      integer,  intent(in)    :: p_i(npair), p_j(npair), pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      real(dp), intent(in)    :: pts(3, npts), zk(npts), q(npts)
      real(dp), intent(inout) :: vmat(nao, nao)

      integer  :: ii, jj, li, lj, ni, nj, a, b, mu, nu, o, n
      real(dp) :: block(NBLK_MAX), abx, aby, abz

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
                        (npts - c)/nch + 1, c, nch, pts, zk, q, block)
      ! On a diagonal pair the block already holds both (mu,nu) and
      ! (nu,mu), so it is added once; off the diagonal the transpose is
      ! this pair's to add too. Adding, not storing: other chunks of the
      ! same pair, and the caller, may be accumulating onto it.
      do b = 0, nj - 1
         nu = sh_ao(jj) + b
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            !$acc atomic update
            vmat(mu, nu) = vmat(mu, nu) + block(1 + a + ni*b)
         end do
      end do
      if (ii /= jj) then
         do b = 0, nj - 1
            nu = sh_ao(jj) + b
            do a = 0, ni - 1
               mu = sh_ao(ii) + a
               !$acc atomic update
               vmat(nu, mu) = vmat(nu, mu) + block(1 + a + ni*b)
            end do
         end do
      end if
   end subroutine pcm_fock_item
#endif
