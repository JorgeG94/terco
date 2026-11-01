!
! The PCM surface-integral GRADIENT driver and its per-item routine, as an
! INCLUDE for the generated module that holds the kernels (trc_pcm_grad_kernels,
! from scripts/gen_1e.py --pcm-grad). Included there for the reason
! trc_1e_grad_driver.inc and trc_pcm_driver.inc both are: stdpar inlines only
! within a translation unit. NC_MAX, NBLK_MAX and PCM_GRAD_PT_CHUNK come from
! the including module.
!
! Combines trc_1e_grad_driver.inc's shell-pair walk over RAW primitives (the
! bra/ket basis derivative needs each primitive's own exponent, not just the
! pair's combined zeta, so this reads sh_e/sh_c directly rather than the pair
! list's precombined pp_* arrays) with trc_pcm_driver.inc's point-chunked work
! split (a cavity has 10^3 to 10^5 points; a thread per shell pair walking
! every point serialises the cavity, as trc_pcm_driver's own header measured
! for the Fock term). Each (shell pair, point chunk) item's own ga/gb/gc are
! lane-accumulated exactly as trc_1e_grad_driver.inc's gp(3,natm,NPRIV) does,
! since many pairs and many chunks can add into the same atom concurrently.
!
#ifndef TRC_PCM_GRAD_DRIVER_INC
#define TRC_PCM_GRAD_DRIVER_INC
   subroutine pcm_grad_driver(nshell, nao, natm, maxnp, sh_l, sh_np, sh_ao, sh_at, &
                              sh_e, sh_c, sh_r, npair, p_i, p_j, &
                              npts, pts, zk, ptatm, q, dmat, grad)
      !! +D:dV/dR of sum_k q_k erf(zeta_k|r-pts_k|)/|r-pts_k|, ADDED into
      !! `grad` -- same convention as trc_1e_gradient: the caller assembles
      !! the rest and adds this term in wherever PCM sits in the sum.
      integer,  intent(in)  :: nshell, nao, natm, maxnp, npair, npts
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell), sh_at(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      real(dp), intent(in)  :: pts(3, npts), zk(npts)
      integer,  intent(in)  :: ptatm(npts)
      real(dp), intent(in)  :: q(npts)
      real(dp), intent(in)  :: dmat(nao, nao)
      real(dp), intent(inout) :: grad(3, natm)

      integer :: ij, c, nch, ga_, gc_, gs
      real(dp) :: gsum
      real(dp), allocatable :: gp(:, :, :)
      !! Lane-private copies, trc_1e_grad_driver.inc's gp(3, natm, NPRIV)
      !! technique: a local allocatable, not a module array, same
      !! nvfortran reason that file's own comment gives.

      nch = min(npts, (npts + PCM_GRAD_PT_CHUNK - 1)/PCM_GRAD_PT_CHUNK)

      allocate (gp(3, natm, PCM_GRAD_NPRIV))
      !$acc enter data create(gp)
      do concurrent(gs = 1:PCM_GRAD_NPRIV, ga_ = 1:natm, gc_ = 1:3)
         gp(gc_, ga_, gs) = 0.0_dp
      end do

      do concurrent(ij = 1:npair, c = 1:nch)
         call pcm_grad_item(p_i(ij), p_j(ij), ij, c, nch, nshell, nao, natm, maxnp, &
                            sh_l, sh_np, sh_ao, sh_at, sh_e, sh_c, sh_r, &
                            npts, pts, zk, ptatm, q, dmat, gp)
      end do

      do concurrent(ga_ = 1:natm, gc_ = 1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, PCM_GRAD_NPRIV
            gsum = gsum + gp(gc_, ga_, gs)
         end do
         grad(gc_, ga_) = grad(gc_, ga_) + gsum
      end do
      !$acc exit data delete(gp)
      deallocate (gp)
   end subroutine pcm_grad_driver

   pure subroutine pcm_grad_item(ii, jj, ij, c, nch, nshell, nao, natm, maxnp, &
                                 sh_l, sh_np, sh_ao, sh_at, sh_e, sh_c, sh_r, &
                                 npts, pts, zk, ptatm, q, dmat, gp)
      !! One (shell pair, share of the points)'s PCM gradient block, into
      !! the lane `mod(ij-1, PCM_GRAD_NPRIV)+1`'s copy of `gp`. Share c of
      !! nch is the points c, c + nch, c + 2 nch, ..., dealt round-robin as
      !! trc_pcm_driver.inc's pcm_fock_item does. `wt` folds the triangular
      !! shell-pair list's double count into one scale, trc_1e_grad_item's
      !! own convention.
      !$acc routine seq
      integer,  value       :: ii, jj, ij, c, nch, nshell, nao, natm, maxnp, npts
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell), sh_at(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(in)  :: pts(3, npts), zk(npts), q(npts)
      integer,  intent(in)  :: ptatm(npts)
      real(dp), intent(in)  :: dmat(nao, nao)
      real(dp), intent(inout) :: gp(3, natm, PCM_GRAD_NPRIV)

      integer  :: li, lj, ni, nj, a, b, mu, nu, k, lane, npt, k0
      real(dp) :: wt, ga(3), gb(3)
      real(dp) :: dblk(NBLK_MAX)

      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      wt = merge(1.0_dp, 2.0_dp, ii == jj)
      lane = mod(ij - 1, PCM_GRAD_NPRIV) + 1
      npt = (npts - c)/nch + 1
      k0 = c

      k = 0
      do b = 0, nj - 1
         nu = sh_ao(jj) + b
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            dblk(k) = dmat(mu, nu)
         end do
      end do

      ! `gc`, the point-derivative term's own-atom output, is a SLICE of
      ! this thread's own lane straight out of `gp`, same reason
      ! trc_1e_grad_driver.inc's one_e_grad_item gives for its gc slice: a
      ! local gc(3, natm) with natm a runtime extent would be a per-thread
      ! dynamic stack allocation that the resident thread count overruns.
      call pcm_grad_dispatch(li, lj, sh_np(ii), sh_np(jj), &
                             sh_e(:, ii), sh_c(:, ii), &
                             sh_e(:, jj), sh_c(:, jj), &
                             sh_r(:, ii), sh_r(:, jj), &
                             npt, k0, nch, pts, zk, ptatm, q, &
                             natm, wt, dblk, ga, gb, gp(1:3, 1:natm, lane))

      ! ga/gb ATOMIC for the same reason trc_1e_grad_driver.inc's own ga/gb
      ! scatter is: two (shell pair, chunk) items sharing both a lane and
      ! an atom run this read-modify-write concurrently.
      !$acc atomic update
      gp(1, sh_at(ii), lane) = gp(1, sh_at(ii), lane) + ga(1)
      !$acc atomic update
      gp(2, sh_at(ii), lane) = gp(2, sh_at(ii), lane) + ga(2)
      !$acc atomic update
      gp(3, sh_at(ii), lane) = gp(3, sh_at(ii), lane) + ga(3)
      !$acc atomic update
      gp(1, sh_at(jj), lane) = gp(1, sh_at(jj), lane) + gb(1)
      !$acc atomic update
      gp(2, sh_at(jj), lane) = gp(2, sh_at(jj), lane) + gb(2)
      !$acc atomic update
      gp(3, sh_at(jj), lane) = gp(3, sh_at(jj), lane) + gb(3)
   end subroutine pcm_grad_item
#endif
