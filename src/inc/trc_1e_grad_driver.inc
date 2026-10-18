!
! The one-electron GRADIENT driver and its per-item routine, as an INCLUDE
! for the module that holds the kernels they dispatch to -- same reason as
! trc_1e_driver.inc: stdpar only inlines within a translation unit, so this
! has to live beside one_e_grad_dispatch rather than call across the module
! boundary from trc_1e_grad_kernels' including module. ONE_E_GRAD_NPRIV and
! NBLK_MAX come from the including module.
!
#ifndef TRC_1E_GRAD_DRIVER_INC
#define TRC_1E_GRAD_DRIVER_INC
   subroutine one_e_grad_driver(nshell, nao, natm, maxnp, sh_l, sh_np, sh_ao, sh_at, &
                                sh_e, sh_c, sh_r, at_z, at_r, &
                                npair, p_i, p_j, wmat, dmat, grad)
      !! -W:dS/dR + D:dT/dR + D:dV/dR, over the same triangular shell-pair
      !! list trc_1e's own driver walks, added into `grad` (already holding
      !! whatever eri%gradient and the nuclear-repulsion term put there).
      integer,  intent(in)  :: nshell, nao, natm, maxnp, npair
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell), sh_at(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(in)  :: at_z(natm), at_r(3, natm)
      real(dp), intent(in)  :: wmat(nao, nao), dmat(nao, nao)
      real(dp), intent(inout) :: grad(3, natm)

      integer :: ij, ga_, gc_, gs
      real(dp) :: gsum
      real(dp), allocatable :: gp(:, :, :)
      !! Lane-private copies, exactly trc_grad.F90's gp(3, natm, GRAD_NPRIV)
      !! technique: multiple shell pairs can add into the SAME atom's
      !! gradient concurrently (unlike trc_1e's smat/tmat/vmat, whose pairs
      !! write disjoint matrix blocks), so the accumulation needs either
      !! this or an atomic; a local allocatable, not a module array, for
      !! the same nvfortran hang/ICE reason trc_grad.F90's own comment
      !! gives for its gp.

      allocate (gp(3, natm, ONE_E_GRAD_NPRIV))
      !$acc enter data create(gp)
      do concurrent(gs=1:ONE_E_GRAD_NPRIV, ga_=1:natm, gc_=1:3)
         gp(gc_, ga_, gs) = 0.0_dp
      end do

      do concurrent(ij=1:npair)
         call one_e_grad_item(p_i(ij), p_j(ij), ij, nshell, nao, natm, maxnp, &
                              sh_l, sh_np, sh_ao, sh_at, sh_e, sh_c, sh_r, &
                              at_z, at_r, wmat, dmat, gp)
      end do

      do concurrent(ga_=1:natm, gc_=1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, ONE_E_GRAD_NPRIV
            gsum = gsum + gp(gc_, ga_, gs)
         end do
         grad(gc_, ga_) = grad(gc_, ga_) + gsum
      end do
      !$acc exit data delete(gp)
      deallocate (gp)
   end subroutine one_e_grad_driver

   pure subroutine one_e_grad_item(ii, jj, ij, nshell, nao, natm, maxnp, sh_l, sh_np, &
                                   sh_ao, sh_at, sh_e, sh_c, sh_r, at_z, at_r, &
                                   wmat, dmat, gp)
      !! One triangular shell pair's full S+T+V gradient contribution,
      !! into the lane `mod(ij-1, ONE_E_GRAD_NPRIV)+1`'s copy of `gp`. `wt`
      !! folds in the factor of 2 an off-diagonal (ii /= jj) pair needs for
      !! BOTH (mu,nu) and (nu,mu) in the full sum_{mu,nu} -- the diagonal
      !! block's own nested AO loop already covers its full square once
      !! (see one_e_item's comment), and a same-atom pair's ga+gb cancel
      !! exactly for S/T (translational invariance) by just being added to
      !! the same atom, with no special case needed.
      !$acc routine seq
      integer,  intent(in) :: ii, jj, ij, nshell, nao, natm, maxnp
      integer,  intent(in) :: sh_l(nshell), sh_np(nshell), sh_ao(nshell), sh_at(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in) :: sh_r(3, nshell)
      real(dp), intent(in) :: at_z(natm), at_r(3, natm)
      real(dp), intent(in) :: wmat(nao, nao), dmat(nao, nao)
      real(dp), intent(inout) :: gp(3, natm, ONE_E_GRAD_NPRIV)

      integer  :: li, lj, ni, nj, a, c, mu, nu, k, lane
      real(dp) :: wt, ga(3), gb(3), gc(3, natm)
      real(dp) :: wblk(NBLK_MAX), dblk(NBLK_MAX)

      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      wt = merge(1.0_dp, 2.0_dp, ii == jj)

      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            wblk(k) = wmat(mu, nu)
            dblk(k) = dmat(mu, nu)
         end do
      end do

      call one_e_grad_dispatch(li, lj, sh_np(ii), sh_np(jj), &
                               sh_e(:, ii), sh_c(:, ii), &
                               sh_e(:, jj), sh_c(:, jj), &
                               sh_r(:, ii), sh_r(:, jj), &
                               natm, at_z, at_r, wt, wblk, dblk, ga, gb, gc)

      lane = mod(ij - 1, ONE_E_GRAD_NPRIV) + 1
      gp(1:3, sh_at(ii), lane) = gp(1:3, sh_at(ii), lane) + ga
      gp(1:3, sh_at(jj), lane) = gp(1:3, sh_at(jj), lane) + gb
      gp(1:3, 1:natm, lane) = gp(1:3, 1:natm, lane) + gc
   end subroutine one_e_grad_item
#endif
