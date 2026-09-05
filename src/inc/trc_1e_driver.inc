!
! The one-electron driver and its per-item routine, as an INCLUDE for the
! module that holds the kernels they dispatch to. They used to live in
! trc_api and call the dispatcher across the module boundary, which is a
! real device call per item -- stdpar only inlines within a translation
! unit. Included after `contains` in the kernel module, they inline; the
! guard keeps a unit from getting two copies. NC_MAX and the NBLK sizes
! come from the including module.
!
#ifndef TRC_1E_DRIVER_INC
#define TRC_1E_DRIVER_INC
   subroutine one_e_driver(nshell, nao, natm, maxnp, sh_l, sh_np, sh_ao, &
                           sh_e, sh_c, sh_r, at_z, at_r, &
                           npair, p_i, p_j, smat, tmat, vmat)
      integer,  intent(in)  :: nshell, nao, natm, maxnp, npair
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(in)  :: at_z(natm), at_r(3, natm)
      real(dp), intent(out) :: smat(nao, nao), tmat(nao, nao), vmat(nao, nao)

      integer :: ij, i, j

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(j=1:nao, i=1:nao)
         smat(i, j) = 0.0_dp
         tmat(i, j) = 0.0_dp
         vmat(i, j) = 0.0_dp
      end do

      ! One thread per SURVIVING shell pair. The list is lower-triangular and
      ! these matrices are symmetric, so each thread writes its block and the
      ! transpose -- which is why the earlier full-square version existed. The
      ! pair list changes the arithmetic: with screening, walking the square
      ! costs a factor of two on top of computing pairs that are known to be
      ! negligible.
      !
      ! The body lives in a `pure` item routine rather than inline, for two
      ! reasons: nvfortran rejects a BLOCK construct inside `do concurrent`,
      ! and declaring the scratch alongside the loop would make it SHARED
      ! across threads -- the failure that silently corrupted the four-centre
      ! per-class kernels once already.
      do concurrent(ij = 1:npair)
         call one_e_item(p_i(ij), p_j(ij), nshell, nao, natm, maxnp, sh_l, &
                         sh_np, sh_ao, sh_e, sh_c, sh_r, at_z, at_r, &
                         smat, tmat, vmat)
      end do
   end subroutine one_e_driver

   pure subroutine one_e_item(ii, jj, nshell, nao, natm, maxnp, sh_l, sh_np, &
                              sh_ao, sh_e, sh_c, sh_r, at_z, at_r, &
                              smat, tmat, vmat)
      !$acc routine seq
      integer,  intent(in)    :: ii, jj, nshell, nao, natm, maxnp
      integer,  intent(in)    :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      real(dp), intent(in)    :: at_z(natm), at_r(3, natm)
      real(dp), intent(inout) :: smat(nao, nao), tmat(nao, nao), vmat(nao, nao)

      integer  :: li, lj, ni, nj, a, c, mu, nu, k
      real(dp) :: sblk(NBLK_MAX), tblk(NBLK_MAX), vblk(NBLK_MAX)

      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      call one_e_dispatch(li, lj, sh_np(ii), sh_np(jj), &
                          sh_e(:, ii), sh_c(:, ii), &
                          sh_e(:, jj), sh_c(:, jj), &
                          sh_r(:, ii), sh_r(:, jj), &
                          natm, at_z, at_r, sblk, tblk, vblk)
      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            smat(mu, nu) = sblk(k)
            tmat(mu, nu) = tblk(k)
            vmat(mu, nu) = vblk(k)
            ! symmetric counterpart; on the diagonal block (ii == jj) this
            ! rewrites the same element with the same value, which is why it
            ! needs no guard
            smat(nu, mu) = sblk(k)
            tmat(nu, mu) = tblk(k)
            vmat(nu, mu) = vblk(k)
         end do
      end do
   end subroutine one_e_item
#endif
