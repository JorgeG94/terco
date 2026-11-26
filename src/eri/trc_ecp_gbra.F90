!
! The ECP gradient's kernels for a g derivative bra: the raised shell of an f
! basis function, against a ket of f or below. The same items as trc_ecp's
! (inc/trc_ecp_block.inc), instantiated at the larger bra sizes in a module of
! their own, so the energy matrix and the gradient of a basis without f keep
! their column sizes and per-thread stack (inc/trc_ecp_block_sizes.inc).
!
! Both public routines launch a kernel for `trc_ecp`'s gradient: every array
! they read is mapped by it already.
!
#include "inc/trc_dc_index.inc"

module trc_ecp_gbra
   use pic_types, only: dp, int64
   use trc_ecp_quadrature, only: ECP_LMAX, ECP_NRS
   use trc_ecp_angular, only: ECP_ORB_LMAX, ECP_BRA_LMAX, CART_POW_Y, CART_POW_Z, &
                              ecp_type1_rad_ang, ecp_type2_facs_ang
   use trc_ecp_grad_types, only: trc_ecp_channels_t, trc_ecp_deriv_basis_t, trc_ecp_grad_tasks_t, &
                                 trc_ecp_k1b_cache_t, trc_ecp_chunk_t
   use trc_ecp_radial, only: ecp_radial_type1_static_facs, ecp_radial_type1_rad_part_levels, &
                             ecp_radial_type2_facs_rad, ecp_radial_type2_lab_levels, ECP_NLEVEL
   implicit none
   private
   public :: ECP_GBRA_FAC1, ECP_GBRA_FAC2, ECP_GBRA_ANG2, ECP_GBRA_COL_BYTES
   public :: trc_ecp_gbra_caches, trc_ecp_gbra_tasks

   integer, parameter :: ECP_BLOCK_LBRA = ECP_BRA_LMAX, ECP_BLOCK_LKET = ECP_ORB_LMAX
   !! This module's instance of the block kernels (trc_ecp_block_sizes.inc).
#include "inc/trc_ecp_block_sizes.inc"

   integer, parameter :: ECP_GBRA_FAC1 = FAC1SIZE_BRA   !! K1b column extents of a g shell
   integer, parameter :: ECP_GBRA_FAC2 = FAC2SIZE_BRA
   integer, parameter :: ECP_GBRA_ANG2 = ANG2SIZE_BRA
   integer(int64), parameter :: ECP_GBRA_COL_BYTES = BRA_COL_BYTES

contains

#include "inc/trc_ecp_block.inc"

   subroutine trc_ecp_gbra_caches(db, ch, cn, cache)
      !! K1b over the chunk's g shell columns, into `cache` (allocated and
      !! mapped, at this module's bra column sizes).
      type(trc_ecp_deriv_basis_t), intent(in) :: db
      type(trc_ecp_channels_t), intent(in) :: ch
      type(trc_ecp_chunk_t), intent(in) :: cn
      type(trc_ecp_k1b_cache_t), intent(inout) :: cache

      call gbra_k1b_kernel(cache%ncol, cn%gcol_sh, cn%gcol_g, db%nshell, db%l, db%np, db%maxnp, db%e, db%c, &
                           db%r, db%natm, db%at_r, ch%ngroup, ch%grp_atom, ch%grp_l, ch%rs, &
                           cache%fac1, cache%fac2, cache%ang2)
   end subroutine trc_ecp_gbra_caches

   subroutine trc_ecp_gbra_tasks(db, ch, tk, density, bra, ket, cn)
      !! The task kernel over the chunk's g-bra tasks, into their slots of
      !! `cn%w` and `cn%conv_fail`: the bra's columns in the g-sized `bra`
      !! caches, the ket's in the f-sized `ket` ones.
      type(trc_ecp_deriv_basis_t), intent(in) :: db
      type(trc_ecp_channels_t), intent(in) :: ch
      type(trc_ecp_grad_tasks_t), intent(in) :: tk
      real(dp), intent(in) :: density(db%nao, db%nao)
      type(trc_ecp_k1b_cache_t), intent(in) :: bra, ket
      type(trc_ecp_chunk_t), intent(inout) :: cn

      call gbra_task_kernel(cn%ntg, cn%nt, cn%tg0, tk%ntask, tk%bra, tk%ket, tk%atm, cn%ci, cn%cj, db%nreal, &
                            ch%natm_ecp, ch%atm_ecp, ch%atm_grp0, ch%atm_grp1, bra%ncol, ket%ncol, &
                            ch%ngroup, ch%grp_l, db%nshell, db%l, db%np, db%maxnp, db%e, db%c, db%r, &
                            db%natm, db%at_r, ch%rs, ch%ws, ch%ur, bra%fac1, bra%fac2, bra%ang2, &
                            ket%fac1, ket%fac2, ket%ang2, db%nao, density, db%nfold, db%parent, db%coef, &
                            NFKET, db%ket_ao, cn%nt + cn%ntg, cn%w, cn%conv_fail)
   end subroutine trc_ecp_gbra_tasks

   subroutine gbra_k1b_kernel(ncol, col_sh, col_g, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                  natm, at_r, ngroup, grp_atom, grp_l, rs, fac1, fac2, ang2)
      !! `trc_ecp_gbra_caches`' kernel.
      integer, intent(in) :: ncol, nshell, maxnp, natm, ngroup
      integer, intent(in) :: col_sh(ncol), col_g(ncol)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: at_r(3, natm)
      integer, intent(in) :: grp_atom(ngroup), grp_l(ngroup)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      real(dp), intent(inout) :: fac1(0:FAC1SIZE_BRA - 1, ncol), fac2(0:FAC2SIZE_BRA - 1, ncol)
      real(dp), intent(inout) :: ang2(0:ANG2SIZE_BRA - 1, ncol)
      integer :: i, g, k

      do concurrent(k=1:ncol) local(i, g)
         i = col_sh(k)
         g = col_g(k)
         call ecp_k1b_item(TRC_DC_INDEX(i), TRC_DC_INDEX(g), nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                           natm, at_r, ngroup, grp_atom, grp_l, rs, fac1(:, k), fac2(:, k), ang2(:, k))
      end do
   end subroutine gbra_k1b_kernel

   subroutine gbra_task_kernel(nt, slot0, task0, ntask, bra, ket, atm, ci, cj, nreal, natm_ecp, &
                                 atm_ecp, atm_grp0, atm_grp1, nbra, nket, ngroup, grp_l, nshell, &
                                 sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, natm, at_r, rs, ws, ur_cache, &
                                 fac1_bra, fac2_bra, ang2_bra, fac1_ket, fac2_ket, ang2_ket, nao, density, &
                                 nfold, parent, coef, nfket_ao, ket_ao, nslot, w, conv_fail)
      !! `trc_ecp_gbra_tasks`' kernel: task `task0 + t` of the list
      !! (bra, ket, atm), t = 1..nt, into slot `slot0 + t` of `ci`, `cj`, `w`
      !! and `conv_fail`.
      integer, intent(in) :: nt, slot0, task0, ntask, nreal, natm_ecp, nbra, nket, ngroup, nshell, maxnp
      integer, intent(in) :: natm, nao, nfold, nfket_ao, nslot
      integer, intent(in) :: bra(ntask), ket(ntask), atm(ntask)
      integer, intent(in) :: ci(nslot), cj(nslot)
      integer, intent(in) :: atm_ecp(natm_ecp), atm_grp0(natm_ecp), atm_grp1(natm_ecp)
      integer, intent(in) :: grp_l(ngroup)
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: at_r(3, natm)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), intent(in) :: ur_cache(0:ECP_NRS - 1, ngroup)
      real(dp), intent(in) :: fac1_bra(0:FAC1SIZE_BRA - 1, nbra), fac2_bra(0:FAC2SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: ang2_bra(0:ANG2SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: fac1_ket(0:FAC1SIZE_KET - 1, nket), fac2_ket(0:FAC2SIZE_KET - 1, nket)
      real(dp), intent(in) :: ang2_ket(0:ANG2SIZE_KET - 1, nket)
      real(dp), intent(in) :: density(nao, nao)
      integer, intent(in) :: parent(3, nfold, nshell - nreal)
      real(dp), intent(in) :: coef(3, nfold, nshell - nreal)
      integer, intent(in) :: ket_ao(nfket_ao, nreal)
      real(dp), intent(inout) :: w(3, nslot)
      integer, intent(inout) :: conv_fail(nslot)
      integer :: t, k, s, i, j, a

      do concurrent(t=1:nt) local(k, s, i, j, a)
         k = task0 + t
         s = slot0 + t
         i = bra(k)
         j = ket(k)
         a = atm(k)
         call ecp_grad_item(w(:, s), conv_fail(s), TRC_DC_INDEX(i), TRC_DC_INDEX(j), nreal, atm_ecp(a), &
                            atm_grp0(a), atm_grp1(a), ci(s), cj(s), nbra, nket, ngroup, grp_l, nshell, &
                            sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, natm, at_r, rs, ws, ur_cache, &
                            fac1_bra, fac2_bra, ang2_bra, fac1_ket, fac2_ket, ang2_ket, &
                            nao, density, nfold, parent, coef, nfket_ao, ket_ao)
      end do
   end subroutine gbra_task_kernel

end module trc_ecp_gbra
