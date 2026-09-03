!
! A density functional, as the integrator sees it: a short list of
! (coefficient, kernel) pairs, a family, and a fraction of exact exchange.
!
! The kernels are the translated ExchCXX routines in the sibling modules.
! This module knows their ids and how to call each one at a point; the
! compositions below are libxc's -- libxc is the reference every energy is
! checked against, so its definition of B3LYP is the one that counts, VWN
! variant included. ExchCXX's own B3LYP uses VWN5 and is `b3lyp5` here.
!
! `xc_eval_point` is the only thing that runs inside a kernel. It takes the
! composition as plain arrays rather than the type, because a derived type
! with allocatable components is exactly what a `do concurrent` body should
! not touch.
!
module trc_xc_functional
   use trc_boys, only: dp
   use trc_xc_slater_exchange, only: slater_exchange_exc_vxc_unpolar
   use trc_xc_vwn, only: vwn_exc_vxc_unpolar
   use trc_xc_vwn_rpa, only: vwn_rpa_exc_vxc_unpolar
   use trc_xc_pbe_x, only: pbe_x_exc_vxc_unpolar
   use trc_xc_pbe_c, only: pbe_c_exc_vxc_unpolar
   use trc_xc_b88, only: b88_exc_vxc_unpolar
   use trc_xc_lyp, only: lyp_exc_vxc_unpolar
   implicit none
   private

   public :: trc_xc_functional_t, xc_eval_point, xc_functional_by_name
   public :: XC_FAMILY_LDA, XC_FAMILY_GGA, XC_MAX_KERNELS
   public :: XC_K_SLATER, XC_K_VWN, XC_K_VWN_RPA, XC_K_PBE_X, XC_K_PBE_C, XC_K_B88, XC_K_LYP

   integer, parameter :: XC_FAMILY_LDA = 1, XC_FAMILY_GGA = 2
   integer, parameter :: XC_MAX_KERNELS = 8

   integer, parameter :: XC_K_SLATER = 1, XC_K_VWN = 2, XC_K_VWN_RPA = 3
   integer, parameter :: XC_K_PBE_X = 4, XC_K_PBE_C = 5, XC_K_B88 = 6, XC_K_LYP = 7

   type :: trc_xc_functional_t
      character(len=16) :: name = ""
      integer :: family = XC_FAMILY_LDA
      integer :: nk = 0
      integer :: kid(XC_MAX_KERNELS) = 0
      real(dp) :: coef(XC_MAX_KERNELS) = 0.0_dp
      real(dp) :: exx = 0.0_dp   !! fraction of exact exchange the caller must add
   end type trc_xc_functional_t

contains

   !
   ! Energy per particle and first derivatives of the composed functional at
   ! one point of a closed-shell density. `sigma` is |grad rho|^2 and is
   ! ignored by an LDA kernel.
   !
   pure subroutine xc_eval_point(nk, kid, coef, rho, sigma, eps, vrho, vsigma)
      !$acc routine seq
      integer, intent(in) :: nk, kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS), rho, sigma
      real(dp), intent(out) :: eps, vrho, vsigma
      real(dp) :: e, vr, vs
      integer :: k

      eps = 0.0_dp; vrho = 0.0_dp; vsigma = 0.0_dp
      do k = 1, nk
         vs = 0.0_dp
         select case (kid(k))
         case (XC_K_SLATER); call slater_exchange_exc_vxc_unpolar(rho, e, vr)
         case (XC_K_VWN); call vwn_exc_vxc_unpolar(rho, e, vr)
         case (XC_K_VWN_RPA); call vwn_rpa_exc_vxc_unpolar(rho, e, vr)
         case (XC_K_PBE_X); call pbe_x_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case (XC_K_PBE_C); call pbe_c_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case (XC_K_B88); call b88_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case (XC_K_LYP); call lyp_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case default; e = 0.0_dp; vr = 0.0_dp
         end select
         eps = eps + coef(k)*e
         vrho = vrho + coef(k)*vr
         vsigma = vsigma + coef(k)*vs
      end do
   end subroutine xc_eval_point

   !
   ! libxc's compositions. Names are case-insensitive.
   !
   function xc_functional_by_name(name, f) result(ok)
      character(len=*), intent(in) :: name
      type(trc_xc_functional_t), intent(out) :: f
      logical :: ok
      character(len=16) :: key
      integer :: i

      key = name
      do i = 1, len_trim(key)
         if (key(i:i) >= 'A' .and. key(i:i) <= 'Z') key(i:i) = achar(iachar(key(i:i)) + 32)
      end do
      f%name = key
      ok = .true.
      select case (trim(key))
      case ("lda_x", "slater")
         call set(f, XC_FAMILY_LDA, [XC_K_SLATER], [1.0_dp])
      case ("svwn", "lda", "svwn5")
         call set(f, XC_FAMILY_LDA, [XC_K_SLATER, XC_K_VWN], [1.0_dp, 1.0_dp])
      case ("svwn_rpa")
         call set(f, XC_FAMILY_LDA, [XC_K_SLATER, XC_K_VWN_RPA], [1.0_dp, 1.0_dp])
      case ("pbe")
         call set(f, XC_FAMILY_GGA, [XC_K_PBE_X, XC_K_PBE_C], [1.0_dp, 1.0_dp])
      case ("blyp")
         call set(f, XC_FAMILY_GGA, [XC_K_B88, XC_K_LYP], [1.0_dp, 1.0_dp])
      case ("pbe0")
         call set(f, XC_FAMILY_GGA, [XC_K_PBE_X, XC_K_PBE_C], [0.75_dp, 1.0_dp])
         f%exx = 0.25_dp
      case ("b3lyp")
         ! libxc HYB_GGA_XC_B3LYP: VWN_RPA in the correlation
         call set(f, XC_FAMILY_GGA, [XC_K_SLATER, XC_K_B88, XC_K_VWN_RPA, XC_K_LYP], &
                  [0.08_dp, 0.72_dp, 0.19_dp, 0.81_dp])
         f%exx = 0.20_dp
      case ("b3lyp5")
         ! libxc HYB_GGA_XC_B3LYP5 and ExchCXX's B3LYP: VWN5
         call set(f, XC_FAMILY_GGA, [XC_K_SLATER, XC_K_B88, XC_K_VWN, XC_K_LYP], &
                  [0.08_dp, 0.72_dp, 0.19_dp, 0.81_dp])
         f%exx = 0.20_dp
      case default
         ok = .false.
      end select
   contains
      subroutine set(f, family, kid, coef)
         type(trc_xc_functional_t), intent(inout) :: f
         integer, intent(in) :: family, kid(:)
         real(dp), intent(in) :: coef(:)
         f%family = family
         f%nk = size(kid)
         f%kid = 0; f%coef = 0.0_dp
         f%kid(1:f%nk) = kid
         f%coef(1:f%nk) = coef
      end subroutine set
   end function xc_functional_by_name

end module trc_xc_functional
