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
   use trc_xc_slater_exchange, only: slater_exchange_exc_vxc_unpolar, slater_exchange_exc_vxc_polar
   use trc_xc_vwn, only: vwn_exc_vxc_unpolar, vwn_exc_vxc_polar
   use trc_xc_vwn_rpa, only: vwn_rpa_exc_vxc_unpolar, vwn_rpa_exc_vxc_polar
   use trc_xc_pbe_x, only: pbe_x_exc_vxc_unpolar, pbe_x_exc_vxc_polar
   use trc_xc_pbe_c, only: pbe_c_exc_vxc_unpolar, pbe_c_exc_vxc_polar
   use trc_xc_b88, only: b88_exc_vxc_unpolar, b88_exc_vxc_polar
   use trc_xc_lyp, only: lyp_exc_vxc_unpolar, lyp_exc_vxc_polar
   use trc_xc_wb97x_v_xc, only: wb97x_v_xc_exc_vxc_unpolar, wb97x_v_xc_exc_vxc_polar
   use trc_xc_wb97x_xc, only: wb97x_xc_exc_vxc_unpolar, wb97x_xc_exc_vxc_polar
   use trc_xc_wb97_xc, only: wb97_xc_exc_vxc_unpolar, wb97_xc_exc_vxc_polar
   implicit none
   private

   public :: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, xc_functional_by_name
   public :: XC_FAMILY_LDA, XC_FAMILY_GGA, XC_MAX_KERNELS
   public :: XC_K_SLATER, XC_K_VWN, XC_K_VWN_RPA, XC_K_PBE_X, XC_K_PBE_C, XC_K_B88, XC_K_LYP
   public :: XC_K_WB97X_V, XC_K_WB97X, XC_K_WB97

   integer, parameter :: XC_FAMILY_LDA = 1, XC_FAMILY_GGA = 2
   integer, parameter :: XC_MAX_KERNELS = 8

   integer, parameter :: XC_K_SLATER = 1, XC_K_VWN = 2, XC_K_VWN_RPA = 3
   integer, parameter :: XC_K_PBE_X = 4, XC_K_PBE_C = 5, XC_K_B88 = 6, XC_K_LYP = 7
   integer, parameter :: XC_K_WB97X_V = 8
   integer, parameter :: XC_K_WB97X = 9, XC_K_WB97 = 10

   type :: trc_xc_functional_t
      character(len=16) :: name = ""
      integer :: family = XC_FAMILY_LDA
      integer :: nk = 0
      integer :: kid(XC_MAX_KERNELS) = 0
      real(dp) :: coef(XC_MAX_KERNELS) = 0.0_dp
      real(dp) :: exx = 0.0_dp   !! fraction of exact exchange the caller must add
      !
      ! Range-separated hybrids (currently only wB97X-V). `exx` keeps its
      ! meaning above unchanged: the fraction of FULL-range exact exchange,
      ! 1/r12. These four are zero for every other functional.
      !
      real(dp) :: exx_lr = 0.0_dp !! extra fraction of LONG-RANGE exact exchange, erf(omega r)/r
      real(dp) :: omega  = 0.0_dp !! range-separation parameter, bohr^-1; 0 = none
      real(dp) :: nlc_b  = 0.0_dp !! VV10 b; 0 = no nonlocal correlation
      real(dp) :: nlc_c  = 0.0_dp !! VV10 C
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
         case (XC_K_WB97X_V); call wb97x_v_xc_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case (XC_K_WB97X); call wb97x_xc_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case (XC_K_WB97); call wb97_xc_exc_vxc_unpolar(rho, sigma, e, vr, vs)
         case default; e = 0.0_dp; vr = 0.0_dp
         end select
         eps = eps + coef(k)*e
         vrho = vrho + coef(k)*vr
         vsigma = vsigma + coef(k)*vs
      end do
   end subroutine xc_eval_point

   !
   ! The same for a spin-polarised density: two densities, three gradient
   ! invariants (aa, ab, bb), and the derivatives with respect to each.
   ! `eps` is still per particle of the TOTAL density.
   !
   pure subroutine xc_eval_point_polar(nk, kid, coef, rho_a, rho_b, s_aa, s_ab, s_bb, &
                                       eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb)
      !$acc routine seq
      integer, intent(in) :: nk, kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS), rho_a, rho_b, s_aa, s_ab, s_bb
      real(dp), intent(out) :: eps, vrho_a, vrho_b, vs_aa, vs_ab, vs_bb
      real(dp) :: e, va, vb, saa, sab, sbb
      integer :: k

      eps = 0.0_dp; vrho_a = 0.0_dp; vrho_b = 0.0_dp
      vs_aa = 0.0_dp; vs_ab = 0.0_dp; vs_bb = 0.0_dp
      do k = 1, nk
         saa = 0.0_dp; sab = 0.0_dp; sbb = 0.0_dp
         select case (kid(k))
         case (XC_K_SLATER); call slater_exchange_exc_vxc_polar(rho_a, rho_b, e, va, vb)
         case (XC_K_VWN); call vwn_exc_vxc_polar(rho_a, rho_b, e, va, vb)
         case (XC_K_VWN_RPA); call vwn_rpa_exc_vxc_polar(rho_a, rho_b, e, va, vb)
         case (XC_K_PBE_X); call pbe_x_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case (XC_K_PBE_C); call pbe_c_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case (XC_K_B88); call b88_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case (XC_K_LYP); call lyp_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case (XC_K_WB97X_V)
            call wb97x_v_xc_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case (XC_K_WB97X)
            call wb97x_xc_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case (XC_K_WB97)
            call wb97_xc_exc_vxc_polar(rho_a, rho_b, s_aa, s_ab, s_bb, e, va, vb, saa, sab, sbb)
         case default; e = 0.0_dp; va = 0.0_dp; vb = 0.0_dp
         end select
         eps = eps + coef(k)*e
         vrho_a = vrho_a + coef(k)*va
         vrho_b = vrho_b + coef(k)*vb
         vs_aa = vs_aa + coef(k)*saa
         vs_ab = vs_ab + coef(k)*sab
         vs_bb = vs_bb + coef(k)*sbb
      end do
   end subroutine xc_eval_point_polar

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
      case ("wb97x-v", "wb97xv", "wb97x_v")
         ! libxc HYB_GGA_XC_WB97X_V. rsh_coeff (omega, alpha, beta) =
         ! (0.3, 1.0, -0.833): alpha is the LONG-RANGE exact-exchange
         ! fraction, hyb = alpha + beta the FULL-range one that multiplies
         ! SR_HFX = K_full - K_erf. pyscf/libxc build
         !   vk = hyb*K_full + (alpha - hyb)*K_erf(omega)
         ! which is exx*K_full + exx_lr*K_erf with exx = hyb, exx_lr = alpha
         ! - hyb. VV10 (nlc_b, nlc_c) from libxc's nlc_coeff, ((6.0, 0.01), 1).
         call set(f, XC_FAMILY_GGA, [XC_K_WB97X_V], [1.0_dp])
         f%exx = 0.167_dp
         f%exx_lr = 0.833_dp
         f%omega = 0.3_dp
         f%nlc_b = 6.0_dp
         f%nlc_c = 0.01_dp
      case ("wb97x")
         ! libxc HYB_GGA_XC_WB97X. rsh_coeff (omega, alpha, beta) =
         ! (0.3, 1.0, -0.842294); hyb = alpha + beta = 0.157706, exx_lr =
         ! alpha - hyb = 0.842294. No VV10 (nlc_coeff returns empty). The
         ! kernel's own translated `beta` parameter is -0.842294 exactly,
         ! the same direct cross-check wb97x-v's gave.
         call set(f, XC_FAMILY_GGA, [XC_K_WB97X], [1.0_dp])
         f%exx = 0.157706_dp
         f%exx_lr = 0.842294_dp
         f%omega = 0.3_dp
      case ("wb97")
         ! libxc HYB_GGA_XC_WB97. rsh_coeff (omega, alpha, beta) =
         ! (0.4, 1.0, -1.0); hyb = alpha + beta = 0, exx_lr = alpha - hyb =
         ! 1 -- a pure range-separated functional, no global HF mixing. No
         ! VV10. The kernel's own translated `beta` is -1.0 exactly.
         call set(f, XC_FAMILY_GGA, [XC_K_WB97], [1.0_dp])
         f%exx = 0.0_dp
         f%exx_lr = 1.0_dp
         f%omega = 0.4_dp
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
