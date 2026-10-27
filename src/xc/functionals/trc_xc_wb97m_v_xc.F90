!
! WB97M_V_XC: range-separated hybrid meta-GGA exchange-correlation kernel.
!
! HAND-WRITTEN, not translated: ExchCXX (~/dev/ExchCXX, pinned commit
! 67be5c6) has no wB97M-V kernel (checked against its upstream HEAD, which
! is the same commit, so there is nothing newer to pull either). The
! functional form and every numeric coefficient below are Mardirossian and
! Head-Gordon's own, J. Chem. Phys. 144, 214110 (2016) ("wB97M-V: A
! combinatorially-optimized, range-separated hybrid, meta-GGA density
! functional with VV10 nonlocal correlation"), Eqs. 3-28 and the "6
! (Final)" column of Table 3. Nothing here is translated or copied from
! libxc (MPL-2.0): libxc/pyscf served only as the NUMERICAL reference
! check_xc.F90 compares against (reference_xc.dat, kernel "wb97m_v_xc").
!
! This kernel is the SEMILOCAL exchange-correlation only -- Eqs. 4-5's
! E_x^mGGA + E_c,ss^mGGA + E_c,os^mGGA. The exact-exchange terms
! (c_x*E_x,sr^exact + E_x,lr^exact) and VV10's E_c,nl are added by the
! caller through trc_xc_functional_t's exx/exx_lr/omega/nlc_b/nlc_c, the
! same split wb97x_v_xc uses for its GGA counterpart.
!
! The functional (per spin sigma in {a,b}, Eqs. 6-28):
!
!   t_sigma  = tau_sigma^UEG / tau_sigma,   tau_sigma^UEG = (3/10)(6 pi^2)^(2/3) rho_sigma^(5/3)
!   w_sigma  = (t_sigma - 1)/(t_sigma + 1)                                    [meta-GGA ingredient]
!   s2_sigma = |grad rho_sigma|^2 / rho_sigma^(8/3)
!   u_sigma  = gamma * s2_sigma / (1 + gamma * s2_sigma)                      [GGA ingredient]
!
! with gamma = 0.004 (exchange), 0.2 (same-spin correlation), 0.006
! (opposite-spin, on s2_ab = (s2_a + s2_b)/2 and t_ab = (t_a + t_b)/2).
! Exchange and same-spin correlation sum their two spins; opposite-spin is
! one Stoll-decomposed (Eq. 20) term built from the TOTAL PW92 correlation
! energy density, e_c,ab = rho*eps_c^PW92(rho_a,rho_b) - rho_a*eps_c^PW92
! (rho_a,0) - rho_b*eps_c^PW92(0,rho_b). Every eps_c^PW92 value and its
! rho-derivative is `pw91_lda_mod_exc_vxc_polar`'s -- the MODIFIED PW92
! high-density correction, not the plain one (`trc_xc_pw91_lda`): checked
! numerically against pyscf's WB97M_V at a w = 0 point, where every other
! term is an exact constant (g = c_*00) and the plain PW92 missed by
! ~3e-7 while the modified one matched to ~1e-13 -- libxc's own wb97m_v.c
! uses LDA_C_PW_MOD internally, not LDA_C_PW, and this kernel now does
! too. That translated ExchCXX kernel already matches libxc's
! LDA_C_PW_MOD at check_xc's tolerance (not in `DISAGREES_WITH_LIBXC`),
! so this kernel reuses it rather than re-deriving Perdew and Wang's
! G(rs,zeta) closed form.
!
! Each inhomogeneity-correction factor is a short power series in w and u,
! Eq. 2's f(x,y) with terco's selectable coefficients from Table 3's final
! column (c_x,ij, c_css,ij, c_cos,ij; every other one of the 135 is zero):
!
!   g_x(w,u)   = c_x00 + c_x10*w + c_x01*u
!   g_css(w,u) = c_css00 + c_css10*w + c_css20*w^2 + c_css43*w^4*u^3 + c_css04*u^4
!   g_cos(w,u) = c_cos00 + c_cos10*w + c_cos20*w^2 + c_cos60*w^6
!                + c_cos21*w^2*u + c_cos61*w^6*u
!
! No g touches sigma_ab (Eqs. 15, 24, 28 use only the diagonal s2_sigma),
! so vsigma_ab is identically zero -- not screened out, just never
! computed. Every analytic derivative below (vrho, vsigma, vtau per
! spin/pair, and the unpolarised reduction two paragraphs down) was cross-
! checked against central finite differences at 1e-6 before being
! transcribed; see the commit message for the script.
!
! `exc_vxc_polar`'s core is the general two-spin formula; `exc_vxc_unpolar`
! calls it at the symmetric point rho_a = rho_b = rho/2, sigma_aa = sigma_ab
! = sigma_bb = sigma/4, tau_a = tau_b = tau/2 (an unpolarised density's
! actual spin-resolution) and reduces the five polarised outputs by the
! chain rule through that substitution:
!
!   eps   = eps_polar                          (same normalisation, no rho_b dependence left)
!   vrho   = (vrho_a + vrho_b)/2
!   vsigma = (vsigma_aa + vsigma_ab + vsigma_bb)/4
!   vtau   = (vtau_a + vtau_b)/2
!
! because d(rho_a)/d(rho) = d(rho_b)/d(rho) = 1/2 and d(sigma_ss)/d(sigma)
! = d(sigma_ab)/d(sigma) = 1/4 along that submanifold, with sigma_aa,
! sigma_ab, sigma_bb independent of rho there (grad rho_a = grad rho_b =
! grad rho / 2 regardless of rho itself). This is cheaper to get right
! than a second from-scratch derivation and was validated the same way
! (finite differences against the general formula, not assumed).
!
! Only `exc`/`exc_vxc` exist, the same as every other translated meta-GGA
! kernel (tools/exchcxx_to_fortran.py's own comment): no fxc, so a
! coupled-perturbed solve or an analytic gradient must refuse this
! functional the way it already refuses every other meta-GGA.
!
module trc_xc_wb97m_v_xc
   use pic_types, only: dp
   use trc_xc_pw91_lda_mod, only: pw91_lda_mod_exc_vxc_polar
   use trc_xc_util, only: m_cbrt_3, m_cbrt_4, m_cbrt_one_ov_pi, m_pi, xc_cbrt, xc_erf
   implicit none
   private

   public :: wb97m_v_xc_exc_unpolar
   public :: wb97m_v_xc_exc_vxc_unpolar
   public :: wb97m_v_xc_exc_polar
   public :: wb97m_v_xc_exc_vxc_polar

   ! libxc's own input sanitisation for a meta-GGA (work_mgga_inc.c's
   ! WORK_MGGA_KERNEL, read -- not translated -- purely to match its
   ! thresholds, since check_xc compares against its numbers): every
   ! meta-GGA floors each spin's rho at `dens_threshold`, tau at
   ! `tau_threshold`, and sigma at `dens_threshold^(4/3)` squared, before
   ! the kernel ever sees them. wB97M-V's own dens_threshold is 1e-13
   ! (hyb_mgga_xc_wb97mv.c's xc_func_info_type; `functionals.c` gives
   ! every meta-GGA the same tau_threshold = 1e-20, and
   ! sigma_threshold = dens_threshold^(4/3)). wB97M-V does NOT carry
   ! XC_FLAGS_ENFORCE_FHC (checked the same way), so the separate
   ! sigma <= 8 rho tau clamp SCAN-family functionals get does not apply
   ! here -- confirming the paper's own functional form (no such bound in
   ! Eqs. 3-28) needed no code to enforce it.
   !
   ! Getting this wrong was a real bug, not a tolerance problem: at
   ! check_xc's worst polarised point, tau_b = 8.49e-21 is BELOW libxc's
   ! 1e-20 floor, so libxc evaluates it at tau_b = 1e-20 while an
   ! unfloored (or `tiny(1.0_dp)`-floored) evaluation here used the true,
   ! smaller value -- amplified by t = tau^UEG/tau's sensitivity at tiny
   ! tau into vsigma_bb/vtau_b differing by up to ~70%. Flooring rho/tau/
   ! sigma at libxc's own thresholds instead of a bare "stay off literal
   ! zero" value brings every check_xc output back under 1e-11 with no
   ! loosened tolerance anywhere (see tools/xc_kernel_discovery.py, which
   ! no longer lists this kernel at all).
   !
   ! dens_tol is the SKIP threshold (a density this small contributes
   ! nothing worth computing, the fast zero path every kernel here has)
   ! and is libxc's own dens_threshold too, so the two codes skip at
   ! exactly the same point.
   real(dp), parameter :: dens_tol = 1.0e-13_dp
   real(dp), parameter :: floor_rho = dens_tol
   real(dp), parameter :: floor_tau = 1.0e-20_dp
   real(dp), parameter :: floor_sigma = dens_tol**(4.0_dp/3.0_dp)*dens_tol**(4.0_dp/3.0_dp)

   ! Table 3's "6 (Final)" column. c_x,00 = 1 - c_x by the paper's one
   ! constraint (Section 5's UEG limit for exchange); c_x is the
   ! short-range exact-exchange fraction trc_xc_functional_t's `exx`
   ! carries, not anything this kernel reads.
   real(dp), parameter :: c_x00 = 0.85_dp
   real(dp), parameter :: c_x10 = 0.259_dp
   real(dp), parameter :: c_x01 = 1.007_dp

   real(dp), parameter :: c_css00 = 0.443_dp
   real(dp), parameter :: c_css10 = -4.535_dp
   real(dp), parameter :: c_css20 = -3.39_dp
   real(dp), parameter :: c_css43 = 4.278_dp
   real(dp), parameter :: c_css04 = -1.437_dp

   real(dp), parameter :: c_cos00 = 1.0_dp
   real(dp), parameter :: c_cos10 = 1.358_dp
   real(dp), parameter :: c_cos20 = 2.924_dp
   real(dp), parameter :: c_cos60 = -1.39_dp
   real(dp), parameter :: c_cos21 = -8.812_dp
   real(dp), parameter :: c_cos61 = 9.142_dp

   ! Becke's 1986 nonlinear exchange parameter and the two correlation ones
   ! fit to Ne/He in 1997 (Section 4's text beside Eqs. 15, 24, 28), and the
   ! range-separation omega wB97X-V's own fit left untouched (Section 5).
   real(dp), parameter :: gamma_x = 0.004_dp
   real(dp), parameter :: gamma_css = 0.2_dp
   real(dp), parameter :: gamma_cos = 0.006_dp
   real(dp), parameter :: omega_sr = 0.3_dp

   ! e_x,sigma^UEG = -c_ueg * rho_sigma^(4/3) (Eq. 9).
   real(dp), parameter :: c_ueg = 1.5_dp*m_cbrt_3*m_cbrt_one_ov_pi/m_cbrt_4
   ! tau_sigma^UEG = c_tau_ueg * rho_sigma^(5/3), the per-spin Thomas-Fermi
   ! kinetic-energy density: c_tau_ueg = (3/10)(6 pi^2)^(2/3) -- the
   ! standard spin-scaled constant (e.g. SCAN's and TPSS's own tau^unif),
   ! consistent with the well-known UNPOLARISED one, (3/10)(3 pi^2)^(2/3),
   ! at rho_a = rho_b = rho/2. Checked against pyscf's WB97M_V at a w = 0
   ! point (tau_sigma set to tau_sigma^UEG exactly) to the kernel's own
   ! tolerance; the module header's finite-difference check alone cannot
   ! catch a wrong overall tau scale, only a wrong rho/tau functional form.
   real(dp), parameter :: c_tau_ueg = 0.3_dp*(6.0_dp*m_pi*m_pi)**(2.0_dp/3.0_dp)

contains

   pure subroutine wb97m_v_xc_exc_unpolar(rho, sigma, tau, eps)
      !$acc routine seq
      real(dp), intent(in) :: rho, sigma, tau
      real(dp), intent(out) :: eps
      real(dp) :: vrho, vsigma, vtau
      call wb97m_v_xc_exc_vxc_unpolar(rho, sigma, tau, eps, vrho, vsigma, vtau)
   end subroutine wb97m_v_xc_exc_unpolar

   pure subroutine wb97m_v_xc_exc_vxc_unpolar(rho, sigma, tau, eps, vrho, vsigma, vtau)
      !$acc routine seq
      real(dp), intent(in) :: rho, sigma, tau
      real(dp), intent(out) :: eps, vrho, vsigma, vtau
      real(dp) :: rho_c, sigma_c, tau_c, rho_s, sigma_s, tau_s, e
      real(dp) :: vra, vrb, vsaa, vsab, vsbb, vta, vtb
      eps = 0.0_dp
      vrho = 0.0_dp
      vsigma = 0.0_dp
      vtau = 0.0_dp
      if (rho <= dens_tol) return
      rho_c = rho
      sigma_c = max(sigma, floor_sigma)
      tau_c = max(tau, floor_tau)
      ! The symmetric point of the polarised formula: an unpolarised
      ! density has rho_a = rho_b and, since grad rho_a = grad rho_b =
      ! grad rho / 2 throughout, sigma_aa = sigma_ab = sigma_bb = sigma/4
      ! and tau_a = tau_b = tau/2 (see the module header for the reduction).
      rho_s = 0.5_dp*rho_c
      sigma_s = 0.25_dp*sigma_c
      tau_s = 0.5_dp*tau_c
      call wb97m_v_core(rho_s, rho_s, sigma_s, sigma_s, sigma_s, tau_s, tau_s, &
         e, vra, vrb, vsaa, vsab, vsbb, vta, vtb)
      eps = e/rho_c
      vrho = 0.5_dp*(vra + vrb)
      vsigma = 0.25_dp*(vsaa + vsab + vsbb)
      vtau = 0.5_dp*(vta + vtb)
   end subroutine wb97m_v_xc_exc_vxc_unpolar

   pure subroutine wb97m_v_xc_exc_polar(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, tau_a, &
      tau_b, eps)
      !$acc routine seq
      real(dp), intent(in) :: rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, tau_a, tau_b
      real(dp), intent(out) :: eps
      real(dp) :: vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b
      call wb97m_v_xc_exc_vxc_polar(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, tau_a, tau_b, &
         eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b)
   end subroutine wb97m_v_xc_exc_polar

   pure subroutine wb97m_v_xc_exc_vxc_polar(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, &
      tau_a, tau_b, eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b)
      !$acc routine seq
      real(dp), intent(in) :: rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, tau_a, tau_b
      real(dp), intent(out) :: eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b
      real(dp) :: rho_a_c, rho_b_c, sigma_aa_c, sigma_bb_c, tau_a_c, tau_b_c, e
      eps = 0.0_dp
      vrho_a = 0.0_dp
      vrho_b = 0.0_dp
      vsigma_aa = 0.0_dp
      vsigma_ab = 0.0_dp
      vsigma_bb = 0.0_dp
      vtau_a = 0.0_dp
      vtau_b = 0.0_dp
      if (rho_a + rho_b <= dens_tol) return
      ! sigma_ab never enters this functional (see the module header), so
      ! it needs no Cauchy-Schwarz clamp -- it is read nowhere below.
      ! floor_rho/floor_tau/floor_sigma are libxc's own per-spin
      ! thresholds (module header), not dens_tol's SUM check: a lopsided
      ! case (rho_a >> rho_b, both above dens_tol's sum but rho_b alone
      ! below it) must still floor rho_b at libxc's per-spin value, not
      ! at the (coarser) total-density skip threshold.
      rho_a_c = max(rho_a, floor_rho)
      rho_b_c = max(rho_b, floor_rho)
      sigma_aa_c = max(sigma_aa, floor_sigma)
      sigma_bb_c = max(sigma_bb, floor_sigma)
      tau_a_c = max(tau_a, floor_tau)
      tau_b_c = max(tau_b, floor_tau)
      call wb97m_v_core(rho_a_c, rho_b_c, sigma_aa_c, sigma_ab, sigma_bb_c, tau_a_c, tau_b_c, &
         e, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b)
      eps = e/(rho_a_c + rho_b_c)
   end subroutine wb97m_v_xc_exc_vxc_polar

   !
   ! The general two-spin formula, Eqs. 4-5 and 12-28 minus exact exchange
   ! and VV10. Inputs are assumed already floored away from zero by the
   ! callers above; `eps` returned here is rho*eps (an energy density, Eq.
   ! 8's convention), not eps itself -- both exported routines divide by
   ! their own total density afterwards.
   !
   pure subroutine wb97m_v_core(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, tau_a, tau_b, &
      eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b)
      !$acc routine seq
      real(dp), intent(in) :: rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb, tau_a, tau_b
      real(dp), intent(out) :: eps, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb, vtau_a, vtau_b
      real(dp) :: ex_a, dexa_drho, dexa_dsigma, dexa_dtau
      real(dp) :: ex_b, dexb_drho, dexb_dsigma, dexb_dtau
      real(dp) :: t_a, t_b, dta_drho, dta_dtau, dtb_drho, dtb_dtau
      real(dp) :: s2_a, s2_b, ds2a_drho, ds2a_dsigma, ds2b_drho, ds2b_dsigma
      real(dp) :: css_a, dcssa_drho, dcssa_dsigma, dcssa_dtau
      real(dp) :: css_b, dcssb_drho, dcssb_dsigma, dcssb_dtau
      real(dp) :: eps_aa, eps_bb, eps_full, e_aa, va1, vb1, e_bb, va2, vb2
      real(dp) :: e_full, va_full, vb_full, e_os
      real(dp) :: dva_os, dvb_os
      real(dp) :: t_ab, w_ab, s2_ab, u_ab, g_cos, dgcos_dw, dgcos_du
      real(dp) :: dwab_drhoa, dwab_dtaua, dwab_drhob, dwab_dtaub
      real(dp) :: duab_drhoa, duab_dsaa, duab_drhob, duab_dsbb
      real(dp) :: dwdt_ab, du_ab_ds2ab

      call exchange_spin(rho_a, sigma_aa, tau_a, ex_a, dexa_drho, dexa_dsigma, dexa_dtau, &
         t_a, dta_drho, dta_dtau, s2_a, ds2a_drho, ds2a_dsigma)
      call exchange_spin(rho_b, sigma_bb, tau_b, ex_b, dexb_drho, dexb_dsigma, dexb_dtau, &
         t_b, dtb_drho, dtb_dtau, s2_b, ds2b_drho, ds2b_dsigma)

      ! Same-spin correlation: e_c,sigma-sigma = rho_sigma*eps_c^PW92(rho_sigma, 0),
      ! the fully spin-polarised PW92 value, from pw91_lda's own polarised
      ! routine evaluated with the other spin's density at zero (Eq. 20's
      ! two subtracted integrals, read off the general two-density kernel
      ! rather than a separate zeta=1 formula).
      call pw91_lda_mod_exc_vxc_polar(rho_a, 0.0_dp, eps_aa, va1, vb1)
      e_aa = rho_a*eps_aa
      call pw91_lda_mod_exc_vxc_polar(0.0_dp, rho_b, eps_bb, va2, vb2)
      e_bb = rho_b*eps_bb
      call css_term(e_aa, va1, t_a, dta_drho, dta_dtau, s2_a, ds2a_drho, ds2a_dsigma, &
         css_a, dcssa_drho, dcssa_dsigma, dcssa_dtau)
      call css_term(e_bb, vb2, t_b, dtb_drho, dtb_dtau, s2_b, ds2b_drho, ds2b_dsigma, &
         css_b, dcssb_drho, dcssb_dsigma, dcssb_dtau)

      ! Opposite-spin correlation: e_c,ab = rho*eps_c^PW92(rho_a,rho_b) - the
      ! two same-spin pieces just computed (Eq. 20).
      call pw91_lda_mod_exc_vxc_polar(rho_a, rho_b, eps_full, va_full, vb_full)
      e_full = (rho_a + rho_b)*eps_full
      e_os = e_full - e_aa - e_bb
      dva_os = va_full - va1
      dvb_os = vb_full - vb2

      t_ab = 0.5_dp*(t_a + t_b)
      w_ab = (t_ab - 1.0_dp)/(t_ab + 1.0_dp)
      s2_ab = 0.5_dp*(s2_a + s2_b)
      u_ab = gamma_cos*s2_ab/(1.0_dp + gamma_cos*s2_ab)
      g_cos = c_cos00 + c_cos10*w_ab + c_cos20*w_ab**2 + c_cos60*w_ab**6 + &
         c_cos21*w_ab**2*u_ab + c_cos61*w_ab**6*u_ab
      dgcos_dw = c_cos10 + 2.0_dp*c_cos20*w_ab + 6.0_dp*c_cos60*w_ab**5 + &
         2.0_dp*c_cos21*w_ab*u_ab + 6.0_dp*c_cos61*w_ab**5*u_ab
      dgcos_du = c_cos21*w_ab**2 + c_cos61*w_ab**6

      dwdt_ab = 2.0_dp/(t_ab + 1.0_dp)**2
      dwab_drhoa = dwdt_ab*0.5_dp*dta_drho
      dwab_dtaua = dwdt_ab*0.5_dp*dta_dtau
      dwab_drhob = dwdt_ab*0.5_dp*dtb_drho
      dwab_dtaub = dwdt_ab*0.5_dp*dtb_dtau

      du_ab_ds2ab = gamma_cos/(1.0_dp + gamma_cos*s2_ab)**2
      duab_drhoa = du_ab_ds2ab*0.5_dp*ds2a_drho
      duab_dsaa = du_ab_ds2ab*0.5_dp*ds2a_dsigma
      duab_drhob = du_ab_ds2ab*0.5_dp*ds2b_drho
      duab_dsbb = du_ab_ds2ab*0.5_dp*ds2b_dsigma

      eps = ex_a + ex_b + css_a + css_b + e_os*g_cos

      vrho_a = dexa_drho + dcssa_drho + dva_os*g_cos + e_os*(dgcos_dw*dwab_drhoa + &
         dgcos_du*duab_drhoa)
      vrho_b = dexb_drho + dcssb_drho + dvb_os*g_cos + e_os*(dgcos_dw*dwab_drhob + &
         dgcos_du*duab_drhob)
      vsigma_aa = dexa_dsigma + dcssa_dsigma + e_os*dgcos_du*duab_dsaa
      vsigma_bb = dexb_dsigma + dcssb_dsigma + e_os*dgcos_du*duab_dsbb
      vsigma_ab = 0.0_dp
      vtau_a = dexa_dtau + dcssa_dtau + e_os*dgcos_dw*dwab_dtaua
      vtau_b = dexb_dtau + dcssb_dtau + e_os*dgcos_dw*dwab_dtaub
   end subroutine wb97m_v_core

   !
   ! One spin's SR-mGGA exchange energy density (Eqs. 9-15) and its rho/
   ! sigma/tau derivatives, plus t_sigma, s2_sigma and their rho/sigma/tau
   ! derivatives for the caller to reuse building w_ab, s2_ab (css/cos
   ! share the same t_sigma and s2_sigma terco's exchange needs; only the
   ! GGA ingredient's gamma differs per use).
   !
   pure subroutine exchange_spin(rho, sigma, tau, ex, dex_drho, dex_dsigma, dex_dtau, &
      t, dt_drho, dt_dtau, s2, ds2_drho, ds2_dsigma)
      !$acc routine seq
      real(dp), intent(in) :: rho, sigma, tau
      real(dp), intent(out) :: ex, dex_drho, dex_dsigma, dex_dtau
      real(dp), intent(out) :: t, dt_drho, dt_dtau, s2, ds2_drho, ds2_dsigma
      real(dp) :: rho13, rho43, kf, a, fx, dfx_da, e_ueg, de_ueg, h, dh_drho
      real(dp) :: tau_ueg, w, dw_dt, dw_drho, dw_dtau, u, du_ds2, du_drho, du_dsigma
      real(dp) :: gx, dgx_dw, dgx_du

      rho13 = xc_cbrt(rho)
      rho43 = rho*rho13  ! rho^(4/3)
      kf = xc_cbrt(6.0_dp*m_pi*m_pi*rho)
      a = omega_sr/kf
      call sr_attenuation(a, fx, dfx_da)
      e_ueg = -c_ueg*rho43
      de_ueg = -c_ueg*(4.0_dp/3.0_dp)*rho13
      h = e_ueg*fx
      dh_drho = de_ueg*fx + e_ueg*dfx_da*(-a/(3.0_dp*rho))

      tau_ueg = c_tau_ueg*rho*rho13*rho13  ! rho^(5/3) = rho * rho^(2/3)
      t = tau_ueg/tau
      dt_drho = (5.0_dp/3.0_dp)*tau_ueg/(rho*tau)
      dt_dtau = -t/tau
      w = (t - 1.0_dp)/(t + 1.0_dp)
      dw_dt = 2.0_dp/(t + 1.0_dp)**2
      dw_drho = dw_dt*dt_drho
      dw_dtau = dw_dt*dt_dtau

      s2 = sigma/(rho43*rho43)  ! rho^(8/3) = (rho^(4/3))^2
      ds2_drho = -(8.0_dp/3.0_dp)*s2/rho
      ds2_dsigma = 1.0_dp/(rho43*rho43)
      u = gamma_x*s2/(1.0_dp + gamma_x*s2)
      du_ds2 = gamma_x/(1.0_dp + gamma_x*s2)**2
      du_drho = du_ds2*ds2_drho
      du_dsigma = du_ds2*ds2_dsigma

      gx = c_x00 + c_x10*w + c_x01*u
      dgx_dw = c_x10
      dgx_du = c_x01

      ex = h*gx
      dex_drho = dh_drho*gx + h*(dgx_dw*dw_drho + dgx_du*du_drho)
      dex_dsigma = h*dgx_du*du_dsigma
      dex_dtau = h*dgx_dw*dw_dtau
   end subroutine exchange_spin

   !
   ! F_x(a) = 1 - (2/3) a [2 sqrt(pi) erf(1/a) - 3a + a^3 + (2a - a^3) exp(-1/a^2)],
   ! Eq. 11, and its derivative dF_x/da (hand-differentiated and checked
   ! against central finite differences; see the module header) -- the
   ! SMALL-a form, exact in double precision up to a ~ 4 (checked against
   ! 50-digit mpmath: 3.3e-12 relative at a = 4, still only 9e-12 at a = 5).
   !
   ! Past that, a^3 and (2a - a^3)*exp(-1/a^2) are two O(a^3) terms whose
   ! difference is F_x itself, O(1/a^2): at a = 450 (check_xc's worst
   ! point, rho_b = 5e-12) the cancellation costs about eleven decimal
   ! digits, which is exactly what turned up as a 70% relative miss in
   ! this kernel's vtau/vsigma against libxc -- not a formula bug (an
   ! arbitrary-precision evaluation of this SAME closed form reproduced
   ! libxc to 1e-10), a genuine double-precision instability in it.
   !
   ! `sr_attenuation_large_a` is the fix: F_x and dF_x/da as a power
   ! series in x = 1/a (small when a is large), from this closed form's
   ! own Taylor expansion around x = 0, carried to x^14/x^15 -- checked
   ! against mpmath from a = 3 to a = 100 and at a = 450, better than
   ! 1e-14 relative everywhere a > 4 is reached from here. The crossover
   ! at a = 4.0 is where both forms already agree to ~1e-12, so switching
   ! there introduces no discontinuity worth worrying about.
   !
   pure subroutine sr_attenuation(a, fx, dfx_da)
      !$acc routine seq
      real(dp), intent(in) :: a
      real(dp), intent(out) :: fx, dfx_da
      real(dp) :: inv_a, erf_term, exp_term, g, gprime

      if (a > 4.0_dp) then
         call sr_attenuation_large_a(a, fx, dfx_da)
         return
      end if
      inv_a = 1.0_dp/a
      erf_term = xc_erf(inv_a)
      exp_term = exp(-inv_a*inv_a)
      g = 2.0_dp*sqrt(m_pi)*erf_term - 3.0_dp*a + a**3 + (2.0_dp*a - a**3)*exp_term
      fx = 1.0_dp - (2.0_dp/3.0_dp)*a*g
      gprime = 2.0_dp*sqrt(m_pi)*erf_term - 6.0_dp*a + 4.0_dp*a**3 + &
         (2.0_dp*a - 4.0_dp*a**3)*exp_term
      dfx_da = -(2.0_dp/3.0_dp)*gprime
   end subroutine sr_attenuation

   !
   ! F_x(1/x) = x^2/9 - x^4/60 + x^6/420 - x^8/3240 + x^10/27720
   !            - x^12/262080 + x^14/2721600
   ! dF_x/da(1/x) = -2x^3/9 + x^5/15 - x^7/70 + x^9/405 - x^11/2772
   !                + x^13/21840 - x^15/194400
   ! (x = 1/a), from sympy's `series` of this module's own closed form
   ! around x = 0 -- not a different functional, the same F_x(a) written
   ! in the variable where it is well behaved.
   !
   pure subroutine sr_attenuation_large_a(a, fx, dfx_da)
      !$acc routine seq
      real(dp), intent(in) :: a
      real(dp), intent(out) :: fx, dfx_da
      real(dp) :: x, x2

      x = 1.0_dp/a
      x2 = x*x
      fx = x2*(1.0_dp/9.0_dp + x2*(-1.0_dp/60.0_dp + x2*(1.0_dp/420.0_dp + &
         x2*(-1.0_dp/3240.0_dp + x2*(1.0_dp/27720.0_dp + x2*(-1.0_dp/262080.0_dp + &
         x2/2721600.0_dp))))))
      dfx_da = x2*x*(-2.0_dp/9.0_dp + x2*(1.0_dp/15.0_dp + x2*(-1.0_dp/70.0_dp + &
         x2*(1.0_dp/405.0_dp + x2*(-1.0_dp/2772.0_dp + x2*(1.0_dp/21840.0_dp - &
         x2/194400.0_dp))))))
   end subroutine sr_attenuation_large_a

   !
   ! One spin's same-spin correlation energy density and derivatives
   ! (Eqs. 21-24): e_css*g_css(w,u), with e_css = rho_sigma*eps_c^PW92
   ! (rho_sigma,0) and its own rho-derivative (dec_drho) handed in by the
   ! caller (from pw91_lda_mod_exc_vxc_polar at the other spin's density
   ! zero), and w, s2 (hence t, its own derivatives) shared with exchange.
   !
   pure subroutine css_term(e_css, dec_drho, t, dt_drho, dt_dtau, s2, ds2_drho, ds2_dsigma, &
      e, de_drho, de_dsigma, de_dtau)
      !$acc routine seq
      real(dp), intent(in) :: e_css, dec_drho
      real(dp), intent(in) :: t, dt_drho, dt_dtau, s2, ds2_drho, ds2_dsigma
      real(dp), intent(out) :: e, de_drho, de_dsigma, de_dtau
      real(dp) :: w, dw_dt, dw_drho, dw_dtau, u, du_ds2, du_drho, du_dsigma
      real(dp) :: g, dg_dw, dg_du

      w = (t - 1.0_dp)/(t + 1.0_dp)
      dw_dt = 2.0_dp/(t + 1.0_dp)**2
      dw_drho = dw_dt*dt_drho
      dw_dtau = dw_dt*dt_dtau

      u = gamma_css*s2/(1.0_dp + gamma_css*s2)
      du_ds2 = gamma_css/(1.0_dp + gamma_css*s2)**2
      du_drho = du_ds2*ds2_drho
      du_dsigma = du_ds2*ds2_dsigma

      g = c_css00 + c_css10*w + c_css20*w**2 + c_css43*w**4*u**3 + c_css04*u**4
      dg_dw = c_css10 + 2.0_dp*c_css20*w + 4.0_dp*c_css43*w**3*u**3
      dg_du = 3.0_dp*c_css43*w**4*u**2 + 4.0_dp*c_css04*u**3

      e = e_css*g
      de_drho = dec_drho*g + e_css*(dg_dw*dw_drho + dg_du*du_drho)
      de_dsigma = e_css*dg_du*du_dsigma
      de_dtau = e_css*dg_dw*dw_dtau
   end subroutine css_term

end module trc_xc_wb97m_v_xc
