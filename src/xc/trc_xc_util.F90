!
! What the translated exchange-correlation kernels need that Fortran does not
! give them: the handful of constants and the few helper functions ExchCXX
! keeps in `impl/builtin/util.hpp` and `constants.hpp`.
!
! The kernels themselves are GENERATED, by tools/exchcxx_to_fortran.py, from
! ExchCXX's Maple-emitted C++ (MIT, Lawrence Berkeley National Laboratory and
! Microsoft). This file is written by hand and is the only piece of the XC
! functional layer that is.
!
! Everything here is `pure` and scalar, and marked `!$acc routine seq`, because
! it is called from inside `do concurrent` over grid points. Nothing here may
! touch a module variable or do I/O.
!
! The cube root and the scaled complementary error function come from
! `bit_repro.f90`, Jorge's bit-reproducible transcendentals (vendored from
! nci/learning_tools/ci_enabled/bitwise_adventures with one change: an
! `!$acc routine seq` beside every `!$omp declare target`, because terco
! offloads through OpenACC and nvfortran wants the routine information in
! that dialect. A directive is a comment to the arithmetic, so the golden
! hash still describes the values). Their contract is no FMA contraction and no
! flush-to-zero; terco's `-fast` GPU build breaks the first half of that, so
! CPU-vs-GPU bit identity is only as good as the flags until that is settled.
!
module trc_xc_util
   use trc_boys, only: dp
   use bit_repro, only: erfc_reprod, exp_reprod, cuberoot
   implicit none
   private

   public :: xc_cbrt, xc_square, pow_3_2, pow_1_4, xc_erfcx
   public :: enforce_polar_sigma_constraints, enforce_fermi_hole_curvature

   ! ExchCXX constants.hpp, digit for digit. The long-double literals there
   ! are rounded to double at the point of use in C++; here they are rounded
   ! once, at compile time, to the same double.
   real(dp), parameter, public :: m_cbrt_2    = 1.259921049894873164767210607278228350570_dp
   real(dp), parameter, public :: m_cbrt_3    = 1.442249570307408382321638310780109588392_dp
   real(dp), parameter, public :: m_cbrt_4    = 1.587401051968199474751705639272308260391_dp
   real(dp), parameter, public :: m_cbrt_6    = 1.817120592832139658891211756327260502428_dp
   real(dp), parameter, public :: m_cbrt_pi   = 1.464591887561523263020142527263790391739_dp
   real(dp), parameter, public :: m_pi        = 3.14159265358979323846e+00_dp
   real(dp), parameter, public :: m_one_ov_pi = 3.18309886183790691216e-01_dp
   real(dp), parameter, public :: m_pi_sq     = 9.869604401089357992305e+00_dp
   real(dp), parameter, public :: m_cbrt_one_ov_pi = 6.82784063255295503581e-01_dp
   real(dp), parameter, public :: m_cbrt_pi_sq     = 2.145029397111025470934e+00_dp
   real(dp), parameter, public :: x2s        = 0.1282782438530421943003109254455883701296_dp
   real(dp), parameter, public :: x_factor_c = 0.9305257363491000250020102180716672510262_dp

contains

   !
   ! Real cube root. Fortran has no cbrt intrinsic, and a pow is a few ulp
   ! off, undefined below zero, and not the same bits on every platform.
   ! `cuberoot` from bit_repro is built from correctly-rounded operations
   ! in a frozen order, so the CPU and the GPU agree bit for bit and CI on
   ! that repository proves it across seven compilers.
   !
   pure function xc_cbrt(x) result(y)
      !$acc routine seq
      real(dp), intent(in) :: x
      real(dp) :: y
      y = cuberoot(x)
   end function xc_cbrt

   pure function xc_square(x) result(y)
      !$acc routine seq
      real(dp), intent(in) :: x
      real(dp) :: y
      y = x*x
   end function xc_square

   pure function pow_3_2(x) result(y)
      !$acc routine seq
      real(dp), intent(in) :: x
      real(dp) :: y
      real(dp) :: s
      s = sqrt(x)
      y = s*s*s
   end function pow_3_2

   pure function pow_1_4(x) result(y)
      !$acc routine seq
      real(dp), intent(in) :: x
      real(dp) :: y
      y = sqrt(sqrt(x))
   end function pow_1_4

   !
   ! exp(x^2) erfc(x). Fortran 2008 has this as `erfc_scaled`, but nvfortran
   ! has no device implementation of it (NVFORTRAN-S-1062 inside a compute
   ! region), and the intrinsic erfc is not the same bits on the CPU and the
   ! GPU anyway. Both factors come from bit_repro instead. Overflow for
   ! large positive x is not a concern: every use is inside a range-separated
   ! exchange kernel with a bounded argument.
   !
   pure function xc_erfcx(x) result(y)
      !$acc routine seq
      real(dp), intent(in) :: x
      real(dp) :: y
      y = exp_reprod(x*x)*erfc_reprod(x)
   end function xc_erfcx

   ! |sigma_ab| may not exceed the mean of the two like-spin gradients.
   pure function enforce_polar_sigma_constraints(sigma_aa, sigma_ab, sigma_bb) result(s)
      !$acc routine seq
      real(dp), intent(in) :: sigma_aa, sigma_ab, sigma_bb
      real(dp) :: s
      real(dp) :: s_ave
      s_ave = 0.5_dp*(sigma_aa + sigma_bb)
      s = sigma_ab
      if (s < -s_ave) s = -s_ave
      if (s > s_ave) s = s_ave
   end function enforce_polar_sigma_constraints

   ! sigma <= 8 rho tau, the meta-GGA bound.
   pure function enforce_fermi_hole_curvature(sigma, rho, tau) result(s)
      !$acc routine seq
      real(dp), intent(in) :: sigma, rho, tau
      real(dp) :: s
      s = min(sigma, 8.0_dp*rho*tau)
   end function enforce_fermi_hole_curvature

end module trc_xc_util
