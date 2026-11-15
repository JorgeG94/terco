module trc_ecp_quadrature
   !!
   !! Numerical primitives for the effective-core-potential integrals: the
   !! adaptive radial grid and the scaled modified spherical Bessel function.
   !!
   !! Ported from libfint's src/cint_ecp_num.f90 (Apache-2.0, itself translated
   !! from PySCF's pyscf/lib/gto/nr_ecp.c -- libcint proper has no ECP code).
   !! Nothing here depends on a basis or a shell, which is what lets it sit in
   !! `core`: it is checked bit for bit against libfint's own values before any
   !! integral is assembled on top of it.
   !!
   !! THE 2047-POINT TABLE IS NOT CARRIED. PySCF's C hard-codes a 2047-point
   !! Gauss-Chebyshev grid and two Bessel-value tables (about 4500 lines); every
   !! one of them is a cache of `ecp_quadrature_gauss_chebyshev`/`ecp_quadrature_sph_ine` below, so the
   !! computation is kept and the cache is not.
   !!
   use pic_types, only: dp
   implicit none
   private

   public :: ECP_LMAX, ECP_LEVEL0, ECP_LEVEL_MAX, ECP_NRS
   public :: ECP_SIM_ZERO, ECP_EXPCUTOFF, ECP_CUTOFF
   public :: ecp_quadrature_gauss_chebyshev, ecp_quadrature_sph_ine
   public :: ecp_quadrature_factorial, ecp_quadrature_factorial2, ecp_quadrature_binom

   integer, parameter :: ECP_LMAX = 5
   !! Highest angular momentum an ECP CHANNEL itself may carry (libfint's
   !! own ceiling; def2's heaviest local channels reach l=3-4).

   real(dp), parameter :: ECP_SIM_ZERO = 1.0e-50_dp
   !! A coefficient below this is treated as absent.
   real(dp), parameter :: ECP_EXPCUTOFF = 39.0_dp
   !! exp(-x) below this is zero to double precision (x=39 gives ~1e-17).
   real(dp), parameter :: ECP_CUTOFF = 460.0_dp
   !! exp(x) above this overflows (x=460 gives ~1e200).

   ! Radial grid refinement levels: 2**(level+1)-1 points at `level`.
   !
   ! Matches PySCF's own nr_ecp.c ceiling (2047 points). An earlier cut-down
   ! version of this port stopped at 9 (511 points) to bound the per-level
   ! buffers this feeds (`trc_ecp_radial.F90`'s radial factors), which are
   ! PER-THREAD scratch in a `do concurrent` item, fixed-size at the worst
   ! case. Validation against pyscf (HI/def2-SVP, I...H at 25 bohr, Xe,
   ! CsH, SnH4) found shell pairs whose value was already stable well
   ! before level 9 but that need one or two further levels to show TWO
   ! consecutive matches at the 1e-12 relative tolerance the convergence
   ! test asks for -- `trc_ecp.F90`'s own ecp_type1_item/ecp_type2_item
   ! separately screen out the pairs that are instead stuck because their
   ! actual contribution has decayed to numerically negligible, which no
   ! amount of refinement would fix. With both effects accounted for, 11
   ! (2047 points) is where every case validated so far converges; each
   ! NRS-sized buffer this feeds is correspondingly about 4x the size the
   ! 511-point cut-down needed (several past 100 KB per thread) -- a real
   ! cost to revisit once this runs on the GPU and gets measured there.
   integer, parameter :: ECP_LEVEL0 = 5
   integer, parameter :: ECP_LEVEL_MAX = 11
   integer, parameter :: ECP_NRS = 2**ECP_LEVEL_MAX - 1
   !! Points in the finest grid: 2047, not 2048 -- see `ecp_quadrature_gauss_chebyshev`.

   ! Highest n tabulated in FACTORIAL/FACTORIAL2 below.
   integer, parameter :: FACTORIAL_NMAX = 23
   integer, parameter :: FACTORIAL2_NMAX = 39

   real(dp), parameter :: FACTORIAL(0:FACTORIAL_NMAX) = [ &
      1.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, 24.0_dp, &
      1.2e+2_dp, 7.2e+2_dp, 5.04e+3_dp, 4.032e+4_dp, 3.6288e+5_dp, &
      3.6288e+6_dp, 3.99168e+7_dp, 4.790016e+8_dp, 6.2270208e+9_dp, 8.71782912e+10_dp, &
      1.307674368e+12_dp, 2.0922789888e+13_dp, 3.55687428096e+14_dp, &
      6.402373705728e+15_dp, 1.21645100408832e+17_dp, &
      2.43290200817664e+18_dp, 5.109094217170944e+19_dp, &
      1.1240007277776077e+21_dp, 2.5852016738884978e+22_dp]

   real(dp), parameter :: FACTORIAL2(0:FACTORIAL2_NMAX) = [ &
      1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 8.0_dp, &
      15.0_dp, 48.0_dp, 105.0_dp, 384.0_dp, 945.0_dp, &
      3840.0_dp, 10395.0_dp, 46080.0_dp, 135135.0_dp, 645120.0_dp, &
      2027025.0_dp, 10321920.0_dp, 34459425.0_dp, 185794560.0_dp, 654729075.0_dp, &
      3715891200.0_dp, 13749310575.0_dp, 81749606400.0_dp, 316234143225.0_dp, &
      1961990553600.0_dp, &
      7905853580625.0_dp, 51011754393600.0_dp, 213458046676875.0_dp, &
      1428329123020800.0_dp, 6190283353629376.0_dp, &
      42849873690624000.0_dp, 1.9189878396251069e+17_dp, &
      1.371195958099968e+18_dp, 6.3326598707628524e+18_dp, &
      4.6620662575398912e+19_dp, 2.2164309547669976e+20_dp, &
      1.6783438527143608e+21_dp, 8.2007945326378929e+21_dp, &
      6.3777066403145712e+22_dp, 3.1983098677287775e+23_dp]

contains

   pure real(dp) function ecp_quadrature_factorial(n) result(f)
      !! n!, and 1 outside the tabulated range (same convention libfint uses).
      !$acc routine seq
      integer, intent(in) :: n
      if (n < 0 .or. n > ubound(FACTORIAL, 1)) then
         f = 1.0_dp
      else
         f = FACTORIAL(n)
      end if
   end function ecp_quadrature_factorial

   pure real(dp) function ecp_quadrature_factorial2(n) result(f)
      !! n!!, with the convention that a negative argument gives 1 -- load
      !! bearing for the type-2 angular factors, which index this at l-1.
      !$acc routine seq
      integer, intent(in) :: n
      if (n < 0 .or. n > ubound(FACTORIAL2, 1)) then
         f = 1.0_dp
      else
         f = FACTORIAL2(n)
      end if
   end function ecp_quadrature_factorial2

   pure real(dp) function ecp_quadrature_binom(n, m) result(b)
      !! Binomial coefficient n choose m.
      !$acc routine seq
      integer, intent(in) :: n, m
      if (m < 0 .or. m > n) then
         b = 0.0_dp
      else
         b = ecp_quadrature_factorial(n)/(ecp_quadrature_factorial(m)*ecp_quadrature_factorial(n - m))
      end if
   end function ecp_quadrature_binom

   pure subroutine ecp_quadrature_gauss_chebyshev(rs, ws, n)
      !! Gauss-Chebyshev radial grid on (0,1], mapped by Perez-Jorda, Becke
      !! and San-Fabian's r = 1 - log2(1+x) so an ECP's short-ranged
      !! integrand gets points where it has support, with no cutoff
      !! parameter.
      !$acc routine seq
      integer, intent(in) :: n
      real(dp), intent(out) :: rs(0:n - 1), ws(0:n - 1)

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp), parameter :: LOG2E = 1.0_dp/log(2.0_dp)
      real(dp) :: step, fac, xinc, x1, x2, x3, x4, xi
      integer :: i

      step = 1.0_dp/real(n + 1, dp)
      fac = 16.0_dp*step/3.0_dp
      xinc = PI*step
      x1 = 0.0_dp
      do i = 0, n - 1
         x1 = x1 + xinc
         x2 = sin(x1)
         x3 = sin(x1*2.0_dp)
         x4 = x2*x2
         xi = real(n - i*2 - 1, dp)*step &
              + (1.0_dp/PI)*(1.0_dp + (2.0_dp/3.0_dp)*x4)*x3
         rs(i) = 1.0_dp - log(1.0_dp + xi)*LOG2E
         ws(i) = fac*x4*x4*LOG2E/(1.0_dp + xi)
      end do
   end subroutine ecp_quadrature_gauss_chebyshev

   pure subroutine ecp_quadrature_sph_ine(out, order, z)
      !! out(l) = i_l(z) exp(-z), l = 0..order: the scaled modified
      !! spherical Bessel function the ECP radial integrand needs (i_l
      !! alone overflows for large z and underflows for small). Three
      !! regimes, boundaries as in PySCF's C: a leading Taylor term below
      !! 1e-7, the terminating asymptotic sum above 16, the ascending
      !! series (bounded iteration; never remotely reached in practice,
      !! kept finite for a GPU thread) in between.
      !$acc routine seq
      integer, intent(in) :: order
      real(dp), intent(in) :: z
      real(dp), intent(out) :: out(0:order)

      integer, parameter :: KMAX = 200
      real(dp) :: z2, ti, s, next, t0
      integer :: i, k

      if (z < 1.0e-7_dp) then
         out(0) = 1.0_dp - z
         do i = 1, order
            out(i) = out(i - 1)*z/real(i*2 + 1, dp)
         end do

      else if (z > 16.0_dp) then
         z2 = -0.5_dp/z
         do i = 0, order
            ti = 0.5_dp/z
            s = ti
            do k = 1, i
               ti = ti*z2
               s = s + ti*ecp_quadrature_factorial(i + k) &
                   /(ecp_quadrature_factorial(k)*ecp_quadrature_factorial(i - k))
            end do
            out(i) = s
         end do

      else
         z2 = 0.5_dp*z*z
         t0 = exp(-z)
         do i = 0, order
            ti = t0
            s = ti
            do k = 1, KMAX
               ! Grouped as `ti*(z2/denom)`, not `ti*z2/denom`: the rounding
               ! differs in the last ULP and this grouping matches PySCF's C.
               ti = ti*(z2/real(k*(k*2 + i*2 + 1), dp))
               next = s + ti
               if (next == s) exit
               s = next
            end do
            t0 = t0*(z/real(i*2 + 3, dp))
            out(i) = s
         end do
      end if
   end subroutine ecp_quadrature_sph_ine

end module trc_ecp_quadrature
