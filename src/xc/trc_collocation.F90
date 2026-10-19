!
! Basis functions on grid points: the value and gradient of every Cartesian
! component of one shell at one point.
!
! This is the piece libcint does not have -- it computes integrals, not
! values -- and it is the piece every XC quantity is built from. Which makes
! its conventions the whole difficulty: the values here must be EXACTLY the
! functions terco's overlap integrals are over, or the density is assembled
! in one basis and the Fock matrix in another and nothing fails loudly.
!
! Two conventions, both inherited from `trc_basis_t`:
!
!   * the coefficients in `sh_c` already carry libcint's primitive
!     normalisation and its `common_fac_sp`, so NOTHING is normalised here.
!     A component is x^i y^j z^k times sum_p c_p exp(-a_p r^2), full stop;
!   * the components of a shell come out in the BUILD's Cartesian order
!     (TRC_CART_ORDER; lexicographic by default -- x descending outermost,
!     y descending inside it), which is the order every AO-indexed quantity
!     the integrals fill uses too. The table below (scripts/gen_cartcomp.py)
!     is generated for the same order and guards against a mismatch.
!
! The gradient is the same expression differentiated:
!
!     d/dx [ x^i y^j z^k R(r) ] = i x^(i-1) y^j z^k R  +  x^(i+1) y^j z^k R'
!
! with R' = sum_p (-2 a_p) c_p exp(-a_p r^2), so the radial part costs one
! extra multiply per primitive and nothing else.
!
module trc_collocation
   use pic_types, only: dp
   implicit none
   private

   public :: shell_collocate, shell_collocate2, NCART_MAX, XC_LMAX, NHESS

   ! The ceiling is generous on purpose: the loop is generic in l and a d
   ! shell costs no more than it would with a generated kernel, since the
   ! exponentials dominate.
   integer, parameter :: XC_LMAX = 4
   integer, parameter :: NCART_MAX = (XC_LMAX + 1)*(XC_LMAX + 2)/2
   integer, parameter :: NHESS = 6   !! independent Cartesian second derivatives: xx, xy, xz, yy, yz, zz
   real(dp), parameter :: COLLOCATE_EXP_MAX = 34.5_dp   !! exp(-34.5) = 1e-15: a primitive past this is skipped

   ! CC_NX/CC_NY/CC_NZ/CC_OFF/CC_NCART: the flat (ix,iy,iz) table
   ! shell_collocate reads below, in the build's Cartesian order
   ! (TRC_CART_ORDER). OUTPUT layout -- see scripts/gen_cartcomp.py, whose
   ! header explains why this is a table and not a `select case`.
#include "trc_collocation_table.inc"

contains

   pure subroutine shell_collocate(l, np, e, c, d, chi, gchi)
      !$acc routine seq
      integer, intent(in) :: l, np
      real(dp), intent(in) :: e(np), c(np)
      real(dp), intent(in) :: d(3)   !! r - A, the point relative to the shell centre
      real(dp), intent(out) :: chi(NCART_MAX), gchi(3, NCART_MAX)
      real(dp) :: r2, rad, drad, ex, ang
      real(dp) :: xp(0:XC_LMAX + 1), yp(0:XC_LMAX + 1), zp(0:XC_LMAX + 1)
      integer :: p, ix, iy, iz, m

      r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
      rad = 0.0_dp
      drad = 0.0_dp
      do p = 1, np
         ! The exponential is the cost of this kernel; a primitive whose
         ! argument is past COLLOCATE_EXP_MAX is below 1e-15 of its
         ! coefficient and is not evaluated (Stocks and Barca, Algorithm 2).
         ! The batch cutoff already keeps every shell here within 1e-10 of
         ! something, so most primitives pass; the diffuse tail of a
         ! contraction is what this skips.
         if (e(p)*r2 > COLLOCATE_EXP_MAX) cycle
         ex = c(p)*exp(-e(p)*r2)
         rad = rad + ex
         drad = drad - 2.0_dp*e(p)*ex
      end do

      xp(0) = 1.0_dp
      yp(0) = 1.0_dp
      zp(0) = 1.0_dp
      do m = 1, l + 1
         xp(m) = xp(m - 1)*d(1)
         yp(m) = yp(m - 1)*d(2)
         zp(m) = zp(m - 1)*d(3)
      end do

      ! (ix,iy,iz) comes from the flat table above rather than from nested
      ! loop counters now, so this loop follows the build's AO output order
      ! (TRC_CART_ORDER) instead of always producing lexicographic's; the
      ! per-point arithmetic itself -- one lookup and the same handful of
      ! multiplies -- is what the nested loop already cost.
      do m = 1, CC_OFF(l + 1) - CC_OFF(l)
         ix = CC_NX(CC_OFF(l) + m)
         iy = CC_NY(CC_OFF(l) + m)
         iz = CC_NZ(CC_OFF(l) + m)
         ang = xp(ix)*yp(iy)*zp(iz)
         chi(m) = ang*rad
         gchi(1, m) = xp(ix + 1)*yp(iy)*zp(iz)*drad
         gchi(2, m) = xp(ix)*yp(iy + 1)*zp(iz)*drad
         gchi(3, m) = xp(ix)*yp(iy)*zp(iz + 1)*drad
         if (ix > 0) gchi(1, m) = gchi(1, m) + real(ix, dp)*xp(ix - 1)*yp(iy)*zp(iz)*rad
         if (iy > 0) gchi(2, m) = gchi(2, m) + real(iy, dp)*xp(ix)*yp(iy - 1)*zp(iz)*rad
         if (iz > 0) gchi(3, m) = gchi(3, m) + real(iz, dp)*xp(ix)*yp(iy)*zp(iz - 1)*rad
      end do
   end subroutine shell_collocate

   !
   ! As `shell_collocate`, plus the Hessian: the GGA exchange-correlation
   ! GRADIENT's second-derivative term (trc_xc_grad.F90's `hchi`) is the only
   ! reader. `shell_collocate` itself is untouched -- every energy-path
   ! caller keeps paying for the value and gradient only -- so this is a
   ! second entry point rather than an optional argument on the first; an
   ! optional dummy on a `!$acc routine seq` procedure is the kind of thing
   ! that works on one compiler and not the next.
   !
   ! d^2/dx_a dx_b [ f(x,y,z) R(r) ] = f_ab R + drad (f_a x_b + f_b x_a)
   !                                   + f delta_ab drad + f x_a x_b ddrad
   !
   ! with f the monomial x^i y^j z^k, R the radial part, drad = dR/dx / x
   ! (shell_collocate's own convention) and ddrad = d(drad)/dx / x, i.e. the
   ! second radial derivative in the same "divide by the coordinate" form
   ! that keeps every piece a plain multiply, no division anywhere.
   !
   pure subroutine shell_collocate2(l, np, e, c, d, chi, gchi, hchi)
      !$acc routine seq
      integer, intent(in) :: l, np
      real(dp), intent(in) :: e(np), c(np)
      real(dp), intent(in) :: d(3)   !! r - A, the point relative to the shell centre
      real(dp), intent(out) :: chi(NCART_MAX), gchi(3, NCART_MAX), hchi(NHESS, NCART_MAX)
      !! hchi components, in order: xx, xy, xz, yy, yz, zz
      real(dp) :: r2, rad, drad, ddrad, ex, ang
      real(dp) :: xp(0:XC_LMAX + 1), yp(0:XC_LMAX + 1), zp(0:XC_LMAX + 1)
      real(dp) :: fx, fy, fz, fxx, fyy, fzz, fxy, fxz, fyz
      integer :: p, ix, iy, iz, m

      r2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
      rad = 0.0_dp
      drad = 0.0_dp
      ddrad = 0.0_dp
      do p = 1, np
         if (e(p)*r2 > COLLOCATE_EXP_MAX) cycle
         ex = c(p)*exp(-e(p)*r2)
         rad = rad + ex
         drad = drad - 2.0_dp*e(p)*ex
         ddrad = ddrad + 4.0_dp*e(p)*e(p)*ex
      end do

      xp(0) = 1.0_dp
      yp(0) = 1.0_dp
      zp(0) = 1.0_dp
      do m = 1, l + 1
         xp(m) = xp(m - 1)*d(1)
         yp(m) = yp(m - 1)*d(2)
         zp(m) = zp(m - 1)*d(3)
      end do

      do m = 1, CC_OFF(l + 1) - CC_OFF(l)
         ix = CC_NX(CC_OFF(l) + m)
         iy = CC_NY(CC_OFF(l) + m)
         iz = CC_NZ(CC_OFF(l) + m)
         ang = xp(ix)*yp(iy)*zp(iz)
         chi(m) = ang*rad

         fx = 0.0_dp
         fy = 0.0_dp
         fz = 0.0_dp
         if (ix > 0) fx = real(ix, dp)*xp(ix - 1)*yp(iy)*zp(iz)
         if (iy > 0) fy = real(iy, dp)*xp(ix)*yp(iy - 1)*zp(iz)
         if (iz > 0) fz = real(iz, dp)*xp(ix)*yp(iy)*zp(iz - 1)

         gchi(1, m) = fx*rad + ang*d(1)*drad
         gchi(2, m) = fy*rad + ang*d(2)*drad
         gchi(3, m) = fz*rad + ang*d(3)*drad

         fxx = 0.0_dp
         fyy = 0.0_dp
         fzz = 0.0_dp
         fxy = 0.0_dp
         fxz = 0.0_dp
         fyz = 0.0_dp
         if (ix > 1) fxx = real(ix*(ix - 1), dp)*xp(ix - 2)*yp(iy)*zp(iz)
         if (iy > 1) fyy = real(iy*(iy - 1), dp)*xp(ix)*yp(iy - 2)*zp(iz)
         if (iz > 1) fzz = real(iz*(iz - 1), dp)*xp(ix)*yp(iy)*zp(iz - 2)
         if (ix > 0 .and. iy > 0) fxy = real(ix*iy, dp)*xp(ix - 1)*yp(iy - 1)*zp(iz)
         if (ix > 0 .and. iz > 0) fxz = real(ix*iz, dp)*xp(ix - 1)*yp(iy)*zp(iz - 1)
         if (iy > 0 .and. iz > 0) fyz = real(iy*iz, dp)*xp(ix)*yp(iy - 1)*zp(iz - 1)

         hchi(1, m) = fxx*rad + 2.0_dp*fx*d(1)*drad + ang*drad + ang*d(1)*d(1)*ddrad
         hchi(2, m) = fxy*rad + drad*(fx*d(2) + fy*d(1)) + ang*d(1)*d(2)*ddrad
         hchi(3, m) = fxz*rad + drad*(fx*d(3) + fz*d(1)) + ang*d(1)*d(3)*ddrad
         hchi(4, m) = fyy*rad + 2.0_dp*fy*d(2)*drad + ang*drad + ang*d(2)*d(2)*ddrad
         hchi(5, m) = fyz*rad + drad*(fy*d(3) + fz*d(2)) + ang*d(2)*d(3)*ddrad
         hchi(6, m) = fzz*rad + 2.0_dp*fz*d(3)*drad + ang*drad + ang*d(3)*d(3)*ddrad
      end do
   end subroutine shell_collocate2

end module trc_collocation
