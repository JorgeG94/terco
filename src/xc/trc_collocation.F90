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
!   * the components of a shell come out in libcint's Cartesian order --
!     x descending outermost, y descending inside it -- which is the order
!     `trc_cart` uses and the order the integrals fill.
!
! The gradient is the same expression differentiated:
!
!     d/dx [ x^i y^j z^k R(r) ] = i x^(i-1) y^j z^k R  +  x^(i+1) y^j z^k R'
!
! with R' = sum_p (-2 a_p) c_p exp(-a_p r^2), so the radial part costs one
! extra multiply per primitive and nothing else.
!
module trc_collocation
   use trc_boys, only: dp
   implicit none
   private

   public :: shell_collocate, NCART_MAX, XC_LMAX

   ! The ceiling is generous on purpose: the loop is generic in l and a d
   ! shell costs no more than it would with a generated kernel, since the
   ! exponentials dominate.
   integer, parameter :: XC_LMAX = 4
   integer, parameter :: NCART_MAX = (XC_LMAX + 1)*(XC_LMAX + 2)/2

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
         ex = c(p)*exp(-e(p)*r2)
         rad = rad + ex
         drad = drad - 2.0_dp*e(p)*ex
      end do

      xp(0) = 1.0_dp; yp(0) = 1.0_dp; zp(0) = 1.0_dp
      do m = 1, l + 1
         xp(m) = xp(m - 1)*d(1)
         yp(m) = yp(m - 1)*d(2)
         zp(m) = zp(m - 1)*d(3)
      end do

      m = 0
      do ix = l, 0, -1
         do iy = l - ix, 0, -1
            iz = l - ix - iy
            m = m + 1
            ang = xp(ix)*yp(iy)*zp(iz)
            chi(m) = ang*rad
            gchi(1, m) = xp(ix + 1)*yp(iy)*zp(iz)*drad
            gchi(2, m) = xp(ix)*yp(iy + 1)*zp(iz)*drad
            gchi(3, m) = xp(ix)*yp(iy)*zp(iz + 1)*drad
            if (ix > 0) gchi(1, m) = gchi(1, m) + real(ix, dp)*xp(ix - 1)*yp(iy)*zp(iz)*rad
            if (iy > 0) gchi(2, m) = gchi(2, m) + real(iy, dp)*xp(ix)*yp(iy - 1)*zp(iz)*rad
            if (iz > 0) gchi(3, m) = gchi(3, m) + real(iz, dp)*xp(ix)*yp(iy)*zp(iz - 1)*rad
         end do
      end do
   end subroutine shell_collocate

end module trc_collocation
