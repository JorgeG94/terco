!
! Radial parts of the effective-core-potential integrals.
!
! Ported from libfint's src/cint_ecp_rad.f90 (Apache-2.0, from PySCF's
! pyscf/lib/gto/nr_ecp.c), with two deliberate departures from a literal
! port, both load-bearing for running this as a `do concurrent` item
! (src/eri/trc_ecp.F90) rather than a per-call allocation:
!
!   * NO GENERAL CONTRACTIONS. `trc_basis_t` never holds one (every BSE
!     coefficient column becomes its own shell at basis-build time --
!     trc_basis_json.F90's own header), so the nc/ic loops PySCF's C and
!     libfint's port carry are dead weight here and are simply not written.
!     `ecp_radial_type2_facs_rad` below accumulates straight into its output
!     instead of through an intermediate per-primitive buffer the way
!     libfint's (which DOES support nc>1) has to.
!
!   * `ecp_radial_rad_part` ALWAYS evaluates the WHOLE finest grid (ECP_NRS
!     points, `trc_ecp_quadrature`'s own ceiling), not a per-level strided
!     subset with an early truncation scan. It is a deterministic,
!     per-point sum of a handful of Gaussians -- no part of it is adaptive
!     -- so evaluating it once per (ECP atom, channel) and letting every
!     shell pair's own adaptive level loop (in trc_ecp.F90) slice a strided
!     view of the SAME cached array is exactly libfint's result, just not
!     recomputed per shell pair. The early-truncation scan libfint's C
!     keeps for speed (skip multiplying by Bessel values once U(r) is
!     already <1e-50) is dropped: omitting it changes nothing past the
!     double-precision tolerance this module is validated to (the dropped
!     tail is already <1e-50 relative), and it would not apply across this
!     one-grid-fits-every-pair cache the way it did per shell pair anyway.
!
! The adaptive LEVEL loop itself -- the actual accuracy control, genuinely
! shell-pair dependent -- is not here: it lives in the item routines
! (trc_ecp.F90) that call `ecp_radial_type1_rad_part`/`ecp_radial_type2_facs_rad` once per
! level, exactly as PySCF's C does.
!
module trc_ecp_radial
   use pic_types, only: dp
   use trc_ecp_quadrature, only: ECP_NRS, ECP_EXPCUTOFF, ECP_CUTOFF, ECP_LMAX, ecp_quadrature_sph_ine
   use trc_ecp_angular, only: ECP_ORB_LMAX, ecp_cache_3dfac, CART_POW_Y, CART_POW_Z
   implicit none
   private

   public :: ecp_radial_rad_part, ecp_radial_type2_facs_rad, ecp_radial_type1_static_facs, ecp_radial_type1_rad_part

contains

   pure integer function cart_len(l) result(n)
      !$acc routine seq
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)/2
   end function cart_len

   !> U(r) on the WHOLE finest grid, summed over one group's rows
   !! (row0:row1), each row applying its OWN r-power to its OWN primitives
   !! -- see trc_ecp_types.F90 for why a channel can be more than one row.
   pure subroutine ecp_radial_rad_part(ur, rs, row0, row1, nrow, maxnp, ecp_rpow, ecp_np, ecp_e, ecp_c)
      !$acc routine seq
      integer, intent(in) :: row0, row1, nrow, maxnp
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      integer, intent(in) :: ecp_rpow(nrow), ecp_np(nrow)
      real(dp), intent(in) :: ecp_e(maxnp, nrow), ecp_c(maxnp, nrow)
      real(dp), intent(out) :: ur(0:ECP_NRS - 1)

      integer :: i, irow, kp, rpow, n
      real(dp) :: r2, s, rr

      do i = 0, ECP_NRS - 1
         r2 = rs(i)*rs(i)
         ur(i) = 0.0_dp
         do irow = row0, row1
            s = 0.0_dp
            do kp = 1, ecp_np(irow)
               s = s + ecp_c(kp, irow)*exp(-ecp_e(kp, irow)*r2)
            end do
            rpow = ecp_rpow(irow)
            rr = 1.0_dp
            do n = 1, rpow
               rr = rr*rs(i)
            end do
            ur(i) = ur(i) + s*rr
         end do
      end do
   end subroutine ecp_radial_rad_part

   !> Radial factors for the type-2 (projected) term, one orbital shell's
   !! primitives against one ECP channel's own angular momentum `lc`:
   !! `facs(i*lilc1+j) `, i = 0..nrs-1 (the level's own grid slice), j the
   !! lambda ladder 0..li+lc. Fixed-size for (ECP_ORB_LMAX, ECP_LMAX); the
   !! caller reads only the first `nrs*(li+lc+1)` entries.
   pure subroutine ecp_radial_type2_facs_rad(facs, li, np, e, c, lc, rca, rs, start, inc, nrs)
      !$acc routine seq
      integer, intent(in) :: li, np, lc, start, inc, nrs
      real(dp), intent(in) :: rca
      real(dp), intent(in) :: e(np), c(np)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      real(dp), intent(out) :: facs(0:ECP_NRS*(ECP_ORB_LMAX + ECP_LMAX + 1) - 1)

      integer :: lilc1, ip, i, j, flat_index
      real(dp) :: ka, ar2, t1, r2, rr
      real(dp) :: bess(0:ECP_ORB_LMAX + ECP_LMAX)

      lilc1 = li + lc + 1
      facs(0:nrs*lilc1 - 1) = 0.0_dp

      do i = 0, nrs - 1
         rr = rs(start + i*inc)
         r2 = (rr - rca)*(rr - rca)
         do ip = 1, np
            ar2 = e(ip)*r2
            if (ar2 > ECP_EXPCUTOFF + 6.0_dp) cycle
            t1 = exp(-ar2)
            ka = 2.0_dp*e(ip)*rca
            call ecp_quadrature_sph_ine(bess, li + lc, ka*rr)
            do j = 0, li + lc
               flat_index = i*lilc1 + j
               facs(flat_index) = facs(flat_index) + c(ip)*bess(j)*t1
            end do
         end do
      end do
   end subroutine ecp_radial_type2_facs_rad

   !> Centre-shift factors for the type-1 term, one block per Cartesian
   !! component of a shell of angular momentum `li` (<= ECP_ORB_LMAX).
   pure subroutine ecp_radial_type1_static_facs(facs, li, ri)
      !$acc routine seq
      integer, intent(in) :: li
      real(dp), intent(in) :: ri(3)
      real(dp), intent(out) :: facs(0:cart_len(ECP_ORB_LMAX)*(ECP_ORB_LMAX + 1)**3 - 1)

      real(dp) :: fac3d(0:3*(ECP_ORB_LMAX + 1)*(ECP_ORB_LMAX + 1) - 1)
      integer :: d1, d2, d3, mi, i, j, k, px, py, pz, ox, oy, oz

      d1 = li + 1
      d2 = d1*d1
      d3 = d2*d1
      ox = 0
      oy = d1*d1
      oz = 2*d1*d1

      call ecp_cache_3dfac(fac3d, li, ri)

      facs(0:cart_len(li)*d3 - 1) = 0.0_dp
      do mi = 0, cart_len(li) - 1
         py = CART_POW_Y(mi)
         pz = CART_POW_Z(mi)
         px = li - py - pz
         do i = 0, px
            do j = 0, py
               do k = 0, pz
                  facs(mi*d3 + i*d2 + j*d1 + k) = &
                     fac3d(ox + px*d1 + i)*fac3d(oy + py*d1 + j)*fac3d(oz + pz*d1 + k)
               end do
            end do
         end do
      end do
   end subroutine ecp_radial_type1_static_facs

   !> Radial integrals for the type-1 (local) term, one primitive pair,
   !! accumulated INTO `rad_all` (the caller zeroes it first, then calls
   !! this once per primitive pair -- unchanged from libfint). `rad_all` is
   !! fixed at (2*ECP_ORB_LMAX+1)^2, matching trc_ecp_angular.F90's
   !! `ecp_type1_rad_ang`, which reads it next.
   !!
   !! The parity stride is not an optimisation: terms where (lab+i) is odd
   !! integrate to zero over the sphere, so stepping by two skips exactly
   !! what the angular factor would multiply by nothing.
   pure subroutine ecp_radial_type1_rad_part(rad_all, lmax, k, aij, ur, rs, start, inc, nrs)
      !$acc routine seq
      integer, intent(in) :: lmax, start, inc, nrs
      real(dp), intent(in) :: k, aij
      real(dp), intent(in) :: ur(0:ECP_NRS - 1), rs(0:ECP_NRS - 1)
      real(dp), intent(inout) :: rad_all(0:(2*ECP_ORB_LMAX + 1)*(2*ECP_ORB_LMAX + 1) - 1)

      integer :: lmax1, lab, i, n
      real(dp) :: kaij, fac, temporary, s, rn
      real(dp) :: rur(0:ECP_NRS - 1), bval(0:ECP_NRS*(2*ECP_ORB_LMAX + 1) - 1)
      real(dp) :: bess(0:2*ECP_ORB_LMAX)

      if (nrs == 0) return
      lmax1 = lmax + 1

      kaij = k/(2.0_dp*aij)
      fac = kaij*kaij*aij
      do n = 0, nrs - 1
         rn = rs(start + n*inc)
         temporary = rn - kaij
         temporary = fac - aij*temporary*temporary
         if (ur(start + n*inc) == 0.0_dp .or. temporary > ECP_CUTOFF &
             .or. temporary < -(ECP_EXPCUTOFF + 6.0_dp + 30.0_dp)) then
            rur(n) = 0.0_dp
            bval(n*lmax1:n*lmax1 + lmax) = 0.0_dp
         else
            rur(n) = ur(start + n*inc)*exp(temporary)
            call ecp_quadrature_sph_ine(bess, lmax, k*rn)
            bval(n*lmax1:n*lmax1 + lmax) = bess(0:lmax)
         end if
      end do

      do lab = 0, lmax
         if (lab > 0) then
            do n = 0, nrs - 1
               rur(n) = rur(n)*rs(start + n*inc)
            end do
         end if
         do i = mod(lab, 2), lmax, 2
            s = rad_all(lab*lmax1 + i)
            do n = 0, nrs - 1
               s = s + rur(n)*bval(n*lmax1 + i)
            end do
            rad_all(lab*lmax1 + i) = s
         end do
      end do
   end subroutine ecp_radial_type1_rad_part

end module trc_ecp_radial
