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
! (trc_ecp.F90) that call `ecp_radial_type1_rad_part_levels`/
! `ecp_radial_type2_lab_levels`, each computing every level's cumulative
! value in one call (`trc_ecp.F90`'s own `pick_converged_level` then picks
! the level the original's per-level loop would have stopped at), exactly
! the SAME arithmetic PySCF's C runs once per level, just not laid out the
! same way in memory -- see this module's own header for why.
!
module trc_ecp_radial
   use pic_types, only: dp
   use trc_ecp_quadrature, only: ECP_NRS, ECP_EXPCUTOFF, ECP_CUTOFF, ECP_LMAX, &
                                 ECP_LEVEL0, ECP_LEVEL_MAX, ecp_quadrature_sph_ine
   use trc_ecp_angular, only: ECP_ORB_LMAX, ecp_cache_3dfac, CART_POW_Y, CART_POW_Z
   implicit none
   private

   public :: ecp_radial_rad_part, ecp_radial_type2_facs_rad, ecp_radial_type1_static_facs
   public :: ecp_radial_type1_rad_part_levels, ecp_radial_type2_lab_levels
   public :: ECP_NLEVEL

   integer, parameter :: ECP_NLEVEL = ECP_LEVEL_MAX - ECP_LEVEL0 + 1
   !! Highest level index + 1 (levels ECP_LEVEL0..ECP_LEVEL_MAX inclusive).

contains

   pure integer function cart_len(l) result(n)
      !$acc routine seq
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)/2
   end function cart_len

   pure subroutine ecp_radial_rad_part(ur, rs, row0, row1, nrow, maxnp, ecp_rpow, ecp_np, ecp_e, ecp_c)
      !! U(r) on the WHOLE finest grid, summed over one group's rows
      !! (row0:row1), each row applying its OWN r-power to its OWN
      !! primitives -- see trc_ecp_types.F90 for why a channel can be more
      !! than one row.
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

   pure subroutine ecp_radial_type2_facs_rad(facs, li, np, e, c, lc, rca, rs, start, inc, nrs)
      !! Radial factors for the type-2 (projected) term, one orbital
      !! shell's primitives against one ECP channel's own angular
      !! momentum `lc`: `facs(i*lilc1+j) `, i = 0..nrs-1 (the level's own
      !! grid slice), j the lambda ladder 0..li+lc. Fixed-size for
      !! (ECP_ORB_LMAX, ECP_LMAX); the caller reads only the first
      !! `nrs*(li+lc+1)` entries.
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

   pure subroutine ecp_radial_type1_static_facs(facs, li, ri)
      !! Centre-shift factors for the type-1 term, one block per Cartesian
      !! component of a shell of angular momentum `li` (<= ECP_ORB_LMAX).
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

   pure subroutine ecp_radial_type1_rad_part_levels(rad_levels, lmax, k, aij, ur, rs, ws)
      !! Radial integrals for the type-1 (local) term, one primitive pair,
      !! ALL levels at once: `rad_levels(:, lvl)` is the cumulative value
      !! through level `ECP_LEVEL0+lvl` (lvl=0..ECP_NLEVEL-1), exactly what
      !! the sequential halving loop this replaces would have held in
      !! `rad_all` after processing that level -- see trc_ecp.F90's own
      !! header for why a `do concurrent` item cannot hold the
      !! ECP_NRS-sized intermediates (`rur`/`bval`) the original per-level
      !! call built, and reads this instead of calling once per level: the
      !! grid loop here is FUSED with the angular-momentum ladder (no
      !! per-point array, `rn**lab` built incrementally per point, same
      !! trick ecp_radial_type2_facs_rad's caller already uses for its own
      !! r^lab ladder) so nothing here scales with ECP_NRS except the loop
      !! trip count.
      !!
      !! The parity stride is not an optimisation: terms where (lab+i) is
      !! odd integrate to zero over the sphere, so stepping by two skips
      !! exactly what the angular factor would multiply by nothing.
      !$acc routine seq
      integer, intent(in) :: lmax
      real(dp), intent(in) :: k, aij
      real(dp), intent(in) :: ur(0:ECP_NRS - 1), rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), intent(out) :: rad_levels(0:(2*ECP_ORB_LMAX + 1)*(2*ECP_ORB_LMAX + 1) - 1, 0:ECP_NLEVEL - 1)

      integer :: lmax1, lab, i, n, lvl, level, nrs0, start, step, nrs
      real(dp) :: kaij, fac, temporary, rn, wtscale, rur0, rpow_n
      real(dp) :: bess(0:2*ECP_ORB_LMAX)

      lmax1 = lmax + 1
      kaij = k/(2.0_dp*aij)
      fac = kaij*kaij*aij

      nrs0 = 2**ECP_LEVEL0 - 1
      step = 2**(ECP_LEVEL_MAX - ECP_LEVEL0)
      start = step - 1
      wtscale = real(step, dp)

      do lvl = 0, ECP_NLEVEL - 1
         nrs = nrs0
         if (lvl == 0) then
            do i = 0, lmax1*lmax1 - 1
               rad_levels(i, lvl) = 0.0_dp
            end do
         else
            do i = 0, lmax1*lmax1 - 1
               rad_levels(i, lvl) = rad_levels(i, lvl - 1)*0.5_dp
            end do
         end if

         do n = 0, nrs - 1
            rn = rs(start + n*step)
            temporary = rn - kaij
            temporary = fac - aij*temporary*temporary
            if (ur(start + n*step) == 0.0_dp .or. temporary > ECP_CUTOFF &
                .or. temporary < -(ECP_EXPCUTOFF + 6.0_dp + 30.0_dp)) cycle
            rur0 = ur(start + n*step)*ws(start + n*step)*wtscale*exp(temporary)
            call ecp_quadrature_sph_ine(bess, lmax, k*rn)
            rpow_n = 1.0_dp
            do lab = 0, lmax
               do i = mod(lab, 2), lmax, 2
                  rad_levels(lab*lmax1 + i, lvl) = rad_levels(lab*lmax1 + i, lvl) + rur0*rpow_n*bess(i)
               end do
               rpow_n = rpow_n*rn
            end do
         end do

         if (lvl < ECP_NLEVEL - 1) then
            level = ECP_LEVEL0 + lvl
            nrs0 = 2**level - 1
            step = 2**(ECP_LEVEL_MAX - level)
            start = (start - 1)/2
            wtscale = wtscale*0.5_dp
         end if
      end do
   end subroutine ecp_radial_type1_rad_part_levels

   pure subroutine ecp_radial_type2_lab_levels(rad_levels, ld, lab, lilc1, ljlc1, ur, ws, rs, radi, radj)
      !! Radial integrals for the type-2 (projected) term, one lambda
      !! value `lab` (0..li+lj), ALL levels at once -- the same
      !! cumulative-snapshot trick as `ecp_radial_type1_rad_part_levels`
      !! above, for the SAME reason: trc_ecp.F90's item cannot hold an
      !! ECP_NRS-sized `ur_level`/`rpow` per level (the original,
      !! per-level call built both). `radi`/`radj` are READ here, not
      !! built: they are `trc_ecp.F90`'s own global per-(shell, channel)
      !! cache, filled once per geometry by `ecp_radial_type2_facs_rad` on
      !! the whole grid, so nothing here scales with ECP_NRS except the
      !! loop trip count, and `lilc1`/`ljlc1` passed in MUST be the same
      !! values the cache was filled with (same shell, same channel) -- a
      !! flat 1-D array indexed by a value shared between writer and
      !! reader, not a declared leading dimension, so there is no
      !! maxnp-style stride trap there. `rad_levels` itself DOES have a
      !! declared leading dimension, `ld` -- the caller's own fixed
      !! ceiling, not lilc1*ljlc1 -- so a caller whose own local array is
      !! sized to its OWN module's larger constant can pass it un-sliced:
      !! slicing it to (0:lilc1*ljlc1-1, :) first would hand this routine a
      !! NON-contiguous section (every "column" strided by the caller's
      !! bigger ld, not by lilc1*ljlc1), which Fortran would then have to
      !! repack through a copy-in/copy-out temporary sized by a RUNTIME
      !! value -- exactly the automatic-array-in-`routine seq` trap this
      !! whole file exists to avoid.
      !$acc routine seq
      integer, intent(in) :: ld, lab, lilc1, ljlc1
      real(dp), intent(in) :: ur(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), rs(0:ECP_NRS - 1)
      real(dp), intent(in) :: radi(0:ECP_NRS*lilc1 - 1), radj(0:ECP_NRS*ljlc1 - 1)
      real(dp), intent(out) :: rad_levels(0:ld - 1, 0:ECP_NLEVEL - 1)

      integer :: d2, i, j, n, pt, lvl, level, nrs0, start, step, nrs
      real(dp) :: wtscale, pointval, rn

      d2 = lilc1*ljlc1
      nrs0 = 2**ECP_LEVEL0 - 1
      step = 2**(ECP_LEVEL_MAX - ECP_LEVEL0)
      start = step - 1
      wtscale = real(step, dp)

      do lvl = 0, ECP_NLEVEL - 1
         nrs = nrs0
         if (lvl == 0) then
            do i = 0, d2 - 1
               rad_levels(i, lvl) = 0.0_dp
            end do
         else
            do i = 0, d2 - 1
               rad_levels(i, lvl) = rad_levels(i, lvl - 1)*0.5_dp
            end do
         end if

         do n = 0, nrs - 1
            ! `radi`/`radj` are K1b's cache, filled on the NATURAL grid
            ! index (start=0, inc=1, nrs=ECP_NRS -- the whole grid, once),
            ! not on this level's densified one: `pt`, not `n`, is the
            ! point this level's n-th sample actually is, and a cached
            ! factor is read at the point it was built for, exactly like
            ! `ur`/`ws`/`rs` themselves. Reading it at `n` instead (this
            ! port's own first version of this routine did) silently reads
            ! a DIFFERENT point's factor at every level but the finest.
            pt = start + n*step
            rn = rs(pt)
            pointval = ur(pt)*ws(pt)*wtscale*rn**lab
            do i = 0, lilc1 - 1
               do j = 0, ljlc1 - 1
                  rad_levels(i*ljlc1 + j, lvl) = rad_levels(i*ljlc1 + j, lvl) &
                     + pointval*radi(pt*lilc1 + i)*radj(pt*ljlc1 + j)
               end do
            end do
         end do

         if (lvl < ECP_NLEVEL - 1) then
            level = ECP_LEVEL0 + lvl
            nrs0 = 2**level - 1
            step = 2**(ECP_LEVEL_MAX - level)
            start = (start - 1)/2
            wtscale = wtscale*0.5_dp
         end if
      end do
   end subroutine ecp_radial_type2_lab_levels

end module trc_ecp_radial
