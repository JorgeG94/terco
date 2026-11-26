!
! The ECP block kernels whose per-thread scratch grows with the shells'
! angular momenta: K1b's item, the (shell pair, ECP atom) block, the gradient
! task item and the two ECP terms. INCLUDEd after `contains` by each module
! that holds a `do concurrent` over them (stdpar only inlines within a
! translation unit), with the sizes of trc_ecp_block_sizes.inc, which says
! why there are two instances. The including module uses what these call:
! ECP_NRS, CART_POW_Y/Z, ecp_type1_rad_ang, ecp_type2_facs_ang and the
! radial routines of trc_ecp_radial.
!
   pure integer function ncart_l(l) result(n)
      !$acc routine seq
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)/2
   end function ncart_l

   pure subroutine ecp_k1b_item(ish, g, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                natm, at_r, ngroup, grp_atom, grp_l, rs, &
                                fac1_slot, fac2_slot, ang2_slot)
      !!
      !! K1b's item: one (orbital shell, ECP channel group). Writes ONLY the
      !! cache this group's own type (local or semi-local) needs -- the other
      !! two output arrays are left untouched, which is fine: K2 never reads a
      !! local group's fac2/ang2 slot or a semi-local group's fac1 slot.
      !! The slots have the instance's bra sizes (a shell's own values fill a
      !! prefix of each).
      !!
      !$acc routine seq
      integer, intent(in) :: ish, g, nshell, maxnp, natm, ngroup
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: at_r(3, natm)
      integer, intent(in) :: grp_atom(ngroup), grp_l(ngroup)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      real(dp), intent(out) :: fac1_slot(0:FAC1SIZE_BRA - 1)
      real(dp), intent(out) :: fac2_slot(0:FAC2SIZE_BRA - 1)
      real(dp), intent(out) :: ang2_slot(0:ANG2SIZE_BRA - 1)

      integer :: li, np, lc
      real(dp) :: rca(3)
      real(dp) :: dca

      li = sh_l(ish)
      np = sh_np(ish)
      lc = grp_l(g)
      rca = at_r(:, grp_atom(g)) - sh_r(:, ish)

      if (lc == -1) then
         call ecp_radial_type1_static_facs(fac1_slot, li, rca)
      else
         dca = sqrt(sum(rca**2))
         call ecp_radial_type2_facs_rad(fac2_slot, li, np, sh_e(1:np, ish), sh_c(1:np, ish), &
                                        lc, dca, rs, 0, 1, ECP_NRS)
         call ecp_type2_facs_ang(ang2_slot, li, lc, rca)
      end if
   end subroutine ecp_k1b_item

   pure subroutine ecp_atom_block(gctr, conv_fail, ish, jsh, atm_id, g0, g1, ci, cj, stride, nbra, nket, &
                                  ngroup, grp_l, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                  natm, at_r, rs, ws, ur_cache, fac1_bra, fac2_bra, ang2_bra, &
                                  fac1_ket, fac2_ket, ang2_ket)
      !!
      !! < shell ish | U of one ECP atom | shell jsh >, lexicographic, into
      !! `gctr`: the sum over that atom's channel groups g0..g1. The two
      !! shells' K1b columns of group g are `ci` and `cj` plus (g - g0)*stride:
      !! stride = nshell in the energy's full (shell, group) caches, 1 in the
      !! gradient's compacted ones. `ci` indexes the bra caches (`*_bra`,
      !! the instance's bra column sizes), `cj` the ket ones: one and the same
      !! arrays where the two sizes agree.
      !!
      !$acc routine seq
      integer, intent(in) :: ish, jsh, atm_id, g0, g1, ci, cj, stride, nbra, nket
      integer, intent(in) :: ngroup
      integer, intent(in) :: grp_l(ngroup)
      integer, intent(in) :: nshell, maxnp
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      integer, intent(in) :: natm
      real(dp), intent(in) :: at_r(3, natm)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), intent(in) :: ur_cache(0:ECP_NRS - 1, ngroup)
      real(dp), intent(in) :: fac1_bra(0:FAC1SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: fac2_bra(0:FAC2SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: ang2_bra(0:ANG2SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: fac1_ket(0:FAC1SIZE_KET - 1, nket)
      real(dp), intent(in) :: fac2_ket(0:FAC2SIZE_KET - 1, nket)
      real(dp), intent(in) :: ang2_ket(0:ANG2SIZE_KET - 1, nket)
      real(dp), intent(out) :: gctr(NFBRA, NFKET)
      integer, intent(out) :: conv_fail

      real(dp) :: rc(3), ri(3), rj(3)
      integer :: li, lj, npi, npj, g, ki, kj, local_fail

      li = sh_l(ish)
      lj = sh_l(jsh)
      npi = sh_np(ish)
      npj = sh_np(jsh)
      ri = sh_r(:, ish)
      rj = sh_r(:, jsh)
      rc = at_r(:, atm_id)
      gctr = 0.0_dp
      conv_fail = 0

      do g = g0, g1
         ki = ci + (g - g0)*stride
         kj = cj + (g - g0)*stride
         if (grp_l(g) == -1) then
            call ecp_type1_item(gctr, li, lj, npi, npj, &
                                sh_e(1:npi, ish), sh_c(1:npi, ish), sh_e(1:npj, jsh), sh_c(1:npj, jsh), &
                                ri, rj, rc, rs, ws, ur_cache(:, g), &
                                fac1_bra(:, ki), fac1_ket(:, kj), local_fail)
         else
            call ecp_type2_item(gctr, li, lj, grp_l(g), &
                                rs, ws, ur_cache(:, g), &
                                fac2_bra(:, ki), fac2_ket(:, kj), &
                                ang2_bra(:, ki), ang2_ket(:, kj), local_fail)
         end if
         conv_fail = max(conv_fail, local_fail)
      end do
   end subroutine ecp_atom_block

   pure subroutine ecp_grad_item(w, conv_fail, ish, jsh, nreal, atm_id, g0, g1, ci, cj, nbra, nket, &
                                 ngroup, grp_l, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                 natm, at_r, rs, ws, ur_cache, fac1_bra, fac2_bra, ang2_bra, &
                                 fac1_ket, fac2_ket, ang2_ket, nao, density, nfold, parent, coef, &
                                 nfket_ao, ket_ao)
      !!
      !! The gradient task kernel's item: derivative bra shell `ish`, real ket
      !! shell `jsh`, one ECP atom. Forms the block in its own scratch and
      !! contracts it with the density through the fold arrays (this file's
      !! header, THE GRADIENT) into its own three numbers. `nfold` and
      !! `nfket_ao` are the fold arrays' function extents, the same for both
      !! instances.
      !!
      !$acc routine seq
      integer, intent(in) :: ish, jsh, nreal, atm_id, g0, g1, ci, cj, nbra, nket
      integer, intent(in) :: nfold, nfket_ao
      integer, intent(in) :: ngroup
      integer, intent(in) :: grp_l(ngroup)
      integer, intent(in) :: nshell, maxnp, nao
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      integer, intent(in) :: natm
      real(dp), intent(in) :: at_r(3, natm)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), intent(in) :: ur_cache(0:ECP_NRS - 1, ngroup)
      real(dp), intent(in) :: fac1_bra(0:FAC1SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: fac2_bra(0:FAC2SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: ang2_bra(0:ANG2SIZE_BRA - 1, nbra)
      real(dp), intent(in) :: fac1_ket(0:FAC1SIZE_KET - 1, nket)
      real(dp), intent(in) :: fac2_ket(0:FAC2SIZE_KET - 1, nket)
      real(dp), intent(in) :: ang2_ket(0:ANG2SIZE_KET - 1, nket)
      real(dp), intent(in) :: density(nao, nao)
      integer, intent(in) :: parent(3, nfold, nshell - nreal)
      real(dp), intent(in) :: coef(3, nfold, nshell - nreal)
      integer, intent(in) :: ket_ao(nfket_ao, nreal)
      real(dp), intent(out) :: w(3)
      integer, intent(out) :: conv_fail

      real(dp) :: gctr(NFBRA, NFKET)
      real(dp) :: acc(3)
      integer :: f, mi, mj, c, nu

      call ecp_atom_block(gctr, conv_fail, ish, jsh, atm_id, g0, g1, ci, cj, 1, nbra, nket, &
                          ngroup, grp_l, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                          natm, at_r, rs, ws, ur_cache, fac1_bra, fac2_bra, ang2_bra, &
                          fac1_ket, fac2_ket, ang2_ket)
      f = ish - nreal
      acc = 0.0_dp
      do mj = 1, ncart_l(sh_l(jsh))
         nu = ket_ao(mj, jsh)
         do mi = 1, ncart_l(sh_l(ish))
            do c = 1, 3
               acc(c) = acc(c) + gctr(mi, mj)*coef(c, mi, f)*density(parent(c, mi, f), nu)
            end do
         end do
      end do
      ! Twice: the ket's derivative gives the same sum for a symmetric density.
      w = 2.0_dp*acc
   end subroutine ecp_grad_item

   pure logical function close_enough(x, y) result(c)
      !$acc routine seq
      real(dp), intent(in) :: x, y
      c = abs(x - y) <= 1.0e-12_dp*max(abs(x), abs(y))
   end function close_enough

   pure subroutine pick_converged_level(ld, d2, rad_levels, picked, failed)
      !! Picks the level a (primitive pair, type 1) or (lambda, type 2)
      !! block converged at, from its ECP_NLEVEL-column snapshot: the same
      !! "two consecutive close levels" test the original per-level loop
      !! ran while the block's own state persisted across levels, now run
      !! once over an already-materialised snapshot. `ld`, the
      !! dummy-bound explicit shape (like `sh_e(maxnp, nshell)` elsewhere
      !! in this file), is the CALLER's own fixed ceiling (D2MAX or
      !! D2MAX2), not the actual `d2`: passing the caller's array
      !! un-sliced, with its own real leading dimension, is what keeps
      !! this a reference rather than a copy-in/copy-out temporary sized
      !! by the runtime `d2` -- see `ecp_radial_type2_lab_levels`'s own
      !! header for the same trap this avoids. `picked` defaults to the
      !! last level and `failed` to true so a caller that never calls this
      !! for an exhausted block still gets the same "use the finest
      !! level, flag non-convergence" fallback the original fell through
      !! to.
      !$acc routine seq
      integer, intent(in) :: ld, d2
      real(dp), intent(in) :: rad_levels(0:ld - 1, 0:ECP_NLEVEL - 1)
      integer, intent(out) :: picked
      logical, intent(out) :: failed

      integer :: lvl, n, converged
      logical :: is_close

      converged = 0
      picked = ECP_NLEVEL - 1
      failed = .true.
      do lvl = 0, ECP_NLEVEL - 1
         is_close = .true.
         do n = 0, d2 - 1
            if (lvl == 0) then
               if (.not. close_enough(0.0_dp, rad_levels(n, lvl))) then
                  is_close = .false.
                  exit
               end if
            else
               if (.not. close_enough(rad_levels(n, lvl - 1), rad_levels(n, lvl))) then
                  is_close = .false.
                  exit
               end if
            end if
         end do
         if (is_close) then
            converged = converged + 1
            if (converged >= 2) then
               picked = lvl
               failed = .false.
               exit
            end if
         else
            converged = 0
         end if
      end do
   end subroutine pick_converged_level

   pure subroutine pick_converged_level_t2(ld, d2, rad_levels, picked, failed)
      !! Same as `pick_converged_level`, for a type-2 lambda block: a
      !! block whose every entry is already below ECP_NEGLIGIBLE is
      !! accepted at whatever level it first reaches that floor, without
      !! the relative test -- see ecp_type2_item's own header for why (a
      !! term already this small can fail a RELATIVE test indefinitely
      !! without ever being wrong to within this port's tolerance,
      !! `ecp_ints_ih_far`'s 25 bohr I-H).
      !$acc routine seq
      integer, intent(in) :: ld, d2
      real(dp), intent(in) :: rad_levels(0:ld - 1, 0:ECP_NLEVEL - 1)
      integer, intent(out) :: picked
      logical, intent(out) :: failed

      integer :: lvl, n, converged
      logical :: is_close

      converged = 0
      picked = ECP_NLEVEL - 1
      failed = .true.
      do lvl = 0, ECP_NLEVEL - 1
         if (maxval(abs(rad_levels(0:d2 - 1, lvl))) < ECP_NEGLIGIBLE) then
            is_close = .true.
         else
            is_close = .true.
            do n = 0, d2 - 1
               if (lvl == 0) then
                  if (.not. close_enough(0.0_dp, rad_levels(n, lvl))) then
                     is_close = .false.
                     exit
                  end if
               else
                  if (.not. close_enough(rad_levels(n, lvl - 1), rad_levels(n, lvl))) then
                     is_close = .false.
                     exit
                  end if
               end if
            end do
         end if
         if (is_close) then
            converged = converged + 1
            if (converged >= 2) then
               picked = lvl
               failed = .false.
               exit
            end if
         else
            converged = 0
         end if
      end do
   end subroutine pick_converged_level_t2

   pure subroutine ecp_type1_item(gctr, li, lj, npi, npj, ei, ci, ej, cj, ri, rj, rc, &
                                  rs, ws, ur, ifac, jfac, conv_fail)
      !!
      !! The type-1 (local) term for one (shell pair, ECP channel). Per-primitive-
      !! pair convergence, exactly as libfint's `ecp_type1_cart`: each (ip,jp)'s
      !! own radial integral is refined to its OWN converged level (picked by
      !! `pick_converged_level` from `ecp_radial_type1_rad_part_levels`'s whole
      !! snapshot, computed and consumed one pair at a time) and the angular
      !! contraction runs only once, after every pair has been picked (or the
      !! level ceiling is hit, which sets `conv_fail`). `ifac`/`jfac` are
      !! K1b's global cache, not rebuilt here.
      !!
      !$acc routine seq
      integer, intent(in) :: li, lj, npi, npj
      real(dp), intent(in) :: ei(npi), ci(npi), ej(npj), cj(npj)
      real(dp), intent(in) :: ri(3), rj(3), rc(3)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(in) :: ifac(0:FAC1SIZE_BRA - 1), jfac(0:FAC1SIZE_KET - 1)
      real(dp), intent(inout) :: gctr(NFBRA, NFKET)
      integer, intent(out) :: conv_fail

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: rca(3), rcb(3), rij(3)
      real(dp) :: cei(NPMAX), cej(NPMAX)
      real(dp) :: rad_levels(0:D2MAX - 1, 0:ECP_NLEVEL - 1)
      real(dp) :: rad_ang(0:D3MAX - 1), rad_ang_all(0:D3MAX - 1)
      integer :: lilj1, d1, d2, d3, di1, dj1
      integer :: ip, jp, mi, mj, n, i1, i2, i3, j1, j2, j3, ix, iy, iz, jx, jy, jz
      integer :: picked
      real(dp) :: r2ca, r2cb, fac
      logical :: failed

      lilj1 = li + lj + 1
      d1 = lilj1
      d2 = d1*d1
      d3 = d2*d1
      di1 = li + 1
      dj1 = lj + 1

      rca = rc - ri
      rcb = rc - rj
      r2ca = sum(rca**2)
      r2cb = sum(rcb**2)
      do ip = 1, npi
         cei(ip) = ci(ip)*exp(-ei(ip)*r2ca)*4.0_dp*PI
      end do
      do jp = 1, npj
         cej(jp) = cj(jp)*exp(-ej(jp)*r2cb)*4.0_dp*PI
      end do

      conv_fail = 0
      rad_ang_all = 0.0_dp
      do ip = 1, npi
         do jp = 1, npj
            ! A primitive pair whose centre-shift prefactor has already
            ! decayed past this floor contributes negligibly to gctr
            ! (fac*rad_ang below) however rad_levels' own magnitude looks:
            ! `ecp_radial_type1_rad_part_levels`'s exponent-shift (complete
            ! the square around kaij, rather than evaluate exp(-aij*rn^2)
            ! directly) trades a huge-but-finite intermediate for one that
            ! would otherwise underflow outright, and that intermediate's
            ! own growth tracks the SAME exponents this prefactor decays
            ! by, so it never overtakes a floor this far below any value
            ! that could matter (measured: a prefactor of 8.4e-251 next to
            ! a "val" of 8.6e+187, net contribution ~1e-63) -- so this pair
            ! never needs the quadrature run on it at all. A tight core
            ! primitive far from the ECP centre (ecp_ints_ih_far's 25 bohr
            ! I-H, ecp_ints_csh's Cs 1s next to the same-atom local channel)
            ! reaches this routinely.
            if (abs(cei(ip)*cej(jp)) < 1.0e-150_dp) cycle
            rij = ei(ip)*rca + ej(jp)*rcb
            call ecp_radial_type1_rad_part_levels(rad_levels, D2MAX, li + lj, &
                                      sqrt(sum(rij**2))*2.0_dp, ei(ip) + ej(jp), ur, rs, ws)
            call pick_converged_level(D2MAX, d2, rad_levels, picked, failed)
            if (failed) conv_fail = 1
            call ecp_type1_rad_ang(rad_ang, li + lj, rij, rad_levels(0:d2 - 1, picked))
            fac = cei(ip)*cej(jp)
            do n = 0, d3 - 1
               rad_ang_all(n) = rad_ang_all(n) + fac*rad_ang(n)
            end do
         end do
      end do

      do mi = 0, ncart_l(li) - 1
         iy = CART_POW_Y(mi)
         iz = CART_POW_Z(mi)
         ix = li - iy - iz
         do mj = 0, ncart_l(lj) - 1
            jy = CART_POW_Y(mj)
            jz = CART_POW_Z(mj)
            jx = lj - jy - jz
            do i1 = 0, ix
               do i2 = 0, iy
                  do i3 = 0, iz
                     do j1 = 0, jx
                        do j2 = 0, jy
                           do j3 = 0, jz
                              gctr(mi + 1, mj + 1) = gctr(mi + 1, mj + 1) &
                                 + ifac(mi*di1*di1*di1 + i1*di1*di1 + i2*di1 + i3) &
                                 *jfac(mj*dj1*dj1*dj1 + j1*dj1*dj1 + j2*dj1 + j3) &
                                 *rad_ang_all((i1 + j1)*d2 + (i2 + j2)*d1 + i3 + j3)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine ecp_type1_item

   pure subroutine ecp_type2_item(gctr, li, lj, lc, &
                                  rs, ws, ur, radi, radj, angi, angj, conv_fail)
      !!
      !! The type-2 (semi-local) term for one (shell pair, ECP channel of
      !! angular momentum lc). No nc/ic loop (trc_basis_t never holds a general
      !! contraction), so the convergence state is indexed by `lab` alone.
      !! `radi`/`radj`/`angi`/`angj` are K1b's global cache -- nothing here
      !! rebuilds a radial or angular factor, or holds one at ECP_NRS size:
      !! each `lab` (0..li+lj) is picked to its OWN converged level one at a
      !! time (`ecp_radial_type2_lab_levels` + `pick_converged_level_t2`), and
      !! only the picked, SLOT_SIZE-scale `rad_all` persists across lab values
      !! -- it has to, since the final contraction below indexes it by i+j,
      !! which repeats across different (i,j) pairs.
      !!
      !$acc routine seq
      integer, intent(in) :: li, lj, lc
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(in) :: radi(0:FAC2SIZE_BRA - 1), radj(0:FAC2SIZE_KET - 1)
      real(dp), intent(in) :: angi(0:ANG2SIZE_BRA - 1), angj(0:ANG2SIZE_KET - 1)
      real(dp), intent(inout) :: gctr(NFBRA, NFKET)
      integer, intent(out) :: conv_fail

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: common_fac
      real(dp) :: rad_levels(0:D2MAX2 - 1, 0:ECP_NLEVEL - 1)
      real(dp) :: rad_all(0:D3MAX2 - 1)
      real(dp) :: buffer(0:BUFMAX - 1)
      !! Written at n*ljlc1 + mj, n < ncart(li)*(2lc+1), mj < lj+lc+1: pyscf's
      !! nfi*(2*ECP_LMAX+1)*(lj+ECP_LMAX+1). Without the LC1KET factor a d-d
      !! pair against a d channel already overran it (#194 review).
      integer :: lilc1, ljlc1, d2, im, mq, dlc
      integer :: lab, i, j, n, mi, mj, picked
      real(dp) :: acc
      logical :: failed

      lilc1 = li + lc + 1
      ljlc1 = lj + lc + 1
      d2 = lilc1*ljlc1
      dlc = 2*lc + 1
      im = ncart_l(li)*dlc
      mq = dlc*ljlc1
      common_fac = 16.0_dp*PI*PI

      conv_fail = 0
      rad_all(0:(li + lj + 1)*d2 - 1) = 0.0_dp
      do lab = 0, li + lj
         call ecp_radial_type2_lab_levels(rad_levels, D2MAX2, lab, lilc1, ljlc1, ur, ws, rs, radi, radj)
         call pick_converged_level_t2(D2MAX2, d2, rad_levels, picked, failed)
         if (failed) conv_fail = 1
         do n = 0, d2 - 1
            rad_all(lab*d2 + n) = rad_levels(n, picked)
         end do
      end do

      do i = 0, li
         do j = 0, lj
            do n = 0, im - 1
               do mj = 0, ljlc1 - 1
                  acc = 0.0_dp
                  do mi = 0, lilc1 - 1
                     acc = acc + rad_all((i + j)*d2 + mi*ljlc1 + mj) &
                           *angi(i*ncart_l(li)*dlc*lilc1 + n*lilc1 + mi)
                  end do
                  buffer(n*ljlc1 + mj) = acc
               end do
            end do
            do mj = 0, ncart_l(lj) - 1
               do mi = 0, ncart_l(li) - 1
                  acc = 0.0_dp
                  do n = 0, mq - 1
                     acc = acc + buffer(mi*mq + n) &
                           *angj(j*ncart_l(lj)*dlc*ljlc1 + mj*mq + n)
                  end do
                  gctr(mi + 1, mj + 1) = gctr(mi + 1, mj + 1) + common_fac*acc
               end do
            end do
         end do
      end do
   end subroutine ecp_type2_item
