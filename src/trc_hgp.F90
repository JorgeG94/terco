!
! Head-Gordon-Pople: Obara-Saika vertical recurrences inside the primitive
! loop, horizontal transfer outside it.
!
! WHY, ALONGSIDE THE MMD PATH
! ---------------------------
! Asadchev & Valeev (arXiv:2307.03452) put matrix-form MD's arithmetic intensity
! at 4.2 FLOP/byte for (ff|ff) and 16.0 for (ii|ii), against a V100 machine
! balance of ~9, and say that at low angular momenta the scheme is memory-bound.
! We work at (dd|dd) and below, i.e. below their weakest case -- and every
! operation-reducing change to the MD kernel here moved nothing, which is what
! memory-bound looks like from the inside.
!
! HGP attacks that directly:
!
!   1. The horizontal transfer has no exponent dependence, so it runs ONCE on
!      contracted integrals instead of once per primitive combination.  For a
!      contracted basis that lifts most of the work out of the hot loop.
!   2. The hot-loop working set is [a0|c0], which for an s,p target is 10x10
!      doubles -- small enough to stay in registers, so intensity rises rather
!      than falls.  (At d it is 35x35 and does not, which is exactly the
!      crossover worth measuring.)
!
! THE HORIZONTAL TRANSFER, IN CLOSED FORM
! ---------------------------------------
! Rather than iterating (a,b+1_i|..) = (a+1_i,b|..) + AB_i (a,b|..) through
! intermediate b levels -- which needs a differently-shaped array at every
! level -- the recurrence is used in its unrolled form:
!
!   (a b|c d) = sum_{b' <= b} sum_{d' <= d} Cb(b') Cd(d') [a+b', 0 | c+d', 0]
!
!   Cb(b') = prod_i binom(b_i, b'_i) AB_i^(b_i - b'_i)
!
! One expression, no intermediate levels, obviously equivalent to the
! recurrence, and cheap because b and d are small: (b_x+1)(b_y+1)(b_z+1) terms,
! at most 4 on average through d.  It runs on contracted data, so the cost is
! per quartet rather than per primitive pair.
!
module trc_hgp
   use trc_boys, only: dp, BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   use trc_cart, only: LCMAX, NCUM, cidx, cnx, cny, cnz, cll, &
                         cdir, cdn1, cdn2, cf2, ncum_of, ncart_of
   implicit none
   private

   public :: hgp_batch, build_pairs_hgp, comp, binom, powi

   integer, parameter :: NPRIM_MAX = 1
   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains

#include "inc/trc_boys_eval.inc"

   !
   ! Primitive-pair data for HGP.  Differs from the MMD builder: no Hermite E
   ! coefficients, but the centres A and B are needed (for PA and AB) and the
   ! coefficient carries K_ab = exp(-mu |AB|^2), which in the MMD path was
   ! hiding inside E_000.
   !
   subroutine build_pairs_hgp(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, cfac, &
                              pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, npp)
      integer,  intent(in)  :: nbas
      integer,  intent(in)  :: sh_l(:), sh_np(:)
      real(dp), intent(in)  :: sh_e(:, :), sh_c(:, :), sh_r(:, :), cfac(:)
      integer,  allocatable, intent(out) :: pp_off(:), pp_n(:)
      real(dp), allocatable, intent(out) :: pp_p(:), pp_r(:, :), pp_ra(:, :), pp_rb(:, :), pp_c(:)
      integer,  intent(out) :: npp

      integer  :: i, j, ki, kj, key, k, d
      real(dp) :: a, b, p, mu, ab2

      allocate (pp_off(nbas*nbas), pp_n(nbas*nbas))
      npp = 0
      do i = 1, nbas
         do j = 1, nbas
            key = (i - 1)*nbas + j
            pp_off(key) = npp
            pp_n(key) = sh_np(i)*sh_np(j)
            npp = npp + pp_n(key)
         end do
      end do

      allocate (pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp))

      do i = 1, nbas
         do j = 1, nbas
            key = (i - 1)*nbas + j
            k = pp_off(key)
            ab2 = 0.0_dp
            do d = 1, 3
               ab2 = ab2 + (sh_r(d, i) - sh_r(d, j))**2
            end do
            do ki = 1, sh_np(i)
               a = sh_e(ki, i)
               do kj = 1, sh_np(j)
                  b = sh_e(kj, j)
                  k = k + 1
                  p = a + b
                  mu = a*b/p
                  pp_p(k) = p
                  do d = 1, 3
                     pp_r(k, d) = (a*sh_r(d, i) + b*sh_r(d, j))/p
                     pp_ra(k, d) = sh_r(d, i)
                     pp_rb(k, d) = sh_r(d, j)
                  end do
                  pp_c(k) = sh_c(ki, i)*cfac(i)*sh_c(kj, j)*cfac(j)*exp(-mu*ab2)
               end do
            end do
         end do
      end do
   end subroutine build_pairs_hgp

   subroutine hgp_batch(lo, hi, nq, nbas, npp, nout, &
                        q_i, q_j, q_k, q_l, q_off, sh_l, &
                        pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, out)
      integer,  intent(in)    :: lo, hi, nq, nbas, npp, nout
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      real(dp), intent(inout) :: out(nout)
      integer :: iq
      do concurrent(iq=lo:hi)
         call hgp_quartet(iq, nq, nbas, npp, nout, q_i, q_j, q_k, q_l, q_off, sh_l, &
                          pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, out)
      end do
   end subroutine hgp_batch

   pure subroutine hgp_quartet(iq, nq, nbas, npp, nout, &
                               q_i, q_j, q_k, q_l, q_off, sh_l, &
                               pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, out)
      !$acc routine seq
      integer,  intent(in)    :: iq, nq, nbas, npp, nout
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      real(dp), intent(inout) :: out(nout)

      integer  :: la, lb, lc, ld, lab, lcd, lt, nca, ncc
      integer  :: keyab, keycd, offab, offcd, nab, ncd
      integer  :: kp, kq, m, d, x, cur, nxt
      real(dp) :: zeta, eta, zpe, rho, tval, pref
      real(dp) :: pa(3), qc(3), wp(3), wq(3), pq(3), wc
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: f(0:BOYS_MMAX)
      ! [a0|c0]^m, two m levels ping-ponged, plus the contracted accumulator.
      ! At LMAX=1 (s,p -- what 6-31G is) NCUM=10, so this is 10x10x2 + 10x10
      ! doubles and stays small; at LMAX=2 it is 35x35 and does not.  That
      ! crossover is the thing worth measuring.
      real(dp) :: v(NCUM*NCUM, 0:1)
      real(dp) :: g(NCUM*NCUM)

      la = sh_l(q_i(iq)); lb = sh_l(q_j(iq))
      lc = sh_l(q_k(iq)); ld = sh_l(q_l(iq))
      lab = la + lb; lcd = lc + ld; lt = lab + lcd
      nca = ncum_of(lab); ncc = ncum_of(lcd)

      keyab = (q_i(iq) - 1)*nbas + q_j(iq)
      keycd = (q_k(iq) - 1)*nbas + q_l(iq)
      offab = pp_off(keyab); nab = pp_n(keyab)
      offcd = pp_off(keycd); ncd = pp_n(keycd)

      do x = 1, nca*ncc
         g(x) = 0.0_dp
      end do

      do kp = offab + 1, offab + nab
         zeta = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            zpe = zeta + eta
            rho = zeta*eta/zpe
            do d = 1, 3
               pq(d) = pp_r(kp, d) - pp_r(kq, d)
               wc = (zeta*pp_r(kp, d) + eta*pp_r(kq, d))/zpe
               pa(d) = pp_r(kp, d) - pp_ra(kp, d)
               qc(d) = pp_r(kq, d) - pp_ra(kq, d)
               wp(d) = wc - pp_r(kp, d)
               wq(d) = wc - pp_r(kq, d)
            end do
            tval = rho*(pq(1)*pq(1) + pq(2)*pq(2) + pq(3)*pq(3))

            call boys_eval(lt, tval, f)

            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
            rz = rho/zeta; re = rho/eta
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

            ! VRR, sweeping m downward so level m reads only level m+1.
            cur = 0
            v(1, cur) = pref*f(lt)
            do m = lt - 1, 0, -1
               nxt = 1 - cur
               call vrr_level(lt - m, nca, ncc, v, nxt, cur, pref*f(m), &
                              pa, qc, wp, wq, oo2z, oo2e, oo2ze, rz, re)
               cur = nxt
            end do

            do x = 1, nca*ncc
               g(x) = g(x) + v(x, cur)
            end do
         end do
      end do

      call hrr_store(iq, nq, nout, la, lb, lc, ld, nca, q_off, &
                     pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1), &
                     pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2), &
                     pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3), &
                     pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1), &
                     pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2), &
                     pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3), g, out)
   end subroutine hgp_quartet

   !
   ! Build every [a0|c0] of total degree <= lim at one m level.
   !
   ! `seed` is [00|00]^m, already carrying the prefactor.  Both loops rely on
   ! the angular-momentum-major ordering: a decrement always lands on a smaller
   ! index, so within a level everything read has already been written.
   !
   pure subroutine vrr_level(lim, nca, ncc, v, dst, src, seed, &
                             pa, qc, wp, wq, oo2z, oo2e, oo2ze, rz, re)
      !$acc routine seq
      integer,  intent(in)    :: lim, nca, ncc, dst, src
      real(dp), intent(inout) :: v(NCUM*NCUM, 0:1)
      real(dp), intent(in)    :: seed, pa(3), qc(3), wp(3), wq(3)
      real(dp), intent(in)    :: oo2z, oo2e, oo2ze, rz, re

      integer  :: ia, ic, d, x, xd1, xd2, xc1, xc2, xac, ai
      real(dp) :: acc

      v(1, dst) = seed

      ! a, at c = 0
      do ia = 2, nca
         if (cll(ia) > lim) exit
         d = cdir(ia); xd1 = cdn1(ia); xd2 = cdn2(ia)
         acc = pa(d)*v(xd1, dst) + wp(d)*v(xd1, src)
         if (xd2 > 0) acc = acc + cf2(ia)*oo2z*(v(xd2, dst) - rz*v(xd2, src))
         v(ia, dst) = acc
      end do

      ! c, for every a
      do ic = 2, ncc
         if (cll(ic) > lim) exit
         d = cdir(ic); xc1 = cdn1(ic); xc2 = cdn2(ic)
         do ia = 1, nca
            if (cll(ia) + cll(ic) > lim) exit
            x = ia + (ic - 1)*nca
            acc = qc(d)*v(ia + (xc1 - 1)*nca, dst) + wq(d)*v(ia + (xc1 - 1)*nca, src)
            if (xc2 > 0) acc = acc + cf2(ic)*oo2e* &
                               (v(ia + (xc2 - 1)*nca, dst) - re*v(ia + (xc2 - 1)*nca, src))
            ! cross term a_d/(2(zeta+eta)) [a-1_d,0|c-1_d,0]^(m+1)
            ai = 0; xac = 0
            if (d == 1) then
               ai = cnx(ia); if (ai > 0) xac = cidx(cnx(ia) - 1, cny(ia), cnz(ia))
            else if (d == 2) then
               ai = cny(ia); if (ai > 0) xac = cidx(cnx(ia), cny(ia) - 1, cnz(ia))
            else
               ai = cnz(ia); if (ai > 0) xac = cidx(cnx(ia), cny(ia), cnz(ia) - 1)
            end if
            if (xac > 0) acc = acc + real(ai, dp)*oo2ze*v(xac + (xc1 - 1)*nca, src)
            v(x, dst) = acc
         end do
      end do
   end subroutine vrr_level

   !
   ! Horizontal transfer in closed form, then store in libcint's buffer order
   ! (first index fastest).
   !
   pure subroutine hrr_store(iq, nq, nout, la, lb, lc, ld, nca, q_off, &
                             abx, aby, abz, cdx, cdy, cdz, g, out)
      !$acc routine seq
      integer,  intent(in)    :: iq, nq, nout, la, lb, lc, ld, nca
      integer,  intent(in)    :: q_off(nq)
      real(dp), intent(in)    :: abx, aby, abz, cdx, cdy, cdz
      real(dp), intent(in)    :: g(NCUM*NCUM)
      real(dp), intent(inout) :: out(nout)

      integer  :: na, nb, nc, nd, ia, ib, ic, id, idx, base
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: px, py, pz, qx, qy, qz, ja, jc
      real(dp) :: acc, cb, cdc

      na = ncart_of(la); nb = ncart_of(lb)
      nc = ncart_of(lc); nd = ncart_of(ld)
      base = q_off(iq)

      idx = 0
      do id = 0, nd - 1
         call comp(ld, id, dx, dy, dz)
      do ic = 0, nc - 1
         call comp(lc, ic, cx, cy, cz)
      do ib = 0, nb - 1
         call comp(lb, ib, bx, by, bz)
      do ia = 0, na - 1
         call comp(la, ia, ax, ay, az)
         idx = idx + 1
         acc = 0.0_dp
         do pz = 0, bz
         do py = 0, by
         do px = 0, bx
            cb = binom(bx, px)*binom(by, py)*binom(bz, pz) &
                 *powi(abx, bx - px)*powi(aby, by - py)*powi(abz, bz - pz)
            if (cb == 0.0_dp) cycle
            ja = cidx(ax + px, ay + py, az + pz)
            do qz = 0, dz
            do qy = 0, dy
            do qx = 0, dx
               cdc = binom(dx, qx)*binom(dy, qy)*binom(dz, qz) &
                     *powi(cdx, dx - qx)*powi(cdy, dy - qy)*powi(cdz, dz - qz)
               if (cdc == 0.0_dp) cycle
               jc = cidx(cx + qx, cy + qy, cz + qz)
               acc = acc + cb*cdc*g(ja + (jc - 1)*nca)
            end do
            end do
            end do
         end do
         end do
         end do
         out(base + idx) = out(base + idx) + acc
      end do
      end do
      end do
      end do
   end subroutine hrr_store

   pure subroutine comp(l, ic, px, py, pz)
      !$acc routine seq
      integer, intent(in)  :: l, ic
      integer, intent(out) :: px, py, pz
      integer :: n, lx, ly
      n = 0
      do lx = l, 0, -1
         do ly = l - lx, 0, -1
            if (n == ic) then
               px = lx; py = ly; pz = l - lx - ly
               return
            end if
            n = n + 1
         end do
      end do
      px = 0; py = 0; pz = 0
   end subroutine comp

   pure real(dp) function binom(n, k)
      !$acc routine seq
      integer, intent(in) :: n, k
      integer :: i
      real(dp) :: r
      r = 1.0_dp
      do i = 1, k
         r = r*real(n - k + i, dp)/real(i, dp)
      end do
      binom = r
   end function binom

   pure real(dp) function powi(x, n)
      !$acc routine seq
      real(dp), intent(in) :: x
      integer,  intent(in) :: n
      integer :: i
      powi = 1.0_dp
      do i = 1, n
         powi = powi*x
      end do
   end function powi

end module trc_hgp
