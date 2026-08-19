!
! Batched 4-centre ERI over a work list of shell quartets.
!
! One source, two targets: gfortran runs the `do concurrent` serially, nvfortran
! with -stdpar=gpu offloads it.  That is deliberate -- a logic bug is found on
! the host in seconds, and what the GPU build then adds is only device concerns.
!
! THE SHAPE
! ---------
! Work item = one shell quartet.  Per-thread state is scalars plus the Boys
! vector; everything large lives in global memory:
!
!   E coefficients  precomputed ONCE per primitive pair (they are geometry, not
!                   kernel work) and read by every quartet that uses that pair.
!                   Recomputing them per quartet would be O(N^4) work on an
!                   O(N^2) quantity.
!   R tensor        a per-quartet slice of scratch, reused across the primitive
!                   loops of that quartet.
!   output          written directly; each quartet owns disjoint slots, so no
!                   atomics and a deterministic result.
!
! LAYOUT
! ------
! Every batched array puts the WORK-ITEM index first and contiguous, so
! consecutive threads touch consecutive addresses.  This is the opposite of the
! natural CPU layout and is the single highest-value decision here.
!
! KNOWN v1 COMPROMISE
! -------------------
! The R scratch is sized for the worst case (LTOT) for every quartet, whatever
! that quartet's actual total degree.  Correct but wasteful; the fix is to sort
! the work list by class and size per class.  Left until it is measured.
!
module trc_batch
   use trc_boys, only: dp, boys_eval, BOYS_MMAX
   use trc_tables, only: LMAX, LTOT, NHERM_MAX, NHERM_PAIR, nherm, tables_init, &
                           hidx, hdir, hm1, hm2, hcf, nherm_of, &
                           ht, hu, hv, hsgn, hshift
   implicit none
   private

   public :: e_off, build_pairs, eri_batch, fock_batch, fock_batch_sym, common_fac_sp, ncart
   public :: cart_pow, one_quartet

   ! E coefficients for one primitive pair: (i, j, t, direction),
   ! i,j <= LMAX and t <= 2*LMAX.
   integer, parameter, public :: NE_PAIR = 3*(LMAX + 1)*(LMAX + 1)*(2*LMAX + 1)

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp   !! 2*pi^(5/2)

contains

   pure integer function ncart(l)
      !$acc routine seq
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! libcint's per-shell angular constant.  See NOTES.md -- sqrt(1/4pi) for s,
   ! sqrt(3/4pi) for p, and 1 from d up, so its Cartesian d functions are
   ! deliberately not unit-normalised.
   !
   pure real(dp) function common_fac_sp(l)
      integer, intent(in) :: l
      select case (l)
      case (0); common_fac_sp = 0.282094791773878143_dp
      case (1); common_fac_sp = 0.488602511902919921_dp
      case default; common_fac_sp = 1.0_dp
      end select
   end function common_fac_sp

   !
   ! Offset of E(i,j,t,d) within one primitive pair's block.  1-based.
   !
   pure integer function e_off(i, j, t, d)
      !$acc routine seq
      integer, intent(in) :: i, j, t, d
      e_off = (d - 1)*(2*LMAX + 1)*(LMAX + 1)*(LMAX + 1) &
              + t*(LMAX + 1)*(LMAX + 1) + j*(LMAX + 1) + i + 1
   end function e_off

   !
   ! Host-side setup: expand every shell pair into its primitive pairs and
   ! store p, the product centre, the combined coefficient, and E.
   !
   ! pp_off/pp_n are indexed by the rectangular pair key (ish-1)*nbas + jsh,
   ! so the kernel finds a pair's primitives with two integer loads.
   !
   subroutine build_pairs(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, &
                          pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, npp)
      integer,  intent(in)  :: nbas
      integer,  intent(in)  :: sh_l(:), sh_np(:)
      real(dp), intent(in)  :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer,  allocatable, intent(out) :: pp_off(:), pp_n(:)
      real(dp), allocatable, intent(out) :: pp_p(:), pp_r(:, :), pp_c(:), pp_e(:, :)
      integer,  intent(out) :: npp

      integer  :: i, j, ki, kj, key, k, d, li, lj, t
      real(dp) :: a, b, p, mu, ci, cj
      real(dp) :: ra(3), rb(3), pc(3)
      real(dp) :: e1(0:LMAX, 0:LMAX, 0:2*LMAX)

      allocate (pp_off(nbas*nbas), pp_n(nbas*nbas))

      ! count first
      npp = 0
      do i = 1, nbas
         do j = 1, nbas
            key = (i - 1)*nbas + j
            pp_off(key) = npp
            pp_n(key) = sh_np(i)*sh_np(j)
            npp = npp + pp_n(key)
         end do
      end do

      allocate (pp_p(npp), pp_r(npp, 3), pp_c(npp), pp_e(npp, NE_PAIR))
      pp_e = 0.0_dp

      do i = 1, nbas
         li = sh_l(i); ra = sh_r(:, i)
         do j = 1, nbas
            lj = sh_l(j); rb = sh_r(:, j)
            key = (i - 1)*nbas + j
            k = pp_off(key)
            do ki = 1, sh_np(i)
               a = sh_e(ki, i); ci = sh_c(ki, i)*common_fac_sp(li)
               do kj = 1, sh_np(j)
                  b = sh_e(kj, j); cj = sh_c(kj, j)*common_fac_sp(lj)
                  k = k + 1
                  p = a + b
                  mu = a*b/p
                  pc = (a*ra + b*rb)/p
                  pp_p(k) = p
                  pp_r(k, :) = pc
                  pp_c(k) = ci*cj
                  do d = 1, 3
                     call hermite_e1d(li, lj, p, mu, ra(d) - rb(d), &
                                      pc(d) - ra(d), pc(d) - rb(d), e1)
                     ! The kernel indexes E by each COMPONENT's own powers,
                     ! not by the shell's l, so every (ii,jj) up to (li,lj) has
                     ! to be stored -- not only the top one.
                     block
                        integer :: ii, jj
                        do ii = 0, li
                           do jj = 0, lj
                              do t = 0, ii + jj
                                 pp_e(k, e_off(ii, jj, t, d)) = e1(ii, jj, t)
                              end do
                           end do
                        end do
                     end block
                  end do
               end do
            end do
         end do
      end do
   end subroutine build_pairs

   !
   ! One Cartesian direction's Hermite expansion coefficients E_t^{ij}.
   ! Host-only (setup), so an assumed-shape result is fine here.
   !
   pure subroutine hermite_e1d(li, lj, p, mu, xab, xpa, xpb, e)
      integer,  intent(in)  :: li, lj
      real(dp), intent(in)  :: p, mu, xab, xpa, xpb
      real(dp), intent(out) :: e(0:, 0:, 0:)
      integer  :: i, j, t
      real(dp) :: oo2p

      e = 0.0_dp
      oo2p = 0.5_dp/p
      e(0, 0, 0) = exp(-mu*xab*xab)

      do i = 0, li - 1
         do t = 0, i
            e(i + 1, 0, t) = xpa*e(i, 0, t)
            if (t > 0) e(i + 1, 0, t) = e(i + 1, 0, t) + oo2p*e(i, 0, t - 1)
            if (t < i) e(i + 1, 0, t) = e(i + 1, 0, t) + real(t + 1, dp)*e(i, 0, t + 1)
         end do
         e(i + 1, 0, i + 1) = oo2p*e(i, 0, i)
      end do

      do i = 0, li
         do j = 0, lj - 1
            do t = 0, i + j
               e(i, j + 1, t) = xpb*e(i, j, t)
               if (t > 0) e(i, j + 1, t) = e(i, j + 1, t) + oo2p*e(i, j, t - 1)
               if (t < i + j) e(i, j + 1, t) = e(i, j + 1, t) + real(t + 1, dp)*e(i, j, t + 1)
            end do
            e(i, j + 1, i + j + 1) = oo2p*e(i, j, i + j)
         end do
      end do
   end subroutine hermite_e1d

   !
   ! THE KERNEL.
   !
   ! Explicit-shape dummies with the extents declared first: an assumed-shape
   ! dummy makes nvfortran walk the descriptor with a per-launch memcpy.
   !
   subroutine eri_batch(lo, hi, nchunk, nq, nbas, npp, nout, &
                        q_i, q_j, q_k, q_l, q_off, &
                        sh_l, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, &
                        rscr, out)
      integer,  intent(in)    :: lo, hi, nchunk, nq, nbas, npp, nout
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_c(npp), pp_e(npp, NE_PAIR)
      real(dp), intent(inout) :: rscr(nchunk, 2*NHERM_MAX)
      real(dp), intent(inout) :: out(nout)

      integer :: iq

      do concurrent(iq=lo:hi)
         call one_quartet(iq, iq - lo + 1, nchunk, nq, nbas, npp, nout, &
                          q_i, q_j, q_k, q_l, q_off, sh_l, &
                          pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, rscr, out)
      end do
   end subroutine eri_batch

   !
   ! One work item.  Marked `acc routine seq`: stdpar only auto-inlines within a
   ! translation unit, and this is called from the parallel loop above.
   !
   pure subroutine one_quartet(iq, ir, nchunk, nq, nbas, npp, nout, &
                          q_i, q_j, q_k, q_l, q_off, sh_l, &
                          pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, rscr, out)
      !$acc routine seq
      integer,  intent(in)    :: iq, ir, nchunk, nq, nbas, npp, nout
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_c(npp), pp_e(npp, NE_PAIR)
      real(dp), intent(inout) :: rscr(nchunk, 2*NHERM_MAX)
      real(dp), intent(inout) :: out(nout)

      integer  :: la, lb, lc, ld, lab, lcd, lt
      integer  :: na, nb, nc, nd
      integer  :: keyab, keycd, offab, offcd, nab, ncd
      integer  :: kp, kq, n, h, deg, base, basen, base1, cur, nxt, rbase
      integer  :: ia, ib, ic, id, idx
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: h1, h2, nh_ab, nh_cd, iab, icd
      real(dp) :: p, q, alpha, tval, pref, acc, ecd, fac
      ! Small thread-locals, fixed extent (never sized by a dummy -- nvfortran
      ! crashes on that in a device routine).  Sized by NHERM_PAIR, the degree
      ! ONE PAIR reaches (35 at LMAX=2), not the quartet's NHERM_MAX (165):
      ! these are indexed by a pair's Hermite index and the wrong bound cost
      ! 2.6 kB of stack frame per thread.
      real(dp) :: tvec(NHERM_PAIR), ecdv(NHERM_PAIR)
      real(dp) :: pqx, pqy, pqz, pqd(3)
      real(dp) :: f(0:BOYS_MMAX)

      la = sh_l(q_i(iq)); lb = sh_l(q_j(iq))
      lc = sh_l(q_k(iq)); ld = sh_l(q_l(iq))
      lab = la + lb; lcd = lc + ld; lt = lab + lcd
      na = ncart(la); nb = ncart(lb); nc = ncart(lc); nd = ncart(ld)

      keyab = (q_i(iq) - 1)*nbas + q_j(iq)
      keycd = (q_k(iq) - 1)*nbas + q_l(iq)
      offab = pp_off(keyab); nab = pp_n(keyab)
      offcd = pp_off(keycd); ncd = pp_n(keycd)

      do kp = offab + 1, offab + nab
         p = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            q = pp_p(kq)
            alpha = p*q/(p + q)
            pqx = pp_r(kp, 1) - pp_r(kq, 1)
            pqy = pp_r(kp, 2) - pp_r(kq, 2)
            pqz = pp_r(kp, 3) - pp_r(kq, 3)
            pqd(1) = pqx; pqd(2) = pqy; pqd(3) = pqz
            tval = alpha*(pqx*pqx + pqy*pqy + pqz*pqz)

            call boys_eval(lt, tval, f)

            !
            ! R recurrence, PING-PONG over two buffers.
            !
            ! Level n is built only from level n+1, so the whole (LTOT+1)-level
            ! table was never needed -- two buffers are.  That cuts the scratch
            ! from NHERM_MAX*(LTOT+1) = 1485 doubles per work item to 2*165,
            ! i.e. 4.5x less footprint and, more to the point, 4.5x less DRAM
            ! traffic.  Measured to matter far more than the contraction did.
            !
            ! fac accumulates (-2 alpha)^n, which seeds R^n_000 at each level.
            !
            fac = 1.0_dp
            do n = 1, lt
               fac = fac*(-2.0_dp*alpha)
            end do
            cur = 0
            rscr(ir, cur*NHERM_MAX + 1) = fac*f(lt)

            do n = lt - 1, 0, -1
               nxt = 1 - cur
               fac = fac/(-2.0_dp*alpha)
               basen = nxt*NHERM_MAX
               base1 = cur*NHERM_MAX
               rscr(ir, basen + 1) = fac*f(n)
               deg = nherm_of(lt - n)
               do h = 2, deg
                  rscr(ir, basen + h) = pqd(hdir(h))*rscr(ir, base1 + hm1(h))
                  if (hm2(h) > 0) &
                     rscr(ir, basen + h) = rscr(ir, basen + h) &
                                           + hcf(h)*rscr(ir, base1 + hm2(h))
               end do
               cur = nxt
            end do
            rbase = cur*NHERM_MAX

            pref = TWO_PI_2_5/(p*q*sqrt(p + q))*pp_c(kp)*pp_c(kq)
            base = q_off(iq)
            nh_ab = nherm_of(lab)
            nh_cd = nherm_of(lcd)

            !
            ! TWO-PASS CONTRACTION.  The naive form evaluates
            !     sum_{tuv} sum_{TUV} E^ab E^cd R
            ! separately for every one of the ncart_ab * ncart_cd components,
            ! costing ncart_ab*ncart_cd*nherm_ab*nherm_cd.  Factoring it as
            !
            !     T(h_ab)      = sum_{h_cd} sgn E^cd(h_cd) R(h_ab + h_cd)
            !     (ab|cd)      = sum_{h_ab} E^ab(h_ab) T(h_ab)
            !
            ! costs nherm_ab*nherm_cd + ncart_ab*nherm_ab per ket component.
            ! For (dd|dd) that is 89.5k against 1.59M -- about 18x.
            !
            ! T is only nherm_ab long (<= 35), so it stays a small thread-local
            ! rather than needing scratch, which is the reason for doing one ket
            ! component at a time instead of forming the whole T matrix.
            !
            ! hshift() carries the index arithmetic, so the hot loop is a pure
            ! gather-multiply-add.
            !
            do icd = 0, nc*nd - 1
               id = icd/nc
               ic = icd - id*nc
               call cart_pow(lc, ic, cx, cy, cz)
               call cart_pow(ld, id, dx, dy, dz)

               ! E^cd over the flat Hermite index, sign folded in
               do h2 = 1, nh_cd
                  ecdv(h2) = 0.0_dp
               end do
               do h2 = 1, nh_cd
                  if (ht(h2) > cx + dx) cycle
                  if (hu(h2) > cy + dy) cycle
                  if (hv(h2) > cz + dz) cycle
                  ecdv(h2) = hsgn(h2) &
                             *pp_e(kq, e_off(cx, dx, ht(h2), 1)) &
                             *pp_e(kq, e_off(cy, dy, hu(h2), 2)) &
                             *pp_e(kq, e_off(cz, dz, hv(h2), 3))
               end do

               ! pass 1: T = R . E^cd
               do h1 = 1, nh_ab
                  tvec(h1) = 0.0_dp
               end do
               do h2 = 1, nh_cd
                  ecd = ecdv(h2)
                  if (ecd == 0.0_dp) cycle
                  do h1 = 1, nh_ab
                     tvec(h1) = tvec(h1) + ecd*rscr(ir, rbase + hshift(h1, h2))
                  end do
               end do

               ! pass 2: (ab|cd) = E^ab . T
               do iab = 0, na*nb - 1
                  ib = iab/na
                  ia = iab - ib*na
                  call cart_pow(la, ia, ax, ay, az)
                  call cart_pow(lb, ib, bx, by, bz)
                  acc = 0.0_dp
                  do h1 = 1, nh_ab
                     if (ht(h1) > ax + bx) cycle
                     if (hu(h1) > ay + by) cycle
                     if (hv(h1) > az + bz) cycle
                     acc = acc + pp_e(kp, e_off(ax, bx, ht(h1), 1)) &
                                 *pp_e(kp, e_off(ay, by, hu(h1), 2)) &
                                 *pp_e(kp, e_off(az, bz, hv(h1), 3))*tvec(h1)
                  end do
                  out(base + 1 + iab + na*nb*icd) = &
                     out(base + 1 + iab + na*nb*icd) + pref*acc
               end do
            end do
         end do
      end do
   end subroutine one_quartet


   !
   ! FUSED FOCK BUILD.  Integrals and digestion in one pass: each thread
   ! computes its quartet into its own slice of `out` and immediately contracts
   ! it, so the integrals never have to exist all at once.
   !
   !   J(mu,nu) += (mu nu|la si) D(la,si)
   !   K(mu,la) += (mu nu|la si) D(nu,si)
   !
   ! The work list here is the full rectangular quartet set, so every ordered
   ! (i,j,k,l) appears exactly once and no permutational prefactors are needed.
   ! That is the honest way to validate the contraction before symmetry is
   ! exploited; exploiting it is a later, separate change with its own check.
   !
   ! ATOMICS.  Unlike the integral pass, whose output slots are disjoint, the
   ! Fock elements are written by many quartets at once, so this needs
   ! `!$acc atomic update` -- verified working and exact inside a
   ! `do concurrent` by probe p01.  The consequence is that the result is NOT
   ! bit-reproducible run to run, because the FP addition order varies.  This is
   ! precisely the cost PLAN.md section 4.1 argues for designing away, and the
   ! structure that would avoid it for J is noted in NOTES.md.
   !
   subroutine fock_batch(lo, hi, nchunk, nq, nbas, npp, nout, nao, &
                         q_i, q_j, q_k, q_l, q_off, &
                         sh_l, ao_off, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, &
                         rscr, out, dmat, jmat, kmat)
      integer,  intent(in)    :: lo, hi, nchunk, nq, nbas, npp, nout, nao
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_c(npp), pp_e(npp, NE_PAIR)
      real(dp), intent(inout) :: rscr(nchunk, 2*NHERM_MAX)
      real(dp), intent(inout) :: out(nout)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: jmat(nao, nao), kmat(nao, nao)

      integer :: iq

      do concurrent(iq=lo:hi)
         call one_quartet(iq, iq - lo + 1, nchunk, nq, nbas, npp, nout, &
                          q_i, q_j, q_k, q_l, q_off, sh_l, &
                          pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, rscr, out)
         call digest_quartet(iq, nq, nbas, nout, nao, &
                             q_i, q_j, q_k, q_l, q_off, sh_l, ao_off, &
                             out, dmat, jmat, kmat)
      end do
   end subroutine fock_batch

   pure subroutine digest_quartet(iq, nq, nbas, nout, nao, &
                                  q_i, q_j, q_k, q_l, q_off, sh_l, ao_off, &
                                  out, dmat, jmat, kmat)
      !$acc routine seq
      integer,  intent(in)    :: iq, nq, nbas, nout, nao
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      real(dp), intent(in)    :: out(nout)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: jmat(nao, nao), kmat(nao, nao)

      integer  :: na, nb, nc, nd, ia, ib, ic, id, idx, base
      integer  :: mu, nu, lam, sig
      real(dp) :: v

      na = ncart(sh_l(q_i(iq))); nb = ncart(sh_l(q_j(iq)))
      nc = ncart(sh_l(q_k(iq))); nd = ncart(sh_l(q_l(iq)))
      base = q_off(iq)

      idx = 0
      do id = 0, nd - 1
         sig = ao_off(q_l(iq)) + id
      do ic = 0, nc - 1
         lam = ao_off(q_k(iq)) + ic
      do ib = 0, nb - 1
         nu = ao_off(q_j(iq)) + ib
      do ia = 0, na - 1
         mu = ao_off(q_i(iq)) + ia
         idx = idx + 1
         v = out(base + idx)
         !$acc atomic update
         jmat(mu, nu) = jmat(mu, nu) + v*dmat(lam, sig)
         !$acc atomic update
         kmat(mu, lam) = kmat(mu, lam) + v*dmat(nu, sig)
      end do
      end do
      end do
      end do
   end subroutine digest_quartet


   !
   ! SYMMETRY-REDUCED FUSED FOCK BUILD.
   !
   ! The work list holds only canonical quartets (i>=j, k>=l, pair_ij>=pair_kl),
   ! so each one stands for up to eight ordered tuples that share a value:
   !
   !   P1 (i,j,k,l)   P2 (j,i,k,l)   P3 (i,j,l,k)   P4 (j,i,l,k)
   !   P5 (k,l,i,j)   P6 (l,k,i,j)   P7 (k,l,j,i)   P8 (l,k,j,i)
   !
   ! Rather than derive a multiplicity prefactor -- which is where these
   ! routines usually go wrong, because the degenerate cases each want a
   ! different factor -- this enumerates the tuples and skips the ones that are
   ! not distinct.  Three booleans decide that: i/=j, k/=l, pair_ij/=pair_kl.
   ! Slightly more code, no factors to get wrong, and it reproduces the full
   ! rectangular result exactly, which is what the test asserts.
   !
   ! For each distinct tuple (a,b,c,d) the contributions are the same two the
   ! unsymmetrised version makes:
   !     J(a,b) += V * D(c,d)
   !     K(a,c) += V * D(b,d)
   !
   subroutine fock_batch_sym(lo, hi, nchunk, nq, nbas, npp, nout, nao, &
                             q_i, q_j, q_k, q_l, q_off, &
                             sh_l, ao_off, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, &
                             rscr, out, dmat, jmat, kmat)
      integer,  intent(in)    :: lo, hi, nchunk, nq, nbas, npp, nout, nao
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_c(npp), pp_e(npp, NE_PAIR)
      real(dp), intent(inout) :: rscr(nchunk, 2*NHERM_MAX)
      real(dp), intent(inout) :: out(nout)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: jmat(nao, nao), kmat(nao, nao)

      integer :: iq

      do concurrent(iq=lo:hi)
         call one_quartet(iq, iq - lo + 1, nchunk, nq, nbas, npp, nout, &
                          q_i, q_j, q_k, q_l, q_off, sh_l, &
                          pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, rscr, out)
         call digest_sym(iq, nq, nbas, nout, nao, &
                         q_i, q_j, q_k, q_l, q_off, sh_l, ao_off, &
                         out, dmat, jmat, kmat)
      end do
   end subroutine fock_batch_sym

   pure subroutine digest_sym(iq, nq, nbas, nout, nao, &
                              q_i, q_j, q_k, q_l, q_off, sh_l, ao_off, &
                              out, dmat, jmat, kmat)
      !$acc routine seq
      integer,  intent(in)    :: iq, nq, nbas, nout, nao
      integer,  intent(in)    :: q_i(nq), q_j(nq), q_k(nq), q_l(nq), q_off(nq)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      real(dp), intent(in)    :: out(nout)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: jmat(nao, nao), kmat(nao, nao)

      integer  :: na, nb, nc, nd, ia, ib, ic, id, idx, base
      integer  :: si, sj, sk, sl, mu, nu, lam, sig
      integer  :: pij, pkl
      logical  :: dij, dkl, dpq
      real(dp) :: v

      si = q_i(iq); sj = q_j(iq); sk = q_k(iq); sl = q_l(iq)
      na = ncart(sh_l(si)); nb = ncart(sh_l(sj))
      nc = ncart(sh_l(sk)); nd = ncart(sh_l(sl))
      base = q_off(iq)

      pij = si*(si - 1)/2 + sj
      pkl = sk*(sk - 1)/2 + sl
      dij = (si /= sj)
      dkl = (sk /= sl)
      dpq = (pij /= pkl)

      idx = 0
      do id = 0, nd - 1
         sig = ao_off(sl) + id
      do ic = 0, nc - 1
         lam = ao_off(sk) + ic
      do ib = 0, nb - 1
         nu = ao_off(sj) + ib
      do ia = 0, na - 1
         mu = ao_off(si) + ia
         idx = idx + 1
         v = out(base + idx)

         ! P1 (i,j,k,l)
         !$acc atomic update
         jmat(mu, nu) = jmat(mu, nu) + v*dmat(lam, sig)
         !$acc atomic update
         kmat(mu, lam) = kmat(mu, lam) + v*dmat(nu, sig)

         if (dij) then                          ! P2 (j,i,k,l)
            !$acc atomic update
            jmat(nu, mu) = jmat(nu, mu) + v*dmat(lam, sig)
            !$acc atomic update
            kmat(nu, lam) = kmat(nu, lam) + v*dmat(mu, sig)
         end if

         if (dkl) then                          ! P3 (i,j,l,k)
            !$acc atomic update
            jmat(mu, nu) = jmat(mu, nu) + v*dmat(sig, lam)
            !$acc atomic update
            kmat(mu, sig) = kmat(mu, sig) + v*dmat(nu, lam)
         end if

         if (dij .and. dkl) then                ! P4 (j,i,l,k)
            !$acc atomic update
            jmat(nu, mu) = jmat(nu, mu) + v*dmat(sig, lam)
            !$acc atomic update
            kmat(nu, sig) = kmat(nu, sig) + v*dmat(mu, lam)
         end if

         if (dpq) then
            ! P5 (k,l,i,j)
            !$acc atomic update
            jmat(lam, sig) = jmat(lam, sig) + v*dmat(mu, nu)
            !$acc atomic update
            kmat(lam, mu) = kmat(lam, mu) + v*dmat(sig, nu)

            if (dkl) then                       ! P6 (l,k,i,j)
               !$acc atomic update
               jmat(sig, lam) = jmat(sig, lam) + v*dmat(mu, nu)
               !$acc atomic update
               kmat(sig, mu) = kmat(sig, mu) + v*dmat(lam, nu)
            end if

            if (dij) then                       ! P7 (k,l,j,i)
               !$acc atomic update
               jmat(lam, sig) = jmat(lam, sig) + v*dmat(nu, mu)
               !$acc atomic update
               kmat(lam, nu) = kmat(lam, nu) + v*dmat(sig, mu)
            end if

            if (dij .and. dkl) then             ! P8 (l,k,j,i)
               !$acc atomic update
               jmat(sig, lam) = jmat(sig, lam) + v*dmat(nu, mu)
               !$acc atomic update
               kmat(sig, nu) = kmat(sig, nu) + v*dmat(lam, mu)
            end if
         end if
      end do
      end do
      end do
      end do
   end subroutine digest_sym

   !
   ! Cartesian powers of component `ic` of a shell of angular momentum l, in
   ! libcint's order (x descending outer, y descending inner).  Computed rather
   ! than looked up: it is a handful of integer ops against a table load, and
   ! it keeps one more array off the device.
   !
   pure subroutine cart_pow(l, ic, px, py, pz)
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
   end subroutine cart_pow

end module trc_batch
