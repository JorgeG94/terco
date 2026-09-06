#!/usr/bin/env python3
"""
Emit one complete kernel per angular-momentum class: unrolled VRR, unrolled
HRR, six-atomic folded Fock, no `select case` anywhere.

WHY, AND WHY IT IS DIFFERENT FROM THE EARLIER ATTEMPT
-----------------------------------------------------
Specialising per class was tried before and lost, for two reasons that no
longer apply:

  * it fragmented the launches -- 81 small kernels on a work list that had not
    been binned, so nothing filled the device.  The bins now group work by
    class already, so a per-class launch is the natural shape rather than an
    imposed one.
  * it only unrolled loops, which lengthened live ranges and pushed registers
    from 146 to 241.

The argument now is different and specific: the single kernel currently holds
the unrolled code for EVERY class behind a `select case`, so ptxas must budget
registers for the worst class on every launch, and occupancy sits at ~11% with
255 registers.  Split per class and each kernel carries only its own
requirement -- an (ss|ss) launch should need a fraction of what (pp|pp) does.

Everything stays inline within its kernel.  The VRR measured 2.4x SLOWER behind
a module boundary and 1.58x faster inlined, so nothing here is allowed to
become a call.
"""

#: Radix for the dispatch key. FIXED, not LMAX+1.
#:
#: The key used to be ((la*(LMAX+1)+lb)*(LMAX+1)+lc)*(LMAX+1)+ld, which encodes
#: LMAX -- so a dispatcher generated at one LMAX and compiled at another still
#: BUILDS and dispatches the wrong kernel. A wrong number, not a build failure.
#: A fixed radix makes the key mean the same thing at every LMAX; 11 covers
#: l <= 10, far past anything that will ever be generated.
CLASS_RADIX = 11
#: Columns per primitive shell the kernels' coefficient scratch allows; must
#: match PS_NCOL_MAX in trc_bins.F90, which splits longer contractions.
PS_NCOL_MAX = 4

import argparse
import re
import importlib.util
import os
import sys


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = os.path.dirname(os.path.abspath(__file__))
gv = _load("gv", os.path.join(HERE, "gen_vrr.py"))
gh = _load("gh", os.path.join(HERE, "gen_hrr.py"))


PROLOGUE = """      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, npp, nao
      integer(kind=8), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      !! Coulomb and exchange scalings. Per call, so they fold into the six
      !! atomic updates and the digestion stays folded -- separating J and K
      !! into two matrices to scale them would give back FOCK6's 22%.
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3)
      real(dp), intent(in)    :: pp_rb(npp, 3), pp_c(npp)
      !! pp_c with the first column's coefficients folded in: what the
      !! scalar kernel multiplies, so a segmented basis pays nothing.
      real(dp), intent(in)    :: pp_cs(npp)
      !! primitive index of each pair's two primitives within their shells
      integer,  intent(in)    :: pp_ki(npp), pp_kj(npp)
      !! GENERAL CONTRACTION: nbas here counts PRIMITIVE shells; sh_l is
      !! theirs, ao_off is the first column's AO offset (the scalar kernel's
      !! only use of it; the blocked one reads the column table), and dsh
      !! is the screen folded to primitive shells. ps_coef holds each
      !! primitive shell's (np x ncol) coefficient matrix column-major from
      !! ps_coff(p)+1; column c of primitive shell p starts at AO
      !! col_ao(ps_soff(p)+c).
      integer,  intent(in)    :: ncoltot, ncoef
      integer,  intent(in)    :: ps_np(nbas), ps_ncol(nbas), ps_soff(nbas), ps_coff(nbas)
      integer,  intent(in)    :: col_ao(ncoltot)
      real(dp), intent(in)    :: ps_coef(ncoef)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao)"""



def emit_boys(m):
    """Straight-line Boys evaluation for a COMPILE-TIME order m.

    `boys_eval` is a cross-module `!$acc routine seq` taking an assumed-shape
    `f(0:)`, called from the innermost primitive loop.  Three costs, all
    avoidable here: NVHPC does not inline `routine seq` across a module
    boundary (the VRR measured 2.4x slower behind exactly that boundary), an
    assumed-shape dummy makes it walk a descriptor, and `m` being a dummy
    leaves the order loop and the Clenshaw loop with runtime bounds.

    Every generated kernel knows its own m, so none of that is necessary: the
    table lookups become literal offsets and the degree-5 Clenshaw recurrence
    unrolls into five FMAs per order.
    """
    o = []
    w = o.append
    w("               if (tval >= BOYS_TMAX) then")
    w("                  btt = sqrt(tval)")
    w("                  f(0) = 0.88622692545275801365_dp*erf(btt)/btt")
    w("                  bet = exp(-tval)")
    for k in range(1, m + 1):
        w(f"                  f({k}) = ({float(2*k - 1)}_dp*f({k - 1})"
          f" - bet)*(0.5_dp/tval)")
    w("               else")
    w("                  bi = int(tval*BOYS_DTINV)")
    w("                  if (bi >= BOYS_NGRID) bi = BOYS_NGRID - 1")
    w("                  bx = 2.0_dp*(tval - real(bi, dp)*BOYS_DT)*BOYS_DTINV"
      " - 1.0_dp")
    w("                  bx2 = 2.0_dp*bx")
    # SYMBOLIC, not a baked literal. This used to emit bi*54, which is
    # (BOYS_MMAX+1)*(BOYS_NCHEB+1) evaluated at MMAX=8 -- so regenerating the
    # Boys table at a different MMAX silently changed the table's stride while
    # every kernel kept indexing with the old one. It compiled, and it gave
    # wrong numbers.
    w("                  bbase = bi*(BOYS_MMAX + 1)*(BOYS_NCHEB + 1)")
    # The TOP order only is interpolated; the rest come by downward
    # recurrence, F_{k} = (2T F_{k+1} + exp(-T)) / (2k+1). One exp and one
    # reciprocal per primitive quartet against six scattered table loads per
    # order: Nsight Compute had the light classes saturated on L1 with the
    # Boys table about 40% of the traffic, and every thread lands on its own
    # interval so those loads never coalesce. Downward is also the stable
    # direction, and it is what gmshpc does.
    w(f"                  bj = bbase + {m}*(BOYS_NCHEB + 1)")
    # Clenshaw, degree 5: c_j sits at boys_table(bj + j + 1)
    w("                  b1 = boys_table(bj + BOYS_NCHEB + 1)")
    w("                  b2 = 0.0_dp")
    for jj in (5, 4, 3, 2):
        w(f"                  b0 = bx2*b1 - b2 + boys_table(bj + {jj});"
          f" b2 = b1; b1 = b0")
    w(f"                  f({m}) = bx*b1 - b2 + boys_table(bj + 1)")
    if m > 0:
        w("                  bet = exp(-tval)")
        w("                  btt = 2.0_dp*tval")
        for k in range(m - 1, -1, -1):
            w(f"                  f({k}) = (btt*f({k + 1}) + bet)*({1.0/(2*k + 1)!r}_dp)")
    w("               end if")
    return "\n".join(o)



def sieve_vrr(vrr_body, hrr_body, fin):
    """Drop VRR statements no HRR target can reach.

    THE SIEVE
    ---------
    The VRR is generated level by level with the classic triangular prune:
    at level m only degrees up to lt-m are built.  That prune is local -- it
    knows the recurrence but not what the recurrence is FOR.  The HRR then
    reads a strict subset of the final buffer: for (ab|cd) it needs
    [e0|f0] with l_e in [la, la+lb] and l_f in [lc, lc+ld], which for (pp|pp)
    is 81 of the 100 entries.  Everything computed only to feed one of the
    other 19 is dead work, and the generator has been emitting it.

    Combining the two DAGs and walking backwards from the HRR's actual reads
    finds it: mark the entries the HRR consumes as live, then sweep the VRR
    statements in reverse, keeping a statement only if its destination is
    live at that point and promoting its sources to live when it is.  This is
    Gill's sieve on the intermediates, applied across the VRR/HRR boundary
    rather than within the VRR alone.

    Measured before implementing: 27.0% of statements dead at (pp|pp), 25.7%
    at (pp|ps), 25.0% at (ps|ps).

    Purely structural -- the surviving statements are bit-for-bit the ones the
    old body already contained, in the same order, so results cannot change.
    """
    # split into statements, keeping Fortran continuations together
    raw, cur = [], []
    for line in vrr_body.split("\n"):
        cur.append(line)
        if not line.rstrip().endswith("&"):
            raw.append("\n".join(cur))
            cur = []
    if cur:
        raw.append("\n".join(cur))

    stmts = []
    for text in raw:
        flat = text.replace("&\n", " ")
        if "=" not in flat or flat.strip().startswith("!"):
            stmts.append((text, None, []))
            continue
        lhs, rhs = flat.split("=", 1)
        m = re.match(r"\s*v\((\d+),(\d+)\)", lhs)
        if not m:
            stmts.append((text, None, []))
            continue
        dst = (int(m.group(1)), int(m.group(2)))
        srcs = [(int(a), int(b))
                for a, b in re.findall(r"v\((\d+),(\d+)\)", rhs)]
        stmts.append((text, dst, srcs))

    live = {(int(i), fin) for i in set(re.findall(r"g\((\d+)\)", hrr_body))}
    keep = [False]*len(stmts)
    for k in range(len(stmts) - 1, -1, -1):
        text, dst, srcs = stmts[k]
        if dst is None:                    # comments ride along with the body
            continue
        if dst in live:
            keep[k] = True
            live.discard(dst)
            live.update(srcs)

    out = [t for k, (t, d, _) in enumerate(stmts)
           if d is None or keep[k]]
    return "\n".join(out)


def _emit_block(la, lb, lc, ld, cidx, vrr_body, hrr_body):
    lab, lcd = la + lb, lc + ld
    lt = lab + lcd
    nca, ncc = gv.ncum(lab), gv.ncum(lcd)
    na, nb = (la + 1)*(la + 2)//2, (lb + 1)*(lb + 2)//2
    nc, nd = (lc + 1)*(lc + 2)//2, (ld + 1)*(ld + 2)//2
    nv = nca*ncc
    tag = f"{la}{lb}{lc}{ld}"
    # Column combinations accumulated at once: the VRR block per
    # combination costs nv doubles of registers, so the chunk shrinks with
    # the class; a quartet with more combinations than this repeats its
    # primitive loops per chunk.
    # Every loop over the chunk has this CONSTANT bound, with zero weights
    # past the live combinations: a runtime bound put g into local memory
    # and cost 80x per primitive quartet. 16 doubles of block per chunk at
    # most on the small classes, one combination at a time on the large.
    # Column pairs accumulated per block on each side: NCAB x NCCD
    # combinations at once, constant bounds, so the block stays in
    # registers; a quartet with more column pairs repeats its primitive
    # loops per block. A single-column quartet takes the scalar path.
    # Sixteen doubles of accumulator at most: (ss|sp) at 4 x 4 x 4 went to
    # 255 registers with 624 bytes of spill, and ran slower than segmented.
    ncab = 4 if nv <= 1 else (2 if nv <= 4 else 1)
    # The block accumulation is written out with literal indices so the
    # accumulator is provably register-resident; a loop over (qab, qcd) is
    # left to the unroller's judgement, which it lost on this kernel.
    acc_zero = "".join(
        f"            do x = 1, {nv}\n               g(x, {a}, {c}) = 0.0_dp\n            end do\n"
        for c in range(1, ncab + 1) for a in range(1, ncab + 1))
    cab_fill = "".join(
        f"               cab({a}) = wta({a})*ps_coef(offa({a}) + ki)*ps_coef(offb({a}) + kj)\n"
        for a in range(1, ncab + 1))
    ccd_fill = "".join(
        f"                  ccd({c}) = wtc({c})*ps_coef(offc({c}) + kk)*ps_coef(offd({c}) + kl)\n"
        for c in range(1, ncab + 1))
    t_zero = "".join(
        f"               do x = 1, {nv}\n                  tc(x, {c}) = 0.0_dp\n               end do\n"
        for c in range(1, ncab + 1))
    t_add = "".join(
        f"                  do x = 1, {nv}\n                     tc(x, {c}) = tc(x, {c}) + ccd({c})*v(x, cur)\n                  end do\n"
        for c in range(1, ncab + 1))
    g_add = "".join(
        f"               do x = 1, {nv}\n                  g(x, {a}, {c}) = g(x, {a}, {c}) + cab({a})*tc(x, {c})\n               end do\n"
        for c in range(1, ncab + 1) for a in range(1, ncab + 1))
    acc_copy = "".join(
        f"            do x = 1, {nv}\n               gl(x, {a}, {c}) = g(x, {a}, {c})\n            end do\n"
        for c in range(1, ncab + 1) for a in range(1, ncab + 1))
    hrr_body = re.sub(r"\bg\(", "g1(", hrr_body)

    return f"""
   !> ({la}{lb}|{lc}{ld}) driver.  The `do concurrent` lives here and the
   !> workspaces live in the item routine below, so they are per THREAD.
   !> Declaring them alongside the loop makes them shared -- the compiler then
   !> emits `implicit copy(v, g, vbuf, f)` per launch and the threads race.
   subroutine pc{tag}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                      npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                      pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, &
                      ndens, dmat, jmat, rank, nranks)
{PROLOGUE}
      integer, intent(in) :: rank, nranks
      integer(kind=8) :: g0, g1, nr, i
      !
      ! LAUNCH GEOMETRY.  `do concurrent` gives nvfortran the whole say, and it
      ! picks 128 threads per block.  Handing the block size back to the
      ! compiler is the point rather than a compromise: the claim this library
      ! is making is that standard Fortran reaches competitive throughput, and
      ! a tuned launch geometry behind a directive would be quietly conceding
      ! it.  Per-class kernels are what recover the occupancy anyway -- ptxas
      ! budgets registers for one (la,lb,lc,ld) instead of for the worst class.
      !
      ! RANKS.  The work list is sorted, so rank r of n taking every n-th item
      ! from r is a static split that balances by construction, and every
      ! rank walks the same list -- the reduction of the result across ranks
      ! happens in the Fock driver, not here.  One rank is rank 0 of 1.
      g0 = sOff(lo) + 1 + rank
      g1 = sOff(hi + 1)
      nr = 0
      if (g1 >= g0) nr = (g1 - g0)/nranks + 1
      do concurrent(i=1:nr)
         call pci{tag}(g0 + (i - 1)*nranks, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                       npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                       pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, ndens, dmat, jmat)
      end do
   end subroutine pc{tag}

   pure subroutine pci{tag}(gt, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                      npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                      pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, &
                      ndens, dmat, jmat)
      !$acc routine seq
      integer(kind=8), intent(in) :: gt
{PROLOGUE}
      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl
      integer(kind=8) :: nsa, u, kx
      real(dp) :: qcut
      integer :: keyab, keycd, offab, offcd, nab, ncd
      integer :: kp, kq, d, x, cur, ia, ib, ic, id, idx, idens
      integer :: mu, nu, lam, sig, mui, nuj, lamk, sigl
      logical :: dij, dkl, dpq
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc
      real(dp) :: pqx, pqy, pqz, pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re, sc, vv
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      real(dp) :: f(0:BOYS_MMAX)
      integer  :: bi, bj, bbase
      real(dp) :: bx, bx2, b0, b1, b2, btt, bet
      real(dp) :: v({nv}, 0:1), g({nv}, {ncab}, {ncab}), gl({nv}, {ncab}, {ncab}), g1({nv}), vbuf({na*nb*nc*nd})
      real(dp) :: wq, w, wab
      integer  :: nca, ncb, nccl, ncdl, npi, npj, npk, npl, ncab, nccd, ab0, cd0, nab_c, ncd_c, qab, qcd
      integer  :: ki, kj, kk, kl, kpl, kql, ia2, ib2, ic2, id2, iabc, icdc
      logical  :: same_ab, same_cd, same_pair
      real(dp) :: cab({ncab}), ccd({ncab}), wta({ncab}), wtc({ncab}), tc({nv}, {ncab})
      integer  :: offa({ncab}), offb({ncab}), offc({ncab}), offd({ncab})
      real(dp) :: jab({na*nb}), jcd({nc*nd}), kac({na*nc})
      real(dp) :: kad({na*nd}), kbc({nb*nc}), kbd({nb*nd})
      real(dp) :: dab({na*nb}), dcd({nc*nd}), dac({na*nc})
      real(dp) :: dad({na*nd}), dbc({nb*nc}), dbd({nb*nd})

      ! locate the segment
         p = lo; q = hi
         do while (p < q)
            mid = (p + q + 1)/2
            if (sOff(mid) < gt) then
               p = mid
            else
               q = mid - 1
            end if
         end do
         seg = p
         t = int(gt - sOff(seg))

         ! KET-UNIFORM DECODE. The bra pair index runs fastest, so the 32
         ! threads of a warp hold 32 bra pairs against ONE ket pair, and the
         ! loads of the inner primitive loop -- the ket side -- are the same
         ! address across the warp: one L1 transaction, broadcast. With the
         ! ket index fastest, as this was, every thread pulled its own ket
         ! record per primitive quartet and Nsight Compute had (ss|ss) at 96%
         ! of L1 throughput with a 28% hit rate. The bra loads now diverge
         ! instead, once per bra primitive rather than once per quartet.
         ! On a symmetric segment the pairs are enumerated column by column,
         ! iab >= icd, with the same closed form inverted.
         if (sD(seg)) then
            nsa = int(sA(seg), 8)
            u = int(t - 1, 8)
            kx = int((real(2*nsa + 1, dp) - sqrt(real(2*nsa + 1, dp)**2 - 8.0_dp*real(u, dp)))/2.0_dp, 8)
            if (kx < 0) kx = 0
            do while (kx > 0)
               if ((kx*(2*nsa + 1) - kx*kx)/2 <= u) exit
               kx = kx - 1
            end do
            do while (((kx + 1)*(2*nsa + 1) - (kx + 1)*(kx + 1))/2 <= u)
               kx = kx + 1
            end do
            icd = int(kx) + 1
            iab = icd + int(u - (kx*(2*nsa + 1) - kx*kx)/2)
         else
            icd = (t - 1)/sA(seg) + 1
            iab = t - (icd - 1)*sA(seg)
         end if

         ! exact Schwarz, per quartet.  The host bin-pair test is only
         ! decade-granular (`int(-log10 Q)` clamped to 9) and therefore admits
         ! up to two decades too much; at RNA3 that was 46% of all quartets.
         qcut = sp_q(sOA(seg) + iab)*sp_q(sOB(seg) + icd)
         if (qcut <= thresh) return

         si = sp_i(sOA(seg) + iab); sj = sp_j(sOA(seg) + iab)
         sk = sp_i(sOB(seg) + icd); sl = sp_j(sOB(seg) + icd)

         ! Density-weighted screen.  What actually enters the Fock matrix is
         ! Q_ab Q_cd times a density block, so a quartet whose integrals are
         ! above threshold still contributes nothing if every density block it
         ! multiplies is negligible.  The J terms carry a factor 4 from the
         ! Coulomb degeneracy, the K terms 1.  For a converged density on an
         ! extended molecule this is where most of the screening comes from --
         ! with a flat model density it correctly fires on nothing.
         if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                      dsh(si, sk), dsh(si, sl), &
                      dsh(sj, sk), dsh(sj, sl)) <= thresh) return

         same_ab = (si == sj); same_cd = (sk == sl)
         same_pair = sD(seg) .and. (iab == icd)

         keyab = (si - 1)*nbas + sj
         keycd = (sk - 1)*nbas + sl
         offab = pp_off(keyab); nab = pp_n(keyab)
         offcd = pp_off(keycd); ncd = pp_n(keycd)
         nca = ps_ncol(si); ncb = ps_ncol(sj); nccl = ps_ncol(sk); ncdl = ps_ncol(sl)
         npi = ps_np(si); npj = ps_np(sj); npk = ps_np(sk); npl = ps_np(sl)

         abx = pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1)
         aby = pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2)
         abz = pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3)
         cdx = pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1)
         cdy = pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2)
         cdz = pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3)

         ! GENERAL CONTRACTION. The column pairs of (ab) and of (cd) are
         ! taken in blocks of NCAB x NCCD; per block the primitive loops run
         ! once and the VRR block is accumulated per combination, weighted by
         ! the two column-pair coefficient products, with constant loop
         ! bounds so it stays in registers. A single-column quartet -- every
         ! quartet of a segmented basis -- takes the scalar path, which is the
         ! kernel as it was with the coefficients read from the tables.
         ncab = nca*ncb; nccd = nccl*ncdl
         do ab0 = 1, ncab, {ncab}
         do cd0 = 1, nccd, {ncab}
         nab_c = min({ncab}, ncab - ab0 + 1)
         ncd_c = min({ncab}, nccd - cd0 + 1)
         if (ncab*nccd == 1) then
            do x = 1, {nv}
               g1(x) = 0.0_dp
            end do
            do kp = offab + 1, offab + nab
               zeta = pp_p(kp)
               kpl = kp - offab - 1
               ki = kpl/npj + 1
               kj = kpl - (ki - 1)*npj + 1
               wab = ps_coef(ps_coff(si) + ki)*ps_coef(ps_coff(sj) + kj)
               do kq = offcd + 1, offcd + ncd
                  eta = pp_p(kq)
                  kql = kq - offcd - 1
                  kk = kql/npl + 1
                  kl = kql - (kk - 1)*npl + 1
                  w = wab*ps_coef(ps_coff(sk) + kk)*ps_coef(ps_coff(sl) + kl)
                  zpe = zeta + eta
                  rho = zeta*eta/zpe
                  pqx = pp_r(kp, 1) - pp_r(kq, 1)
                  pqy = pp_r(kp, 2) - pp_r(kq, 2)
                  pqz = pp_r(kp, 3) - pp_r(kq, 3)
                  pax = pp_r(kp, 1) - pp_ra(kp, 1)
                  pay = pp_r(kp, 2) - pp_ra(kp, 2)
                  paz = pp_r(kp, 3) - pp_ra(kp, 3)
                  qcx = pp_r(kq, 1) - pp_ra(kq, 1)
                  qcy = pp_r(kq, 2) - pp_ra(kq, 2)
                  qcz = pp_r(kq, 3) - pp_ra(kq, 3)
                  wc = (zeta*pp_r(kp, 1) + eta*pp_r(kq, 1))/zpe
                  wpx = wc - pp_r(kp, 1); wqx = wc - pp_r(kq, 1)
                  wc = (zeta*pp_r(kp, 2) + eta*pp_r(kq, 2))/zpe
                  wpy = wc - pp_r(kp, 2); wqy = wc - pp_r(kq, 2)
                  wc = (zeta*pp_r(kp, 3) + eta*pp_r(kq, 3))/zpe
                  wpz = wc - pp_r(kp, 3); wqz = wc - pp_r(kq, 3)
                  tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)

   {emit_boys(lt)}

                  oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
                  rz = rho/zeta; re = rho/eta
                  pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

   {vrr_body}
                  do x = 1, {nv}
                     g1(x) = g1(x) + w*v(x, cur)
                  end do
               end do
            end do
            do x = 1, {nv}
               gl(x, 1, 1) = g1(x)
            end do
         else
{acc_zero}            ! The column-pair coefficient offsets of this block, decoded once
            ! per block and not per primitive: the integer divisions were
            ! most of the primitive loop. A slot past the live combinations
            ! reads column 1 with a zero weight, so there is no branch.
            do qab = 1, {ncab}
               iabc = min(ab0 + qab - 1, ncab)
               ia2 = mod(iabc - 1, nca) + 1; ib2 = (iabc - 1)/nca + 1
               offa(qab) = ps_coff(si) + (ia2 - 1)*npi
               offb(qab) = ps_coff(sj) + (ib2 - 1)*npj
               wta(qab) = merge(1.0_dp, 0.0_dp, ab0 + qab - 1 <= ncab)
               icdc = min(cd0 + qab - 1, nccd)
               ic2 = mod(icdc - 1, nccl) + 1; id2 = (icdc - 1)/nccl + 1
               offc(qab) = ps_coff(sk) + (ic2 - 1)*npk
               offd(qab) = ps_coff(sl) + (id2 - 1)*npl
               wtc(qab) = merge(1.0_dp, 0.0_dp, cd0 + qab - 1 <= nccd)
            end do
            do kp = offab + 1, offab + nab
               zeta = pp_p(kp)
               ki = pp_ki(kp); kj = pp_kj(kp)
{cab_fill}{t_zero}               do kq = offcd + 1, offcd + ncd
                  eta = pp_p(kq)
                  kk = pp_ki(kq); kl = pp_kj(kq)
{ccd_fill}                  zpe = zeta + eta
                  rho = zeta*eta/zpe
                  pqx = pp_r(kp, 1) - pp_r(kq, 1)
                  pqy = pp_r(kp, 2) - pp_r(kq, 2)
                  pqz = pp_r(kp, 3) - pp_r(kq, 3)
                  pax = pp_r(kp, 1) - pp_ra(kp, 1)
                  pay = pp_r(kp, 2) - pp_ra(kp, 2)
                  paz = pp_r(kp, 3) - pp_ra(kp, 3)
                  qcx = pp_r(kq, 1) - pp_ra(kq, 1)
                  qcy = pp_r(kq, 2) - pp_ra(kq, 2)
                  qcz = pp_r(kq, 3) - pp_ra(kq, 3)
                  wc = (zeta*pp_r(kp, 1) + eta*pp_r(kq, 1))/zpe
                  wpx = wc - pp_r(kp, 1); wqx = wc - pp_r(kq, 1)
                  wc = (zeta*pp_r(kp, 2) + eta*pp_r(kq, 2))/zpe
                  wpy = wc - pp_r(kp, 2); wqy = wc - pp_r(kq, 2)
                  wc = (zeta*pp_r(kp, 3) + eta*pp_r(kq, 3))/zpe
                  wpz = wc - pp_r(kp, 3); wqz = wc - pp_r(kq, 3)
                  tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)

   {emit_boys(lt)}

                  oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
                  rz = rho/zeta; re = rho/eta
                  pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

   {vrr_body}
                  ! cd side first: NCAB FMAs per primitive quartet, and the
                  ! NCAB x NCAB block only once per bra primitive.
{t_add}               end do
{g_add}            end do
            ! To local memory once: the per-combination read below has a
            ! runtime index and must not touch the register copy.
{acc_copy}         end if

         do qcd = 1, ncd_c
         do qab = 1, nab_c
         iabc = ab0 + qab - 1
         icdc = cd0 + qcd - 1
         ia2 = mod(iabc - 1, nca) + 1; ib2 = (iabc - 1)/nca + 1
         ic2 = mod(icdc - 1, nccl) + 1; id2 = (icdc - 1)/nccl + 1
         ! Only the canonical contracted quartets, as the contracted path
         ! enumerated them: the column pairs of one primitive shell in one
         ! order, and the two column pairs of one primitive-shell pair in one
         ! order.
         if (same_ab .and. ia2 < ib2) cycle
         if (same_cd .and. ic2 < id2) cycle
         if (same_pair .and. iabc < icdc) cycle
         dij = .not. (same_ab .and. ia2 == ib2)
         dkl = .not. (same_cd .and. ic2 == id2)
         dpq = .not. (same_pair .and. iabc == icdc)
         mui = col_ao(ps_soff(si) + ia2); nuj = col_ao(ps_soff(sj) + ib2)
         lamk = col_ao(ps_soff(sk) + ic2); sigl = col_ao(ps_soff(sl) + id2)
         do x = 1, {nv}
            g1(x) = gl(x, qab, qcd)
         end do

         ! --- HRR ---
{hrr_body}

#ifdef TRC_NO_DIGEST
         !
         ! Evaluation only, for like-for-like comparison against published
         ! numbers that were measured with digestion removed (the libERI paper
         ! edits QUICK's source to strip it, so its Tables 2 and 3 are ERI
         ! evaluation alone).  Every integral is still formed -- the whole VRR
         ! and HRR run and every component of vbuf is read, so nothing is
         ! dead-code eliminated -- but the density loads and the atomic
         ! scatters into the Fock matrix are gone.  The guard is opaque to the
         ! compiler and never fires, so jmat is untouched and the ANSWER IS
         ! DELIBERATELY WRONG.  Timing only.
         !
         sc = 0.0_dp
         do idx = 1, {na*nb*nc*nd}
            sc = sc + vbuf(idx)
         end do
         if (sc == huge(1.0_dp)) jmat(1, 1, 1) = sc
#else
         !
         ! BLOCK-ACCUMULATED DIGESTION.
         !
         ! The previous form did six atomic updates and six scattered `dmat`
         ! loads PER CARTESIAN COMPONENT.  For (pp|pp) that is 81 components x
         ! 6 = 486 of each, per quartet, even though only six small BLOCKS of
         ! jmat and six of dmat are ever touched.
         !
         ! gpu4pyscf's rys_contract_jk.cu does it the other way round: pull the
         ! density blocks into registers once (`load_dm`), contract every
         ! component against them there (`dot_dm`), and write each output block
         ! out once.  Same arithmetic, an order of magnitude less traffic --
         ! 54 loads and 54 atomics for (pp|pp) instead of 486.
         !
         ! The degeneracy factors are per-quartet, so they collapse into one
         ! scalar applied as the components are consumed.
         !
         wq = 1.0_dp
         if (.not. dij) wq = wq*0.5_dp
         if (.not. dkl) wq = wq*0.5_dp
         if (.not. dpq) wq = wq*0.5_dp

         !
         ! BATCHED OVER DENSITIES.
         !
         ! The integral is formed once, in `vbuf`, and contracted against every
         ! density in the batch. In the coupled-perturbed equations that is the
         ! difference between one integral pass and a hundred: the dynamic
         ! polarizabilities need nine perturbations times twelve imaginary
         ! frequencies, each a Fock build on a different response density.
         !
         ! The density loop is OUTSIDE the block accumulators, not inside, and
         ! that is the whole trick. The six blocks are zeroed, filled and
         ! written per density, so REGISTER PRESSURE DOES NOT GROW WITH THE
         ! BATCH -- holding N sets of blocks at once would have cost 54 more
         ! doubles per density at (pp|pp) and 216 at (dd|dd), on a kernel
         ! already spilling. Cost is `eval + ndens*digest`, and the ceiling is
         ! one over the digestion fraction.
         !
         do idens = 1, ndens

         do ib = 0, {nb - 1}
            do ia = 0, {na - 1}
               dab(1 + ia + {na}*ib) = dmat(idens, mui + ia, nuj + ib)
               jab(1 + ia + {na}*ib) = 0.0_dp
            end do
         end do
         do id = 0, {nd - 1}
            do ic = 0, {nc - 1}
               dcd(1 + ic + {nc}*id) = dmat(idens, lamk + ic, sigl + id)
               jcd(1 + ic + {nc}*id) = 0.0_dp
            end do
         end do
         do ic = 0, {nc - 1}
            do ia = 0, {na - 1}
               dac(1 + ia + {na}*ic) = dmat(idens, mui + ia, lamk + ic)
               kac(1 + ia + {na}*ic) = 0.0_dp
            end do
         end do
         do id = 0, {nd - 1}
            do ia = 0, {na - 1}
               dad(1 + ia + {na}*id) = dmat(idens, mui + ia, sigl + id)
               kad(1 + ia + {na}*id) = 0.0_dp
            end do
         end do
         do ic = 0, {nc - 1}
            do ib = 0, {nb - 1}
               dbc(1 + ib + {nb}*ic) = dmat(idens, nuj + ib, lamk + ic)
               kbc(1 + ib + {nb}*ic) = 0.0_dp
            end do
         end do
         do id = 0, {nd - 1}
            do ib = 0, {nb - 1}
               dbd(1 + ib + {nb}*id) = dmat(idens, nuj + ib, sigl + id)
               kbd(1 + ib + {nb}*id) = 0.0_dp
            end do
         end do

         idx = 0
         do id = 0, {nd - 1}
         do ic = 0, {nc - 1}
         do ib = 0, {nb - 1}
         do ia = 0, {na - 1}
            idx = idx + 1
            sc = wq*vbuf(idx)
            jab(1 + ia + {na}*ib) = jab(1 + ia + {na}*ib) &
                                    + 4.0_dp*jfac*sc*dcd(1 + ic + {nc}*id)
            jcd(1 + ic + {nc}*id) = jcd(1 + ic + {nc}*id) &
                                    + 4.0_dp*jfac*sc*dab(1 + ia + {na}*ib)
            kac(1 + ia + {na}*ic) = kac(1 + ia + {na}*ic) &
                                    - kfac*sc*dbd(1 + ib + {nb}*id)
            kad(1 + ia + {na}*id) = kad(1 + ia + {na}*id) &
                                    - kfac*sc*dbc(1 + ib + {nb}*ic)
            kbc(1 + ib + {nb}*ic) = kbc(1 + ib + {nb}*ic) &
                                    - kfac*sc*dad(1 + ia + {na}*id)
            kbd(1 + ib + {nb}*id) = kbd(1 + ib + {nb}*id) &
                                    - kfac*sc*dac(1 + ia + {na}*ic)
         end do
         end do
         end do
         end do

         do ib = 0, {nb - 1}
            do ia = 0, {na - 1}
               !$acc atomic update
               jmat(idens, mui + ia, nuj + ib) = jmat(idens, mui + ia, nuj + ib) &
                                          + jab(1 + ia + {na}*ib)
            end do
         end do
         do id = 0, {nd - 1}
            do ic = 0, {nc - 1}
               !$acc atomic update
               jmat(idens, lamk + ic, sigl + id) = jmat(idens, lamk + ic, sigl + id) &
                                            + jcd(1 + ic + {nc}*id)
            end do
         end do
         do ic = 0, {nc - 1}
            do ia = 0, {na - 1}
               !$acc atomic update
               jmat(idens, mui + ia, lamk + ic) = jmat(idens, mui + ia, lamk + ic) &
                                           + kac(1 + ia + {na}*ic)
            end do
         end do
         do id = 0, {nd - 1}
            do ia = 0, {na - 1}
               !$acc atomic update
               jmat(idens, mui + ia, sigl + id) = jmat(idens, mui + ia, sigl + id) &
                                           + kad(1 + ia + {na}*id)
            end do
         end do
         do ic = 0, {nc - 1}
            do ib = 0, {nb - 1}
               !$acc atomic update
               jmat(idens, nuj + ib, lamk + ic) = jmat(idens, nuj + ib, lamk + ic) &
                                           + kbc(1 + ib + {nb}*ic)
            end do
         end do
         do id = 0, {nd - 1}
            do ib = 0, {nb - 1}
               !$acc atomic update
               jmat(idens, nuj + ib, sigl + id) = jmat(idens, nuj + ib, sigl + id) &
                                           + kbd(1 + ib + {nb}*id)
            end do
         end do

         end do   ! idens
#endif
         end do   ! qab
         end do   ! qcd
         end do   ! cd0
         end do   ! ab0
   end subroutine pci{tag}
"""


def emit_kernel(la, lb, lc, ld, cidx, vrr_body, hrr_body):
    """Two kernels per class from one template.

    The BLOCK kernel takes the multi-column segments: its VRR block is
    accumulated per column combination.  The SCALAR kernel takes the
    single-column segments and is the kernel as it was before general
    contraction -- one accumulator, the pair coefficients folded into
    pp_cs, no column decode -- so a segmented basis, and every segmented
    part of a general one, runs exactly the code it ran before.  They are
    separate routines so ptxas budgets registers for each on its own and no
    warp carries both paths."""
    tag = f"{la}{lb}{lc}{ld}"
    txt = _emit_block(la, lb, lc, ld, cidx, vrr_body, hrr_body)

    m1 = "         if (ncab*nccd == 1) then\n"
    i1 = txt.index(m1)
    i2 = txt.index("\n         else\n            do x = 1, ", i1) + 1
    m3 = "         end if\n\n         do qcd = 1, ncd_c\n"
    i3 = txt.index(m3)
    scalar_prims = txt[i1 + len(m1):i2]
    block = txt[:i1] + txt[i2 + len("         else\n"):i3] + txt[i3 + len("         end if\n"):]

    # --- scalar item: head, one accumulator, folded coefficients, no decode
    ihead = txt.index("      ! locate the segment")
    igc = txt.index("         ! GENERAL CONTRACTION.")
    head = txt[ihead:igc]
    prims = scalar_prims
    for drop in ("               kpl = kp - offab - 1\n",
                 "               ki = kpl/npj + 1\n",
                 "               kj = kpl - (ki - 1)*npj + 1\n",
                 "               wab = ps_coef(ps_coff(si) + ki)*ps_coef(ps_coff(sj) + kj)\n",
                 "                  kql = kq - offcd - 1\n",
                 "                  kk = kql/npl + 1\n",
                 "                  kl = kql - (kk - 1)*npl + 1\n",
                 "                  w = wab*ps_coef(ps_coff(sk) + kk)*ps_coef(ps_coff(sl) + kl)\n"):
        assert prims.count(drop) == 1, drop
        prims = prims.replace(drop, "")
    prims = prims.replace("*pp_c(kp)*pp_c(kq)", "*pp_cs(kp)*pp_cs(kq)")
    prims = prims.replace("g1(x) = g1(x) + w*v(x, cur)", "g1(x) = g1(x) + v(x, cur)")
    igl = prims.rfind("            do x = 1, ")   # the gl copy, last loop
    assert igl > 0 and "gl(x, 1, 1) = g1(x)" in prims[igl:]
    prims = prims[:igl]
    ihrr = txt.index("         ! --- HRR ---\n")
    idig = txt.index("#ifdef TRC_NO_DIGEST")
    iend = txt.index("         end do   ! qab\n")
    hrr_and_digest = txt[ihrr:iend]
    decode = """         ! Single column on every side: the canonical enumeration and the
         ! degeneracy weights are the pair-level ones.
         dij = .not. same_ab
         dkl = .not. same_cd
         dpq = .not. same_pair
         mui = ao_off(si); nuj = ao_off(sj); lamk = ao_off(sk); sigl = ao_off(sl)

"""
    idecl0 = txt.index("      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl\n")
    decls = txt[idecl0:ihead]
    decls = re.sub(r"      real\(dp\) :: v\((\d+), 0:1\), g\(.*?\), gl\(.*?\), g1\((\d+)\), vbuf\((\d+)\)\n",
                   r"      real(dp) :: v(\1, 0:1), g1(\2), vbuf(\3)\n", decls)
    decls = re.sub(r"      real\(dp\) :: cab\(\d+\), ccd\(\d+\).*\n", "", decls)
    decls = re.sub(r"      integer  :: offa\(.*\n", "", decls)
    decls = decls.replace("      real(dp) :: wq, w, wab\n", "      real(dp) :: wq\n")
    decls = decls.replace("      integer  :: nca, ncb, nccl, ncdl, npi, npj, npk, npl, ncab, nccd, ab0, cd0, nab_c, ncd_c, qab, qcd\n", "")
    decls = decls.replace("      integer  :: ki, kj, kk, kl, kpl, kql, ia2, ib2, ic2, id2, iabc, icdc\n", "")
    head = head.replace("         nca = ps_ncol(si); ncb = ps_ncol(sj); nccl = ps_ncol(sk); ncdl = ps_ncol(sl)\n", "")
    head = head.replace("         npi = ps_np(si); npj = ps_np(sj); npk = ps_np(sk); npl = ps_np(sl)\n", "")

    isig = txt.index("   pure subroutine pci" + tag + "(")
    sig = txt[isig:idecl0].replace("pci" + tag, "pcsi" + tag)
    scalar_item = (sig + decls + head + prims + decode + hrr_and_digest
                   + "   end subroutine pcsi" + tag + "\n")

    # --- scalar driver: the block driver with the names swapped
    idrv = txt.index("   subroutine pc" + tag + "(")
    idrv_end = txt.index("   end subroutine pc" + tag + "\n") + len("   end subroutine pc" + tag + "\n")
    drv = txt[idrv:idrv_end].replace("pci" + tag, "pcsi" + tag).replace("pc" + tag, "pcs" + tag)
    drv = drv.replace("driver.  The `do concurrent`", "SCALAR driver.  The `do concurrent`")

    return block + "\n" + drv + "\n" + scalar_item


def _radix(txt):
    """Fill the radix placeholder. A placeholder rather than an f-string
    field because this block is emitted from two templates, one f-string
    and one plain, and a brace works in only one of them."""
    return txt.replace("__RADIX__", str(CLASS_RADIX))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=1)
    ap.add_argument("--split", metavar="DIR", default=None,
                    help="emit one module per class into DIR, plus a "
                         "dispatcher. At LMAX=2 the single-file form is "
                         "142k lines in one module and nvfortran's fort2 "
                         "grows past a gigabyte chewing on it; (22|22) alone "
                         "is 20k lines. The kernels are independent, and the "
                         "inlining rule only binds DEVICE code -- pci must "
                         "share a file with its pc driver, which per-class "
                         "files satisfy, while pc_dispatch is host code "
                         "launching kernels and may cross modules freely.")
    ap.add_argument("-o", "--output", default="src/trc_pc_kernels.F90")
    args = ap.parse_args()
    L = args.lmax

    lst, idx_v, info = gv.build_tables(2*L)
    idx_h = gh.cart_cum(4*L)

    pieces = []
    out = [f"""!
! One kernel per angular-momentum class: unrolled VRR, unrolled HRR, six-atomic
! folded Fock, and no `select case` in sight.
!
! GENERATED by scripts/gen_perclass.py --lmax {L} -- do not edit.
!
! The shared kernel holds the unrolled code for every class behind a select
! case, so ptxas budgets registers for the worst class on every launch --
! 255 registers and ~11% occupancy.  Each kernel here carries only its own
! requirement.
!
! Everything is inline within its kernel.  The VRR measured 2.4x slower behind
! a module boundary and 1.58x faster included, so nothing here may become a
! call.
!
module trc_pc_kernels
   use trc_boys, only: dp, boys_eval, BOYS_MMAX, boys_table, &
                        BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, &
                        BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   implicit none
   private
   public :: pc_dispatch, CLASS_RADIX

   !> Radix of the dispatch key, FIXED so it does not encode LMAX.
   !> The caller in trc_binkernel imports this rather than repeating the
   !> literal, because a key formed two ways that disagree dispatches the
   !> wrong kernel and still builds.
   integer, parameter :: CLASS_RADIX = __RADIX__

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains
"""]

    names = []
    for la in range(L + 1):
        for lb in range(L + 1):
            for lc in range(L + 1):
                for ld in range(L + 1):
                    # Gill's sieve needs to know what the HRR actually reads,
                    # so the HRR is generated first and its target set drives
                    # the VRR.
                    hrr_pre = gh.emit_class(la, lb, lc, ld, idx_h)
                    nca_ = gv.ncum(la + lb)
                    reads = {int(x) for x in re.findall(r"g\((\d+)\)", hrr_pre)}
                    targets = {((r - 1) % nca_ + 1, (r - 1)//nca_ + 1, 0)
                               for r in reads}
                    Pl = [None] + gv.cart_list(2*L)
                    Il = {q: i for i, q in enumerate(Pl) if q is not None}
                    ref, _f0, _n0 = gv.emit_class(la + lb, lc + ld, lst, idx_v, info)
                    vrr, _fin, nlive = gv.emit_class_sieved(
                        la + lb, lc + ld, targets, Pl, Il, nca_)
                    print(f"  ({la}{lb}|{lc}{ld}): VRR statements "
                          f"{len([x for x in ref.split(chr(10)) if '=' in x])}"
                          f" -> {nlive}")
                    vrr = "\n".join("      " + ln if ln.strip() else ln
                                    for ln in vrr.split("\n"))
                    vrr += f"\n               cur = {(la + lb + lc + ld) % 2}"
                    hrr = gh.emit_class(la, lb, lc, ld, idx_h)
                    key = ((la*CLASS_RADIX + lb)*CLASS_RADIX + lc)*CLASS_RADIX + ld
                    names.append((key, f"{la}{lb}{lc}{ld}"))
                    body = emit_kernel(la, lb, lc, ld, idx_h, vrr, hrr)
                    pieces.append((f"{la}{lb}{lc}{ld}", body))
                    out.append(body)

    out.append(f"""
   subroutine pc_dispatch(key, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                          pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, general, &
                          ndens, dmat, jmat, rank, nranks)
      integer, intent(in) :: key
{PROLOGUE}
      !! .true. sends the segments to the blocked kernel: every quartet in
      !! them has a multi-column side.  .false. is the scalar kernel.
      logical, intent(in) :: general
      integer, intent(in) :: rank, nranks
      if (general) then
      select case (key)""")
    for key, tag in names:
        out.append(f"""      case ({key}); call pc{tag}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                          pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, ndens, dmat, jmat, rank, nranks)""")
    out.append("      end select\n      else\n      select case (key)")
    for key, tag in names:
        out.append(f"""      case ({key}); call pcs{tag}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                          pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, pp_cs, pp_ki, pp_kj, ncoltot, ncoef, ps_np, ps_ncol, ps_soff, ps_coff, col_ao, ps_coef, ndens, dmat, jmat, rank, nranks)""")
    out.append("""      end select
      end if
   end subroutine pc_dispatch

end module trc_pc_kernels
""")

    txt = "\n".join(out)
    if args.split:
        import os
        os.makedirs(args.split, exist_ok=True)
        use_lines, names = [], []
        for tag, body in pieces:
            mod = f"trc_pc_k{tag}"
            names.append((tag, mod))
            fn = os.path.join(args.split, mod + ".F90")
            with open(fn, "w") as fh:
                fh.write(f"""!
! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) kernel.
!
! GENERATED by scripts/gen_perclass.py --lmax {L} --split -- do not edit.
!
! One module per class so no single compilation unit is enormous. The item
! routine stays in the same file as its driver, which is what the inlining
! rule requires; the dispatcher is host code and may cross modules.
!
module {mod}
   use trc_boys, only: dp, boys_eval, BOYS_MMAX, boys_table, &
                         BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, &
                         BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   implicit none
   private
   public :: pc{tag}, pcs{tag}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains
{body}
end module {mod}
""")
            use_lines.append(f"   use {mod}, only: pc{tag}, pcs{tag}")
        disp = [f"""!
! Kernel dispatcher for the per-class modules.
!
! GENERATED by scripts/gen_perclass.py --lmax {L} --split -- do not edit.
!
module trc_pc_kernels
   use trc_boys, only: dp
   use trc_tables, only: LMAX
""" + "\n".join(use_lines) + """
   implicit none
   private
   public :: pc_dispatch, CLASS_RADIX

   !> Radix of the dispatch key, FIXED so it does not encode LMAX.
   !> The caller in trc_binkernel imports this rather than repeating the
   !> literal, because a key formed two ways that disagree dispatches the
   !> wrong kernel and still builds.
   integer, parameter :: CLASS_RADIX = __RADIX__

contains
"""]
        disp.append(txt[txt.index("   subroutine pc_dispatch"):
                        txt.index("end module trc_pc_kernels")])
        disp.append("end module trc_pc_kernels\n")
        with open(os.path.join(args.split, "trc_pc_kernels.F90"), "w") as fh:
            fh.write(_radix("\n".join(disp)))
        print(f"wrote {len(names)} class modules + dispatcher into "
              f"{args.split}")
        return

    open(args.output, "w").write(_radix(txt))
    print(f"wrote {args.output}: {len(names)} kernels, {txt.count(chr(10))} lines")


if __name__ == "__main__":
    main()
