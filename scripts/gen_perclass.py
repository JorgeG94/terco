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


PROLOGUE = """      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, nao
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      !! Primitive base of the bra bin / ket bin, into qp_* below -- the
      !! bin-interleaved second copy of pp_p/pp_r/pp_ra/pp_c, laid out so
      !! that lane iab's primitive k sits next to lane iab+1's.
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: pcut
      !! Cutoff for the primitive-quartet prescreen: a primitive quartet
      !! whose prefactor cannot reach this is skipped before the Boys
      !! function. It is a fraction of `thresh`, and a runtime argument
      !! rather than a literal so the fraction can be measured against the
      !! energy it costs instead of guessed once.
      !!
      !! Tested against the UNATTENUATED prefactor (the primitive loop below
      !! computes `pref` before `theta`, deliberately): erf(omega*r12)/r12 is
      !! pointwise smaller than 1/r12, so the Coulomb prefactor bounds the
      !! attenuated integral too, just loosely at omega > 0. Still a valid
      !! prescreen, never a missed quartet -- see trc_screen.F90 for the
      !! Schwarz bound's identical reasoning.
      real(dp), intent(in)    :: jfac, kfac
      !! Coulomb and exchange scalings. Per call, so they fold into the six
      !! atomic updates and the digestion stays folded -- separating J and K
      !! into two matrices to scale them would give back FOCK6's 22%.
      real(dp), intent(in)    :: omega
      !! Range-separation parameter, bohr^-1. Zero (the default every caller
      !! but wB97X-V's long-range exchange pass uses) is today's Coulomb
      !! operator; omega > 0 attenuates it to erf(omega*r12)/r12. A plain
      !! scalar, not optional, so it costs one comparison per primitive
      !! quartet and nothing else -- see the primitive loop below, where
      !! omega == 0 makes every quantity it touches bit-identical to the
      !! unattenuated value.
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: nqp
      !! Bin-interleaved primitive-pair data: qp_z/qp_c/qp_px/qp_py/qp_pz
      !! are what the kernel reads for the bra/ket primitives, at
      !! sPA(seg) + (k-1)*sA(seg) + iab (bra) or sPB(seg) + (k-1)*sNB(seg)
      !! + icd (ket).
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      !! Per pair slot, aligned with sp_i/sp_j/sp_q: the A centre, A - B,
      !! and the real kept primitive-pair count pp_n(key) gave before.
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: sq_n(npair)
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
    # LARGE T, NO TRANSCENDENTALS.
    #
    # Past BOYS_TMAX = 30, erf(sqrt(T)) is 1 to within 1e-14 and exp(-T) is
    # below 1e-13, so the asymptotic form F_m(T) = (2m-1)!!/(2T)^m
    # sqrt(pi/(4T)) is exact to the last digit that matters -- the prescreen
    # margin already puts 3e-7 into the Fock matrix. This branch used to
    # call erf AND exp, and a warp pays for it if ONE of its lanes takes it.
    # Deep contraction guarantees that: the tight primitives of a cc-pVDZ s
    # shell push T far past 30 while the diffuse ones stay under it, so the
    # classes that hurt most were the ones diverging into two transcendentals
    # per primitive quartet. What is left is a reciprocal square root.
    w("               if (tval >= BOYS_TMAX) then")
    w("                  btt = 1.0_dp/tval")
    w("                  f(0) = 0.88622692545275801365_dp*sqrt(btt)")
    w("                  bet = 0.0_dp")
    # bet is zero here, and 1/T is already in hand, so this is the
    # asymptotic (2k-1)!!/(2T)^k ladder with no divisions of its own.
    for k in range(1, m + 1):
        w(f"                  f({k}) = {float(2*k - 1)}_dp*f({k - 1})*(0.5_dp*btt)")
    w("               else")
    w("                  bi = int(tval*BOYS_DTINV)")
    w("                  bi = min(bi, BOYS_NGRID - 1)")
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


def _unroll_combos(txt, ncab):
    """The per-combination loop, written out at literal column pairs.

    The accumulator is indexed by (x, qab, qcd) and x alone has a constant
    bound; with qab and qcd runtime the whole block goes to local memory,
    which is why this used to copy it there once and read the copy. At a
    literal column pair the block stays in registers and the copy is gone,
    which is what pays for a wider chunk. The three canonical filters were
    `cycle` statements and become a flag, since there is no loop left."""
    head = "         do qcd = 1, ncd_c\n         do qab = 1, nab_c\n"
    tail = "         end do   ! qab\n         end do   ! qcd\n"
    i = txt.index(head)
    j = txt.index(tail, i)
    body = txt[i + len(head):j]
    out = ["         ! === combinations, unrolled at literal column pairs ===\n"]
    for qcd in range(1, ncab + 1):
        for qab in range(1, ncab + 1):
            b = body
            for old, new in (
                ("iabc = ab0 + qab - 1", f"iabc = ab0 + {qab} - 1"),
                ("icdc = cd0 + qcd - 1", f"icdc = cd0 + {qcd} - 1"),
                ("         if (.not. ok) cycle\n", "         if (ok) then\n"),
                ("g1(x) = g(x, qab, qcd)", f"g1(x) = g(x, {qab}, {qcd})")):
                assert b.count(old) == 1, old
                b = b.replace(old, new)
            out.append(f"         if ({qab} <= nab_c .and. {qcd} <= ncd_c) then\n"
                       + b + "         end if\n         end if\n")
    return txt[:i] + "".join(out) + txt[j + len(tail):]


#: ---------------------------------------------------------------------
#: SCAFFOLDING SHARED WITH THE ROTATED-AXIS EMITTER.
#:
#: `scripts/rotaxis/emit_docon.py` calls these same functions rather than
#: keeping a second copy: the segment locate/decode, the exact Schwarz and
#: density-weighted screens, the same_ab/same_cd/same_pair canonicalisation,
#: the pair-table offset lookup, the scalar-path degeneracy/AO-offset decode,
#: and the block-accumulated, density-batched digestion do not depend on WHICH integral core fills
#: `vbuf` -- HGP's VRR/HRR or rotaxis's frame-and-three-levels -- only on the
#: pair tables and the four shells' Cartesian counts. Factoring them out is
#: what makes "byte-identical before and after" a meaningful check: this
#: file's own emission is unchanged (same text, now produced by a function
#: call instead of an inline literal), and the rotaxis emitter gets the
#: identical text rather than a hand-copied approximation of it.
#: ---------------------------------------------------------------------

def emit_segment_decode(bra_fastest=True):
    """Locate the segment `gt` falls in and decode (iab, icd) within it.

    `bra_fastest` (the default, and the ONLY value HGP's own kernels ever
    pass -- byte-identity depends on that) makes 32 consecutive threads hold
    32 consecutive bra pairs against one ket pair: the right mapping when
    the PRIMITIVE loop nested inside is ket-outside/bra-inside, i.e. HGP's,
    because then the innermost loads -- the ket side -- are the same address
    across the warp (broadcast). See the KET-UNIFORM DECODE comment this
    emits for the measurement behind it.

    `bra_fastest=False` mirrors it for rotaxis's OWN primitive nesting,
    which FRAME.md fixes the other way round (ket outside, bra inside, so
    the level-1 accumulators can be per-bra-primitive and component-
    independent -- README_TERCO.md section 2.4). With the same thread
    mapping as HGP, rotaxis's innermost loop -- bra primitives -- would be
    the one that diverges instead of broadcasting, undoing the whole point.
    Swapping which pair index is uniform across the warp fixes that: 32
    consecutive threads then hold 32 consecutive KET pairs against one BRA
    pair, so the bra-primitive loop's addresses (`offab`, fixed for the
    warp) are identical across every lane.

    The symmetric-segment closed form (`sD(seg)`) enumerates the SAME set of
    (iab, icd) pairs with iab >= icd either way -- verified by exhaustive
    count against the bra-fastest form for every segment size up to several
    thousand pairs and spot-checked by round-trip up to 65000 -- just row by
    iab (mirrored) instead of row by icd. The non-symmetric form swaps which
    of `sA` (bra pairs) / `sNB` (ket pairs) is the modulus."""
    o = ["""      ! locate the segment
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
"""]
    if bra_fastest:
        o.append("""
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
            nsa = int(sA(seg), int64)
            u = int(t - 1, int64)
            kx = int((real(2*nsa + 1, dp) - sqrt(real(2*nsa + 1, dp)**2 - 8.0_dp*real(u, dp)))/2.0_dp, int64)
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
""")
    else:
        o.append("""
         ! BRA-UNIFORM DECODE (mirrored from HGP's KET-UNIFORM form above --
         ! see emit_segment_decode's docstring). The ket pair index runs
         ! fastest, so the 32 threads of a warp hold 32 ket pairs against ONE
         ! bra pair, and rotaxis's innermost primitive loop -- the bra side,
         ! FRAME.md's ket-outside/bra-inside nesting -- reads the same
         ! address across the warp: one L1 transaction, broadcast. The ket
         ! loads diverge instead, once per ket primitive rather than once per
         ! quartet -- the outer loop, run far fewer times.
         ! On a symmetric segment the pairs are enumerated row by row in
         ! iab (mirrored from HGP's row-by-icd form), iab >= icd throughout,
         ! with the standard triangular-number inversion.
         if (sD(seg)) then
            nsa = int(sA(seg), int64)
            u = int(t - 1, int64)
            kx = int((sqrt(8.0_dp*real(u, dp) + 1.0_dp) - 1.0_dp)/2.0_dp, int64)
            if (kx < 0) kx = 0
            do while (kx > 0)
               if (kx*(kx + 1)/2 <= u) exit
               kx = kx - 1
            end do
            do while ((kx + 1)*(kx + 2)/2 <= u)
               kx = kx + 1
            end do
            iab = int(kx) + 1
            icd = int(u - kx*(kx + 1)/2) + 1
         else
            iab = (t - 1)/sNB(seg) + 1
            icd = t - (iab - 1)*sNB(seg)
         end if
""")
    return "".join(o)


def emit_screens():
    """Exact Schwarz and density-weighted screens, and the same_ab/same_cd/
    same_pair canonicalisation flags. Fixed text, shared for the same reason
    as `emit_segment_decode`."""
    return """
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
"""


def emit_pair_offsets(with_columns):
    """Pair-table offsets for the bra and ket shell pairs.

    `with_columns` additionally decodes the per-shell column/primitive
    counts the general-contraction (block) kernel needs; the scalar HGP path
    and every rotaxis kernel -- one column, no block accumulation -- pass
    False."""
    o = """
         keyab = (si - 1)*nbas + sj
         keycd = (sk - 1)*nbas + sl
         offab = pp_off(keyab); nab = pp_n(keyab)
         offcd = pp_off(keycd); ncd = pp_n(keycd)
"""
    if with_columns:
        o += "         nca = ps_ncol(si); ncb = ps_ncol(sj); nccl = ps_ncol(sk); ncdl = ps_ncol(sl)\n"
        o += "         npi = ps_np(si); npj = ps_np(sj); npk = ps_np(sk); npl = ps_np(sl)\n"
    return o


def emit_scalar_decode():
    """Degeneracy flags and AO offsets on the single-column (scalar) path:
    the canonical enumeration and the degeneracy weights are the pair-level
    ones directly, with no column decode. Identical for HGP and rotaxis."""
    return """         ! Single column on every side: the canonical enumeration and the
         ! degeneracy weights are the pair-level ones.
         dij = .not. same_ab
         dkl = .not. same_cd
         dpq = .not. same_pair
         mui = ao_off(si); nuj = ao_off(sj); lamk = ao_off(sk); sigl = ao_off(sl)

"""


#: First line of every emitted digestion; gen_rys.py and the block-form
#: splice below cut the kernel text there.
DIGEST_MARKER = "         ! BLOCK-ACCUMULATED DIGESTION.\n"


def _do(var, n, ind):
    """Open the loop `do var = 0, n - 1` at indent `ind`; for a single-trip
    loop (n == 1, an s shell) assign `var = 0` instead and open nothing."""
    return f"{ind}do {var} = 0, {n - 1}\n" if n > 1 else f"{ind}{var} = 0\n"


def _end(n, ind):
    """Close what `_do(var, n, ind)` opened: nothing for n == 1."""
    return f"{ind}end do\n" if n > 1 else ""


def _jblocks_load(na, nb, nc, nd):
    """The general branch's Coulomb gathers under `jfac /= 0`, zero-filled
    otherwise. Emitted without single-trip loops, as `_jonly_digestion` is:
    with them nvfortran 26.5's fort2 segfaulted on classes with two s shells
    at some build-path lengths (Gadi, trc_pc_rys0201/0210/1021)."""
    i9, i12, i15 = " " * 9, " " * 12, " " * 15
    o = f"{i9}if (jfac /= 0.0_dp) then\n"
    o += _do("ib", nb, i9) + _do("ia", na, i12)
    o += (f"{i15}dab(1 + ia + {na}*ib) = dmat(idens, mui + ia, nuj + ib)\n"
          f"{i15}jab(1 + ia + {na}*ib) = 0.0_dp\n")
    o += _end(na, i12) + _end(nb, i9)
    o += _do("id", nd, i9) + _do("ic", nc, i12)
    o += (f"{i15}dcd(1 + ic + {nc}*id) = dmat(idens, lamk + ic, sigl + id)\n"
          f"{i15}jcd(1 + ic + {nc}*id) = 0.0_dp\n")
    o += _end(nc, i12) + _end(nd, i9)
    o += f"{i9}else\n"
    o += _do("ib", nb, i9) + _do("ia", na, i12)
    o += (f"{i15}dab(1 + ia + {na}*ib) = 0.0_dp\n"
          f"{i15}jab(1 + ia + {na}*ib) = 0.0_dp\n")
    o += _end(na, i12) + _end(nb, i9)
    o += _do("id", nd, i9) + _do("ic", nc, i12)
    o += (f"{i15}dcd(1 + ic + {nc}*id) = 0.0_dp\n"
          f"{i15}jcd(1 + ic + {nc}*id) = 0.0_dp\n")
    o += _end(nc, i12) + _end(nd, i9)
    o += f"{i9}end if\n"
    return o


def _jblocks_store(na, nb, nc, nd):
    """The general branch's Coulomb atomics, skipped when jfac = 0; no
    single-trip loops (see _jblocks_load)."""
    i9, i12, i15 = " " * 9, " " * 12, " " * 15
    o = f"{i9}if (jfac /= 0.0_dp) then\n"
    o += _do("ib", nb, i9) + _do("ia", na, i12)
    o += (f"{i15}!$acc atomic update\n"
          f"{i15}jmat(idens, mui + ia, nuj + ib) = jmat(idens, mui + ia, nuj + ib) &\n"
          f"{' ' * 42}+ jab(1 + ia + {na}*ib)\n")
    o += _end(na, i12) + _end(nb, i9)
    o += _do("id", nd, i9) + _do("ic", nc, i12)
    o += (f"{i15}!$acc atomic update\n"
          f"{i15}jmat(idens, lamk + ic, sigl + id) = jmat(idens, lamk + ic, sigl + id) &\n"
          f"{' ' * 44}+ jcd(1 + ic + {nc}*id)\n")
    o += _end(nc, i12) + _end(nd, i9)
    o += f"{i9}end if\n"
    return o


def _jonly_digestion(na, nb, nc, nd):
    """The kfac = 0 branch of emit_digestion, from its density loop through
    the J scatters. Its statements are the general branch's J statements,
    but it emits no single-trip loop: an s shell's `do ia = 0, 0` becomes
    `ia = 0`, so every index and every floating-point operation is the one
    the loop would have run.

    Why: with those loops, nvfortran 26.5's fort2 fails on the classes with
    two or more s shells -- signal 11 in acc_remove_code, or the internal
    error "ili_hash_func: invalid opc 65535" -- while compiling the
    `routine seq` item routine. Whether a given file fails depends on
    fort2's heap layout, not on its code: the same file compiles or not as
    the length of its absolute path changes (a build-directory name is
    enough; a Perlmutter LMAX=3 MPI build lost (10|10) that way), and the
    include list, the mpif90 wrapper and TRC_LMAX move the failing lengths.
    With LMAX=3 flags, (00|00), (00|10), (01|00), (10|00) and (10|10) each
    failed for some path lengths between 45 and 200; at LMAX=2, (00|01),
    (01|00), (00|10) and (10|00) failed between 60 and 130. With no
    single-trip loop in this branch, no class failed at any of those
    lengths and configurations, and 26.9 compiles either form. The general
    branch keeps its loops: without this branch, it failed in none of them
    either."""
    i9, i12, i15 = " " * 9, " " * 12, " " * 15
    o = f"{i9}do idens = 1, ndens\n\n"
    o += _do("ib", nb, i9) + _do("ia", na, i12)
    o += (f"{i15}dab(1 + ia + {na}*ib) = dmat(idens, mui + ia, nuj + ib)\n"
          f"{i15}jab(1 + ia + {na}*ib) = 0.0_dp\n")
    o += _end(na, i12) + _end(nb, i9)
    o += _do("id", nd, i9) + _do("ic", nc, i12)
    o += (f"{i15}dcd(1 + ic + {nc}*id) = dmat(idens, lamk + ic, sigl + id)\n"
          f"{i15}jcd(1 + ic + {nc}*id) = 0.0_dp\n")
    o += _end(nc, i12) + _end(nd, i9)
    o += """
         ! Spelled differently from the general branch, on purpose, and with
         ! the same arithmetic: the component index in closed form and the
         ! scaled integral inline rather than in `idx` and `sc`. A second
         ! counted loop fails to translate the Rys launch under nvfortran
         ! 26.5 ("Missing increment value", trc_pc_rys1102).
"""
    o += _do("id", nd, i9) + _do("ic", nc, i9) + _do("ib", nb, i9) + _do("ia", na, i9)
    vb = f"(wq*vbuf(1 + ia + {na}*ib + {na*nb}*ic + {na*nb*nc}*id))"
    o += (f"{i12}jab(1 + ia + {na}*ib) = jab(1 + ia + {na}*ib) &\n"
          f"{' ' * 36}+ 4.0_dp*jfac*{vb} &\n"
          f"{' ' * 36}*dcd(1 + ic + {nc}*id)\n"
          f"{i12}jcd(1 + ic + {nc}*id) = jcd(1 + ic + {nc}*id) &\n"
          f"{' ' * 36}+ 4.0_dp*jfac*{vb} &\n"
          f"{' ' * 36}*dab(1 + ia + {na}*ib)\n")
    o += _end(na, i9) + _end(nb, i9) + _end(nc, i9) + _end(nd, i9)
    o += "\n"
    o += _do("ib", nb, i9) + _do("ia", na, i12)
    o += (f"{i15}!$acc atomic update\n"
          f"{i15}jmat(idens, mui + ia, nuj + ib) = jmat(idens, mui + ia, nuj + ib) &\n"
          f"{' ' * 42}+ jab(1 + ia + {na}*ib)\n")
    o += _end(na, i12) + _end(nb, i9)
    o += _do("id", nd, i9) + _do("ic", nc, i12)
    o += (f"{i15}!$acc atomic update\n"
          f"{i15}jmat(idens, lamk + ic, sigl + id) = jmat(idens, lamk + ic, sigl + id) &\n"
          f"{' ' * 44}+ jcd(1 + ic + {nc}*id)\n")
    o += _end(nc, i12) + _end(nd, i9)
    o += "\n"
    return o


def emit_digestion(na, nb, nc, nd):
    """Block-accumulated, density-batched Fock digestion from `vbuf`.
    Its second line is DIGEST_MARKER, which the splicing generators search
    for. Parameterised only by the four shells' Cartesian counts, so it is the same text for every kernel -- HGP or
    rotaxis -- that fills `vbuf(1:na*nb*nc*nd)` in the ia-fastest order this
    reads (FRAME.md section 6 confirms rotaxis's derivation order already
    matches it, so no remapping is needed to reuse this)."""
    return f"""         !
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
         ! J-ONLY DIGESTION, for kfac = 0 (a pure LDA/GGA functional, or the
         ! Coulomb-only pass of the response equations). The exchange half of
         ! the general digestion below is then four density-block gathers and
         ! four multiply-adds per component that all come out exactly
         ! 0.0_dp, plus four accumulator blocks held live for nothing. This
         ! copy drops them; its J statements are the general branch's,
         ! character for character, so every jmat element is bit-identical
         ! to what that branch would have produced. A separate branch and not
         ! a test inside the component loop: guarding the K updates per
         ! component was measured 22% slower on a hybrid (w20/6-31G* B3LYP,
         ! V100), because the test sits on the branch that is always taken.
         ! kfac is uniform over the launch, so no warp diverges here.
         !
         if (kfac == 0.0_dp) then
{_jonly_digestion(na, nb, nc, nd)}         end do   ! idens
         else

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

         !
         ! jfac = 0 whenever the MD J engine forms J (always for a g basis):
         ! the Coulomb blocks are then zero-filled rather than gathered and
         ! their atomics skipped below. The component loop is left alone (a
         ! test inside it costs more than it saves, see the J-only note).
         ! jfac is uniform over the launch, so no warp diverges.
         !
{_jblocks_load(na, nb, nc, nd)}         do ic = 0, {nc - 1}
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

{_jblocks_store(na, nb, nc, nd)}         do ic = 0, {nc - 1}
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
         end if
"""


def emit_scalar_driver(name, item):
    """The scalar launch driver: `do concurrent` over the quartets of one
    segment run, ranks split by static stride, the item routine called once
    per iteration. Fixed text parameterised only by the driver's and item
    routine's own names, so HGP and rotaxis kernels share it exactly.

    HGP's own scalar driver is still produced by editing the general-
    contraction (block) driver's text (see `emit_kernel` below) rather than
    by calling this function: that derivation is what `--split`'s byte-
    identity check was verified against, and changing it over risked that
    invariant for no behavioural gain, since the two are already the same
    text by construction. `scripts/rotaxis/emit_docon.py` -- which has no
    block-form driver to derive from -- calls this directly."""
    return f"""   subroutine {name}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                      npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                      nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                      sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                      ndens, dmat, jmat, rank, nranks)
{PROLOGUE}
      integer, intent(in) :: rank, nranks
      integer(int64) :: g0, g1, nr, i
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
         call {item}(g0 + (i - 1)*nranks, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                       npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                       nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                       sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat)
      end do
   end subroutine {name}
"""


def emit_item_signature(name):
    """The item routine's signature: same argument list as the driver plus
    the flat quartet index `gt`, `pure` and `!$acc routine seq` so it inlines
    into its caller (the reason it lives in the same module and file)."""
    return f"""   pure subroutine {name}(gt, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                      npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                      nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                      sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                      ndens, dmat, jmat)
      !$acc routine seq
      integer(int64), intent(in) :: gt
{PROLOGUE}
"""


def _emit_block(la, lb, lc, ld, cidx, vrr_body, hrr_body, unroll=True):
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
    # REGISTERS, NOT FLOPS, SET THIS. The accumulator costs nv*ncab^2
    # doubles per thread beside the VRR block's 2*nv and g1's nv; ptxas has
    # 255 registers, 127 doubles, and vbuf and the six digestion blocks want
    # their share. Budget the accumulator at 64 doubles and take the largest
    # power-of-two chunk that fits, capped at PS_NCOL_MAX columns a side.
    # The old rule (4 at nv<=1, 2 at nv<=4, 1 otherwise) left every class
    # with a p on both sides at one combination per pass, which is the
    # segmented cost with extra bookkeeping -- and those are exactly the
    # classes cc-pVDZ spends its time in.
    ncab = 4
    while ncab > 1 and nv*ncab*ncab > 64:
        ncab //= 2
    # The block accumulation is written out with literal indices so the
    # accumulator is provably register-resident; a loop over (qab, qcd) is
    # left to the unroller's judgement, which it lost on this kernel.
    acc_zero = "".join(
        f"            do x = 1, {nv}\n               g(x, {a}, {c}) = 0.0_dp\n            end do\n"
        for c in range(1, ncab + 1) for a in range(1, ncab + 1))
    cab_fill = "".join(
        f"               cab({a}) = wta({a})*ps_coef(offa({a}) + ki)*ps_coef(offb({a}) + kj)\n"
        for a in range(1, ncab + 1))
    # The prescreen's bound needs the largest coefficient the chunk can
    # apply, since the pair factor on this path carries none of them.
    cab_fill += ("               cabmax = " + ("abs(cab(1))" if ncab == 1 else
                 "max(" + ", ".join(f"abs(cab({a}))" for a in range(1, ncab + 1)) + ")") + "\n")
    ccd_fill = "".join(
        f"                  ccd({c}) = wtc({c})*ps_coef(offc({c}) + kk)*ps_coef(offd({c}) + kl)\n"
        for c in range(1, ncab + 1))
    ccd_fill += ("                  ccdmax = " + ("abs(ccd(1))" if ncab == 1 else
                 "max(" + ", ".join(f"abs(ccd({c}))" for c in range(1, ncab + 1)) + ")") + "\n")
    # The column-pair weight of every combination, formed once per
    # primitive quartet, and the accumulation straight into the block. The
    # ket-side partial this used to carry saved ncab^2 - ncab FMAs per
    # quartet and cost nv*ncab registers, which at ncab = 2 is the wrong
    # side of the trade: registers are what stop the chunk being wider.
    w_fill = "".join(
        f"                  w2({a}, {c}) = cab({a})*ccd({c})\n"
        for c in range(1, ncab + 1) for a in range(1, ncab + 1))
    g_add = "".join(
        f"                  do x = 1, {nv}\n                     g(x, {a}, {c}) = g(x, {a}, {c}) + w2({a}, {c})*v(x, cur)\n                  end do\n"
        for c in range(1, ncab + 1) for a in range(1, ncab + 1))
    hrr_body = re.sub(r"\bg\(", "g1(", hrr_body)

    txt = f"""
   subroutine pc{tag}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                      npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                      nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                      sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                      ndens, dmat, jmat, rank, nranks)
      !! ({la}{lb}|{lc}{ld}) driver.  The `do concurrent` lives here and the
      !! workspaces live in the item routine below, so they are per THREAD.
      !! Declaring them alongside the loop makes them shared -- the compiler then
      !! emits `implicit copy(v, g, vbuf, f)` per launch and the threads race.
{PROLOGUE}
      integer, intent(in) :: rank, nranks
      integer(int64) :: g0, g1, nr, i
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
         call pci{tag}(g0 + (i - 1)*nranks, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                       npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                       nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                       sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat)
      end do
   end subroutine pc{tag}

   pure subroutine pci{tag}(gt, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                      npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                      nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                      sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                      ndens, dmat, jmat)
      !$acc routine seq
      integer(int64), intent(in) :: gt
{PROLOGUE}
      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl
      integer(int64) :: nsa, u, kx
      real(dp) :: qcut, bnd
      integer :: keyab, keycd, offab, offcd, nab, ncd
      integer :: kp, kq, d, x, cur, ia, ib, ic, id, idx, idens
      integer :: mu, nu, lam, sig, mui, nuj, lamk, sigl
      logical :: dij, dkl, dpq, ok
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc
      real(dp) :: theta, rho_theta
      real(dp) :: pqx, pqy, pqz, pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re, sc, vv
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      real(dp) :: rpx, rpy, rpz, rzeta, reta, rzpe
      real(dp) :: f(0:BOYS_MMAX)
      integer  :: bi, bj, bbase
      real(dp) :: bx, bx2, b0, b1, b2, btt, bet
      real(dp) :: v({nv}, 0:1), g({nv}, {ncab}, {ncab}), g1({nv}), vbuf({na*nb*nc*nd})
      real(dp) :: wq, w, wab, cabmax, ccdmax, qchunk
      integer  :: nca, ncb, nccl, ncdl, npi, npj, npk, npl, ncab, nccd, ab0, cd0, nab_c, ncd_c, qab, qcd
      integer  :: ki, kj, kk, kl, kpl, kql, ia2, ib2, ic2, id2, iabc, icdc
      logical  :: same_ab, same_cd, same_pair
      real(dp) :: cab({ncab}), ccd({ncab}), wta({ncab}), wtc({ncab}), w2({ncab}, {ncab})
      integer  :: offa({ncab}), offb({ncab}), offc({ncab}), offd({ncab})
      real(dp) :: jab({na*nb}), jcd({nc*nd}), kac({na*nc})
      real(dp) :: kad({na*nd}), kbc({nb*nc}), kbd({nb*nd})
      real(dp) :: dab({na*nb}), dcd({nc*nd}), dac({na*nc})
      real(dp) :: dad({na*nd}), dbc({nb*nc}), dbd({nb*nd})

{emit_segment_decode()}{emit_screens()}{emit_pair_offsets(True)}
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
               ki = pp_ki(kp); kj = pp_kj(kp)
               wab = ps_coef(ps_coff(si) + ki)*ps_coef(ps_coff(sj) + kj)
               ! Bra-only quantities, out of the ket loop. They depend on
               ! kp alone and were being reloaded and recomputed once per
               ! KET primitive, which on a deeply contracted shell pair is
               ! 144 times over for values that never change. The compiler
               ! does not lift them, presumably because it cannot prove the
               ! loads invariant. Measured on (ps|ss) of the silica slice in
               ! cc-pVDZ: 3.21 s -> 2.29 s, with every other class in the
               ! same profile flat.
               rpx = pp_r(kp, 1); rpy = pp_r(kp, 2); rpz = pp_r(kp, 3)
               pax = rpx - pp_ra(kp, 1)
               pay = rpy - pp_ra(kp, 2)
               paz = rpz - pp_ra(kp, 3)
               ! 1/zeta with them: the bra exponent does not change across
               ! the ket loop, so its reciprocal is one division per bra
               ! primitive rather than three per primitive QUARTET.
               rzeta = 1.0_dp/zeta
               oo2z = 0.5_dp*rzeta
               do kq = offcd + 1, offcd + ncd
                  eta = pp_p(kq)
                  kk = pp_ki(kq); kl = pp_kj(kq)
                  w = wab*ps_coef(ps_coff(sk) + kk)*ps_coef(ps_coff(sl) + kl)
                  zpe = zeta + eta
                  ! PRIMITIVE-QUARTET PRESCREEN. The prefactor bounds the
                  ! primitive (ss|ss) integral, and with the normalisation
                  ! in the coefficients it bounds the higher ones to within
                  ! the polynomial factors the cutoff's three decades of
                  ! margin cover. Tested before the Boys function and the
                  ! VRR, which is nearly all of a primitive quartet's cost;
                  ! on a generally contracted basis two thirds of the
                  ! quartets that survive the pair pruning die here.
                  pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)
                  if (abs(pref*w) <= pcut) cycle
                  ! One reciprocal each for eta and zeta+eta, then multiply.
                  ! This loop used to issue ten double-precision divisions
                  ! per primitive quartet -- rho, three for the W centre,
                  ! three halves, two rho ratios and the prefactor -- against
                  ! a recurrence that for the light classes is five lines.
                  ! Division has no fast reciprocal in double precision, so
                  ! that was the arithmetic, not the recurrence.
                  reta = 1.0_dp/eta
                  rzpe = 1.0_dp/zpe
                  rho = zeta*eta*rzpe
                  ! RANGE SEPARATION. erf(omega*r12)/r12 is the Coulomb
                  ! integral with rho replaced everywhere by
                  ! rho_theta = theta*rho, theta = omega^2/(rho + omega^2) --
                  ! the W centre (wpx/wqx below), the Boys argument tval, and
                  ! the oo2ze/rz/re VRR coefficients -- and the prefactor
                  ! picks up a sqrt(theta) factor. oo2z and oo2e come from
                  ! zeta and eta alone and are not touched: the kernels do
                  ! NOT use rho consistently, and those two are the
                  ! quantities that only equal a rho expression through
                  ! rho*(zeta+eta) = zeta*eta, not through rho itself.
                  !
                  ! omega == 0 sets theta = 1.0_dp exactly, so
                  ! rho_theta equals rho bit-for-bit and every "theta*" below
                  ! is an exact multiply by one: the attenuated path costs
                  ! this one branch and is otherwise identical to the
                  ! Coulomb result, to the last bit.
                  if (omega > 0.0_dp) then
                     theta = omega*omega/(rho + omega*omega)
                  else
                     theta = 1.0_dp
                  end if
                  rho_theta = theta*rho
                  pqx = rpx - pp_r(kq, 1)
                  pqy = rpy - pp_r(kq, 2)
                  pqz = rpz - pp_r(kq, 3)
                  qcx = pp_r(kq, 1) - pp_ra(kq, 1)
                  qcy = pp_r(kq, 2) - pp_ra(kq, 2)
                  qcz = pp_r(kq, 3) - pp_ra(kq, 3)
                  wc = (zeta*rpx + eta*pp_r(kq, 1))*rzpe
                  wpx = wc - pp_r(kp, 1); wqx = wc - pp_r(kq, 1)
                  wc = (zeta*rpy + eta*pp_r(kq, 2))*rzpe
                  wpy = wc - pp_r(kp, 2); wqy = wc - pp_r(kq, 2)
                  wc = (zeta*rpz + eta*pp_r(kq, 3))*rzpe
                  wpz = wc - pp_r(kp, 3); wqz = wc - pp_r(kq, 3)
                  wpx = theta*wpx; wpy = theta*wpy; wpz = theta*wpz
                  wqx = theta*wqx; wqy = theta*wqy; wqz = theta*wqz
                  tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)

   {emit_boys(lt)}

                  pref = pref*sqrt(theta)
                  oo2e = 0.5_dp*reta; oo2ze = theta*(0.5_dp*rzpe)
                  rz = rho_theta*rzeta; re = rho_theta*reta

   {vrr_body}
                  do x = 1, {nv}
                     g1(x) = g1(x) + w*v(x, cur)
                  end do
               end do
            end do
            do x = 1, {nv}
               g(x, 1, 1) = g1(x)
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
            ! Nothing in this chunk can contribute: skip it whole, which
            ! is the only place the shared VRR can be skipped at all.
            if (nqc > 1) then
               qchunk = 0.0_dp
               do qab = 1, nab_c
               do qcd = 1, ncd_c
                  qchunk = max(qchunk, &
                     cb{tag}(ab0 + qab - 1, cd0 + qcd - 1, nca, nccl, si, sj, sk, sl, nbas, ncoltot, ps_soff, col_sh, nshc, nqc, q_col, dsh_c))
               end do
               end do
               if (qchunk <= thresh) cycle
            end if
            do kp = offab + 1, offab + nab
               zeta = pp_p(kp)
               ki = pp_ki(kp); kj = pp_kj(kp)
{cab_fill}               ! Bra-only quantities, out of the ket loop. They depend on
               ! kp alone and were being reloaded and recomputed once per
               ! KET primitive, which on a deeply contracted shell pair is
               ! 144 times over for values that never change. The compiler
               ! does not lift them, presumably because it cannot prove the
               ! loads invariant. Measured on (ps|ss) of the silica slice in
               ! cc-pVDZ: 3.21 s -> 2.29 s, with every other class in the
               ! same profile flat.
               rpx = pp_r(kp, 1); rpy = pp_r(kp, 2); rpz = pp_r(kp, 3)
               pax = rpx - pp_ra(kp, 1)
               pay = rpy - pp_ra(kp, 2)
               paz = rpz - pp_ra(kp, 3)
               ! 1/zeta with them: the bra exponent does not change across
               ! the ket loop, so its reciprocal is one division per bra
               ! primitive rather than three per primitive QUARTET.
               rzeta = 1.0_dp/zeta
               oo2z = 0.5_dp*rzeta
               do kq = offcd + 1, offcd + ncd
                  eta = pp_p(kq)
                  kk = pp_ki(kq); kl = pp_kj(kq)
{ccd_fill}                  zpe = zeta + eta
                  !
                  ! PRIMITIVE-QUARTET PRESCREEN. The pair factor on this
                  ! path carries no contraction coefficients -- they are
                  ! applied per column combination -- so the bound takes the
                  ! largest of them on each side. Without it this kernel
                  ! evaluated every primitive quartet in the pair list while
                  ! the scalar one skipped two thirds of them on cc-pVDZ,
                  ! which is most of why the primitive-shell view measured
                  ! slower than the path it was meant to beat.
                  !
                  pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)
                  if (abs(pref)*cabmax*ccdmax <= pcut) cycle
                  ! One reciprocal each for eta and zeta+eta, then multiply.
                  ! This loop used to issue ten double-precision divisions
                  ! per primitive quartet -- rho, three for the W centre,
                  ! three halves, two rho ratios and the prefactor -- against
                  ! a recurrence that for the light classes is five lines.
                  ! Division has no fast reciprocal in double precision, so
                  ! that was the arithmetic, not the recurrence.
                  reta = 1.0_dp/eta
                  rzpe = 1.0_dp/zpe
                  rho = zeta*eta*rzpe
                  ! RANGE SEPARATION. erf(omega*r12)/r12 is the Coulomb
                  ! integral with rho replaced everywhere by
                  ! rho_theta = theta*rho, theta = omega^2/(rho + omega^2) --
                  ! the W centre (wpx/wqx below), the Boys argument tval, and
                  ! the oo2ze/rz/re VRR coefficients -- and the prefactor
                  ! picks up a sqrt(theta) factor. oo2z and oo2e come from
                  ! zeta and eta alone and are not touched: the kernels do
                  ! NOT use rho consistently, and those two are the
                  ! quantities that only equal a rho expression through
                  ! rho*(zeta+eta) = zeta*eta, not through rho itself.
                  !
                  ! omega == 0 sets theta = 1.0_dp exactly, so
                  ! rho_theta equals rho bit-for-bit and every "theta*" below
                  ! is an exact multiply by one: the attenuated path costs
                  ! this one branch and is otherwise identical to the
                  ! Coulomb result, to the last bit.
                  if (omega > 0.0_dp) then
                     theta = omega*omega/(rho + omega*omega)
                  else
                     theta = 1.0_dp
                  end if
                  rho_theta = theta*rho
                  pqx = rpx - pp_r(kq, 1)
                  pqy = rpy - pp_r(kq, 2)
                  pqz = rpz - pp_r(kq, 3)
                  qcx = pp_r(kq, 1) - pp_ra(kq, 1)
                  qcy = pp_r(kq, 2) - pp_ra(kq, 2)
                  qcz = pp_r(kq, 3) - pp_ra(kq, 3)
                  wc = (zeta*rpx + eta*pp_r(kq, 1))*rzpe
                  wpx = wc - pp_r(kp, 1); wqx = wc - pp_r(kq, 1)
                  wc = (zeta*rpy + eta*pp_r(kq, 2))*rzpe
                  wpy = wc - pp_r(kp, 2); wqy = wc - pp_r(kq, 2)
                  wc = (zeta*rpz + eta*pp_r(kq, 3))*rzpe
                  wpz = wc - pp_r(kp, 3); wqz = wc - pp_r(kq, 3)
                  wpx = theta*wpx; wpy = theta*wpy; wpz = theta*wpz
                  wqx = theta*wqx; wqy = theta*wqy; wqz = theta*wqz
                  tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)

   {emit_boys(lt)}

                  pref = pref*sqrt(theta)
                  oo2e = 0.5_dp*reta; oo2ze = theta*(0.5_dp*rzpe)
                  rz = rho_theta*rzeta; re = rho_theta*reta

   {vrr_body}
{w_fill}{g_add}               end do
            end do
         end if

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
         ok = .not. (same_ab .and. ia2 < ib2)
         if (same_cd .and. ic2 < id2) ok = .false.
         if (same_pair .and. iabc < icdc) ok = .false.
         ! ... and only those whose own column bound survives. The quartet
         ! test upstream used the merged maximum over the columns, which is
         ! the right admission test for the quartet and much too generous
         ! for any one combination of it.
         if (nqc > 1) then
            if (cb{tag}(iabc, icdc, nca, nccl, si, sj, sk, sl, nbas, ncoltot, ps_soff, col_sh, nshc, nqc, q_col, dsh_c) <= thresh) ok = .false.
         end if
         if (.not. ok) cycle
         dij = .not. (same_ab .and. ia2 == ib2)
         dkl = .not. (same_cd .and. ic2 == id2)
         dpq = .not. (same_pair .and. iabc == icdc)
         mui = col_ao(ps_soff(si) + ia2); nuj = col_ao(ps_soff(sj) + ib2)
         lamk = col_ao(ps_soff(sk) + ic2); sigl = col_ao(ps_soff(sl) + id2)
         do x = 1, {nv}
            g1(x) = g(x, qab, qcd)
         end do

         ! --- HRR ---
{hrr_body}

{emit_digestion(na, nb, nc, nd)}         end do   ! qab
         end do   ! qcd
         end do   ! cd0
         end do   ! ab0
   end subroutine pci{tag}

   pure real(dp) function cb{tag}(iab_, icd_, nca, nccl, si, sj, sk, sl, &
                                  nbas, ncoltot, ps_soff, col_sh, nshc, nqc, q_col, dsh_c) result(qb)
      !! Schwarz times the density blocks for ONE column combination, at
      !! CONTRACTED resolution: the same test the quartet gets, sharpened from
      !! the merged maximum over the columns to the columns actually being
      !! evaluated. A module procedure rather than a contained one because
      !! `!$acc routine` cannot capture host-subprogram data, and in the same
      !! module as its caller so it inlines.
      !$acc routine seq
      integer,  intent(in) :: iab_, icd_, nca, nccl, si, sj, sk, sl, nbas, ncoltot, nshc, nqc
      integer,  intent(in) :: ps_soff(nbas), col_sh(ncoltot)
      real(dp), intent(in) :: q_col(nqc), dsh_c(nshc, nshc)
      integer :: a_, b_, c_, d_, ja, jb, jc, jd
      a_ = mod(iab_ - 1, nca) + 1; b_ = (iab_ - 1)/nca + 1
      c_ = mod(icd_ - 1, nccl) + 1; d_ = (icd_ - 1)/nccl + 1
      ja = col_sh(ps_soff(si) + a_); jb = col_sh(ps_soff(sj) + b_)
      jc = col_sh(ps_soff(sk) + c_); jd = col_sh(ps_soff(sl) + d_)
      qb = q_col(max(ja, jb)*(max(ja, jb) - 1)/2 + min(ja, jb)) &
           *q_col(max(jc, jd)*(max(jc, jd) - 1)/2 + min(jc, jd))
      qb = qb*max(4.0_dp*dsh_c(ja, jb), 4.0_dp*dsh_c(jc, jd), &
                  dsh_c(ja, jc), dsh_c(ja, jd), dsh_c(jb, jc), dsh_c(jb, jd))
   end function cb{tag}
"""
    if unroll and ncab > 1:
        txt = _unroll_combos(txt, ncab)
    return txt


#: WARP-ROLE DECOMPOSITION, per class: which shell's Cartesian components the
#: work of one quartet is split over.  Measured, not guessed -- see
#: `--role-report` for the numbers behind each entry.
#:
#: THE MECHANISM.  A (dp|dp) thread carries v(400,0:1), g1(400) and vbuf(324):
#: ~1500 doubles against 127 a thread can hold, so ptxas spills 15.9 kB per
#: thread and Nsight Compute has the kernel 94% stalled on local-memory
#: latency with the FP64 pipe 2.6% busy.  Capping registers is monotonically
#: worse and a hand-written CUDA port of the SAME decomposition lands on the
#: same 255 registers: the wall is the working set, not the compiler.
#:
#: What moves it is giving each thread a SLICE of the quartet.  Thread
#: (quartet, role) owns the HRR outputs whose component of the split shell is
#: `role`, runs only the VRR dependency closure those outputs need, and
#: accumulates only its g1 entries.  Redundant arithmetic across roles --
#: 2.9x for (dp|dp) -- is cheaper than the local memory it replaces: 15.9 kB
#: -> 1.7 kB spill frame, local traffic per primitive quartet 318 -> 64
#: sectors, FP64 pipe 2.6% -> 41%, 6.7x on the race harness.  Pure
#: `do concurrent (i, role)`: no thread hierarchy is named, and nvfortran
#: was measured to keep `role` uniform across a warp in either index order
#: (it picks the fast thread index from the access pattern, not the source).
#:
#: THE COST.  Each role loads the density blocks and does its slice of the
#: digestion, and the Fock blocks that do not involve the split shell get one
#: partial atomic update per role instead of one -- (dp|dp): 342 atomics per
#: quartet instead of 117.  Screening is per thread but before the split, so
#: a screened quartet costs one early return per role.
#: Measured on the race harness (scalar vs role, same work, screening off,
#: ns per primitive quartet, speedup at 1/2/3 primitives per shell):
#:   2121 c: 2.60 / 3.83 / 4.48      2111 c: 1.69 / 2.37 / 2.52
#:   2221 c: 1.22 / 2.16 / 2.47      2211 c: 1.12 / 1.13 / 1.13
#:   2120 c: 1.02 / 1.21 / 1.26
#: Left scalar, the split LOSES: 2110 c 0.57-0.90, 2010 c 0.42-0.60,
#: 2011 c 0.35-0.65 (a: 0.18-0.22), 2111 a 0.65-1.52 (c is better).
#: The small classes never spilled, so the redundancy buys nothing.
#: Splitting over the ket's first shell (c) won everywhere it was tried: the
#: HRR closure of one ket component is the smallest slice of the VRR.
#: Extended to every canonical class whose scalar working set spills. Only 21
#: of the 81 generated classes are ever dispatched (bra >= ket, la >= lb,
#: lc >= ld), and across those the split's sign tracks the scalar working set
#: -- 2*v + g1 + vbuf doubles -- almost exactly:
#:     win  2221 2748   2211 1374   2121 1524   2120 708   2111 762
#:     lose 2011  354   2010  138   2110  294
#: so the threshold sits between ~350 and ~700 doubles, which is where the
#: kernel starts spilling to local memory and the redundancy finally buys
#: something. Two unsplit classes above it are added here; every other
#: unsplit class is at or below 381 doubles and is left alone.
#:
#: (dd|dd) IS ABOVE THE THRESHOLD AND IS STILL LEFT SCALAR, on compile time.
#: Its working set is 4971 doubles, the largest of any class, and splitting it
#: six ways is structurally the most promising of the lot -- accumulator 961
#: down to 310. But six roles over a 2179-statement VRR generate 52563 lines,
#: a quarter of everything this script emits, and take the full build from 20
#: minutes to 45. On cc-pVDZ that buys about 0.4 s of a 68 s Fock build,
#: because (dd|dd) is 0.8% of integral time when there is one d shell per
#: heavy atom. Revisit for a basis with more d functions, or for f, where the
#: class is worth more and the trade may reverse.
#: (dd|ps) IS NOT HERE, and must not be put back without checking cc80.
#: Its role kernel makes nvfortran 26.5 die with a signal 11 in fort2 --
#: "TERMINATED by signal 11", no diagnostic -- when compiled for cc80. The
#: SAME file compiles for cc70, which is how it passed review on a V100 and
#: broke an A100 build. Every other role class here compiles for both. If a
#: later compiler fixes it the class is worth about 1.8 s of a 68 s Fock
#: build on the silica slice.
ROLE_SPLITS = {"2121": "c", "2111": "c", "2221": "c", "2211": "c", "2120": "c",
               "2220": "c"}


def _vrr_stmts(vrr):
    """(dst, srcs) per VRR assignment, in program order; comments dropped."""
    raw, cur = [], []
    for line in vrr.split("\n"):
        cur.append(line)
        if not line.rstrip().endswith("&"):
            raw.append("\n".join(cur))
            cur = []
    out = []
    for text in raw:
        flat = text.replace("&\n", " ")
        m = re.match(r"\s*v\((\d+),(\d+)\)\s*=", flat)
        if not m:
            continue
        dst = (int(m.group(1)), int(m.group(2)))
        srcs = [(int(a), int(b)) for a, b in
                re.findall(r"v\((\d+),(\d+)\)", flat.split("=", 1)[1])]
        out.append((text, dst, srcs))
    return out


def _hrr_stmts(hrr):
    """vbuf index -> (statement text, set of g(...) indices it reads)."""
    out = {}
    for m in re.finditer(r"([ \t]*vbuf\((\d+)\)\s*=.*?)(?=\n[ \t]*vbuf\(|\Z)",
                         hrr, re.S):
        out[int(m.group(2))] = (m.group(1).rstrip("\n"),
                                {int(x) for x in re.findall(r"g\((\d+)\)", m.group(1))})
    return out


def role_plan(la, lb, lc, ld, vrr, hrr, shell):
    """Per-role work for splitting over the components of `shell` (a|b|c|d).

    Returns (stmts, roles) where stmts is the parsed VRR and each role is a
    dict: keep (VRR statement ids, in order), vmap (v index -> compact),
    g1 (sorted final-buffer indices, 1-based), gmap, vb (vbuf indices in
    digestion order), and the component list `comps`."""
    stmts = _vrr_stmts(vrr)
    cur = (la + lb + lc + ld) % 2
    defs, deps = {}, []
    for k, (_t, dst, srcs) in enumerate(stmts):
        deps.append({defs[s] for s in srcs})
        defs[dst] = k
    hs = _hrr_stmts(hrr)
    na, nb, nc, nd = [(l + 1)*(l + 2)//2 for l in (la, lb, lc, ld)]
    assert len(hs) == na*nb*nc*nd
    nX = {"a": na, "b": nb, "c": nc, "d": nd}[shell]

    def closure(t):
        seen, stack = set(), list(t)
        while stack:
            k = stack.pop()
            if k in seen:
                continue
            seen.add(k)
            stack.extend(deps[k])
        return seen
    roles = []
    for comp in range(nX):
        vb = [1 + ia + na*ib + na*nb*ic + na*nb*nc*id
              for id in range(nd) for ic in range(nc)
              for ib in range(nb) for ia in range(na)
              if {"a": ia, "b": ib, "c": ic, "d": id}[shell] == comp]
        g = sorted(set().union(*(hs[i][1] for i in vb)))
        keep = sorted(closure([defs[(i, cur)] for i in g]))
        vidx = sorted({i for k in keep for i in
                       [stmts[k][1][0]] + [s[0] for s in stmts[k][2]]})
        roles.append(dict(comps=[comp], keep=keep,
                          vmap={i: k + 1 for k, i in enumerate(vidx)},
                          g1=g, gmap={i: k + 1 for k, i in enumerate(g)},
                          vb=vb, bmap={i: k + 1 for k, i in enumerate(vb)}))
    return stmts, roles, hs


def role_report(la, lb, lc, ld, vrr, hrr):
    """One line per candidate split shell: roles, redundancy, sizes."""
    na, nb, nc, nd = [(l + 1)*(l + 2)//2 for l in (la, lb, lc, ld)]
    tot = len(_vrr_stmts(vrr))
    hs = _hrr_stmts(hrr)
    allg = len(set().union(*(v[1] for v in hs.values())))
    cols = []
    for shell, nX in (("a", na), ("b", nb), ("c", nc), ("d", nd)):
        if nX == 1:
            continue
        _s, roles, _h = role_plan(la, lb, lc, ld, vrr, hrr, shell)
        R = sum(len(r["keep"]) for r in roles)/tot
        cols.append(f"{shell}: {nX} roles R={R:.2f} closure<={max(len(r['keep']) for r in roles)}"
                    f" acc={max(len(r['g1']) for r in roles)}")
    return (f"  ({la}{lb}|{lc}{ld}) VRR {tot} g1 {allg} vbuf {na*nb*nc*nd}: "
            + " | ".join(cols))


def emit_role_kernel(la, lb, lc, ld, vrr, hrr, pieces, shell):
    """The role kernel of one class: driver pcr<tag> + item pcri<tag>.

    Built from the SAME pieces as the scalar kernel (`_scalar_pieces`): the
    head -- segment decode, Schwarz and density screens, pair data -- is
    shared by every role and runs once per thread; a `select case (role)`
    then holds, per role, the primitive loops with that role's VRR closure
    and compacted arrays, its slice of the HRR, and the digestion restricted
    to its components of the split shell."""
    tag = f"{la}{lb}{lc}{ld}"
    na, nb, nc, nd = [(l + 1)*(l + 2)//2 for l in (la, lb, lc, ld)]
    nv = gv.ncum(la + lb)*gv.ncum(lc + ld)
    stmts, roles, hs = role_plan(la, lb, lc, ld, vrr, hrr, shell)
    nrole = len(roles)
    ivar = {"a": "ia", "b": "ib", "c": "ic", "d": "id"}[shell]
    nX = {"a": na, "b": nb, "c": nc, "d": nd}[shell]
    nvr = max(len(r["vmap"]) for r in roles)
    ngr = max(len(r["g1"]) for r in roles)
    nbr = max(len(r["vb"]) for r in roles)

    # --- driver: the scalar driver with a role index on the loop
    drv = pieces["drv"].replace("pcsi" + tag, "pcri" + tag).replace("pcs" + tag, "pcr" + tag)
    drv = drv.replace("SCALAR driver.", "WARP-ROLE driver (split over shell "
                      + shell + f", {nrole} roles).")
    drv = drv.replace("      integer(int64) :: g0, g1, nr, i\n",
                      "      integer(int64) :: g0, g1, nr, i, it\n      integer :: role\n")
    old = ("      do concurrent(i=1:nr)\n"
           "         call pcri" + tag + "(g0 + (i - 1)*nranks, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &\n"
           "                       npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &\n"
           "                       nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &\n"
           "                       sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat)\n"
           "      end do\n")
    assert drv.count(old) == 1
    drv = drv.replace(old, f"""      !
      ! ONE index, decoded.  Thread i holds quartet `it` and role `role` with
      ! 32 consecutive i sharing a role: warp = 32 quartets x ONE role, so the
      ! role bodies never share a warp.  A two-index `do concurrent (i, role)`
      ! was tried first and nvfortran made `role` the fast thread index --
      ! 5.3 active lanes of 32 and the kernel slower than the scalar one.  The
      ! single index leans only on what every kernel here already leans on:
      ! consecutive iterations are consecutive threads (the ket-uniform
      ! decode in the item routine depends on the same thing).
      !
      ! `local(it, role)` is redundant for the standard -- an unnamed scalar
      ! assigned in a `do concurrent` is already per-iteration -- but the
      ! OpenMP port rewrites this construct as `!$omp parallel do`, whose
      ! default is SHARED.  Naming them here is what makes that rewrite a
      ! `private(it, role)` instead of a data race; `tools/dc_locality_lint.py`
      ! enforces it.
      !
      do concurrent(i=1:((nr + 31)/32)*32*{nrole}) local(it, role)
         role = int(mod((i - 1)/32, {nrole}_int64)) + 1
         it = ((i - 1)/(32*{nrole}))*32 + mod(i - 1, 32_int64) + 1
         if (it <= nr) then
            call pcri{tag}(g0 + (it - 1)*nranks, role, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                          nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat)
         end if
      end do
""")

    # --- item: signature and declarations
    sig = pieces["sig"].replace("pcsi" + tag + "(gt, ", "pcri" + tag + "(gt, role, ")
    assert "pcri" + tag in sig
    sig = sig.replace("      integer(int64), intent(in) :: gt\n",
                      "      integer(int64), intent(in) :: gt\n      integer, intent(in) :: role\n")
    decls = pieces["decls"]
    decls = re.sub(r"      real\(dp\) :: v\(\d+, 0:1\), g1\(\d+\), vbuf\(\d+\)\n",
                   f"      real(dp) :: v({nvr}, 0:1), g1({ngr}), vbuf({nbr})\n", decls)
    decls += "      integer  :: gi\n"
    for r, role in enumerate(roles):
        decls += (f"      integer, parameter :: rc{r}({len(role['comps'])}) = [" +
                  ", ".join(str(c) for c in role["comps"]) + "]\n")

    prims, hd = pieces["prims"], pieces["hrr_and_digest"]
    ihrr_end = hd.index(DIGEST_MARKER) - len("         !\n")
    digest = hd[ihrr_end:]
    ivrr = prims.index(vrr)
    acc_old = (f"                  do x = 1, {nv}\n"
               "                     g1(x) = g1(x) + v(x, cur)\n"
               "                  end do\n")
    assert prims.count(acc_old) == 1
    zero_old = f"            do x = 1, {nv}\n               g1(x) = 0.0_dp\n            end do\n"
    assert prims.count(zero_old) == 1
    cases = []
    for r, role in enumerate(roles):
        vmap, gmap, bmap = role["vmap"], role["gmap"], role["bmap"]

        def cv(t):
            return re.sub(r"\bv\((\d+),(\d+)\)",
                          lambda m: f"v({vmap[int(m.group(1))]},{m.group(2)})", t)
        body = "\n".join(cv(stmts[k][0]) for k in role["keep"])
        body += f"\n               cur = {(la + lb + lc + ld) % 2}"
        acc = "".join(f"                  g1({gmap[i]}) = g1({gmap[i]}) + v({vmap[i]}, cur)\n"
                      for i in role["g1"])
        pr = prims[:ivrr] + body + prims[ivrr + len(vrr):]
        pr = pr.replace(acc_old, acc)
        pr = pr.replace(zero_old, f"            do x = 1, {len(role['g1'])}\n"
                                  "               g1(x) = 0.0_dp\n            end do\n")
        hr = "\n".join(
            re.sub(r"\bvbuf\((\d+)\)", lambda m: f"vbuf({bmap[int(m.group(1))]})",
                   re.sub(r"\bg\((\d+)\)", lambda m: f"g1({gmap[int(m.group(1))]})",
                          hs[i][0]))
            for i in role["vb"])
        dg = digest.replace(f"         do idx = 1, {na*nb*nc*nd}\n",
                            f"         do idx = 1, {len(role['vb'])}\n")
        # The J-only branch indexes vbuf in closed form; a role's vbuf holds
        # only its one component of the split shell, so that shell drops out
        # of the strides.
        assert len(role["comps"]) == 1
        jl_old = f"(wq*vbuf(1 + ia + {na}*ib + {na*nb}*ic + {na*nb*nc}*id))"
        assert dg.count(jl_old) == 2, (tag, shell)
        dims = [("ia", na), ("ib", nb), ("ic", nc), ("id", nd)]
        terms, stride = ["1"], 1
        for v, n in dims:
            if v == ivar:
                continue
            terms.append(v if stride == 1 else f"{stride}*{v}")
            stride *= n
        dg = dg.replace(jl_old, "(wq*vbuf(" + " + ".join(terms) + "))")
        n_loops = dg.count(f"do {ivar} = 0, {nX - 1}")
        # 10, plus one: the jfac = 0 zero-fill of the Coulomb blocks repeats
        # the gather loop of whichever of dab/dcd holds the split shell.
        assert n_loops == 11, (tag, shell, n_loops)
        dg = re.sub(r"^([ \t]*)do " + ivar + f" = 0, {nX - 1}\n",
                    lambda m: f"{m.group(1)}do gi = 1, {len(role['comps'])}\n"
                              f"{m.group(1)}   {ivar} = rc{r}(gi)\n", dg, flags=re.M)
        cases.append(f"         case ({r + 1})\n" + pr + hr + "\n" + dg)
    item = (sig + decls + pieces["head"] + pieces["decode"]
            + "         select case (role)\n" + "".join(cases)
            + "         end select\n"
            + "   end subroutine pcri" + tag + "\n")
    return drv + "\n" + item


def emit_kernel(la, lb, lc, ld, cidx, vrr_body, hrr_body):
    """The scalar kernel for one class.

    One accumulator, the pair coefficients folded into pp_cs, no column
    decode. It is still SLICED OUT of the block template below, because that
    template is where the VRR, HRR and digestion are written and there is no
    reason to keep two copies of them -- but the block kernel itself is no
    longer emitted.

    It used to be. The blocked kernel accumulated its VRR block per column
    combination so the columns of a general contraction could share it, and
    was the answer to cc-pVDZ costing five times 6-31G*. `trc_decontract`
    is a better answer: it splits the contraction before anything is built,
    which leaves every shell single-column and the blocked kernel with
    nothing to share. Emitting it cost 162 subroutines of the 324 here and
    roughly half the compile time of this library for code that production
    could no longer reach."""
    tag = f"{la}{lb}{lc}{ld}"
    txt = _emit_block(la, lb, lc, ld, cidx, vrr_body, hrr_body)
    # The scalar kernel is sliced out of the LOOP form: its markers are
    # statements rather than the unrolled form's literal blocks, and the
    # two differ only after the combination head, which the scalar kernel
    # replaces wholesale anyway.
    txt_l = _emit_block(la, lb, lc, ld, cidx, vrr_body, hrr_body, unroll=False)

    m1 = "         if (ncab*nccd == 1) then\n"
    i1 = txt.index(m1)
    i2 = txt.index("\n         else\n            do x = 1, ", i1) + 1
    if "! === combinations" in txt:
        m3 = "         end if\n\n         ! === combinations"
    else:
        m3 = "         end if\n\n         do qcd = 1, ncd_c\n"
    i3 = txt.index(m3)
    scalar_prims = txt_l[i1 + len(m1):i2]

    # --- scalar item: head, one accumulator, folded coefficients, no decode
    ihead = txt_l.index("      ! locate the segment")
    igc = txt_l.index("         ! GENERAL CONTRACTION.")
    head = txt_l[ihead:igc]
    prims = scalar_prims
    for drop in ("               ki = pp_ki(kp); kj = pp_kj(kp)\n",
                 "               wab = ps_coef(ps_coff(si) + ki)*ps_coef(ps_coff(sj) + kj)\n",
                 "                  kk = pp_ki(kq); kl = pp_kj(kq)\n",
                 "                  w = wab*ps_coef(ps_coff(sk) + kk)*ps_coef(ps_coff(sl) + kl)\n"):
        assert prims.count(drop) == 1, drop
        prims = prims.replace(drop, "")
    prims = prims.replace("*pp_c(kp)*pp_c(kq)", "*pp_cs(kp)*pp_cs(kq)")
    # pp_cs carries the coefficients, so the prescreen weight is already in pref
    assert prims.count("if (abs(pref*w) <= pcut) cycle") == 1
    # The pairs are sorted by |c|/p descending (build_pairs_hgp), so the
    # bound 2 pi^2.5 (|c_p|/zeta)(|c_q|/eta)/sqrt(zeta) falls along the ket
    # loop: EXIT at the first failure, and skip a bra primitive whose best
    # ket already fails. Exact pref is still tested for the rare quartet
    # the loose bound admits.
    prims = prims.replace("if (abs(pref*w) <= pcut) cycle",
                          "if (bnd*abs(pp_cs(kq))/eta <= pcut) exit\n"
                          "                  if (abs(pref) <= pcut) cycle")
    prims = prims.replace("               zeta = pp_p(kp)\n",
                          "               zeta = pp_p(kp)\n"
                          "               bnd = TWO_PI_2_5*abs(pp_cs(kp))/(zeta*sqrt(zeta))\n"
                          "               if (bnd*abs(pp_cs(offcd + 1))/pp_p(offcd + 1) <= pcut) cycle\n", 1)
    prims = prims.replace("g1(x) = g1(x) + w*v(x, cur)", "g1(x) = g1(x) + v(x, cur)")
    igl = prims.rfind("            do x = 1, ")   # the gl copy, last loop
    assert igl > 0 and "g(x, 1, 1) = g1(x)" in prims[igl:]
    prims = prims[:igl]
    ihrr = txt_l.index("         ! --- HRR ---\n")
    idig = txt_l.index(DIGEST_MARKER) - len("         !\n")
    iend = txt_l.index("         end do   ! qab\n")
    hrr_and_digest = txt_l[ihrr:iend]
    decode = emit_scalar_decode()
    idecl0 = txt_l.index("      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl\n")
    decls = txt_l[idecl0:ihead]
    decls = re.sub(r"      real\(dp\) :: v\((\d+), 0:1\), g\(.*?\), g1\((\d+)\), vbuf\((\d+)\)\n",
                   r"      real(dp) :: v(\1, 0:1), g1(\2), vbuf(\3)\n", decls)
    decls = re.sub(r"      real\(dp\) :: cab\(\d+\), ccd\(\d+\).*\n", "", decls)
    decls = decls.replace("      logical :: dij, dkl, dpq, ok\n", "      logical :: dij, dkl, dpq\n")
    decls = re.sub(r"      integer  :: offa\(.*\n", "", decls)
    decls = decls.replace("      real(dp) :: wq, w, wab\n", "      real(dp) :: wq\n")
    decls = decls.replace("      integer  :: nca, ncb, nccl, ncdl, npi, npj, npk, npl, ncab, nccd, ab0, cd0, nab_c, ncd_c, qab, qcd\n", "")
    decls = decls.replace("      integer  :: ki, kj, kk, kl, kpl, kql, ia2, ib2, ic2, id2, iabc, icdc\n", "")
    head = head.replace("         nca = ps_ncol(si); ncb = ps_ncol(sj); nccl = ps_ncol(sk); ncdl = ps_ncol(sl)\n", "")
    head = head.replace("         npi = ps_np(si); npj = ps_np(sj); npk = ps_np(sk); npl = ps_np(sl)\n", "")

    #
    # BIN-INTERLEAVED PRIMITIVE-PAIR LAYOUT.
    #
    # Replace the pp_off/pp_n shell-pair lookup and the pp_p/pp_r/pp_ra/pp_cs
    # reads of the primitive loop with the bin-interleaved qp_*/sq_* tables
    # (trc_bins%build_pair_quartets): `slotab`/`slotcd` are the pair's own
    # slot, aligned with sp_i/sp_j/sp_q, and `kp`/`kq` become addresses into
    # qp_* rather than flat pp_* indices. Same screens, same cycle/exit, same
    # order of operations -- only where each value comes from changes.
    #
    decls = decls.replace("      integer :: keyab, keycd, offab, offcd, nab, ncd\n",
                          "      integer :: slotab, slotcd, nab, ncd\n"
                          "      integer :: kp, kq, k1, k2\n")
    assert decls.count("      integer :: kp, kq, d, x, cur, ia, ib, ic, id, idx, idens\n") == 1
    decls = decls.replace("      integer :: kp, kq, d, x, cur, ia, ib, ic, id, idx, idens\n",
                          "      integer :: d, x, cur, ia, ib, ic, id, idx, idens\n")

    head_old = ("         keyab = (si - 1)*nbas + sj\n"
                "         keycd = (sk - 1)*nbas + sl\n"
                "         offab = pp_off(keyab); nab = pp_n(keyab)\n"
                "         offcd = pp_off(keycd); ncd = pp_n(keycd)\n"
                "\n"
                "         abx = pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1)\n"
                "         aby = pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2)\n"
                "         abz = pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3)\n"
                "         cdx = pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1)\n"
                "         cdy = pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2)\n"
                "         cdz = pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3)\n")
    assert head.count(head_old) == 1
    head = head.replace(head_old,
                        "         slotab = sOA(seg) + iab; slotcd = sOB(seg) + icd\n"
                        "         nab = sq_n(slotab); ncd = sq_n(slotcd)\n"
                        "\n"
                        "         abx = sq_abx(slotab); aby = sq_aby(slotab); abz = sq_abz(slotab)\n"
                        "         cdx = sq_abx(slotcd); cdy = sq_aby(slotcd); cdz = sq_abz(slotcd)\n")

    prims_old = ("            do kp = offab + 1, offab + nab\n"
                 "               zeta = pp_p(kp)\n"
                 "               bnd = TWO_PI_2_5*abs(pp_cs(kp))/(zeta*sqrt(zeta))\n"
                 "               if (bnd*abs(pp_cs(offcd + 1))/pp_p(offcd + 1) <= pcut) cycle\n")
    assert prims.count(prims_old) == 1
    prims = prims.replace(prims_old,
                          "            do k1 = 1, nab\n"
                          "               kp = sPA(seg) + (k1 - 1)*sA(seg) + iab\n"
                          "               zeta = qp_z(kp)\n"
                          "               bnd = TWO_PI_2_5*abs(qp_c(kp))/(zeta*sqrt(zeta))\n"
                          "               if (bnd*abs(qp_c(sPB(seg) + icd))/qp_z(sPB(seg) + icd) <= pcut) cycle\n")
    prims = prims.replace("               rpx = pp_r(kp, 1); rpy = pp_r(kp, 2); rpz = pp_r(kp, 3)\n"
                          "               pax = rpx - pp_ra(kp, 1)\n"
                          "               pay = rpy - pp_ra(kp, 2)\n"
                          "               paz = rpz - pp_ra(kp, 3)\n",
                          "               rpx = qp_px(kp); rpy = qp_py(kp); rpz = qp_pz(kp)\n"
                          "               pax = rpx - sq_ax(slotab)\n"
                          "               pay = rpy - sq_ay(slotab)\n"
                          "               paz = rpz - sq_az(slotab)\n")
    prims = prims.replace("               do kq = offcd + 1, offcd + ncd\n"
                          "                  eta = pp_p(kq)\n",
                          "               do k2 = 1, ncd\n"
                          "                  kq = sPB(seg) + (k2 - 1)*sNB(seg) + icd\n"
                          "                  eta = qp_z(kq)\n")
    prims = prims.replace("                  pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_cs(kp)*pp_cs(kq)\n"
                          "                  if (bnd*abs(pp_cs(kq))/eta <= pcut) exit\n",
                          "                  pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*qp_c(kp)*qp_c(kq)\n"
                          "                  if (bnd*abs(qp_c(kq))/eta <= pcut) exit\n")
    prims = prims.replace("                  pqx = rpx - pp_r(kq, 1)\n"
                          "                  pqy = rpy - pp_r(kq, 2)\n"
                          "                  pqz = rpz - pp_r(kq, 3)\n"
                          "                  qcx = pp_r(kq, 1) - pp_ra(kq, 1)\n"
                          "                  qcy = pp_r(kq, 2) - pp_ra(kq, 2)\n"
                          "                  qcz = pp_r(kq, 3) - pp_ra(kq, 3)\n"
                          "                  wc = (zeta*rpx + eta*pp_r(kq, 1))*rzpe\n"
                          "                  wpx = wc - pp_r(kp, 1); wqx = wc - pp_r(kq, 1)\n"
                          "                  wc = (zeta*rpy + eta*pp_r(kq, 2))*rzpe\n"
                          "                  wpy = wc - pp_r(kp, 2); wqy = wc - pp_r(kq, 2)\n"
                          "                  wc = (zeta*rpz + eta*pp_r(kq, 3))*rzpe\n"
                          "                  wpz = wc - pp_r(kp, 3); wqz = wc - pp_r(kq, 3)\n",
                          "                  pqx = rpx - qp_px(kq)\n"
                          "                  pqy = rpy - qp_py(kq)\n"
                          "                  pqz = rpz - qp_pz(kq)\n"
                          "                  qcx = qp_px(kq) - sq_ax(slotcd)\n"
                          "                  qcy = qp_py(kq) - sq_ay(slotcd)\n"
                          "                  qcz = qp_pz(kq) - sq_az(slotcd)\n"
                          "                  wc = (zeta*rpx + eta*qp_px(kq))*rzpe\n"
                          "                  wpx = wc - qp_px(kp); wqx = wc - qp_px(kq)\n"
                          "                  wc = (zeta*rpy + eta*qp_py(kq))*rzpe\n"
                          "                  wpy = wc - qp_py(kp); wqy = wc - qp_py(kq)\n"
                          "                  wc = (zeta*rpz + eta*qp_pz(kq))*rzpe\n"
                          "                  wpz = wc - qp_pz(kp); wqz = wc - qp_pz(kq)\n")

    isig = txt_l.index("   pure subroutine pci" + tag + "(")
    sig = txt_l[isig:idecl0].replace("pci" + tag, "pcsi" + tag)
    scalar_item = (sig + decls + head + prims + decode + hrr_and_digest
                   + "   end subroutine pcsi" + tag + "\n")

    # --- scalar driver: the block driver with the names swapped
    idrv = txt_l.index("   subroutine pc" + tag + "(")
    idrv_end = txt_l.index("   end subroutine pc" + tag + "\n") + len("   end subroutine pc" + tag + "\n")
    drv = txt_l[idrv:idrv_end].replace("pci" + tag, "pcsi" + tag).replace("pc" + tag, "pcs" + tag)
    drv = drv.replace("driver.  The `do concurrent`", "SCALAR driver.  The `do concurrent`")

    pieces = dict(drv=drv, sig=sig, decls=decls, head=head, prims=prims,
                  decode=decode, hrr_and_digest=hrr_and_digest)
    return drv + "\n" + scalar_item, pieces


def _radix(txt):
    """Fill the radix placeholder. A placeholder rather than an f-string
    field because this block is emitted from two templates, one f-string
    and one plain, and a brace works in only one of them."""
    return txt.replace("__RADIX__", str(CLASS_RADIX))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=1)
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="Cartesian order of the AO slots the digestion "
                         "scatters into (mu = mui + ia, etc). This is "
                         "OUTPUT layout only: the VRR sieve's targets "
                         "(gv.cart_list) and the HRR's internal index into "
                         "`g` (gh.cart_cum) stay lexicographic regardless -- "
                         "see the comment where gh.emit_class is called "
                         "below.")
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
    ap.add_argument("--role-report", action="store_true",
                    help="print the warp-role redundancy table for every "
                         "class with a d shell and exit; nothing is written")
    ap.add_argument("--classes", default=None,
                    help="comma-separated four-digit class tags (e.g. "
                         "3330,3333) to emit instead of the full "
                         "(--lmax + 1)^4 sweep; requires --split, since a "
                         "partial set has no single-file form and no "
                         "dispatcher (see the --split handling below). "
                         "Paired with --skeleton to produce one or a few "
                         "classes' splicer template without the full-sweep "
                         "cost of --lmax 3 or 4.")
    ap.add_argument("--skeleton", action="store_true",
                    help="skip the VRR/HRR generation (gen_vrr.py/gen_hrr.py) "
                         "for every selected class and emit a placeholder "
                         "integral body instead. gen_rys.py and "
                         "gen_rys_grad.py only ever read the scalar kernel "
                         "up to the primitive loop's `tval` line, and from "
                         "the digestion onward (see either generator's "
                         "docstring for THE CUT POINTS); the HGP VRR/HRR "
                         "algebra in between is spliced away unread. A "
                         "skeleton module keeps every marker and array size "
                         "those splicers key on -- the same emit_kernel, "
                         "the same digestion, the same declarations -- and "
                         "skips only the part neither of them looks at. At "
                         "(ff|ff) the real VRR/HRR text alone is hundreds "
                         "of thousands of lines; the skeleton's is a "
                         "handful.")
    args = ap.parse_args()
    if args.classes and not args.role_report:
        assert args.split, "--classes writes one module per class; pass --split DIR too"
    L = args.lmax
    order = args.order
    order_val = 0 if order == "lexicographic" else 1
    # Suppressed when default: the shipped set's header must stay
    # byte-identical to what it was before --order existed.
    order_arg = "" if order == "lexicographic" else f" --order {order}"

    lst, idx_v, info = gv.build_tables(2*L)
    idx_h = gh.cart_cum(4*L)

    pieces = []
    out = [f"""!
! One kernel per angular-momentum class: unrolled VRR, unrolled HRR, six-atomic
! folded Fock, and no `select case` in sight.
!
! GENERATED by scripts/gen_perclass.py --lmax {L}{order_arg} -- do not edit.
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
! TRC_CART_ORDER GUARD. Every kernel's digestion below scatters into an AO
! slot fixed at generation time by --order (0 = lexicographic, the default;
! 1 = gamess); TRC_CART_ORDER selects it at compile time, and an undefined
! TRC_CART_ORDER reads as 0 in a preprocessor #if, so the default build needs
! no new -D flag. Same idea as the TRC_LMAX guard docs/GENERATORS.md
! describes, for a different compile-time constant.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
module trc_pc_kernels
   use pic_types, only: dp, int64
   use trc_boys, only: boys_eval, BOYS_MMAX, boys_table, &
                        BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, &
                        BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   implicit none
   private
   public :: pc_dispatch, CLASS_RADIX

   integer, parameter :: CLASS_RADIX = __RADIX__
   !! Radix of the dispatch key, FIXED so it does not encode LMAX.
   !! The caller in trc_binkernel imports this rather than repeating the
   !! literal, because a key formed two ways that disagree dispatches the
   !! wrong kernel and still builds.

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains
"""]

    if args.classes:
        combos = []
        for t in (s.strip() for s in args.classes.split(",")):
            if not t:
                continue
            assert len(t) == 4 and t.isdigit(), f"bad --classes tag {t!r}"
            la, lb, lc, ld = (int(c) for c in t)
            assert max(la, lb, lc, ld) <= L, (
                f"--classes {t} needs --lmax >= {max(la, lb, lc, ld)}")
            combos.append((la, lb, lc, ld))
    else:
        combos = [(la, lb, lc, ld)
                  for la in range(L + 1) for lb in range(L + 1)
                  for lc in range(L + 1) for ld in range(L + 1)]

    names = []
    for la, lb, lc, ld in combos:
        tag = f"{la}{lb}{lc}{ld}"
        if args.skeleton:
            # See --skeleton's help: gen_rys.py/gen_rys_grad.py never read
            # the HGP VRR/HRR algebra at all, only the markers and sizes
            # around it, so the generation gv.emit_class_sieved/
            # gh.emit_class would otherwise do here (THE expensive step --
            # hundreds of thousands of lines at (ff|ff)) is skipped
            # outright. `cur` is the one piece of this range a downstream
            # statement (the accumulate loop right after it) still reads;
            # everything else here is dropped by the splice unread.
            vrr = (f"      ! HGP VRR omitted -- --skeleton (scripts/gen_perclass.py)\n"
                   f"               cur = {(la + lb + lc + ld) % 2}")
            hrr = "         ! HGP HRR omitted -- --skeleton (scripts/gen_perclass.py)\n"
            print(f"  ({la}{lb}|{lc}{ld}): skeleton, HGP body omitted")
        else:
            # Gill's sieve needs to know what the HRR actually reads, so
            # the HRR is generated first and its target set drives the VRR.
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
            # OUTPUT: this is the digestion's own AO slot order, so (unlike
            # hrr_pre above, whose reads set is order-invariant -- see the
            # module docstring) this one gets the requested order.
            hrr = gh.emit_class(la, lb, lc, ld, idx_h, order)
        key = ((la*CLASS_RADIX + lb)*CLASS_RADIX + lc)*CLASS_RADIX + ld
        names.append((key, tag))
        if args.role_report:
            if not args.skeleton and 2 in (la, lb, lc, ld):
                print(role_report(la, lb, lc, ld, vrr, hrr))
            continue
        body, parts = emit_kernel(la, lb, lc, ld, idx_h, vrr, hrr)
        # The role-split kernels are an HGP-only redundant-warp optimisation
        # (pcr{tag}/pcri{tag}) that gen_rys.py/gen_rys_grad.py never read --
        # see gen_rys.py's THE CUT POINTS on why it scopes its search to
        # pcsi{tag} alone. A skeleton module has no real HRR for
        # emit_role_kernel to split, so it is skipped here rather than fed
        # the placeholder.
        if tag in ROLE_SPLITS and not args.skeleton:
            body += "\n" + emit_role_kernel(la, lb, lc, ld, vrr, hrr,
                                            parts, ROLE_SPLITS[tag])
        pieces.append((tag, body))
        out.append(body)

    if args.role_report:
        return
    out.append(f"""
   subroutine pc_dispatch(key, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                          nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                          ndens, dmat, jmat, rank, nranks)
      integer, intent(in) :: key
{PROLOGUE}
      integer, intent(in) :: rank, nranks
      select case (key)""")
    for key, tag in names:
        kern = "pcr" if tag in ROLE_SPLITS else "pcs"
        out.append(f"""      case ({key}); call {kern}{tag}(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                          npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, omega, dsh, nbas, nao, sh_l, ao_off, &
                          nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, ndens, dmat, jmat, rank, nranks)""")
    out.append("""      end select
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
! GENERATED by scripts/gen_perclass.py --lmax {L}{order_arg} --split -- do not edit.
!
! One module per class so no single compilation unit is enormous. The item
! routine stays in the same file as its driver, which is what the inlining
! rule requires; the dispatcher is host code and may cross modules.
!
! TRC_CART_ORDER GUARD -- see this generator's non-split header (above, in
! gen_perclass.py) for the full explanation; same guard, one per class file.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
module {mod}
   use pic_types, only: dp, int64
   use trc_boys, only: boys_eval, BOYS_MMAX, boys_table, &
                         BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, &
                         BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   implicit none
   private
   public :: pcs{tag}{", pcr" + tag if tag in ROLE_SPLITS else ""}

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp

contains
{body}
end module {mod}
""")
            use_lines.append(f"   use {mod}, only: pcs{tag}"
                             + (f", pcr{tag}" if tag in ROLE_SPLITS else ""))
        if args.classes:
            # A --classes run is a requested subset, not the full
            # (--lmax + 1)^4 sweep pc_dispatch's `select case` is built to
            # cover -- writing it here would dispatch every OTHER class to
            # nothing. gen_rys.py/gen_rys_grad.py (the only consumers this
            # mode exists for) read trc_pc_k{tag}.F90 directly and never
            # look at the dispatcher at all, so it is skipped outright
            # rather than written partial and misleading.
            print(f"wrote {len(names)} class module(s) "
                  f"(skeleton={args.skeleton}) into {args.split}; "
                  "no dispatcher -- --classes is a partial set")
            return
        disp = [f"""!
! Kernel dispatcher for the per-class modules.
!
! GENERATED by scripts/gen_perclass.py --lmax {L}{order_arg} --split -- do not edit.
!
! No TRC_CART_ORDER guard here: this file only maps a class key to a call,
! and that mapping never depends on Cartesian order -- the guard lives in
! each class module instead, where the order-dependent digestion actually is.
!
module trc_pc_kernels
   use pic_types, only: dp, int64
   use trc_tables, only: LMAX
""" + "\n".join(use_lines) + """
   implicit none
   private
   public :: pc_dispatch, CLASS_RADIX

   integer, parameter :: CLASS_RADIX = __RADIX__
   !! Radix of the dispatch key, FIXED so it does not encode LMAX.
   !! The caller in trc_binkernel imports this rather than repeating the
   !! literal, because a key formed two ways that disagree dispatches the
   !! wrong kernel and still builds.

contains
"""]
        disp.append(txt[txt.index("   subroutine pc_dispatch"):
                        txt.index("end module trc_pc_kernels")])
        disp.append("end module trc_pc_kernels\n")
        with open(os.path.join(args.split, "trc_pc_kernels.F90"), "w") as fh:
            fh.write(_radix("\n".join(disp)))
        print(f"wrote {len(names)} class modules + dispatcher into {args.split}")
        return

    open(args.output, "w").write(_radix(txt))
    print(f"wrote {args.output}: {len(names)} kernels, {txt.count(chr(10))} lines")


if __name__ == "__main__":
    main()
