#!/usr/bin/env python3
"""Emit a Rys kernel for one class, by splicing terco's own scalar kernel.

WHY SPLICE RATHER THAN GENERATE WHOLE
--------------------------------------
The comparison this exists to serve is HGP against Rys with everything else
held fixed. Reusing the scalar HGP kernel's prologue (segment decode,
screening, the bin-interleaved primitive-pair loads, the primitive-quartet
prescreen) and its folded six-atomic digestion VERBATIM is what makes that
true: only the middle -- the integral itself -- differs, so a timing
difference is the scheme and nothing else.

This reads `src/generated/trc_pc_k{tag}.F90`, the file `gen_perclass.py
--split` already produced for the scalar (decontracted, single-column)
kernel `pcs{tag}`/`pcsi{tag}`. General contraction is not this file's
concern at all: `trc_decontract` splits it away before any kernel runs, so
every class file -- HGP or Rys -- only ever needs the scalar form.

THE CUT POINTS
---------------
`pcsi{tag}` falls into three parts, found by text search rather than by
counting lines, so a reflow of the surrounding generator does not silently
misalign the splice:

  1. HEAD, up to and including `tval = rho*(pqx*pqx + ...)`. Segment decode,
     screens, the pair-table offsets, `abx..cdz`, and the bra/ket primitive
     loop through everything two shell-pair centres and one primitive
     quartet need. Identical for HGP and Rys, kept verbatim.
  2. The INTEGRAL itself -- HGP's Boys evaluation, VRR and the accumulation
     into its pre-HRR block. Replaced by the Rys root lookup, the 2D
     recurrence and the horizontal transfer, run per primitive quartet
     (HGP contracts over primitives and only then transforms once; Rys
     builds the shifted class per primitive quartet instead). This writes
     the final Cartesian block straight into `vbuf`, so there is no
     separate HRR pass afterward.
  3. TAIL: the degeneracy/AO-offset decode (`! Single column on every
     side:`), the digestion (from `! BLOCK-ACCUMULATED DIGESTION.`), and
     `end subroutine`. Scheme-independent, kept verbatim. HGP's own
     "--- HRR ---" text, which sits between the decode and the digestion in
     the source file, is dropped: Rys's vbuf is already the transferred
     block by the time the decode runs.

THE ALGEBRA
-----------
From gamess-libERI's int*_rysgen, reduced: writing P-A and Q-C through PQ
makes the c00/cp00 coefficients collapse to one term each, so the kernel
needs only pax/qcx/pqx -- all of which the prologue already has, since HGP
needs them too.  `tval` IS the Rys argument and `pref` already carries
2*pi^(5/2)/(zeta*eta*sqrt(zpe))*(coefficients). This algebra was validated
against libcint through (dd|dd) on explore/rys-4c (not carried over here);
what checks it in this tree is `test/check_rys_vs_hgp.F90`, comparing the
Fock matrix this produces against HGP's on the same density.

RANGE SEPARATION
----------------
`omega` and `theta` arrive already computed in the spliced HEAD (see
scripts/gen_perclass.py's primitive loop): `tval` is already theta*T, and
`pref` is scaled by sqrt(theta) once the roots/weights are in hand. The one
thing this generator still has to get right itself is `u2` -- see the
comment where it is built (`unrolled_energy_middle` and this module's
`main`): t^2 -> theta*t^2 inside u2's own formula, NOT rho -> theta*rho the
way HGP's VRR does, confirmed by comparing against pyscf's
`mol.with_range_coulomb` (test/check_erf_fock_vs_pyscf.F90).

THE CHUNKED FOCK PASS
---------------------
A class in FOCK_CHUNK_PLAN (the f classes the pair bins launch, and the d
classes FOCK_RYS_DEFAULT routes from HGP to Rys) keeps the
HEAD but replaces the integral AND the digestion: it is cut into chunks of
output components, each its own primitive/root loop and digestion, so one
thread never holds the whole Cartesian block. See the section comment
above FOCK_CHUNK_PLAN for the design and the measurements.
"""
import argparse
import importlib.util
import itertools
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


co = _load("co", os.path.join(HERE, "cartorder.py"))
ri = _load("ri", os.path.join(HERE, "rys_inline.py"))
# ---------------------------------------------------------------------------
# CHUNKED FOCK PASS
# ---------------------------------------------------------------------------
#
# WHY. The f classes are Rys only, and the array form below holds one
# quartet's WHOLE Cartesian block in one thread: vbuf(na*nb*nc*nd), the
# gx/gy/gz and hx/hy/hz tables indexed by runtime loop variables, and the
# twelve density/Fock blocks of the digestion. That is far past the
# 255-register file, so it lives in local memory, and with ~20k resident
# threads it streams through DRAM. Nsight Compute, w20 cc-pVTZ, 1e-13, V100,
# one SCF launch of (fp|dp): 137.5 ms, 255 registers, 7.5 KB of stack, 103 GB
# of local loads and 38 GB of local stores against 3.2 GB of global loads,
# 93 GB of DRAM traffic, 3.5% SM throughput. The heaviest classes were worse
# for a second reason: (ff|ff) has 10000 components and few quartets, so its
# 97 KB-stack threads were also too FEW to fill the GPU.
#
# WHAT. The item routine is ordered, as the gradient's CHUNKED RHF COMPACT
# PASS (gen_rys_grad.py) is,
#
#     chunk of output components -> primitive quartet -> root
#
# The four shells split into INSIDE and OUTER (FOCK_CHUNK_PLAN). A chunk is
# every Cartesian component of the inside shells for one component of each
# outer shell, and its integrals accumulate in named scalars vq_0.. across
# the primitive and root loops. Per root, the chunk builds the 2D-table
# corner its components read and the horizontal transfer, all named scalars
# with literal indices: over an inside shell with the literal binomials, over
# an outer shell with RUNTIME weights -- wb_s = C(eb, s) ab^(eb-s), zero past
# eb, set once per chunk from the outer component's exponent eb (the cn
# tables) -- so that h(..) = sum_s wb_s g(.. + s ..) has literal indices
# either way and EVERY chunk runs the same code. The source is one chunk
# long, not one class long. Then the chunk is digested on its own: each Fock
# element it touches gets one partial sum (the chunk's components that land
# on it times the density element each pairs with, read straight from dmat)
# and one atomic.
#
# The outer sets supported are, per side, nothing, the second shell, or both
# (_FOCK_OUTER_SETS): an outer FIRST shell alone would pick a 2D-table row
# by a runtime index, and a runtime index into per-thread scratch is exactly
# what puts it back in local memory.
#
# Chunk -> primitive, not the other way round: with the chunk loop inside
# the primitive loop the compiler hoists every chunk's state out of it and
# rebuilds the whole block in registers (measured on the gradient, see
# gen_rys_grad.py).
#
# LOOP OR THREAD. 'loop' runs every chunk in one thread per quartet;
# 'thread' gives each (quartet, chunk) its own `do concurrent` iteration,
# quartet fastest so a warp shares one chunk, and redoes the decode and the
# screens per chunk. 'thread' multiplies the parallelism by the chunk count,
# which is what the heavy classes lacked: (ff|ff) 502 -> 225 ms over the run
# in loop form, 15 ms in thread form.
#
# The price is recomputation -- the primitive setup, the roots and the
# 2D-table corner per chunk -- and atomics: a Fock element several chunks
# touch is written once per chunk instead of once per quartet. Which shape
# wins is not predictable from a register model (the gradient's live-set
# rule picks badly here), so FOCK_CHUNK_PLAN is measured.
#
# The density loop runs per chunk, around the digestion only, so any ndens
# (UHF, the response batches) reuses the chunk's integrals the way the array
# form reuses vbuf. kfac = 0 skips the exchange half per chunk: one uniform
# test per chunk, not per component.
#
# A first version kept every chunk's indices literal, gen_rys_grad.py's
# select case over chunks. It was no faster where it compiled (within 5% of
# the plan below on every class it was built for) and it did not compile in
# budget: 570-583 s per object at 900-1080 components, and (ff|fd) at 1000
# had not finished after 25 minutes. The runtime-outer form takes at most
# 12 s per object.

#: The measured chunk shape per class: outer shells (0..3 = a, b, c, d; ()
#: is one chunk of the whole class) and the form. w20 (the first 60 atoms of
#: w172), cc-pVTZ, 1e-13, V100, every one of the nine shapes in both forms;
#: the per-class time over the whole 11-iteration run, array form -> this
#: plan, in ms (the full table is in the commit that introduced this):
#:   (fp|dp) 1496 -> 385, (fd|dp) 1120 -> 327, (fp|ps) 1104 -> 231,
#:   (fp|pp) 1056 -> 218, (fp|ds) 966 -> 254, (fs|dp) 943 -> 266,
#:   (fd|ps) 835 -> 226, (fp|dd) 760 -> 208, (fd|pp) 757 -> 170,
#:   (fd|ds) 650 -> 194, (fd|fp) 588 -> 155, (fs|pp) 580 -> 117,
#:   (fp|fp) 578 -> 101, (ff|dp) 518 -> 141, (ff|ff) 502 -> 15,
#:   (ff|pp) 481 -> 77, (ff|fp) 476 -> 63, (fd|dd) 476 -> 159,
#:   (fs|dd) 409 -> 129, (fp|fs) 386 -> 108, (fd|fs) 329 -> 90,
#:   (ff|fd) 326 -> 46, (ff|ps) 287 -> 95, (ff|dd) 281 -> 65,
#:   (fd|fd) 276 -> 68, (ff|ds) 271 -> 90, (fs|fs) 266 -> 48,
#:   (fs|ds) 177 -> 146, (ff|fs) 141 -> 47, (ff|ss) 122 -> 39.
#:
#: The d classes listed are the ones FOCK_RYS_DEFAULT routes away from HGP,
#: raced the same way against HGP (pcs/pcr) on six systems; see there.
#: A class not listed keeps the array form: the d and p classes that stay
#: on HGP by default (their Rys kernels exist for TRC_RYS_CLASSES and are
#: exactly as they were), the tags the pair bins never launch (la < lb,
#: lc < ld or (la, lb) < (lc, ld); trc_bins.F90's ANGULAR-MOMENTUM
#: CANONICAL), and (fs|ss), (fs|ps), (fp|ss), (fd|ss), which measured the
#: same in every shape (104, 153, 81, 64 ms).
FOCK_CHUNK_PLAN = {
    "2111": ((1,), "thread"),
    "2120": ((1,), "thread"),
    "2211": ((1,), "thread"),
    "2121": ((1, 3), "thread"),
    "2210": ((2, 3), "thread"),
    "2221": ((2, 3), "thread"),
    "2222": ((2, 3), "thread"),
    "3020": ((), "loop"),
    "3030": ((), "loop"),
    "3300": ((), "loop"),
    "3011": ((3,), "thread"),
    "3021": ((3,), "thread"),
    "3022": ((3,), "thread"),
    "3310": ((1,), "loop"),
    "3320": ((1,), "loop"),
    "3110": ((2, 3), "loop"),
    "3120": ((2, 3), "loop"),
    "3210": ((2, 3), "loop"),
    "3111": ((2, 3), "thread"),
    "3121": ((2, 3), "thread"),
    "3122": ((2, 3), "thread"),
    "3130": ((2, 3), "thread"),
    "3131": ((2, 3), "thread"),
    "3211": ((2, 3), "thread"),
    "3220": ((2, 3), "thread"),
    "3221": ((2, 3), "thread"),
    "3222": ((2, 3), "thread"),
    "3230": ((2, 3), "thread"),
    "3231": ((2, 3), "thread"),
    "3232": ((2, 3), "thread"),
    "3311": ((1, 3), "thread"),
    "3321": ((1, 3), "thread"),
    "3330": ((1, 3), "thread"),
    "3322": ((0, 1), "thread"),
    "3331": ((0, 1), "thread"),
    "3332": ((0, 1), "thread"),
    "3333": ((0, 1), "thread"),
    # g classes up to la+lb+lc+ld = 9 (TERCO_RYS_G_LTOT). w20 cc-pVQZ, 1e-12,
    # V100, ms per Fock build (4th iteration), all nine shapes raced: best,
    # runner-up, and the MD exchange kernel it replaces. Sum 2.62 s against
    # 6.90 s in MD and 8.92 s in the array form (GPU4PySCF: 2.34 s).
    "4000": ((3,), "thread"),  # (gs|ss) 18; next (1, 3) thread 18; MD 130
    "4010": ((1,), "thread"),  # (gs|ps) 33; next (1, 3) thread 33; MD 282
    "4011": ((3,), "thread"),  # (gs|pp) 35; next (1, 3) thread 35; MD 173
    "4020": ((1,), "thread"),  # (gs|ds) 56; next (1, 3) thread 56; MD 258
    "4021": ((1, 3), "thread"),  # (gs|dp) 103; next (3,) thread 103; MD 322
    "4022": ((1, 3), "thread"),  # (gs|dd) 71; next (3,) thread 71; MD 225
    "4030": ((2, 3), "loop"),  # (gs|fs) 58; next (1, 2, 3) thread 67; MD 239
    "4031": ((1, 2, 3), "thread"),  # (gs|fp) 125; next (2, 3) thread 125; MD 327
    "4032": ((1, 2, 3), "thread"),  # (gs|fd) 137; next (2, 3) thread 137; MD 318
    "4040": ((2, 3), "loop"),  # (gs|gs) 13; next (1, 2, 3) thread 13; MD 80
    "4100": ((3,), "thread"),  # (gp|ss) 21; next () loop 21; MD 181
    "4110": ((2, 3), "loop"),  # (gp|ps) 76; next (1,) thread 77; MD 320
    "4111": ((2, 3), "loop"),  # (gp|pp) 75; next (2, 3) thread 76; MD 194
    "4120": ((2, 3), "loop"),  # (gp|ds) 91; next (2, 3) thread 94; MD 298
    "4121": ((2, 3), "thread"),  # (gp|dp) 184; next (2, 3) loop 189; MD 330
    "4122": ((2, 3), "thread"),  # (gp|dd) 114; next (2, 3) loop 115; MD 150
    "4130": ((2, 3), "thread"),  # (gp|fs) 93; next (2, 3) loop 94; MD 255
    "4131": ((2, 3), "thread"),  # (gp|fp) 172; next (2, 3) loop 184; MD 260
    "4140": ((2, 3), "loop"),  # (gp|gs) 35; next (2, 3) thread 35; MD 81
    "4200": ((3,), "thread"),  # (gd|ss) 26; next () loop 26; MD 192
    "4210": ((1,), "thread"),  # (gd|ps) 86; next (1, 3) thread 86; MD 350
    "4211": ((1, 3), "thread"),  # (gd|pp) 92; next (2, 3) loop 108; MD 187
    "4220": ((2, 3), "loop"),  # (gd|ds) 121; next (2, 3) thread 127; MD 295
    "4221": ((2, 3), "loop"),  # (gd|dp) 260; next (2, 3) thread 276; MD 271
    "4230": ((2, 3), "loop"),  # (gd|fs) 132; next (2, 3) thread 141; MD 209
    "4300": ((1, 3), "thread"),  # (gf|ss) 38; next (1,) thread 38; MD 192
    "4310": ((1, 3), "thread"),  # (gf|ps) 90; next (1,) thread 91; MD 307
    "4311": ((1, 3), "thread"),  # (gf|pp) 84; next (1, 2, 3) thread 138; MD 141
    "4320": ((1, 3), "thread"),  # (gf|ds) 134; next (1,) thread 135; MD 219
    "4400": ((1, 3), "thread"),  # (gg|ss) 12; next (1,) thread 12; MD 48
    "4410": ((1,), "thread"),  # (gg|ps) 32; next (1, 3) thread 32; MD 70
}

#: The classes up to d whose Fock build goes to the chunked Rys kernel
#: rather than HGP when TRC_RYS_CLASSES is unset (gen_rys_dispatch.py's
#: rys_default_class; trc_rys_select.F90). HGP's heavy d kernels have the
#: f kernels' array-form disease: 255 registers plus local-memory stack,
#: e.g. 3.1 and 1.9 KB per thread in the pcr2221 and pcr2121 warp-role
#: kernels that run (cuobjdump -res-usage).
#:
#: Each of the 21 canonical s/p/d classes was raced (the s/p ones through
#: Rys kernels built for the race only), HGP against the array Rys form and
#: every chunk shape in both forms, on six systems at 1e-13 on a V100: w20
#: (the first 60 atoms of w172) in cc-pVTZ, cc-pVDZ and def2-SVP,
#: gly10/6-31G* (Cartesian), w60 (180 atoms)/def2-SVP and gly3/def2-SVPD.
#: A class is here only if its plan shape beat HGP on all six (repeat runs
#: vary by under 0.5% per class), and on NA-040 (89 atoms)/def2-SVPD too.
#: Per class, HGP -> chunked Rys, ms over the whole run (w20 cc-pVTZ,
#: gly10, w60 def2-SVP, gly3 def2-SVPD) or over three iterations (NA-040):
#:            w20 TZ       gly10       w60 SVP     gly3 SVPD      NA-040
#:   (dd|dp) 1469 -> 452  356 -> 122  296 -> 114  601 -> 212  15576 -> 5073
#:   (dp|dp) 1204 -> 539  321 -> 169  587 -> 310  433 -> 211  10488 -> 5595
#:   (dp|ds)  677 -> 595  278 -> 176  481 -> 442  256 -> 171   5530 -> 4512
#:   (dd|pp)  571 -> 212  202 ->  77  327 -> 161  232 ->  84   4626 -> 1699
#:   (dd|ps)  486 -> 290  170 -> 111  353 -> 241  142 ->  95   3088 -> 1862
#:   (dd|dd)  564 -> 133  176 ->  30   67 ->  15  254 ->  59   5541 -> 1547
#:   (dp|pp)  427 -> 381  186 -> 158  574 -> 561  180 -> 149   3331 -> 2819
#: HGP stays faster on every system for (ss|ss), (ps|ss), (ps|ps), (pp|ss),
#: (pp|ps), (ds|ss), (ds|ps), (dp|ss), (dp|ps), (ds|ds) and (dd|ss), by up
#: to 80% against the best Rys shape; (pp|pp), (dd|ds) and (ds|pp) win on
#: some systems and lose or tie on others, so they stay on HGP too.
FOCK_RYS_DEFAULT = ("2111", "2120", "2121", "2210", "2211", "2221", "2222")

#: Canonical classes above d that keep the array form because no chunk shape
#: beat it (measured, see the FOCK_CHUNK_PLAN comment): (fs|ss), (fs|ps),
#: (fp|ss), (fd|ss).
FOCK_ARRAY_MEASURED = ("3000", "3010", "3100", "3200")

#: CHUNKING IS REQUIRED above d. A canonical class with an f or g shell must
#: be in FOCK_CHUNK_PLAN or FOCK_ARRAY_MEASURED, or main() refuses to
#: generate it: the array form holds the whole Cartesian block per thread,
#: which spills to local memory at these sizes. The g classes below
#: TERCO_RYS_G_LTOT ran 8.92 s array vs 3.17 s chunked ((2,3) thread) on w20
#: cc-pVQZ, V100 -- 2.2x slower than MD in array form, 2.2x faster chunked.
#: A new class gets its shape by racing the shapes (--chunk-classes,
#: --chunk-outer, --chunk-form), not by default.
#:
#: Live-set ceiling, in doubles, for chunk_plan_fock: only the fallback for
#: a class forced with --chunk-classes and not in FOCK_CHUNK_PLAN.
CHUNK_LIVE_MAX = 160


#: Fock pairs the digestion writes, as (shell, shell, the density pair it
#: multiplies, the scale): J_ab with D_cd, J_cd with D_ab, and the four
#: exchange blocks with the opposite pair. Shells are 0..3 = a, b, c, d.
_FOCK_PAIRS = (((0, 1), (2, 3), "cj"), ((2, 3), (0, 1), "cj"),
               ((0, 2), (1, 3), "ck"), ((0, 3), (1, 2), "ck"),
               ((1, 2), (0, 3), "ck"), ((1, 3), (0, 2), "ck"))
_AO = ("mui", "nuj", "lamk", "sigl")
_DN = ("a", "b", "c", "d")


def _ao(s, x):
    if x == "o":
        return f"{_AO[s]} + jo{s}"
    return _AO[s] if x == 0 else f"{_AO[s]} + {x}"


def _dname(p, q, x, y):
    return f"d{_DN[p]}{_DN[q]}_{x}_{y}"


def _sum_lines(terms, pad):
    """`terms` joined by ' + ', wrapped at four terms per line."""
    rows = [" + ".join(terms[i:i + 4]) for i in range(0, len(terms), 4)]
    return (" &\n" + pad + "   + ").join(rows)


def _chunk_digest(chunk, vnames, pad, dnames, which):
    """One chunk's digestion for one density: for each Fock element the chunk
    touches in the pairs `which` (indices into _FOCK_PAIRS), one partial
    sum over the chunk's components that land on it and one atomic. The
    density elements are loaded once into named scalars first."""
    out = []
    need = {}
    for w in which:
        (p, q), (r, s), _ = _FOCK_PAIRS[w]
        for comp in chunk:
            need[(r, s, comp[r], comp[s])] = None
    for (r, s, x, y) in need:
        nmd = _dname(r, s, x, y)
        dnames.add(nmd)
        out.append(f"{pad}{nmd} = dmat(idens, {_ao(r, x)}, {_ao(s, y)})\n")
    for w in which:
        (p, q), (r, s), scale = _FOCK_PAIRS[w]
        groups = {}
        for comp, vn in zip(chunk, vnames):
            groups.setdefault((comp[p], comp[q]), []).append(
                f"{vn}*{_dname(r, s, comp[r], comp[s])}")
        sign = "+" if scale == "cj" else "-"
        for (x, y), terms in groups.items():
            out.append(f"{pad}sc = {_sum_lines(terms, pad)}\n")
            el = f"jmat(idens, {_ao(p, x)}, {_ao(q, y)})"
            out.append(f"{pad}!$acc atomic update\n"
                       f"{pad}{el} = {el} {sign} {scale}*sc\n")
    return "".join(out)


#: The outer sets the chunked pass supports: per side, nothing, the second
#: shell, or both -- never the first shell alone (see the section comment).
_FOCK_OUTER_SETS = [tuple(sorted(br + kt)) for br in ((), (1,), (0, 1)) for kt in ((), (3,), (2, 3))]


def _chunk_keys(comps, outer, axis):
    """Per axis, the (bra key, ket key) of every inside component: the
    exponents of the inside shells on each side, () where a side is all
    outer. Returns the inside component list and the per-component keys."""
    inside = [s for s in range(4) if s not in outer]
    dims = [len(c) for c in comps]
    inner = []
    for ic in itertools.product(*[range(dims[x]) for x in reversed(inside)]):
        full = dict(zip(reversed(inside), ic))
        inner.append(tuple(full.get(x, "o") for x in range(4)))
    keys = []
    for comp in inner:
        e = [None if comp[x] == "o" else comps[x][comp[x]][axis] for x in range(4)]
        bk = tuple(v for v in e[:2] if v is not None)
        kk = tuple(v for v in e[2:] if v is not None)
        keys.append((bk, kk))
    return inner, keys


def _chunk_axis(tag, outer, keys, ax, pad):
    """One axis of a chunk, per root: the 2D-table corner,
    the bra transfer (literal binomials where both bra shells are inside,
    runtime weights wb/wbb where they are not), then the ket transfer the
    same way. Returns (text, names, {(bk, kk): scalar})."""
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    bo = (1 in outer) + (0 in outer)
    ko = (3 in outer) + (2 in outer)
    g = lambda i, k: ri.nm(f"qg{ax}", i, k)
    out, names = [], []
    hkeys = sorted(set(keys))

    # The k range the ket stage reads, per bra key.
    def kneed(kk):
        if ko == 0:
            ec, ed = kk
            return [ec + u for u in range(ed + 1)]
        if ko == 1:
            return [kk[0] + u for u in range(ld + 1)]
        return list(range(nk + 1))
    tneed = sorted({(bk, k) for (bk, kk) in hkeys for k in kneed(kk)})

    def ineed(bk):
        if bo == 0:
            ea, eb = bk
            return [ea + s for s in range(eb + 1)]
        if bo == 1:
            return [bk[0] + s for s in range(lb + 1)]
        return list(range(ni + 1))
    imax = max(i for (bk, k) in tneed for i in ineed(bk))
    kmax = max(k for (bk, k) in tneed)
    g00 = "f00" if ax == "z" else "1.0_dp"
    out.append(ri.rys2d_scalar_assignments(imax, kmax, f"c00{ax}", f"cp0{ax}",
                                           "b00r", "b10r", "bp01r", g00, f"qg{ax}", len(pad)))
    names += ri.rys2d_scalar_names(imax, kmax, f"qg{ax}")

    tname = {}
    for (bk, k) in tneed:
        if bo == 0:
            ea, eb = bk
            if eb == 0:
                tname[(bk, k)] = g(ea, k)
                continue
            terms = []
            for t in range(eb + 1):
                coef, p = ri.binom(eb, t), eb - t
                pre = "" if coef == 1 else f"{float(coef)}_dp*"
                mul = "" if p == 0 else (f"ab{ax}*" if p == 1 else f"ab{ax}{p}*")
                terms.append(f"{pre}{mul}{g(ea + t, k)}")
        elif bo == 1:
            terms = [f"wb{ax}_{t}*{g(bk[0] + t, k)}" for t in range(lb + 1)]
        else:
            terms = [f"wbb{ax}_{i}*{g(i, k)}" for i in range(ni + 1)]
        nmv = ri.nm(f"qt{ax}", *bk, k)
        names.append(nmv)
        tname[(bk, k)] = nmv
        out.append(f"{pad}{nmv} = {_sum_lines(terms, pad)}")
    hmap = {}
    for (bk, kk) in hkeys:
        if ko == 0:
            ec, ed = kk
            if ed == 0:
                hmap[(bk, kk)] = tname[(bk, ec)]
                continue
            terms = []
            for u in range(ed + 1):
                coef, p = ri.binom(ed, u), ed - u
                pre = "" if coef == 1 else f"{float(coef)}_dp*"
                mul = "" if p == 0 else (f"cd{ax}*" if p == 1 else f"cd{ax}{p}*")
                terms.append(f"{pre}{mul}{tname[(bk, ec + u)]}")
        elif ko == 1:
            terms = [f"wd{ax}_{u}*{tname[(bk, kk[0] + u)]}" for u in range(ld + 1)]
        else:
            terms = [f"wdd{ax}_{k}*{tname[(bk, k)]}" for k in range(nk + 1)]
        nmv = ri.nm(f"qh{ax}", *bk, "k", *kk)
        names.append(nmv)
        hmap[(bk, kk)] = nmv
        out.append(f"{pad}{nmv} = {_sum_lines(terms, pad)}")
    return "\n".join(out) + "\n", names, hmap


def chunk_live_fock(tag, comps, outer):
    """The live-set model, in doubles: the chunk's integrals plus, per
    axis, the 2D-table corner, the transfer intermediates and the slice."""
    inner, _ = _chunk_keys(comps, outer, 0)
    live = len(inner)
    for axis in range(3):
        _, keys = _chunk_keys(comps, outer, axis)
        txt, nms, _ = _chunk_axis(tag, outer, keys, "x", "")
        live += len(nms)
    return live


def chunk_plan_fock(tag, comps, live_max=CHUNK_LIVE_MAX):
    """The supported outer set with the fewest chunks under `live_max`."""
    dims = [len(c) for c in comps]
    best = None
    for outer in _FOCK_OUTER_SETS:
        nch = 1
        for x in outer:
            nch *= dims[x]
        live = chunk_live_fock(tag, comps, outer)
        if live > live_max:
            continue
        key = (nch, live)
        if best is None or key < best[0]:
            best = (key, outer)
    return None if best is None else best[1]


def _weight_lines(tag, outer, pad):
    """The runtime transfer weights of one chunk, per axis: wb_s = C(eb, s)
    ab^(eb-s) (0 past eb) for an outer b, and wbb_i = sum_s [ea + s = i] wb_s
    when a is outer too; the same on the ket side. Per chunk, outside the
    primitive loop: they depend on the quartet and the chunk only."""
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    out, names = [], []
    for (first, second, lsec, ltot, w, W, ab) in ((0, 1, lb, ni, "wb", "wbb", "ab"),
                                                   (2, 3, ld, nk, "wd", "wdd", "cd")):
        if second not in outer:
            continue
        for ax in "xyz":
            out.append(f"{pad}select case (ex{second}{ax})")
            for e in range(lsec + 1):
                out.append(f"{pad}case ({e})")
                for t in range(lsec + 1):
                    if t > e:
                        val = "0.0_dp"
                    else:
                        coef, p = ri.binom(e, t), e - t
                        f = [] if coef == 1 else [f"{float(coef)}_dp"]
                        if p == 1:
                            f.append(f"{ab}{ax}")
                        elif p > 1:
                            f.append(f"{ab}{ax}{p}")
                        val = "*".join(f) if f else "1.0_dp"
                    out.append(f"{pad}   {w}{ax}_{t} = {val}")
            out.append(f"{pad}end select")
            names += [f"{w}{ax}_{t}" for t in range(lsec + 1)]
            if first in outer:
                for i in range(ltot + 1):
                    out.append(f"{pad}{W}{ax}_{i} = 0.0_dp")
                for t in range(lsec + 1):
                    for i in range(t, ltot + 1):
                        out.append(f"{pad}if (ex{first}{ax} == {i - t}) {W}{ax}_{i} = {W}{ax}_{i} + {w}{ax}_{t}")
                names += [f"{W}{ax}_{i}" for i in range(ltot + 1)]
    return "\n".join(out) + ("\n" if out else ""), names


def chunked_fock_middle(tag, comps, outer, nroots, pre_loop, prim_hdr, thread=False):
    """The CHUNKED FOCK PASS's item body after the decode (see the section
    comment above): the chunk loop (or, for `thread`, the one chunk `jch`
    this iteration owns), the outer components and their transfer weights,
    the primitive loop (`prim_hdr`, the spliced text from `do k1` through
    the tval line) and roots, and the chunk's digestion. Returns (text, real
    names, integer names, chunk count, the outer shells' cn tables)."""
    la, lb, lc, ld = (int(c) for c in tag)
    dims = [len(c) for c in comps]
    nch = 1
    for x in outer:
        nch *= dims[x]
    inner, _ = _chunk_keys(comps, outer, 0)
    vn = [f"vq_{n}" for n in range(len(inner))]
    names = set(vn) | {"cj", "ck"}
    dnames = set()
    pw = []
    for ax in "xyz":
        for p in range(2, lb + 1):
            pw.append(f"ab{ax}{p} = " + ("*".join([f"ab{ax}"]*p)))
            names.add(f"ab{ax}{p}")
        for p in range(2, ld + 1):
            pw.append(f"cd{ax}{p} = " + ("*".join([f"cd{ax}"]*p)))
            names.add(f"cd{ax}{p}")
    body = [pre_loop,
            "         ! The six Fock scatters' scales, with the degeneracy weight in.\n"
            "         cj = 4.0_dp*jfac*wq\n"
            "         ck = kfac*wq\n"]
    if pw:
        body.append("".join(f"         {x}\n" for x in pw))
    if thread:
        body.append(f"         ! CHUNKED FOCK PASS (scripts/gen_rys.py): chunk jch of {nch}.\n")
    else:
        body.append(f"         ! CHUNKED FOCK PASS (scripts/gen_rys.py): {nch} chunks.\n"
                    f"         do jch = 1, {nch}\n")
    # Outer component indices (0-based) and their exponents, per chunk.
    rem = "jch - 1"
    for n, x in enumerate(outer):
        if n == len(outer) - 1:
            body.append(f"            jo{x} = {rem}\n")
        else:
            body.append(f"            jo{x} = mod({rem}, {dims[x]})\n")
            rem = f"({rem})/{dims[x]}"
        for a, ax in enumerate("xyz"):
            body.append(f"            ex{x}{ax} = cn{ax}{x}_{tag}(jo{x} + 1)\n")
    wtxt, wn = _weight_lines(tag, outer, " "*12)
    body.append(wtxt)
    names.update(wn)
    body.append("".join(f"            {v} = 0.0_dp\n" for v in vn))
    body.append(prim_hdr)
    body.append(ri.roots_inline_body(nroots, "tval", 18))
    body.append(f"""\
                  pref = pref*sqrt(theta)
                  do jrt = 1, {nroots}
                     ! RANGE SEPARATION -- see unrolled_energy_middle's own
                     ! comment (scripts/gen_rys.py): t^2 -> theta*t^2 in u2.
                     u2 = (theta*ru(jrt))/(1.0_dp - theta*ru(jrt))*rho
                     dinv = 1.0_dp/(zeta*eta + u2*(zeta + eta))
                     b00r  = 0.5_dp*u2*dinv
                     b10r  = 0.5_dp*(eta + u2)*dinv
                     bp01r = 0.5_dp*(zeta + u2)*dinv
                     f00 = pref*rw(jrt)
                     c00x = pax - u2*eta*dinv*pqx
                     c00y = pay - u2*eta*dinv*pqy
                     c00z = paz - u2*eta*dinv*pqz
                     cp0x = qcx + u2*zeta*dinv*pqx
                     cp0y = qcy + u2*zeta*dinv*pqy
                     cp0z = qcz + u2*zeta*dinv*pqz
""")
    pad = " "*21
    hmaps, keyss = [], []
    for axis, ax in enumerate("xyz"):
        _, keys = _chunk_keys(comps, outer, axis)
        txt, nms, hmap = _chunk_axis(tag, outer, keys, ax, pad)
        body.append(txt)
        names.update(nms)
        hmaps.append(hmap)
        keyss.append(keys)
    for m in range(len(inner)):
        h = [hmaps[axis][keyss[axis][m]] for axis in range(3)]
        body.append(f"{pad}{vn[m]} = {vn[m]} + {h[0]}*{h[1]}*{h[2]}\n")
    body.append("""\
                  end do
               end do
            end do
            do idens = 1, ndens
""")
    pad = " "*15
    body.append(_chunk_digest(inner, vn, pad, dnames, (0, 1)))
    body.append(f"{pad}! kfac is uniform over the launch: no warp diverges here.\n"
                f"{pad}if (kfac /= 0.0_dp) then\n")
    body.append(_chunk_digest(inner, vn, pad + "   ", dnames, (2, 3, 4, 5)))
    body.append(f"{pad}end if\n            end do\n")
    if not thread:
        body.append("         end do\n")
    tabs = ""
    for x in outer:
        for cnm, a in (("cnx", 0), ("cny", 1), ("cnz", 2)):
            vals = ", ".join(str(c[a]) for c in comps[x])
            tabs += f"   integer, parameter :: {cnm}{x}_{tag}({dims[x]}) = [{vals}]\n"
    ints = [f"jo{x}" for x in outer] + [f"ex{x}{ax}" for x in outer for ax in "xyz"]
    return "".join(body), sorted(names) + sorted(dnames), ints, nch, tabs


def _parse_unroll_max(s):
    """`--unroll-max-comp` takes an int, or 'inf'/'all' for no ceiling --
    the measurement sweep's top end, where every class is unrolled."""
    if s.strip().lower() in ("inf", "all"):
        return math.inf
    return int(s)


def unrolled_energy_middle(tag, la, lb, lc, ld, ni, nk, nroots, dims, order):
    """The unrolled mode's energy middle (see rys_inline.py's UNROLLED MODE
    section): the roots lookup once, then each root as its own block --
    nroots is a compile-time constant per class -- with the 2D recurrence
    and HRR writing named scalars instead of gx/gy/gz/hx/hy/hz arrays, and
    the (aa, bb, cc, dd2) component loop unrolled Python-side so every
    `hx(cnx0(aa), ...)` becomes a literal-indexed named scalar read straight
    off this generator's own `co.components` tables -- no cnx*_tag lookup
    table is even needed for an unrolled class."""
    comps = [co.components(l, order) for l in (la, lb, lc, ld)]
    d0, d1, d2, d3 = dims
    pad = " "*21
    out = [ri.roots_inline_body(nroots, "tval", 18),
          "                  pref = pref*sqrt(theta)\n"]
    for rr in range(1, nroots + 1):
        out.append(f"""\
                  ! RANGE SEPARATION. ru(n) is already the root of the
                  ! quadrature at the ATTENUATED argument tval = theta*T (see
                  ! the primitive loop above), so the Coulomb u2 = t^2/(1-t^2)
                  ! * rho does NOT simply take rho -> theta*rho here the way
                  ! the HGP VRR's oo2ze/rz/re do: u2 is t^2 reparametrised
                  ! into the Laplace variable this root stands for, and the
                  ! substitution that keeps that reparametrisation consistent
                  ! with the SAME theta is t^2 -> theta*t^2 inside u2's own
                  ! formula, not an extra factor of theta outside it (checked
                  ! against pyscf's mol.with_range_coulomb, not derived and
                  ! trusted: the rho*theta analogue of HGP's substitution
                  ! passed every omega == 0 regression and still disagreed
                  ! with pyscf by O(1) once omega > 0). theta*ru({rr}) replaces
                  ! every bare ru({rr}) in the ORIGINAL expression, in the
                  ! ORIGINAL grouping -- omega == 0 makes theta = 1.0_dp
                  ! exactly, reproducing the unattenuated u2 bit-for-bit.
                  u2 = (theta*ru({rr}))/(1.0_dp - theta*ru({rr}))*rho
                  dinv = 1.0_dp/(zeta*eta + u2*(zeta + eta))
                  b00r  = 0.5_dp*u2*dinv
                  b10r  = 0.5_dp*(eta + u2)*dinv
                  bp01r = 0.5_dp*(zeta + u2)*dinv
                  f00 = pref*rw({rr})
                  c00x = pax - u2*eta*dinv*pqx
                  c00y = pay - u2*eta*dinv*pqy
                  c00z = paz - u2*eta*dinv*pqz
                  cp0x = qcx + u2*zeta*dinv*pqx
                  cp0y = qcy + u2*zeta*dinv*pqy
                  cp0z = qcz + u2*zeta*dinv*pqz
""")
        out.append(ri.rys2d_scalar_assignments(ni, nk, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", 18))
        out.append(ri.rys2d_scalar_assignments(ni, nk, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", 18))
        out.append(ri.rys2d_scalar_assignments(ni, nk, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00",    "gz", 18))
        out.append(ri.hrr_scalar_assignments(ni, nk, la, lb, lc, ld, "gx", "abx", "cdx", "hx", 18))
        out.append(ri.hrr_scalar_assignments(ni, nk, la, lb, lc, ld, "gy", "aby", "cdy", "hy", 18))
        out.append(ri.hrr_scalar_assignments(ni, nk, la, lb, lc, ld, "gz", "abz", "cdz", "hz", 18))
        idx = 0
        acc = []
        for dd2 in range(1, d3 + 1):
            for cc in range(1, d2 + 1):
                for bb in range(1, d1 + 1):
                    for aa in range(1, d0 + 1):
                        idx += 1
                        eix, eiy, eiz = comps[0][aa - 1]
                        ejx, ejy, ejz = comps[1][bb - 1]
                        ekx, eky, ekz = comps[2][cc - 1]
                        elx, ely, elz = comps[3][dd2 - 1]
                        hxn = ri.nm("hx", eix, ejx, ekx, elx)
                        hyn = ri.nm("hy", eiy, ejy, eky, ely)
                        hzn = ri.nm("hz", eiz, ejz, ekz, elz)
                        acc.append(f"{pad}vbuf({idx}) = vbuf({idx}) + {hxn}*{hyn}*{hzn}")
        out.append("\n".join(acc) + "\n")
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("class_tag", help="four digits, e.g. 1111 or 2222")
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="must match the --order the scalar HGP file in "
                         "--src-dir was generated with; the digestion it is "
                         "spliced from already scatters into this order, "
                         "and the component tables below must agree with it")
    ap.add_argument("--src-dir", default="src/generated",
                    help="where trc_pc_k{tag}.F90 (the scalar HGP kernel) "
                         "already lives")
    ap.add_argument("--unroll-max-comp", default="0",
                    help="unroll the component and root loops (see "
                         "rys_inline.py's UNROLLED MODE) for a class whose "
                         "na*nb*nc*nd is at most this; 'inf'/'all' unrolls "
                         "every class. Default 0: no class is unrolled, "
                         "i.e. today's array/loop form, unchanged")
    ap.add_argument("--chunk-classes", default="auto",
                    help="classes that get the CHUNKED FOCK PASS: 'auto' "
                         "(FOCK_CHUNK_PLAN's), 'none', or a comma-separated "
                         "list of tags (for measurement)")
    ap.add_argument("--chunk-outer", default=None,
                    help="for measurement: the outer shells, e.g. '2,3', or "
                         "'none' for one chunk, overriding FOCK_CHUNK_PLAN "
                         "(must be one of _FOCK_OUTER_SETS)")
    ap.add_argument("--chunk-form", choices=("loop", "thread"), default=None,
                    help="for measurement: override FOCK_CHUNK_PLAN's form")
    ap.add_argument("--chunk-live-max", type=int, default=CHUNK_LIVE_MAX,
                    help="live-set ceiling chunk_plan_fock picks the shape "
                         "under, for a class FOCK_CHUNK_PLAN does not list")
    ap.add_argument("-o", "--output", required=True,
                    help="required, and never a path under src/generated: "
                         "unlike the Rys root tables, these per-class "
                         "kernels are not committed -- CMake regenerates "
                         "them at configure time into the build directory, "
                         "same as a non-default TERCO_LMAX's HGP kernels")
    args = ap.parse_args()
    tag = args.class_tag
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    nroots = (la + lb + lc + ld)//2 + 1
    dims = [len(co.components(l, args.order)) for l in (la, lb, lc, ld)]
    nv = dims[0]*dims[1]*dims[2]*dims[3]
    comps = [co.components(l, args.order) for l in (la, lb, lc, ld)]
    # chunk_outer: the CHUNKED FOCK PASS's outer shells, None for the array
    # form. () is a real plan (one chunk), so test it with `is None`.
    if args.chunk_classes == "auto":
        chunked = tag in FOCK_CHUNK_PLAN
    elif args.chunk_classes == "none":
        chunked = False
    else:
        chunked = tag in args.chunk_classes.split(",")
    canonical = la >= lb and lc >= ld and (la, lb) >= (lc, ld)
    if (args.chunk_classes == "auto" and canonical and max(la, lb, lc, ld) > 2
            and tag not in FOCK_CHUNK_PLAN and tag not in FOCK_ARRAY_MEASURED):
        raise SystemExit(f"gen_rys.py: ({tag}) has a shell above d and no measured chunk plan; "
                         "race the shapes and add it to FOCK_CHUNK_PLAN (see CHUNKING IS REQUIRED)")
    chunk_outer, chunk_form = None, None
    if chunked:
        chunk_outer, chunk_form = FOCK_CHUNK_PLAN.get(tag, (None, "loop"))
        if args.chunk_outer is not None:
            chunk_outer = () if args.chunk_outer == "none" else \
                tuple(sorted(int(x) for x in args.chunk_outer.split(",")))
        if chunk_outer is None:
            chunk_outer = chunk_plan_fock(tag, comps, args.chunk_live_max)
        assert chunk_outer in _FOCK_OUTER_SETS, chunk_outer
        if args.chunk_form is not None:
            chunk_form = args.chunk_form
    # The chunked form is never unrolled: it already has literal indices.
    unroll = chunk_outer is None and nv <= _parse_unroll_max(args.unroll_max_comp)

    src_path = os.path.join(args.src_dir, f"trc_pc_k{tag}.F90")
    src = open(src_path).read()

    # Scope every search to pcsi{tag} itself (not pcs{tag}, its driver, and
    # not any pcr{tag}/pcri{tag} warp-role kernel the same file may also
    # carry for HGP -- ROLE_SPLITS in gen_perclass.py -- which would also
    # contain a `tval = rho*(pqx...` line and silently splice the wrong one).
    i_item = src.index(f"   pure subroutine pcsi{tag}(")
    i_end = src.index(f"end subroutine pcsi{tag}\n") + len(f"end subroutine pcsi{tag}\n")
    whole = src[:i_end]
    assert i_item < i_end

    tval_line = "tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tval = whole.index(tval_line) + len(tval_line)
    head = whole[:i_tval]

    decode_marker = "         ! Single column on every side:"
    hrr_marker = "         ! --- HRR ---\n"
    dig_marker = "         !\n         ! BLOCK-ACCUMULATED DIGESTION.\n"
    i_decode = whole.index(decode_marker)
    i_hrr = whole.index(hrr_marker)
    i_dig = whole.index(dig_marker)
    decode = whole[i_decode:i_hrr]
    tail = whole[i_dig:]

    # The pre-HRR VRR accumulator's size, read back from its own
    # declaration rather than recomputed -- it is whatever gen_perclass.py
    # actually emitted, and asserting it here catches a mismatch as a clear
    # error instead of a wrong array bound three edits later.
    decl_re = re.compile(
        r"      real\(dp\) :: v\((\d+), 0:1\), g1\(\1\), vbuf\((\d+)\)\n")
    m = decl_re.search(head)
    assert m, "could not find the HGP v/g1/vbuf declaration to replace"
    nv_vrr, nv_in_file = int(m.group(1)), int(m.group(2))
    assert nv_in_file == nv, (
        f"vbuf({nv_in_file}) in {src_path} disagrees with "
        f"na*nb*nc*nd = {nv} computed here for --order {args.order}; "
        "the --order argument likely does not match the source file")

    decl_new = (
        ("" if chunk_outer is not None else f"      real(dp) :: vbuf({nv})\n")
        + f"      real(dp) :: ru({nroots}), rw({nroots})\n"
        f"      real(dp) :: u2, dinv, b00r, b10r, bp01r, f00\n"
        f"      real(dp) :: c00x, c00y, c00z, cp0x, cp0y, cp0z\n"
        + ri.roots_decls())
    if chunk_outer is not None:
        decl_new += "      integer  :: jch, jrt\n@@CHUNK_DECLS@@"
    elif unroll:
        decl_new += ri.unrolled_table_decls(ni, nk, la, lb, lc, ld, 6)
    else:
        decl_new += (
            f"      real(dp) :: gx(0:{ni}, 0:{nk}), gy(0:{ni}, 0:{nk}), gz(0:{ni}, 0:{nk})\n"
            f"      real(dp) :: hx(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})\n"
            f"      real(dp) :: hy(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})\n"
            f"      real(dp) :: hz(0:{la}, 0:{lb}, 0:{lc}, 0:{ld})\n"
            f"      integer  :: rr, aa, bb, cc, dd2\n"
            + ri.rys2d_hrr_decls(la, lb, ld, nk))
    head = head[:m.start()] + decl_new + head[m.end():]

    # Boys scratch, now unused.
    for pat in (r"      real\(dp\) :: f\(0:BOYS_MMAX\)\n",
                r"      integer  :: bi, bj, bbase\n",
                r"      real\(dp\) :: bx, bx2, b0, b1, b2, btt, bet\n"):
        head, n1 = re.subn(pat, "", head, count=1)
        assert n1 == 1, pat

    # HGP zeroes its pre-HRR accumulator once before the primitive loop;
    # Rys zeroes `vbuf` there instead, since nothing defers the transfer.
    zero_re = re.compile(r"(            do x = 1, )" + str(nv_vrr) +
                         r"(\n               )g1\(x\)( = 0\.0_dp\n            end do\n)")
    if chunk_outer is None:
        head, n1 = zero_re.subn(lambda g: f"{g.group(1)}{nv}{g.group(2)}vbuf(x){g.group(3)}", head, 1)
    else:
        # The chunked pass zeroes its own accumulators per chunk; the split
        # point is where the chunk loop goes, and the digestion blocks are
        # replaced by its named scalars.
        head, n1 = zero_re.subn("@@CHUNK_SPLIT@@", head, 1)
        for pat in (r"      real\(dp\) :: jab\(\d+\), jcd\(\d+\), kac\(\d+\)\n",
                    r"      real\(dp\) :: kad\(\d+\), kbc\(\d+\), kbd\(\d+\)\n",
                    r"      real\(dp\) :: dab\(\d+\), dcd\(\d+\), dac\(\d+\)\n",
                    r"      real\(dp\) :: dad\(\d+\), dbc\(\d+\), dbd\(\d+\)\n"):
            head, n2 = re.subn(pat, "", head, count=1)
            assert n2 == 1, pat
    assert n1 == 1, "could not find the g1 zero loop to retarget at vbuf"

    # No Boys function; the root lookup replaces it, inlined (see
    # rys_inline.py) rather than called.
    head, n1 = re.subn(r"   use trc_boys, only:.*?BOYS_DT, BOYS_DTINV\n",
                       ri.roots_table_use(nroots),
                       head, count=1, flags=re.S)
    assert n1 == 1

    def ren(t):
        t = t.replace(f"trc_pc_k{tag}", f"trc_pc_rys{tag}")
        t = t.replace(f"pcsi{tag}", f"rysi{tag}")
        t = t.replace(f"pcs{tag}", f"rys{tag}")
        return t

    head = re.sub(r"public :: pcs" + tag + r".*\n", f"public :: rys{tag}\n", head, count=1)
    head, decode, tail = ren(head), ren(decode), ren(tail)

    # The file-top comment still describes the HGP file it was read from;
    # replace it with one that says what this file actually is.
    old_header_re = re.compile(r"\A!\n.*?#endif\n", re.S)
    assert old_header_re.search(head)
    order_val = 0 if args.order == "lexicographic" else 1
    chunk_note = ""
    if chunk_outer is not None:
        chunk_note = ("\n!\n! CHUNKED FOCK PASS (gen_rys.py): the integral and the digestion are\n"
                      "! replaced, chunk by chunk; the prologue is still spliced.")
    new_header = f"""!
! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) kernel, Rys quadrature.
!
! GENERATED by scripts/gen_rys.py {tag} -- do not edit.
!
! Spliced from the scalar HGP kernel in src/generated/trc_pc_k{tag}.F90: the
! prologue (segment decode, screens, the primitive-pair loads, the
! primitive-quartet prescreen) and the folded six-atomic digestion are that
! file's text, verbatim; only the integral itself differs. See this
! generator's module docstring for the cut points.{chunk_note}
!
! TRC_CART_ORDER GUARD -- see gen_perclass.py's non-split header for the
! full explanation; same guard, carried over from the file this was spliced
! from, since the digestion below scatters into the same AO slot order.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""
    head = old_header_re.sub(new_header, head, count=1)

    # The Rys integral itself: roots, the 2D recurrence per axis, the
    # horizontal transfer per axis, and the outer product straight into
    # vbuf -- no deferred HRR, see the module docstring. Every piece is
    # inlined straight-line code (scripts/rys_inline.py): no CALL survives
    # into the SASS, which is this generator's whole point -- see that
    # module's docstring for why a call boundary, even within the same
    # module, could not be trusted to disappear on its own.
    if unroll:
        middle = "\n" + unrolled_energy_middle(tag, la, lb, lc, ld, ni, nk, nroots, dims, args.order) + \
            "               end do\n            end do\n"
    else:
        middle = f"""
{ri.roots_inline_body(nroots, "tval", 18)}
                  pref = pref*sqrt(theta)
                  do rr = 1, {nroots}
                     ! RANGE SEPARATION -- see unrolled_energy_middle's own
                     ! comment (scripts/gen_rys.py) for the derivation: u2's
                     ! substitution is t^2 -> theta*t^2 inside its own
                     ! formula, not rho -> theta*rho the way HGP's VRR does.
                     u2 = (theta*ru(rr))/(1.0_dp - theta*ru(rr))*rho
                     dinv = 1.0_dp/(zeta*eta + u2*(zeta + eta))
                     b00r  = 0.5_dp*u2*dinv
                     b10r  = 0.5_dp*(eta + u2)*dinv
                     bp01r = 0.5_dp*(zeta + u2)*dinv
                     f00 = pref*rw(rr)
                     c00x = pax - u2*eta*dinv*pqx
                     c00y = pay - u2*eta*dinv*pqy
                     c00z = paz - u2*eta*dinv*pqz
                     cp0x = qcx + u2*zeta*dinv*pqx
                     cp0y = qcy + u2*zeta*dinv*pqy
                     cp0z = qcz + u2*zeta*dinv*pqz

{ri.rys2d_inline_body(ni, nk, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", 21)}
{ri.rys2d_inline_body(ni, nk, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", 21)}
{ri.rys2d_inline_body(ni, nk, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00",    "gz", 21)}
{ri.hrr_inline_body(ni, nk, la, lb, lc, ld, "gx", "abx", "cdx", "hx", 21)}
{ri.hrr_inline_body(ni, nk, la, lb, lc, ld, "gy", "aby", "cdy", "hy", 21)}
{ri.hrr_inline_body(ni, nk, la, lb, lc, ld, "gz", "abz", "cdz", "hz", 21)}

                     idx = 0
                     do dd2 = 1, {dims[3]}
                        do cc = 1, {dims[2]}
                           do bb = 1, {dims[1]}
                              do aa = 1, {dims[0]}
                                 idx = idx + 1
                                 vbuf(idx) = vbuf(idx) &
                                    + hx(cnx0_{tag}(aa), cnx1_{tag}(bb), cnx2_{tag}(cc), cnx3_{tag}(dd2)) &
                                    * hy(cny0_{tag}(aa), cny1_{tag}(bb), cny2_{tag}(cc), cny3_{tag}(dd2)) &
                                    * hz(cnz0_{tag}(aa), cnz1_{tag}(bb), cnz2_{tag}(cc), cnz3_{tag}(dd2))
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
"""

    # Component tables: aa/bb/cc/dd2 above index into these to pick up each
    # shell's (lx,ly,lz) triple in the BUILD'S OWN Cartesian order -- the
    # same order gen_perclass.py's digestion scattered into when it wrote
    # the tail this file reused verbatim. An unrolled class needs none of
    # this: the component triples are read straight off `co.components` at
    # GENERATION time instead (see unrolled_energy_middle).
    if chunk_outer is not None:
        pre, prim_hdr = head.split("@@CHUNK_SPLIT@@")
        pre_loop = decode + """\
         ! Degeneracy weight, as the block digestion forms it.
         wq = 1.0_dp
         if (.not. dij) wq = wq*0.5_dp
         if (.not. dkl) wq = wq*0.5_dp
         if (.not. dpq) wq = wq*0.5_dp
"""
        thread = chunk_form == "thread"
        middle, cnames, cints, nch, ctabs = chunked_fock_middle(
            tag, comps, chunk_outer, nroots, pre_loop, prim_hdr, thread)
        cdecl = ri.scalar_decl_lines(cnames, 6)
        if cints:
            cdecl += f"      integer  :: {', '.join(cints)}\n"
        head = pre.replace("@@CHUNK_DECLS@@", cdecl)
        if thread:
            # One thread per (quartet, chunk), quartet fastest, so a warp's
            # lanes share one chunk and its select case does not diverge.
            for a, b in (("      do concurrent(i=1:nr)\n", f"      do concurrent(i=1:nr*{nch})\n"),
                         (f"call rysi{tag}(g0 + (i - 1)*nranks, ",
                          f"call rysi{tag}(g0 + mod(i - 1, nr)*nranks, int((i - 1)/nr) + 1, "),
                         (f"pure subroutine rysi{tag}(gt, ", f"pure subroutine rysi{tag}(gt, jch, "),
                         ("      integer  :: jch, jrt\n", "      integer  :: jrt\n")):
                assert head.count(a) == 1, a
                head = head.replace(a, b)
            a = "      !$acc routine seq\n      integer(int64), intent(in) :: gt\n"
            assert head.count(a) == 1
            head = head.replace(a, a + "      integer, intent(in) :: jch\n")
        decode, tail = "", f"   end subroutine rysi{tag}\n"

    tabs = ctabs if chunk_outer is not None else ""
    if not unroll and chunk_outer is None:
        for cnm, ax in (("cnx", 0), ("cny", 1), ("cnz", 2)):
            for which, l in enumerate((la, lb, lc, ld)):
                comps = co.components(l, args.order)
                vals = ", ".join(str(c[ax]) for c in comps)
                tabs += f"   integer, parameter :: {cnm}{which}_{tag}({len(comps)}) = [{vals}]\n"

    # `ren` renamed the module's own "contains" line nowhere -- it is bare
    # text -- so insert the component tables ahead of the module's original
    # `contains` (the one right after its TWO_PI_2_5 parameter, in the
    # specification part, before any subroutine; the tables are `parameter`
    # arrays and may not follow one). No helper subroutines follow them any
    # more -- the roots lookup, 2D recurrence and HRR are inlined straight
    # into `middle` below instead of called, so `contains` goes right back
    # to introducing rys{tag}/rysi{tag}, same as the HGP file this was
    # spliced from.
    helpers = f"""
{tabs}
contains
"""
    assert head.count("\ncontains\n") == 1
    head = head.replace("\ncontains\n", helpers, 1)

    out = head + middle + decode + tail + f"\nend module trc_pc_rys{tag}\n"
    dest = args.output
    with open(dest, "w") as fh:
        fh.write(out)
    print(f"  {dest}: ({la}{lb}|{lc}{ld}), {nroots} roots, vbuf({nv}), {len(out.splitlines())} lines")


if __name__ == "__main__":
    main()
