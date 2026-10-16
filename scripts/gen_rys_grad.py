#!/usr/bin/env python3
"""Emit a Rys GRADIENT kernel for one class, by splicing terco's own scalar
HGP kernel -- the same technique gen_rys.py uses for the energy kernel, with
the same reasons (see that file's docstring for why splicing rather than
generating whole).

WHAT IS DIFFERENT FROM THE ENERGY SPLICE
-----------------------------------------
The energy kernel defers the density fold to AFTER the primitive/root loop:
it sums the bare integral into `vbuf` and only then multiplies by the density
blocks, once per quartet. The gradient cannot do that -- the derivative
weight `fix`/`fjx`/`fkx`/`flx` depends on the PRIMITIVE's own exponent, so it
has to be formed, and folded against the density, INSIDE the primitive/root
loop. That means the density blocks have to be loaded BEFORE that loop
starts rather than after it ends, so the HEAD/DECODE/TAIL split gen_rys.py
uses does not carry over unchanged:

  - HEAD is still the segment decode, screens, bin-interleaved primitive
    loads and primitive-quartet prescreen, verbatim -- but with the density
    blocks, the degeneracy weight `wq`, and a precomputed `gam(idx)` per
    output Cartesian component folded in right where the HGP file already
    has every piece they need (si/sj/sk/sl, same_ab/same_cd/same_pair, the
    AO offsets) and right before the primitive loop starts. This is new
    text, not copied from HGP's tail -- the energy digestion's "Single
    column..." decode and block digestion are not reused at
    all.
  - MIDDLE replaces the Rys integral build exactly like the energy splice,
    except it builds ONE extended horizontal-transfer table per axis,
    `h(0:la+1, 0:lb, 0:lc+1, 0:ld)` -- i and k (the bra- and ket-combined
    totals) raised by one each, j and l never raised -- and folds the
    derivative directly into 12 running sums (A, B, C, D x, y, z) instead of
    an output integral block.
  - TAIL is new too: four atomic updates into `grad(3, natm)`, through the
    shell -> atom map `sh_at`, in place of HGP's six-block Fock digestion.

THE ALGEBRA (GPU4PySCF's, not re-derived)
------------------------------------------
Checked against gpu4pyscf's actual production kernel
(gvhf-rys/rys_contract_jk_ip1.cu, the BUILD_4C_GXYZ/fix/fjx/fkx/flx block),
not only its description:

  - Raise ONLY the bra- and ket-combined 2D-integral indices (i, k), never
    j or l. The j derivative comes from the horizontal-transfer identity
    instead of a third raised table:
        fix = 2*alpha_i*h(i+1,j,k,l) [- i*h(i-1,j,k,l) if i>0]
        fkx = 2*alpha_k*h(i,j,k+1,l) [- k*h(i,j,k-1,l) if k>0]
        fjx = 2*alpha_j*(h(i+1,j,k,l) + (Ax-Bx)*h(i,j,k,l)) [- j*h(i,j-1,k,l) if j>0]
        flx = -fix - fjx - fkx                      (translational invariance)
  - nroots = (L + 1)//2 + 1 with L = la+lb+lc+ld (integer division), NOT
    (L+2)//2+1: every derivative TERM raises only one of i or k, never both
    in the same product, so the polynomial degree actually summed is L+1,
    not L+2, even though the table's nominal range covers both raised by
    one. Confirmed against gpu4pyscf's own bounds (gvhf-rys/...ip1.cu
    around its `order = li+lj+lk+ll; nroots = (order+1)/2+1`). This means
    the existing n = 1..5 root tables cover every class up to (dd|dd); no
    new table is needed.
  - Gamma, from the standard RHF/UHF two-electron gradient identity (not
    copied from the Fock digestion's internal 4x/divide-by-8 bookkeeping,
    which is a Fock-MATRIX convention that does not transfer to a gradient
    SCALAR):
        RHF (ndens=1, dmat(:,:,1) = total D):
            Gamma = jfac*D_ab*D_cd - kfac/4*(D_ac*D_bd + D_ad*D_bc)
        UHF (ndens=2, dmat(:,:,1)=Da, dmat(:,:,2)=Db):
            Gamma = jfac*(Da+Db)_ab*(Da+Db)_cd
                    - kfac*(Da_ac*Da_bd+Da_ad*Da_bc + Db_ac*Db_bd+Db_ad*Db_bc)
    jfac/kfac default to 1 (pure HF); a hybrid functional's exchange
    fraction multiplies kfac, same meaning as eri_fock's k_scale.
  - A quartet whose four shells all sit on the same atom contributes
    nothing (A=B=C=D already satisfies translational invariance on its
    own) and is skipped before the primitive loop.

Never hand-edit the output; change this generator and regenerate.
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


def _parse_unroll_max(s):
    """Same parsing as gen_rys.py's --unroll-max-comp: an int, or
    'inf'/'all' for no ceiling."""
    if s.strip().lower() in ("inf", "all"):
        return math.inf
    return int(s)


def unrolled_gam_block(d0, d1, d2, d3):
    """The unrolled mode's twin of the density-block/Gamma fill (see the
    module docstring's inject text, just below `kfac_eff`): same algebra,
    but every `dab`/`dcd`/`dac`/`dad`/`dbc`/`dbd`/`gam` array element
    becomes a named scalar (`nm`), and the `do idens = 1, ndens` loop --
    ndens is a RUNTIME value, but only ever 1 (RHF) or 2 (UHF) -- is
    unrolled into the idens=1 term (always present) plus an `if (ndens ==
    2)` for the idens=2 term, rather than keeping a 2-element array indexed
    by a loop variable. `dab`/`dcd` already summed over idens in the array
    form (they hold the TOTAL density), so they keep that shape; `dac`/
    `dad`/`dbc`/`dbd` keep one spin apart (`dac1_*`/`dac2_*`, ...), since
    Gamma's K term needs each spin's own product, not their sum."""
    lines = []
    for ib in range(d1):
        for ia in range(d0):
            v = ri.nm("dab", ia, ib)
            lines.append(f"         {v} = dmat(1, mui + {ia}, nuj + {ib})")
            lines.append(f"         if (ndens == 2) {v} = {v} + dmat(2, mui + {ia}, nuj + {ib})")
    for idd in range(d3):
        for ic in range(d2):
            v = ri.nm("dcd", ic, idd)
            lines.append(f"         {v} = dmat(1, lamk + {ic}, sigl + {idd})")
            lines.append(f"         if (ndens == 2) {v} = {v} + dmat(2, lamk + {ic}, sigl + {idd})")
    for ic in range(d2):
        for ia in range(d0):
            v1, v2 = ri.nm("dac1", ia, ic), ri.nm("dac2", ia, ic)
            lines.append(f"         {v1} = dmat(1, mui + {ia}, lamk + {ic})")
            lines.append(f"         if (ndens == 2) {v2} = dmat(2, mui + {ia}, lamk + {ic})")
    for idd in range(d3):
        for ia in range(d0):
            v1, v2 = ri.nm("dad1", ia, idd), ri.nm("dad2", ia, idd)
            lines.append(f"         {v1} = dmat(1, mui + {ia}, sigl + {idd})")
            lines.append(f"         if (ndens == 2) {v2} = dmat(2, mui + {ia}, sigl + {idd})")
    for ic in range(d2):
        for ib in range(d1):
            v1, v2 = ri.nm("dbc1", ib, ic), ri.nm("dbc2", ib, ic)
            lines.append(f"         {v1} = dmat(1, nuj + {ib}, lamk + {ic})")
            lines.append(f"         if (ndens == 2) {v2} = dmat(2, nuj + {ib}, lamk + {ic})")
    for idd in range(d3):
        for ib in range(d1):
            v1, v2 = ri.nm("dbd1", ib, idd), ri.nm("dbd2", ib, idd)
            lines.append(f"         {v1} = dmat(1, nuj + {ib}, sigl + {idd})")
            lines.append(f"         if (ndens == 2) {v2} = dmat(2, nuj + {ib}, sigl + {idd})")
    for idd in range(d3):
        for ic in range(d2):
            for ib in range(d1):
                for ia in range(d0):
                    g = ri.nm("gam", ia, ib, ic, idd)
                    dabn, dcdn = ri.nm("dab", ia, ib), ri.nm("dcd", ic, idd)
                    dac1, dac2 = ri.nm("dac1", ia, ic), ri.nm("dac2", ia, ic)
                    dad1, dad2 = ri.nm("dad1", ia, idd), ri.nm("dad2", ia, idd)
                    dbc1, dbc2 = ri.nm("dbc1", ib, ic), ri.nm("dbc2", ib, ic)
                    dbd1, dbd2 = ri.nm("dbd1", ib, idd), ri.nm("dbd2", ib, idd)
                    lines.append(f"         {g} = 8.0_dp*wq*jfac*{dabn}*{dcdn} &")
                    lines.append(f"              - 8.0_dp*wq*kfac_eff*({dac1}*{dbd1} + {dad1}*{dbc1})")
                    lines.append(f"         if (ndens == 2) {g} = {g} - 8.0_dp*wq*kfac_eff*({dac2}*{dbd2} + {dad2}*{dbc2})")
    return "\n".join(lines) + "\n"


def unrolled_gam_decls(d0, d1, d2, d3, indent=6):
    """Declarations for `unrolled_gam_block`'s named scalars."""
    names = ([ri.nm("dab", ia, ib) for ib in range(d1) for ia in range(d0)]
             + [ri.nm("dcd", ic, idd) for idd in range(d3) for ic in range(d2)]
             + [ri.nm(p, ia, ic) for p in ("dac1", "dac2") for ic in range(d2) for ia in range(d0)]
             + [ri.nm(p, ia, idd) for p in ("dad1", "dad2") for idd in range(d3) for ia in range(d0)]
             + [ri.nm(p, ib, ic) for p in ("dbc1", "dbc2") for ic in range(d2) for ib in range(d1)]
             + [ri.nm(p, ib, idd) for p in ("dbd1", "dbd2") for idd in range(d3) for ib in range(d1)]
             + [ri.nm("gam", ia, ib, ic, idd)
                for idd in range(d3) for ic in range(d2) for ib in range(d1) for ia in range(d0)])
    return ri.scalar_decl_lines(names, indent)


def _axis_derivative(e_i, e_j, e_k, e_l, hprefix, ab_name):
    """One axis's vix/fix/fkx/fjx (GPU4PySCF's BUILD_4C_GXYZ identity, see
    this generator's module docstring), with the exponents e_i/e_j/e_k
    literal: `real(e, dp)` folds to its value, and a `- e*h(e-1)` term is
    simply never emitted when that exponent is 0, rather than guarded by a
    runtime `if` the way the loop form (gen_rys_grad.py's main `middle`)
    still is."""
    h = lambda a, b, c, d: ri.nm(hprefix, a, b, c, d)
    fi = f"2.0_dp*alpha_i*{h(e_i + 1, e_j, e_k, e_l)}"
    if e_i > 0:
        fi += f" - {float(e_i)}_dp*{h(e_i - 1, e_j, e_k, e_l)}"
    fk = f"2.0_dp*alpha_k*{h(e_i, e_j, e_k + 1, e_l)}"
    if e_k > 0:
        fk += f" - {float(e_k)}_dp*{h(e_i, e_j, e_k - 1, e_l)}"
    fj = f"2.0_dp*alpha_j*({h(e_i + 1, e_j, e_k, e_l)} + {ab_name}*{h(e_i, e_j, e_k, e_l)})"
    if e_j > 0:
        fj += f" - {float(e_j)}_dp*{h(e_i, e_j - 1, e_k, e_l)}"
    return h(e_i, e_j, e_k, e_l), fi, fk, fj


def unrolled_gradient_middle(tag, la, lb, lc, ld, nig, nkg, nroots, dims, order):
    """The unrolled mode's gradient middle (see rys_inline.py's UNROLLED
    MODE section and gen_rys.py's unrolled_energy_middle, its energy twin):
    the roots lookup once, each root as its own block, the raised 2D
    recurrence and HRR writing named scalars, and the (aa, bb, cc, dd2)
    derivative fold unrolled Python-side -- every `hx(eix, ...)` read
    becomes a literal-indexed named scalar, and a derivative term that
    GPU4PySCF's `- e*h(e-1)` would multiply by a literal exponent of 0 is
    simply never emitted, rather than guarded by a runtime `if`."""
    comps = [co.components(l, order) for l in (la, lb, lc, ld)]
    d0, d1, d2, d3 = dims
    pad = " "*21
    out = [ri.roots_inline_body(nroots, "tval", 18)]
    for rr in range(1, nroots + 1):
        out.append(f"""\
                  u2 = ru({rr})/(1.0_dp - ru({rr}))*rho
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
        out.append(ri.rys2d_scalar_assignments(nig, nkg, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", 18))
        out.append(ri.rys2d_scalar_assignments(nig, nkg, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", 18))
        out.append(ri.rys2d_scalar_assignments(nig, nkg, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00",    "gz", 18))
        out.append(ri.hrr_scalar_assignments(nig, nkg, la + 1, lb, lc + 1, ld, "gx", "abx", "cdx", "hx", 18))
        out.append(ri.hrr_scalar_assignments(nig, nkg, la + 1, lb, lc + 1, ld, "gy", "aby", "cdy", "hy", 18))
        out.append(ri.hrr_scalar_assignments(nig, nkg, la + 1, lb, lc + 1, ld, "gz", "abz", "cdz", "hz", 18))

        for dd2 in range(1, d3 + 1):
            for cc in range(1, d2 + 1):
                for bb in range(1, d1 + 1):
                    for aa in range(1, d0 + 1):
                        eix, eiy, eiz = comps[0][aa - 1]
                        ejx, ejy, ejz = comps[1][bb - 1]
                        ekx, eky, ekz = comps[2][cc - 1]
                        elx, ely, elz = comps[3][dd2 - 1]
                        gnm = ri.nm("gam", aa - 1, bb - 1, cc - 1, dd2 - 1)
                        blk = []

                        vix, fix_e, fkx_e, fjx_e = _axis_derivative(eix, ejx, ekx, elx, "hx", "abx")
                        viy, fiy_e, fky_e, fjy_e = _axis_derivative(eiy, ejy, eky, ely, "hy", "aby")
                        viz, fiz_e, fkz_e, fjz_e = _axis_derivative(eiz, ejz, ekz, elz, "hz", "abz")
                        blk.append(f"{pad}vix = {vix}; viy = {viy}; viz = {viz}")
                        blk.append(f"{pad}pyz = viy*viz; pxz = vix*viz; pxy = vix*viy")
                        blk.append(f"{pad}fix = {fix_e}")
                        blk.append(f"{pad}fkx = {fkx_e}")
                        blk.append(f"{pad}fjx = {fjx_e}")
                        blk.append(f"{pad}flx = -fix - fjx - fkx")
                        blk.append(f"{pad}fiy = {fiy_e}")
                        blk.append(f"{pad}fky = {fky_e}")
                        blk.append(f"{pad}fjy = {fjy_e}")
                        blk.append(f"{pad}fly = -fiy - fjy - fky")
                        blk.append(f"{pad}fiz = {fiz_e}")
                        blk.append(f"{pad}fkz = {fkz_e}")
                        blk.append(f"{pad}fjz = {fjz_e}")
                        blk.append(f"{pad}flz = -fiz - fjz - fkz")
                        blk.append(f"{pad}vax = vax + fix*pyz*{gnm}; vay = vay + fiy*pxz*{gnm}; vaz = vaz + fiz*pxy*{gnm}")
                        blk.append(f"{pad}vbx = vbx + fjx*pyz*{gnm}; vby = vby + fjy*pxz*{gnm}; vbz = vbz + fjz*pxy*{gnm}")
                        blk.append(f"{pad}vcx = vcx + fkx*pyz*{gnm}; vcy = vcy + fky*pxz*{gnm}; vcz = vcz + fkz*pxy*{gnm}")
                        blk.append(f"{pad}vdx = vdx + flx*pyz*{gnm}; vdy = vdy + fly*pxz*{gnm}; vdz = vdz + flz*pxy*{gnm}")
                        out.append("\n".join(blk) + "\n")
    return "".join(out)


def _same_atom_check(term):
    """The same-atom skip, exactly as injected into the fused item routine
    (see the `inject` string in `main`), with `return` swapped for
    whatever control-flow word fits the caller -- `cycle` inside the
    screen pass's chunk loop, `return` (the default, unused here) inside
    a `pure` item routine."""
    return f"""\
         ! A quartet whose four shells sit on one atom is translation-
         ! invariant on its own: dGrad/dR = 0 for every R trivially, so skip
         ! it before the primitive loop rather than compute and cancel it.
         if (sh_at(si) == sh_at(sj) .and. sh_at(si) == sh_at(sk) &
             .and. sh_at(si) == sh_at(sl)) {term}
"""


def _return_to(s, term, n):
    """Swap every `...) return` for `...) term` in a decode/screen
    fragment sliced out of the fused item routine -- `cycle` turns a
    `return`-from-item screen into a `cycle`-this-work-item screen inside
    the chunk loop below. `n` is the expected count, asserted so a future
    change to the spliced text cannot silently change how many screens
    the chunk loop actually applies."""
    s2, cnt = re.subn(r"\) return\n", f") {term}\n", s)
    assert cnt == n, (n, cnt, s)
    return s2


# ---------------------------------------------------------------------------
# CHUNKED RHF COMPACT PASS
# ---------------------------------------------------------------------------
#
# WHY. The heavy classes (three or more non-s shells) were bound by DRAM
# traffic to their OWN per-thread stack, not by arithmetic: one thread held
# the whole raised HRR table hx/hy/hz, the whole Gamma(nv) and the density
# blocks across the primitive loop, far beyond the 255-register file, and
# with ~20k resident threads a V100 has ~300 B of L2 per thread, so every
# spilled element streamed through DRAM. Nsight Compute, w150 6-31G*,
# 1e-13: (pp|pp) 640 GB/s of DRAM (71% of peak), (dp|pp) 400 GB/s, 75-80%
# of stall cycles on L1TEX long scoreboard, 2848 B / 3616 B of stack.
#
# WHAT. For RHF (ndens = 1, known at compile time in this routine) the
# compact pass of every class chunk_auto selects -- the heavy ones, and the
# light ones too, which gain from the RHF and D trims alone -- becomes
# rysgi{tag}_crhf, ordered
#
#     chunk of output components -> primitive quartet -> root
#
# A chunk is every Cartesian component of the shells in `inside` for one
# fixed component of each other shell. Per chunk: its Gamma values, read
# straight from dmat (no materialised gam(nv), no density blocks held);
# then the primitive loop, the roots, and per root ONLY the 2D-table corner
# and the HRR slice that chunk's derivative terms read (hrr_slice_
# assignments in rys_inline.py), all named scalars with literal indices.
# The chunk loop is a runtime `do` with a `select case` on the chunk (each
# case keeps its indices literal), and the root loop is a runtime `do`:
# unrolling both took (pp|pp) past nvfortran's 590 s per-object budget
# (27 chunks x 3 roots: no object after 590 s), while the runtime form
# compiles in ~2 min and keeps the stack at 112-960 B (only ru/rw, read with
# the runtime root index, and a little spill). D comes from translational
# invariance once per quartet, vd = -(va + vb + vc), instead of a fourth set
# of accumulators.
#
# The price is recomputation: the primitive-quartet setup and the roots are
# redone once per chunk. So chunks are as FEW as possible while their live
# set stays small; see chunk_plan. UHF (ndens = 2) keeps the old routine.
#
# Order matters. With the chunk loop INSIDE the primitive loop, nvfortran
# hoists every Gamma (loop-invariant there) out of the primitive loops and
# rebuilds the materialised gam(nv) in registers: (pp|pp) went to 5936 B
# of stack and 0.83 s, slower than the 0.63 s it started from.
#
#: Live-set ceiling, in doubles, for one chunk: Gamma values + the 2D-table
#: box + the HRR slice (chunk_live). Measured, w150 6-31G*, 1e-13, V100,
#: compact pass, (chunks x components per chunk, live) -> time:
#:   (pp|pp) 27x3 (50) 0.467 s; 9x9 (75) 0.234 s; 3x27 (115) 0.129 s;
#:           1x81 (225) 0.226 s. Old form 0.635 s.
#:   (dp|pp) 27x6 (69) 0.738 s; 9x18 (101-106) 0.37 s; 6x27 (147) 0.33 s;
#:           3x54 (170) 0.463 s. Old form 1.55 s.
#:   (ds|pp) 3x18 (91) 0.251 s; 1x54 (168) 0.306 s.
#:   (dp|ps) 3x18 (96) 0.430 s; 1x54 (165) 0.433 s.
#:   (dd|ds) 36x6 (77) 0.268 s; 6x36 (155) 0.115 s.
#:   (dd|pp) 18x18 (110) 0.152 s; 9x36 (156) 0.136 s.
#: Fewer chunks win up to a live set of ~156; at 165-170 the gain is gone
#: or reversed, and at 225 the chunk spills like the old form.
CHUNK_LIVE_MAX = 160

#: Largest class (na*nb*nc*nd) chunk_auto chunks: compile time, not speed.
#: The chunked routine is ~10 lines of fold per component; one nvfortran
#: 26.5 object took 187-206 s at 324 components ((dp|dp), (dd|pp)) and did
#: not finish in 590 s at 648 ((dd|dp), 12.4 K lines), the per-object
#: budget of terco's builds. (dd|dp) and (dd|dd) keep the old form.
CHUNK_MAX_COMP = 324


def _chunk_h_needs(comps, chunk, axis):
    """The h(i, j, k, l) elements one axis of the derivative fold reads for
    the components in `chunk` -- exactly the terms _axis_derivative emits:
    h(e), h(ei+1), h(ei-1) if ei > 0, h(ek+1), h(ek-1) if ek > 0, h(ej-1)
    if ej > 0."""
    need = set()
    for (a, b, c, d) in chunk:
        ei, ej, ek, el = (comps[0][a][axis], comps[1][b][axis],
                          comps[2][c][axis], comps[3][d][axis])
        need.add((ei, ej, ek, el))
        need.add((ei + 1, ej, ek, el))
        need.add((ei, ej, ek + 1, el))
        if ei > 0:
            need.add((ei - 1, ej, ek, el))
        if ek > 0:
            need.add((ei, ej, ek - 1, el))
        if ej > 0:
            need.add((ei, ej - 1, ek, el))
    return need


def _g_box(hneed):
    """The 2D-table corner g(0:imax, 0:kmax) an HRR slice reads. The
    recurrence fills a corner from itself alone, so the corner is all that
    has to be built."""
    return (max(i + j for (i, j, k, l) in hneed),
            max(k + l for (i, j, k, l) in hneed))


def _chunks(comps, inside):
    """The chunks for a given `inside` set: one list of (a, b, c, d) per
    combination of the OTHER shells' components."""
    dims = [len(c) for c in comps]
    outer = [s for s in range(4) if s not in inside]
    res = []
    for oc in itertools.product(*[range(dims[s]) for s in reversed(outer)]):
        fixed = dict(zip(reversed(outer), oc))
        chunk = []
        for ic in itertools.product(*[range(dims[s]) for s in reversed(sorted(inside))]):
            full = {**fixed, **dict(zip(reversed(sorted(inside)), ic))}
            chunk.append((full[0], full[1], full[2], full[3]))
        res.append(chunk)
    return res


def chunk_live(comps, chunk):
    """Doubles one chunk keeps live per root: its Gamma values, the 2D-table
    box and the HRR slice, summed over the three axes. A model, calibrated
    against the measurements at CHUNK_LIVE_MAX, not a register count."""
    live = len(chunk)
    for axis in range(3):
        hneed = _chunk_h_needs(comps, chunk, axis)
        imax, kmax = _g_box(hneed)
        live += (imax + 1)*(kmax + 1) + len(hneed)
    return live


def chunk_auto(tag):
    """Whether a class's RHF compact pass is chunked by default: every class
    the pair bins can reach (la >= lb, lc >= ld, (la, lb) >= (lc, ld) -- see
    trc_bins.F90's ANGULAR-MOMENTUM CANONICAL; the other 60 tags are built
    but never launched, so chunking them would only cost compile time) up
    to CHUNK_MAX_COMP components. Measured on every one of them, w150
    6-31G*, 1e-13, V100, compact pass: none got slower; the heavy classes
    gained 50-80% ((pp|pp) 0.635 -> 0.127 s, (dp|pp) 1.56 -> 0.33 s, (dp|ds)
    1.47 -> 0.29 s, (dd|ps) 1.05 -> 0.19 s), and even the one-chunk light
    classes gained from the RHF specialisation and D-by-invariance alone
    ((ps|ss) 1.19 -> 1.05 s, (ps|ps) 0.54 -> 0.40 s, (ss|ss) 0.72 -> 0.69 s)."""
    la, lb, lc, ld = (int(c) for c in tag)
    canonical = la >= lb and lc >= ld and (la, lb) >= (lc, ld)
    nv = 1
    for l in (la, lb, lc, ld):
        nv *= (l + 1)*(l + 2)//2
    return canonical and nv <= CHUNK_MAX_COMP


def chunk_plan(comps, live_max=CHUNK_LIVE_MAX):
    """The `inside` set with the FEWEST chunks whose worst chunk stays at or
    under `live_max`; ties go to the smaller live set. Fewest chunks, since
    every chunk redoes the primitive setup and the roots."""
    best = None
    for r in range(1, 5):
        for inside in itertools.combinations(range(4), r):
            chunks = _chunks(comps, inside)
            live = max(chunk_live(comps, ch) for ch in chunks)
            if live > live_max:
                continue
            key = (len(chunks), live)
            if best is None or key < best[0]:
                best = (key, inside)
    return None if best is None else best[1]


def _gamma_from_dmat(a, b, c, d):
    """Gamma for one component, RHF, read straight from dmat:
    cj*D_ab*D_cd - ck*(D_ac*D_bd + D_ad*D_bc), with cj = 8 wq jfac and
    ck = 8 wq kfac/4 (see the item routine's comments on the 8 and the 1/4)."""
    D = lambda o0, x, o1, y: f"dmat(1, {o0} + {x}, {o1} + {y})"
    return (f"cj*{D('mui', a, 'nuj', b)}*{D('lamk', c, 'sigl', d)} &\n"
            f"                    - ck*({D('mui', a, 'lamk', c)}*{D('nuj', b, 'sigl', d)} &\n"
            f"                    + {D('mui', a, 'sigl', d)}*{D('nuj', b, 'lamk', c)})")


def _chunk_root_body(comps, lb, ld, chunk, pad, names):
    """One root's work for one chunk: the 2D-table corner per axis, the HRR
    slice, and the derivative fold of every component in the chunk into
    the nine A/B/C accumulators."""
    out, hmaps = [], []
    for axis, ax in enumerate("xyz"):
        hneed = _chunk_h_needs(comps, chunk, axis)
        imax, kmax = _g_box(hneed)
        g00 = "f00" if ax == "z" else "1.0_dp"
        out.append(ri.rys2d_scalar_assignments(imax, kmax, f"c00{ax}", f"cp0{ax}",
                                               "b00r", "b10r", "bp01r", g00,
                                               f"qg{ax}", len(pad)))
        names.update(ri.rys2d_scalar_names(imax, kmax, f"qg{ax}"))
        txt, hn, hmap = ri.hrr_slice_assignments(hneed, f"qg{ax}", f"qt{ax}", f"qh{ax}",
                                                 f"ab{ax}", f"cd{ax}", len(pad))
        out.append(txt)
        names.update(hn)
        hmaps.append(hmap)
    for n, (a, b, c, d) in enumerate(chunk):
        v, fi, fj, fk = [], [], [], []
        for axis, ax in enumerate("xyz"):
            ei, ej, ek, el = (comps[0][a][axis], comps[1][b][axis],
                              comps[2][c][axis], comps[3][d][axis])
            H = lambda i, j, k, l: hmaps[axis][(i, j, k, l)]
            v.append(H(ei, ej, ek, el))
            e = f"ai2*{H(ei + 1, ej, ek, el)}"
            if ei > 0:
                e += f" - {float(ei)}_dp*{H(ei - 1, ej, ek, el)}"
            fi.append(e)
            e = f"aj2*({H(ei + 1, ej, ek, el)} + ab{ax}*{H(ei, ej, ek, el)})"
            if ej > 0:
                e += f" - {float(ej)}_dp*{H(ei, ej - 1, ek, el)}"
            fj.append(e)
            e = f"ak2*{H(ei, ej, ek + 1, el)}"
            if ek > 0:
                e += f" - {float(ek)}_dp*{H(ei, ej, ek - 1, el)}"
            fk.append(e)
        out.append(f"{pad}gyz = {v[1]}*{v[2]}*gm_{n}; gxz = {v[0]}*{v[2]}*gm_{n}; "
                   f"gxy = {v[0]}*{v[1]}*gm_{n}\n")
        for axis, (ax, gg) in enumerate((("x", "gyz"), ("y", "gxz"), ("z", "gxy"))):
            out.append(f"{pad}va{ax} = va{ax} + ({fi[axis]})*{gg}\n")
            out.append(f"{pad}vb{ax} = vb{ax} + ({fj[axis]})*{gg}\n")
            out.append(f"{pad}vc{ax} = vc{ax} + ({fk[axis]})*{gg}\n")
    return "".join(out)


#: Declarations of the item routine that the chunked routine never reads:
#: the old full tables, density blocks and materialised Gamma, in both the
#: unrolled (named-scalar) and the loop (array) forms.
_CHUNK_DEAD_DECL = re.compile(
    r"^(gx|gy|gz|hx|hy|hz|ht|dab|dcd|dac1|dac2|dad1|dad2|dbc1|dbc2|dbd1|dbd2|"
    r"gam|dac|dad|dbc|dbd)(_|\()")


def _strip_dead_decls(text):
    """Drop `real(dp) ::` statements (continuations included) whose every
    name is in _CHUNK_DEAD_DECL, and the dead names from mixed ones."""
    keep, stmt = [], []
    for ln in text.split("\n"):
        stmt.append(ln)
        if ln.rstrip().endswith("&"):
            continue
        src, stmt = stmt, []
        joined = " ".join(x.strip().rstrip("&").strip() for x in src)
        m = re.match(r"real\(dp\)\s*::\s*(.*)", joined)
        if m:
            nms = [x.strip() for x in re.split(r",(?![^(]*\))", m.group(1)) if x.strip()]
            dead = [x for x in nms if _CHUNK_DEAD_DECL.match(x)]
            if dead:
                live = [x for x in nms if not _CHUNK_DEAD_DECL.match(x)]
                if live:
                    keep.append("      real(dp) :: " + ", ".join(live))
                continue
        keep.extend(src)
    return "\n".join(keep)


def chunked_rhf_item_text(tag, item, comps, inside):
    """rysgi{tag}_crhf: the CHUNKED RHF COMPACT PASS (see the section
    comment above), built from the compact item routine `item` -- its
    signature, declarations, survivor setup, same-atom skip, AO offsets,
    degeneracy weight and primitive-loop header are reused verbatim; the
    density fold, the integral build and the D accumulators are replaced."""
    la, lb, lc, ld = (int(c) for c in tag)
    nroots = (la + lb + lc + ld + 1)//2 + 1

    kfac_line = "         kfac_eff = merge(0.25_dp, 0.5_dp, ndens == 1)*kfac\n"
    i_kfac = item.index(kfac_line)
    i_k1 = item.index("            do k1 = 1, nab\n")
    tval_line = "tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tv = item.index(tval_line, i_k1) + len(tval_line)
    i_at = item.index("\n         !$acc atomic update\n")
    pre = _strip_dead_decls(item[:i_kfac])
    prim_hdr = item[i_k1:i_tv]
    tail = item[i_at:]

    names = {"cj", "ck", "ai2", "aj2", "ak2", "gyz", "gxz", "gxy"}
    chunks = _chunks(comps, inside)
    nmax = max(len(ch) for ch in chunks)
    names.update(f"gm_{n}" for n in range(nmax))
    body = [
        "         ! RHF: ndens is 1 in this routine, so kfac_eff is kfac/4\n"
        "         ! (see the item routine above for the 1/4 and the 8).\n"
        "         kfac_eff = 0.25_dp*kfac\n"
        "         cj = 8.0_dp*wq*jfac\n"
        "         ck = 8.0_dp*wq*kfac_eff\n"]
    if lb == 2 or ld == 2:
        body.append("         abx2 = abx*abx; aby2 = aby*aby; abz2 = abz*abz\n"
                    "         cdx2 = cdx*cdx; cdy2 = cdy*cdy; cdz2 = cdz*cdz\n")
        names.update(("abx2", "aby2", "abz2", "cdx2", "cdy2", "cdz2"))
    body.append("         vax = 0.0_dp; vay = 0.0_dp; vaz = 0.0_dp\n"
                "         vbx = 0.0_dp; vby = 0.0_dp; vbz = 0.0_dp\n"
                "         vcx = 0.0_dp; vcy = 0.0_dp; vcz = 0.0_dp\n")
    body.append(f"         do jch = 1, {len(chunks)}\n")
    body.append("            ! This chunk's Gamma, once, before its primitive loop.\n")
    body.append("            select case (jch)\n")
    for n, ch in enumerate(chunks):
        body.append(f"            case ({n + 1})\n")
        for m, (a, b, c, d) in enumerate(ch):
            body.append(f"               gm_{m} = {_gamma_from_dmat(a, b, c, d)}\n")
    body.append("            end select\n")
    body.append(prim_hdr)
    body.append("                  ai2 = 2.0_dp*alpha_i; aj2 = 2.0_dp*alpha_j; ak2 = 2.0_dp*alpha_k\n")
    body.append(ri.roots_inline_body(nroots, "tval", 18))
    body.append(f"""\
                  do jrt = 1, {nroots}
                     u2 = ru(jrt)/(1.0_dp - ru(jrt))*rho
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
                     select case (jch)
""")
    for n, ch in enumerate(chunks):
        body.append(f"                     case ({n + 1})\n")
        body.append(_chunk_root_body(comps, lb, ld, ch, " "*24, names))
    body.append("                     end select\n"
                "                  end do\n"
                "               end do\n"
                "            end do\n"
                "         end do\n")
    body.append("\n         ! D from translational invariance, once per quartet, in place of\n"
                "         ! a fourth set of accumulators: flx = -fix - fjx - fkx per term.\n"
                "         vdx = -(vax + vbx + vcx); vdy = -(vay + vby + vcy); vdz = -(vaz + vbz + vcz)\n")

    decls = ri.scalar_decl_lines(sorted(names), 6) + "      integer  :: jch, jrt\n"
    i_si = pre.index("si = sp_i(sOA(seg) + iab)")
    i_si = pre.rfind("\n", 0, i_si)
    pre = pre[:i_si] + "\n" + decls + pre[i_si:]
    pre = pre.replace(f"subroutine rysgi{tag}_compact(", f"subroutine rysgi{tag}_crhf(", 1)
    tail = tail.replace(f"end subroutine rysgi{tag}_compact", f"end subroutine rysgi{tag}_crhf")
    return pre + "".join(body) + tail


def screen_and_compact_text(tag, head, middle, comps=None, chunk_inside=None):
    """ORDERED SURVIVOR COMPACTION (ps|ss) prototype, generalised to every
    class: split the fused rysg{tag}/rysgi{tag} this module already built
    (in `head`/`middle`) into the two kernels below. The fused pair itself
    is only the source text they are sliced from; main() does not emit it.
    The two kernels are

      - rysg{tag}_screen / rysgi{tag}_screen_chunk: the decode and the
        three quartet-level screens (exact Schwarz, density-weighted,
        same-atom) ONLY, run RUNLEN=32 consecutive work items per `do
        concurrent` iteration (one warp's worth), appending survivors'
        (seg, iab, icd) to a dense list with ONE atomic-capture per
        chunk, sized to the chunk's own local survivor count, written in
        item order. The decode is bra-fastest, so RUNLEN consecutive
        items usually share one ket pair -- the chunked append keeps
        those survivors contiguous, instead of scattering them the way a
        per-item atomic-capture append would. RUNLEN=32 measured against
        64/128/256 on (ps|ss), w150 6-31G*, V100: 32 matches the warp
        width and wins (3.97s fused -> 3.31s per-item append -> 2.15s
        RUNLEN=32 at thresh=1e-13; 1.62 -> 0.92s at 1e-10).
      - rysg{tag}_compact / rysgi{tag}_compact: the SAME integral body
        the fused kernel runs (same_ab/same_cd/same_pair, the AO offsets,
        the density fold, Gamma, the Rys root/HRR/derivative build, the
        atomic grad updates) over that dense survivor list, by slicing
        the fused item routine's own text right after its decode and
        screens -- nothing here is retyped or re-derived by hand. The
        same-atom check is spliced in TWICE (once in the screen pass as
        `cycle`, once more where the fused text already has it before the
        AO offsets) rather than edited out of the second copy: it is
        dead the second time (the screen pass already removed every
        same-atom survivor) but removing it would mean hand-editing
        spliced text instead of just slicing around it.

    With `chunk_inside` set (a chunk_plan result; `comps` is the four
    shells' Cartesian components), the compact driver also gets the
    CHUNKED RHF COMPACT PASS, rysgi{tag}_crhf, and runs it whenever
    ndens == 1; the per-quartet rysgi{tag}_compact stays for ndens == 2.

    Every piece below is sliced out of `head`/`middle` by anchoring on
    literal text that does not depend on (la, lb, lc, ld) -- the decode,
    the screens and the primitive-pair loop headers are byte-identical
    across every class (only array sizes and the integral build differ,
    and those live entirely in the parts this function does not touch).
    """
    decl_end_marker = "real(dp), intent(inout) :: grad(3, natm)\n"
    decl_start = head.rfind(decl_end_marker) + len(decl_end_marker)
    assert decl_start > len(decl_end_marker)

    decode_start_marker = "      ! locate the segment\n"
    decode_start = head.index(decode_start_marker, decl_start)
    item_decls = head[decl_start:decode_start]

    schwarz_return = "if (qcut <= thresh) return\n"
    decode_end = head.index(schwarz_return, decode_start) + len(schwarz_return)
    decode_block = head[decode_start:decode_end]

    si_marker = "si = sp_i(sOA(seg) + iab); sj = sp_j(sOA(seg) + iab)\n"
    si_start = head.index(si_marker, decode_end)
    density_comment = "         ! Density-weighted screen."
    density_start = head.index(density_comment, si_start)
    si_only_block = head[si_start:density_start]

    density_return = "dsh(sj, sk), dsh(sj, sl)) <= thresh) return\n"
    density_end = head.index(density_return, density_start) + len(density_return)
    density_block = head[density_start:density_end]

    # same_ab/same_cd/same_pair, slotab/slotcd/nab/ncd, abx..cdz, the
    # injected same-atom+mui+wq+kfac_eff+gam block, and the primitive-pair
    # loop headers through the tval line -- unchanged, verbatim, for the
    # compact pass only (the screen pass stops at the same-atom check).
    rest_after_density = head[density_end:]

    decode_block_c = _return_to(decode_block, "cycle", 1)
    density_block_c = _return_to(density_block, "cycle", 1)
    same_atom_c = _same_atom_check("cycle")

    decl_line_old = "      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl\n"
    decl_line_new = "      integer :: p, q, mid, t, si, sj, sk, sl\n"
    compact_decls, ndc = re.subn(re.escape(decl_line_old), decl_line_new, item_decls, count=1)
    assert ndc == 1, "compact decl rewrite: expected decode-local decl line not found"

    end_marker_old = f"      end subroutine rysgi{tag}\n"
    end_marker_new = f"      end subroutine rysgi{tag}_compact\n"
    assert middle.count(end_marker_old) == 1
    compact_middle = middle.replace(end_marker_old, end_marker_new, 1)

    # LANE-PRIVATE GRADIENT COPIES: the compact pass adds into its own
    # copy gp(:, :, slot) of the gradient, not into grad itself -- see the
    # compact driver's comment below for why and what it measured.
    compact_middle, ngp = re.subn(r"grad\((\d), sh_at\((s[ijkl])\)\)",
                                  r"gp(\1, sh_at(\2), slot)", compact_middle)
    assert ngp == 24, f"compact gp rewrite: expected 24 grad references, found {ngp}"
    # Only classes of total angular momentum <= 3 spread over the copies;
    # heavier ones all add into copy 1 (the old single-grad behaviour).
    # Their compact pass is 20-460 ns per quartet against ~2 ns for
    # (ss|ss), so 12 atomics are not their bound, and the computed slot
    # alone cost them 2-7% (V100, w150 1e-13: (dp|ps) 1.33 -> 1.37 s,
    # (dp|pp) 1.55 -> 1.60 s, (dd|dd) 0.073 -> 0.078 s; back to
    # 1.35/1.57/0.073 s with slot = 1). L = 3 still gains: (ds|ps) -4%,
    # (dp|ss) -2.5%; (pp|ps) is flat.
    if sum(int(c) for c in tag) <= 3:
        slot_expr = "int(mod(i - 1, int(npriv, int64))) + 1"
    else:
        slot_expr = "1"

    screen = f"""
   !
   ! ORDERED SURVIVOR COMPACTION -- screen pass, class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}).
   ! See screen_and_compact_text's docstring (scripts/gen_rys_grad.py) for
   ! the technique and the measurement it is based on.
   !
   subroutine rysg{tag}_screen(lo, hi, w0, w1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                             out_seg, out_iab, out_icd, nmax, ctr)
      !! Screens work items [w0, w1], which must lie inside segments
      !! [lo, hi] (the decode's segment search runs over [lo, hi]); the
      !! caller windows a segment too large for one survivor list.
      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, nmax
      integer(int64), intent(in) :: w0, w1
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_at(nbas)
      integer,  intent(inout) :: out_seg(nmax), out_iab(nmax), out_icd(nmax)
      integer(int64), intent(inout) :: ctr(1)
      integer, parameter :: RUNLEN = 32
      !! Run length a chunk screens per `do concurrent` iteration: one
      !! warp's worth. Measured against 64/128/256 on (ps|ss), w150
      !! 6-31G*, V100 -- 32 matches the warp width and wins; see
      !! scripts/gen_rys_grad.py's screen_and_compact_text docstring.
      integer(int64) :: nr, nchunks, c, clen
      nr = 0
      if (w1 >= w0) nr = w1 - w0 + 1
      nchunks = (nr + RUNLEN - 1)/RUNLEN
      do concurrent(c=1:nchunks)
         call rysgi{tag}_screen_chunk(w0 + (c - 1)*RUNLEN, &
                               min(int(RUNLEN, int64), nr - (c - 1)*RUNLEN), &
                               lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                               npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                               out_seg, out_iab, out_icd, nmax, ctr)
      end do
   end subroutine rysg{tag}_screen

   pure subroutine rysgi{tag}_screen_chunk(gt0, clen, lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                             out_seg, out_iab, out_icd, nmax, ctr)
      !$acc routine seq
      integer(int64), intent(in) :: gt0, clen
      integer,  intent(in)    :: lo, hi, nseg, npair, nbas, nmax
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_at(nbas)
      integer,  intent(inout) :: out_seg(nmax), out_iab(nmax), out_icd(nmax)
      integer(int64), intent(inout) :: ctr(1)
      integer, parameter :: RUNLEN = 32
      integer :: p, q, mid, seg, t, iab, icd, si, sj, sk, sl
      integer(int64) :: nsa, u, kx, gt, slot0
      real(dp) :: qcut
      integer :: loc_seg(RUNLEN), loc_iab(RUNLEN), loc_icd(RUNLEN)
      integer :: nlocal, j

      nlocal = 0
      do j = 0, int(clen) - 1
         gt = gt0 + j

{decode_block_c}{si_only_block}{density_block_c}{same_atom_c}
         nlocal = nlocal + 1
         loc_seg(nlocal) = seg
         loc_iab(nlocal) = iab
         loc_icd(nlocal) = icd
      end do

      if (nlocal > 0) then
         !$acc atomic capture
         slot0 = ctr(1)
         ctr(1) = ctr(1) + nlocal
         !$acc end atomic
         do j = 1, nlocal
            if (slot0 + (j - 1) < int(nmax, int64)) then
               out_seg(slot0 + j) = loc_seg(j)
               out_iab(slot0 + j) = loc_iab(j)
               out_icd(slot0 + j) = loc_icd(j)
            end if
         end do
      end if
   end subroutine rysgi{tag}_screen_chunk
"""

    def compact_loop(item, indent):
        p = " "*indent
        return (f"{p}do concurrent(i=1:nout)\n"
                f"{p}   call {item}(out_seg(i), out_iab(i), out_icd(i), nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &\n"
                f"{p}                          npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &\n"
                f"{p}                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &\n"
                f"{p}                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &\n"
                f"{p}                          ndens, dmat, npriv, gp, {slot_expr})\n"
                f"{p}end do\n")

    if chunk_inside is None:
        compact_loops = compact_loop(f"rysgi{tag}_compact", 6)
    else:
        compact_loops = (
            "      ! RHF runs the CHUNKED RHF COMPACT PASS (gen_rys_grad.py); UHF keeps\n"
            "      ! the per-quartet form below.\n"
            "      if (ndens == 1) then\n"
            + compact_loop(f"rysgi{tag}_crhf", 9)
            + "      else\n"
            + compact_loop(f"rysgi{tag}_compact", 9)
            + "      end if\n")

    compact = f"""
   !
   ! ORDERED SURVIVOR COMPACTION -- compact pass, class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}).
   ! The integral body, sliced out of the fused per-work-item kernel the
   ! generator builds internally (and does not emit) to start right after
   ! the decode and the three quartet-level screens: seg/iab/icd come in
   ! directly instead of being decoded from a work-item index, and the
   ! Schwarz/density screens are not repeated (the screen pass already
   ! applied them). See scripts/gen_rys_grad.py's
   ! screen_and_compact_text docstring.
   !
   subroutine rysg{tag}_compact(nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             ndens, dmat, npriv, gp)
      integer(int64), intent(in) :: nout
      integer,  intent(in)    :: out_seg(*), out_iab(*), out_icd(*)
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac
      integer,  intent(in)    :: ao_off(nbas), sh_at(nbas), natm
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: sq_n(npair)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      integer,  intent(in)    :: npriv
      real(dp), intent(inout) :: gp(3, natm, npriv)
      !! LANE-PRIVATE GRADIENT COPIES: survivor i adds into copy
      !! mod(i-1, npriv)+1, so the 32 lanes of a warp (consecutive i) add
      !! into 32 different copies instead of all into the same grad entry.
      !! The survivor list is ket-uniform (every lane of a warp, and every
      !! warp of a segment column, shares one ket pair, and runs of bra
      !! pairs share atom A), so with a single grad(3, natm) the 12 atomics
      !! per quartet are same-address and serialise at L2: on (ss|ss) at
      !! w150 6-31G* 1e-13 they held 55% of all stall samples. grad_bins
      !! sums the copies into grad once per eri%gradient, on the device. V100,
      !! compact pass: (ss|ss) 1.34 -> 0.71 s, (ps|ss) 1.77 -> 1.20 s with
      !! 32 copies; 8 copies give 0.72/1.21 s, 64 the same as 32, and
      !! copy-major against atom-major order makes no difference. Classes
      !! with la+lb+lc+ld > 3 pass slot = 1: see gen_rys_grad.py.
      integer(int64) :: i
{compact_loops}   end subroutine rysg{tag}_compact

   pure subroutine rysgi{tag}_compact(seg, iab, icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             ndens, dmat, npriv, gp, slot)
      !$acc routine seq
      integer,  intent(in)    :: seg, iab, icd
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac
      integer,  intent(in)    :: ao_off(nbas), sh_at(nbas), natm
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: sq_n(npair)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      integer,  intent(in)    :: npriv, slot
      real(dp), intent(inout) :: gp(3, natm, npriv)
{compact_decls}
{si_only_block}{rest_after_density}{compact_middle}"""

    if chunk_inside is not None:
        i_item = compact.index(f"   pure subroutine rysgi{tag}_compact(")
        compact += "\n" + chunked_rhf_item_text(tag, compact[i_item:], comps, chunk_inside)
    return screen + compact


#: `pcsi{tag}`'s primitive-quartet head now computes the wB97X-V long-range
#: exchange substitution (rho_theta = theta*rho, theta from `omega`) -- see
#: scripts/gen_perclass.py. The gradient has no range-separated operator and
#: no `omega` dummy argument (`eri%gradient` error-stops if ever given one),
#: so the splice below removes that substitution instead of inheriting an
#: undeclared `omega` into the spliced text: `rho_theta` collapses back to
#: plain `rho`, and the two attenuation multiplies on wpx/wqx (omega == 0's
#: theta == 1.0_dp, so removing them is the identical specialisation, not an
#: approximation) go with it.
_RANGE_SEP_BLOCK = """                  ! RANGE SEPARATION. erf(omega*r12)/r12 is the Coulomb
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
"""
_RANGE_SEP_WP_WQ = ("                  wpx = theta*wpx; wpy = theta*wpy; wpz = theta*wpz\n"
                    "                  wqx = theta*wqx; wqy = theta*wqy; wqz = theta*wqz\n")


def _degauge_range_separation(text):
    """Specialise the spliced head to omega == 0 (see `_RANGE_SEP_BLOCK`)."""
    n = text.count(_RANGE_SEP_BLOCK)
    assert n == 1, n
    text = text.replace(_RANGE_SEP_BLOCK, "                  rho_theta = rho\n")
    n = text.count(_RANGE_SEP_WP_WQ)
    assert n == 1, n
    text = text.replace(_RANGE_SEP_WP_WQ, "")
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("class_tag", help="four digits, e.g. 1111 or 2222")
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="must match the --order the scalar HGP file in "
                         "--src-dir was generated with")
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
                    help="classes whose RHF compact pass uses the CHUNKED "
                         "RHF COMPACT PASS: 'auto' (chunk_auto's rule), "
                         "'none', 'all', or a comma-separated list of tags")
    ap.add_argument("--chunk-live-max", type=int, default=CHUNK_LIVE_MAX,
                    help="live-set ceiling chunk_plan picks the chunk shape under")
    ap.add_argument("-o", "--output", required=True,
                    help="never a path under src/generated: regenerated at "
                         "configure time into the build directory")
    args = ap.parse_args()
    tag = args.class_tag
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    nig, nkg = ni + 1, nk + 1       # the gradient-extended table's dims
    L = la + lb + lc + ld
    nroots = (L + 1)//2 + 1
    if nroots > 5:
        raise SystemExit(f"class {tag}: nroots={nroots} exceeds the 5 "
                          "committed root tables")
    dims = [len(co.components(l, args.order)) for l in (la, lb, lc, ld)]
    nv = dims[0]*dims[1]*dims[2]*dims[3]
    unroll = nv <= _parse_unroll_max(args.unroll_max_comp)

    src_path = os.path.join(args.src_dir, f"trc_pc_k{tag}.F90")
    src = open(src_path).read()

    i_item = src.index(f"   pure subroutine pcsi{tag}(")
    i_end = src.index(f"end subroutine pcsi{tag}\n") + len(f"end subroutine pcsi{tag}\n")
    whole = src[:i_end]
    assert i_item < i_end

    tval_line = "tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tval = whole.index(tval_line) + len(tval_line)
    head = whole[:i_tval]
    head = _degauge_range_separation(head)

    # --- 1. thread qp_a (the bra/ket primitive's own exponent) through the
    #     signatures and the one call between them, exactly like qp_z. ---
    head, n1 = re.subn(r"qp_z\(nqp\), qp_c\(nqp\)",
                       "qp_z(nqp), qp_c(nqp), qp_a(nqp)", head)
    assert n1 == 2, n1
    head, n1 = re.subn(r"nqp, qp_z, qp_c, qp_px, qp_py, qp_pz",
                       "nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz", head)
    assert n1 == 3, n1

    # --- 2. thread sh_at/natm/grad in place of jmat, next to sh_l/ao_off. ---
    head, n1 = re.subn(r"sh_l, ao_off,", "sh_l, ao_off, sh_at, natm,", head)
    assert n1 == 3, n1
    head, n1 = re.subn(r"sh_l\(nbas\), ao_off\(nbas\)",
                       "sh_l(nbas), ao_off(nbas)\n"
                       "      integer,  intent(in)    :: sh_at(nbas), natm",
                       head)
    assert n1 == 2, n1
    head, n1 = re.subn(r"dmat, jmat", "dmat, grad", head)
    assert n1 == 3, n1
    head, n1 = re.subn(r"real\(dp\), intent\(inout\) :: jmat\(ndens, nao, nao\)",
                       "real(dp), intent(inout) :: grad(3, natm)", head)
    assert n1 == 2, n1

    # --- 3. alpha_i/alpha_k: the individual primitive exponent, out of the
    #     SAME bin-interleaved slot qp_z/qp_px read -- not a new primitive
    #     loop, just one more field of the pair already being visited. ---
    head, n1 = re.subn(
        r"(rpx = qp_px\(kp\); rpy = qp_py\(kp\); rpz = qp_pz\(kp\)\n)",
        r"\1               alpha_i = qp_a(kp)\n"
        r"               alpha_j = zeta - alpha_i\n", head)
    assert n1 == 1, n1
    head, n1 = re.subn(
        r"(eta = qp_z\(kq\)\n)",
        r"\1                  alpha_k = qp_a(kq)\n", head)
    assert n1 == 1, n1

    # --- 4. drop the pre-HRR VRR declaration and the Boys scratch, exactly
    #     like gen_rys.py -- the gradient needs neither. ---
    decl_re = re.compile(
        r"      real\(dp\) :: v\((\d+), 0:1\), g1\(\1\), vbuf\((\d+)\)\n")
    m = decl_re.search(head)
    assert m, "could not find the HGP v/g1/vbuf declaration to replace"
    nv_vrr = int(m.group(1))
    nv_in_file = int(m.group(2))
    assert nv_in_file == nv, (
        f"vbuf({nv_in_file}) disagrees with na*nb*nc*nd = {nv}; "
        "the --order argument likely does not match the source file")

    decl_new = (
        f"      real(dp) :: ru({nroots}), rw({nroots})\n"
        f"      real(dp) :: u2, dinv, b00r, b10r, bp01r, f00\n"
        f"      real(dp) :: c00x, c00y, c00z, cp0x, cp0y, cp0z\n"
        f"      real(dp) :: alpha_i, alpha_j, alpha_k, kfac_eff\n"
        f"      real(dp) :: vix, viy, viz, pyz, pxz, pxy\n"
        f"      real(dp) :: fix, fiy, fiz, fjx, fjy, fjz\n"
        f"      real(dp) :: fkx, fky, fkz, flx, fly, flz\n"
        f"      real(dp) :: vax, vay, vaz, vbx, vby, vbz\n"
        f"      real(dp) :: vcx, vcy, vcz, vdx, vdy, vdz\n"
        + ri.roots_decls())
    if unroll:
        decl_new += (ri.unrolled_table_decls(nig, nkg, la + 1, lb, lc + 1, ld, 6)
                     + unrolled_gam_decls(*dims, 6))
    else:
        decl_new += (
            f"      real(dp) :: gx(0:{nig}, 0:{nkg}), gy(0:{nig}, 0:{nkg}), gz(0:{nig}, 0:{nkg})\n"
            f"      real(dp) :: hx(0:{la + 1}, 0:{lb}, 0:{lc + 1}, 0:{ld})\n"
            f"      real(dp) :: hy(0:{la + 1}, 0:{lb}, 0:{lc + 1}, 0:{ld})\n"
            f"      real(dp) :: hz(0:{la + 1}, 0:{lb}, 0:{lc + 1}, 0:{ld})\n"
            f"      integer  :: rr, aa, bb, cc, dd2\n"
            f"      integer  :: eix, eiy, eiz, ejx, ejy, ejz\n"
            f"      integer  :: ekx, eky, ekz, elx, ely, elz\n"
            f"      real(dp) :: dab({dims[0]*dims[1]}), dcd({dims[2]*dims[3]})\n"
            f"      real(dp) :: dac(2, {dims[0]*dims[2]}), dad(2, {dims[0]*dims[3]})\n"
            f"      real(dp) :: dbc(2, {dims[1]*dims[2]}), dbd(2, {dims[1]*dims[3]})\n"
            f"      real(dp) :: gam({nv})\n"
            + ri.rys2d_hrr_decls(la + 1, lb, ld, nkg))
    head = head[:m.start()] + decl_new + head[m.end():]

    for pat in (r"      real\(dp\) :: f\(0:BOYS_MMAX\)\n",
                r"      integer  :: bi, bj, bbase\n",
                r"      real\(dp\) :: bx, bx2, b0, b1, b2, btt, bet\n"):
        head, n1 = re.subn(pat, "", head, count=1)
        assert n1 == 1, pat

    # The six Fock-digestion block declarations (jab/jcd/kac/kad/kbc/kbd and
    # dab/dcd/dac/dad/dbc/dbd) are dead for the gradient -- `dab`/`dcd` are
    # redeclared above with the shape this file needs, and the rest are
    # never read at all once `tail` is replaced. Drop all four lines.
    for pat in (r"      real\(dp\) :: jab\(\d+\), jcd\(\d+\), kac\(\d+\)\n",
                r"      real\(dp\) :: kad\(\d+\), kbc\(\d+\), kbd\(\d+\)\n",
                r"      real\(dp\) :: dab\(\d+\), dcd\(\d+\), dac\(\d+\)\n",
                r"      real\(dp\) :: dad\(\d+\), dbc\(\d+\), dbd\(\d+\)\n"):
        head, n1 = re.subn(pat, "", head, count=1)
        assert n1 == 1, pat

    # No Boys function; the root lookup replaces it, inlined (see
    # rys_inline.py) rather than called, same as the energy splice.
    head, n1 = re.subn(r"   use trc_boys, only:.*?BOYS_DT, BOYS_DTINV\n",
                       ri.roots_table_use(nroots),
                       head, count=1, flags=re.S)
    assert n1 == 1

    # --- 5. the zero-the-VRR-accumulator loop becomes: same-atom skip, the
    #     degeneracy weight, the AO offsets, the density blocks and the
    #     per-component Gamma, all computed ONCE per quartet before the
    #     primitive loop -- see the module docstring for why this has to
    #     move here instead of staying in a deferred digestion. ---
    zero_re = re.compile(r"            do x = 1, " + str(nv_vrr) +
                         r"\n               g1\(x\) = 0\.0_dp\n            end do\n")
    m2 = zero_re.search(head)
    assert m2, "could not find the g1 zero loop to replace"

    d0, d1, d2, d3 = dims
    inject = f"""
         ! A quartet whose four shells sit on one atom is translation-
         ! invariant on its own: dGrad/dR = 0 for every R trivially, so skip
         ! it before the primitive loop rather than compute and cancel it.
         if (sh_at(si) == sh_at(sj) .and. sh_at(si) == sh_at(sk) &
             .and. sh_at(si) == sh_at(sl)) return

         mui = ao_off(si); nuj = ao_off(sj); lamk = ao_off(sk); sigl = ao_off(sl)
         wq = 1.0_dp
         if (same_ab) wq = wq*0.5_dp
         if (same_cd) wq = wq*0.5_dp
         if (same_pair) wq = wq*0.5_dp
         ! The 8 below is NOT part of wq: terco visits one representative
         ! quartet per (ab|cd)/(cd|ab)/(ba|cd)/... 8-fold permutation class
         ! (wq already halves it where two of those labels coincide), but
         ! dE2/dR = sum Gamma (ab|cd)/dR runs over the FULL unrestricted
         ! index sum, where a generic (non-degenerate) quartet's integral
         ! appears in all 8 equivalent terms. Found empirically against the
         ! finite-difference oracle (check_rys_grad_fd.F90): without it the
         ! analytic gradient was a uniform 1/8 of the correct one on every
         ! component, for every system and direction tried.
         !
         ! kfac_eff: RHF's Gamma = D_ab D_cd - kfac/4 (D_ac D_bd + D_ad D_bc)
         ! using the TOTAL density: the 1/4 is standard (Szabo & Ostlund),
         ! checked against the FD oracle directly. UHF's K term is NOT the
         ! same coefficient reused per spin -- it is 1/2, not 1/4 or 1 --
         ! found the same way the 8 above was: Da = Db is the RHF solution
         ! with D = 2 Da, so eri%gradient(Da, Da) must equal eri%gradient
         ! (2 Da) exactly; only kfac_eff = 1/2 per spin reproduces that
         ! identity (1/4 and 1 both failed it before the FD oracle confirmed
         ! the full density was actually wrong, not just this constant).
         kfac_eff = merge(0.25_dp, 0.5_dp, ndens == 1)*kfac

"""
    if unroll:
        inject += unrolled_gam_block(d0, d1, d2, d3)
    else:
        inject += f"""\
         do ib = 0, {d1 - 1}
            do ia = 0, {d0 - 1}
               dab(1 + ia + {d0}*ib) = 0.0_dp
               do idens = 1, ndens
                  dab(1 + ia + {d0}*ib) = dab(1 + ia + {d0}*ib) + dmat(idens, mui + ia, nuj + ib)
               end do
            end do
         end do
         do id = 0, {d3 - 1}
            do ic = 0, {d2 - 1}
               dcd(1 + ic + {d2}*id) = 0.0_dp
               do idens = 1, ndens
                  dcd(1 + ic + {d2}*id) = dcd(1 + ic + {d2}*id) + dmat(idens, lamk + ic, sigl + id)
               end do
            end do
         end do
         do idens = 1, ndens
            do ic = 0, {d2 - 1}
               do ia = 0, {d0 - 1}
                  dac(idens, 1 + ia + {d0}*ic) = dmat(idens, mui + ia, lamk + ic)
               end do
            end do
            do id = 0, {d3 - 1}
               do ia = 0, {d0 - 1}
                  dad(idens, 1 + ia + {d0}*id) = dmat(idens, mui + ia, sigl + id)
               end do
            end do
            do ic = 0, {d2 - 1}
               do ib = 0, {d1 - 1}
                  dbc(idens, 1 + ib + {d1}*ic) = dmat(idens, nuj + ib, lamk + ic)
               end do
            end do
            do id = 0, {d3 - 1}
               do ib = 0, {d1 - 1}
                  dbd(idens, 1 + ib + {d1}*id) = dmat(idens, nuj + ib, sigl + id)
               end do
            end do
         end do

         idx = 0
         do dd2 = 1, {d3}
            do cc = 1, {d2}
               do bb = 1, {d1}
                  do aa = 1, {d0}
                     idx = idx + 1
                     gam(idx) = 8.0_dp*wq*jfac*dab(aa + {d0}*(bb - 1))*dcd(cc + {d2}*(dd2 - 1))
                     do idens = 1, ndens
                        gam(idx) = gam(idx) - 8.0_dp*wq*kfac_eff*( &
                           dac(idens, aa + {d0}*(cc - 1))*dbd(idens, bb + {d1}*(dd2 - 1)) &
                           + dad(idens, aa + {d0}*(dd2 - 1))*dbc(idens, bb + {d1}*(cc - 1)))
                     end do
                  end do
               end do
            end do
         end do
"""
    inject += """
         vax = 0.0_dp; vay = 0.0_dp; vaz = 0.0_dp
         vbx = 0.0_dp; vby = 0.0_dp; vbz = 0.0_dp
         vcx = 0.0_dp; vcy = 0.0_dp; vcz = 0.0_dp
         vdx = 0.0_dp; vdy = 0.0_dp; vdz = 0.0_dp
"""
    head = head[:m2.start()] + inject + head[m2.end():]

    def ren(t):
        t = t.replace(f"trc_pc_k{tag}", f"trc_pc_rysg{tag}")
        t = t.replace(f"pcsi{tag}", f"rysgi{tag}")
        t = t.replace(f"pcs{tag}", f"rysg{tag}")
        return t

    head = re.sub(r"public :: pcs" + tag + r".*\n",
                  f"public :: rysg{tag}_screen, rysg{tag}_compact\n",
                  head, count=1)
    head = ren(head)

    old_header_re = re.compile(r"\A!\n.*?#endif\n", re.S)
    assert old_header_re.search(head)
    order_val = 0 if args.order == "lexicographic" else 1
    new_header = f"""!
! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) GRADIENT kernel, Rys quadrature.
!
! GENERATED by scripts/gen_rys_grad.py {tag} -- do not edit.
!
! Spliced from the scalar HGP kernel in src/generated/trc_pc_k{tag}.F90: the
! segment decode, screens and primitive-pair loads are that file's text,
! verbatim. The density fold and the digestion are NOT HGP's -- see this
! generator's module docstring for why they could not be.
!
! TRC_CART_ORDER GUARD -- see gen_perclass.py's non-split header for the
! full explanation; same guard, carried over from the file this was spliced
! from.
!
#if TRC_CART_ORDER != {order_val}
#error "This file was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
"""
    head = old_header_re.sub(new_header, head, count=1)

    # --- the Rys gradient integral: roots, the 2D recurrence per axis, the
    #     extended horizontal transfer, then the derivative fold straight
    #     into the 12 running sums -- no output integral block is ever
    #     materialised. ---
    if unroll:
        middle = ("\n" + unrolled_gradient_middle(tag, la, lb, lc, ld, nig, nkg, nroots, dims, args.order)
                  + "               end do\n            end do\n")
    else:
        middle = f"""
{ri.roots_inline_body(nroots, "tval", 18)}
                  do rr = 1, {nroots}
                     u2 = ru(rr)/(1.0_dp - ru(rr))*rho
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

{ri.rys2d_inline_body(nig, nkg, "c00x", "cp0x", "b00r", "b10r", "bp01r", "1.0_dp", "gx", 21)}
{ri.rys2d_inline_body(nig, nkg, "c00y", "cp0y", "b00r", "b10r", "bp01r", "1.0_dp", "gy", 21)}
{ri.rys2d_inline_body(nig, nkg, "c00z", "cp0z", "b00r", "b10r", "bp01r", "f00",    "gz", 21)}
{ri.hrr_inline_body(nig, nkg, la + 1, lb, lc + 1, ld, "gx", "abx", "cdx", "hx", 21)}
{ri.hrr_inline_body(nig, nkg, la + 1, lb, lc + 1, ld, "gy", "aby", "cdy", "hy", 21)}
{ri.hrr_inline_body(nig, nkg, la + 1, lb, lc + 1, ld, "gz", "abz", "cdz", "hz", 21)}

                     idx = 0
                     do dd2 = 1, {d3}
                        do cc = 1, {d2}
                           do bb = 1, {d1}
                              do aa = 1, {d0}
                                 idx = idx + 1
                                 eix = cnx0_{tag}(aa); eiy = cny0_{tag}(aa); eiz = cnz0_{tag}(aa)
                                 ejx = cnx1_{tag}(bb); ejy = cny1_{tag}(bb); ejz = cnz1_{tag}(bb)
                                 ekx = cnx2_{tag}(cc); eky = cny2_{tag}(cc); ekz = cnz2_{tag}(cc)
                                 elx = cnx3_{tag}(dd2); ely = cny3_{tag}(dd2); elz = cnz3_{tag}(dd2)

                                 vix = hx(eix, ejx, ekx, elx); viy = hy(eiy, ejy, eky, ely)
                                 viz = hz(eiz, ejz, ekz, elz)
                                 pyz = viy*viz; pxz = vix*viz; pxy = vix*viy

                                 fix = 2.0_dp*alpha_i*hx(eix + 1, ejx, ekx, elx)
                                 if (eix > 0) fix = fix - real(eix, dp)*hx(eix - 1, ejx, ekx, elx)
                                 fkx = 2.0_dp*alpha_k*hx(eix, ejx, ekx + 1, elx)
                                 if (ekx > 0) fkx = fkx - real(ekx, dp)*hx(eix, ejx, ekx - 1, elx)
                                 fjx = 2.0_dp*alpha_j*(hx(eix + 1, ejx, ekx, elx) + abx*hx(eix, ejx, ekx, elx))
                                 if (ejx > 0) fjx = fjx - real(ejx, dp)*hx(eix, ejx - 1, ekx, elx)
                                 flx = -fix - fjx - fkx

                                 fiy = 2.0_dp*alpha_i*hy(eiy + 1, ejy, eky, ely)
                                 if (eiy > 0) fiy = fiy - real(eiy, dp)*hy(eiy - 1, ejy, eky, ely)
                                 fky = 2.0_dp*alpha_k*hy(eiy, ejy, eky + 1, ely)
                                 if (eky > 0) fky = fky - real(eky, dp)*hy(eiy, ejy, eky - 1, ely)
                                 fjy = 2.0_dp*alpha_j*(hy(eiy + 1, ejy, eky, ely) + aby*hy(eiy, ejy, eky, ely))
                                 if (ejy > 0) fjy = fjy - real(ejy, dp)*hy(eiy, ejy - 1, eky, ely)
                                 fly = -fiy - fjy - fky

                                 fiz = 2.0_dp*alpha_i*hz(eiz + 1, ejz, ekz, elz)
                                 if (eiz > 0) fiz = fiz - real(eiz, dp)*hz(eiz - 1, ejz, ekz, elz)
                                 fkz = 2.0_dp*alpha_k*hz(eiz, ejz, ekz + 1, elz)
                                 if (ekz > 0) fkz = fkz - real(ekz, dp)*hz(eiz, ejz, ekz - 1, elz)
                                 fjz = 2.0_dp*alpha_j*(hz(eiz + 1, ejz, ekz, elz) + abz*hz(eiz, ejz, ekz, elz))
                                 if (ejz > 0) fjz = fjz - real(ejz, dp)*hz(eiz, ejz - 1, ekz, elz)
                                 flz = -fiz - fjz - fkz

                                 vax = vax + fix*pyz*gam(idx); vay = vay + fiy*pxz*gam(idx); vaz = vaz + fiz*pxy*gam(idx)
                                 vbx = vbx + fjx*pyz*gam(idx); vby = vby + fjy*pxz*gam(idx); vbz = vbz + fjz*pxy*gam(idx)
                                 vcx = vcx + fkx*pyz*gam(idx); vcy = vcy + fky*pxz*gam(idx); vcz = vcz + fkz*pxy*gam(idx)
                                 vdx = vdx + flx*pyz*gam(idx); vdy = vdy + fly*pxz*gam(idx); vdz = vdz + flz*pxy*gam(idx)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
"""

    middle += f"""
         !$acc atomic update
         grad(1, sh_at(si)) = grad(1, sh_at(si)) + vax
         !$acc atomic update
         grad(2, sh_at(si)) = grad(2, sh_at(si)) + vay
         !$acc atomic update
         grad(3, sh_at(si)) = grad(3, sh_at(si)) + vaz
         !$acc atomic update
         grad(1, sh_at(sj)) = grad(1, sh_at(sj)) + vbx
         !$acc atomic update
         grad(2, sh_at(sj)) = grad(2, sh_at(sj)) + vby
         !$acc atomic update
         grad(3, sh_at(sj)) = grad(3, sh_at(sj)) + vbz
         !$acc atomic update
         grad(1, sh_at(sk)) = grad(1, sh_at(sk)) + vcx
         !$acc atomic update
         grad(2, sh_at(sk)) = grad(2, sh_at(sk)) + vcy
         !$acc atomic update
         grad(3, sh_at(sk)) = grad(3, sh_at(sk)) + vcz
         !$acc atomic update
         grad(1, sh_at(sl)) = grad(1, sh_at(sl)) + vdx
         !$acc atomic update
         grad(2, sh_at(sl)) = grad(2, sh_at(sl)) + vdy
         !$acc atomic update
         grad(3, sh_at(sl)) = grad(3, sh_at(sl)) + vdz
      end subroutine rysgi{tag}
"""

    # As in gen_rys.py: an unrolled class needs no runtime component lookup
    # table at all -- the triples come straight off `co.components` at
    # GENERATION time (see unrolled_gradient_middle).
    tabs = ""
    if not unroll:
        for cnm, ax in (("cnx", 0), ("cny", 1), ("cnz", 2)):
            for which, l in enumerate((la, lb, lc, ld)):
                comps = co.components(l, args.order)
                vals = ", ".join(str(c[ax]) for c in comps)
                tabs += f"   integer, parameter :: {cnm}{which}_{tag}({len(comps)}) = [{vals}]\n"

    # No helper subroutines follow the component tables any more -- the
    # roots lookup, 2D recurrence and HRR are inlined straight into `middle`
    # above instead of called (see scripts/rys_inline.py).
    helpers = f"""
{tabs}
contains
"""

    assert head.count("\ncontains\n") == 1
    head = head.replace("\ncontains\n", helpers, 1)

    comps = [co.components(l, args.order) for l in (la, lb, lc, ld)]
    sel = args.chunk_classes
    if sel == "auto":
        chunked = chunk_auto(tag)
    elif sel in ("none", "all"):
        chunked = sel == "all"
    else:
        chunked = tag in sel.split(",")
    inside = chunk_plan(comps, args.chunk_live_max) if chunked else None
    sc_text = screen_and_compact_text(tag, head, middle, comps, inside)

    # Emit the module prelude and the screen/compact pair only. The fused
    # rysg{tag}/rysgi{tag} in head + middle is just the text sc_text was
    # sliced from: nothing calls it, and compiling it would roughly double
    # each file's compile time for nothing.
    prelude = head[:head.index("\ncontains\n") + len("\ncontains\n")]
    out = prelude + sc_text + f"\nend module trc_pc_rysg{tag}\n"
    with open(args.output, "w") as fh:
        fh.write(out)
    print(f"  {args.output}: ({la}{lb}|{lc}{ld}), {nroots} roots, "
          f"{len(out.splitlines())} lines"
          + (f", chunked inside={inside}" if inside is not None else ""))


if __name__ == "__main__":
    main()
