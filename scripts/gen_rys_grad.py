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
    (dd|dd) needs 5 roots and (ff|ff) 7, both inside the committed n =
    1..9 tables the energy kernels already use.
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
    fraction multiplies kfac, same meaning as eri_fock's k_scale. omega
    defaults to 0 (the ordinary Coulomb operator); omega > 0 switches to
    erf(omega*r12)/r12, the same operator eri_fock's own omega builds --
    see `_RANGE_SEP_BLOCK`'s own comment for the Rys quadrature's
    (different-looking, algebraically equivalent) form of the same HGP
    oo2ze/rz/re substitution.
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
gr = _load("gr", os.path.join(HERE, "gen_rys.py"))


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
                  u2 = theta*ru({rr})/(1.0_dp - theta*ru({rr}))*rho
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


def _gamma_from_dmat(a, b, c, d, uhf=False):
    """Gamma for one component, read straight from dmat. RHF:
    cj*D_ab*D_cd - ck*(D_ac*D_bd + D_ad*D_bc), with cj = 8 wq jfac and
    ck = 8 wq kfac/4 (see the item routine's comments on the 8 and the 1/4).
    UHF: the Coulomb term on Da + Db, the exchange term summed per spin,
    ck = 8 wq kfac/2."""
    D = lambda sp, o0, x, o1, y: f"dmat({sp}, {o0} + {x}, {o1} + {y})"
    if not uhf:
        return (f"cj*{D(1, 'mui', a, 'nuj', b)}*{D(1, 'lamk', c, 'sigl', d)} &\n"
                f"                    - ck*({D(1, 'mui', a, 'lamk', c)}*{D(1, 'nuj', b, 'sigl', d)} &\n"
                f"                    + {D(1, 'mui', a, 'sigl', d)}*{D(1, 'nuj', b, 'lamk', c)})")
    return (f"cj*({D(1, 'mui', a, 'nuj', b)} + {D(2, 'mui', a, 'nuj', b)}) &\n"
            f"                    *({D(1, 'lamk', c, 'sigl', d)} + {D(2, 'lamk', c, 'sigl', d)}) &\n"
            f"                    - ck*({D(1, 'mui', a, 'lamk', c)}*{D(1, 'nuj', b, 'sigl', d)} &\n"
            f"                    + {D(1, 'mui', a, 'sigl', d)}*{D(1, 'nuj', b, 'lamk', c)} &\n"
            f"                    + {D(2, 'mui', a, 'lamk', c)}*{D(2, 'nuj', b, 'sigl', d)} &\n"
            f"                    + {D(2, 'mui', a, 'sigl', d)}*{D(2, 'nuj', b, 'lamk', c)})")


def _chunk_root_body(comps, lb, ld, chunk, pad, names, cols=None):
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
        if cols is not None:
            out.append(_pair_fold(pad, v, fi, fj, fk, f"gm_{n}", cols))
            continue
        out.append(f"{pad}gyz = {v[1]}*{v[2]}*gm_{n}; gxz = {v[0]}*{v[2]}*gm_{n}; "
                   f"gxy = {v[0]}*{v[1]}*gm_{n}\n")
        for axis, (ax, gg) in enumerate((("x", "gyz"), ("y", "gxz"), ("z", "gxy"))):
            out.append(f"{pad}va{ax} = va{ax} + ({fi[axis]})*{gg}\n")
            out.append(f"{pad}vb{ax} = vb{ax} + ({fj[axis]})*{gg}\n")
            out.append(f"{pad}vc{ax} = vc{ax} + ({fk[axis]})*{gg}\n")
    return "".join(out)


#: Scalars every pair-table item routine of the chunked forms adds to its
#: declarations on top of the per-column sums (see _pair_fold).
_PAIR_FOLD_NAMES = ("z_yz", "z_xz", "z_xy", "z_g1", "z_g2", "z_g3",
                    "z_tix", "z_tiy", "z_tiz", "z_tjx", "z_tjy", "z_tjz",
                    "z_tkx", "z_tky", "z_tkz")


def _pair_fold(pad, v, fi, fj, fk, gm, cols):
    """One component's derivative fold for the pair-table passes: the factors
    once, then each column's nine sums times its own Gamma `gm`_c."""
    out = [f"{pad}z_yz = {v[1]}*{v[2]}; z_xz = {v[0]}*{v[2]}; z_xy = {v[0]}*{v[1]}\n"]
    for axis, ax in enumerate("xyz"):
        out.append(f"{pad}z_ti{ax} = {fi[axis]}; z_tj{ax} = {fj[axis]}; z_tk{ax} = {fk[axis]}\n")
    for c in cols:
        out.append(f"{pad}z_g1 = z_yz*{gm}_{c}; z_g2 = z_xz*{gm}_{c}; z_g3 = z_xy*{gm}_{c}\n")
        for ax, gg in (("x", "z_g1"), ("y", "z_g2"), ("z", "z_g3")):
            out.append(f"{pad}va{ax}_{c} = va{ax}_{c} + z_ti{ax}*{gg}; "
                       f"vb{ax}_{c} = vb{ax}_{c} + z_tj{ax}*{gg}; "
                       f"vc{ax}_{c} = vc{ax}_{c} + z_tk{ax}*{gg}\n")
    return "".join(out)


def _pair_gamma_block(pad, comps, cols):
    """The Gamma of every component in `comps` for every column, as ONE loop
    over the column's pairs per column. `comps` is a list of (gm, dl, dr):
    the scalar's name stem and the density-element functions `dl(p, q)` /
    `dr(p, q)` of that component (shell p against shell q, left/right pool
    density of the pair).

    Each sum gm_c still accumulates in pair order with the same statement as
    before; only the loop around it is shared. One loop per component and
    column (115 of them in one (21|20) kernel) took nvfortran 2.2-3x longer
    to compile than the fused form, for the same arithmetic."""
    out = "".join(f"{pad}{gm}_{c} = 0.0_dp\n" for gm, _, _ in comps for c in cols)
    for c in cols:
        out += (f"{pad}do ip = plo_{c}, phi_{c}\n"
                f"{pad}   jl = pl(ip); jr = pr(ip)\n")
        for gm, dl, dr in comps:
            out += (f"{pad}   {gm}_{c} = {gm}_{c} + cw*(pj(ip)*({dl(0, 1)}*{dr(2, 3)} + {dr(0, 1)}*{dl(2, 3)}) &\n"
                    f"{pad}      - pk(ip)*({dl(0, 2)}*{dr(1, 3)} + {dr(0, 2)}*{dl(1, 3)} &\n"
                    f"{pad}      + {dl(0, 3)}*{dr(1, 2)} + {dr(0, 3)}*{dl(1, 2)}))\n")
        out += f"{pad}end do\n"
    return out


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


def chunked_rhf_item_text(tag, item, comps, inside, uhf=False):
    """rysgi{tag}_crhf: the CHUNKED RHF COMPACT PASS (see the section
    comment above), built from the compact item routine `item` -- its
    signature, declarations, survivor setup, same-atom skip, AO offsets,
    degeneracy weight and primitive-loop header are reused verbatim; the
    density fold, the integral build and the D accumulators are replaced."""
    la, lb, lc, ld = (int(c) for c in tag)
    nroots = (la + lb + lc + ld + 1)//2 + 1
    rname = f"rysgi{tag}_cuhf" if uhf else f"rysgi{tag}_crhf"

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
    kfac_text = ("         ! UHF: ndens is 2 in this routine, so kfac_eff is kfac/2 per spin\n"
                 "         ! (see the item routine above for the 1/2 and the 8).\n"
                 "         kfac_eff = 0.5_dp*kfac\n") if uhf else (
                 "         ! RHF: ndens is 1 in this routine, so kfac_eff is kfac/4\n"
                 "         ! (see the item routine above for the 1/4 and the 8).\n"
                 "         kfac_eff = 0.25_dp*kfac\n")
    body = [kfac_text +
            "         cj = 8.0_dp*wq*jfac\n"
            "         ck = 8.0_dp*wq*kfac_eff\n"]
    # Every power the horizontal transfer reads, (A-B)^p and (C-D)^p up to
    # p = lb, ld (j and l are never raised): squares at d, cubes at f.
    for pw in range(2, max(lb, ld) + 1):
        prev = {2: "{0}"}.get(pw, "{0}" + str(pw - 1))
        body.append("         " + "; ".join(
            f"{v}{pw} = {prev.format(v)}*{v}" for v in ("abx", "aby", "abz")) + "\n")
        body.append("         " + "; ".join(
            f"{v}{pw} = {prev.format(v)}*{v}" for v in ("cdx", "cdy", "cdz")) + "\n")
        names.update(f"{v}{pw}" for v in ("abx", "aby", "abz", "cdx", "cdy", "cdz"))
    body.append("         vax = 0.0_dp; vay = 0.0_dp; vaz = 0.0_dp\n"
                "         vbx = 0.0_dp; vby = 0.0_dp; vbz = 0.0_dp\n"
                "         vcx = 0.0_dp; vcy = 0.0_dp; vcz = 0.0_dp\n")
    body.append(f"         do jch = 1, {len(chunks)}\n")
    body.append("            ! This chunk's Gamma, once, before its primitive loop.\n")
    body.append("            select case (jch)\n")
    for n, ch in enumerate(chunks):
        body.append(f"            case ({n + 1})\n")
        for m, (a, b, c, d) in enumerate(ch):
            body.append(f"               gm_{m} = {_gamma_from_dmat(a, b, c, d, uhf)}\n")
    body.append("            end select\n")
    body.append(prim_hdr)
    body.append("                  ai2 = 2.0_dp*alpha_i; aj2 = 2.0_dp*alpha_j; ak2 = 2.0_dp*alpha_k\n")
    body.append(ri.roots_inline_body(nroots, "tval", 18))
    body.append(f"""\
                  do jrt = 1, {nroots}
                     u2 = theta*ru(jrt)/(1.0_dp - theta*ru(jrt))*rho
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
    pre = pre.replace(f"subroutine rysgi{tag}_compact(", f"subroutine {rname}(", 1)
    tail = tail.replace(f"end subroutine rysgi{tag}_compact", f"end subroutine {rname}")
    return pre + "".join(body) + tail


# ---------------------------------------------------------------------------
# RUNTIME-OUTER CHUNKED RHF PASS
# ---------------------------------------------------------------------------
#
# The chunked pass above keeps every chunk's indices literal (a `select
# case` over chunks), so its source grows with the class: ~10 lines of fold
# per component, which is why CHUNK_MAX_COMP stops it at 324 components and
# every heavier class, (dd|dp) and up and nearly every f class, kept the
# array form and its local-memory traffic. This is gen_rys.py's CHUNKED
# FOCK PASS design carried over to the gradient: the OUTER shells' component
# comes in at run time (jo0..jo3, the cn tables), and the horizontal
# transfer over an outer shell uses run-time weights wb_s = C(eb, s)
# ab^(eb-s) against literal 2D-table indices, so one chunk's source serves
# every chunk. See gen_rys.py's section comment for the design and its
# measurements; what the gradient adds:
#
#   - Per axis and component it needs h at (a+1, b), (a-1, b) and (a, b-1)
#     on the bra and (c+1, d), (c-1, d) on the ket (l by translational
#     invariance). A shift of a or c is a LITERAL offset into the 2D table
#     even when that shell is outer, h(ea+1, eb) = sum_s w_s g(ea+1+s), so
#     it costs no new weights; the factor ea in front of h(ea-1) is the run-
#     time exponent, and the g(-1) term it would read only ever carries ea =
#     0. Only b-1 with b outer needs a second weight set, wbm_s = C(eb-1, s)
#     ab^(eb-1-s) (and wbbm with a outer too).
#   - Gamma is per chunk, read from dmat at the chunk's outer components
#     before the primitive loop, as the literal pass does.
#   - 'loop' runs every chunk in one thread per quartet; 'thread' gives each
#     (quartet, chunk) its own `do concurrent` iteration and atomics, as in
#     gen_rys.py.
#
#: The measured shape per class, (outer shells, form); a class not listed
#: keeps the literal chunked pass (chunk_auto) or the array form. Every
#: canonical class with an f shell or more than 81 components was raced
#: over all nine outer sets in both forms (none of more than 360 components
#: per chunk: (ff|ff) in one chunk is 137K lines and did not compile in
#: 2.5 h), the literal pass and the array form. w20 (the first 60 atoms of
#: w172), cc-pVTZ, RHF, 1e-13, V100; ms per class over one gradient: the
#: winner, the runner-up, and the literal pass (or the array form where the
#: literal pass is too big). The literal pass won (dd|ds), (dd|ps), (dp|ds),
#: (dp|dp), (fs|dp), (fp|ds), (fp|fs) and (ff|ps); the array form won none.
#: Over the 42 classes: 1368 ms, against 2947 ms for the literal pass where
#: it compiles and the array form elsewhere.
#:
#: Checked on a second system, NA-040 (89 atoms, C34H40N8O7), def2-TZVP,
#: raced the same way (17 shapes; not the array form or the all-outer sets,
#: which lost by 2-9x on w20): this plan takes 56.40 s over the 42 classes,
#: the per-class best over both races 56.26 s. Only (dd|dd) and (fp|dd)
#: would move (23:thread -> 23:loop, -0.25% there, +2.6 ms on w20), so the
#: plan is left as measured on w20.
#: Classes where the LITERAL pass won the race but its source is large
#: enough that a literal UHF twin doubles the object's compile time: they
#: keep the literal pass for RHF and take a run-time chunked UHF routine
#: (their best run-time shape over both races) instead. With a literal UHF
#: routine, nvfortran 26.5 took 639 s for (fp|fs), 538 s for (dp|dp) and
#: 512 s for (ff|ps), past the ~590 s per-object budget; moving the RHF side
#: too would cost NA-040/def2-TZVP 1.0 s (literal 2011/1697/662 ms against
#: run-time 2616/2052/745 ms).
GRAD_UHF_RT_PLAN = {
    "3130": ((0, 1), "loop"),
    "2121": ((2, 3), "loop"),
    "3310": ((1,), "loop"),
}

GRAD_CHUNK_PLAN = {
    "2111": ((2, 3), "loop"),  # (21|11) 86; next literal 91; literal 91
    "2211": ((2, 3), "loop"),  # (22|11) 44; next 13:loop 49; literal 49
    "2221": ((2, 3), "loop"),  # (22|21) 80; next 23:thread 82; array 308
    "2222": ((2, 3), "thread"),  # (22|22) 20; next 13:thread 21; array 62
    "3000": ((), "loop"),  # (30|00) 12; next :thread 12; literal 14
    "3010": ((), "thread"),  # (30|10) 32; next :loop 32; literal 37
    "3011": ((3,), "loop"),  # (30|11) 19; next 3:thread 23; literal 23
    "3020": ((), "thread"),  # (30|20) 23; next :loop 24; literal 30
    "3022": ((3,), "loop"),  # (30|22) 20; next 3:thread 20; array 66
    "3030": ((), "loop"),  # (30|30) 6; next literal 6; literal 6
    "3100": ((), "loop"),  # (31|00) 19; next :thread 19; literal 21
    "3110": ((1,), "loop"),  # (31|10) 40; next literal 42; literal 42
    "3111": ((2, 3), "loop"),  # (31|11) 35; next 23:thread 38; literal 43
    "3121": ((2, 3), "loop"),  # (31|21) 72; next 23:thread 75; array 300
    "3122": ((2, 3), "thread"),  # (31|22) 35; next 23:loop 37; array 113
    "3131": ((2, 3), "thread"),  # (31|31) 15; next 01:thread 19; array 70
    "3200": ((), "loop"),  # (32|00) 13; next :thread 13; literal 16
    "3210": ((1,), "loop"),  # (32|10) 33; next literal 36; literal 36
    "3211": ((1, 3), "thread"),  # (32|11) 40; next 13:loop 41; array 170
    "3220": ((1,), "loop"),  # (32|20) 52; next 1:thread 56; array 151
    "3221": ((2, 3), "thread"),  # (32|21) 76; next 01:loop 81; array 211
    "3222": ((0, 1), "thread"),  # (32|22) 37; next 23:thread 39; array 102
    "3230": ((0, 1), "loop"),  # (32|30) 28; next 1:thread 31; array 58
    "3231": ((0, 1), "thread"),  # (32|31) 40; next 23:thread 40; array 114
    "3232": ((2, 3), "thread"),  # (32|32) 15; next 23:loop 17; array 39
    "3300": ((1,), "loop"),  # (33|00) 5; next literal 5; literal 5
    "3311": ((1, 3), "thread"),  # (33|11) 14; next 13:loop 21; array 64
    "3320": ((1,), "loop"),  # (33|20) 17; next 1:thread 20; array 49
    "3321": ((0, 1), "thread"),  # (33|21) 33; next 13:loop 40; array 89
    "3322": ((0, 1), "thread"),  # (33|22) 15; next 01:loop 16; array 39
    "3330": ((0, 1), "thread"),  # (33|30) 14; next 1:thread 14; array 35
    "3331": ((0, 1), "thread"),  # (33|31) 14; next 01:loop 20; array 62
    "3332": ((0, 1), "thread"),  # (33|32) 16; next 23:loop 19; array 39
    "3333": ((0, 1), "thread"),  # (33|33) 3; next 23:thread 3; array 48
}


def _rt_shell_exp(comps, comp, s, axis):
    """The exponent of shell s on `axis` for an inner component: a literal
    int when s is inside, the run-time name ex{s}{ax} when it is outer."""
    if comp[s] == "o":
        return f"ex{s}{'xyz'[axis]}"
    return comps[s][comp[s]][axis]


def _rt_weight_lines(tag, outer, pad):
    """The run-time transfer weights of one chunk, per axis, outside the
    primitive loop: gen_rys.py's wb/wbb and wd/wdd, plus wbm (and wbbm)
    for the b-1 term when b is outer."""
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    out, names = [], []

    def wset(first, second, lsec, ltot, w, W, ab, minus):
        for ax in "xyz":
            n = lsec - minus
            out.append(f"{pad}select case (ex{second}{ax})")
            for e in range(lsec + 1):
                out.append(f"{pad}case ({e})")
                ee = e - minus
                for t in range(n + 1):
                    if ee < 0 or t > ee:
                        val = "0.0_dp"
                    else:
                        coef, p = ri.binom(ee, t), ee - t
                        f = [] if coef == 1 else [f"{float(coef)}_dp"]
                        if p == 1:
                            f.append(f"{ab}{ax}")
                        elif p > 1:
                            f.append(f"{ab}{ax}{p}")
                        val = "*".join(f) if f else "1.0_dp"
                    out.append(f"{pad}   {w}{ax}_{t} = {val}")
            out.append(f"{pad}end select")
            names.extend(f"{w}{ax}_{t}" for t in range(n + 1))
            if first in outer:
                for i in range(ltot - minus + 1):
                    out.append(f"{pad}{W}{ax}_{i} = 0.0_dp")
                for t in range(n + 1):
                    for i in range(t, ltot - minus + 1):
                        out.append(f"{pad}if (ex{first}{ax} == {i - t}) {W}{ax}_{i} = {W}{ax}_{i} + {w}{ax}_{t}")
                names.extend(f"{W}{ax}_{i}" for i in range(ltot - minus + 1))

    if 1 in outer:
        wset(0, 1, lb, ni, "wb", "wbb", "ab", 0)
        if lb > 0:
            wset(0, 1, lb, ni, "wbm", "wbbm", "ab", 1)
    if 3 in outer:
        wset(2, 3, ld, nk, "wd", "wdd", "cd", 0)
    return "\n".join(out) + ("\n" if out else ""), names


def _rt_chunk_axis(tag, comps, outer, inner, axis, pad):
    """One axis of a chunk, per root: the 2D-table corner, then the bra and
    ket transfers to every h the derivative fold reads. Returns (text,
    names, H) with H(m, sa, jm, sc) the scalar holding component m's h at
    a shifted by sa, b lowered by jm and c shifted by sc."""
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    ax = "xyz"[axis]
    bo = (0 in outer, 1 in outer)
    ko = (2 in outer, 3 in outer)
    g = lambda i, k: ri.nm(f"qg{ax}", i, k)
    lit = lambda v: isinstance(v, int)

    def key(m):
        e = [_rt_shell_exp(comps, inner[m], s, axis) for s in range(4)]
        bk = tuple(v for s, v in enumerate(e[:2]) if lit(v))
        kk = tuple(v for s, v in enumerate(e[2:]) if lit(v))
        return bk, kk, e

    # Which (bk, sa, jm) x (kk, sc) values each component needs.
    need = {}
    for m in range(len(inner)):
        bk, kk, e = key(m)
        ea, eb, ec = e[0], e[1], e[2]
        variants = [(0, 0, 0), (1, 0, 0), (0, 0, 1)]
        if not lit(ea) or ea > 0:
            variants.append((-1, 0, 0))
        if not lit(eb) or eb > 0:
            variants.append((0, 1, 0))
        if not lit(ec) or ec > 0:
            variants.append((0, 0, -1))
        need[m] = (bk, kk, variants)

    # Bra transfer terms: list of (weight text or None, g i index) for (bk, sa, jm).
    def bra_terms(bk, sa, jm):
        if not bo[1]:                      # b inside (so is a, by the outer sets)
            ea, eb = bk
            ebm = eb - jm
            res = []
            for t in range(ebm + 1):
                coef, p = ri.binom(ebm, t), ebm - t
                f = [] if coef == 1 else [f"{float(coef)}_dp"]
                if p == 1:
                    f.append(f"ab{ax}")
                elif p > 1:
                    f.append(f"ab{ax}{p}")
                res.append(("*".join(f) if f else None, ea + sa + t))
            return res
        w = "wbm" if jm else "wb"
        if not bo[0]:                      # b outer, a inside
            (ea,) = bk
            return [(f"{w}{ax}_{t}", ea + sa + t) for t in range(lb - jm + 1)
                    if ea + sa + t >= 0]
        W = "wbbm" if jm else "wbb"
        return [(f"{W}{ax}_{i}", i + sa) for i in range(ni - jm + 1) if i + sa >= 0]

    def ket_terms(kk, sc):
        if not ko[1]:
            ec, ed = kk
            res = []
            for u in range(ed + 1):
                coef, p = ri.binom(ed, u), ed - u
                f = [] if coef == 1 else [f"{float(coef)}_dp"]
                if p == 1:
                    f.append(f"cd{ax}")
                elif p > 1:
                    f.append(f"cd{ax}{p}")
                res.append(("*".join(f) if f else None, ec + sc + u))
            return res
        if not ko[0]:
            (ec,) = kk
            return [(f"wd{ax}_{u}", ec + sc + u) for u in range(ld + 1) if ec + sc + u >= 0]
        return [(f"wdd{ax}_{k}", k + sc) for k in range(nk + 1) if k + sc >= 0]

    bvars = sorted({(bk, sa, jm) for (bk, kk, vs) in need.values() for (sa, jm, sc) in vs})
    kvars = {}
    for (bk, kk, vs) in need.values():
        for (sa, jm, sc) in vs:
            kvars.setdefault((bk, sa, jm), set()).add((kk, sc))
    # the k indices each bra target is needed at
    tneed = {}
    for bv, ks in kvars.items():
        tneed[bv] = sorted({k for (kk, sc) in ks for (_, k) in ket_terms(kk, sc)})
    imax = max(i for bv in bvars for (_, i) in bra_terms(*bv))
    kmax = max(k for ks in tneed.values() for k in ks)
    g00 = "f00" if ax == "z" else "1.0_dp"
    out = [ri.rys2d_scalar_assignments(imax, kmax, f"c00{ax}", f"cp0{ax}",
                                       "b00r", "b10r", "bp01r", g00, f"qg{ax}", len(pad))]
    names = list(ri.rys2d_scalar_names(imax, kmax, f"qg{ax}"))

    def tag_of(bk, sa, jm):
        return "_".join(str(x) for x in bk) + f"s{sa + 1}j{jm}"

    def wsum(terms, val):
        parts = [(val(i) if w is None else f"{w}*{val(i)}") for (w, i) in terms]
        return gr._sum_lines(parts, pad) if parts else "0.0_dp"

    tname = {}
    for bv in bvars:
        terms = bra_terms(*bv)
        for k in tneed.get(bv, []):
            if len(terms) == 1 and terms[0][0] is None:
                tname[(bv, k)] = g(terms[0][1], k)
                continue
            nmv = f"qt{ax}_{tag_of(*bv)}_{k}"
            names.append(nmv)
            tname[(bv, k)] = nmv
            out.append(f"{pad}{nmv} = {wsum(terms, lambda i: g(i, k))}")
    hname = {}
    for bv, ks in kvars.items():
        for (kk, sc) in sorted(ks):
            terms = ket_terms(kk, sc)
            if len(terms) == 1 and terms[0][0] is None:
                hname[(bv, kk, sc)] = tname[(bv, terms[0][1])]
                continue
            nmv = f"qh{ax}_{tag_of(*bv)}_k" + "_".join(str(x) for x in kk) + f"s{sc + 1}"
            names.append(nmv)
            hname[(bv, kk, sc)] = nmv
            out.append(f"{pad}{nmv} = {wsum(terms, lambda k: tname[(bv, k)])}")

    def H(m, sa, jm, sc):
        bk, kk, _ = need[m]
        return hname[((bk, sa, jm), kk, sc)]
    return "\n".join(out) + "\n", names, H, [key(m)[2] for m in range(len(inner))]


def _rt_gamma(comp, uhf=False):
    """Gamma for one inner component, read from dmat, the outer shells'
    components at their run-time offsets jo{s}. RHF: cj D_ab D_cd - ck (D_ac
    D_bd + D_ad D_bc) on the total density. UHF: the Coulomb term on Da + Db,
    the exchange term summed per spin (the per-quartet routine's Gamma)."""
    ao = ("mui", "nuj", "lamk", "sigl")
    o = lambda s: f"{ao[s]} + jo{s}" if comp[s] == "o" else f"{ao[s]} + {comp[s]}"
    D = lambda sp, p, q: f"dmat({sp}, {o(p)}, {o(q)})"
    if not uhf:
        return (f"cj*{D(1, 0, 1)}*{D(1, 2, 3)} &\n"
                f"                    - ck*({D(1, 0, 2)}*{D(1, 1, 3)} &\n"
                f"                    + {D(1, 0, 3)}*{D(1, 1, 2)})")
    return (f"cj*({D(1, 0, 1)} + {D(2, 0, 1)}) &\n"
            f"                    *({D(1, 2, 3)} + {D(2, 2, 3)}) &\n"
            f"                    - ck*({D(1, 0, 2)}*{D(1, 1, 3)} + {D(1, 0, 3)}*{D(1, 1, 2)} &\n"
            f"                    + {D(2, 0, 2)}*{D(2, 1, 3)} + {D(2, 0, 3)}*{D(2, 1, 2)})")


def _parse_rt_plan(text):
    """--rt-plan: 'none', or '<outer digits>:<loop|thread>'."""
    if text == "none":
        return None
    outer, form = text.split(":")
    assert form in ("loop", "thread"), form
    return (tuple(int(c) for c in outer), form)


def rt_chunked_rhf_item_text(tag, item, comps, outer, thread, uhf=False):
    """rysgi{tag}_crhf (or, with `uhf`, rysgi{tag}_cuhf: ndens = 2, Gamma
    from both spin densities) in the RUNTIME-OUTER form (see the section
    comment above), from the compact item routine `item`, as
    chunked_rhf_item_text builds the literal form. Returns (routine text,
    chunk count, cn tables)."""
    rname = f"rysgi{tag}_cuhf" if uhf else f"rysgi{tag}_crhf"
    la, lb, lc, ld = (int(c) for c in tag)
    nroots = (la + lb + lc + ld + 1)//2 + 1
    dims = [len(c) for c in comps]
    nch = 1
    for x in outer:
        nch *= dims[x]

    kfac_line = "         kfac_eff = merge(0.25_dp, 0.5_dp, ndens == 1)*kfac\n"
    i_kfac = item.index(kfac_line)
    i_k1 = item.index("            do k1 = 1, nab\n")
    tval_line = "tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tv = item.index(tval_line, i_k1) + len(tval_line)
    i_at = item.index("\n         !$acc atomic update\n")
    pre = _strip_dead_decls(item[:i_kfac])
    prim_hdr = item[i_k1:i_tv]
    tail = item[i_at:]

    inner, _ = gr._chunk_keys(comps, outer, 0)
    names = {"cj", "ck", "ai2", "aj2", "ak2", "gyz", "gxz", "gxy"}
    names.update(f"gm_{n}" for n in range(len(inner)))
    kfac_text = ("         ! UHF: ndens is 2 in this routine, so kfac_eff is kfac/2 per spin\n"
                 "         ! (see the item routine above for the 1/2 and the 8).\n"
                 "         kfac_eff = 0.5_dp*kfac\n") if uhf else (
                 "         ! RHF: ndens is 1 in this routine, so kfac_eff is kfac/4\n"
                 "         ! (see the item routine above for the 1/4 and the 8).\n"
                 "         kfac_eff = 0.25_dp*kfac\n")
    body = [kfac_text +
            "         cj = 8.0_dp*wq*jfac\n"
            "         ck = 8.0_dp*wq*kfac_eff\n"]
    for pw in range(2, max(lb, ld) + 1):
        prev = {2: "{0}"}.get(pw, "{0}" + str(pw - 1))
        body.append("         " + "; ".join(
            f"{v}{pw} = {prev.format(v)}*{v}" for v in ("abx", "aby", "abz")) + "\n")
        body.append("         " + "; ".join(
            f"{v}{pw} = {prev.format(v)}*{v}" for v in ("cdx", "cdy", "cdz")) + "\n")
        names.update(f"{v}{pw}" for v in ("abx", "aby", "abz", "cdx", "cdy", "cdz"))
    body.append("         vax = 0.0_dp; vay = 0.0_dp; vaz = 0.0_dp\n"
                "         vbx = 0.0_dp; vby = 0.0_dp; vbz = 0.0_dp\n"
                "         vcx = 0.0_dp; vcy = 0.0_dp; vcz = 0.0_dp\n")
    if thread:
        body.append(f"         ! RUNTIME-OUTER CHUNKED RHF PASS: chunk jch of {nch}.\n")
    else:
        body.append(f"         ! RUNTIME-OUTER CHUNKED RHF PASS: {nch} chunks.\n"
                    f"         do jch = 1, {nch}\n")
    rem = "jch - 1"
    for n, x in enumerate(outer):
        if n == len(outer) - 1:
            body.append(f"            jo{x} = {rem}\n")
        else:
            body.append(f"            jo{x} = mod({rem}, {dims[x]})\n")
            rem = f"({rem})/{dims[x]}"
        for ax in "xyz":
            body.append(f"            ex{x}{ax} = cn{ax}{x}_{tag}(jo{x} + 1)\n")
    wtxt, wn = _rt_weight_lines(tag, outer, " "*12)
    body.append(wtxt)
    names.update(wn)
    for m, comp in enumerate(inner):
        body.append(f"            gm_{m} = {_rt_gamma(comp, uhf)}\n")
    body.append(prim_hdr)
    body.append("                  ai2 = 2.0_dp*alpha_i; aj2 = 2.0_dp*alpha_j; ak2 = 2.0_dp*alpha_k\n")
    body.append(ri.roots_inline_body(nroots, "tval", 18))
    body.append(f"""\
                  do jrt = 1, {nroots}
                     u2 = theta*ru(jrt)/(1.0_dp - theta*ru(jrt))*rho
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
    Hs, Es = [], []
    for axis in range(3):
        txt, nms, H, E = _rt_chunk_axis(tag, comps, outer, inner, axis, pad)
        body.append(txt)
        names.update(nms)
        Hs.append(H)
        Es.append(E)

    def fac(e):
        return f"{float(e)}_dp" if isinstance(e, int) else f"real({e}, dp)"

    def nonzero(e):
        return not isinstance(e, int) or e > 0

    for m in range(len(inner)):
        v, fi, fj, fk = [], [], [], []
        for axis, ax in enumerate("xyz"):
            H = lambda sa, jm, sc: Hs[axis](m, sa, jm, sc)
            ea, eb, ec = Es[axis][m][0], Es[axis][m][1], Es[axis][m][2]
            v.append(H(0, 0, 0))
            e = f"ai2*{H(1, 0, 0)}"
            if nonzero(ea):
                e += f" - {fac(ea)}*{H(-1, 0, 0)}"
            fi.append(e)
            e = f"aj2*({H(1, 0, 0)} + ab{ax}*{H(0, 0, 0)})"
            if nonzero(eb):
                e += f" - {fac(eb)}*{H(0, 1, 0)}"
            fj.append(e)
            e = f"ak2*{H(0, 0, 1)}"
            if nonzero(ec):
                e += f" - {fac(ec)}*{H(0, 0, -1)}"
            fk.append(e)
        body.append(f"{pad}gyz = {v[1]}*{v[2]}*gm_{m}; gxz = {v[0]}*{v[2]}*gm_{m}; "
                    f"gxy = {v[0]}*{v[1]}*gm_{m}\n")
        for axis, (ax, gg) in enumerate((("x", "gyz"), ("y", "gxz"), ("z", "gxy"))):
            body.append(f"{pad}va{ax} = va{ax} + ({fi[axis]})*{gg}\n")
            body.append(f"{pad}vb{ax} = vb{ax} + ({fj[axis]})*{gg}\n")
            body.append(f"{pad}vc{ax} = vc{ax} + ({fk[axis]})*{gg}\n")
    body.append("                  end do\n"
                "               end do\n"
                "            end do\n")
    if not thread:
        body.append("         end do\n")
    body.append("\n         ! D from translational invariance, once per item, in place of\n"
                "         ! a fourth set of accumulators: flx = -fix - fjx - fkx per term.\n"
                "         vdx = -(vax + vbx + vcx); vdy = -(vay + vby + vcy); vdz = -(vaz + vbz + vcz)\n")

    ints = [f"jo{x}" for x in outer] + [f"ex{x}{ax}" for x in outer for ax in "xyz"]
    decls = ri.scalar_decl_lines(sorted(names), 6) + "      integer  :: jrt\n"
    if not thread:
        decls += "      integer  :: jch\n"
    if ints:
        decls += f"      integer  :: {', '.join(ints)}\n"
    i_si = pre.index("si = sp_i(sOA(seg) + iab)")
    i_si = pre.rfind("\n", 0, i_si)
    pre = pre[:i_si] + "\n" + decls + pre[i_si:]
    pre = pre.replace(f"subroutine rysgi{tag}_compact(", f"subroutine {rname}(", 1)
    if thread:
        a = f"subroutine {rname}(seg, "
        assert pre.count(a) == 1, a
        pre = pre.replace(a, f"subroutine {rname}(jch, seg, ", 1)
        a = "      !$acc routine seq\n"
        assert pre.count(a) == 1
        pre = pre.replace(a, a + "      integer, intent(in) :: jch\n", 1)
    tail = tail.replace(f"end subroutine rysgi{tag}_compact", f"end subroutine {rname}")
    tabs = ""
    for x in outer:
        for cnm, a in (("cnx", 0), ("cny", 1), ("cnz", 2)):
            vals = ", ".join(str(c[a]) for c in comps[x])
            tabs += f"   integer, parameter :: {cnm}{x}_{tag}({dims[x]}) = [{vals}]\n"
    return pre + "".join(body) + tail, nch, tabs


def screen_and_compact_text(tag, head, middle, comps=None, chunk_inside=None, rt_plan=None):
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
    With `rt_plan` set instead (a GRAD_CHUNK_PLAN entry), rysgi{tag}_crhf
    is the RUNTIME-OUTER form. Returns (text, the cn tables that form
    needs in the module's specification part).

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
    # `dscr` belongs to the Fock screen (see density_block below).
    assert item_decls.count("      real(dp) :: qcut, bnd, dscr\n") == 1
    item_decls = item_decls.replace("      real(dp) :: qcut, bnd, dscr\n", "      real(dp) :: qcut, bnd\n")

    schwarz_return = "if (qcut <= thresh) return\n"
    decode_end = head.index(schwarz_return, decode_start) + len(schwarz_return)
    decode_block = head[decode_start:decode_end]

    si_marker = "si = sp_i(sOA(seg) + iab); sj = sp_j(sOA(seg) + iab)\n"
    si_start = head.index(si_marker, decode_end)
    density_comment = "         ! Density-weighted screen."
    density_start = head.index(density_comment, si_start)
    si_only_block = head[si_start:density_start]

    density_return = "if (qcut*dscr <= thresh) return\n"
    density_end = head.index(density_return, density_start) + len(density_return)
    # The Fock kernel's screen drops the blocks of a zero factor (jfac,
    # kfac); a gradient has no such factors and contracts all six blocks,
    # so it keeps the six-block test instead of the spliced one.
    density_block = """         ! Density-weighted screen.  What actually enters the Fock matrix is
         ! Q_ab Q_cd times a density block, so a quartet whose integrals are
         ! above threshold still contributes nothing if every density block it
         ! multiplies is negligible.  The J terms carry a factor 4 from the
         ! Coulomb degeneracy, the K terms 1.  For a converged density on an
         ! extended molecule this is where most of the screening comes from --
         ! with a flat model density it correctly fires on nothing.
         if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                      dsh(si, sk), dsh(si, sl), &
                      dsh(sj, sk), dsh(sj, sl)) <= thresh) return
"""

    # same_ab/same_cd/same_pair, slotab/slotcd/nab/ncd, abx..cdz, the
    # injected same-atom+mui+wq+kfac_eff+gam block, and the primitive-pair
    # loop headers through the tval line -- unchanged, verbatim, for the
    # compact pass only (the screen pass stops at the same-atom check).
    rest_after_density = head[density_end:]

    # The primitive-slot expressions must index on the BIN stride sSA/sSB,
    # not on the segment's extent sA/sNB: see the sSA/sSB declaration below
    # for what the latter costs on a sub-range segment. Asserted rather than
    # assumed, because this text is gen_perclass.py's.
    assert "(k1 - 1)*sSA(seg)" in rest_after_density
    assert "(k2 - 1)*sSB(seg)" in rest_after_density
    assert "(k1 - 1)*sA(seg)" not in rest_after_density
    assert "(k2 - 1)*sNB(seg)" not in rest_after_density
    # The screen pass reaches the pairs through sOA/sOB alone and never the
    # interleaved primitive tables, so a sub-range costs it nothing.
    assert "sSA" not in decode_block and "sSB" not in decode_block

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

    def compact_loop(item, indent, nch=0):
        p = " "*indent
        if nch:
            # One iteration per (survivor, chunk), survivor fastest, so a
            # warp's lanes share one chunk (gen_rys.py's 'thread' form).
            return (f"{p}do concurrent(i=1:nout*{nch})\n"
                    f"{p}   call {item}(int((i - 1)/nout) + 1, out_seg(mod(i - 1, nout) + 1), &\n"
                    f"{p}                          out_iab(mod(i - 1, nout) + 1), out_icd(mod(i - 1, nout) + 1), &\n"
                    f"{p}                          nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &\n"
                    f"{p}                          npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, omega, nbas, nao, ao_off, sh_at, natm, &\n"
                    f"{p}                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &\n"
                    f"{p}                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &\n"
                    f"{p}                          ndens, dmat, npriv, gp, {slot_expr})\n"
                    f"{p}end do\n")
        return (f"{p}do concurrent(i=1:nout)\n"
                f"{p}   call {item}(out_seg(i), out_iab(i), out_icd(i), nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &\n"
                f"{p}                          npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, omega, nbas, nao, ao_off, sh_at, natm, &\n"
                f"{p}                          nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &\n"
                f"{p}                          sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &\n"
                f"{p}                          ndens, dmat, npriv, gp, {slot_expr})\n"
                f"{p}end do\n")

    rt_nch, rt_thread = 0, False
    if rt_plan is not None:
        rt_outer, rt_form = rt_plan
        rt_thread = rt_form == "thread"
        rt_nch = 1
        for x in rt_outer:
            rt_nch *= len(comps[x])
    if chunk_inside is None and rt_plan is None:
        compact_loops = compact_loop(f"rysgi{tag}_compact", 6)
    elif rt_plan is not None:
        compact_loops = (
            "      ! RUNTIME-OUTER CHUNKED PASS (gen_rys_grad.py): crhf for RHF, cuhf\n"
            "      ! for UHF (ndens = 2).\n"
            "      if (ndens == 1) then\n"
            + compact_loop(f"rysgi{tag}_crhf", 9, rt_nch if rt_thread else 0)
            + "      else\n"
            + compact_loop(f"rysgi{tag}_cuhf", 9, rt_nch if rt_thread else 0)
            + "      end if\n")
    else:
        compact_loops = (
            "      ! CHUNKED COMPACT PASS (gen_rys_grad.py): crhf for RHF, cuhf for\n"
            "      ! UHF (ndens = 2).\n"
            "      if (ndens == 1) then\n"
            + compact_loop(f"rysgi{tag}_crhf", 9)
            + "      else\n"
            + compact_loop(f"rysgi{tag}_cuhf", 9)
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
   subroutine rysg{tag}_compact(nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &
                             npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, omega, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             ndens, dmat, npriv, gp)
      integer(int64), intent(in) :: nout
      integer,  intent(in)    :: out_seg(*), out_iab(*), out_icd(*)
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      integer,  intent(in)    :: sSA(nseg), sSB(nseg)
      !! Pairs in the whole bra / ket BIN: the stride of the bin-interleaved
      !! primitive tables, which sA/sNB are NOT once a segment covers only
      !! part of its bin. Every segment of a skeleton gradient plan
      !! (bins_plan_build_skeleton) is such a sub-range, and so is an
      !! exchange tile. Indexing on sA/sNB there reads another pair's
      !! primitives: measured, a symmetrized skeleton gradient wrong by
      !! 16.6 (water), 3.4 (B2H6), 2.3e+05 (SF6) and 1.8e+03 (benzene)
      !! against components of order one.
      logical,  intent(in)    :: sD(nseg)
      real(dp), intent(in)    :: sW(nseg)
      !! The segment's orbit weight, folded into the degeneracy factor `wq`
      !! below. One for an ordinary plan -- grad_bins hands an array of ones
      !! rather than branch in here, as binkernel_fock_bins does -- and a
      !! multiple of one half for a skeleton plan, where it stands in for the
      !! images of the bra representatives the plan dropped.
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac, omega
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

   pure subroutine rysgi{tag}_compact(seg, iab, icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &
                             npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, omega, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             ndens, dmat, npriv, gp, slot)
      !$acc routine seq
      integer,  intent(in)    :: seg, iab, icd
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      integer,  intent(in)    :: sSA(nseg), sSB(nseg)
      !! Pairs in the whole bra / ket BIN: the stride of the bin-interleaved
      !! primitive tables, which sA/sNB are NOT once a segment covers only
      !! part of its bin. Every segment of a skeleton gradient plan
      !! (bins_plan_build_skeleton) is such a sub-range, and so is an
      !! exchange tile. Indexing on sA/sNB there reads another pair's
      !! primitives: measured, a symmetrized skeleton gradient wrong by
      !! 16.6 (water), 3.4 (B2H6), 2.3e+05 (SF6) and 1.8e+03 (benzene)
      !! against components of order one.
      logical,  intent(in)    :: sD(nseg)
      real(dp), intent(in)    :: sW(nseg)
      !! The segment's orbit weight, folded into the degeneracy factor `wq`
      !! below. One for an ordinary plan -- grad_bins hands an array of ones
      !! rather than branch in here, as binkernel_fock_bins does -- and a
      !! multiple of one half for a skeleton plan, where it stands in for the
      !! images of the bra representatives the plan dropped.
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac, omega
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

    tabs = ""
    if chunk_inside is not None:
        i_item = compact.index(f"   pure subroutine rysgi{tag}_compact(")
        item = compact[i_item:]
        # As for the run-time form below: both from the compact routine
        # alone, which is then not emitted (nothing calls it any more).
        uhf_rt = GRAD_UHF_RT_PLAN.get(tag)
        if uhf_rt is None:
            txt_u = chunked_rhf_item_text(tag, item, comps, chunk_inside, uhf=True)
        else:
            txt_u, _, tabs = rt_chunked_rhf_item_text(tag, item, comps, uhf_rt[0],
                                                      uhf_rt[1] == "thread", uhf=True)
            assert uhf_rt[1] == "loop", "a thread-form UHF override needs its own launch"
        compact = (compact[:i_item] + chunked_rhf_item_text(tag, item, comps, chunk_inside)
                   + "\n" + txt_u)
    elif rt_plan is not None:
        i_item = compact.index(f"   pure subroutine rysgi{tag}_compact(")
        # Both from the compact routine alone: slicing after the first is
        # appended would carry a copy of it into the second one's tail.
        item = compact[i_item:]
        txt, _, tabs = rt_chunked_rhf_item_text(tag, item, comps, rt_plan[0], rt_thread)
        txt_u, _, _ = rt_chunked_rhf_item_text(tag, item, comps, rt_plan[0], rt_thread, uhf=True)
        # The per-quartet routine is only their source text now: RHF and UHF
        # both run a chunked pass, so it is not emitted (nothing calls it).
        compact = compact[:i_item] + txt + "\n" + txt_u
    return screen + compact, tabs


#: `pcsi{tag}`'s primitive-quartet head already computes the wB97X-V
#: long-range exchange substitution (rho_theta = theta*rho, theta from
#: `omega`) -- see scripts/gen_perclass.py, whose comment is reproduced
#: below verbatim since the splice keeps it. `omega` is already a dummy
#: argument of the spliced head (gen_perclass.py's PROLOGUE declares it
#: unconditionally), so nothing needs inheriting from this splice point
#: on down -- but what the gradient's per-root Rys coefficients (`u2`
#: onward, below `pcsi{tag}`'s tval line) need is NOT "read rho_theta
#: instead of rho" (that distorts the (zeta*eta+u2*(zeta+eta)) denominator
#: u2 sits in, measured wrong by 10-70% against a finite difference of
#: the attenuated exchange energy -- two different wrong substitutions
#: were tried and measured before this one). `tval` IS the one place
#: X = rho_theta*PQ^2 belongs (in `head`, unchanged): that argument is
#: what the root/weight TABLE is evaluated at. What ALSO needs omega,
#: checked against gpu4pyscf's own range-separated root routine
#: (gpu4pyscf/lib/gvhf-rys/rys_roots.cu's `rys_roots_rs`: roots/weights
#: looked up at the attenuated X, then the ROOT itself scaled by theta
#: and the WEIGHT by sqrt(theta), on top of that lookup) is the ROOT
#: `ru(n)` itself, wherever `u2` reads it -- `theta*ru(n)` in place of
#: `ru(n)`, with `rho` in u2's own formula left plain. The weight
#: `rw(n)` gets its sqrt(theta) through `pref` instead (linear, so
#: either side of the `f00 = pref*rw(n)` product is equivalent) -- see
#: below.
#:
#: ONE thing the splice does NOT inherit: HGP's `pref = pref*sqrt(theta)`
#: sits in its TAIL, right after the Boys evaluation this generator never
#: splices in (the Rys roots replace it) -- so it is added explicitly
#: here, in the same place relative to `rho_theta`, before `pref` is
#: first read (the per-root `f00 = pref*rw(rr)`).
#:
#: omega == 0 keeps every one of these bit-identical to the unattenuated
#: path: theta == 1.0_dp exactly, so rho_theta == rho, pref*sqrt(theta) ==
#: pref and theta*ru(n) == ru(n), each an exact multiply by one.
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


def _thread_range_separation(text):
    """Keep the spliced head's range-separation substitution (see
    `_RANGE_SEP_BLOCK`'s own comment) and add the one line HGP's tail has
    that this splice does not reach: `pref = pref*sqrt(theta)`, right
    after `rho_theta` is formed."""
    n = text.count(_RANGE_SEP_BLOCK)
    assert n == 1, n
    return text.replace(_RANGE_SEP_BLOCK,
                        _RANGE_SEP_BLOCK + "                  pref = pref*sqrt(theta)\n")


#: Columns one pair-table sweep carries. The widest pair kernel keeps
#: gam(nv, NCB) and acc(12, NCB) per thread, so this bounds its private
#: memory; the dispatcher (gen_rys_grad_dispatch.py) runs another sweep per
#: further block of NCB columns.
PAIRS_NCB = 4

#: Largest gm count (Gamma values live per chunk and column) a class may carry
#: at each width. A first guard against the steepest compile times; it does
#: not separate the slow kernels from the fast ones by itself, which is what
#: PAIRS_NARROW_ONLY below is for. Measured 2026-11-17 (nvfortran 26.5,
#: -stdpar=gpu -acc -gpu=cc70,mem:separate -fast -O3, TERCO_LMAX=3, each
#: kernel compiled alone, 16 at a time on a box with load average 10-25):
#: every runtime-outer ("rt") kernel took 20-90 s whatever its width, and
#: every literal one (a class with no rt plan) grew with the class far more
#: than with nc*gm, e.g. (21|10) 27/68/193 s at widths 1/2/4 but (30|21)
#: 190/526 s at 1/2 for a smaller nc*gm.
PAIRS_GM_BUDGET = {4: 72, 2: 80}

#: Classes whose wider-than-1 pair kernels took more than 300 s to compile
#: (seconds measured as above): they are generated at width 1 alone, and the
#: dispatcher walks their column blocks one column at a time. The pair
#: kernels compile serially per file, so a file much past five minutes holds
#: the whole build, and no caller in src/ uses more than one column yet, so
#: compile time wins over the sweep sharing.
#:
#: RUN-TIME COST, not yet measured: for these classes an ncol-column
#: gradient_pairs call recomputes every integral derivative ncol times, not
#: ceil(ncol/2) times (w2 was the widest they had) or, for (31|30), not
#: ceil(ncol/4). Only the classes in the table pay it. The one routine that
#: takes ncol > 1 is eri%gradient_pairs (eri_gradient_pairs in
#: src/eri/trc_eri.F90), meant for excited-state gradients with one column
#: per root; nothing in src/ or app/ calls it yet, only
#: test/check_grad_pairs.F90 does.
#:
#: Not fixed by this table: (21|21), (31|30) and (33|10) at width 1 still
#: take 500-600 s each to compile (width 1 must exist for the dispatcher).
#: Bringing them under five minutes needs a different kernel shape.
#:
#: Seconds, nvfortran 26.5, LMAX=3, kernels alone ("-": no such width):
#:
#:   class   w1      w2      w4
#:   2121    >600    >600    -
#:   2220    264     >600    -
#:   3021    190     526     -
#:   3120    258     517     -
#:   3130    568     >600    >600
#:   3310    509     >600    -
#:
#: Widths under 300 s that stay, for reference (w1/w2/w4): (11|11) 72/147/-,
#: (20|11) 20/57/206, (21|10) 27/68/193, (21|20) 97/222/-, (22|10) 125/261/-.
PAIRS_NARROW_ONLY = ("2121", "2220", "3021", "3120", "3130", "3310")

#: Column widths a pair kernel is generated for. The twelve running sums of
#: every column are unrolled into scalar accumulators so they stay in
#: registers (a runtime column loop forces them into local memory: measured
#: 5.9 GB of local loads against 0.14 GB for the RHF/UHF kernel on one
#: class, 4x the time), which means the width is a compile-time constant.
#: A block of ncb columns runs the narrowest variant that holds it; the
#: columns past ncb of a wider variant are skipped.
PAIRS_VARIANTS = (1, 2, PAIRS_NCB)


def _dim(loopvar, dims):
    """Extent of the density-block loop variable ia/ib/ic/id."""
    return dims["abcd".index(loopvar[1])]


def _pairs_header(t, tag, nc, need):
    """Rewrite the header of the per-quartet compact routine text `t` into a
    pair-table item routine's: no jfac/kfac, the pool and the column table
    in place of ndens and the density, and gp with a column index."""
    def sub1(text, old, new):
        need(text.count(old) == 1, f"expected exactly one {old!r}, found {text.count(old)}")
        return text.replace(old, new, 1)

    t = sub1(t, "pcut, jfac, kfac, omega, nbas,", "pcut, omega, nbas,")
    t = sub1(t, "real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac, omega\n",
             "real(dp), intent(in)    :: sp_q(npair), pcut, omega\n")
    t = sub1(t, "ndens, dmat, npriv, gp, slot)",
             "npool, dmat, ncb, pstart, pl, pr, pj, pk, npriv, gp, slot)")
    return sub1(t, """      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      integer,  intent(in)    :: npriv, slot
      real(dp), intent(inout) :: gp(3, natm, npriv)
""", """      integer,  intent(in)    :: npool, ncb
      real(dp), intent(in)    :: dmat(npool, nao, nao)
      integer,  intent(in)    :: pstart(ncb + 1), pl(*), pr(*)
      real(dp), intent(in)    :: pj(*), pk(*)
      !! The column pair table (see pairs_modules in gen_rys_grad.py).
      integer,  intent(in)    :: npriv, slot
      real(dp), intent(inout) :: gp(3, natm, npriv, ncb)
""")


def _pair_names(ngm, cols):
    """Scalar names the chunked pair passes declare: Gamma per component and
    column, the per-column sums, the fold scratch and the 8 wq."""
    names = {"cw", "ai2", "aj2", "ak2"}
    names.update(_PAIR_FOLD_NAMES)
    names.update(f"gm_{n}_{c}" for n in range(ngm) for c in cols)
    names.update(f"v{a}{ax}_{c}" for a in "abcd" for ax in "xyz" for c in cols)
    return names


def _pair_column_setup(cols):
    """Per-column pair ranges (empty past ncb) and zeroed sums."""
    out = ""
    for c in cols:
        out += (f"         plo_{c} = 1; phi_{c} = 0\n"
                f"         if (ncb >= {c}) then\n"
                f"            plo_{c} = pstart({c}); phi_{c} = pstart({c + 1}) - 1\n"
                f"         end if\n")
        for a in "abc":
            out += "         " + "; ".join(f"v{a}{ax}_{c} = 0.0_dp" for ax in "xyz") + "\n"
    return out


def _pair_column_tail(cols):
    """D from translational invariance, then each column's atomics into its
    own slice of the lane-private copies."""
    out = "\n         ! D from translational invariance, once per column.\n"
    for c in cols:
        out += "         " + "; ".join(
            f"vd{ax}_{c} = -(va{ax}_{c} + vb{ax}_{c} + vc{ax}_{c})" for ax in "xyz") + "\n"
    for c in cols:
        out += f"         if (ncb >= {c}) then\n"
        for sh, a in zip("ijkl", "abcd"):
            for k, ax in enumerate("xyz"):
                out += ("            !$acc atomic update\n"
                        f"            gp({k + 1}, sh_at(s{sh}), slot, {c}) = "
                        f"gp({k + 1}, sh_at(s{sh}), slot, {c}) + v{a}{ax}_{c}\n")
        out += "         end if\n"
    return out


def _pair_item_frame(tag, item, nc, rname, need):
    """The pieces of a chunked pair routine taken from the per-quartet
    compact routine `item`: (header, primitive-loop head, tail)."""
    kfac_line = "         kfac_eff = merge(0.25_dp, 0.5_dp, ndens == 1)*kfac\n"
    i_kfac = item.index(kfac_line)
    i_k1 = item.index("            do k1 = 1, nab\n")
    tval_line = "tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tv = item.index(tval_line, i_k1) + len(tval_line)
    i_at = item.index("\n         !$acc atomic update\n")
    pre = _strip_dead_decls(item[:i_kfac])
    pre = _pairs_header(pre[pre.index(f"   pure subroutine rysgi{tag}_compact("):], tag, nc, need)
    tail = item[i_at:]
    tail = tail[tail.index(f"end subroutine rysgi{tag}_compact"):]
    return pre, item[i_k1:i_tv], tail


def _pair_item_finish(tag, pre, decls, rname, tail, thread=False):
    i_si = pre.index("si = sp_i(sOA(seg) + iab)")
    i_si = pre.rfind("\n", 0, i_si)
    pre = pre[:i_si] + "\n" + decls + pre[i_si:]
    pre = pre.replace(f"subroutine rysgi{tag}_compact(", f"subroutine {rname}(", 1)
    if thread:
        a = f"subroutine {rname}(seg, "
        assert pre.count(a) == 1, a
        pre = pre.replace(a, f"subroutine {rname}(jch, seg, ", 1)
        a = "      !$acc routine seq\n"
        assert pre.count(a) == 1
        pre = pre.replace(a, a + "      integer, intent(in) :: jch\n", 1)
    return pre, tail.replace(f"end subroutine rysgi{tag}_compact", f"end subroutine {rname}")


def _rt_pair_dens(comp):
    ao = ("mui", "nuj", "lamk", "sigl")
    o = lambda s: f"{ao[s]} + jo{s}" if comp[s] == "o" else f"{ao[s]} + {comp[s]}"
    return (lambda p, q: f"dmat(jl, {o(p)}, {o(q)})"), (lambda p, q: f"dmat(jr, {o(p)}, {o(q)})")


def _literal_pair_dens(comp):
    ao = ("mui", "nuj", "lamk", "sigl")
    return ((lambda p, q: f"dmat(jl, {ao[p]} + {comp[p]}, {ao[q]} + {comp[q]})"),
            (lambda p, q: f"dmat(jr, {ao[p]} + {comp[p]}, {ao[q]} + {comp[q]})"))


def _pair_power_lines(lb, ld, names):
    out = ""
    for pw in range(2, max(lb, ld) + 1):
        prev = {2: "{0}"}.get(pw, "{0}" + str(pw - 1))
        out += "         " + "; ".join(f"{v}{pw} = {prev.format(v)}*{v}" for v in ("abx", "aby", "abz")) + "\n"
        out += "         " + "; ".join(f"{v}{pw} = {prev.format(v)}*{v}" for v in ("cdx", "cdy", "cdz")) + "\n"
        names.update(f"{v}{pw}" for v in ("abx", "aby", "abz", "cdx", "cdy", "cdz"))
    return out


_ROOT_LOOP_HEAD = """\
                  do jrt = 1, {nroots}
                     u2 = theta*ru(jrt)/(1.0_dp - theta*ru(jrt))*rho
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
"""


def chunked_pairs_item_text(tag, item, comps, inside, nc):
    """rysgi{tag}_pairs{nc} in the LITERAL chunked form of chunked_rhf_item_text:
    each chunk forms Gamma for its components and every column, then folds
    the shared derivative factors into nc sets of nine scalar sums."""
    def need(cond, msg):
        if not cond:
            raise SystemExit(f"gen_rys_grad.py chunked_pairs_item_text({tag}): {msg}")

    rname = f"rysgi{tag}_pairs{nc}"
    la, lb, lc, ld = (int(c) for c in tag)
    nroots = (la + lb + lc + ld + 1)//2 + 1
    cols = range(1, nc + 1)
    pre, prim_hdr, tail = _pair_item_frame(tag, item, nc, rname, need)
    chunks = _chunks(comps, inside)
    nmax = max(len(ch) for ch in chunks)
    names = _pair_names(nmax, cols)
    body = ["         cw = 8.0_dp*wq\n", _pair_power_lines(lb, ld, names), _pair_column_setup(cols)]
    body.append(f"         do jch = 1, {len(chunks)}\n")
    body.append("            ! This chunk's Gamma per column, once, before its primitive loop.\n")
    body.append("            select case (jch)\n")
    for n, ch in enumerate(chunks):
        body.append(f"            case ({n + 1})\n")
        body.append(_pair_gamma_block(
            " "*15, [("gm_%d" % m, *_literal_pair_dens(comp)) for m, comp in enumerate(ch)], cols))
    body.append("            end select\n")
    body.append(prim_hdr)
    body.append("                  ai2 = 2.0_dp*alpha_i; aj2 = 2.0_dp*alpha_j; ak2 = 2.0_dp*alpha_k\n")
    body.append(ri.roots_inline_body(nroots, "tval", 18))
    body.append(_ROOT_LOOP_HEAD.format(nroots=nroots))
    body.append("                     select case (jch)\n")
    for n, ch in enumerate(chunks):
        body.append(f"                     case ({n + 1})\n")
        body.append(_chunk_root_body(comps, lb, ld, ch, " "*24, names, cols=cols))
    body.append("                     end select\n"
                "                  end do\n"
                "               end do\n"
                "            end do\n"
                "         end do\n")
    body.append(_pair_column_tail(cols))
    ints = ["ip", "jl", "jr"] + [f"plo_{c}" for c in cols] + [f"phi_{c}" for c in cols]
    decls = (ri.scalar_decl_lines(sorted(names), 6) + "      integer  :: jch, jrt\n"
             + f"      integer  :: {', '.join(ints)}\n")
    pre, tail = _pair_item_finish(tag, pre, decls, rname, tail)
    return pre + "".join(body) + "      " + tail


def rt_pairs_item_text(tag, item, comps, outer, thread, nc):
    """rysgi{tag}_pairs{nc} in the RUNTIME-OUTER chunked form of
    rt_chunked_rhf_item_text (same integral build; see chunked_pairs_item_text
    for what differs). Column c > ncb has no pairs and is skipped in the tail."""
    def need(cond, msg):
        if not cond:
            raise SystemExit(f"gen_rys_grad.py rt_pairs_item_text({tag}): {msg}")

    rname = f"rysgi{tag}_pairs{nc}"
    la, lb, lc, ld = (int(c) for c in tag)
    nroots = (la + lb + lc + ld + 1)//2 + 1
    dims = [len(c) for c in comps]
    nch = 1
    for x in outer:
        nch *= dims[x]
    cols = range(1, nc + 1)
    pre, prim_hdr, tail = _pair_item_frame(tag, item, nc, rname, need)
    inner, _ = gr._chunk_keys(comps, outer, 0)
    names = _pair_names(len(inner), cols)
    body = ["         cw = 8.0_dp*wq\n", _pair_power_lines(lb, ld, names), _pair_column_setup(cols)]
    if thread:
        body.append(f"         ! RUNTIME-OUTER CHUNKED PAIR PASS: chunk jch of {nch}.\n")
    else:
        body.append(f"         ! RUNTIME-OUTER CHUNKED PAIR PASS: {nch} chunks.\n"
                    f"         do jch = 1, {nch}\n")
    rem = "jch - 1"
    for n, x in enumerate(outer):
        if n == len(outer) - 1:
            body.append(f"            jo{x} = {rem}\n")
        else:
            body.append(f"            jo{x} = mod({rem}, {dims[x]})\n")
            rem = f"({rem})/{dims[x]}"
        for ax in "xyz":
            body.append(f"            ex{x}{ax} = cn{ax}{x}_{tag}(jo{x} + 1)\n")
    wtxt, wn = _rt_weight_lines(tag, outer, " "*12)
    body.append(wtxt)
    names.update(wn)
    body.append(_pair_gamma_block(
        " "*12, [("gm_%d" % m, *_rt_pair_dens(comp)) for m, comp in enumerate(inner)], cols))
    body.append(prim_hdr)
    body.append("                  ai2 = 2.0_dp*alpha_i; aj2 = 2.0_dp*alpha_j; ak2 = 2.0_dp*alpha_k\n")
    body.append(ri.roots_inline_body(nroots, "tval", 18))
    body.append(_ROOT_LOOP_HEAD.format(nroots=nroots))
    pad = " "*21
    Hs, Es = [], []
    for axis in range(3):
        txt, nms, H, E = _rt_chunk_axis(tag, comps, outer, inner, axis, pad)
        body.append(txt)
        names.update(nms)
        Hs.append(H)
        Es.append(E)

    def fac(e):
        return f"{float(e)}_dp" if isinstance(e, int) else f"real({e}, dp)"

    def nonzero(e):
        return not isinstance(e, int) or e > 0

    for m in range(len(inner)):
        v, fi, fj, fk = [], [], [], []
        for axis, ax in enumerate("xyz"):
            H = lambda sa, jm, sc: Hs[axis](m, sa, jm, sc)
            ea, eb, ec = Es[axis][m][0], Es[axis][m][1], Es[axis][m][2]
            v.append(H(0, 0, 0))
            e = f"ai2*{H(1, 0, 0)}"
            if nonzero(ea):
                e += f" - {fac(ea)}*{H(-1, 0, 0)}"
            fi.append(e)
            e = f"aj2*({H(1, 0, 0)} + ab{ax}*{H(0, 0, 0)})"
            if nonzero(eb):
                e += f" - {fac(eb)}*{H(0, 1, 0)}"
            fj.append(e)
            e = f"ak2*{H(0, 0, 1)}"
            if nonzero(ec):
                e += f" - {fac(ec)}*{H(0, 0, -1)}"
            fk.append(e)
        body.append(_pair_fold(pad, v, fi, fj, fk, f"gm_{m}", cols))
    body.append("                  end do\n"
                "               end do\n"
                "            end do\n")
    if not thread:
        body.append("         end do\n")
    body.append(_pair_column_tail(cols))
    ints = ([f"jo{x}" for x in outer] + [f"ex{x}{ax}" for x in outer for ax in "xyz"]
            + ["ip", "jl", "jr"] + [f"plo_{c}" for c in cols] + [f"phi_{c}" for c in cols])
    decls = ri.scalar_decl_lines(sorted(names), 6) + "      integer  :: jrt\n"
    if not thread:
        decls += "      integer  :: jch\n"
    decls += f"      integer  :: {', '.join(ints)}\n"
    pre, tail = _pair_item_finish(tag, pre, decls, rname, tail, thread)
    return pre + "".join(body) + "      " + tail


def _pairs_variant(tag, item, dims, nc):
    """One pair item routine, rysgi{tag}_pairs{nc}, for blocks of at most nc
    columns: pairs_modules' docstring says what it computes."""
    d0, d1, d2, d3 = dims
    nv = d0*d1*d2*d3

    def need(cond, msg):
        # Not an assert: `python -O` must not strip a check on spliced text.
        if not cond:
            raise SystemExit(f"gen_rys_grad.py pairs_modules({tag}): {msg}")

    def sub1(text, old, new):
        need(text.count(old) == 1, f"expected exactly one {old!r}, found {text.count(old)}")
        return text.replace(old, new, 1)

    i_item = item.index(f"   pure subroutine rysgi{tag}_compact(")
    t = item[i_item:]
    need(t.count(f"rysgi{tag}_compact") == 2, "the compact routine's name should occur twice (header, end)")
    t = t.replace(f"rysgi{tag}_compact", f"rysgi{tag}_pairs{nc}")
    t = sub1(t, "pcut, jfac, kfac, omega, nbas,", "pcut, omega, nbas,")
    t = sub1(t, "real(dp), intent(in)    :: sp_q(npair), pcut, jfac, kfac, omega\n",
             "real(dp), intent(in)    :: sp_q(npair), pcut, omega\n")
    t = sub1(t, "ndens, dmat, npriv, gp, slot)",
             "npool, dmat, ncb, pstart, pl, pr, pj, pk, npriv, gp, slot)")
    t = sub1(t, """      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      integer,  intent(in)    :: npriv, slot
      real(dp), intent(inout) :: gp(3, natm, npriv)
""", f"""      integer,  intent(in)    :: npool, ncb
      real(dp), intent(in)    :: dmat(npool, nao, nao)
      integer,  intent(in)    :: pstart(ncb + 1), pl(*), pr(*)
      real(dp), intent(in)    :: pj(*), pk(*)
      !! The column pair table (see pairs_modules in gen_rys_grad.py).
      integer,  intent(in)    :: npriv, slot
      real(dp), intent(inout) :: gp(3, natm, npriv, ncb)
""")
    t = sub1(t, "      real(dp) :: vax, vay, vaz, vbx, vby, vbz\n"
                "      real(dp) :: vcx, vcy, vcz, vdx, vdy, vdz\n",
             f"      ! {nc} columns wide (PAIRS_VARIANTS in gen_rys_grad.py): acc stays in\n"
             f"      ! registers, gam is nv*{nc} per-thread doubles.\n"
             f"      real(dp) :: acc(12, {nc})\n"
             "      integer  :: col, ip, jl, jr\n"
             f"      real(dp) :: lab({d0*d1}), lcd({d2*d3}), lac({d0*d2}), lad({d0*d3}), lbc({d1*d2}), lbd({d1*d3})\n"
             f"      real(dp) :: rab({d0*d1}), rcd({d2*d3}), rac({d0*d2}), rad({d0*d3}), rbc({d1*d2}), rbd({d1*d3})\n")
    for dl in (f"      real(dp) :: dab({d0*d1}), dcd({d2*d3})\n",
               f"      real(dp) :: dac(2, {d0*d2}), dad(2, {d0*d3})\n",
               f"      real(dp) :: dbc(2, {d1*d2}), dbd(2, {d1*d3})\n"):
        t = sub1(t, dl, "")
    t = sub1(t, f"      real(dp) :: gam({nv})\n", f"      real(dp) :: gam({nv}, {nc})\n")
    # Gamma: from the kfac_eff line through the end of the gam fill.
    a = t.index("         kfac_eff = merge(")
    e = t.index("         vdx = 0.0_dp; vdy = 0.0_dp; vdz = 0.0_dp\n")
    e = t.index("\n", e) + 1
    # The span being replaced may hold only the density-block / Gamma fill
    # and the sum zeroing: a line a later generator edit puts here would
    # otherwise be dropped silently.
    dropped = {"kfac_eff", "dab", "dcd", "dac", "dad", "dbc", "dbd", "gam",
               "vax", "vay", "vaz", "vbx", "vby", "vbz", "vcx", "vcy", "vcz", "vdx", "vdy", "vdz",
               "idx"}
    cont = False
    for ln in t[a:e].splitlines():
        code = ln.strip()
        was_cont, cont = cont, code.endswith("&")
        if was_cont or not code or code.startswith("!"):
            continue
        ok = (re.match(r"(do \w+ = .*|end do)$", code) is not None
              or all(m.split("(")[0].strip() in dropped
                     for m in re.findall(r"(?:^|;)\s*(\w+(?:\([^=]*\))?)\s*=[^=]", code)) and "=" in code)
        need(ok, f"unexpected line in the Gamma span being replaced: {code!r}")
    def load_blocks(pre, dens):
        """Density blocks of pool density `dens` the quartet touches, into
        the thread-private arrays {pre}ab ... {pre}bd (the same six blocks
        the RHF/UHF kernel caches), so Gamma below reads small local arrays
        and not nv * 8 strided pool entries per pair."""
        out = ""
        for nm, (x0, n0, x1, n1, o0, o1) in (
                ("ab", ("mui", d0, "nuj", d1, "ia", "ib")),
                ("cd", ("lamk", d2, "sigl", d3, "ic", "id")),
                ("ac", ("mui", d0, "lamk", d2, "ia", "ic")),
                ("ad", ("mui", d0, "sigl", d3, "ia", "id")),
                ("bc", ("nuj", d1, "lamk", d2, "ib", "ic")),
                ("bd", ("nuj", d1, "sigl", d3, "ib", "id"))):
            out += (f"               do {o1} = 0, {_dim(o1, dims) - 1}\n"
                    f"                  do {o0} = 0, {_dim(o0, dims) - 1}\n"
                    f"                     {pre}{nm}(1 + {o0} + {n0}*{o1}) = dmat({dens}, {x0} + {o0}, {x1} + {o1})\n"
                    f"                  end do\n"
                    f"               end do\n")
        return out

    def gam_loop(body):
        return (f"               idx = 0\n"
                f"               do dd2 = 1, {d3}\n"
                f"                  do cc = 1, {d2}\n"
                f"                     do bb = 1, {d1}\n"
                f"                        do aa = 1, {d0}\n"
                f"                           idx = idx + 1\n"
                f"                           gam(idx, col) = gam(idx, col) + 8.0_dp*wq*( &\n"
                f"{body}"
                f"                              )\n"
                f"                        end do\n"
                f"                     end do\n"
                f"                  end do\n"
                f"               end do\n")

    def blk(pre, nm, aa, bb, n0):
        return f"{pre}{nm}({aa} + {n0}*({bb} - 1))"

    cross = gam_loop(
        f"                              pj(ip)*({blk('l','ab','aa','bb',d0)}*{blk('r','cd','cc','dd2',d2)} &\n"
        f"                                     + {blk('r','ab','aa','bb',d0)}*{blk('l','cd','cc','dd2',d2)}) &\n"
        f"                              - pk(ip)*({blk('l','ac','aa','cc',d0)}*{blk('r','bd','bb','dd2',d1)} &\n"
        f"                                      + {blk('r','ac','aa','cc',d0)}*{blk('l','bd','bb','dd2',d1)} &\n"
        f"                                      + {blk('l','ad','aa','dd2',d0)}*{blk('r','bc','bb','cc',d1)} &\n"
        f"                                      + {blk('r','ad','aa','dd2',d0)}*{blk('l','bc','bb','cc',d1)})")
    # A column past ncb (a block narrower than this variant) gets no pairs:
    # its gam is zero and its tail is skipped.
    gam = f"""         do col = 1, {nc}
            do idx = 1, {nv}
               gam(idx, col) = 0.0_dp
            end do
            if (col > ncb) cycle
            do ip = pstart(col), pstart(col + 1) - 1
               jl = pl(ip); jr = pr(ip)
{load_blocks('l', 'jl')}{load_blocks('r', 'jr')}{cross}            end do
         end do

         acc = 0.0_dp
"""
    t = t[:a] + gam + t[e:]
    # Accumulation.
    acc_old = ("                                 vax = vax + fix*pyz*gam(idx); vay = vay + fiy*pxz*gam(idx); vaz = vaz + fiz*pxy*gam(idx)\n"
               "                                 vbx = vbx + fjx*pyz*gam(idx); vby = vby + fjy*pxz*gam(idx); vbz = vbz + fjz*pxy*gam(idx)\n"
               "                                 vcx = vcx + fkx*pyz*gam(idx); vcy = vcy + fky*pxz*gam(idx); vcz = vcz + fkz*pxy*gam(idx)\n"
               "                                 vdx = vdx + flx*pyz*gam(idx); vdy = vdy + fly*pxz*gam(idx); vdz = vdz + flz*pxy*gam(idx)\n")
    # One statement per column and sum, constant column indices: the sums are
    # scalarised into registers.
    acc_new = ""
    for c in range(1, nc + 1):
        for k, (f_, p_) in enumerate((("fix", "pyz"), ("fiy", "pxz"), ("fiz", "pxy"),
                                      ("fjx", "pyz"), ("fjy", "pxz"), ("fjz", "pxy"),
                                      ("fkx", "pyz"), ("fky", "pxz"), ("fkz", "pxy"),
                                      ("flx", "pyz"), ("fly", "pxz"), ("flz", "pxy")), 1):
            acc_new += (f"                                 acc({k}, {c}) = acc({k}, {c}) "
                        f"+ {f_}*{p_}*gam(idx, {c})\n")
    t = sub1(t, acc_old, acc_new)
    # Tail.
    a = t.index("         !$acc atomic update\n         gp(1, sh_at(si), slot)")
    e = t.index(f"      end subroutine rysgi{tag}_pairs")
    tail = ""
    for c in range(1, nc + 1):
        tail += f"         if (ncb >= {c}) then\n"
        for k, sh in enumerate("ijkl"):
            for ax in range(3):
                tail += ("            !$acc atomic update\n"
                         f"            gp({ax + 1}, sh_at(s{sh}), slot, {c}) = gp({ax + 1}, sh_at(s{sh}), slot, {c}) "
                         f"+ acc({3*k + ax + 1}, {c})\n")
        tail += "         end if\n"
    t = t[:a] + tail + t[e:]
    need("ndens" not in t.replace("npool", ""), "pair item still reads ndens")
    code_only = "\n".join(ln for ln in t.splitlines() if not ln.strip().startswith("!"))
    bad = [ln for ln in code_only.splitlines()
           if "jfac" in ln or "kfac" in ln.replace("kfac_eff", "")]
    need(not bad, f"pair item still reads jfac/kfac: {bad[:2]}")
    t = t.replace("      real(dp) :: alpha_i, alpha_j, alpha_k, kfac_eff\n",
                  "      real(dp) :: alpha_i, alpha_j, alpha_k\n")
    return t



def _wide_widths(comps, rt_plan, inside):
    """pairs_variants() before PAIRS_NARROW_ONLY: the widths the Gamma budget allows."""
    if rt_plan is not None:
        inner, _ = gr._chunk_keys(comps, rt_plan[0], 0)
        gm = len(inner)
    elif inside is not None:
        gm = max(len(ch) for ch in _chunks(comps, inside))
    else:
        return (1,)
    return tuple(nc for nc in PAIRS_VARIANTS if nc == 1 or nc*gm <= PAIRS_GM_BUDGET[nc])


def pairs_variants(tag, comps, rt_plan, inside):
    """The column widths class `tag`'s pair kernel is generated for: every
    width of PAIRS_VARIANTS whose Gamma count per sweep fits PAIRS_GM_BUDGET
    (nc * gm per chunk), unless the class is in PAIRS_NARROW_ONLY; width 1
    always."""
    if tag in PAIRS_NARROW_ONLY:
        return (1,)
    return _wide_widths(comps, rt_plan, inside)


def pair_widths(tag, comps, rt_plan, inside):
    """The widths class `tag`'s pair kernel is generated for: none for an
    orientation the pair bins never launch, else pairs_variants()."""
    la, lb, lc, ld = (int(c) for c in tag)
    if not (la >= lb and lc >= ld and (la, lb) >= (lc, ld)):
        # The pair bins never present this orientation (trc_bins.F90's
        # ANGULAR-MOMENTUM CANONICAL), so no kernel is built for it; the
        # dispatcher still names the class.
        return ()
    return pairs_variants(tag, comps, rt_plan, inside)


def pair_module_name(tag, nc):
    """Module holding the width-`nc` pair kernel of class `tag`."""
    return f"trc_pc_rysgp{tag}w{nc}"


def pair_driver_name(tag, nc):
    return f"rysgp{tag}w{nc}_pairs_compact"


def pairs_modules(tag, item, dims, comps, rt_plan, inside, prelude, tabs):
    """One module per width in pair_widths(), as {nc: text}. Each holds the
    pair item routine rysgi{tag}_pairs{nc} and its driver
    rysgp{tag}w{nc}_pairs_compact, which owns the `do concurrent` (that body
    may call only same-module routines). The block walk and the choice of
    the narrowest width that holds a block live in the dispatcher
    (gen_rys_grad_dispatch.py); a width is its own file so the widths
    compile in parallel, and apart from the class's other kernels.

    `prelude` is the class module's text up to `contains` (the root-table
    uses and the constants the item routines read), `tabs` the per-shell
    component tables of the array form.

    The item routines are rewritten from the per-quartet array-form
    `rysgi{tag}_compact` text `item` (everything from its header on).

    Same decode, screens, primitive loops and derivative build; what changes
    is Gamma and where it goes. Column c owns the pairs pstart(c) ..
    pstart(c+1)-1 of a table (pl, pr, pj, pk): pair p adds
        8 wq [ pj (Dl_ab Dr_cd + Dr_ab Dl_cd)
             - pk (Dl_ac Dr_bd + Dr_ac Dl_bd + Dl_ad Dr_bc + Dr_ad Dl_bc) ]
    to gam(., c). The host folds everything else into pj and pk (the
    Coulomb/exchange weights, the factor for counting both orderings of a
    cross pair, and the sign an antisymmetric pair's exchange takes), so one
    kernel serves symmetric and antisymmetric densities alike. The twelve
    running sums become acc(12, nc) for a variant nc columns wide, one
    scalarised set per column, and the tail adds column c into
    gp(:, :, slot, c): every column rides the SAME integral derivatives.
    The dispatcher runs the narrowest variant holding a block of columns.
    """
    slot_expr = ("int(mod(i - 1, int(npriv, int64))) + 1"
                 if sum(int(c) for c in tag) <= 3 else "1")
    rt_thread = rt_plan is not None and rt_plan[1] == "thread"
    rt_nch = 1
    if rt_plan is not None:
        for x in rt_plan[0]:
            rt_nch *= len(comps[x])
    seg_args = ("out_seg(mod(i - 1, nout) + 1), out_iab(mod(i - 1, nout) + 1), out_icd(mod(i - 1, nout) + 1)"
                if rt_thread else "out_seg(i), out_iab(i), out_icd(i)")
    lead = "int((i - 1)/nout) + 1, " if rt_thread else ""
    args = (seg_args + ", nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &\n"
            "                                npair, sp_i, sp_j, sp_q, pcut, omega, nbas, nao, ao_off, sh_at, natm, &\n"
            "                                nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &\n"
            "                                sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &\n"
            "                                npool, dmat, ncb, pstart, pl, pr, pj, pk, npriv, "
            "gp, " + slot_expr + ")")
    out = {}
    for nc in pair_widths(tag, comps, rt_plan, inside):
        if rt_plan is not None:
            item_nc = rt_pairs_item_text(tag, item, comps, rt_plan[0], rt_thread, nc)
        elif inside is not None:
            item_nc = chunked_pairs_item_text(tag, item, comps, inside, nc)
        else:
            item_nc = _pairs_variant(tag, item, dims, nc)
        drv = f"""
   subroutine {pair_driver_name(tag, nc)}(nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sW, &
                             npair, sp_i, sp_j, sp_q, pcut, omega, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             npool, dmat, ncb, pstart, pl, pr, pj, pk, npriv, gp)
      !! The pair-table compact pass of class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) for a block of at most
      !! {nc} columns: the survivors of the screen pass, each contracted against
      !! every column of the block at once. The caller passes the block as its
      !! own sub-table (pstart and gp start at the block's first column).
      integer(int64), intent(in) :: nout
      integer,  intent(in)    :: out_seg(*), out_iab(*), out_icd(*)
      integer,  intent(in)    :: nseg, npair, nbas, nao
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      integer,  intent(in)    :: sSA(nseg), sSB(nseg)
      logical,  intent(in)    :: sD(nseg)
      real(dp), intent(in)    :: sW(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), pcut, omega
      integer,  intent(in)    :: ao_off(nbas), sh_at(nbas), natm
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: sq_n(npair)
      integer,  intent(in)    :: npool, ncb
      real(dp), intent(in)    :: dmat(npool, nao, nao)
      integer,  intent(in)    :: pstart(ncb + 1), pl(*), pr(*)
      real(dp), intent(in)    :: pj(*), pk(*)
      integer,  intent(in)    :: npriv
      real(dp), intent(inout) :: gp(3, natm, npriv, ncb)
      integer(int64) :: i
      do concurrent(i=1:nout{'*%d' % rt_nch if rt_thread else ''})
         call rysgi{tag}_pairs{nc}({lead}{args}
      end do
   end subroutine {pair_driver_name(tag, nc)}

"""
        mod = pair_module_name(tag, nc)
        pre = prelude
        for old, new_ in (
                (f"module trc_pc_rysg{tag}\n", f"module {mod}\n"),
                (f"public :: rysg{tag}_screen, rysg{tag}_compact\n",
                 f"public :: {pair_driver_name(tag, nc)}\n"),
                (f"! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) GRADIENT kernel, Rys quadrature.",
                 f"! Class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}) GRADIENT kernel, Rys quadrature: the pair-table\n"
                 f"! pass, {nc} columns wide (gen_rys_grad.py pairs_modules).")):
            # Not an assert: `python -O` must not strip a check on spliced text.
            if pre.count(old) != 1:
                raise SystemExit(f"gen_rys_grad.py pairs_modules({tag}): expected exactly one "
                                 f"{old!r} in the module prelude, found {pre.count(old)}")
            pre = pre.replace(old, new_, 1)
        # Only the component tables this module reads and the prelude lacks.
        seen, keep = set(), []
        for ln in tabs.splitlines():
            nm = ln.split("::")[1].split("(")[0].strip()
            if nm + "(" in pre or nm in seen:
                continue
            seen.add(nm)
            keep.append(ln + "\n")
        out[nc] = (pre + "\n" + "".join(keep) + "\ncontains\n" + drv + item_nc
                   + f"\nend module {mod}\n")
    return out


def add_plan_options(ap):
    """The options that decide a class's plan, and with it which pair widths
    it has. The generator, the CMake --pair-widths query and the dispatcher
    all take them from here, so the three cannot disagree on a default; a
    full build must pass the same values to all three. (--rt-plan, which
    overrides one class, is the generator's alone.)"""
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="must match the --order the scalar HGP file in "
                         "--src-dir was generated with")
    ap.add_argument("--chunk-classes", default="auto",
                    help="classes whose RHF compact pass uses the CHUNKED "
                         "RHF COMPACT PASS: 'auto' (chunk_auto's rule), "
                         "'none', 'all', or a comma-separated list of tags")
    ap.add_argument("--chunk-live-max", type=int, default=CHUNK_LIVE_MAX,
                    help="live-set ceiling chunk_plan picks the chunk shape under")


def default_plan_args(order="lexicographic"):
    """A namespace holding every shared plan option at its default, but `order`."""
    ap = argparse.ArgumentParser()
    add_plan_options(ap)
    return ap.parse_args(["--order", order])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("class_tag", nargs="?", help="four digits, e.g. 1111 or 2222")
    add_plan_options(ap)
    # Not in add_plan_options: it is per class, so the dispatcher and the
    # --pair-widths query (one run for many classes) cannot take it.
    ap.add_argument("--rt-plan", default=None,
                    help="override GRAD_CHUNK_PLAN for this class: 'none', or the outer "
                         "shells and form, e.g. '23:thread' or ':loop' (one chunk)")
    ap.add_argument("--src-dir", default="src/generated",
                    help="where trc_pc_k{tag}.F90 (the scalar HGP kernel) "
                         "already lives")
    ap.add_argument("--unroll-max-comp", default="0",
                    help="unroll the component and root loops (see "
                         "rys_inline.py's UNROLLED MODE) for a class whose "
                         "na*nb*nc*nd is at most this; 'inf'/'all' unrolls "
                         "every class. Default 0: no class is unrolled, "
                         "i.e. today's array/loop form, unchanged")
    ap.add_argument("-o", "--output",
                    help="never a path under src/generated: regenerated at "
                         "configure time into the build directory. The pair "
                         "kernels (trc_pc_rysgp{tag}w{nc}.F90, one per width) "
                         "are written beside it")
    ap.add_argument("--pair-widths", metavar="CLASSES",
                    help="instead of generating (no class_tag, no -o), print one 'tag w1,w2,..' line per "
                         "class of this comma-separated list: the pair kernel widths "
                         "a run for that class writes (CMake reads it at configure "
                         "time to name the files). Empty after the tag: none")
    args = ap.parse_args()
    if args.pair_widths is not None:
        if args.class_tag or args.output:
            ap.error("--pair-widths takes no class_tag and no -o/--output")
        if args.rt_plan is not None:
            ap.error("--pair-widths cannot take --rt-plan (it overrides one class; "
                     "the dispatcher could not follow it)")
        if not args.pair_widths:
            ap.error("--pair-widths needs a comma-separated list of class tags")
        tags = args.pair_widths.split(",")
        for t in tags:
            if not re.fullmatch(r"[0-3]{4}", t):
                ap.error(f"--pair-widths: {t!r} is not a class tag (four digits, l <= 3)")
        for t in tags:
            print(t, ",".join(str(w) for w in class_pair_widths(t, args)))
        return
    if not args.class_tag or not args.output:
        ap.error("class_tag and -o/--output are required")
    run(args)


def class_plan(tag, args, generic_only=False):
    """(component tables, runtime-outer plan, inside chunk plan) of class
    `tag` under the chunk options of `args`: the choice run() and
    class_pair_widths() must agree on."""
    la, lb, lc, ld = (int(c) for c in tag)
    comps = [co.components(l, args.order) for l in (la, lb, lc, ld)]
    sel = args.chunk_classes
    if sel == "auto":
        chunked = chunk_auto(tag)
    elif sel in ("none", "all"):
        chunked = sel == "all"
    else:
        chunked = tag in sel.split(",")
    rt_arg = getattr(args, "rt_plan", None)
    rt_plan = GRAD_CHUNK_PLAN.get(tag) if rt_arg is None else _parse_rt_plan(rt_arg)
    if generic_only:
        chunked, rt_plan = False, None
    if rt_plan is not None:
        chunked = False
    inside = chunk_plan(comps, args.chunk_live_max) if chunked else None
    return comps, rt_plan, inside


def class_pair_widths(tag, args):
    """The widths whose pair modules run() writes for class `tag`; what
    CMake (--pair-widths) and the dispatcher list."""
    comps, rt_plan, inside = class_plan(tag, args)
    return pair_widths(tag, comps, rt_plan, inside)


def run(args, generic_only=False):
    """Emit one class's kernel module. With `generic_only` return instead
    (compact text, component tables) of the per-quartet ARRAY-form routine,
    unrolled or chunked for nobody: the text `pairs_modules` rewrites."""
    tag = args.class_tag
    la, lb, lc, ld = (int(c) for c in tag)
    ni, nk = la + lb, lc + ld
    nig, nkg = ni + 1, nk + 1       # the gradient-extended table's dims
    L = la + lb + lc + ld
    nroots = (L + 1)//2 + 1
    if nroots > 9:
        raise SystemExit(f"class {tag}: nroots={nroots} exceeds the 9 "
                          "committed root tables")
    dims = [len(co.components(l, args.order)) for l in (la, lb, lc, ld)]
    nv = dims[0]*dims[1]*dims[2]*dims[3]
    unroll = nv <= _parse_unroll_max(args.unroll_max_comp) and not generic_only

    src_path = os.path.join(args.src_dir, f"trc_pc_k{tag}.F90")
    src = open(src_path).read()

    i_item = src.index(f"   pure subroutine pcsi{tag}(")
    i_end = src.index(f"end subroutine pcsi{tag}\n") + len(f"end subroutine pcsi{tag}\n")
    whole = src[:i_end]
    assert i_item < i_end

    tval_line = "tval = rho_theta*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tval = whole.index(tval_line) + len(tval_line)
    head = whole[:i_tval]
    head = _thread_range_separation(head)

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
         ! The segment's orbit weight, one load and one multiply. It rides on
         ! the integral exactly as it does on the Fock side (W2), so one launch still covers
         ! a whole class run whatever weights it holds, and jfac/kfac stay
         ! unweighted launch scalars. One on every plan but a skeleton one.
         wq = wq*sW(seg)
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
                     u2 = theta*ru(rr)/(1.0_dp - theta*ru(rr))*rho
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

    comps, rt_plan, inside = class_plan(tag, args, generic_only)
    sc_text, rt_tabs = screen_and_compact_text(tag, head, middle, comps, inside, rt_plan)
    if generic_only:
        return sc_text, tabs

    # Emit the module prelude and the screen/compact pair only. The fused
    # rysg{tag}/rysgi{tag} in head + middle is just the text sc_text was
    # sliced from: nothing calls it, and compiling it would roughly double
    # each file's compile time for nothing.
    # The pair-table kernel, from the array-form routine of its own run: the
    # class's normal text may be unrolled or chunked, which the rewrite
    # cannot start from.
    gen_text, gen_tabs = run(args, generic_only=True)
    prelude = head[:head.index("\ncontains\n")]
    pair_mods = pairs_modules(tag, gen_text, dims, comps, rt_plan, inside, prelude, gen_tabs)
    rt_tabs = "".join(ln + "\n" for ln in rt_tabs.splitlines()
                      if ln.split("::")[1].split("(")[0].strip() + "(" not in prelude)
    prelude += "\n" + rt_tabs + "\ncontains\n"
    out = prelude + sc_text + f"\nend module trc_pc_rysg{tag}\n"
    # The pair kernels are one file per width, beside the class's own file:
    # nvfortran compiles them in parallel, and one width's compile time no
    # longer holds up the ordinary gradient kernels of its class.
    outdir = os.path.dirname(os.path.abspath(args.output))
    for nc, text in pair_mods.items():
        with open(os.path.join(outdir, pair_module_name(tag, nc) + ".F90"), "w") as fh:
            fh.write(text)
    with open(args.output, "w") as fh:
        fh.write(out)
    print(f"  {args.output}: ({la}{lb}|{lc}{ld}), {nroots} roots, "
          f"{len(out.splitlines())} lines, pair files "
          + (", ".join(pair_module_name(tag, nc) + ".F90" for nc in pair_mods) or "none")
          + (f", chunked inside={inside}" if inside is not None else "")
          + (f", runtime-outer {rt_plan}" if rt_plan is not None else ""))


def _check_narrow_only():
    """PAIRS_NARROW_ONLY names real, launched classes that would otherwise get
    a wider width: a stale entry would silently do nothing."""
    for t in PAIRS_NARROW_ONLY:
        if not re.fullmatch(r"[0-3]{4}", t):
            raise SystemExit(f"gen_rys_grad.py: PAIRS_NARROW_ONLY entry {t!r} is not a class tag")
        la, lb, lc, ld = (int(c) for c in t)
        if not (la >= lb and lc >= ld and (la, lb) >= (lc, ld)):
            raise SystemExit(f"gen_rys_grad.py: PAIRS_NARROW_ONLY entry {t}: the pair bins never launch it")
        comps, rt_plan, inside = class_plan(t, default_plan_args())
        if max(_wide_widths(comps, rt_plan, inside)) < 2:
            raise SystemExit(f"gen_rys_grad.py: PAIRS_NARROW_ONLY entry {t} has no width above 1 to drop")


_check_narrow_only()


if __name__ == "__main__":
    main()
