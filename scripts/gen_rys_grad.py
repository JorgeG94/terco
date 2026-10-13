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
    column..." decode and `#ifdef TRC_NO_DIGEST` block are not reused at
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


def screen_and_compact_text(tag, head, middle):
    """ORDERED SURVIVOR COMPACTION (ps|ss) prototype, generalised to every
    class: split the fused rysg{tag}/rysgi{tag} this module already built
    (in `head`/`middle`) into

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

    screen = f"""
   !
   ! ORDERED SURVIVOR COMPACTION -- screen pass, class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}).
   ! See screen_and_compact_text's docstring (scripts/gen_rys_grad.py) for
   ! the technique and the measurement it is based on.
   !
   subroutine rysg{tag}_screen(lo, hi, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                             out_seg, out_iab, out_icd, nmax, ctr)
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
      !: Run length a chunk screens per `do concurrent` iteration: one
      !: warp's worth. Measured against 64/128/256 on (ps|ss), w150
      !: 6-31G*, V100 -- 32 matches the warp width and wins; see this
      !: file's module-level screen_and_compact_text docstring.
      integer, parameter :: RUNLEN = 32
      integer(int64) :: g0, g1, nr, nchunks, c, clen
      g0 = sOff(lo) + 1
      g1 = sOff(hi + 1)
      nr = 0
      if (g1 >= g0) nr = g1 - g0 + 1
      nchunks = (nr + RUNLEN - 1)/RUNLEN
      do concurrent(c=1:nchunks)
         call rysgi{tag}_screen_chunk(g0 + (c - 1)*RUNLEN, &
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

    compact = f"""
   !
   ! ORDERED SURVIVOR COMPACTION -- compact pass, class ({tag[0]}{tag[1]}|{tag[2]}{tag[3]}).
   ! The same integral body rysgi{tag} above runs, sliced to start right
   ! after the decode and the three quartet-level screens: seg/iab/icd
   ! come in directly instead of being decoded from a work-item index,
   ! and the Schwarz/density screens are not repeated (the screen pass
   ! already applied them). See screen_and_compact_text's docstring.
   !
   subroutine rysg{tag}_compact(nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             ndens, dmat, grad)
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
      real(dp), intent(inout) :: grad(3, natm)
      integer(int64) :: i
      do concurrent(i=1:nout)
         call rysgi{tag}_compact(out_seg(i), out_iab(i), out_icd(i), nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &
                                nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                                sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                ndens, dmat, grad)
      end do
   end subroutine rysg{tag}_compact

   pure subroutine rysgi{tag}_compact(seg, iab, icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                             npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &
                             nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                             sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                             ndens, dmat, grad)
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
      real(dp), intent(inout) :: grad(3, natm)
{compact_decls}
{si_only_block}{rest_after_density}{compact_middle}"""

    return screen + compact


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

    tval_line = "tval = rho*(pqx*pqx + pqy*pqy + pqz*pqz)\n"
    i_tval = whole.index(tval_line) + len(tval_line)
    head = whole[:i_tval]

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
            + ri.rys2d_hrr_decls(la + 1, lb, nkg))
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
                       f"   use trc_boys, only: dp\n" + ri.roots_table_use(nroots),
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
                  f"public :: rysg{tag}\n"
                  f"   public :: rysg{tag}_screen, rysg{tag}_compact\n",
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
    # above instead of called (see scripts/rys_inline.py), so `contains`
    # goes right back to introducing rysg{tag}/rysgi{tag}.
    helpers = f"""
{tabs}
contains
"""

    assert head.count("\ncontains\n") == 1
    head = head.replace("\ncontains\n", helpers, 1)

    sc_text = screen_and_compact_text(tag, head, middle)

    out = head + middle + sc_text + f"\nend module trc_pc_rysg{tag}\n"
    with open(args.output, "w") as fh:
        fh.write(out)
    print(f"  {args.output}: ({la}{lb}|{lc}{ld}), {nroots} roots, "
          f"{len(out.splitlines())} lines")


if __name__ == "__main__":
    main()
