#!/usr/bin/env python3
"""Numerical reference check for terco's rotated-axis MD algebra.

Evaluates the two-electron Coulomb integral (ab|cd) by walking the SAME
three factorisation levels a generated Fortran kernel would (see
scripts/rotaxis/factor.py and libfint's scripts/rotaxis_mmd/emit.py), using
Python floats instead of Fortran statements, starting from the same
per-primitive-pair quantities terco's own kernels are handed (the pp_p,
pp_r, pp_ra, pp_rb, pp_c arrays built by src/trc_hgp.F90's
build_pairs_hgp), and compares the lab-frame, terco-normalised result
against pyscf's int2e_cart on random contracted shells and geometries,
including the degenerate frames a GPU kernel must not choke on.

Run as:  python3 scripts/rotaxis/check_numeric.py
Exits nonzero if any quartet's scaled error exceeds the tolerance.

WHERE THE NUMBERS COME FROM
----------------------------
terco's pp_* convention (read from src/trc_hgp.F90, src/trc_api.F90,
src/generated/trc_pc_k1111.F90 -- see FRAME.md for the derivation in full):

  pp_p(k)      p = a + b                    (bra) / q = c + d        (ket)
  pp_r(k,:)    P = (a*A + b*B)/p             (bra) / Q = (c*C+d*D)/q (ket)
  pp_ra(k,:)   A                             (bra) / C               (ket)
  pp_rb(k,:)   B                             (bra) / D               (ket)
  pp_c(k)      c_i*c_j*exp(-mu*R_AB^2)       (bra) / c_k*c_l*exp(-nu*R_CD^2) (ket)

  where c_i is terco's sh_c: the shell's contraction coefficient already
  multiplied by the primitive's own normalisation AND by libcint's
  common_fac_sp(l) (0.282094791773878143 for l=0, 0.488602511902919921 for
  l=1) -- see src/trc_api.F90:223-230,284 and docs/GENERATORS.md section
  "Shared conventions". This is verified empirically below by
  `verify_pp_c_convention`, run once at start-up: pyscf's own _env
  coefficients are the primitive+contraction-normalised value alone (no
  common_fac_sp), and multiplying by common_fac_sp(l) once per shell
  reproduces mol.intor('int2e_cart') exactly for a single-primitive
  (ss|ss) test built by hand, independent of the rest of this file's
  machinery.

  The per-quartet prefactor terco's HGP kernel applies (trc_hgp.F90:306,
  trc_batch.F90:386, trc_binkernel.F90:900) is

      pref = TWO_PI_2_5/(zeta*eta*sqrt(zeta+eta)) * pp_c(kp) * pp_c(kq)

  which is exactly 2 pi^(5/2)/(p q sqrt(p+q)) K_ab K_cd c_i c_j c_k c_l --
  the textbook MD prefactor. Splitting this the way libfint's rotaxis
  kernels want it (kab carries 1/p, kcd carries 2 pi^(5/2)/q, B_n carries
  1/sqrt(p+q)):

      kab = pp_c(kp) / pp_p(kp)
      kcd = TWO_PI_2_5 * pp_c(kq) / pp_p(kq)

  reproduces the same product bit for bit (both are just a regrouping of
  the same four factors).

BOYS FUNCTION
-------------
F_n(t) = gamma(n+1/2) / (2 t^(n+1/2)) * P(n+1/2, t), P the regularised
lower incomplete gamma function (scipy.special.gammainc), with the t=0
limit 1/(2n+1) handled directly. Every F_n is evaluated from this closed
form directly (no downward or upward recursion between orders), so there
is no recursion-cancellation failure mode to inherit from libfint's
gamma_inc_like -- this checks the ALGEBRA, not a Boys approximation.
scipy's gammainc is accurate across the whole range tested here
(t = 0 .. 1e8), which the adversarial cases below confirm.
"""
import math
import os
import sys

import numpy as np
from scipy.special import gammainc, gammaln

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rotaxis import derive, factor  # noqa: E402

from pyscf import gto  # noqa: E402

TWO_PI_2_5 = 34.986836655249725
COMMON_FAC_SP = {0: 0.282094791773878143, 1: 0.488602511902919921}
KIND_L = {"s": 0, "p": 1}

TOL = 1.0e-10  # scaled to the largest |integral| in the quartet


# ---------------------------------------------------------------------
# Boys function, closed form per order (see module docstring).
# ---------------------------------------------------------------------
def boys(n, t):
    if t < 1.0e-14:
        return 1.0 / (2 * n + 1)
    a = n + 0.5
    return math.exp(gammaln(a)) / (2.0 * t ** a) * gammainc(a, t)


def boys_array(L, t):
    return [boys(n, t) for n in range(L + 1)]


# ---------------------------------------------------------------------
# The local frame -- transcribed from libfint's cint_rotaxis_2e.f90
# `local_frame`, `bra_pairs`, `ket_pairs`. Operates on terco's pp_*
# per-primitive-pair convention directly.
# ---------------------------------------------------------------------
def local_frame(rA, rB, rC, rD):
    rA, rB, rC, rD = (np.asarray(x, dtype=float) for x in (rA, rB, rC, rD))
    ab = rB - rA
    rab = math.sqrt(float(ab @ ab))
    if rab > 0.0:
        zh = ab / rab
    else:
        zh = np.array([0.0, 0.0, 1.0])
    cd = rD - rC
    rcd = math.sqrt(float(cd @ cd))
    w = cd - (cd @ zh) * zh
    wn = math.sqrt(float(w @ w))
    if wn > 1.0e-12 * max(rcd, 1.0):
        xh = w / wn
    else:
        # CD parallel to AB (or a point): any x perpendicular to z will do.
        # Take the lab axis least aligned with z and project -- exactly
        # libfint's fallback (cint_rotaxis_2e.f90:local_frame).
        kmin = int(np.argmin(np.abs(zh)))
        w = np.zeros(3)
        w[kmin] = 1.0
        w = w - (w @ zh) * zh
        xh = w / math.sqrt(float(w @ w))
    yh = np.cross(zh, xh)
    rot = np.column_stack([xh, yh, zh])  # columns = local axes in lab coords

    cl = np.array([xh @ (rC - rA), yh @ (rC - rA), zh @ (rC - rA)])
    dl = np.array([xh @ (rD - rA), yh @ (rD - rA), zh @ (rD - rA)])
    Rs = dl[0] - cl[0]
    Rc = dl[2] - cl[2]
    cy = cl[1]
    return rot, rab, rcd, (Rs, Rc, cy), cl


def bra_pairs(prims_a, prims_b, rab):
    """prims_* = [(exponent, sh_c_coefficient), ...] for the two bra shells,
    already carrying terco's normalisation (primitive norm * contraction
    coefficient * common_fac_sp(l)). Returns a list of dicts, one per
    primitive pair, holding exactly what pp_p/pp_c would give a kernel:
    p, ip=1/(2p), za, zb, and the rotaxis "kab" weight pp_c/pp_p."""
    out = []
    for a, ca in prims_a:
        for b, cb in prims_b:
            p = a + b
            mu = a * b / p
            c = ca * cb * math.exp(-mu * rab * rab)  # this IS terco's pp_c(k)
            za = b * rab / p
            zb = za - rab
            out.append(dict(p=p, ip=0.5 / p, za=za, zb=zb, kab=c / p))
    return out


def eval_factorised(f, bras, kets, gc):
    """Numerically evaluate every component of the class, walking
    Factorised's three levels exactly as emit_kernel (libfint's
    scripts/rotaxis_mmd/emit.py) would: ket primitives outside, bra
    primitives inside, one S accumulator per distinct bra monomial, one r
    accumulator per distinct ket polynomial, then the quartet assembly.
    `bras`/`kets` are lists of per-primitive-pair dicts as built by
    `bra_pairs`/the ket loop in `evaluate_quartet`. Single contraction
    column only (NTT = 1 for every s/p class, no L shells)."""
    d = f.d
    L = d.L
    Rs, Rc, cy = gc
    NS, NR = len(f.s_list), len(f.r_list)
    ncomp = len(f.components)

    r_acc = [0.0] * NR
    for ket in kets:
        q, iq, yc, yd, ax, az, kcd = (ket[k] for k in ("q", "iq", "yc", "yd", "ax", "az", "kcd"))
        s_acc = [0.0] * NS
        for bra in bras:
            p, ip, za, zb, kab = (bra[k] for k in ("p", "ip", "za", "zb", "kab"))
            zq = za - az
            alpha = p * q / (p + q)
            T = alpha * (ax * ax + cy * cy + zq * zq)
            Fm = boys_array(L, T)
            w = 1.0 / math.sqrt(p + q)
            m2a = -2.0 * alpha
            B = []
            for n in range(L + 1):
                B.append(w * Fm[n])
                w *= m2a
            gvals = {}
            for g in f.g_list:
                e_ip, e_za, e_zb, e_ea = g
                gvals[g] = (ip ** e_ip) * (za ** e_za) * (zb ** e_zb)
            zqpow = [1.0] * (f.max_zq_power + 1)
            for kk in range(1, f.max_zq_power + 1):
                zqpow[kk] = zqpow[kk - 1] * zq
            for si, (g, n, kz, tt) in enumerate(f.s_list):
                s_acc[si] += kab * gvals[g] * B[n] * zqpow[kz]
        kmvals = {}
        for km in f.km_list:
            e_iq, e_yc, e_yd, e_ax = km
            kmvals[km] = (iq ** e_iq) * (yc ** e_yc) * (yd ** e_yd) * (ax ** e_ax)
        for ri, (terms, tt) in enumerate(f.r_list):
            acc = 0.0
            for (km_id, s_id), c in terms:
                acc += float(c) * kmvals[f.km_list[km_id]] * s_acc[s_id]
            r_acc[ri] += kcd * acc

    cmvals = {}
    for cm in f.cm_list:
        e_Rs, e_Rc, e_cy = cm
        cmvals[cm] = (Rs ** e_Rs) * (Rc ** e_Rc) * (cy ** e_cy)
    comps = [0.0] * ncomp
    for ci, terms in enumerate(f.assembly):
        acc = 0.0
        for scale, cm_id, r_id in terms:
            acc += float(scale) * cmvals[f.cm_list[cm_id]] * r_acc[r_id]
        comps[ci] = acc
    return comps


def rotate_axis(arr, axis, rot):
    """out[...,io,...] = sum_ii rot[io,ii] * arr[...,ii,...] along `axis`."""
    arr2 = np.moveaxis(arr, axis, -1)
    arr2 = arr2 @ rot.T
    return np.moveaxis(arr2, -1, axis)


def evaluate_quartet(f, kinds, centers, prims):
    """centers = (A,B,C,D) lab coordinates; prims = 4 lists of
    (exponent, terco_sh_c) for the four contracted shells, in slot order.
    Returns the lab-frame, terco-normalised block as a flat list in
    derive.py's all_components() order (index a fastest, then b, c, d)."""
    A, B, C, D = centers
    rot, rab, rcd, gc, cl = local_frame(A, B, C, D)
    Rs, Rc, cy = gc

    bras = bra_pairs(prims[0], prims[1], rab)
    kets = []
    for c, cc in prims[2]:
        for dd, cd in prims[3]:
            q = c + dd
            nu = c * dd / q
            coef = cc * cd * math.exp(-nu * rcd * rcd)  # terco's pp_c(k), ket side
            yc, yd = c / q, dd / q
            ax = cl[0] + yd * Rs
            az = cl[2] + yd * Rc
            kets.append(dict(q=q, iq=0.5 / q, yc=yc, yd=yd, ax=ax, az=az,
                              kcd=TWO_PI_2_5 * coef / q))

    local = eval_factorised(f, bras, kets, gc)

    na = len(derive.KINDS[kinds[0]])
    nb = len(derive.KINDS[kinds[1]])
    nc = len(derive.KINDS[kinds[2]])
    nd = len(derive.KINDS[kinds[3]])
    arr = np.array(local, dtype=float).reshape(nd, nc, nb, na)
    axis_of_slot = {0: 3, 1: 2, 2: 1, 3: 0}
    for slot in range(4):
        if kinds[slot] == "p":
            arr = rotate_axis(arr, axis_of_slot[slot], rot)
    return arr  # shape (nd, nc, nb, na)


def swapped_component_map(kinds):
    """index i (0-based, in Derivation(kc,kd,ka,kb)'s own all_components()
    order) -> (ia, ib, ic, id) in TERCO's own (a,b,c,d) slot-index order.

    Derivation(kc,kd,ka,kb).all_components() iterates its own four slots
    slowest-to-fastest as (kb, ka, kd, kc) -- see derive.py's nested loop,
    "for cd: for cc: for cb: for ca:" applied to kinds (kc,kd,ka,kb), so its
    slot 0 (fastest) is kc, slot 1 is kd, slot 2 is ka, slot 3 (slowest) is
    kb. Building this map once by direct construction, rather than by
    reshape/transpose of the flat result, is what makes the swap trivial to
    get right: every entry says explicitly which terco (ia,ib,ic,id) it is,
    with no axis-order arithmetic to get backwards."""
    ka, kb, kc, kd = kinds
    NA, NB, NC, ND = (len(derive.KINDS[k]) for k in (ka, kb, kc, kd))
    mapping = []
    for id_ in range(NB):       # slowest: swapped slot 3 = kb = terco's b
        for ic_ in range(NA):   # swapped slot 2 = ka = terco's a
            for ib_ in range(ND):   # swapped slot 1 = kd = terco's d
                for ia_ in range(NC):  # fastest: swapped slot 0 = kc = terco's c
                    # (ia_ indexes terco's c-slot, ib_ terco's d-slot,
                    #  ic_ terco's a-slot, id_ terco's b-slot)
                    mapping.append((ic_, id_, ia_, ib_))  # (terco a,b,c,d)
    return mapping


def evaluate_quartet_swapped(f_swapped, kinds, centers, prims):
    """The role-swapped evaluation (coordinator review, milestone-1 follow-
    up): derive Derivation(kc,kd,ka,kb) instead of Derivation(ka,kb,kc,kd),
    so its "bra" (level 1, the inner accumulation loop) is physically
    terco's KET pair (C,D) and its "ket" (level 2, outer) is physically
    terco's BRA pair (A,B). The frame is built from (C,D,A,B) instead of
    (A,B,C,D) -- z along C->D, x from AB perpendicular to CD, local
    coordinates of A (not C) at the frame's origin -- exactly `local_frame`
    called with its arguments in the other order, since (ab|cd) = (cd|ab)
    is the same integral computed from the other pair's point of view.

    `f_swapped` must be `factor.Factorised(derive.Derivation(kinds[2],
    kinds[3], kinds[0], kinds[1]))` -- the caller's, so it is built once per
    class rather than once per geometry.

    Returns the SAME shape and convention as `evaluate_quartet`: (nd, nc,
    nb, na), lab-frame, rotated, ready to `.ravel()` against pyscf. Nothing
    here reshapes or transposes the flat result -- `swapped_component_map`
    says explicitly which terco slot each computed value belongs to, and
    values are scattered directly into a terco-shaped array before the
    (unmodified) rotation step runs, so a transpose-order mistake cannot
    hide in an axis permutation."""
    A, B, C, D = centers
    rot2, rcd_dist, rab_dist, gc2, cl2 = local_frame(C, D, A, B)
    Rs2, Rc2, cy2 = gc2

    # derivation-bra = physically (C, D): the pair whose OWN separation is
    # rcd_dist (local_frame's "rab" return when called C, D first).
    bras2 = bra_pairs(prims[2], prims[3], rcd_dist)
    # derivation-ket = physically (A, B): cl2/Rs2/Rc2 are A's local
    # coordinates and (B - A)'s local offset in the NEW frame, and the
    # Gaussian-product exponential needs (A, B)'s own separation, rab_dist
    # (local_frame's "rcd" return here).
    kets2 = []
    for a, ca in prims[0]:
        for b, cb in prims[1]:
            q = a + b
            mu = a * b / q
            coef = ca * cb * math.exp(-mu * rab_dist * rab_dist)
            yc, yd = a / q, b / q
            ax = cl2[0] + yd * Rs2
            az = cl2[2] + yd * Rc2
            kets2.append(dict(q=q, iq=0.5 / q, yc=yc, yd=yd, ax=ax, az=az,
                              kcd=TWO_PI_2_5 * coef / q))

    local_swapped = eval_factorised(f_swapped, bras2, kets2, gc2)

    na = len(derive.KINDS[kinds[0]])
    nb = len(derive.KINDS[kinds[1]])
    nc = len(derive.KINDS[kinds[2]])
    nd = len(derive.KINDS[kinds[3]])
    mapping = swapped_component_map(kinds)
    assert len(mapping) == len(local_swapped) == na * nb * nc * nd
    arr = np.zeros((nd, nc, nb, na), dtype=float)
    for val, (ia, ib, ic, id_) in zip(local_swapped, mapping):
        arr[id_, ic, ib, ia] = val

    axis_of_slot = {0: 3, 1: 2, 2: 1, 3: 0}
    for slot in range(4):
        if kinds[slot] == "p":
            arr = rotate_axis(arr, axis_of_slot[slot], rot2)
    return arr  # shape (nd, nc, nb, na), identical convention to evaluate_quartet


# ---------------------------------------------------------------------
# pyscf reference
# ---------------------------------------------------------------------
ELEMENTS = ["H", "He", "Li", "Be"]


def make_mol(kinds, centers, prims_raw):
    """prims_raw[slot] = [(exp, coef), ...] RAW (un-normalised) input
    coefficients; pyscf normalises them the same way libcint does. Returns
    the Mole plus, per slot, the list of (exponent, terco_sh_c) pairs read
    back out of mol._env -- i.e. exactly what terco's basis_build would
    hand the kernel (primitive-normalised, contraction-normalised,
    common_fac_sp folded in)."""
    atom = []
    basis = {}
    for slot in range(4):
        elem = ELEMENTS[slot]
        atom.append([elem, tuple(centers[slot])])
        l = KIND_L[kinds[slot]]
        basis[elem] = [[l] + [[e, c] for e, c in prims_raw[slot]]]
    mol = gto.M(atom=atom, basis=basis, unit="Bohr", cart=True,
                spin=None, charge=0, verbose=0, dump_input=False)
    # spin=None triggers auto; some element combinations have odd electron
    # counts, so force closed-shell bookkeeping off by using UHF-agnostic
    # electron count -- irrelevant to int2e_cart, which never reads nelec.
    prims_terco = []
    for slot in range(4):
        b = mol._bas[slot]
        l = b[gto.ANG_OF]
        nprim = b[gto.NPRIM_OF]
        pe, pc = b[gto.PTR_EXP], b[gto.PTR_COEFF]
        exps = mol._env[pe:pe + nprim]
        coefs = mol._env[pc:pc + nprim] * COMMON_FAC_SP[l]
        prims_terco.append(list(zip(exps.tolist(), coefs.tolist())))
    return mol, prims_terco


def pyscf_block(mol, eri, slots=(0, 1, 2, 3)):
    ao_loc = mol.ao_loc_nr()
    a0, a1 = ao_loc[slots[0]], ao_loc[slots[0] + 1]
    b0, b1 = ao_loc[slots[1]], ao_loc[slots[1] + 1]
    c0, c1 = ao_loc[slots[2]], ao_loc[slots[2] + 1]
    d0, d1 = ao_loc[slots[3]], ao_loc[slots[3] + 1]
    return eri[a0:a1, b0:b1, c0:c1, d0:d1]


def flatten_pyscf_order(block):
    """block has numpy axes (a,b,c,d); flatten in a-fastest order to match
    derive.py's all_components()."""
    na, nb, nc, nd = block.shape
    out = np.empty(na * nb * nc * nd)
    idx = 0
    for id_ in range(nd):
        for ic in range(nc):
            for ib in range(nb):
                for ia in range(na):
                    out[idx] = block[ia, ib, ic, id_]
                    idx += 1
    return out


# ---------------------------------------------------------------------
# The check itself
# ---------------------------------------------------------------------
def verify_pp_c_convention():
    """Standalone confirmation, independent of everything else in this
    file, that terco's sh_c = (pyscf's _env coefficient) * common_fac_sp(l):
    a single s-primitive on one centre, four copies, by hand against the
    textbook MD formula (T=0 exactly since all four centres coincide, so
    F_0=1 and there is nothing to get subtly wrong in the Boys function).
    Returns (ok, computed, pyscf_reference)."""
    mol = gto.M(atom=[["He", (0.0, 0.0, 0.0)]], basis={"He": [[0, [1.3, 1.0]]]},
                unit="Bohr", cart=True, verbose=0)
    b = mol._bas[0]
    a = mol._env[b[gto.PTR_EXP]]
    c_raw = mol._env[b[gto.PTR_COEFF]]
    c = c_raw * COMMON_FAC_SP[0]
    p = q = 2.0 * a
    pref = TWO_PI_2_5 / (p * q * math.sqrt(p + q))
    computed = pref * c ** 4  # F_0(0) = 1
    ref = mol.intor("int2e_cart")[0, 0, 0, 0]
    ok = abs(computed - ref) <= 1e-13 * max(abs(ref), 1e-300)
    return ok, computed, ref


def random_prims(rng, n=2):
    return [(rng.uniform(0.3, 3.0), rng.uniform(0.3, 1.3)) for _ in range(n)]


ALL_CLASSES = [(a, b, c, d)
               for a in "sp" for b in "sp" for c in "sp" for d in "sp"]


def geometries(rng):
    """Named (centers, label) tuples: two generic random geometries plus
    every adversarial case the task calls out."""
    cases = []
    for i in range(2):
        pts = [np.array([rng.uniform(-1.5, 1.5) for _ in range(3)]) for _ in range(4)]
        cases.append((pts, f"generic-{i}"))

    # A == B: R_AB = 0, the frame's z axis is undefined (falls back to lab z)
    A = np.array([0.3, -0.2, 0.5])
    C = np.array([1.1, 0.4, -0.6])
    D = np.array([-0.7, 1.3, 0.2])
    cases.append(([A, A.copy(), C, D], "A==B"))

    # C == D: analogous on the ket side
    B = np.array([-0.9, 0.6, 1.2])
    cases.append(([A, B, C, C.copy()], "C==D"))

    # AB parallel to CD: the y axis (AB x CD) is undefined; also covers the
    # x-axis fallback in local_frame (wn ~ 0).
    A2 = np.array([0.0, 0.0, 0.0])
    B2 = np.array([0.0, 0.0, 1.7])
    C2 = np.array([0.4, -0.3, -0.5])
    D2 = C2 + np.array([0.0, 0.0, 0.9])   # D - C along the same z direction
    cases.append(([A2, B2, C2, D2], "AB-parallel-CD"))

    # all four centres coincide: both degenerate branches at once, T == 0
    P0 = np.array([0.25, -0.4, 0.1])
    cases.append(([P0, P0.copy(), P0.copy(), P0.copy()], "all-coincide"))

    # a very small Boys argument without being exactly zero: two nearly
    # identical, very diffuse primitive pairs a hair apart
    E1 = np.array([0.0, 0.0, 0.0])
    E2 = np.array([1.0e-9, 0.0, 0.0])
    cases.append(([E1, E1.copy(), E2, E2.copy()], "tiny-T"))

    # a very large Boys argument (T ~ 1e7-1e8): far-apart, tight primitives
    F1 = np.array([0.0, 0.0, 0.0])
    F2 = np.array([0.0, 0.0, 0.5])
    F3 = np.array([6000.0, 0.0, 0.0])
    F4 = np.array([6000.0, 0.0, 0.5])
    cases.append(([F1, F2, F3, F4], "huge-T"))

    return cases


def prims_for_case(rng, label):
    """Tight, high-exponent primitives for huge-T so alpha*(PQ)^2 is
    actually large; otherwise the generic random range."""
    if label == "huge-T":
        return [random_prims(rng, n=1) for _ in range(4)]
    return [random_prims(rng, n=2) for _ in range(4)]


def run():
    ok, computed, ref = verify_pp_c_convention()
    print(f"pp_c convention check (single-primitive ssss, by hand): "
          f"computed={computed:.15e} pyscf={ref:.15e} ok={ok}")
    if not ok:
        print("FAILED: terco's common_fac_sp folding does not reproduce "
              "pyscf -- stop, this is not a tolerance question.")
        return 1
    print()

    rng = np.random.default_rng(20260928)
    worst_overall = 0.0
    worst_case = None
    n_checked = 0
    failures = []

    print(f"{'class':6} {'case':16} {'worst scaled err':>18} {'max|int|':>12}")
    for kinds in ALL_CLASSES:
        d = derive.Derivation(*kinds)
        f = factor.Factorised(d)
        for centers, label in geometries(rng):
            prims_raw = prims_for_case(rng, label)
            mol, prims_terco = make_mol(kinds, centers, prims_raw)
            eri = mol.intor("int2e_cart")
            ref_block = pyscf_block(mol, eri)
            ref_flat = flatten_pyscf_order(ref_block)

            local_lab = evaluate_quartet(f, kinds, centers, prims_terco)
            # arr has axes (d,c,b,a); ravel (C order) gives a fastest -- the
            # same order flatten_pyscf_order builds.
            got_flat = local_lab.ravel()

            scale = max(np.max(np.abs(ref_flat)), 1e-300)
            err = np.max(np.abs(got_flat - ref_flat)) / scale
            n_checked += got_flat.size
            tag = "".join(kinds)
            print(f"{tag:6} {label:16} {err:18.3e} {scale:12.3e}")
            if err > worst_overall:
                worst_overall = err
                worst_case = (tag, label)
            if err > TOL:
                failures.append((tag, label, err))

    n_geom = len(geometries(np.random.default_rng(0)))  # count only, throwaway rng
    print()
    print(f"checked {n_checked} scalar components over "
          f"{len(ALL_CLASSES)} classes x {n_geom} geometries per class")
    print(f"worst scaled error: {worst_overall:.3e}  at {worst_case}")
    if failures:
        print(f"FAILED: {len(failures)} (class, geometry) pairs exceeded "
              f"tolerance {TOL:.1e}:")
        for tag, label, err in failures:
            print(f"  {tag:6} {label:16} {err:.3e}")
        return 1
    print(f"PASSED: every (class, geometry) pair scaled error <= {TOL:.1e}")
    return 0


def run_swapped():
    """The role-swapped evaluation (coordinator review, milestone-1 follow-
    up): Derivation(kc,kd,ka,kb) instead of Derivation(ka,kb,kc,kd), frame
    from (C,D,A,B). Same 16 classes, same 8 geometries (the existing
    degenerate cases already swap roles correctly here too: "C==D" is the
    frame's undefined-z case for THIS evaluator, exactly as "A==B" is for
    the unswapped one, since local_frame is called with its first two
    arguments swapped and does not know or care which physical pair they
    came from), same tolerance. A separate function, run only after `run()`
    has already passed, so the original check is never weakened or
    reordered by this one."""
    rng = np.random.default_rng(20260928)  # same seed: same geometries/prims
    worst_overall = 0.0
    worst_case = None
    n_checked = 0
    failures = []

    print(f"{'class':6} {'case':16} {'worst scaled err':>18} {'max|int|':>12}")
    for kinds in ALL_CLASSES:
        ka, kb, kc, kd = kinds
        d_swapped = derive.Derivation(kc, kd, ka, kb)
        f_swapped = factor.Factorised(d_swapped)
        for centers, label in geometries(rng):
            prims_raw = prims_for_case(rng, label)
            mol, prims_terco = make_mol(kinds, centers, prims_raw)
            eri = mol.intor("int2e_cart")
            ref_block = pyscf_block(mol, eri)
            ref_flat = flatten_pyscf_order(ref_block)

            local_lab = evaluate_quartet_swapped(f_swapped, kinds, centers, prims_terco)
            got_flat = local_lab.ravel()

            scale = max(np.max(np.abs(ref_flat)), 1e-300)
            err = np.max(np.abs(got_flat - ref_flat)) / scale
            n_checked += got_flat.size
            tag = "".join(kinds)
            print(f"{tag:6} {label:16} {err:18.3e} {scale:12.3e}")
            if err > worst_overall:
                worst_overall = err
                worst_case = (tag, label)
            if err > TOL:
                failures.append((tag, label, err))

    n_geom = len(geometries(np.random.default_rng(0)))
    print()
    print(f"checked {n_checked} scalar components over "
          f"{len(ALL_CLASSES)} classes x {n_geom} geometries per class (SWAPPED)")
    print(f"worst scaled error: {worst_overall:.3e}  at {worst_case}")
    if failures:
        print(f"FAILED (swapped): {len(failures)} (class, geometry) pairs exceeded "
              f"tolerance {TOL:.1e}:")
        for tag, label, err in failures:
            print(f"  {tag:6} {label:16} {err:.3e}")
        return 1
    print(f"PASSED (swapped): every (class, geometry) pair scaled error <= {TOL:.1e}")
    return 0


if __name__ == "__main__":
    rc = run()
    if rc == 0:
        print()
        print("=" * 70)
        print("ROLE-SWAPPED EVALUATION (coordinator review, milestone-1 follow-up)")
        print("=" * 70)
        rc = run_swapped()
    sys.exit(rc)
