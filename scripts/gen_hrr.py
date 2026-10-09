#!/usr/bin/env python3
"""
Emit the unrolled horizontal transfer + digestion, as `case` blocks to include
directly inside the hot loop.

WHY
---
The generic HRR is the worst-shaped loop in the library.  Per Cartesian
component it calls `comp` to decode the component, then runs a six-deep sum in
which every term recomputes three `binom` calls and three `powi` calls and a
`cidx` table lookup -- and `binom`/`powi` were cross-module until recently, so
they were real calls in the innermost loop.  The ablation put the whole
digestion block at 31% of runtime.

Every one of those is known at generation time.  The binomial products are
literal constants, the `cidx` lookups are literal indices, and the powers of
AB and CD are short products of scalars already in registers.  What is left is
a straight-line dot product per component:

    vv = c1*abx*cdy*g(17) + c2*g(23) + ...

INCLUDE, NOT A MODULE
---------------------
The VRR unrolling measured 2.4x SLOWER as a separate module and 1.58x faster as
an include, because NVHPC will not inline an `!$acc routine seq` across a module
boundary and `g` then has to cross the call in memory.  Same reasoning applies
here, so this is an include from the start.

THE HRR IS A SIEVED DAG, THE SAME WAY THE VRR ALREADY IS
----------------------------------------------------------
`emit_class` above turns the straight-line dot product per component into
literal numbers, but it still computes every one of the na*nb*nc*nd outputs
as its OWN independent closed-form sum: two outputs that share a partial
bra or ket transfer recompute it from scratch. gen_vrr.py's sieve already
builds the VRR as a DAG of shared intermediates reachable from the targets
it actually needs, pruned to the cheapest route to each one
(docs/GENERATORS.md, "Gill's generation-step sieve"); `emit_class_dag`
below is that same construct applied one step further out, where it had
not yet been: build the HRR back from the targets to the VRR buffer
(the (e0|f0) entries), one node per elementary bra or ket transfer, shared
by every output that reaches it, and emit the nodes in the order that
keeps the fewest of them live at once (its own docstring has the measured
counts). `emit_class` is kept only as this rewrite's numeric oracle
(scripts/check_hrr_dag_numeric.py evaluates both against each other on
random buffer values) -- no generator's command line reaches it any more.
"""

import argparse

import cartorder


def cart(l, order="lexicographic"):
    """(nx,ny,nz) for one shell, in `order`.

    Thin wrapper: the order itself is the single definition in cartorder.py.
    Called with NO explicit order (the default) at every INTERNAL use below --
    `cart_cum`'s index into `g`, the VRR intermediate array -- so those stay
    pinned to lexicographic regardless of what `--order` asks for. Only
    `emit_class`'s own OUTPUT enumeration (the digestion's AO slot order)
    takes `order` explicitly. See docs/GENERATORS.md's rule on OUTPUT vs
    INTERNAL enumerations.
    """
    return cartorder.components(l, order)


def cart_cum(lmax):
    out = []
    for l in range(lmax + 1):
        out.extend(cart(l))
    return {p: i + 1 for i, p in enumerate(out)}


def ncum(l):
    return (l + 1) * (l + 2) * (l + 3) // 6


def binom(n, k):
    r = 1
    for i in range(1, k + 1):
        r = r * (n - k + i) // i
    return r


def powexpr(var, n):
    """`abx**n` written out; empty string means the factor is 1."""
    if n == 0:
        return ""
    return "*".join([var] * n)


def _axis(t):
    """First nonzero axis (0=x,1=y,2=z) of a Cartesian exponent tuple, or
    None if t is (0,0,0). The direction choice in build_class_dag: fixed,
    not value-dependent, so every recursion path that can reach the same
    (a,b,c,d) node reaches it by reducing the SAME axis first -- which is
    what makes memoisation on that tuple actually share work instead of
    building two unequal routes to equal content."""
    for k in range(3):
        if t[k] > 0:
            return k
    return None


def _dec(t, k):
    t = list(t)
    t[k] -= 1
    return tuple(t)


def _inc(t, k):
    t = list(t)
    t[k] += 1
    return tuple(t)


def _make_hrr_rec(la, lb, lc, ld, idx, on_new_node):
    """Factory for the memoised HRR recursion, shared by every caller that
    needs the DAG: `build_class_dag` (structure only, for the numeric check
    and for a role's closure), `emit_class_dag` (Fortran text, with its own
    emission order), and `gen_perclass.py`'s role adaptation.

    Bra: (a,b|c,d) = (a+1_k, b-1_k|c,d) + AB_k*(a, b-1_k|c,d), k the first
    nonzero axis of b. Ket the same on (c,d) with CD, once b is (0,0,0).
    Terminal, b=d=(0,0,0): the VRR's own output, g(e,f) with e=a, f=c
    (both now totalling the class's full lab/lcd on one centre each).

    A node is a 4-tuple (a,b,c,d); recursion is memoised on that literal
    tuple, which is the DAG's only sharing axis -- see the module
    docstring's example `(a, b-1_k|c,d)`. This is a strict subset of the
    sharing a fully separated bra/ket factorisation could give (nodes do
    not vary freely over ALL of c,d while a,b transfer, only over the
    targets actually asked for), chosen because the alternative multiplies
    live temporaries by the far side's component count -- exactly the
    register pressure these kernels are already walled by (see
    docs/GENERATORS.md and the sieve's own history).

    `on_new_node(nid, r1, coefvar, r2)` fires once, the FIRST time a node is
    built -- never on a memo hit -- which is what lets the caller choose
    emission order: a caller that walks targets in a fixed order and emits
    each new node right when it is first demanded gets a DAG interleaved
    with its consumers, rather than every node built before any is used.
    See `emit_class_dag`'s docstring for why that interleaving matters.

    Returns `rec(a, b, c, d) -> ref`, where a ref is ('g', slot) for a VRR
    leaf or ('h', id) for a DAG node; slot is 1-based into the contracted
    VRR buffer `g`, exactly as `emit_class` addresses it."""
    nca = ncum(la + lb)
    memo = {}
    counter = [0]

    def rec(a, b, c, d):
        key = (a, b, c, d)
        hit = memo.get(key)
        if hit is not None:
            return hit
        if b == (0, 0, 0) and d == (0, 0, 0):
            ref = ("g", idx[a] + (idx[c] - 1)*nca)
        elif b != (0, 0, 0):
            k = _axis(b)
            r1 = rec(_inc(a, k), _dec(b, k), c, d)
            r2 = rec(a, _dec(b, k), c, d)
            counter[0] += 1
            nid = counter[0]
            coef = ("abx", "aby", "abz")[k]
            on_new_node(nid, r1, coef, r2)
            ref = ("h", nid)
        else:
            k = _axis(d)
            r1 = rec(a, b, _inc(c, k), _dec(d, k))
            r2 = rec(a, b, c, _dec(d, k))
            counter[0] += 1
            nid = counter[0]
            coef = ("cdx", "cdy", "cdz")[k]
            on_new_node(nid, r1, coef, r2)
            ref = ("h", nid)
        memo[key] = ref
        return ref

    return rec


def build_class_dag(la, lb, lc, ld, idx):
    """The HRR as a DAG with shared intermediates, for one (la,lb,lc,ld).
    Structure only (see `_make_hrr_rec`): builds every node reachable from
    the na*nb*nc*nd actual output targets, in dependency (topological)
    order -- children always precede the parent, since a node is appended
    only after both the recursive calls that build its children return.

    Returns (nodes, rec): `nodes` is [(id, ref1, coefvar, ref2), ...];
    `rec(a, b, c, d)` gives any tuple's own reference, target or not (a
    memo hit after this call, so cheap)."""
    nodes = []
    rec = _make_hrr_rec(la, lb, lc, ld, idx,
                        lambda nid, r1, coef, r2: nodes.append((nid, r1, coef, r2)))
    for a in cart(la):
        for b in cart(lb):
            for c in cart(lc):
                for d in cart(ld):
                    rec(a, b, c, d)
    return nodes, rec


def _ref_text(ref):
    kind, val = ref
    return f"g({val})" if kind == "g" else f"h({val})"


def emit_class_dag(la, lb, lc, ld, idx, order="lexicographic"):
    """DAG form of `emit_class`: identical semantics, slots and component
    order (see that function's docstring for the OUTPUT/INTERNAL split --
    the same rule applies here), but the na*nb*nc*nd outputs share whatever
    the HRR recursion found reusable instead of each being an independent
    closed-form sum.

    EMISSION ORDER IS INTERLEAVED, NOT "ALL NODES THEN ALL OUTPUTS". A
    node's live range runs from its own statement to the LAST statement
    that reads it; building the whole DAG before writing any vbuf makes
    every node live until the class's final output, which measured as the
    dominant cost -- peak simultaneously-live h(...) entries for (dd|dd)
    was 1406 that way. Emitting each output's CURRENTLY-MISSING nodes right
    before that output, in the existing OUTPUT loop order (which already
    groups every (iia,iib) under a fixed (iic,iid), i.e. under a fixed
    ket target -- exactly the DAG's own sharing axis, see
    `_make_hrr_rec`), lets most nodes die within the block that created
    them: measured peak 637 for (dd|dd), 12 for (dp|dp) (was 326), 6 for
    (pp|pp) (was 83). Still not free of risk on (dd|dd) -- see
    docs/GENERATORS.md's measured numbers -- but a straightforward, large
    win everywhere else.

    Returns (text, nnodes). `text` is `h(k) = ...` node definitions
    interleaved with `vbuf(slot) = ...` output writes, each output now a
    single reference, not a sum."""
    na_, nb_, nc_ = len(cart(la)), len(cart(lb)), len(cart(lc))
    o = []
    w = o.append
    nnodes = [0]

    def on_new(nid, r1, coef, r2):
        nnodes[0] += 1
        w(f"         h({nid}) = {_ref_text(r1)} + {coef}*{_ref_text(r2)}")

    rec = _make_hrr_rec(la, lb, lc, ld, idx, on_new)

    for iid, (dx, dy, dz) in enumerate(cart(ld, order)):
        for iic, (cx, cy, cz) in enumerate(cart(lc, order)):
            for iib, (bx, by, bz) in enumerate(cart(lb, order)):
                for iia, (ax, ay, az) in enumerate(cart(la, order)):
                    slot = 1 + iia + na_*(iib + nb_*(iic + nc_*iid))
                    ref = rec((ax, ay, az), (bx, by, bz), (cx, cy, cz), (dx, dy, dz))
                    w(f"         vbuf({slot}) = {_ref_text(ref)}")
    return "\n".join(o), nnodes[0]


def emit_class(la, lb, lc, ld, idx, order="lexicographic"):
    """Straight-line HRR + digestion for one (la,lb,lc,ld).

    `order` governs ONLY the four outer loops below: `iia`/`iib`/`iic`/`iid`
    become the AO offsets `mu = mui + iia` etc at the call site
    (trc_binkernel.F90's `#include "trc_hrr_cases.inc"`), so this is OUTPUT.
    `idx`, built by `cart_cum` above with no order argument, addresses `g` --
    the VRR's intermediate array -- and stays lexicographic regardless.
    """
    lab, lcd = la + lb, lc + ld
    nca = ncum(lab)
    o = []
    w = o.append

    for iid, (dx, dy, dz) in enumerate(cart(ld, order)):
        for iic, (cx, cy, cz) in enumerate(cart(lc, order)):
            for iib, (bx, by, bz) in enumerate(cart(lb, order)):
                for iia, (ax, ay, az) in enumerate(cart(la, order)):
                    terms = []
                    for px in range(bx + 1):
                        for py in range(by + 1):
                            for pz in range(bz + 1):
                                cb = binom(bx, px)*binom(by, py)*binom(bz, pz)
                                ja = idx[(ax + px, ay + py, az + pz)]
                                fb = [powexpr("abx", bx - px),
                                      powexpr("aby", by - py),
                                      powexpr("abz", bz - pz)]
                                for qx in range(dx + 1):
                                    for qy in range(dy + 1):
                                        for qz in range(dz + 1):
                                            cd = binom(dx, qx)*binom(dy, qy)*binom(dz, qz)
                                            jc = idx[(cx + qx, cy + qy, cz + qz)]
                                            fd = [powexpr("cdx", dx - qx),
                                                  powexpr("cdy", dy - qy),
                                                  powexpr("cdz", dz - qz)]
                                            fac = [f for f in fb + fd if f]
                                            c = cb*cd
                                            head = "" if c == 1 else f"{float(c)}_dp*"
                                            body = "*".join(fac)
                                            if body:
                                                body += "*"
                                            terms.append(
                                                f"{head}{body}g({ja + (jc - 1)*nca})")
                    #
                    # ONLY the dot product is unrolled.  An earlier version
                    # emitted the 40-line digestion block per component too,
                    # which for (pp|pp) is 81 copies of it -- that measured 17%
                    # SLOWER than the generic HRR (5.76 s against 4.91 s at 75
                    # waters).  Code bloat without work removal: the ablation
                    # had already shown the atomics are free, so replicating
                    # them buys nothing and costs instruction cache.
                    #
                    # The scatter stays a single loop over vbuf below.
                    #
                    na_, nb_, nc_ = len(cart(la)), len(cart(lb)), len(cart(lc))
                    slot = 1 + iia + na_*(iib + nb_*(iic + nc_*iid))
                    expr = terms[0]
                    for t in terms[1:]:
                        expr += f" &\n            + {t}"
                    w(f"         vbuf({slot}) = {expr}")
    return "\n".join(o)


# The eight-permutation digestion, unchanged from the generic path.  dij/dkl/dpq
# are per-quartet runtime values (they depend on the shell indices, not on the
# class), so they stay conditional; the ablation showed the atomics themselves
# are free.
DIGEST = """         !$acc atomic update
         jmat(mu, nu) = jmat(mu, nu) + vv*dmat(lam, sig)
         !$acc atomic update
         kmat(mu, lam) = kmat(mu, lam) + vv*dmat(nu, sig)
         if (dij) then
            !$acc atomic update
            jmat(nu, mu) = jmat(nu, mu) + vv*dmat(lam, sig)
            !$acc atomic update
            kmat(nu, lam) = kmat(nu, lam) + vv*dmat(mu, sig)
         end if
         if (dkl) then
            !$acc atomic update
            jmat(mu, nu) = jmat(mu, nu) + vv*dmat(sig, lam)
            !$acc atomic update
            kmat(mu, sig) = kmat(mu, sig) + vv*dmat(nu, lam)
         end if
         if (dij .and. dkl) then
            !$acc atomic update
            jmat(nu, mu) = jmat(nu, mu) + vv*dmat(sig, lam)
            !$acc atomic update
            kmat(nu, sig) = kmat(nu, sig) + vv*dmat(mu, lam)
         end if
         if (dpq) then
            !$acc atomic update
            jmat(lam, sig) = jmat(lam, sig) + vv*dmat(mu, nu)
            !$acc atomic update
            kmat(lam, mu) = kmat(lam, mu) + vv*dmat(sig, nu)
            if (dkl) then
               !$acc atomic update
               jmat(sig, lam) = jmat(sig, lam) + vv*dmat(mu, nu)
               !$acc atomic update
               kmat(sig, mu) = kmat(sig, mu) + vv*dmat(lam, nu)
            end if
            if (dij) then
               !$acc atomic update
               jmat(lam, sig) = jmat(lam, sig) + vv*dmat(nu, mu)
               !$acc atomic update
               kmat(lam, nu) = kmat(lam, nu) + vv*dmat(sig, mu)
            end if
            if (dij .and. dkl) then
               !$acc atomic update
               jmat(sig, lam) = jmat(sig, lam) + vv*dmat(nu, mu)
               !$acc atomic update
               kmat(sig, nu) = kmat(sig, nu) + vv*dmat(lam, mu)
            end if
         end if"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=1)
    ap.add_argument("--order", choices=("lexicographic", "gamess"),
                    default="lexicographic",
                    help="Cartesian order of the AO slots the digestion "
                         "below scatters into. INTERNAL indices into `g` "
                         "(the VRR's intermediate array) are unaffected: "
                         "see cart_cum and emit_class's own docstring.")
    ap.add_argument("-o", "--output", default="src/trc_hrr_cases.inc")
    ap.add_argument("--hmax-output", default="src/trc_hrr_hmax.inc",
                    help="Where to write the TRC_HRR_HMAX macro the shared "
                         "kernel declares `h(...)` with. The shared kernel "
                         "is one procedure for every class (a `select "
                         "case`), so its `h` buffer must be sized for the "
                         "worst class, unlike a per-class kernel's own "
                         "local declaration.")
    args = ap.parse_args()
    L = args.lmax
    order = args.order
    order_val = 0 if order == "lexicographic" else 1
    idx = cart_cum(4*L)

    # Suppressed when default, like every generator here now: the shipped
    # set's header must stay byte-identical to what it was before --order
    # existed, or "regenerates byte for byte with default options"
    # (docs/GENERATORS.md) stops being true.
    order_arg = "" if order == "lexicographic" else f" --order {order}"

    out = [f"""!
! Unrolled horizontal transfer + digestion, as `case` blocks to include INSIDE
! the hot loop's own procedure.
!
! GENERATED by scripts/gen_hrr.py --lmax {L}{order_arg} -- do not edit.
!
!
! STALE-INCLUDE GUARD. The `case` labels below were computed with the angular
! momentum limit this file was generated for, and the dispatch key is built
! from TRC_LMAX at compile time. If the two disagree the keys land on the
! wrong labels -- or on none, since there is no `case default` -- and the
! result is silently wrong for some classes and right for others, which is
! far worse than a failure.
!
! That is not hypothetical: this file was generated once for lmax 1 and
! included unchanged in an LMAX=2 build, where the key base went from 2 to 3
! and (pp|pp) asked for a label that did not exist. The folded path never
! noticed because it uses the per-class kernels; only the enumerated
! `fock_nosym` reaches this code, and it was wrong for a month.
!
#if TRC_LMAX != {L}
#error "This .inc was generated for a different TRC_LMAX. Regenerate it."
#endif
!
! TRC_CART_ORDER GUARD, same idea, a different compile-time constant: the AO
! slot each component scatters to below was fixed at generation time by
! --order (0 = lexicographic, the default; 1 = gamess), and TRC_CART_ORDER
! selects it at compile time. An undefined TRC_CART_ORDER reads as 0 in a
! preprocessor #if, so the default build needs no new -D flag.
!
#if TRC_CART_ORDER != {order_val}
#error "This .inc was generated for a different TRC_CART_ORDER. Regenerate it."
#endif
! Every output shares DAG nodes `h(...)` with the other outputs that need
! them (see emit_class_dag below), the same sieved-DAG treatment gen_vrr.py
! already gives the VRR, so each component is a short chain of two-term
! sums over `g` and `h` instead of its own closed-form dot product. `h`
! must be declared by the including file, sized TRC_HRR_HMAX (see
! trc_hrr_hmax.inc, generated alongside this file).
!
! An include rather than a module for the reason the VRR established by
! measurement: the same generated code ran 2.4x slower behind a module boundary
! and 1.58x faster inlined, because NVHPC does not inline `acc routine seq`
! across modules.
!
! Expects in scope: g, h, vv, mu/nu/lam/sig, mui/nuj/lamk/sigl, abx/aby/abz,
! cdx/cdy/cdz, dij/dkl/dpq, jmat, kmat, dmat.
!
      select case (((la*(LMAX + 1) + lb)*(LMAX + 1) + lc)*(LMAX + 1) + ld)"""]

    n = 0
    hmax = 0
    for la in range(L + 1):
        for lb in range(L + 1):
            for lc in range(L + 1):
                for ld in range(L + 1):
                    key = ((la*(L + 1) + lb)*(L + 1) + lc)*(L + 1) + ld
                    out.append(f"      case ({key})")
                    body, nnodes = emit_class_dag(la, lb, lc, ld, idx, order)
                    hmax = max(hmax, nnodes)
                    out.append(body)
                    n += 1
    out.append("      end select")

    txt = "\n".join(out)
    with open(args.output, "w") as fh:
        fh.write(txt + "\n")
    print(f"wrote {args.output}: {n} classes, {txt.count(chr(10))} lines")

    hmax = max(hmax, 1)   # a zero-size array is legal but needless care
    hmax_txt = f"""!
! TRC_HRR_HMAX: the largest number of shared HRR DAG nodes any one class in
! this file needs.
! The shared kernel (trc_binkernel.F90) is one procedure for every class,
! so its `h` buffer is sized once, for the worst case, unlike a per-class
! kernel's own exactly-sized local declaration.
!
! GENERATED by scripts/gen_hrr.py --lmax {L}{order_arg} -- do not edit.
!
#if TRC_LMAX != {L}
#error "This .inc was generated for a different TRC_LMAX. Regenerate it."
#endif
#define TRC_HRR_HMAX {hmax}
"""
    with open(args.hmax_output, "w") as fh:
        fh.write(hmax_txt)
    print(f"wrote {args.hmax_output}: TRC_HRR_HMAX={hmax}")


if __name__ == "__main__":
    main()
