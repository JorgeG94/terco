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
"""

import argparse


def cart(l):
    """(nx,ny,nz) for one shell, libcint order."""
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


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


def emit_class(la, lb, lc, ld, idx):
    """Straight-line HRR + digestion for one (la,lb,lc,ld)."""
    lab, lcd = la + lb, lc + ld
    nca = ncum(lab)
    o = []
    w = o.append

    for iid, (dx, dy, dz) in enumerate(cart(ld)):
        for iic, (cx, cy, cz) in enumerate(cart(lc)):
            for iib, (bx, by, bz) in enumerate(cart(lb)):
                for iia, (ax, ay, az) in enumerate(cart(la)):
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
    ap.add_argument("-o", "--output", default="src/trc_hrr_cases.inc")
    args = ap.parse_args()
    L = args.lmax
    idx = cart_cum(4*L)

    out = [f"""!
! Unrolled horizontal transfer + digestion, as `case` blocks to include INSIDE
! the hot loop's own procedure.
!
! GENERATED by scripts/gen_hrr.py --lmax {L} -- do not edit.
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
! The generic form calls `comp` per component and then runs a six-deep sum
! recomputing three `binom`, three `powi` and a `cidx` lookup per term.  All of
! it is known here, so each component becomes one straight-line dot product
! over `g` with literal coefficients and literal indices.
!
! An include rather than a module for the reason the VRR established by
! measurement: the same generated code ran 2.4x slower behind a module boundary
! and 1.58x faster inlined, because NVHPC does not inline `acc routine seq`
! across modules.
!
! Expects in scope: g, vv, mu/nu/lam/sig, mui/nuj/lamk/sigl, abx/aby/abz,
! cdx/cdy/cdz, dij/dkl/dpq, jmat, kmat, dmat.
!
      select case (((la*(LMAX + 1) + lb)*(LMAX + 1) + lc)*(LMAX + 1) + ld)"""]

    n = 0
    for la in range(L + 1):
        for lb in range(L + 1):
            for lc in range(L + 1):
                for ld in range(L + 1):
                    key = ((la*(L + 1) + lb)*(L + 1) + lc)*(L + 1) + ld
                    out.append(f"      case ({key})")
                    out.append(emit_class(la, lb, lc, ld, idx))
                    n += 1
    out.append("      end select")

    txt = "\n".join(out)
    with open(args.output, "w") as fh:
        fh.write(txt + "\n")
    print(f"wrote {args.output}: {n} classes, {txt.count(chr(10))} lines")


if __name__ == "__main__":
    main()
