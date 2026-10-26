"""Hermite index sets and the McMurchie-Davidson E coefficients.

A product of two Cartesian Gaussians on A and B with exponents a and b is a
Gaussian on P = (aA + bB)/p, p = a + b, times K_AB = exp(-ab/p |AB|^2), and
it expands exactly in Hermite Gaussians Lambda_tuv on P:

    G_a G_b = K_AB sum_tuv E^ab_tuv Lambda_tuv,
    E^ab_tuv = E^{ax bx}_t E^{ay by}_u E^{az bz}_v.

The one-dimensional coefficients come from (Helgaker, Jorgensen and Olsen,
eqs 9.5.6-9.5.7, with K_AB taken out: terco's primitive-pair coefficient
already carries it, see build_pairs_hgp)

    E^{00}_0      = 1
    E^{i+1,j}_t   = 1/(2p) E^{ij}_{t-1} + X_PA E^{ij}_t + (t+1) E^{ij}_{t+1}
    E^{i,j+1}_t   = 1/(2p) E^{ij}_{t-1} + X_PB E^{ij}_t + (t+1) E^{ij}_{t+1}

and E^{ij}_t = 0 outside 0 <= t <= i + j.

Statements are kept in one small linear form, shared with `rtree`, so the
same list can be emitted as Fortran or executed in Python by the reference:

    lhs = sum_k  const_k * factor_k * operand_k

with `factor` a named scalar of the primitive pair (or None) and `operand`
the name of an earlier statement (or None for a bare constant).
"""

from dataclasses import dataclass
from functools import lru_cache

DIRS = "xyz"


def nherm(l):
    """Hermite functions of total degree <= l."""
    return (l + 1) * (l + 2) * (l + 3) // 6


@lru_cache(maxsize=None)
def herm(l):
    """(t, u, v) with t + u + v <= l, degree-major, t then u descending.

    Same order as trc_tables's hidx, so a flat Hermite index means the same
    thing to the generated code and to the hand-written tables. Degree-major
    is what lets "every index of degree <= d" be the first nherm(d) entries.
    """
    out = []
    for n in range(l + 1):
        for t in range(n, -1, -1):
            for u in range(n - t, -1, -1):
                out.append((t, u, n - t - u))
    return tuple(out)


@lru_cache(maxsize=None)
def herm_index(l):
    """(t, u, v) -> 0-based position in herm(l)."""
    return {r: k for k, r in enumerate(herm(l))}


@dataclass(frozen=True)
class Term:
    const: float
    factor: str | None
    operand: str | None


@dataclass(frozen=True)
class Stmt:
    lhs: str
    terms: tuple


def e_name(d, i, j, t):
    return f"e{DIRS[d]}{i}{j}_{t}"


def e1d_statements(la, lb, d, needed=None):
    """E^{ij}_t in direction d for i <= la, j <= lb, as straight-line code.

    Built i first along j = 0, then j up at every i. `needed` (a set of
    (i, j, t)) prunes to what a consumer reads, walking the dependencies
    back; None keeps everything. E^{00}_0 = 1 is folded as a constant.
    """
    xa = f"p{DIRS[d]}a"   # X_PA for this direction
    xb = f"p{DIRS[d]}b"   # X_PB
    deps = {}
    order = []

    def rule(i, j, t):
        # returns list of (const, factor, (i,j,t)) for E^{ij}_t
        if j > 0:
            src, x = (i, j - 1), xb
        else:
            src, x = (i - 1, 0), xa
        si, sj = src
        terms = []
        if 0 <= t - 1 <= si + sj:
            terms.append((1.0, "oo2p", (si, sj, t - 1)))
        if 0 <= t <= si + sj:
            terms.append((1.0, x, (si, sj, t)))
        if 0 <= t + 1 <= si + sj:
            terms.append((float(t + 1), None, (si, sj, t + 1)))
        return terms

    for i in range(la + 1):
        for j in range(lb + 1):
            for t in range(i + j + 1):
                if i == 0 and j == 0:
                    continue
                deps[(i, j, t)] = rule(i, j, t)
                order.append((i, j, t))

    if needed is None:
        keep = set(order)
    else:
        keep = set()
        stack = [n for n in needed if n != (0, 0, 0)]
        while stack:
            n = stack.pop()
            if n in keep:
                continue
            keep.add(n)
            for _, _, src in deps[n]:
                if src != (0, 0, 0):
                    stack.append(src)

    stmts = []
    for n in order:
        if n not in keep:
            continue
        terms = []
        for c, f, src in deps[n]:
            if src == (0, 0, 0):
                terms.append(Term(c, f, None))       # E^{00}_0 = 1
            else:
                terms.append(Term(c, f, e_name(d, *src)))
        stmts.append(Stmt(e_name(d, *n), tuple(terms)))
    return stmts


def e_operand(d, i, j, t):
    """Name of E^{ij}_t, or None when it is the constant 1."""
    if i == 0 and j == 0 and t == 0:
        return None
    return e_name(d, i, j, t)
