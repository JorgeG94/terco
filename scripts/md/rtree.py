"""The Hermite Coulomb integrals R^(m)_tuv as a pruned DAG, and its batches.

For a primitive bra pair P (exponent p) and ket pair Q (exponent q),
alpha = pq/(p + q), the two-centre Hermite integrals are

    R^(m)_000  = (-2 alpha)^m F_m(alpha |PQ|^2)
    R^(m)_r    = (r_d - 1) R^(m+1)_{r - 2 e_d} + X_PQ,d R^(m+1)_{r - e_d}

for any direction d with r_d > 0, and [p|q] = (-1)^|q| R^(0)_{p+q}. A class
whose bra and ket reach Hermite degrees lp and lq needs every R^(0)_r with
|r| <= L = lp + lq (Palethorpe and Barca, JCTC 21, 9388 (2025), sec 3.2).

Which direction each node reduces along is free, and it decides how many
intermediates exist. This module follows that paper:

  * MILD: reduce along the smallest non-zero component of r (Ryu et al.),
    ties broken round-robin. A node on a face or an edge of the degree
    tetrahedron then steps to the nearest edge, and a unit component costs
    one multiply and one intermediate.
  * Gill's sieve on the resulting DAG: a node all of whose consumers can be
    formed along another direction from nodes still present is deleted and
    its consumers rerouted, to a fixed point. It never adds a node.
  * Batching (their Algorithm 4): the targets are split by the direction of
    their largest component, then by the range of that component's value,
    into 1, 3, 6, 9, ... batches. Each batch gets its own DAG; only the
    fundamental (000)^(m) and some axial nodes are shared between batches,
    so the duplicated work is small while each kernel's live set shrinks.

The statement form is `hermite.Stmt`, so the reference can execute exactly
what the Fortran emitter writes.
"""

from collections import deque

from .hermite import DIRS, Stmt, Term, herm


def r_name(r, m):
    t, u, v = r
    return f"r{t}_{u}_{v}_{m}"


def leaf_name(m):
    """R^(m)_000: the prefactor times (-2 alpha)^m F_m, computed by the caller."""
    return f"rb{m}"


def _sub(r, d, k):
    s = list(r)
    s[d] -= k
    return tuple(s)


def parents(r, m, d):
    """Nodes R^(m)_r reads when reduced along d."""
    if r[d] < 1:
        raise ValueError(f"cannot reduce {r} along {DIRS[d]}")
    ps = [(_sub(r, d, 1), m + 1)]
    if r[d] >= 2:
        ps.append((_sub(r, d, 2), m + 1))
    return ps


def mild_dir(r):
    """The minimum non-zero component's direction, round-robin on ties.

    The rotation is keyed on |r| so that ties at one degree do not all fall
    the same way; any fixed rule is correct, this one only spreads them.
    """
    nz = [d for d in range(3) if r[d] > 0]
    lo = min(r[d] for d in nz)
    cands = [d for d in nz if r[d] == lo]
    s = sum(r) % 3
    for k in range(3):
        d = (s + k) % 3
        if d in cands:
            return d
    raise AssertionError


class RDag:
    """Nodes (r, m) with a reduction direction each; leaves are r = 000."""

    def __init__(self, targets):
        self.targets = [(r, 0) for r in targets]
        self.tset = set(self.targets)
        self.dir = {}
        self.nodes = set()
        q = deque(self.targets)
        while q:
            n = q.popleft()
            if n in self.nodes:
                continue
            self.nodes.add(n)
            r, m = n
            if r == (0, 0, 0):
                continue
            d = mild_dir(r)
            self.dir[n] = d
            for p in parents(r, m, d):
                if p not in self.nodes:
                    q.append(p)

    def node_parents(self, n):
        r, m = n
        if r == (0, 0, 0):
            return []
        return parents(r, m, self.dir[n])

    def consumers(self):
        cons = {n: [] for n in self.nodes}
        for n in self.nodes:
            for p in self.node_parents(n):
                cons[p].append(n)
        return cons

    def sieve(self):
        """Gill's sieve to a fixed point; returns the number of nodes removed.

        Removing a node can make another removable, so a single pass leaves
        most of the win behind (gen_vrr.py's experience, same algorithm).
        The visiting order is breadth-first from the targets, as in the
        paper's Algorithm 3; it moves the result by a percent or two only.
        """
        removed = 0
        changed = True
        while changed:
            changed = False
            cons = self.consumers()
            order = sorted(self.nodes - self.tset,
                           key=lambda n: (n[1], -sum(n[0]), n[0]))
            for i in order:
                if i not in self.nodes:
                    continue
                cs = cons[i]
                if not cs:
                    self._drop(i, cons)
                    removed += 1
                    changed = True
                    continue
                reroute = {}
                ok = True
                for c in cs:
                    r, m = c
                    alt = None
                    for d in range(3):
                        if d == self.dir[c] or r[d] == 0:
                            continue
                        ps = parents(r, m, d)
                        if all(p in self.nodes and p != i for p in ps):
                            alt = d
                            break
                    if alt is None:
                        ok = False
                        break
                    reroute[c] = alt
                if not ok:
                    continue
                for c, d in reroute.items():
                    for p in self.node_parents(c):
                        cons[p].remove(c)
                    self.dir[c] = d
                    for p in self.node_parents(c):
                        cons[p].append(c)
                self._drop(i, cons)
                removed += 1
                changed = True
        return removed

    def _drop(self, i, cons):
        for p in self.node_parents(i):
            cons[p].remove(i)
        self.nodes.discard(i)
        self.dir.pop(i, None)
        cons.pop(i, None)

    def counts(self):
        """(intermediates, flops): nodes with m > 0, and arithmetic per R set.

        A reduction with r_d = 1 is one multiply; with r_d >= 2 a multiply
        by the integer, a multiply by X_PQ and an add. Leaves cost one
        multiply each (the prefactor power times F_m).
        """
        inter = sum(1 for n in self.nodes if n[1] > 0)
        flops = 0
        for n in self.nodes:
            r, m = n
            if r == (0, 0, 0):
                flops += 1
            elif r[self.dir[n]] == 1:
                flops += 1
            else:
                flops += 3
        return inter, flops

    def mmax(self):
        return max(m for _, m in self.nodes)

    def statements(self):
        """Straight-line code, highest m first; level-0 nodes are targets.

        Every level m reads only level m + 1, so emitting by level descending
        is a topological order, and level m + 1 is dead once level m is done.
        Leaves are names the caller fills (`leaf_name`).
        """
        out = []
        for m in range(self.mmax(), -1, -1):
            lvl = sorted((n for n in self.nodes if n[1] == m and n[0] != (0, 0, 0)),
                         key=lambda n: (sum(n[0]), n[0]))
            for n in lvl:
                r, _ = n
                d = self.dir[n]
                terms = []
                for k, (pr, pm) in enumerate(parents(r, m, d)):
                    name = leaf_name(pm) if pr == (0, 0, 0) else r_name(pr, pm)
                    if k == 0:
                        terms.append(Term(1.0, f"pq{DIRS[d]}", name))
                    else:
                        terms.append(Term(float(r[d] - 1), None, name))
                out.append(Stmt(r_name(r, m), tuple(terms)))
        return out

    def leaves(self):
        return sorted(m for r, m in self.nodes if r == (0, 0, 0))


def _group_dir(r):
    """Direction of the largest component; rx = ry = rz goes to x.

    Ties between two components are broken by the reverse of mild_dir's
    rotation, as the paper does, so the X, Y and Z groups come out the
    same size.
    """
    hi = max(r)
    cands = [d for d in range(3) if r[d] == hi]
    if len(cands) == 3:
        return 0
    s = sum(r) % 3
    for k in range(3):
        d = (s - k) % 3
        if d in cands:
            return d
    raise AssertionError


def batches(L, nbatch):
    """Split {r : |r| <= L} into `nbatch` target sets (1 or a multiple of 3).

    A batch is a direction D and a range [a, b] of the largest component,
    D_{a,b}; the ranges of one direction are chosen so the largest batch is
    as small as possible (their stated criterion: each batch then digests a
    similar number of targets).
    """
    allr = list(herm(L))
    if nbatch == 1:
        return [allr]
    if nbatch % 3:
        raise ValueError("batch count must be 1 or a multiple of 3")
    nsub = nbatch // 3
    groups = [[], [], []]
    for r in allr:
        groups[_group_dir(r)].append(r)
    out = []
    for d in range(3):
        g = groups[d]
        if not g:
            continue
        hist = [0] * (L + 1)
        for r in g:
            hist[max(r)] += 1
        cuts = _balanced_cuts(hist, nsub)
        lo = 0
        for hi in cuts + [L]:
            b = [r for r in g if lo <= max(r) <= hi]
            if b:
                out.append(b)
            lo = hi + 1
    return out


def _balanced_cuts(hist, nsub):
    """Upper ends of nsub - 1 contiguous ranges minimising the largest sum."""
    n = len(hist)
    nsub = min(nsub, n)
    if nsub <= 1:
        return []
    best = None
    best_cuts = None

    def rec(start, left, cuts, worst):
        nonlocal best, best_cuts
        if left == 1:
            w = max(worst, sum(hist[start:]))
            if best is None or w < best:
                best, best_cuts = w, list(cuts)
            return
        for hi in range(start, n - left + 1):
            s = sum(hist[start:hi + 1])
            if best is not None and max(worst, s) >= best:
                continue
            rec(hi + 1, left - 1, cuts + [hi], max(worst, s))

    rec(0, nsub, [], 0)
    return best_cuts


def batch_dags(L, nbatch, sieve=True):
    out = []
    for tg in batches(L, nbatch):
        g = RDag(tg)
        if sieve:
            g.sieve()
        out.append(g)
    return out
