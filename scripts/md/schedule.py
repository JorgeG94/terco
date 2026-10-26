"""Statement order for an R DAG that keeps few values live at once.

rtree.RDag.statements() emits level by level (m descending), which is a
valid topological order but holds every level-m value until the whole of
level m - 1 is done. When the consumer scatters each m = 0 target as soon
as it is formed, a depth-first order from the targets lets most
intermediates die early. `peak_live` measures either order.
"""

from . import rtree


def dfs_order(g, targets=None):
    """Statements in depth-first order from each target (targets in order)."""
    by = {st.lhs: st for st in g.statements()}
    done = set()
    out = []

    def name(n):
        r, m = n
        return rtree.leaf_name(m) if r == (0, 0, 0) else rtree.r_name(r, m)

    def visit(n):
        nm = name(n)
        if nm in done or n[0] == (0, 0, 0):
            return
        for p in g.node_parents(n):
            visit(p)
        done.add(nm)
        out.append(by[nm])

    for t in (targets or g.targets):
        visit(t)
    return out


def peak_live(stmts, final):
    """Largest number of computed values still needed later, for an order.

    `final` are names consumed right after being produced (the targets,
    scattered at once): they are live only for their own statement.
    """
    last = {}
    for k, st in enumerate(stmts):
        for t in st.terms:
            if t.operand is not None:
                last[t.operand] = k
    live, peak = set(), 0
    for k, st in enumerate(stmts):
        if st.lhs not in final:
            live.add(st.lhs)
        live = {v for v in live if last.get(v, -1) > k or v == st.lhs and st.lhs not in final}
        peak = max(peak, len(live))
    return peak
