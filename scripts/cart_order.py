#!/usr/bin/env python3
"""
The one definition of within-shell Cartesian component order.

WHY A SEPARATE MODULE
----------------------
gen_vrr.py, gen_hrr.py, gen_1e.py and gen_df.py each carried their own copy of
the same function -- `cart_list(lmax)` or `cart(l)` -- textually identical,
emitting (lx,ly,lz) in libcint's order: x power descending outer, y power
descending inner. Four copies of one definition is four places that could
silently drift, and it is exactly what stands in the way of a build-wide
switch to GAMESS order (d = XX YY ZZ XY XZ YZ): that switch needs ONE place
to flip, not four to find and keep in sync.

This module is that one place. The generators call into it; none of them
enumerates a shell's components itself any more.

`order` IS AN EXPLICIT PARAMETER, NEVER A MODULE-LEVEL GLOBAL. A future tool
that compares two orderings in the same process -- or simply two generators
run back to back with different settings -- must not be able to contaminate
each other through shared mutable state. Every entry point here takes the
order it means and nothing is remembered between calls except the read-only
cache below.

Only "libcint" is implemented. Adding "gamess" is one more table entry in
`_ORDERS`; nothing else changes.
"""

from functools import lru_cache


def _libcint(l):
    """(lx,ly,lz) for one shell of angular momentum l, libcint order: x
    power descending outer, y power descending inner (z is whatever is
    left)."""
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


#: One entry per supported order. This is a fixed table of pure functions,
#: not a "current order" -- it carries no state and adding "gamess" later is
#: adding a key here, nothing more.
_ORDERS = {
    "libcint": _libcint,
}


def components(l, order="libcint"):
    """(lx,ly,lz) tuples for one shell of angular momentum l, in `order`."""
    try:
        fn = _ORDERS[order]
    except KeyError:
        raise ValueError(f"unknown Cartesian order {order!r}; "
                          f"have {sorted(_ORDERS)}") from None
    return fn(l)


@lru_cache(maxsize=None)
def _table(l, order):
    """0-based position of each (lx,ly,lz) within its shell, in `order`. A
    LOOKUP built from components() -- never a second, closed-form encoding
    of the order."""
    return {p: i for i, p in enumerate(components(l, order))}


def index(lx, ly, lz, order="libcint"):
    """0-based position of (lx,ly,lz) within its shell, in `order`."""
    return _table(lx + ly + lz, order)[(lx, ly, lz)]
