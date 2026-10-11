"""AO permutation between LEXICOGRAPHIC order (libcint's, and what pyscf
always speaks -- it has no GAMESS mode) and terco's build's own order.

The Fortran-side counterpart is src/trc_cart_perm.F90, generated from the
same cartorder.py tables; this module is the Python side, for the *_ref.py
scripts that compare a terco probe file against a live pyscf array. Neither
reimplements the permutation independently of cartorder.py -- both are
lookups into it -- so the two cannot silently disagree.

The build's order is never guessed from a probe file: it is passed in
explicitly, because ctest sets the TERCO_CART_ORDER environment variable
from the CMake cache entry of the same name (see test/CMakeLists.txt), and
"no guessing from file contents" is the point -- a probe file's own bytes
carry no order marker, deliberately, since every terco writer/reader on the
Fortran side already agrees on the order by construction.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import cartorder  # noqa: E402


def build_order():
    """This build's Cartesian order, from the TERCO_CART_ORDER environment
    variable ctest sets. Defaults to "lexicographic" only when the variable
    is entirely unset (a script run outside ctest, e.g. by hand against a
    default build) -- an empty or unrecognised value is an error, not a
    silent fallback, since a wrong guess here would make a real ordering bug
    look like a passing test."""
    v = os.environ.get("TERCO_CART_ORDER")
    if v is None:
        return "lexicographic"
    v = v.strip()
    if v in ("", "libcint"):
        return "lexicographic"
    if v not in ("lexicographic", "gamess"):
        raise ValueError(f"TERCO_CART_ORDER={v!r} is not lexicographic or "
                         "gamess")
    return v


def ao_permutation(sh_l, order=None):
    """The AO permutation for a basis of shells with angular momenta `sh_l`
    (in shell order), between LEXICOGRAPHIC order and `order` (default:
    `build_order()`).

    Returns an int array `perm` of length nao with `perm[i] = j` meaning:
    this build's AO slot `i` is the SAME physical function as lexicographic
    slot `j`. Apply as:

      build_vec  = lex_vec[perm]                        # a vector
      build_mat  = lex_mat[np.ix_(perm, perm)]           # a matrix, both AO axes
      build_mat3 = lex_mat3[:, perm][:, :, perm]         # AO axes 1, 2 of a
                                                          # (k, nao, nao) tensor

    and, since `perm` is a permutation of range(nao), its own argsort is the
    inverse (lexicographic order from build order):

      inv = np.argsort(perm)
      lex_vec = build_vec[inv]
    """
    if order is None:
        order = build_order()
    perm = []
    off = 0
    for l in sh_l:
        lex = cartorder.components(l, "lexicographic")
        build = cartorder.components(l, order)
        lex_index = {c: i for i, c in enumerate(lex)}
        for c in build:
            perm.append(off + lex_index[c])
        off += len(lex)
    return perm
