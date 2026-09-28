# Vendored from libfint

Four files are copied verbatim from libfint, branch
`generator/cartesian-order` at commit `d324ef7` ("Derive the GAMESS rescale
from a closed form, and assert the tables against it"):

| here | libfint |
|---|---|
| `scripts/cartorder.py` | `scripts/cartorder.py` |
| `scripts/rotaxis/poly.py` | `scripts/rotaxis_mmd/poly.py` |
| `scripts/rotaxis/derive.py` | `scripts/rotaxis_mmd/derive.py` |
| `scripts/rotaxis/factor.py` | `scripts/rotaxis_mmd/factor.py` |

`cartorder.py` is the one definition of the within-shell Cartesian order,
shared with libfint rather than restated: libcint's order, GAMESS's (read out
of GAMESS's own int2a.src tables, through g), and `gamess_norm_sq`, the
closed form for GAMESS's per-component normalisation. terco's own generators
import it too, so one module decides the order for every kernel terco emits.
It sits in `scripts/` because `derive.py` imports it by absolute name.

The other three are the target-agnostic half of libfint's rotated-axis
McMurchie-Davidson generator: the symbolic derivation, the component order
and the three-level factorisation. Nothing in them knows what language is
emitted. libfint's own `emit.py` (CPU Fortran, per-quartet ABI) is
deliberately not vendored; terco's is `emit_docon.py`.

**Keep them verbatim.** To take an upstream change, copy the files again and
update the commit above. To check for drift:

    B=generator/cartesian-order   # or main, once merged
    git -C ../libfint show $B:scripts/cartorder.py | diff - scripts/cartorder.py
    for f in poly derive factor; do
        git -C ../libfint show $B:scripts/rotaxis_mmd/$f.py \
            | diff - scripts/rotaxis/$f.py
    done

A change terco needs in any of them should go upstream first, not be made
in these copies.
