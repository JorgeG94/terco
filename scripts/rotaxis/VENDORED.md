# Vendored from libfint

Four files are copied verbatim from libfint, branch
`generator/cartesian-order` at commit `f39b3cc` ("Fix check_numeric's
energy half, and run the generators in CI"). Their content there is the
same as at `e2c20cf` ("Rename the default order to 'lexicographic'").

That branch was force-pushed on 2026-09-30 to take libfint's
README_TERCO.md out of its history; the code did not change. terco's own
commit messages cite the pre-rewrite hashes, which map as follows:

| cited in terco | now |
|---|---|
| `c181abb` Define the Cartesian component order once | `35f8690` |
| `07cd704` Transcribe GAMESS's f and g orders | `607f482` |
| `f477d4b` Decode the GAMESS rescale | `c643701` |
| `d324ef7` Derive the GAMESS rescale from a closed form | `49eb0d9` |
| `59f9e38` Rename the default order to "lexicographic" | `e2c20cf` |

The guide libfint wrote for this work is kept in terco at
`docs/LIBFINT_GENERATORS.md`, since it is no longer in libfint's git.

| here | libfint |
|---|---|
| `scripts/cartorder.py` | `scripts/cartorder.py` |
| `scripts/rotaxis/poly.py` | `scripts/rotaxis_mmd/poly.py` |
| `scripts/rotaxis/derive.py` | `scripts/rotaxis_mmd/derive.py` |
| `scripts/rotaxis/factor.py` | `scripts/rotaxis_mmd/factor.py` |

`cartorder.py` is the one definition of the within-shell Cartesian order,
shared with libfint rather than restated: the lexicographic order
(libcint's, and terco's default), GAMESS's (read out
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
