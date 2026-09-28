# Vendored from libfint

`poly.py`, `derive.py` and `factor.py` are copied verbatim from libfint's
`scripts/rotaxis_mmd/` at commit `bb9c33a` ("Measure the gradient crossover:
it is contraction depth, not angular momentum"). They are the
target-agnostic half of libfint's rotated-axis McMurchie-Davidson generator:
the symbolic derivation, the component order and the three-level
factorisation. Nothing in them knows what language is emitted.

Everything terco-specific lives beside them in files libfint does not have.
libfint's own `emit.py` (CPU Fortran, per-quartet ABI) is deliberately not
vendored.

**Keep them verbatim.** To take an upstream change, copy the files again and
update the commit above. To check for drift:

    for f in poly derive factor; do
        diff scripts/rotaxis/$f.py ../libfint/scripts/rotaxis_mmd/$f.py
    done

A change terco needs in the algebra -- for example the
component order in `derive.cart_components` -- should go upstream first, or
be made by wrapping, not by editing these copies in place.
