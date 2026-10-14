#!/usr/bin/env python3
"""Shared helper: which ExchCXX built-in kernels does the translator handle.

Both tools/gen_check_xc.py (the libxc comparison test's plumbing) and
tools/gen_xc_registry.py (the functional registry and docs/functionals.md)
need the same two facts -- which of ExchCXX's kernel headers
tools/exchcxx_to_fortran.py can translate, and why it refuses the rest --
so this asks the translator itself rather than keeping a second list that
could drift from the first.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import exchcxx_to_fortran as x2f  # noqa: E402

# Headers under kernels/ that are not themselves a kernel: the screening
# wrapper template exchcxx_to_fortran.py's docstring describes, emitted from
# a hand-written template rather than translated.
NOT_A_KERNEL = {"screening_interface"}

# Kernels the translator accepts -- the Maple arithmetic is straight-line
# and every function it calls is known -- but whose translated output does
# NOT agree with libxc at check_xc's tolerances (1e-11 v, 1e-10 fxc), on the
# FULL grid test/xc_ref.py probes, not only at a density floor the way
# hjs_pbe_x's relatives were first suspected to need (raising the floor a
# decade did not fix it; the disagreement is not confined to one corner).
# Measured worst relative errors (unpolarised / polarised, vxc / fxc):
#   ityh_x              2e-11 / 3e-10   9e-08 / 1e-04
#   ityh_x_015          9e-12 / 2e-10   2e-07 / 5e-05
#   ityh_x_033          5e-11 / 6e-10   8e-07 / 7e-04
#   hjs_pbe_x           8e-06 / 4e-10   7e-08 / 5e-08  (vxc fine, fxc not)
#   lrcwpbe_hjs_pbe_x   4e-09 / 8e-09   1e-06 / 1e-05
#   lrcwpbeh_hjs_pbe_x  9e-10 / 5e-09   9e-08 / 2e-06
#   wb97x_d3_xc         (one polarised vxc point, 1.66e-11 vs the 1e-11
#                        tolerance, at its own omega-attenuation switch --
#                        not used by any registered functional: wB97XD3
#                        is already `missing` for its dispersion term)
# terco's own four evaluation routines agree with each other exactly
# (the "self" column is 0) for every one of these, so the discrepancy is
# between ExchCXX's Maple derivation and libxc's C, not a translation slip.
# Not loosened: excluded from check_xc (which only asserts kernels that DO
# agree) and from the registry's notion of an available kernel, each with
# this reason. See docs/functionals.md.
DISAGREES_WITH_LIBXC = {
    "ityh_x": "translated kernel disagrees with libxc (fxc off by up to 1.3e-4)",
    "ityh_x_015": "translated kernel disagrees with libxc (fxc off by up to 4.6e-5)",
    "ityh_x_033": "translated kernel disagrees with libxc (fxc off by up to 6.6e-4)",
    "hjs_pbe_x": "translated kernel disagrees with libxc (fxc off by up to 5.3e-8)",
    "lrcwpbe_hjs_pbe_x": "translated kernel disagrees with libxc (fxc off by up to 1.4e-5)",
    "lrcwpbeh_hjs_pbe_x": "translated kernel disagrees with libxc (fxc off by up to 1.9e-6)",
    "wb97x_d3_xc": "translated kernel disagrees with libxc (vxc off by 1.66e-11 at its own "
                   "omega-attenuation switch, just past the 1e-11 tolerance)",
    "vwn3": "translated kernel disagrees with libxc on the GPU (polarised eps off by "
            "1.41e-11 at rho ~ 1e-8, zeta ~ 1; agrees on the CPU), and its UKS potential "
            "missed pyscf on terco's grid -- under investigation",
}


def discover(exchcxx: Path) -> tuple[list[tuple[str, str]], dict[str, str]]:
    """Return (translated [(stem, family), ...], refused {stem: reason}).

    A kernel is only really translated once `emit_module` has walked every
    statement of every one of its eight routines: `parse_kernel` alone does
    not call `translate_expr` on the body, so an unsupported function deep
    in a statement (the scaled exponential integral in the HSE/LC-wPBE
    kernels, say) would otherwise slip past as "translated" here and then
    fail for real when tools/exchcxx_to_fortran.py is actually run.
    """
    kdir = exchcxx / x2f.KERNEL_DIR
    translated = []
    refused = {}
    for path in sorted(kdir.glob("*.hpp")):
        stem = path.stem
        if stem in NOT_A_KERNEL:
            continue
        try:
            k = x2f.parse_kernel(path)
            x2f.emit_module(k, path.resolve())
            translated.append((stem, k.family))
        except x2f.Unsupported as e:
            refused[stem] = str(e)
    return translated, refused
