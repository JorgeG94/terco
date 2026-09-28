#!/usr/bin/env python3
"""
Emit the rotated-axis MD kernels for the 16 ordered s/p classes: one module
per class (`trc_ra_kXXXX.F90`) plus the dispatcher (`trc_ra_kernels.F90`).

WHY A SEPARATE GENERATOR FROM gen_perclass.py
----------------------------------------------
The integral core is a different algorithm (rotated-axis MD, following
scripts/rotaxis/FRAME.md, rather than Obara-Saika VRR + HGP HRR), so it needs
its own derivation (scripts/rotaxis/{derive,factor}.py) and its own emission
(scripts/rotaxis/emit_docon.py). Everything BESIDES the integral core --
segment locate/decode, the exact Schwarz and density-weighted screens, the
same_ab/same_cd/same_pair logic, the scalar-path AO decode, the block-
accumulated density-batched digestion (and its TRC_NO_DIGEST branch), the
scalar launch driver and the item routine signature -- is the SAME text
gen_perclass.py already emits, called from the same functions; see
emit_docon.py's module docstring.

This generator is fixed at s/p (kind rank <= 1): FRAME.md section 5 notes d
is not implemented (libfint's own rotaxis path cuts d because it loses to
Rys there), so there is no `--lmax` here the way gen_perclass.py has one --
the class set is always the 16 ordered {s,p}^4 combinations, at every
TRC_LMAX >= 1 the rest of the library supports.

USAGE
    gen_rotaxis.py --split src/generated
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from rotaxis import emit_docon  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", metavar="DIR", required=True,
                     help="emit one module per class into DIR, plus "
                          "trc_ra_kernels.F90, matching gen_perclass.py "
                          "--split's own layout.")
    ap.add_argument("--scalarize", action="store_true",
                     help="name sacc/racc/kmv/cmv/gb/bn/bz as individual "
                          "scalars instead of small arrays -- the register-"
                          "promotion comparison of the milestone-1 review. "
                          "Every index into these is a compile-time literal "
                          "either way; this only changes whether ptxas is "
                          "TOLD to keep the array out of local memory or has "
                          "to work it out on its own. Measure before trusting "
                          "either form (see the review's ptxinfo table).")
    ap.add_argument("--bra-uniform", action="store_true",
                     help="use the mirrored (bra-uniform, ket fastest) "
                          "segment decode instead of HGP's own (bra fastest) "
                          "one -- coordinator review item B. NOT recommended: "
                          "measured net negative on a 20-water 6-31G cluster "
                          "(psps improved, psss/ppss/ppps got worse, some "
                          "substantially), see emit_class_module's comment "
                          "for the likely mechanism. Exists so the experiment "
                          "is reproducible, not because it should be turned "
                          "on.")
    ap.add_argument("--legacy-bra-inner", action="store_true",
                     help="derive Derivation(ka,kb,kc,kd) directly (physical "
                          "BRA primitives innermost) instead of the DEFAULT "
                          "Derivation(kc,kd,ka,kb) (physical KET primitives "
                          "innermost, matching HGP's own nesting -- FRAME.md "
                          "section 10, coordinator review). MEASURED, not "
                          "assumed: on a 20-water 6-31G cluster the 75-water "
                          "Fock build is 11.8% faster with the ket-inner "
                          "default than with this flag (1.14x HGP vs 1.29x), "
                          "reproducibly to within 0.1% over repeats, and no "
                          "class regressed by more than the fixed ~10% "
                          "decision-rule bound -- see emit_class_module's "
                          "comment and the commit message for the full "
                          "per-class table. This flag exists to reproduce "
                          "the earlier (milestone-1) default for comparison, "
                          "not because it measures better.")
    args = ap.parse_args()

    # The output directory is not recorded, as in gen_perclass.py: the header
    # names the options that shape the code, so regenerating anywhere else
    # reproduces the shipped files byte for byte.
    cmd = "scripts/gen_rotaxis.py --split"
    if args.scalarize:
        cmd += " --scalarize"
    if args.bra_uniform:
        cmd += " --bra-uniform"
    if args.legacy_bra_inner:
        cmd += " --legacy-bra-inner"
    os.makedirs(args.split, exist_ok=True)
    files = emit_docon.emit_files(cmd, scalarize=args.scalarize,
                                  bra_fastest=not args.bra_uniform,
                                  ket_inner=not args.legacy_bra_inner)
    for name, text in files.items():
        with open(os.path.join(args.split, name), "w") as fh:
            fh.write(text)
    print(f"wrote {len(files)} rotaxis files ({len(files) - 1} classes + "
          f"dispatcher) into {args.split}")


if __name__ == "__main__":
    main()
