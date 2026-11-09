#!/usr/bin/env python3
"""Write the basis sets terco ships (basis_sets/) from the Basis Set Exchange.

Every file is basis_set_exchange's own JSON export (get_basis(name,
fmt="json")), under its own file name (transform_basis_name: lower case, "*"
as "_st_", "/" as "_sl_"), which is the first spelling the CLI's -b lookup
tries (app_resolve_basis). Nothing is converted or trimmed, so the files are
reproducible from the pinned BSE version alone:

    python3 -m pip install basis_set_exchange==0.12
    tools/fetch_basis.py            # write basis_sets/
    tools/fetch_basis.py --check    # exit 1 if any file differs from BSE

test/basis_sets/ is separate on purpose: it holds the tests' frozen inputs,
some of them converted from pyscf's bundled data with a different contraction
layout, and it does not change with this list.
"""
import argparse
import os
import sys

BSE_VERSION = "0.12"

# The families v0.1.0 ships. QZ and the 6-311G(3df,3pd)-style sets reach g
# (or f) functions, so they need a TERCO_LMAX=4 (or 3) build; the CLI refuses
# them in a lower build with the ceiling in the message.
BASES = [
    # Pople
    "STO-3G", "3-21G",
    "6-31G", "6-31G*", "6-31G**", "6-31+G", "6-31+G*", "6-31+G**",
    "6-31++G", "6-31++G*", "6-31++G**", "6-31G(2df,p)",
    "6-311G", "6-311G*", "6-311G**", "6-311+G*", "6-311+G**",
    "6-311++G", "6-311++G**", "6-311++G(2d,2p)", "6-311++G(3df,3pd)",
    # Karlsruhe def2 (the ECPs ride inside each file; def2-ECP alone too)
    "def2-SV(P)", "def2-SVP", "def2-SVPD",
    "def2-TZVP", "def2-TZVPD", "def2-TZVPP", "def2-TZVPPD",
    "def2-QZVP", "def2-QZVPD", "def2-QZVPP", "def2-QZVPPD",
    "def2-ECP",
    # Dunning
    "cc-pVDZ", "cc-pVTZ", "cc-pVQZ",
    "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ",
    # Fitting sets for --df (J/K) and --rimp2 (RI)
    "def2-universal-JFIT", "def2-universal-JKFIT",
    "def2-SVP-RIFIT", "def2-SVPD-RIFIT", "def2-TZVP-RIFIT", "def2-TZVPD-RIFIT",
    "def2-TZVPP-RIFIT", "def2-TZVPPD-RIFIT", "def2-QZVPP-RIFIT", "def2-QZVPPD-RIFIT",
    "cc-pVDZ-RIFIT", "cc-pVTZ-RIFIT", "cc-pVQZ-RIFIT",
    "aug-cc-pVDZ-RIFIT", "aug-cc-pVTZ-RIFIT", "aug-cc-pVQZ-RIFIT",
    "cc-pVTZ-JKFIT", "cc-pVQZ-JKFIT",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare instead of writing; exit 1 on any difference")
    ap.add_argument("-o", "--out", default=os.path.join(os.path.dirname(__file__), "..", "basis_sets"))
    a = ap.parse_args()

    import basis_set_exchange as bse
    from basis_set_exchange.misc import transform_basis_name
    if bse.version() != BSE_VERSION:
        sys.exit(f"fetch_basis: basis_set_exchange {bse.version()}, this list is pinned to {BSE_VERSION}")

    out = os.path.normpath(a.out)
    os.makedirs(out, exist_ok=True)
    bad = []
    for name in BASES:
        text = bse.get_basis(name, fmt="json")
        path = os.path.join(out, transform_basis_name(name) + ".json")
        if a.check:
            try:
                with open(path) as f:
                    same = f.read() == text
            except FileNotFoundError:
                same = False
            if not same:
                bad.append(path)
        else:
            with open(path, "w") as f:
                f.write(text)
    if a.check:
        listed = {transform_basis_name(n) + ".json" for n in BASES}
        extra = sorted(set(f for f in os.listdir(out) if f.endswith(".json")) - listed)
        for p in bad:
            print(f"fetch_basis: {p} differs from basis_set_exchange {BSE_VERSION}")
        for f in extra:
            print(f"fetch_basis: {os.path.join(out, f)} is not in the list")
        sys.exit(1 if bad or extra else 0)
    print(f"fetch_basis: {len(BASES)} basis sets from basis_set_exchange {BSE_VERSION} in {out}")


if __name__ == "__main__":
    main()
