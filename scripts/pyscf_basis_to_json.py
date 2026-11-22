#!/usr/bin/env python3
"""pyscf's bundled basis data -> the MolSSI BSE "complete" JSON terco's
trc_basis_json.F90 reads, for a basis BSE itself does not have handy (or
that basis_set_exchange, which is not installed here, would otherwise
have produced). Every existing file under test/basis_sets/ came from a
real BSE download; this script exists only because that path is closed in
this environment and pyscf 2.14 carries the same data offline
(`pyscf.gto.basis.load`).

Schema (byte-identical keys/nesting to every committed basis_sets/*.json,
read one -- e.g. cc-pvdz.json -- before trusting this comment):

  { "molssi_bse_schema": {"schema_type": "complete", "schema_version": "0.1"},
    "revision_description": ..., "revision_date": ...,
    "elements": {"<Z>": {"electron_shells": [ {
        "function_type": "gto" | "gto_spherical",
        "region": "",
        "angular_momentum": [l],
        "exponents": ["1.234000E+01", ...],
        "coefficients": [["5.678000E-01", ...], ...]   # one row per column
      }, ... ], "references": [...] } } }

ONE JSON electron_shell PER ENTRY pyscf.gto.basis.load() returns, not
regrouped to match however the official BSE download happens to combine
or split a general contraction -- trc_basis_json.F90 reads "one shell per
coefficient column" regardless of how the columns are bundled into JSON
shell objects (see its own header), and a general contraction split into
two segmented pieces integrates to the exact same functions as the same
split bundled into one shell with the unused primitives' coefficients
zero-padded (a product against a zero is zero either way the terms are
grouped). Verified empirically for cc-pVDZ/O below, not just argued: see
this script's own round-trip check and the report it prints.

angular_momentum is always length 1 here: pyscf's own load() never
returns the length>1 "SP-combined" shape bse_json_basis's docstring
describes (that shape is specific to how some download formats, e.g.
6-31G's GAMESS text, bundle an s and a p on one exponent set; pyscf's
internal tables keep s and p as separate entries already).

function_type follows the existing files' own convention (cc-pvdz.json,
def2-svp.json): "gto" for s/p (no cartesian/spherical distinction exists
below l=2), "gto_spherical" for l>=2, since cc-pVTZ and def2-TZVP are
spherical-harmonic families in their usual definition, same as cc-pVDZ
and def2-SVP already committed. This tag is read by trc_basis_from_json
to set b%spherical, but check_rys_vs_hgp/check_rys_grad_vs_pyscf-style
tests that call eri%build/eri%fock/eri%gradient directly never apply that
transform -- bas%nao is the CARTESIAN count regardless of the tag, which
is what lets those tests build the same pyscf cart=True comparison this
converter's own round-trip check below relies on.

Usage:
    python3 pyscf_basis_to_json.py <pyscf-basis-name> <out.json> \
        --elements H,C,N,O [--description TEXT]

Restrict --elements to what a test actually needs: these files are
committed, and a full-periodic-table dump would bloat the repo for
elements no test exercises.
"""
import argparse
import json
from datetime import date

import pyscf
from pyscf import gto
from pyscf.data.elements import charge


#: Digits after the decimal point of every number written. 6 matches the
#: BSE downloads' own style (cc-pvdz.json's "1.301000E+01"); --digits 17
#: writes pyscf's doubles exactly, which is what a comparison against pyscf
#: needs when it must see the identical exponents and coefficients
#: (cc-pvqz-pyscf17.json is such a file).
DIGITS = 6


def _fmt(x: float) -> str:
    """Scientific, upper-case E, DIGITS digits after the decimal point."""
    return "%.*E" % (DIGITS, x)


def pick_primitives(entries, nprim):
    """`nprim` uncontracted single-primitive shells of every angular
    momentum named in `nprim` ({l: n}): the n exponents of that l spread
    evenly (by rank, tightest first) over all the exponents pyscf lists for
    it, each with coefficient 1. A cut-down basis that still has every
    shell type, so a test can reach h functions without paying for a full
    cc-pV5Z. pyscf renormalises a lone primitive, as terco does."""
    out = []
    for l in sorted(nprim):
        exps = sorted({p[0] for e in entries if e[0] == l for p in e[1:]}, reverse=True)
        if not exps:
            continue
        n = min(nprim[l], len(exps))
        idx = sorted({round(i * (len(exps) - 1) / max(n - 1, 1)) for i in range(n)}) if n > 1 else [len(exps) // 2]
        for i in idx:
            out.append([l, [exps[i], 1.0]])
    return out


def shells_for_element(basis_name: str, symbol: str, nprim=None) -> list:
    entries = gto.basis.load(basis_name, symbol)
    if nprim:
        entries = pick_primitives(entries, nprim)
    shells = []
    for entry in entries:
        l = entry[0]
        prims = entry[1:]
        ncol = len(prims[0]) - 1
        exps = [p[0] for p in prims]
        cols = [[p[1 + c] for p in prims] for c in range(ncol)]
        shells.append({
            "function_type": "gto" if l < 2 else "gto_spherical",
            "region": "",
            "angular_momentum": [l],
            "exponents": [_fmt(e) for e in exps],
            "coefficients": [[_fmt(c) for c in col] for col in cols],
        })
    return shells


def build_json(basis_name: str, symbols: list, description: str, nprim=None) -> dict:
    elements = {}
    for sym in symbols:
        z = charge(sym)
        elements[str(z)] = {
            "electron_shells": shells_for_element(basis_name, sym, nprim),
            "references": [{
                "reference_description": description,
                "reference_keys": [],
            }],
        }
    return {
        "molssi_bse_schema": {"schema_type": "complete", "schema_version": "0.1"},
        "revision_description": description,
        "revision_date": date.today().isoformat(),
        "elements": elements,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("basis", help="pyscf basis name, e.g. cc-pvtz, def2-tzvp")
    ap.add_argument("out", help="output JSON path")
    ap.add_argument("--elements", required=True,
                     help="comma-separated element symbols, e.g. H,C,N,O")
    ap.add_argument("--description", default=None,
                     help="revision_description text (default: generated)")
    ap.add_argument("--digits", type=int, default=6,
                     help="digits after the decimal point (default 6; 17 is exact)")
    ap.add_argument("--primitives", default=None,
                     help="cut-down uncontracted basis: l:count pairs, e.g. s:3,p:2,d:2,f:2,g:1,h:1 "
                          "(exponents picked from the named pyscf basis)")
    args = ap.parse_args()
    global DIGITS
    DIGITS = args.digits
    nprim = None
    if args.primitives:
        nprim = {"spdfghik".index(item.split(":")[0]): int(item.split(":")[1])
                 for item in args.primitives.split(",")}

    symbols = [s.strip() for s in args.elements.split(",") if s.strip()]
    desc = args.description or (
        f"{args.basis}, converted from pyscf {pyscf.__version__}'s "
        f"bundled basis data (pyscf.gto.basis.load), for elements {','.join(symbols)}"
    )
    data = build_json(args.basis, symbols, desc, nprim)
    with open(args.out, "w") as f:
        json.dump(data, f, indent=4)
        f.write("\n")
    print(f"wrote {args.out}: {', '.join(symbols)}")


if __name__ == "__main__":
    main()
