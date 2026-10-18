#!/usr/bin/env python3
"""Compare one check's own numbers against its committed golden.

Standard library only. ctest runs this as terco/golden/<case>, after the
check that writes <out_dir>/<case>.out (test/golden_io.F90). See "Two
tiers" in test/pyscf_refs/README.md for when goldens and pyscf references
are each regenerated.

    python3 test/goldens/compare_goldens.py <case> <out_dir>
    python3 test/goldens/compare_goldens.py --report <out_dir>

--report prints every quantity of every case in <out_dir> against its
golden, then the largest |diff| per quantity type, and always exits 0. CI
runs it after each job's tests, so the spread between compilers is in
every job's log, not only a failing one's.

Every key in the output must have a golden, and must agree with it to the
tolerance of its quantity (TOLERANCES below). A golden key the output does
not have is reported, not failed: some cases only run at a higher
TERCO_LMAX. Exit 0 with "SKIP" when the check skipped itself.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: (key prefix, absolute tolerance, why), first match wins. Set from the
#: cross-compiler spread measured in CI -- see "Two tiers" in
#: test/pyscf_refs/README.md.
TOLERANCES = [
    ("E_", 1.0e-10, "total energies, Eh"),
    ("TR_D", 1.0e-10, "tr(D G) of a fixed-density Fock build, Eh"),
    ("NORM_", 1.0e-10, "Frobenius norm of a fixed-density Fock build"),
    ("WSUM_", 1.0e-10, "weighted element sum of a fixed-density Fock build"),
]


def tolerance(key: str) -> float:
    for prefix, tol, _ in TOLERANCES:
        if key.startswith(prefix):
            return tol
    raise SystemExit(f"compare_goldens: no tolerance for key {key!r}; add its prefix to TOLERANCES")


def read_values(path: Path) -> tuple[dict[str, float], bool]:
    """(values, complete) of a golden or an output file."""
    values: dict[str, float] = {}
    complete = False
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line == "END":
            complete = True
            continue
        key, val = line.split()
        if key in values:
            raise SystemExit(f"compare_goldens: {path.name}: duplicate key {key}")
        values[key] = float(val)
    return values, complete


def report(out_dir: Path) -> int:
    worst: dict[str, tuple[float, str]] = {}
    for out in sorted(out_dir.glob("*.out")):
        case = out.stem
        got, complete = read_values(out)
        gold = HERE / f"{case}.gold"
        if not (got and complete and gold.exists()):
            continue
        want, _ = read_values(gold)
        for key, val in got.items():
            if key not in want:
                continue
            prefix = next((p for p, _, _ in TOLERANCES if key.startswith(p)), "?")
            diff = abs(val - want[key])
            print(f"golden-spread {case:34s} {key:44s} {diff:9.2e}")
            if diff >= worst.get(prefix, (-1.0, ""))[0]:
                worst[prefix] = (diff, f"{case} {key}")
    for prefix, (diff, where) in sorted(worst.items()):
        print(f"golden-spread-max {prefix:6s} {diff:9.2e}  ({where})")
    return 0


def main() -> int:
    if len(sys.argv) == 3 and sys.argv[1] == "--report":
        return report(Path(sys.argv[2]))
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    case, out_dir = sys.argv[1], Path(sys.argv[2])
    out = out_dir / f"{case}.out"
    gold = HERE / f"{case}.gold"
    if not out.exists():
        print(f"SKIP: {out} not written (the check did not run, or skipped itself)")
        return 0
    got, complete = read_values(out)
    if not got and not complete:
        print(f"SKIP: the check skipped itself; {out.name} has no values")
        return 0
    if not complete:
        print(f"FAIL: {out} is incomplete (no END) -- the check stopped part way")
        return 1
    if not gold.exists():
        print(f"FAIL: no golden {gold.relative_to(HERE.parent.parent)} -- run "
              f"test/goldens/regen_goldens.py (see 'Two tiers' in test/pyscf_refs/README.md)")
        return 1
    want, _ = read_values(gold)
    nfail = 0
    print(f"{'key':48s} {'golden':>24s} {'this build':>24s} {'|diff|':>9s} {'tol':>8s}")
    for key, val in got.items():
        if key not in want:
            print(f"{key:48s} {'(none)':>24s} {val:24.16e}   FAIL: no golden for this key")
            nfail += 1
            continue
        tol = tolerance(key)
        diff = abs(val - want[key])
        ok = diff <= tol
        nfail += not ok
        print(f"{key:48s} {want[key]:24.16e} {val:24.16e} {diff:9.2e} {tol:8.0e}"
              f"{'' if ok else '   FAIL'}")
    for key in want:
        if key not in got:
            print(f"{key:48s} not run in this build")
    if nfail:
        print(f"compare_goldens {case}: {nfail} quantity(ies) off their golden -- FAIL")
        return 1
    print(f"compare_goldens {case}: {len(got)} quantities match their golden -- PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
