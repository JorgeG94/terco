#!/usr/bin/env python3
"""Regenerate test/goldens/*.gold from terco's own runs. Never runs pyscf.

    python3 test/goldens/regen_goldens.py --build build [--build build-lmax3]

Runs every golden producer in each build directory through ctest (which
points TERCO_GOLDEN_OUT at <build>/test/golden_out), then writes each
completed output as test/goldens/<case>.gold. A case the first build skips
(it needs a higher TERCO_LMAX) is taken from the next build that runs it;
a case two builds both run is taken from the first, and the largest
difference between them is printed.

Regenerate only for an intentional numerical change, and only once the
pyscf anchors still pass -- see "Two tiers" in test/pyscf_refs/README.md.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
from compare_goldens import read_values  # noqa: E402


def producers(build: Path) -> list[str]:
    """ctest names of every golden producer in `build`, from its golden tests."""
    out = subprocess.run(["ctest", "--test-dir", str(build), "-N", "-R", "^terco/golden/"],
                         check=True, capture_output=True, text=True).stdout
    return [line.split("terco/golden/", 1)[1].strip() for line in out.splitlines()
            if "terco/golden/" in line]


def run(build: Path, jobs: int) -> dict[str, tuple[dict[str, float], str]]:
    cases = producers(build)
    if not cases:
        raise SystemExit(f"{build}: no terco/golden/* tests -- configure with Python3 available")
    # The golden tests require their producers' fixtures, so running the
    # golden tests runs the producers first; the comparisons themselves may
    # fail here (the goldens are what is being replaced), which is fine.
    subprocess.run(["ctest", "--test-dir", str(build), "-R", "^terco/golden/", "--parallel", str(jobs)],
                   capture_output=True, text=True)
    out_dir = build / "test" / "golden_out"
    found = {}
    for case in cases:
        f = out_dir / f"{case}.out"
        if not f.exists():
            continue
        values, complete = read_values(f)
        if complete and values:
            found[case] = (values, f.read_text().splitlines()[0])
    return found


def compiler_of(build: Path) -> str:
    cache = (build / "CMakeCache.txt").read_text().splitlines()
    fc = next((l.split("=", 1)[1] for l in cache if l.startswith("CMAKE_Fortran_COMPILER:")), "?")
    return Path(fc).name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="append", required=True, type=Path,
                    help="a configured and built terco build directory (repeatable)")
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    taken: dict[str, tuple[dict[str, float], str]] = {}
    for build in a.build:
        build = build.resolve()
        for case, (values, _) in sorted(run(build, a.jobs).items()):
            if case in taken:
                first, src = taken[case]
                common = set(first) & set(values)
                worst = max((abs(first[k] - values[k]) for k in common), default=0.0)
                print(f"  {case:32s} also in {build.name}: max |diff| vs {src} {worst:.2e}")
                continue
            taken[case] = (values, build.name)
            src = f"{compiler_of(build)} build {build.name} at {commit}"
            lines = ["# terco golden v1", f"# case: {case}",
                     f"# from: {src}, by test/goldens/regen_goldens.py"]
            lines += [f"{k} {v:.16e}" for k, v in values.items()]
            lines.append("END")
            (HERE / f"{case}.gold").write_text("\n".join(lines) + "\n")
            print(f"  wrote test/goldens/{case}.gold  ({len(values)} quantities, {build.name})")
    print(f"regen_goldens: {len(taken)} golden(s) written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
