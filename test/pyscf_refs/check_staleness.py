#!/usr/bin/env python3
"""Cheap, pyscf-free staleness guard for every committed pyscf reference.

Standard library only (hashlib, pathlib, re, sys) -- this runs in CI, where
pyscf and numpy are never installed (see this directory's README.md). It
does not re-derive any number; it only recomputes each reference's own
"hash:"/"files:" header (test/pyscf_refs/generate.py's RefWriter.header()
writes both: the hash of exactly the files named on the "files:" line,
which always includes generate.py itself) and compares. A mismatch means
some input generate.py read -- geometry, basis JSON, the Fortran source
that builds a fixed density or grid, or the generator script -- changed
since the reference was written.

Usage:
    python3 test/pyscf_refs/check_staleness.py          # check every *.ref
                                                          # and reference_xc.dat
Exit code 0: every reference is fresh. Exit code 1: at least one is stale
or malformed, with a message naming it and the fix (rerun generate.py).

Registered as its own ctest (terco/check_pyscf_refs_fresh, test/CMakeLists.txt)
so a stale reference fails the suite instead of silently comparing against
a number that no longer describes its own inputs.
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent     # test/pyscf_refs
TEST = HERE.parent
ROOT = TEST.parent

HASH_RE = re.compile(r"^#\s*hash:\s*([0-9a-f]{64})\s*$")
FILES_RE = re.compile(r"^#\s*files:\s*(.*\S)\s*$")


def sha256_of(*paths: Path) -> str:
    h = hashlib.sha256()
    for p in paths:
        h.update(p.read_bytes())
    return h.hexdigest()


def check_one(path: Path) -> str | None:
    """Returns None if `path` is fresh, else an error message."""
    text = path.read_text(errors="replace")
    stored_hash = None
    files = None
    for line in text.splitlines():
        if not line.startswith("#"):
            # Headers are always the leading comment block; once a non-'#'
            # line shows up, both tags (if present) have already been seen.
            break
        m = HASH_RE.match(line)
        if m:
            stored_hash = m.group(1)
        m = FILES_RE.match(line)
        if m:
            files = m.group(1).split()
    if stored_hash is None or files is None:
        return f"{path.relative_to(ROOT)}: no 'hash:'/'files:' header -- not written by generate.py?"
    missing = [f for f in files if not (ROOT / f).exists()]
    if missing:
        return (f"{path.relative_to(ROOT)}: input file(s) no longer exist: {missing} "
                f"-- rerun test/pyscf_refs/generate.py")
    recomputed = sha256_of(*[ROOT / f for f in files])
    if recomputed != stored_hash:
        return (f"{path.relative_to(ROOT)}: STALE (inputs changed since this reference was "
                f"written: {', '.join(files)}) -- rerun test/pyscf_refs/generate.py")
    return None


#: df_vs_pyscf.ref (~18 MB) is local-only -- .gitignore excludes it, and
#: check_df_vs_pyscf.F90 SKIPs rather than fails when it is absent (see
#: that program's own comment). Optional here in exactly the same sense:
#: `sorted(HERE.glob("*.ref"))` below simply does not find it when it does
#: not exist, so nothing special is needed to SKIP it -- but its absence
#: is worth a line in the output, not silence, and if a developer HAS
#: generated it locally, it is checked like every other reference.
OPTIONAL_REFS = {"df_vs_pyscf.ref"}


def main() -> int:
    targets = sorted(HERE.glob("*.ref"))
    for name in sorted(OPTIONAL_REFS):
        if not (HERE / name).exists():
            print(f"  (absent, optional) {name}")
    xc_dat = TEST / "reference_xc.dat"
    if xc_dat.exists():
        targets.append(xc_dat)
    if not targets:
        print("check_staleness: no references found under test/pyscf_refs/ -- "
              "nothing to check (run test/pyscf_refs/generate.py first)")
        return 0
    bad = []
    for t in targets:
        err = check_one(t)
        status = "stale" if err else "fresh"
        print(f"  {status:6s} {t.relative_to(ROOT)}")
        if err:
            bad.append(err)
    if bad:
        print("\ncheck_staleness: FAILED")
        for b in bad:
            print("  " + b)
        return 1
    print(f"\ncheck_staleness: all {len(targets)} reference(s) fresh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
