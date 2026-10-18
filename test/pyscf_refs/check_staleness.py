#!/usr/bin/env python3
"""Cheap, pyscf-free staleness guard for every committed pyscf reference.

Standard library only (ast, hashlib, pathlib, re, sys) -- this runs in CI,
where pyscf and numpy are never installed (see this directory's README.md).
It does not re-derive any number; it only recomputes each reference's own
"hash:" header line and compares. generate.py writes that line with
refhash.case_hash, and this script recomputes it with the same function:
the hash of the reference's own CASES entry in generate.py (plus every
definition of generate.py that entry reaches, comments and docstrings
ignored) and of the files on its "files:" line. A mismatch means the code
that wrote this reference, or an input it read -- geometry, basis JSON, the
Fortran source that builds a fixed density or grid -- changed since it was
written. Editing a different case's generator does not restale it.

Usage:
    python3 test/pyscf_refs/check_staleness.py          # check every *.ref
                                                          # and reference_xc.dat
Exit code 0: every reference is fresh. Exit code 1: at least one is stale
or malformed, with a message naming it and the fix (rerun generate.py for
that case).

Registered as its own ctest (terco/check_pyscf_refs_fresh, test/CMakeLists.txt)
so a stale reference fails the suite instead of silently comparing against
a number that no longer describes its own inputs.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent     # test/pyscf_refs
TEST = HERE.parent
ROOT = TEST.parent

sys.path.insert(0, str(HERE))
import refhash  # noqa: E402

HASH_RE = re.compile(r"^#\s*hash:\s*([0-9a-f]{64})\s*$")
FILES_RE = re.compile(r"^#\s*files:(.*)$")
SCHEME_RE = re.compile(r"^#\s*hash-scheme:\s*(\S+)\s*$")
CASE_RE = re.compile(r'^#\s*generator:.*CASES\["([^"]+)"\]')


def check_one(path: Path, gen: refhash.Generator) -> str | None:
    """Returns None if `path` is fresh, else an error message."""
    rel = path.relative_to(ROOT)
    text = path.read_text(errors="replace")
    stored_hash = files = scheme = case = None
    for line in text.splitlines():
        if not line.startswith("#"):
            # Headers are always the leading comment block; once a non-'#'
            # line shows up, every tag (if present) has already been seen.
            break
        for regex, kind in ((HASH_RE, "hash"), (FILES_RE, "files"),
                            (SCHEME_RE, "scheme"), (CASE_RE, "case")):
            m = regex.match(line)
            if m:
                if kind == "hash":
                    stored_hash = m.group(1)
                elif kind == "files":
                    files = m.group(1).split()
                elif kind == "scheme":
                    scheme = m.group(1)
                else:
                    case = m.group(1)
    if stored_hash is None or files is None:
        return f"{rel}: no 'hash:'/'files:' header -- not written by generate.py?"
    if scheme != refhash.SCHEME or case is None:
        return (f"{rel}: header is not hash scheme {refhash.SCHEME} with a CASES key "
                f"-- rerun test/pyscf_refs/generate.py for it")
    if case not in gen.cases:
        return f"{rel}: generate.py has no CASES entry {case!r} any more"
    missing = [f for f in files if not (ROOT / f).exists()]
    if missing:
        return (f"{rel}: input file(s) no longer exist: {missing} "
                f"-- rerun test/pyscf_refs/generate.py {case}")
    recomputed = refhash.case_hash(gen, case, [ROOT / f for f in files])
    if recomputed != stored_hash:
        what = ", ".join(["generate.py " + gen.function(case), *files])
        return (f"{rel}: STALE (its code or inputs changed since it was written: {what}) "
                f"-- rerun test/pyscf_refs/generate.py {case}")
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
    gen = refhash.Generator.from_path(HERE / "generate.py")
    bad = []
    for t in targets:
        err = check_one(t, gen)
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
