#!/usr/bin/env python3
"""One-time migration of every committed reference to the per-case hash.

Until this change each reference's "hash:" line covered the whole of
generate.py, so editing one generator restaled all of them. refhash.py now
hashes only the code that wrote each reference. The numbers themselves were
produced by the current generators, so they do not need regenerating --
only the header does. This script:

  1. checks every reference is fresh under the OLD scheme (sha256 of the
     bytes of the files on its "files:" line, generate.py read from git at
     --old-rev) and refuses to touch any that is not;
  2. rewrites its generator/hash/files header lines to exactly what
     generate.py now writes (refhash.case_hash, generate.py dropped from
     "files:"), leaving every other line byte-identical.

Standard library only; pyscf is never imported or run.

    python3 test/pyscf_refs/migrate_ref_headers.py [--old-rev ffcf3a6]
"""
from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEST = HERE.parent
ROOT = TEST.parent
sys.path.insert(0, str(HERE))
import refhash  # noqa: E402

GENERATE_REL = "test/pyscf_refs/generate.py"

#: References whose CASES key is not their file stem.
CASE_OF = {
    "rimp2_ccpvdz.ref": "rimp2",
    "pcm_cavity_water.ref": "pcm_cavity",
    "reference_xc.dat": "xc_kernels",
    "cphf_water_6-31gs.ref": "cphf",
    "cphf_water_cc-pvdz.ref": "cphf",
    "cphf_formamide_6-31gs.ref": "cphf",
    "cphf_formamide_cc-pvdz.ref": "cphf",
    "cphf_uhf_oh_6-31gs.ref": "cphf_uhf",
    "full_grad_ccpvdz.ref": "full_grad",
}


def case_of(path: Path, gen: refhash.Generator) -> str:
    case = CASE_OF.get(path.name, path.stem)
    if case not in gen.cases:
        raise SystemExit(f"{path.name}: no CASES entry {case!r}")
    return case


def migrate(path: Path, gen: refhash.Generator, old_generate: bytes) -> bool:
    lines = path.read_text().split("\n")
    head = 0
    while head < len(lines) and lines[head].startswith("#"):
        head += 1
    idx = {}
    for i in range(head):
        for tag in ("generator", "hash", "files"):
            if re.match(rf"^#\s*{tag}:", lines[i]):
                idx[tag] = i
    if "hash-scheme" in "".join(lines[:head]):
        print(f"  already  {path.relative_to(ROOT)}")
        return False
    if "hash" not in idx or "files" not in idx:
        raise SystemExit(f"{path.name}: no hash/files header")
    files = lines[idx["files"]].split(":", 1)[1].split()
    h = hashlib.sha256()
    for f in files:
        h.update(old_generate if f == GENERATE_REL else (ROOT / f).read_bytes())
    stored = lines[idx["hash"]].split(":", 1)[1].strip()
    if h.hexdigest() != stored:
        raise SystemExit(f"{path.name}: stale under the old scheme -- regenerate it instead")
    case = case_of(path, gen)
    inputs = [f for f in files if f != GENERATE_REL]
    new_hash = refhash.case_hash(gen, case, [ROOT / f for f in inputs])
    generator = (f'# generator: {GENERATE_REL} CASES["{case}"] -> {gen.function(case)}')
    hashing = ["# hash-scheme: " + refhash.SCHEME, "# hash: " + new_hash,
               "# files: " + " ".join(inputs)]
    # A .ref has a generator line, which keeps its place; reference_xc.dat
    # has none, and gets it just ahead of the hash, as generate.py writes it.
    out = []
    for i, line in enumerate(lines):
        if i == idx.get("generator"):
            out.append(generator)
        elif i == idx["hash"]:
            if "generator" not in idx:
                out.append(generator)
            out.extend(hashing)
        elif i == idx["files"]:
            continue
        else:
            out.append(line)
    path.write_text("\n".join(out))
    print(f"  migrated {path.relative_to(ROOT)}  ({case})")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-rev", default="HEAD",
                    help="git revision whose generate.py wrote the references (default HEAD)")
    a = ap.parse_args()
    old_generate = subprocess.run(["git", "show", f"{a.old_rev}:{GENERATE_REL}"], cwd=ROOT,
                                  check=True, capture_output=True).stdout
    gen = refhash.Generator.from_path(HERE / "generate.py")
    targets = sorted(HERE.glob("*.ref")) + [TEST / "reference_xc.dat"]
    n = sum(migrate(t, gen, old_generate) for t in targets)
    print(f"migrate_ref_headers: {n} of {len(targets)} reference(s) migrated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
