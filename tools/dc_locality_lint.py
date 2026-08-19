#!/usr/bin/env python3
"""Require every `do concurrent` to name the scalars it assigns.

    python3 tools/dc_locality_lint.py [paths...]      # default: src

WHY THIS IS A LINT AND NOT A STYLE PREFERENCE
---------------------------------------------
`do concurrent` privatises a variable that is assigned before it is read in an
iteration, so a loop that omits `local(...)` is perfectly correct as written.
That is precisely what makes the omission dangerous: it is invisible until the
construct is translated.

`tools/dc_to_omp.py` rewrites these loops as `!$omp parallel do` (host) or
`!$omp target teams distribute parallel do` (offload), and OpenMP's default is
**shared**. A scalar the standard privatised for you silently becomes one
variable torn between every thread.

This is not hypothetical. `fock_all_nosym` had exactly one such loop, and after
conversion it disagreed with the folded digestion on a symmetric density --
where the two are the same operator by construction. Everything else passed.

So: name them. The cost is a few words per loop; the alternative is a race that
only appears in a build nobody runs by default.

WHAT IT CANNOT SEE
------------------
Assignments through a called procedure's `intent(out)` dummy, and anything
built by string manipulation. It is a text lint over one construct, not an
analyser -- it catches the case that has actually bitten, cheaply, in CI.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

DC_RE = re.compile(r'^(\s*)do concurrent\s*\(([^)]*)\)(.*)$', re.IGNORECASE)
OPEN_RE = re.compile(r'(do|do concurrent|if\b.*\bthen)\b', re.IGNORECASE)
CLOSE_RE = re.compile(r'end\s*(do|if)\b', re.IGNORECASE)
ASSIGN_RE = re.compile(r'^\s*([A-Za-z]\w*)\s*=[^=]')
DOVAR_RE = re.compile(r'^\s*do\s+([A-Za-z]\w*)\s*=', re.IGNORECASE)


def check_file(path: Path) -> list[str]:
    problems: list[str] = []
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        m = DC_RE.match(line)
        if not m:
            continue
        header, tail = m.group(2), m.group(3).lower()
        named = set()
        for kw in ("local", "local_init", "shared", "reduce"):
            for grp in re.findall(kw + r'\s*\(([^)]*)\)', tail):
                named |= {t.strip().lower().split(':')[-1] for t in grp.split(',')}
        indices = {v.lower() for v in re.findall(r'(\w+)\s*=', header)}

        depth, body = 1, []
        for j in range(i + 1, len(lines)):
            stripped = lines[j].split('!')[0].strip()
            if OPEN_RE.match(stripped):
                depth += 1
            if CLOSE_RE.match(stripped):
                depth -= 1
            if depth == 0:
                break
            body.append(lines[j])

        assigned = set()
        for b in body:
            code = b.split('!')[0]
            for a in ASSIGN_RE.findall(code):
                assigned.add(a.lower())
            for a in DOVAR_RE.findall(code):
                assigned.add(a.lower())

        missing = sorted(assigned - indices - named)
        if missing:
            problems.append(
                f"{path}:{i + 1}: `do concurrent` assigns "
                f"{', '.join(missing)} without naming them: add "
                f"local({', '.join(missing)})")
    return problems


def main(argv: list[str]) -> int:
    paths = [Path(p) for p in (argv[1:] or ["src"])]
    files: list[Path] = []
    for p in paths:
        files.extend(sorted(p.rglob("*.F90")) if p.is_dir() else [p])

    problems: list[str] = []
    for f in files:
        problems.extend(check_file(f))

    for p in problems:
        print(f"  {p}")
    print(f"\n  {len(files)} files checked, {len(problems)} problem(s)")
    if problems:
        print("  See the docstring in this file for why this is not optional.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
