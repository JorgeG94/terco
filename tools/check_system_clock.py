#!/usr/bin/env python3
"""Keep `system_clock` out of terco: time with pic_timer's `timer_type`.

    python3 tools/check_system_clock.py [paths...]   # default: src app test

AGENTS.md: time goes through pic_timer (`start`, `stop`, `get_elapsed_time`),
not `system_clock`. This fails, listing file:line, on any `system_clock` in a
Fortran source (.F90 .f90 .F .f). Comments and string literals are ignored;
`src/generated` is skipped, since generators own it.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

EXTS = {'.f90', '.F90', '.f', '.F'}
SKIP_DIRS = {'generated', 'build', '.git'}
USE_RE = re.compile(r'\bsystem_clock\b', re.IGNORECASE)


def code_part(line: str) -> str:
    """The line with string literals blanked and any trailing comment cut."""
    out = []
    quote = None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
            out.append(' ')
        elif ch in '\'"':
            quote = ch
            out.append(' ')
        elif ch == '!':
            break
        else:
            out.append(ch)
    return ''.join(out)


def sources(paths):
    for p in map(Path, paths):
        if p.is_file():
            yield p
            continue
        for f in sorted(p.rglob('*')):
            if f.suffix in EXTS and not (SKIP_DIRS & set(f.parts)):
                yield f


def main(argv: list[str]) -> int:
    paths = argv or ['src', 'app', 'test']
    hits = []
    for f in sources(paths):
        for n, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if USE_RE.search(code_part(line)):
                hits.append(f'{f}:{n}: {line.strip()}')
    for h in hits:
        print(h)
    if hits:
        print(f'check_system_clock: {len(hits)} use(s) of system_clock; time with pic_timer\'s timer_type instead')
        return 1
    print('check_system_clock: no system_clock')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
