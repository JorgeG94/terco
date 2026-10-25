#!/usr/bin/env python3
"""Keep `system_clock` out of terco: time with pic_timer's `timer_type`.

    python3 tools/check_system_clock.py [paths...]          # check (default: src app test)
    python3 tools/check_system_clock.py --mark [paths...]   # tag every untagged use
    python3 tools/check_system_clock.py --list [paths...]   # print every use

A use of `system_clock` is allowed only while it carries a `TODO(fix)` comment
on the statement's line or on the comment line directly above it. The check
fails on any use that does not, so the existing ones are a shrinking backlog
and no new one can be added. Fixing one means replacing it with `timer_type`
(`start`, `stop`, `get_elapsed_time`) and deleting its marker; a marker with
no `system_clock` left under it is reported as stale.

`--mark` is how the backlog was created: it inserts
`! TODO(fix): system_clock -> pic_timer timer_type` above each untagged
statement, at the statement's indentation, and is idempotent.

It is a text scan over Fortran sources (.F90 .f90 .F .f): comments and string
literals are ignored, continuation lines are followed back to the start of the
statement. `src/generated` is skipped, since generators own it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

EXTS = {'.f90', '.F90', '.f', '.F'}
SKIP_DIRS = {'generated', 'build', '.git'}
USE_RE = re.compile(r'\bsystem_clock\b', re.IGNORECASE)
MARKER = 'TODO(fix)'
MARK_TEXT = '! TODO(fix): system_clock -> pic_timer timer_type'


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


def statement_start(lines, i):
    """Index of the first line of the statement that line i continues."""
    while i > 0 and code_part(lines[i - 1]).rstrip().endswith('&'):
        i -= 1
    return i


def is_marked(lines, i):
    if MARKER in lines[i]:
        return True
    s = statement_start(lines, i)
    return s > 0 and lines[s - 1].lstrip().startswith('!') and MARKER in lines[s - 1]


def scan(path):
    """(uses, stale): untagged-or-tagged use line indices, and orphan markers."""
    lines = path.read_text().split('\n')
    uses = [i for i, ln in enumerate(lines) if USE_RE.search(code_part(ln))]
    stale = []
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith('!') and MARKER in ln and 'system_clock' in ln:
            nxt = i + 1
            while nxt < len(lines) and lines[nxt].lstrip().startswith('!'):
                nxt += 1
            end = nxt
            while end < len(lines) and code_part(lines[end]).rstrip().endswith('&'):
                end += 1
            if not any(USE_RE.search(code_part(lines[j]))
                       for j in range(nxt, min(end + 1, len(lines)))):
                stale.append(i)
    return lines, uses, stale


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('paths', nargs='*', default=['src', 'app', 'test'])
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--mark', action='store_true', help='tag every untagged use')
    g.add_argument('--list', action='store_true', help='print every use')
    args = ap.parse_args()

    bad, total, tagged = [], 0, 0
    for path in sources(args.paths):
        lines, uses, stale = scan(path)
        todo = [i for i in uses if not is_marked(lines, i)]
        total += len(uses)
        tagged += len(uses) - len(todo)
        if args.list:
            for i in uses:
                flag = 'TODO' if is_marked(lines, i) else 'NEW '
                print(f'{flag} {path}:{i + 1}: {lines[i].strip()}')
        if args.mark and todo:
            done = set()
            for i in reversed(todo):
                s = statement_start(lines, i)
                if s in done:
                    continue
                done.add(s)
                indent = re.match(r'\s*', lines[s]).group(0)
                lines.insert(s, indent + MARK_TEXT)
            path.write_text('\n'.join(lines))
            print(f'marked {len(done)} statement(s) in {path}')
        elif not args.mark and not args.list:
            bad += [f'{path}:{i + 1}: untagged system_clock: {lines[i].strip()}' for i in todo]
            bad += [f'{path}:{i + 1}: stale {MARKER} marker, no system_clock follows' for i in stale]

    if args.mark or args.list:
        print(f'{total} use(s), {tagged} already tagged')
        return 0
    if bad:
        print('\n'.join(bad))
        print('\nTime with pic_timer (timer_type: start, stop, get_elapsed_time), not '
              'system_clock. Existing uses are tagged TODO(fix); new ones are not allowed.')
        return 1
    print(f'ok: no new system_clock ({total} tagged TODO(fix) remaining)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
