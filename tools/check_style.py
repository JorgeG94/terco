#!/usr/bin/env python3
"""Style rules from AGENTS.md that a grep can enforce.

    python3 tools/check_style.py [paths...]   # default: src app test

Fails, listing file:line and the rule, on a use in a Fortran source
(.F90 .f90 .F .f) of:

  system_clock   time goes through pic_timer's `timer_type`.
  bare MPI       `use mpi`, `use mpi_f08`, `include 'mpif.h'` or a call to an
                 `MPI_*(` routine: MPI goes through pic-mpi (`pic_mpi_lib`:
                 `comm_t`, `allreduce`, `bcast`, ...). The constants it
                 re-exports (`MPI_SUM`) are not calls and are fine.
  Doxygen        `!>`, `!<`, `@param`, `@brief`, `@return`: docs are FORD only
                 (`!!` after the entity).
  literal kinds  `real(8)`, `real(kind=8)` and `d` exponents (`1.0d0`): kinds
                 come from pic_types (`dp`, `1.0_dp`).
  console output a bare `print` or `write(*, ...)` in src/ and app/: output goes
                 through pic_logger's `global_logger`. Tests are exempt.

Comments and string literals are ignored for every rule but Doxygen, which
reads the comment. `src/generated` is skipped, since generators own it.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

EXTS = {'.f90', '.F90', '.f', '.F'}
SKIP_DIRS = {'generated', 'build', '.git'}
I = re.IGNORECASE

# (rule name, regex, what to do instead). Applied to the code part of a line.
CODE_RULES = [
    ('system_clock', re.compile(r'\bsystem_clock\b', I), "time with pic_timer's timer_type"),
    ('bare MPI', re.compile(r"\buse\s*,?\s*(?:intrinsic\s*)?(?:::)?\s*mpi(?:_f08)?\b|\bMPI_[A-Za-z0-9_]+\s*\(", I),
     'go through pic-mpi (pic_mpi_lib), not the MPI library directly'),
    ('literal kind', re.compile(r'\breal\s*\(\s*(?:kind\s*=\s*)?[48]\s*\)|(?<![\w.])\d+\.?\d*[dD][+-]?\d+\b'),
     'use pic_types kinds (dp, 1.0_dp)'),
]
# Only under these top-level directories; tests print and write freely.
LIB_RULES = [
    ('console output', re.compile(r"^\s*print\s*[*'\"(]|\bwrite\s*\(\s*\*\s*,", I),
     "log through pic_logger's global_logger"),
]
INCLUDE_MPIF = re.compile(r"^\s*#?\s*include\s*[<'\"]\s*mpif", I)
DOXYGEN = re.compile(r'^[<>]|@(?:param|brief|return)\b|\\brief\b')


def split_line(line: str) -> tuple[str, str]:
    """(code, comment): the code with the inside of string literals blanked
    (the quotes stay, so `print '(a)'` still reads as a print of a string), and
    the text after the first `!` that is not inside a string literal."""
    out = []
    quote = None
    for i, ch in enumerate(line):
        if quote:
            if ch == quote:
                quote = None
                out.append(ch)
            else:
                out.append(' ')
        elif ch in '\'"':
            quote = ch
            out.append(ch)
        elif ch == '!':
            return ''.join(out), line[i + 1:]
        else:
            out.append(ch)
    return ''.join(out), ''


def sources(paths):
    for p in map(Path, paths):
        if p.is_file():
            yield p
            continue
        for f in sorted(p.rglob('*')):
            if f.suffix in EXTS and not (SKIP_DIRS & set(f.parts)):
                yield f


def check_line(path: Path, line: str) -> list[tuple[str, str]]:
    """(rule, advice) for every rule this line breaks."""
    code, comment = split_line(line)
    hits = [(name, advice) for name, rx, advice in CODE_RULES if rx.search(code)]
    # The file name is inside a string literal, which `code` blanks.
    if INCLUDE_MPIF.match(line):
        hits.append(('bare MPI', 'go through pic-mpi (pic_mpi_lib), not the MPI library directly'))
    if path.parts and path.parts[0] in ('src', 'app'):
        hits += [(name, advice) for name, rx, advice in LIB_RULES if rx.search(code)]
    # `!$acc` / `!$omp` directives and `!!` FORD blocks are not Doxygen.
    if DOXYGEN.search(comment):
        hits.append(('Doxygen', 'FORD only: a `!!` block after the entity, never `!>`, `!<` or @param'))
    return hits


def main(argv: list[str]) -> int:
    paths = argv or ['src', 'app', 'test']
    count = 0
    for f in sources(paths):
        for n, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            for name, advice in check_line(f, line):
                print(f'{f}:{n}: {name}: {line.strip()}')
                count += 1
    if count:
        print(f'check_style: {count} violation(s); see the header of tools/check_style.py for the rules')
        return 1
    print('check_style: clean')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
