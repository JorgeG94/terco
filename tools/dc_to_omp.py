#!/usr/bin/env python3
"""Rewrite `do concurrent` blocks as `!$omp target teams distribute parallel do`.

    do concurrent (j=1:ny, i=1:nx) local(a, b)
        ...
    end do

becomes

    !$omp target teams distribute parallel do collapse(2) private(a, b)
    do j = 1, ny
    do i = 1, nx
        ...
    end do
    end do
    !$omp end target teams distribute parallel do

Rules:
  - Preserve header order: leftmost induction → outermost loop
  - `local(...)` → `private(...)` (only clause present in this codebase)
  - Drop `collapse(N)` when N == 1
  - Match the `end do` of the block by depth-tracking forward
  - Skip occurrences inside comments / string literals
  - Join `&` continuations on the header
  - Strip the `pure` attribute from any procedure that ends up containing an
    `!$omp target` region, or an `atomic` directive.  An executable target region launches a device
    kernel + moves data — side effects a `pure` procedure may not have, so
    NVHPC/ifx reject `pure` + `!$omp target`.  (`do concurrent` itself is
    side-effect-free and is fine in a pure procedure, which is why the
    OpenACC/`-stdpar` source keeps `pure`; only this OpenMP-target variant
    needs it removed.)  Declarative directives like `!$omp declare target`
    are pure-compatible and are left alone.

Targets (--target):
  gpu (default)  `!$omp target teams distribute parallel do` (offload).
  cpu            `!$omp parallel do` (host multicore) — drops `target teams
                 distribute`, whose host-fallback path is flaky on some
                 compilers (notably GNU). The `pure`-strip cascade keys on the
                 host `!$omp parallel do` opening directive instead of
                 `!$omp target` (a host parallel region is likewise rejected
                 inside a PURE procedure). Build with
                 -DRAKALI_PARALLEL_BACKEND=openmp -DRAKALI_ENABLE_GPU=OFF.

Usage:
    python tools/dc_to_omp.py                  # dry-run, scans src/ (gpu)
    python tools/dc_to_omp.py --write          # apply in place (gpu)
    python tools/dc_to_omp.py --target cpu --write src   # host multicore
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from functools import partial
from pathlib import Path

from _parallel import pmap, resolve_jobs
from _progress import Progress, track

DC_RE = re.compile(r"\bdo\s+concurrent\b", re.IGNORECASE)
DO_RE = re.compile(r"^\s*(?:[a-zA-Z_][a-zA-Z0-9_]*\s*:\s*)?do\b", re.IGNORECASE)
END_DO_RE = re.compile(r"^\s*end\s*do\b", re.IGNORECASE)
KNOWN_CLAUSES = ("local_init", "local", "shared", "default", "reduce")

# An *opening* `!$omp target` directive (not `!$omp end target`, and not a
# declarative `!$omp declare target`).
OMP_TARGET_RE = re.compile(r"^\s*!\$omp\s+target\b", re.IGNORECASE)
# The CPU-target opening worksharing directive (`!$omp parallel do`, not its
# `!$omp end parallel do`). Like a target region, a host parallel region is a
# side-effecting construct that the compiler rejects inside a PURE procedure.
OMP_PARALLEL_DO_RE = re.compile(r"^\s*!\$omp\s+parallel\s+do\b", re.IGNORECASE)

#
# An ATOMIC also costs a procedure its purity, and for a reason unrelated to
# launching anything:
#
#     Error: OpenMP directives other than SIMD or DECLARE TARGET at (1)
#            may not appear in PURE procedures
#
# OpenMP forbids an atomic inside a PURE procedure; OpenACC permits it. So a
# codebase whose scatter lives in a pure routine -- which is what
# `do concurrent` requires of anything it calls -- has to shed `pure` there
# too before the atomics can be translated. Matches either sentinel, since
# this may run before or after acc_to_omp.py.
#
OMP_ATOMIC_RE = re.compile(r"^\s*!\$(omp|acc)\s+atomic\b", re.IGNORECASE)
# A procedure-opening statement (subroutine/function with a name), e.g.
# `pure subroutine foo(`, `pure real(wp) function bar(`.
PROC_OPEN_RE = re.compile(
    r"^\s*[^!]*\b(?:subroutine|function)\b\s+[a-zA-Z]\w*", re.IGNORECASE
)
PROC_END_RE = re.compile(r"^\s*end\s*(?:subroutine|function)\b", re.IGNORECASE)
# The `pure` prefix token plus trailing blanks (leaves indentation intact;
# `\bpure\b` does not match inside `impure`).
PURE_RE = re.compile(r"\bpure\b[ \t]*", re.IGNORECASE)


def strip_strings(line: str) -> str:
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(" " if ch != quote else ch)
            if ch == quote:
                quote = None
        else:
            if ch in ("'", '"'):
                quote = ch
            out.append(ch)
    return "".join(out)


def code_part(line: str) -> str:
    s = strip_strings(line)
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def split_top_level(s: str, sep: str = ",") -> list[str]:
    out, buf, depth = [], [], 0
    for ch in s:
        if ch == "(":
            depth += 1; buf.append(ch)
        elif ch == ")":
            depth -= 1; buf.append(ch)
        elif ch == sep and depth == 0:
            out.append("".join(buf).strip()); buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf).strip())
    return out


@dataclass
class Header:
    indent: str
    indices: list[tuple[str, str, str, str | None]]  # (var, lo, hi, step)
    privates: list[str]
    has_mask: bool


def parse_triplet(item: str) -> tuple[str, str, str, str | None] | None:
    if "=" not in item:
        return None
    name, _, rng = item.partition("=")
    name = name.strip()
    parts = split_top_level(rng, ":")
    if len(parts) < 2:
        return None
    lo = parts[0].strip()
    hi = parts[1].strip()
    step = parts[2].strip() if len(parts) >= 3 else None
    return (name, lo, hi, step)


def find_dc_block(lines: list[str], start: int) -> tuple[Header, int, int] | None:
    """Locate a `do concurrent` block starting at or after `start`.
    Returns (header, header_first_line, header_last_line) or None."""
    i = start
    while i < len(lines):
        c = code_part(lines[i])
        m = DC_RE.search(c)
        if not m:
            i += 1
            continue
        # Join continuations on the header.
        joined_parts = []
        first = i
        j = i
        while j < len(lines):
            piece = code_part(lines[j]).rstrip()
            if piece.endswith("&"):
                joined_parts.append(piece[:-1])
                j += 1
                if j < len(lines):
                    nxt = lines[j].lstrip()
                    if nxt.startswith("&"):
                        lines[j] = lines[j].replace("&", " ", 1)
            else:
                joined_parts.append(piece)
                break
        joined = " ".join(joined_parts)
        last = j

        # Parse header.
        m2 = DC_RE.search(joined)
        rest = joined[m2.end():].lstrip()
        if not rest.startswith("("):
            return None
        depth, end = 0, -1
        for k, ch in enumerate(rest):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    end = k
                    break
        if end < 0:
            return None
        inside = rest[1:end]
        trailing = rest[end + 1:].strip()

        indices: list[tuple[str, str, str, str | None]] = []
        has_mask = False
        for it in split_top_level(inside, ","):
            t = parse_triplet(it)
            if t is not None:
                indices.append(t)
            elif it.strip():
                has_mask = True

        privates: list[str] = []
        t = trailing
        while t:
            tlow = t.lower()
            matched = None
            for c in KNOWN_CLAUSES:
                if tlow.startswith(c) and (len(t) == len(c) or t[len(c)] in "( "):
                    matched = c
                    break
            if not matched:
                break
            after = t[len(matched):].lstrip()
            if after.startswith("("):
                d, jj = 0, -1
                for k, ch in enumerate(after):
                    if ch == "(":
                        d += 1
                    elif ch == ")":
                        d -= 1
                        if d == 0:
                            jj = k
                            break
                if jj < 0:
                    break
                payload = after[1:jj]
                if matched == "local":
                    privates.extend(v.strip() for v in split_top_level(payload, ","))
                t = after[jj + 1:].lstrip()
            else:
                t = after

        if has_mask:
            sys.stderr.write(
                f"  skip: header has mask at line {first + 1}; manual handling needed\n"
            )
            return None
        if not indices:
            return None

        indent = lines[first][: len(lines[first]) - len(lines[first].lstrip())]
        return (Header(indent=indent, indices=indices, privates=privates, has_mask=False),
                first, last)
    return None


def find_matching_end_do(lines: list[str], header_last: int) -> int | None:
    """Scan forward from line after header to find the matching `end do`.
    Tracks `do` depth (any `do` opens, `end do` closes)."""
    depth = 1
    i = header_last + 1
    while i < len(lines):
        c = code_part(lines[i])
        if END_DO_RE.match(c):
            depth -= 1
            if depth == 0:
                return i
        elif DO_RE.match(c):
            depth += 1
        i += 1
    return None


# Worksharing directive emitted per target. CPU drops `target teams distribute`
# (the host-fallback path that is flaky on some compilers, notably GNU) for a
# plain host `parallel do`.
DIRECTIVE_KW = {
    "gpu": "target teams distribute parallel do",
    "cpu": "parallel do",
}


def emit_replacement(h: Header, target: str = "gpu") -> tuple[list[str], list[str]]:
    """Return (header_lines, footer_lines)."""
    n = len(h.indices)
    kw = DIRECTIVE_KW[target]
    parts = [f"!$omp {kw}"]
    if n > 1:
        parts.append(f"collapse({n})")
    if h.privates:
        parts.append(f"private({', '.join(h.privates)})")
    directive = h.indent + " ".join(parts)
    do_lines = []
    for (var, lo, hi, step) in h.indices:
        rng = f"{lo}, {hi}" + (f", {step}" if step else "")
        do_lines.append(f"{h.indent}do {var} = {rng}")
    end_lines = [f"{h.indent}end do" for _ in h.indices]
    end_lines.append(f"{h.indent}!$omp end {kw}")
    return ([directive] + do_lines, end_lines)


def transform_lines(lines: list[str], target: str = "gpu") -> tuple[list[str], int]:
    """Transform all `do concurrent` blocks. Returns (new_lines, n_changed)."""
    out = list(lines)
    changed = 0
    i = 0
    while i < len(out):
        block = find_dc_block(out, i)
        if block is None:
            break
        h, first, last = block
        end_idx = find_matching_end_do(out, last)
        if end_idx is None:
            sys.stderr.write(f"  warn: no matching end do for header at {first + 1}\n")
            i = last + 1
            continue
        new_header, new_footer = emit_replacement(h, target)
        out = (out[:first] + new_header + out[last + 1:end_idx] + new_footer + out[end_idx + 1:])
        changed += 1
        i = first + len(new_header) + (end_idx - last - 1) + len(new_footer)
    return (out, changed)


PROC_NAME_RE = re.compile(r"\b(?:subroutine|function)\b\s+([a-zA-Z]\w*)",
                          re.IGNORECASE)


def _parse_procedures(lines: list[str]) -> list[dict]:
    """List every subroutine/function as {open, end, name}.

    Stack-matches procedure-opening statements to their `end subroutine` /
    `end function` so nested (contained) procedures pair correctly.
    """
    procs: list[dict] = []
    stack: list[tuple[int, str]] = []
    for idx, raw in enumerate(lines):
        c = code_part(raw)
        if PROC_END_RE.match(c):
            if stack:
                open_idx, name = stack.pop()
                procs.append({"open": open_idx, "end": idx, "name": name})
        elif PROC_OPEN_RE.match(c):
            m = PROC_NAME_RE.search(c)
            stack.append((idx, m.group(1).lower() if m else ""))
    return procs


def seed_target_strips(lines: list[str], procs: list[dict],
                       open_re: re.Pattern = OMP_TARGET_RE) -> set[int]:
    """Open-line indices of procedures that DIRECTLY contain an opening
    worksharing directive (`!$omp target` on gpu, `!$omp parallel do` on cpu) —
    those launch a kernel / spawn a team (a side effect) and so can't stay
    `pure` — OR an ATOMIC, which OpenMP forbids in a pure procedure outright.
    The nearest enclosing procedure of each such directive is the one that
    loses purity, and the caller cascade takes it from there."""
    seed: set[int] = set()
    for idx, raw in enumerate(lines):
        if not (open_re.match(raw) or OMP_ATOMIC_RE.match(raw)):
            continue
        enclosing = None
        for p in procs:
            if p["open"] < idx < p["end"]:
                if enclosing is None or (p["end"] - p["open"]) < (
                    enclosing["end"] - enclosing["open"]):
                    enclosing = p
        if enclosing is not None:
            seed.add(enclosing["open"])
    return seed


IDENT_RE = re.compile(r"[A-Za-z_]\w*")


def _convert_task(path: Path, target: str = "gpu") -> tuple:
    """Phase-1 unit of work: convert one file and seed its local impure names.

    Self-contained (reads its own file, shares nothing), so it is safe to run in
    a worker process. Returns the `unit` tuple plus the names this file
    contributes to the global impure set.
    """
    open_re = OMP_PARALLEL_DO_RE if target == "cpu" else OMP_TARGET_RE
    text = path.read_text()
    lines, n_loops = transform_lines(text.splitlines(), target)
    procs = _parse_procedures(lines)
    strip_idx = seed_target_strips(lines, procs, open_re)
    # Seed with every procedure that is non-pure after conversion: the
    # target-stripped ones AND any procedure that already lacks `pure`
    # (legitimately impure, or stripped by a prior run — re-running must still
    # cascade to ITS pure callers). In valid Fortran a pure procedure never
    # calls an already-impure one, so the already-impure names never trigger a
    # spurious strip.
    local_impure = {p["name"] for p in procs
                    if p["open"] in strip_idx
                    or not PURE_RE.search(lines[p["open"]])}
    return (path, lines, n_loops, procs, strip_idx,
            text.endswith("\n"), local_impure)


def _write_task(unit: tuple, write: bool = False) -> tuple:
    """Phase-3 unit of work: apply this file's pure-strips and write it back."""
    path, lines, n_loops, _procs, strip_idx, trailing_nl = unit
    n_pure = 0
    for idx in strip_idx:
        if PURE_RE.search(lines[idx]):
            lines[idx] = PURE_RE.sub("", lines[idx], count=1)
            n_pure += 1
    if (n_loops or n_pure) and write:
        path.write_text("\n".join(lines) + ("\n" if trailing_nl else ""))
    return path, n_loops, n_pure


def walk(roots: list[Path]) -> list[Path]:
    files: list[Path] = []
    for r in roots:
        if r.is_file():
            files.append(r); continue
        for p in sorted(r.rglob("*.F90")):
            files.append(p)
        for p in sorted(r.rglob("*.f90")):
            files.append(p)
    return files


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="apply changes in place (default: dry-run)")
    ap.add_argument(
        "--target", choices=["gpu", "cpu"], default="gpu",
        help="gpu (default): rewrite do concurrent as !$omp target teams "
             "distribute parallel do. cpu: !$omp parallel do (host multicore).")
    ap.add_argument(
        "-j", "--jobs", "-t", dest="jobs", type=int, default=1, metavar="N",
        help="run the per-file phases across N worker processes (0 = one per "
             "available core). The whole-program pure cascade stays serial; "
             "output is order-independent either way.")
    ap.add_argument("paths", nargs="*", default=["src"])
    args = ap.parse_args()

    roots = [Path(p) for p in args.paths]
    for r in roots:
        if not r.exists():
            print(f"error: {r} does not exist", file=sys.stderr)
            return 1
    jobs = resolve_jobs(args.jobs)

    # ---- Phase 1: per-file DC→target conversion + seed the impure set ----
    # `procs` are parsed from the POST-conversion lines so body spans line up.
    # `impure` is the GLOBAL set of procedure names that lose `pure` — seeded
    # with every procedure that directly contains an `!$omp target`.
    units = []  # [path, lines, n_loops, procs, strip_idx, trailing_nl]
    impure: set[str] = set()
    for (path, lines, n_loops, procs, strip_idx, trailing_nl,
         local_impure) in pmap(partial(_convert_task, target=args.target),
                               walk(roots), jobs, "[1/3] convert"):
        impure |= local_impure
        units.append([path, lines, n_loops, procs, strip_idx, trailing_nl])

    # ---- Phase 2: WHOLE-PROGRAM cascade to a fixpoint ----
    # A `pure` procedure that references an impure procedure (in ANY file —
    # a pure procedure may only call pure procedures) must lose `pure` too.
    # Iterate across all files until the impure set stops growing, so the
    # impurity propagates up cross-module call chains (e.g. ml_seed_tracer →
    # ml_init_tracer lives in one file, but its caller in rki_multilayer_state
    # is in another).
    #
    # Run as a worklist over a reverse index rather than re-sweeping every
    # file. The old form rebuilt each procedure body and ran one regex per
    # impure name over it, on every sweep — O(sweeps x procs x |impure|), and
    # |impure| grows into the thousands, which is what made this the slow
    # stage. A `\b`-anchored case-insensitive search for a Fortran identifier
    # hits iff that identifier occurs as a token in the body, so a lowercase
    # token set is an exact substitute for the regex; the fixpoint is
    # unchanged, it is just reached by propagating along the call graph
    # instead of rediscovering it each pass.
    candidates = []            # [unit_idx, open_idx, name_lc] still strippable
    referenced_by: dict[str, list[int]] = {}
    for u_i, (_p, lines, _nl, procs, strip_idx, _t) in enumerate(
            track(units, "[2/3] index")):
        for p in procs:
            if p["open"] in strip_idx:
                continue
            if not PURE_RE.search(lines[p["open"]]):
                continue  # already non-pure / nothing to strip
            body = "\n".join(code_part(lines[k])
                             for k in range(p["open"] + 1, p["end"]))
            c_i = len(candidates)
            candidates.append([u_i, p["open"], p["name"].lower()])
            for tok in {m.group(0).lower() for m in IDENT_RE.finditer(body)}:
                referenced_by.setdefault(tok, []).append(c_i)

    cascade = Progress("[2/3] pure cascade", total=None)
    cascade.draw(force=True)
    stripped = [False] * len(candidates)
    impure_lc = {c.lower() for c in impure if c}
    work = list(impure_lc)
    while work:
        callee = work.pop()
        for c_i in referenced_by.get(callee, ()):
            if stripped[c_i]:
                continue
            u_i, open_idx, name_lc = candidates[c_i]
            stripped[c_i] = True
            units[u_i][4].add(open_idx)
            cascade.advance(suffix=f"{len(impure_lc)} impure")
            if name_lc and name_lc not in impure_lc:
                impure_lc.add(name_lc)
                work.append(name_lc)
    cascade.done(f"{len(impure_lc)} impure procedures")

    # ---- Phase 3: apply pure-strips + write ----
    total = 0
    total_pure = 0
    files_changed = 0
    for path, n_loops, n_pure in pmap(partial(_write_task, write=args.write),
                                      units, jobs, "[3/3] write"):
        if n_loops == 0 and n_pure == 0:
            continue
        total += n_loops
        total_pure += n_pure
        files_changed += 1
        verb = "rewrote" if args.write else "would rewrite"
        extra = f"  (-{n_pure} pure)" if n_pure else ""
        print(f"  {verb} {n_loops:3d}  {path}{extra}")
    mode = "WRITE" if args.write else "DRY RUN"
    print(f"\n[{mode}] {total} loops, {total_pure} pure attrs removed "
          f"across {files_changed} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
