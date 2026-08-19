#!/usr/bin/env python3
"""Convert OpenACC directives to OpenMP target equivalents.

Directive-block aware: joins `&`-continued directives into one logical block,
applies head + clause transformations uniformly across the block, and re-emits.

Head substitutions:
  !$acc enter data ...          →  !$omp target enter data ...
  !$acc exit data ...           →  !$omp target exit data ...
  !$acc data ...                →  !$omp target data ...
  !$acc end data                →  !$omp end target data
  !$acc update ...              →  !$omp target update ...
  !$acc parallel loop ...       →  !$omp target teams distribute parallel do ...
  !$acc host_data use_device(   →  !$omp target data use_device_addr(
  !$acc end host_data           →  !$omp end target data
  !$acc declare create(X)       →  !$omp declare target(X)
  !$acc set device_num(X)       →  call omp_set_default_device(X)   (warn: needs omp_lib)
  !$acc routine seq             →  !$omp declare target
  !$acc atomic [update]         →  LEFT ALONE. OpenMP forbids an atomic in a
                                   PURE procedure and OpenACC allows it, and a
                                   procedure called from `do concurrent` has to
                                   be pure. Inert under OpenMP, live under
                                   OpenACC, correct either way.
  !$acc&  (continuation)        →  !$omp&

Clause translations (applied across the directive block — head + continuations):
  copyin(  → map(to:
  copyout( → map(from:
  copy(    → map(tofrom:
  create(  → map(alloc:
  delete(  → map(delete:
  self(    → from(    (only on update directives)
  device(  → to(      (only on update directives)

Stripped (OpenACC-only, no OpenMP equivalent — emit as-is would fail):
  present(...)        whole clause + args
  if_present          keyword

Left untouched:
  do concurrent loops
  !$omp directives that already exist

Targets (--target):
  gpu (default)  Head substitutions above — `!$omp target ...` offload, with
                 copyin/out/create translated to map().
  cpu            Host multicore. The compute loop becomes `!$omp parallel do`
                 (no `target teams distribute` — that path is flaky on some host
                 fallbacks, notably GNU). Every device-data / offload directive
                 (enter/exit data, data, update, host_data, declare, set
                 device_num, routine seq) is a host no-op and is neutralised to
                 an inert comment. Data/locality clauses on the compute loop
                 (present, copyin, ...) are stripped (invalid on parallel do).
                 Build the result with -DRAKALI_PARALLEL_BACKEND=openmp
                 -DRAKALI_ENABLE_GPU=OFF.

Usage:
    python tools/acc_to_omp.py             # dry-run on src/ (gpu)
    python tools/acc_to_omp.py --write     # apply (gpu)
    python tools/acc_to_omp.py --target cpu --write src tests benchmarks
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from functools import partial
from pathlib import Path

from _parallel import pmap, resolve_jobs

ACC_DIR = re.compile(r"^(\s*)!\$acc\b", re.IGNORECASE)
ACC_CONT = re.compile(r"^(\s*)!\$acc&", re.IGNORECASE)
OMP_DIR = re.compile(r"^(\s*)!\$omp\b", re.IGNORECASE)
OMP_CONT = re.compile(r"^(\s*)!\$omp&", re.IGNORECASE)
OMP_TARGET_HEAD = re.compile(
    r"!\$omp\s+(target|declare\s+target|end\s+target)\b", re.IGNORECASE,
)

# Head replacements. Order matters: longer/more-specific first.
HEAD_PATTERNS: list[tuple[re.Pattern, str, str]] = [
    #
    # ATOMICS ARE RECOGNISED AND DELIBERATELY LEFT ALONE.
    #
    # `!$omp atomic` and `!$acc atomic` mean the same thing, so translating
    # them looks obviously right. It is not:
    #
    #     Error: OpenMP directives other than SIMD or DECLARE TARGET at (1)
    #            may not appear in PURE procedures
    #
    # OpenMP forbids an atomic inside a PURE procedure and OpenACC permits it.
    # A procedure called from `do concurrent` MUST be pure, so any code whose
    # scatter lives in a pure routine -- which is the whole point of the
    # construct -- cannot have its atomics translated.
    #
    # Left untranslated the directive is an inert comment under an OpenMP
    # compiler, which is correct wherever `do concurrent` runs serially, and
    # still live under an OpenACC one, which is correct wherever it does not.
    # Neutralising it to a comment instead would be wrong under
    # `-stdpar=multicore`, where the loop really is threaded.
    #
    # Listed first so `!$acc atomic update` is never mistaken for an `update`
    # directive on the way past.
    #
    (re.compile(r"!\$acc\s+atomic\b", re.IGNORECASE),
     None, "atomic"),
    (re.compile(r"!\$acc\s+routine\s+seq\b", re.IGNORECASE),
     "!$omp declare target", "routine_seq"),
    (re.compile(r"!\$acc\s+declare\s+create\s*\(", re.IGNORECASE),
     "!$omp declare target(", "declare"),
    (re.compile(r"!\$acc\s+set\s+device_num\s*\(", re.IGNORECASE),
     "call omp_set_default_device(", "set_device"),
    (re.compile(r"!\$acc\s+enter\s+data\b", re.IGNORECASE),
     "!$omp target enter data", "enter_data"),
    (re.compile(r"!\$acc\s+exit\s+data\b", re.IGNORECASE),
     "!$omp target exit data", "exit_data"),
    (re.compile(r"!\$acc\s+end\s+data\b", re.IGNORECASE),
     "!$omp end target data", "end_data"),
    (re.compile(r"!\$acc\s+data\b", re.IGNORECASE),
     "!$omp target data", "data"),
    (re.compile(r"!\$acc\s+update\b", re.IGNORECASE),
     "!$omp target update", "update"),
    (re.compile(r"!\$acc\s+end\s+parallel\s+loop\b", re.IGNORECASE),
     "!$omp end target teams distribute parallel do", "end_parallel_loop"),
    (re.compile(r"!\$acc\s+parallel\s+loop\b", re.IGNORECASE),
     "!$omp target teams distribute parallel do", "parallel_loop"),
    (re.compile(r"!\$acc\s+host_data\s+use_device\s*\(", re.IGNORECASE),
     "!$omp target data use_device_addr(", "host_data"),
    (re.compile(r"!\$acc\s+end\s+host_data\b", re.IGNORECASE),
     "!$omp end target data", "end_host_data"),
]

# Clause translations applied to the body of the directive block.
# `self(` and `device(` are only translated for update directives.
CLAUSE_MAP_DATA = {
    "copyin":             "map(to:",
    "copyout":            "map(from:",
    "copy":               "map(tofrom:",
    "create":             "map(alloc:",
    "delete":             "map(delete:",
    # `present_or_*` forms behave like `*` in OpenACC 3.0+ — "reuse if
    # already mapped, otherwise copy/alloc". OpenMP map() has the same
    # reference-count semantic, so the translation is the same target.
    "present_or_copyin":  "map(to:",
    "present_or_copyout": "map(from:",
    "present_or_copy":    "map(tofrom:",
    "present_or_create":  "map(alloc:",
    # Short-form aliases supported by NVHPC.
    "pcopyin":            "map(to:",
    "pcopyout":           "map(from:",
    "pcopy":              "map(tofrom:",
    "pcreate":            "map(alloc:",
}
CLAUSE_MAP_UPDATE = {
    "self":   "from(",
    "device": "to(",
}
CLAUSE_RE_DATA = re.compile(
    r"\b(" + "|".join(sorted(CLAUSE_MAP_DATA, key=len, reverse=True)) + r")\s*\(",
    re.IGNORECASE,
)
CLAUSE_RE_UPDATE = re.compile(
    r"\b(" + "|".join(sorted(CLAUSE_MAP_UPDATE, key=len, reverse=True)) + r")\s*\(",
    re.IGNORECASE,
)

# OpenACC-only clauses to strip.
PRESENT_RE = re.compile(r"\bpresent\s*\(", re.IGNORECASE)
IF_PRESENT_RE = re.compile(r"\s+if_present\b", re.IGNORECASE)
# `default(present)` is an OpenACC-only clause: it asserts every referenced
# array is already device-resident.  OpenMP has no equivalent (it is NOT
# `defaultmap(present)`, which means something else), and both ifx and
# amdflang reject it — the Frontier/Intel OpenMP build failure.
# PRESENT_RE above cannot catch it: that matches `present` followed by `(`,
# whereas here `present` is the ARGUMENT.  Needs its own strip.
DEFAULT_PRESENT_RE = re.compile(r"\s+default\s*\(\s*present\s*\)", re.IGNORECASE)

# OpenACC-only LOOP-mapping clauses. Invalid on an OpenMP loop construct
# (`target teams distribute parallel do` on GPU, `parallel do` on host) — the
# parallelism is expressed by the construct itself, not by clauses, so these
# are dropped. Parenthesised forms (gang(n)/vector(n)/vector_length(n)/
# num_gangs(n)/num_workers(n)/tile(...)/device_type(...)) are removed with
# balanced parens; the bare keywords (gang/vector/worker/independent/auto/seq)
# only at paren-depth 0, so a variable named e.g. `gang` inside a
# reduction(...)/private(...) clause is never touched. `vector_length` is
# ordered before `vector` so the longer name wins.
ACC_LOOP_PAREN_RE = re.compile(
    r"\b(num_gangs|num_workers|vector_length|gang|vector|worker|tile|device_type)\s*\(",
    re.IGNORECASE,
)
ACC_LOOP_BARE_RE = re.compile(
    r"(?i)(\s+)(gang|vector|worker|independent|auto|seq)\b(?!\s*\()"
)


def strip_acc_loop_clauses(joined: str) -> tuple[str, int]:
    """Remove OpenACC-only loop-mapping clauses from a converted loop directive
    block. Returns (joined, count). Parenthesised forms go via balanced-paren
    removal; bare keywords are removed only where the paren depth is 0 (so they
    are never mistaken for a variable inside another clause's parens)."""
    joined, n_paren = strip_balanced_clause(joined, ACC_LOOP_PAREN_RE)
    # `default(present)`: OpenACC-only, and NOT reachable by the `present(`
    # strip elsewhere (there `present` is the argument, not the clause name).
    # Left in place it produces OpenMP that ifx and amdflang both reject.
    joined, n_defpres = DEFAULT_PRESENT_RE.subn("", joined)
    n_bare = 0
    while True:
        hit = None
        for m in ACC_LOOP_BARE_RE.finditer(joined):
            kw = m.start(2)
            if joined.count("(", 0, kw) == joined.count(")", 0, kw):  # depth 0
                hit = m
                break
        if hit is None:
            break
        joined = joined[:hit.start()] + joined[hit.end():]
        n_bare += 1
    return joined, n_paren + n_bare + n_defpres

# --- CPU (multicore host) target ---------------------------------------------
# On the host there is no separate device memory, so every data-movement /
# offload-management directive is a no-op and is neutralised to an inert comment.
# Only the compute loop survives, rewritten to a host `!$omp parallel do`.
CPU_DROP_KINDS = {
    "enter_data", "exit_data", "end_data", "data", "update",
    "host_data", "end_host_data", "declare", "set_device", "routine_seq",
}
# Short human label per kind for the inert drop-comment.
CPU_DROP_LABEL = {
    "enter_data": "enter data", "exit_data": "exit data", "end_data": "end data",
    "data": "data", "update": "update", "host_data": "host_data",
    "end_host_data": "end host_data", "declare": "declare",
    "set_device": "set device_num", "routine_seq": "routine seq",
}
# Data/locality clauses valid on `!$acc parallel loop` but NOT on a host
# `!$omp parallel do` (no map/copy semantics) — stripped whole, balanced parens.
CPU_STRIP_RE = re.compile(
    r"\b(" + "|".join(sorted(set(CLAUSE_MAP_DATA) | {"present"},
                             key=len, reverse=True)) + r")\s*\(",
    re.IGNORECASE,
)


def collect_block(lines: list[str], i: int, kind: str) -> tuple[list[str], int]:
    """Collect a directive block of `kind` ('acc' or 'omp') starting at lines[i],
    following &-continuations. Return (block, last_index)."""
    block = [lines[i]]
    j = i
    head_re = ACC_DIR if kind == "acc" else OMP_DIR
    cont_re = ACC_CONT if kind == "acc" else OMP_CONT
    while block[-1].rstrip().rstrip("\n").endswith("&") and j + 1 < len(lines):
        nxt = lines[j + 1]
        if head_re.match(nxt) or cont_re.match(nxt):
            block.append(nxt)
            j += 1
        else:
            break
    return block, j


def strip_balanced_clause(joined: str, clause_re: re.Pattern) -> tuple[str, int]:
    """Remove every `clause(...)` matched by `clause_re`, honouring nested
    parens. Drops a dangling `&` that was only holding the removed clause when
    nothing else follows on that continuation fragment. Returns (joined, count)."""
    count = 0
    while True:
        m = clause_re.search(joined)
        if not m:
            break
        depth, end = 0, -1
        for k in range(m.end() - 1, len(joined)):
            ch = joined[k]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    end = k
                    break
        if end < 0:
            break  # malformed — leave alone
        before = joined[:m.start()].rstrip(" \t")
        after = joined[end + 1:]
        if before.endswith("&") and after.lstrip(" \t").startswith("\n"):
            before = before[:-1].rstrip(" \t")
        joined = before + (" " if not after.startswith(("\n", " ")) else "") + after
        count += 1
    return joined, count


def normalize_continuations(joined: str) -> str:
    """Tidy a converted directive block's continuation lines.

    Drops sentinel-only `!$omp` / `!$omp&` lines (left behind when a strip,
    e.g. `present(...)`, empties a continuation line down to its sentinel),
    then removes a now-dangling trailing `&` on the final directive line.
    A directive block's last line must never end with `&`."""
    parts = joined.splitlines(keepends=True)
    parts = [p for p in parts
             if not re.fullmatch(r"[ \t]*!\$omp&?[ \t]*\n?", p, re.IGNORECASE)]
    if parts:
        last = parts[-1]
        nl = "\n" if last.endswith("\n") else ""
        body = (last[:-1] if nl else last).rstrip()
        if body.endswith("&"):
            parts[-1] = body[:-1].rstrip() + nl
    return "".join(parts)


def cleanup_omp_block(block: list[str], stats: Counter) -> list[str]:
    """Idempotent cleanup pass for already-converted !$omp target blocks:
    translates any leftover OpenACC-style clauses and strips OpenACC-only
    clauses. Safe to run on fully-converted code (no-op)."""
    joined = "".join(block)

    if not OMP_TARGET_HEAD.search(joined):
        return block  # only touch target / declare target / end target blocks

    # Strip present(...)
    while True:
        m = PRESENT_RE.search(joined)
        if not m:
            break
        depth, end = 0, -1
        for k in range(m.end() - 1, len(joined)):
            ch = joined[k]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    end = k
                    break
        if end < 0:
            break
        before = joined[:m.start()].rstrip(" \t")
        after = joined[end + 1:]
        if before.endswith("&") and after.lstrip(" \t").startswith("\n"):
            before = before[:-1].rstrip(" \t")
        joined = before + (" " if not after.startswith(("\n", " ")) else "") + after
        stats["present_stripped"] += 1

    # Strip if_present
    joined, n_ifp = IF_PRESENT_RE.subn("", joined)
    stats["if_present_stripped"] += n_ifp

    # Strip default(present) — OpenACC-only, invalid OpenMP (see the regex).
    joined, n_dp = DEFAULT_PRESENT_RE.subn("", joined)
    stats["default_present_stripped"] += n_dp

    # Translate data-style clauses (copy/copyin/copyout/create/delete) anywhere
    # in the block — including continuation lines.
    joined, n_data = CLAUSE_RE_DATA.subn(
        lambda m: CLAUSE_MAP_DATA[m.group(1).lower()] + " ", joined,
    )
    stats["clause_data"] += n_data

    # Translate update-style clauses only on omp target update blocks.
    if re.search(r"!\$omp\s+target\s+update\b", joined, re.IGNORECASE):
        joined, n_upd = CLAUSE_RE_UPDATE.subn(
            lambda m: CLAUSE_MAP_UPDATE[m.group(1).lower()], joined,
        )
        stats["clause_update"] += n_upd

    # Drop degenerate `!$omp& &` empty continuations.
    joined = re.sub(
        r"^\s*!\$omp&\s*&?\s*\n", "", joined,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    # Drop sentinel-only continuations + dangling trailing `&` (e.g. left when a
    # present() strip empties a continuation line).
    joined = normalize_continuations(joined)

    new_lines = joined.splitlines(keepends=True)
    if block and not block[-1].endswith("\n") and new_lines and new_lines[-1].endswith("\n"):
        new_lines[-1] = new_lines[-1].rstrip("\n")
    return new_lines


def transform_cpu_block(block: list[str], joined: str, head_kind: str | None,
                        stats: Counter) -> list[str]:
    """CPU (multicore host) variant of a directive block.

    Compute loops become `!$omp parallel do`; all device-data / offload
    directives are host no-ops, neutralised to an inert comment (a single line,
    collapsing any continuations). Unrecognised acc directives are left intact."""
    if head_kind in CPU_DROP_KINDS:
        indent = re.match(r"^(\s*)", block[0]).group(1)
        stats["cpu_dropped"] += 1
        return [f"{indent}! [acc->omp-cpu] dropped device directive "
                f"(host no-op): {CPU_DROP_LABEL[head_kind]}\n"]

    if head_kind == "atomic":
        # Untranslated on purpose -- see HEAD_PATTERNS. Not dropped either:
        # dropping it would race under `-stdpar=multicore`.
        stats["atomic_left_as_is"] += 1
        return block

    if head_kind == "end_parallel_loop":
        joined = re.sub(r"!\$acc\s+end\s+parallel\s+loop\b", "!$omp end parallel do",
                        joined, count=1, flags=re.IGNORECASE)
        stats["end_parallel_loop"] += 1
        return joined.splitlines(keepends=True)

    if head_kind == "parallel_loop":
        # Strip data/locality clauses that a host parallel do can't carry
        # (present, present_or_*, copyin/out, create, ...) + the if_present kw.
        joined, n_strip = strip_balanced_clause(joined, CPU_STRIP_RE)
        stats["cpu_clause_stripped"] += n_strip
        joined, n_ifp = IF_PRESENT_RE.subn("", joined)
        stats["if_present_stripped"] += n_ifp
        # Drop OpenACC-only loop-mapping clauses (gang/vector/worker/...) —
        # invalid on a host `!$omp parallel do` too.
        joined, n_loop = strip_acc_loop_clauses(joined)
        stats["acc_loop_clause_stripped"] += n_loop
        joined = re.sub(r"!\$acc\s+parallel\s+loop\b", "!$omp parallel do",
                        joined, count=1, flags=re.IGNORECASE)
        stats["parallel_loop"] += 1
        # Continuation sentinels (both !$acc& and the !$acc <clause> space form).
        joined, n_cont = re.subn(r"(?im)^([ \t]*)!\$acc\b", r"\1!$omp", joined)
        stats["continuation"] += n_cont
        joined = normalize_continuations(joined)
        new_lines = joined.splitlines(keepends=True)
        if block and not block[-1].endswith("\n") and new_lines and new_lines[-1].endswith("\n"):
            new_lines[-1] = new_lines[-1].rstrip("\n")
        return new_lines

    return block  # unrecognised acc directive — leave untouched


def transform_block(block: list[str], stats: Counter, warnings: list[str],
                    path: Path, lineno: int, target: str = "gpu") -> list[str]:
    """Apply head + clause + strip transformations to a directive block."""
    joined = "".join(block)

    # 0. Drop CUDA-graph-capture wrappers: `!$acc kernels [async(...)]`,
    #    `!$acc end kernels`, `!$acc wait(...)`.  These wrap the do-concurrent
    #    loops so stdpar can capture them into a CUDA graph; under OpenMP
    #    offload they're meaningless, and left live an async `!$acc kernels`
    #    region around an `!$omp target` launch corrupts the device stream /
    #    data attach (null component descriptors, cuStreamSynchronize errors).
    #    Neutralise to an inert plain comment (traceable, no `!$acc` sentinel).
    #    Applies to both targets (CPU has no CUDA graph either).
    m_wrap = re.match(r"^(\s*)!\$acc\s+(end\s+kernels|kernels|wait)\b(.*)$",
                      block[0], re.IGNORECASE)
    if m_wrap:
        indent = m_wrap.group(1)
        detail = (m_wrap.group(2) + " " + m_wrap.group(3)).strip()
        detail = re.sub(r"\s+", " ", detail)
        stats["graph_wrapper_dropped"] += 1
        return [f"{indent}! [acc->omp] dropped CUDA-graph wrapper: {detail}\n"]

    # 0b. Leave any directive carrying an OpenACC `async(...)` clause UNTRANSLATED.
    #    OpenACC async queues (e.g. the SI CG launch-batching `!$acc parallel loop
    #    ... async(cg_q)`) have no OpenMP-target equivalent; translating the head
    #    while keeping `async(cg_q)` emits an invalid clause on an `!$omp target
    #    teams distribute parallel do` and breaks the AMD/Intel offload build.
    #    Left as-is, the `!$acc` directive is inert under an OpenMP compiler, so
    #    the loop runs correctly (synchronously). `!$acc kernels/wait async(...)`
    #    is already neutralised by step 0 above; this catches the parallel-loop
    #    form. (The async path will be unified so this duplication goes away.)
    if re.search(r"\basync\s*\(", joined, re.IGNORECASE):
        stats["async_skipped"] += 1
        return block

    # CPU target diverges entirely from here: device-data directives vanish and
    # the compute loop becomes a host worksharing region. Detect the head kind
    # (search only) and hand off.
    if target == "cpu":
        head_kind = next((name for pat, _, name in HEAD_PATTERNS
                          if pat.search(joined)), None)
        return transform_cpu_block(block, joined, head_kind, stats)

    # 1. Strip `present(...)` clause (with balanced parens).
    joined, n_pres = strip_balanced_clause(joined, PRESENT_RE)
    stats["present_stripped"] += n_pres

    # 2. Strip `if_present` keyword.
    joined, n_ifp = IF_PRESENT_RE.subn("", joined)
    stats["if_present_stripped"] += n_ifp

    # 3. Determine directive type from the first line (post-strip but pre-head-sub),
    #    so we know whether to apply update-specific clause translations.
    first_line = joined.splitlines()[0] if joined else ""
    is_update = bool(re.search(r"!\$acc\s+update\b", first_line, re.IGNORECASE))

    # 4. Apply head substitution (first match wins — patterns are ordered).
    head_kind = None
    for pat, repl, name in HEAD_PATTERNS:
        if pat.search(joined):
            if repl is None:
                # Recognised, deliberately untranslated -- see HEAD_PATTERNS.
                stats[name + "_left_as_is"] += 1
                return block
            joined = pat.sub(repl, joined, count=1)
            head_kind = name
            stats[name] += 1
            if name == "set_device":
                warnings.append(
                    f"  {path}:{lineno}  !$acc set device_num → omp_set_default_device "
                    f"(ensure surrounding scope has `use omp_lib`)"
                )
            break

    # 5. Translate clauses (data-style on data/enter/exit/parallel/declare;
    #    update-style on update directives).
    if head_kind in ("enter_data", "exit_data", "data", "parallel_loop"):
        joined, n = CLAUSE_RE_DATA.subn(
            lambda m: CLAUSE_MAP_DATA[m.group(1).lower()] + " ", joined,
        )
        stats["clause_data"] += n
    elif head_kind == "update":
        joined, n = CLAUSE_RE_UPDATE.subn(
            lambda m: CLAUSE_MAP_UPDATE[m.group(1).lower()] + "", joined,
        )
        stats["clause_update"] += n

    # 5b. Drop OpenACC-only loop-mapping clauses (gang/vector/worker/num_gangs/
    #     vector_length/...). They are invalid on the OpenMP combined loop
    #     construct; the parallelism is carried by `teams distribute parallel
    #     do`. Without this, a `!$acc parallel loop gang vector` translates to
    #     `!$omp target teams distribute parallel do gang vector` (ifx #5082).
    if head_kind == "parallel_loop":
        joined, n_loop = strip_acc_loop_clauses(joined)
        stats["acc_loop_clause_stripped"] += n_loop

    # 6. Convert continuation sentinels to !$omp. The head sentinel was already
    #    replaced in step 4, so any remaining line-leading !$acc is a
    #    &-continuation of this same block — covering BOTH the `!$acc&` form and
    #    the `!$acc   <clause>` space form this codebase actually uses (the latter
    #    was previously left as a stray !$acc under an !$omp head — invalid
    #    OpenMP). Guarded on head_kind so an unrecognised acc directive is left
    #    fully intact rather than half-converted.
    if head_kind is not None:
        joined, n_cont = re.subn(r"(?im)^([ \t]*)!\$acc\b", r"\1!$omp", joined)
        stats["continuation"] += n_cont

    # 7. Tidy continuations emptied by the strips above (sentinel-only !$omp
    #    lines + any resulting dangling trailing `&`).
    joined = normalize_continuations(joined)

    # Re-split into lines, restoring trailing newlines that the original had.
    new_lines = joined.splitlines(keepends=True)
    # If the last original line had no trailing newline, mirror that.
    if block and not block[-1].endswith("\n") and new_lines and new_lines[-1].endswith("\n"):
        new_lines[-1] = new_lines[-1].rstrip("\n")
    return new_lines


DECLARE_TARGET_RE = re.compile(r"^\s*!\$omp\s+declare\s+target\s*$", re.IGNORECASE)
# Skippable between duplicates: blank lines and PLAIN comments (`!` not `!$`),
# including `!!` FORD docstrings. Any sentinel directive (`!$...`) other than the
# bare declare target itself ends the run.
_SKIP_RE = re.compile(r"^\s*(!(?!\$).*)?$")


BARE_TARGET_DATA_RE = re.compile(r"^(\s*)!\$omp\s+target\s+data\s*$", re.IGNORECASE)
TARGET_DATA_HEAD_RE = re.compile(r"^\s*!\$omp\s+target\s+data\b", re.IGNORECASE)
END_TARGET_DATA_RE = re.compile(r"^(\s*)!\$omp\s+end\s+target\s+data\s*$", re.IGNORECASE)


def drop_bare_target_data(lines: list[str]) -> tuple[list[str], int]:
    """Drop `!$omp target data` regions that ended up with NO clauses, plus
    their matching `!$omp end target data`.

    An `!$acc data present(...)` is a pure ASSERTION that the arrays are
    already device-resident — it moves nothing.  Stripping `present(...)`
    (invalid in OpenMP) therefore leaves a bare `!$omp target data`, which is
    not merely useless but ILL-FORMED: OpenMP requires at least one of
    `map` / `use_device_ptr` / `use_device_addr` on `target data`.  amdflang
    rejects it outright.  The correct OpenMP rendering of the original intent
    is nothing at all.

    Matching end is found by depth-counting `target data` opens against
    `end target data` closes, so a bare region nested inside a clause-carrying
    one (e.g. from `host_data use_device_addr`) drops only its own end.

    Mirrors what the CPU target already does for these kinds, and the
    `[acc->omp]` drop-comment idiom used elsewhere.
    """
    out, n = [], 0
    depth_stack: list[bool] = []          # True => this open was dropped
    for line in lines:
        if TARGET_DATA_HEAD_RE.match(line):
            bare = BARE_TARGET_DATA_RE.match(line)
            depth_stack.append(bool(bare))
            if bare:
                out.append("%s! [acc->omp] dropped clause-less target data: "
                           "acc data present(...) asserts residency, moves nothing "
                           "(bare target data is invalid OpenMP)\n" % bare.group(1))
                n += 1
                continue
            out.append(line)
            continue
        m_end = END_TARGET_DATA_RE.match(line)
        if m_end and depth_stack:
            was_dropped = depth_stack.pop()
            if was_dropped:
                out.append("%s! [acc->omp] dropped clause-less end target data\n"
                           % m_end.group(1))
                n += 1
                continue
        out.append(line)
    return out, n


def dedup_declare_target(lines: list[str]) -> tuple[list[str], int]:
    """Collapse a run of identical bare `!$omp declare target` directives to one.

    A device helper's source carries BOTH `!$acc routine seq` (the OpenACC/stdpar
    build) AND a hand-written `!$omp declare target` (a direct OpenMP build of
    main, where the `!$acc` line is inert). Converting the former yields a
    SECOND, adjacent `!$omp declare target`, which ifx rejects with #8682
    ("Multiple DECLARE TARGET"). Keep the first; drop the rest of the run.

    Blank / plain-comment lines between duplicates are preserved and do not end
    the run; any real Fortran statement or other `!$` directive resets it, so
    declare targets in different procedures are never merged. Returns
    (lines, n_removed)."""
    out: list[str] = []
    removed = 0
    seen = False
    for ln in lines:
        if DECLARE_TARGET_RE.match(ln):
            if seen:
                removed += 1
                continue
            seen = True
            out.append(ln)
        elif _SKIP_RE.match(ln.rstrip("\n")):
            out.append(ln)
        else:
            seen = False
            out.append(ln)
    return out, removed


def process_file(path: Path, write: bool, stats: Counter, warnings: list[str],
                 target: str = "gpu") -> int:
    text = path.read_text()
    lines = text.splitlines(keepends=True)
    out: list[str] = []
    changed = 0
    i = 0
    while i < len(lines):
        line = lines[i]
        if ACC_DIR.match(line):
            block, last = collect_block(lines, i, "acc")
            new_block = transform_block(block, stats, warnings, path, i + 1, target)
            if "".join(new_block) != "".join(block):
                changed += len(block)
            out.extend(new_block)
            i = last + 1
        elif OMP_DIR.match(line) and OMP_TARGET_HEAD.search(line):
            # Idempotent cleanup of already-converted !$omp target blocks.
            block, last = collect_block(lines, i, "omp")
            new_block = cleanup_omp_block(block, stats)
            if "".join(new_block) != "".join(block):
                changed += len(block)
            out.extend(new_block)
            i = last + 1
        else:
            out.append(line)
            i += 1
    # Collapse any duplicate `!$omp declare target` produced when the source
    # carries both `!$acc routine seq` and a hand-written `!$omp declare target`
    # (ifx #8682).
    out, n_bare_td = drop_bare_target_data(out)
    if n_bare_td:
        stats["bare_target_data_dropped"] += n_bare_td
        changed += n_bare_td
    out, n_dedup = dedup_declare_target(out)
    if n_dedup:
        stats["declare_target_deduped"] += n_dedup
        changed += n_dedup
    if changed == 0:
        return 0
    if write:
        path.write_text("".join(out))
    return changed


def _file_task(path: Path, write: bool = False,
               target: str = "gpu") -> tuple[Path, int, Counter, list[str]]:
    """Whole-file unit of work, safe to run in a worker process.

    `process_file` accumulates into caller-owned `stats`/`warnings`, so give it
    private ones and hand them back for the parent to merge — a shared Counter
    cannot cross a process boundary.
    """
    stats: Counter = Counter()
    warnings: list[str] = []
    n = process_file(path, write, stats, warnings, target)
    return path, n, stats, warnings


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
    ap.add_argument("--write", action="store_true", help="apply in place (default: dry-run)")
    ap.add_argument(
        "--target", choices=["gpu", "cpu"], default="gpu",
        help="gpu (default): !$omp target teams distribute parallel do + map() "
             "data directives. cpu: !$omp parallel do (host multicore); all "
             "device-data directives become inert no-op comments.")
    ap.add_argument(
        "-j", "--jobs", "-t", dest="jobs", type=int, default=1, metavar="N",
        help="translate files across N worker processes (0 = one per available "
             "core). Files are independent; output order is unchanged.")
    ap.add_argument("paths", nargs="*", default=["src"])
    args = ap.parse_args()

    roots = [Path(p) for p in args.paths]
    for r in roots:
        if not r.exists():
            print(f"error: {r} does not exist", file=sys.stderr)
            return 1

    stats: Counter = Counter()
    warnings: list[str] = []
    files_changed = 0
    jobs = resolve_jobs(args.jobs)
    results = pmap(partial(_file_task, write=args.write, target=args.target),
                   walk(roots), jobs, "acc_to_omp")
    for path, n, fstats, fwarn in results:
        stats.update(fstats)
        warnings.extend(fwarn)
        if n:
            files_changed += 1
            verb = "rewrote" if args.write else "would rewrite"
            print(f"  {verb} {n:3d}  {path}")

    print()
    print("=== conversion counts ===")
    order = ("atomic_left_as_is",
             "enter_data", "exit_data", "data", "end_data", "update",
             "parallel_loop", "end_parallel_loop", "host_data", "end_host_data",
             "declare", "set_device", "routine_seq",
             "clause_data", "clause_update", "continuation",
             "present_stripped", "if_present_stripped",
             "graph_wrapper_dropped", "async_skipped",
             "cpu_dropped", "cpu_clause_stripped", "acc_loop_clause_stripped")
    # Anything transformed but not listed above would otherwise vanish from
    # this report while still changing the source -- which is how 513 atomic
    # conversions once went unreported. Name the gap rather than hide it.
    unlisted = sorted(k for k in stats if k not in order)
    if unlisted:
        print("  (untallied keys, please add them to `order`: "
              + ", ".join(unlisted) + ")")

    for k in order:
        v = stats.get(k, 0)
        if v:
            print(f"  {k:<20} {v:5d}")
    if warnings:
        print()
        print("=== warnings ===")
        for w in warnings:
            print(w)
    mode = "WRITE" if args.write else "DRY RUN"
    print(f"\n[{mode}] {sum(stats.values())} transformations, {files_changed} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
