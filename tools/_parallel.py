"""Optional multi-process fan-out for the per-file translation passes.

Processes, not threads, on purpose: these passes are pure-Python `re` work and
CPython's `re` does not release the GIL, so a thread pool would serialise and
buy nothing. A Frontier compute node has 64 physical cores / 128 hardware
threads, which a process pool can actually use.

`pmap` preserves INPUT ORDER regardless of completion order (it consumes
`Executor.map`, not `as_completed`). That matters: these tools emit a report and
a fixpoint seeded from the per-file results, and a source-to-source translator
whose output depends on scheduling would be untrustworthy.

Only whole-file, independent work belongs here — each task must read and write
its own file and share no mutable state with its siblings.
"""
from __future__ import annotations

import os
from typing import Callable, Iterable, List, Optional, Sequence, TypeVar

from _progress import Progress

A = TypeVar("A")
B = TypeVar("B")

# Below this many items the pool setup + pickling costs more than it saves.
_MIN_PARALLEL_ITEMS = 8


def default_jobs() -> int:
    """Cores available to THIS process (respects cpuset/SLURM pinning)."""
    try:
        return len(os.sched_getaffinity(0))  # Linux: honours SLURM binding
    except AttributeError:
        return os.cpu_count() or 1


def resolve_jobs(requested: Optional[int]) -> int:
    """Map a --jobs value onto a worker count. 0 or None -> auto-detect."""
    if not requested:
        return default_jobs()
    if requested < 0:
        return 1
    return requested


def pmap(fn: Callable[[A], B], items: Iterable[A], jobs: int = 1,
         label: str = "work") -> List[B]:
    """Apply `fn` to every item, optionally across processes, with a bar.

    `fn` must be picklable (a module-level function, or a `functools.partial`
    of one) and must not rely on state shared with other tasks.
    """
    seq: Sequence[A] = list(items)
    p = Progress(label, len(seq))
    p.draw(force=True)
    out: List[B] = []

    if jobs <= 1 or len(seq) < _MIN_PARALLEL_ITEMS:
        for it in seq:
            p.clear()
            out.append(fn(it))
            p.advance()
        p.done()
        return out

    # Chunk so each worker gets a decent batch of files but the bar still moves
    # often enough to be informative.
    chunksize = max(1, len(seq) // (jobs * 8))
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for r in ex.map(fn, seq, chunksize=chunksize):
            p.clear()
            out.append(r)
            p.advance()
    p.done(f"({jobs} procs)")
    return out
