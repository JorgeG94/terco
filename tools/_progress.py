"""Tiny stdlib-only progress bar for the source-translation tools.

Used by `dc_audit.py`, `acc_to_omp.py` and `dc_to_omp.py`, which each sweep the
whole tree (~1500 files) and, in `dc_to_omp`'s case, run an unbounded
whole-program fixpoint — long enough that a silent terminal is indistinguishable
from a hang.

Deliberately dependency-free (no tqdm): the repo rule is to use only what the
loaded module environment already provides.

Rendering goes to stderr so it never contaminates a redirected stdout report,
and is suppressed entirely when stderr is not a TTY, keeping CI logs clean.
"""
from __future__ import annotations

import os
import sys
import time
from typing import Iterable, Iterator, Optional, TypeVar

T = TypeVar("T")

_BAR_W = 28
_MIN_REDRAW = 0.05  # s — cap redraws at ~20 Hz so rendering is never the cost


def _enabled(stream) -> bool:
    """Render only to an interactive terminal, and honour NO_COLOR-style opt-out."""
    if os.environ.get("RAKALI_NO_PROGRESS"):
        return False
    try:
        return stream.isatty()
    except Exception:
        return False


def _fmt_eta(seconds: float) -> str:
    if seconds < 0 or seconds != seconds or seconds == float("inf"):
        return "--:--"
    m, s = divmod(int(seconds), 60)
    if m >= 60:
        h, m = divmod(m, 60)
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:02d}:{s:02d}"


class Progress:
    """Render a bar that coexists with normal `print()` output.

    The line is cleared before each item is handed to the caller, so anything
    the loop body prints lands on a clean line; the bar is then redrawn beneath
    it. That keeps the existing per-file "rewrote N path" reporting intact.
    """

    def __init__(self, label: str, total: Optional[int] = None, stream=None):
        self.label = label
        self.total = total
        self.stream = stream if stream is not None else sys.stderr
        self.on = _enabled(self.stream)
        self.n = 0
        self.start = time.monotonic()
        self._last_draw = 0.0
        self._dirty = False

    # -- line handling -----------------------------------------------------
    def clear(self) -> None:
        if self.on and self._dirty:
            self.stream.write("\r\x1b[K")
            self.stream.flush()
            self._dirty = False

    def draw(self, force: bool = False, suffix: str = "") -> None:
        if not self.on:
            return
        now = time.monotonic()
        if not force and (now - self._last_draw) < _MIN_REDRAW:
            return
        self._last_draw = now
        elapsed = now - self.start
        if self.total:
            frac = min(1.0, self.n / self.total)
            filled = int(_BAR_W * frac)
            bar = "#" * filled + "-" * (_BAR_W - filled)
            rate = self.n / elapsed if elapsed > 0 else 0.0
            eta = (self.total - self.n) / rate if rate > 0 else float("inf")
            msg = (f"\r\x1b[K  {self.label} [{bar}] {self.n}/{self.total} "
                   f"({100 * frac:4.1f}%) eta {_fmt_eta(eta)}")
        else:
            # Unknown total (fixpoint sweeps): spinner + running count.
            spin = "|/-\\"[int(now * 8) % 4]
            msg = f"\r\x1b[K  {self.label} {spin} {self.n} in {_fmt_eta(elapsed)}"
        if suffix:
            msg += f"  {suffix}"
        self.stream.write(msg)
        self.stream.flush()
        self._dirty = True

    def advance(self, k: int = 1, suffix: str = "") -> None:
        self.n += k
        self.draw(suffix=suffix)

    def done(self, note: str = "") -> None:
        if not self.on:
            return
        self.clear()
        elapsed = time.monotonic() - self.start
        tail = f"  {note}" if note else ""
        self.stream.write(f"  {self.label}: {self.n} in {_fmt_eta(elapsed)}{tail}\n")
        self.stream.flush()


def track(it: Iterable[T], label: str, total: Optional[int] = None,
          stream=None) -> Iterator[T]:
    """Wrap an iterable in a progress bar.

    `total` is inferred with len() when the iterable is sized.
    """
    if total is None:
        try:
            total = len(it)  # type: ignore[arg-type]
        except TypeError:
            total = None
    p = Progress(label, total, stream)
    p.draw(force=True)
    for item in it:
        p.clear()
        yield item
        p.advance()
    p.done()
