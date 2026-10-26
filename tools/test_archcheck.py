"""Unit tests for tools/archcheck.py (run with pytest)."""
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import archcheck  # noqa: E402


def test_declared_names():
    f = archcheck.declared_names
    assert f("real(dp), intent(in) :: dens(n, m), thr = 1.0") == ["dens", "thr"]
    assert f("type(foo_t), allocatable :: ctx") == ["ctx"]
    assert f("integer :: i, j") == ["i", "j"]
    assert f("type :: thing_t") == []
    assert f("call foo(a, b)") == []


def test_bad_words_matches_whole_words_only():
    assert archcheck.bad_words("dens_a") == ["dens"]
    assert archcheck.bad_words("density") == []
    assert archcheck.bad_words("ctx") == ["ctx"]
    assert archcheck.bad_words("idx_map") == ["idx"]
    assert archcheck.bad_words("matrix") == []


def test_scan_procedures(tmp_path):
    src = tmp_path / "m.F90"
    src.write_text(textwrap.dedent("""\
        module trc_demo
           implicit none
           integer :: tmp_global
        contains
           subroutine demo_run(a, b, &
                               c)
              integer, intent(in) :: a, b, c
              real :: ctx
              ctx = a + b + c
           end subroutine demo_run
           pure function demo_twice(x) result(y)
              !$acc routine seq
              real, intent(in) :: x
              real :: y
              y = 2*x
           end function demo_twice
        end module trc_demo
    """))
    procs = {p["name"]: p for p in archcheck.scan_procedures(src)}
    run = procs["demo_run"]
    assert run["args"] == ["a", "b", "c"]
    assert (run["end"] - run["start"] + 1) == 6
    assert ("ctx", 8) in run["names"]
    assert procs["demo_twice"]["device"]
    assert ("tmp_global", 3) in procs["trc_demo"]["names"]


def test_scan_skips_plain_interface_bodies(tmp_path):
    # A plain interface block declares external procedures (C bindings,
    # LAPACK): their names and argument counts are not the module's to set.
    src = tmp_path / "m.F90"
    src.write_text(textwrap.dedent("""\
        module trc_demo
           implicit none
           interface
              function cDgemm(h, ta, tb, m, n, k, al, a, lda, b, ldb, be, c, ldc) &
                 result(st) bind(c, name="cDgemm")
                 integer :: h, ta, tb, m, n, k, lda, ldb, ldc, st
                 real :: al, a, b, be, c
              end function cDgemm
           end interface
        contains
           subroutine demo_run(a)
              integer, intent(in) :: a
           end subroutine demo_run
        end module trc_demo
    """))
    names = {p["name"] for p in archcheck.scan_procedures(src)}
    assert "cdgemm" not in names
    assert "demo_run" in names
