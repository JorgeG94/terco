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


def test_kernel_driver(tmp_path):
    # A thin launch is a kernel driver; one with real host logic (more than
    # KERNEL_HOST_STMTS statements outside its do concurrent) is not, and a
    # procedure with no do concurrent never is.
    host = "\n".join(f"      n = n + {k}" for k in range(archcheck.KERNEL_HOST_STMTS + 1))
    src = tmp_path / "k.F90"
    src.write_text(textwrap.dedent("""\
        module trc_k
           implicit none
        contains
           subroutine thin(n, a)
              integer, intent(in) :: n
              real, intent(inout) :: a(n)
              integer :: i
              !$acc enter data copyin(a)
              do concurrent(i=1:n)
                 a(i) = 2*a(i)
              end do
              !$acc exit data copyout(a)
           end subroutine thin
           subroutine busy(n, a)
              integer, intent(inout) :: n
              real, intent(inout) :: a(n)
              integer :: i
        """) + host + "\n" + textwrap.dedent("""\
              do concurrent(i=1:n)
                 a(i) = 2*a(i)
              end do
           end subroutine busy
           subroutine plain(n)
              integer, intent(inout) :: n
              n = n + 1
           end subroutine plain
        end module trc_k
    """))
    procs = {p["name"]: p for p in archcheck.scan_procedures(src)}
    assert archcheck.is_kernel_driver(procs["thin"])
    assert procs["thin"]["host_exec"] == 0
    assert not archcheck.is_kernel_driver(procs["busy"])
    assert not archcheck.is_kernel_driver(procs["plain"])
