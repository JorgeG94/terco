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


def test_pointer_names():
    f = archcheck.pointer_names
    assert f("type(trc_linalg_t), pointer :: lp") == ["lp"]
    assert f("real(dp), pointer, contiguous :: fxc(:) => null()") == ["fxc"]
    assert f("procedure(cb_t), pointer :: cb => null()") == ["cb"]
    assert f("pointer :: a, b(:)") == ["a", "b"]
    assert f("type(json_value), pointer :: root") == []
    assert f("type(json_value), pointer, intent(in) :: root") == []
    assert f("real(dp), allocatable :: pointer_free(:)") == []
    assert f("logical :: has_pointer = .false.") == []
    assert f("call c_f_pointer(handle, box)") == []


def test_bad_words_matches_whole_words_only():
    assert archcheck.bad_words("dens_a") == ["dens"]
    assert archcheck.bad_words("density") == []
    assert archcheck.bad_words("ctx") == ["ctx"]
    assert archcheck.bad_words("idx_map") == ["idx"]
    assert archcheck.bad_words("matrix") == []


def test_cli_layer_dependencies(tmp_path, monkeypatch):
    for layer in ("core", "cli"):
        (tmp_path / "src" / layer).mkdir(parents=True)
    (tmp_path / "app").mkdir()
    (tmp_path / "src/core/trc_core.F90").write_text(
        "module trc_core\nuse trc_cli\nend module trc_core\n")
    (tmp_path / "src/cli/trc_cli.F90").write_text(
        "module trc_cli\nuse trc_core\nend module trc_cli\n")
    monkeypatch.setattr(archcheck, "ROOT", tmp_path)
    monkeypatch.setattr(archcheck.modmap, "ROOT", tmp_path)

    violations = [v for v in archcheck.check() if v.rule == "layer-up"]
    assert [(v.scope, v.detail) for v in violations] == [("trc_core", "trc_cli")]


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


def test_procedure_length_excludes_only_initial_declarations(tmp_path):
    src = tmp_path / "length.F90"
    src.write_text(textwrap.dedent("""\
        module trc_demo
        contains
           subroutine run(a)
              use pic_types, only: dp
              real(dp), intent(in) :: &
                 ! comment inside a continued declaration
                 a
              integer :: i
              !$acc routine seq
              ! a comment
              block
                 integer :: j
                 j = 1
              end block
              i = 1
           end subroutine run
        end module trc_demo
    """))
    proc = next(p for p in archcheck.scan_procedures(src) if p["name"] == "run")
    assert archcheck.procedure_length(src.read_text().splitlines(), proc) == 9


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


def test_bind_c_args_flagged_free(tmp_path):
    # bind(c) dummy names are the C ABI's parameter names: scanned, but
    # marked so short-name leaves them alone.
    src = tmp_path / "c.F90"
    src.write_text(textwrap.dedent("""\
        module trc_c
           use iso_c_binding
           implicit none
        contains
           function trc_c_get(buf, max_iter) result(status) bind(c, name="trc_c_get")
              character(kind=c_char), intent(out) :: buf(*)
              integer(c_int), value :: max_iter
              integer(c_int) :: status
              status = max_iter
           end function trc_c_get
        end module trc_c
    """))
    procs = {p["name"]: p for p in archcheck.scan_procedures(src)}
    assert procs["trc_c_get"]["bind_c"]
    assert set(procs["trc_c_get"]["args"]) == {"buf", "max_iter"}
