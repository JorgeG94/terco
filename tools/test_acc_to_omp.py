"""Unit tests for the ACC->OpenMP translator's loop-clause handling.

Run with: pytest tools/test_acc_to_omp.py

Focus: OpenACC-only loop-mapping clauses (gang/vector/worker/num_gangs/
vector_length/...) must be dropped from the converted loop directive — they are
invalid on the OpenMP combined-loop construct (GPU `target teams distribute
parallel do`, host `parallel do`).  Regression guard for the ifx #5082 the
naive head-swap produced (`!$omp target teams distribute parallel do gang
vector`).
"""
import importlib.util
import pathlib
from collections import Counter

_SPEC = importlib.util.spec_from_file_location(
    "acc_to_omp", pathlib.Path(__file__).parent / "acc_to_omp.py")
acc_to_omp = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(acc_to_omp)


def _xform(line, target):
    return "".join(acc_to_omp.transform_block(
        [line], Counter(), [], pathlib.Path("t.F90"), 1, target))


def test_gang_vector_dropped_gpu():
    out = _xform("      !$acc parallel loop gang vector\n", "gpu")
    assert "!$omp target teams distribute parallel do" in out
    assert "gang" not in out and "vector" not in out


def test_gang_vector_dropped_cpu():
    out = _xform("      !$acc parallel loop gang vector\n", "cpu")
    assert "!$omp parallel do" in out
    assert "gang" not in out and "vector" not in out


def test_paren_loop_clauses_dropped():
    out = _xform("      !$acc parallel loop gang vector_length(128) num_gangs(8)\n", "gpu")
    assert "vector_length" not in out and "num_gangs" not in out and "gang" not in out
    assert "!$omp target teams distribute parallel do" in out


def test_collapse_and_reduction_preserved():
    out = _xform("      !$acc parallel loop collapse(2) reduction(+:rr)\n", "gpu")
    assert "collapse(2)" in out
    assert "reduction(+:rr)" in out


def test_bare_keyword_inside_paren_not_stripped():
    # A variable named like an ACC loop keyword inside another clause's parens
    # (paren depth > 0) must survive.
    joined, _ = acc_to_omp.strip_acc_loop_clauses(
        "!$omp target teams distribute parallel do reduction(max:gang) private(vector)")
    assert "reduction(max:gang)" in joined
    assert "private(vector)" in joined


def test_present_still_stripped_with_gang():
    out = _xform("      !$acc parallel loop present(ms) gang vector\n", "gpu")
    assert "present(" not in out
    assert "gang" not in out and "vector" not in out


def test_async_parallel_loop_left_untranslated():
    # OpenACC async() queues have no OpenMP-target equivalent; a `parallel loop`
    # carrying async(...) must be left as inert `!$acc` (compiles + runs
    # synchronously under an OMP compiler), NOT head-swapped to an `!$omp`
    # construct with an invalid `async` clause. Both targets.
    src = "      !$acc parallel loop collapse(2) reduction(+:rr) async(cg_q)\n"
    for target in ("gpu", "cpu"):
        out = _xform(src, target)
        assert "!$omp" not in out
        assert "async(cg_q)" in out and out.lstrip().startswith("!$acc")


def test_async_kernels_still_dropped():
    # The `!$acc kernels async(...)` CUDA-graph wrapper stays neutralised by the
    # earlier graph-wrapper drop — the async-skip must not resurrect it.
    out = _xform("         !$acc kernels async(cg_q)\n", "gpu")
    assert "!$acc" not in out and "!$omp" not in out
    assert "dropped CUDA-graph wrapper" in out


def test_duplicate_declare_target_deduped(tmp_path):
    # Source carries `!$acc routine seq` (OpenACC build) AND a hand-written
    # `!$omp declare target` (direct OpenMP build of main). Converting the
    # former must NOT leave two declare targets (ifx #8682).
    src = ("   pure subroutine rescale(x)\n"
           "      !$acc routine seq\n"
           "      !$omp declare target\n"
           "      !! a docstring between, must not break the dedup\n"
           "      real :: x\n"
           "   end subroutine rescale\n")
    p = tmp_path / "m.F90"
    p.write_text(src)
    acc_to_omp.process_file(p, True, Counter(), [], "gpu")
    out = p.read_text()
    assert out.count("!$omp declare target") == 1
    assert "!$acc routine seq" not in out


def test_declare_targets_in_separate_procedures_kept(tmp_path):
    # Two different procedures each get one declare target — the run must reset
    # at the procedure boundary, so BOTH survive.
    src = ("   pure subroutine a(x)\n"
           "      !$acc routine seq\n"
           "      real :: x\n"
           "   end subroutine a\n"
           "   pure subroutine b(y)\n"
           "      !$acc routine seq\n"
           "      real :: y\n"
           "   end subroutine b\n")
    p = tmp_path / "m.F90"
    p.write_text(src)
    acc_to_omp.process_file(p, True, Counter(), [], "gpu")
    out = p.read_text()
    assert out.count("!$omp declare target") == 2
