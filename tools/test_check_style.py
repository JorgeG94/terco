"""Unit tests for tools/check_style.py (run with pytest)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_style as cs  # noqa: E402

LIB = Path('src/x.F90')
TEST = Path('test/x.F90')


def rules(line, path=LIB):
    return [name for name, _ in cs.check_line(path, line)]


def test_system_clock():
    assert rules('      call system_clock(t0, rate)') == ['system_clock']
    assert rules('      call SYSTEM_CLOCK(t0)') == ['system_clock']
    assert rules('      ! was timed with system_clock') == []
    assert rules('      msg = "system_clock"') == []


def test_bare_mpi():
    assert rules('   use mpi_f08') == ['bare MPI']
    assert rules('   use mpi') == ['bare MPI']
    assert rules('   use, intrinsic :: mpi') == ['bare MPI']
    assert rules("   include 'mpif.h'") == ['bare MPI']
    assert rules('   call MPI_Allreduce(a, b, n, MPI_DOUBLE, MPI_SUM, comm, ierr)') == ['bare MPI']
    assert rules('   ierr = mpi_barrier(comm)') == ['bare MPI']


def test_pic_mpi_is_fine():
    assert rules('   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM') == []
    assert rules('   call allreduce(comm, x, MPI_SUM)') == []
    assert rules('   ! the handle comes from pic-mpi, not from `use mpi_f08`') == []


def test_doxygen():
    assert rules('   real(dp) :: x   !< the thing') == ['Doxygen']
    assert rules('   !> the module') == ['Doxygen']
    assert rules('   ! @param n the count') == ['Doxygen']
    assert rules('   ! @brief the thing') == ['Doxygen']
    assert rules('   !! a FORD block') == []
    assert rules('   !$acc enter data copyin(a)') == []
    assert rules('   !$omp parallel do') == []
    assert rules('   ! plain comment about a > b') == []
    assert rules('   s = "!> not a comment"') == []


def test_literal_kinds():
    assert rules('   real(8) :: x') == ['literal kind']
    assert rules('   real(kind=8) :: x') == ['literal kind']
    assert rules('   x = 1.0d0') == ['literal kind']
    assert rules('   x = 2.5D-3*y') == ['literal kind']
    assert rules('   real(dp) :: x') == []
    assert rules('   x = 1.0_dp') == []
    assert rules('   integer :: nd0') == []


def test_console_output_only_in_the_library():
    assert rules('   print *, x') == ['console output']
    assert rules("   print '(a)', msg") == ['console output']
    assert rules('   write(*, "(a)") msg') == ['console output']
    assert rules('   write(buf, "(i0)") n') == []
    assert rules("   call global_logger%info('x')") == []
    assert rules('   print *, x', TEST) == []
    assert rules('   print *, x', Path('app/main.F90')) == ['console output']


def test_the_tree_is_clean(tmp_path, monkeypatch):
    f = tmp_path / 'ok.F90'
    f.write_text('module m\n   use pic_mpi_lib, only: comm_t\ncontains\nend module m\n')
    assert cs.main([str(f)]) == 0
    g = tmp_path / 'bad.F90'
    g.write_text('program p\n   call system_clock(t)\nend program p\n')
    assert cs.main([str(g)]) == 1
