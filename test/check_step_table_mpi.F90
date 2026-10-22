module check_step_table_mpi_suite
   !! step_table_t%reduce over the ranks of a run, on however many ranks it
   !! is launched with: a step's time becomes the slowest rank's, its
   !! operations and bytes the sum over the ranks, and a replicated step's
   !! work (one every rank runs in full, like RI-MP2's metric Cholesky)
   !! counts once. Without it RI-MP2's table was rank 0's share, so its
   !! FLOP/s stayed flat when ranks were added.
   use pic_types, only: dp, int64
   use pic_mpi_lib, only: comm_t
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_step_table, only: step_table_t
   implicit none
   private
   public :: collect_step_table_mpi, suite_comm

   type(comm_t) :: suite_comm
   !! The run's communicator, set by the program before the suite runs.

   integer, parameter :: ST_SHARED = 1, ST_OWN = 2, ST_BYTES = 3

contains

   subroutine collect_step_table_mpi(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("time is the slowest rank's", test_time), &
                  new_unittest("distributed work is summed", test_summed), &
                  new_unittest("replicated work counts once", test_replicated), &
                  new_unittest("bytes are summed", test_bytes), &
                  new_unittest("rank count is recorded", test_nranks) &
                  ]
   end subroutine collect_step_table_mpi

   subroutine filled(tbl)
      !! Rank r ran the replicated step for 1 s and 1000 FLOP, its own step
      !! for 1 + r s and 100 (r + 1) FLOP, and moved 8 (r + 1) bytes.
      type(step_table_t), intent(out) :: tbl
      integer :: r
      r = suite_comm%rank()
      call tbl%init("t", [character(len=8) :: "shared", "own", "copy"], bytes_steps=[ST_BYTES], &
                    replicated_steps=[ST_SHARED])
      call tbl%add_time(ST_SHARED, 1.0_dp)
      call tbl%add_flop(ST_SHARED, 1000_int64)
      call tbl%add_time(ST_OWN, real(1 + r, dp))
      call tbl%add_flop(ST_OWN, 100_int64*(r + 1))
      call tbl%add_bytes(ST_BYTES, 8_int64*(r + 1))
      call tbl%reduce(suite_comm)
   end subroutine filled

   subroutine test_time(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      call filled(tbl)
      call check(error, tbl%t(ST_OWN), real(suite_comm%size(), dp))
      if (allocated(error)) return
      call check(error, tbl%t(ST_SHARED), 1.0_dp)
   end subroutine test_time

   subroutine test_summed(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      integer(int64) :: n
      call filled(tbl)
      n = suite_comm%size()
      call check(error, tbl%fl(ST_OWN)%get(), 50_int64*n*(n + 1))
   end subroutine test_summed

   subroutine test_replicated(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      call filled(tbl)
      call check(error, tbl%fl(ST_SHARED)%get(), 1000_int64)
   end subroutine test_replicated

   subroutine test_bytes(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      integer(int64) :: n
      call filled(tbl)
      n = suite_comm%size()
      call check(error, tbl%bytes(ST_BYTES), 4_int64*n*(n + 1))
   end subroutine test_bytes

   subroutine test_nranks(error)
      type(error_type), allocatable, intent(out) :: error
      type(step_table_t) :: tbl
      call filled(tbl)
      call check(error, tbl%nranks_reduced, suite_comm%size())
   end subroutine test_nranks

end module check_step_table_mpi_suite

program check_step_table_mpi
   use, intrinsic :: iso_fortran_env, only: error_unit
   use pic_mpi_lib, only: comm_world, pic_mpi_init, pic_mpi_finalize, allreduce, MPI_MAX
   use testdrive, only: run_testsuite
   use check_step_table_mpi_suite, only: collect_step_table_mpi, suite_comm
   implicit none
   integer :: stat
   call pic_mpi_init()
   suite_comm = comm_world()
   stat = 0
   ! Serial: every case is a collective, so every rank must run the cases in
   ! the same order, one at a time.
   call run_testsuite(collect_step_table_mpi, error_unit, stat, parallel=.false.)
   call allreduce(suite_comm, stat, op=MPI_MAX)
   call pic_mpi_finalize()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed on some rank"
      error stop 1
   end if
end program check_step_table_mpi
