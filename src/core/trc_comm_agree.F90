!
! One agreement on a failure across the ranks of a communicator.
!
! A collective routine whose ranks can fail on their own (a refusal that only
! rank 0 decides, an allocation, a build) must not let the failing rank return
! while the others sit in the next collective call. Every rank calls
! comm_agree_failure with its own status and message BEFORE that call; all of
! them come back with the same pair, the one of the lowest-numbered failing
! rank, or with status 0 untouched when no rank failed.
!
! No character broadcast is needed: the failing rank with the lowest number
! is found by a MIN reduction of rank numbers, and then only that rank
! contributes its text, one character code per integer, to a SUM reduction.
!
module trc_comm_agree
   use pic_types, only: int32
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM, MPI_MIN
   implicit none
   private
   public :: comm_agree_failure, COMM_AGREE_MSG_LEN

   integer, parameter :: COMM_AGREE_MSG_LEN = 1024
   !! Longest message that survives the agreement; a longer one is cut.

contains

   subroutine comm_agree_failure(comm, status, message)
      !! Collective over `comm`. `status` 0 means this rank succeeded. On
      !! return `status` and `message` are those of the lowest-numbered rank
      !! that failed (the same on every rank), or unchanged and 0 when none did.
      type(comm_t), intent(in) :: comm
      integer(int32), intent(inout) :: status
      character(len=:), allocatable, intent(inout) :: message
      integer(int32) :: first, nranks
      integer(int32), allocatable :: codes(:)
      integer :: i, n
      nranks = int(comm%size(), int32)
      first = nranks
      if (status /= 0_int32) first = int(comm%rank(), int32)
      call allreduce(comm, first, op=MPI_MIN)
      if (first >= nranks) return      ! nobody failed
      allocate (codes(1 + COMM_AGREE_MSG_LEN))
      codes = 0_int32
      if (int(comm%rank(), int32) == first) then
         codes(1) = status
         n = 0
         if (allocated(message)) n = min(len(message), COMM_AGREE_MSG_LEN)
         do i = 1, n
            codes(1 + i) = int(iachar(message(i:i)), int32)
         end do
      end if
      call allreduce(comm, codes, op=MPI_SUM)
      status = codes(1)
      n = 0
      do i = 1, COMM_AGREE_MSG_LEN
         if (codes(1 + i) == 0_int32) exit
         n = i
      end do
      message = repeat(" ", n)
      do i = 1, n
         message(i:i) = achar(codes(1 + i))
      end do
      deallocate (codes)
   end subroutine comm_agree_failure

end module trc_comm_agree
