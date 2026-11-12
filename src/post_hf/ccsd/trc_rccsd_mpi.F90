!
! The collective calls of multi-rank RI-CCSD(T), in one place (module trc_rccsd_driver and trc_rccsd_triples
! describe the scheme they serve). Every routine is collective over `comm`: all of its ranks call it, in the same
! order. None is called without a communicator of more than one rank, so a one-rank run makes no MPI call.
!
! Host buffers throughout. The amplitudes are sent once per run (no^2 nv^2 doubles), and a host broadcast is the
! same code on every MPI; RI-MP2's ring, which moves a block per step, is where handing MPI device buffers pays
! (trc_mpi_device).
module trc_rccsd_mpi
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_MIN
   implicit none
   private

   public :: rccsd_mpi_agree_text, rccsd_mpi_bcast_amps, rccsd_mpi_bcast_orbitals

contains

   function rccsd_mpi_agree_text(comm, text, failed) result(any_failed)
      !! True on every rank when `failed` is true on at least one. Then `text` on EVERY rank becomes the text of the
      !! lowest-numbered failing rank (MIN of the failing rank index, then a broadcast of its text), so a reason is never
      !! lost to a rank that merely heard of the failure. A failing rank must have a non-blank text. Collective; the
      !! two calls inside are made on every rank in the same order.
      type(comm_t), intent(in) :: comm
      character(len=*), intent(inout) :: text
      logical, intent(in) :: failed
      logical :: any_failed

      integer(int32) :: first
      real(dp), allocatable :: codes(:)
      integer :: i

      first = merge(int(comm%rank(), int32), int(comm%size(), int32), failed)
      call allreduce(comm, first, op=MPI_MIN)
      any_failed = first < int(comm%size(), int32)
      if (.not. any_failed) return
      allocate (codes(len(text)))
      do i = 1, len(text)
         codes(i) = real(iachar(text(i:i)), dp)
      end do
      call bcast(comm, codes, int(len(text), int32), first)
      do i = 1, len(text)
         text(i:i) = achar(nint(codes(i)))
      end do
      deallocate (codes)
   end function rccsd_mpi_agree_text

   subroutine rccsd_mpi_bcast_amps(comm, t1, t2, root)
      !! Rank `root`'s host t1(no, nv) and t2(no, no, nv, nv) to every rank's. t2 goes one b-slice (no^2 nv doubles)
      !! at a time through a flat buffer: the array is rank 4 and pic's broadcast takes a rank-1 buffer, and a slice is
      !! far below the 2^31 element limit of an MPI count at any size this method reaches.
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: root

      real(dp), allocatable :: buffer(:)
      integer :: b, n1, n2

      n1 = size(t1)
      allocate (buffer(max(n1, size(t2, 1)*size(t2, 2)*size(t2, 3))))
      buffer(1:n1) = reshape(t1, [n1])
      call bcast(comm, buffer, int(n1, int32), int(root, int32))
      t1 = reshape(buffer(1:n1), shape(t1))
      n2 = size(t2, 1)*size(t2, 2)*size(t2, 3)
      do b = 1, size(t2, 4)
         buffer(1:n2) = reshape(t2(:, :, :, b), [n2])
         call bcast(comm, buffer, int(n2, int32), int(root, int32))
         t2(:, :, :, b) = reshape(buffer(1:n2), shape(t2(:, :, :, b)))
      end do
      deallocate (buffer)
   end subroutine rccsd_mpi_bcast_amps

   subroutine rccsd_mpi_bcast_orbitals(comm, cmo, eps, root)
      !! Rank `root`'s orbitals and orbital energies to every rank's: a GUARD. The collective SCF of the CLI and the
      !! API already gives every rank the same ones; this protects a caller who hands ranks different orbitals (SCFs
      !! run separately agree to the threshold but can differ in the sign of a column, and amplitudes from one rank's
      !! orbitals are only amplitudes in that rank's integrals: the (T) share of another rank came out wrong by 5e-5 Eh
      !! on water). One nmo^2 message.
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: cmo(:, :), eps(:)
      integer, intent(in) :: root

      real(dp), allocatable :: buffer(:)
      integer :: j

      allocate (buffer(max(size(cmo, 1), size(eps))))
      buffer(1:size(eps)) = eps
      call bcast(comm, buffer, int(size(eps), int32), int(root, int32))
      eps = buffer(1:size(eps))
      do j = 1, size(cmo, 2)
         buffer(1:size(cmo, 1)) = cmo(:, j)
         call bcast(comm, buffer, int(size(cmo, 1), int32), int(root, int32))
         cmo(:, j) = buffer(1:size(cmo, 1))
      end do
      deallocate (buffer)
   end subroutine rccsd_mpi_bcast_orbitals

end module trc_rccsd_mpi
