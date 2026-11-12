!
! The collective calls of multi-rank RI-CCSD(T), in one place (module trc_rccsd_driver and trc_rccsd_triples
! describe the scheme they serve). Every routine is collective over `comm`: all of its ranks call it, in the same
! order. None is called without a communicator of more than one rank, so a one-rank run makes no MPI call.
!
! The per-iteration transfers (the amplitudes down, the partial ladders up) are in trc_rccsd_xfer; this module has the
! small collectives and the host broadcast that ends a run.
module trc_rccsd_mpi
   use pic_types, only: dp, int32
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_MIN
   implicit none
   private

   public :: rccsd_mpi_agree_text, rccsd_mpi_bcast_amps, rccsd_mpi_bcast_orbitals
   public :: rccsd_mpi_control, rccsd_mpi_gather_values

   integer, parameter, public :: RCCSD_MPI_TEXT_LEN = 300  !! length of the texts the ranks agree on (the result's message)
   integer, parameter :: BCAST_CHUNK = 33554432  !! doubles per broadcast call: 256 MiB, far below the 2^31-1 count limit

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
      !! Rank `root`'s host t1(no, nv) and t2(no, no, nv, nv) to every rank's, each as a flat array in chunks of
      !! BCAST_CHUNK doubles (an MPI count is a 32-bit integer; no reshaped copy). Collective. Sent a b-slice at a time
      !! (130 calls of 0.65 MB at adenine cc-pVDZ) it took 0.26 s for 84 MB on two V100s, in flat chunks 26 ms.
      type(comm_t), intent(in) :: comm
      real(dp), intent(inout) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: root

      call bcast_flat(comm, size(t1), t1, root)
      call bcast_flat(comm, size(t2), t2, root)
   end subroutine rccsd_mpi_bcast_amps

   subroutine bcast_flat(comm, n, x, root)
      !! x(1:n) of rank `root` to every rank, in chunks
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, root
      real(dp), intent(inout) :: x(n)

      integer :: k, m

      k = 1
      do while (k <= n)
         m = min(BCAST_CHUNK, n - k + 1)
         call bcast(comm, x(k:k + m - 1), int(m, int32), int(root, int32))
         k = k + m
      end do
   end subroutine bcast_flat

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

   function rccsd_mpi_control(comm, word) result(received)
      !! The first rank's control word to every rank (one integer): the iteration number to run, or 0 to stop.
      !! Collective; `word` is read on rank 0 only.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: word
      integer :: received

      integer(int32) :: w

      w = int(word, int32)
      call bcast(comm, w, 1_int32, 0_int32)
      received = int(w)
   end function rccsd_mpi_control

   function rccsd_mpi_gather_values(comm, value) result(values)
      !! A few numbers per rank, values(:, r + 1) of rank r, on every rank: each rank writes its own column of a zero
      !! array and the arrays are summed (pic-mpi has no gather of reals). Collective; for the timings the report shows
      !! once per run.
      type(comm_t), intent(in) :: comm
      real(dp), intent(in) :: value(:)
      real(dp) :: values(size(value), comm%size())

      real(dp) :: flat(size(value)*comm%size())

      flat = 0.0_dp
      flat(comm%rank()*size(value) + 1:(comm%rank() + 1)*size(value)) = value
      call allreduce(comm, flat)
      values = reshape(flat, shape(values))
   end function rccsd_mpi_gather_values

end module trc_rccsd_mpi
