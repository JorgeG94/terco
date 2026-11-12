!
! The two transfers of a CCSD iteration over ranks (trc_rccsd_driver, RANKS): the amplitudes from the first rank to
! each other rank, and each other rank's partial ladder to the first rank. Point to point, o^2 v^2 doubles at most,
! never an all-ranks reduction: a follower sends its partial once, the first rank receives it into a staging array
! and adds it to its own with a kernel on the device, nothing comes back, and a rank whose share is empty sends
! nothing.
!
! PAIR SYMMETRY HALVES THE DATA. When t2(i,j,a,b) = t2(j,i,b,a) (every production run: `tau_symmetric`, which the
! control word carries) only the unique half is sent: with p = (i,a), q = (j,b), the N (N + 1) / 2 elements p <= q,
! N = o v, at the packed position p + q (q - 1) / 2. The ladder of the half form's one-contraction path has the
! same symmetry exactly (h = N + N^T, a sum that is the same in either order). trc_rccsd_diis has the same pairing
! but keeps its pack/unpack kernels private and averages the two halves, which a lossless transfer must not do, so the
! pair is written here. The unpacked amplitudes are symmetric to the bit; the ladder of the half form
! symmetrises tau when it gathers it, so a rank whose t2 was symmetric only to round-off computes the same
! ladder from the unpacked one. The general (non-symmetric) case sends the whole arrays, from the arrays themselves.
!
! TWO PATHS, ONE SWITCH. `trc_mpi_device_direct()` (trc_mpi_device: the build's default, TERCO_MPI_DEVICE=0/1 in the
! environment overrides it) decides, read once by `init`. Direct: the device array goes to MPI inside
! `host_data use_device`, by pic-mpi's specific names and in the chunks RI-MP2's ring uses (trc_rimp2_driver:
! comm_isend_real_dp_array_n, comm_recv_real_dp_array_n, MSG_CHUNK doubles), so nothing passes through the host.
! Host: `update self` before the send and `update device` after the receive, the generic calls. A build without a
! device has only the second, and its `update` and `host_data` lines are no-ops.
!
! The t1 (o v doubles) always goes through the host, in one message.
module trc_rccsd_xfer
   use pic_types, only: dp, int32, int64
   use pic_mpi_lib, only: comm_t, isend, recv, request_t, waitall, comm_isend_real_dp_array_n, comm_recv_real_dp_array_n
   use trc_mpi_device, only: trc_mpi_device_direct
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   implicit none
   private

   integer(int64), parameter :: MSG_CHUNK = 2_int64**30  !! doubles per message, as RI-MP2's ring (MPI counts are 32-bit)
   integer, parameter :: TAG_T2 = 31, TAG_T1 = 32, TAG_LADDER = 33

   public :: trc_rccsd_xfer_t
   public :: rccsd_xfer_amps_down, rccsd_xfer_ladder_up, rccsd_xfer_ladder_in
   public :: rccsd_xfer_amps_bytes, rccsd_xfer_ladder_bytes

   type :: trc_rccsd_xfer_t
      !! The staging array of the packed half, on the host and the device
      integer :: no = 0
      integer :: nv = 0
      integer :: nhalf = 0               !! N (N + 1) / 2, N = o v
      logical :: direct = .false.        !! MPI is handed device addresses (module header)
      real(dp), allocatable :: stage(:)  !! nhalf doubles
   contains
      procedure :: init => xfer_init
      procedure :: release => xfer_release
   end type trc_rccsd_xfer_t

contains

   subroutine xfer_init(this, no, nv, ok)
      !! Allocates the staging array for (no, nv) and reads the path switch; `ok` false when it does not fit a default
      !! integer index or the allocation failed (nothing is left allocated)
      class(trc_rccsd_xfer_t), intent(inout) :: this
      integer, intent(in) :: no, nv
      logical, intent(out) :: ok

      integer(int64) :: n
      integer :: stat

      call this%release()
      ok = .false.
      n = int(no, int64)*nv
      if (n*(n + 1)/2 > int(huge(0_int32), int64)) return
      this%no = no
      this%nv = nv
      this%nhalf = int(n*(n + 1)/2)
      this%direct = trc_mpi_device_direct()
      allocate (this%stage(this%nhalf), stat=stat)
      if (stat /= 0) then
         call this%release()
         return
      end if
      !$acc enter data create(this%stage)
      ok = .true.
   end subroutine xfer_init

   subroutine xfer_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_xfer_t), intent(inout) :: this

      if (allocated(this%stage)) then
         !$acc exit data delete(this%stage)
         deallocate (this%stage)
      end if
      this%no = 0
      this%nv = 0
      this%nhalf = 0
   end subroutine xfer_release

   pure function rccsd_xfer_amps_bytes(no, nv, symmetric) result(bytes)
      !! Bytes of one copy of the amplitudes as sent (t1 whole, t2 whole or its unique half)
      integer, intent(in) :: no, nv
      logical, intent(in) :: symmetric
      integer(int64) :: bytes

      integer(int64) :: n

      n = int(no, int64)*nv
      bytes = 8_int64*(n + merge(n*(n + 1)/2, n*n, symmetric))
   end function rccsd_xfer_amps_bytes

   pure function rccsd_xfer_ladder_bytes(no, nv, symmetric) result(bytes)
      !! Bytes of one partial ladder as sent
      integer, intent(in) :: no, nv
      logical, intent(in) :: symmetric
      integer(int64) :: bytes

      integer(int64) :: n

      n = int(no, int64)*nv
      bytes = 8_int64*merge(n*(n + 1)/2, n*n, symmetric)
   end function rccsd_xfer_ladder_bytes

   ! ---------------------------------------------------------------------
   ! the transfers
   ! ---------------------------------------------------------------------

   subroutine rccsd_xfer_amps_down(comm, lead, xf, amps, symmetric, active)
      !! Collective over the ranks with `active` true. The first rank's amplitudes (on its device) to the device of every
      !! other active rank, in turn; the active others receive. A rank whose share of the ladder is empty is not active:
      !! it has no use for the amplitudes, is not sent them and does not call this. `symmetric` and `active` (one entry
      !! per rank, rank r at r + 1) must be the same on every rank that calls.
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: lead, symmetric
      logical, intent(in) :: active(:)
      type(trc_rccsd_xfer_t), intent(inout) :: xf
      type(trc_rccsd_amps_t), intent(inout) :: amps

      integer :: no, nv, f

      no = xf%no
      nv = xf%nv
      if (lead) then
         !$acc update self(amps%t1)
         if (symmetric) then
            call kernel_pack(no, nv, xf%nhalf, amps%t2, xf%stage)
            if (.not. xf%direct) then
               !$acc update self(xf%stage)
            end if
         else if (.not. xf%direct) then
            !$acc update self(amps%t2)
         end if
         do f = 1, comm%size() - 1
            if (.not. active(f + 1)) cycle
            call send_host(comm, size(amps%t1), amps%t1, f, TAG_T1)
            if (symmetric) then
               call send_array(comm, xf%direct, xf%nhalf, xf%stage, f, TAG_T2)
            else
               call send_array(comm, xf%direct, size(amps%t2), amps%t2, f, TAG_T2)
            end if
         end do
         return
      end if
      call recv_host(comm, size(amps%t1), amps%t1, 0, TAG_T1)
      !$acc update device(amps%t1)
      if (symmetric) then
         call recv_array(comm, xf%direct, xf%nhalf, xf%stage, 0, TAG_T2)
         call kernel_unpack(no, nv, xf%nhalf, xf%stage, amps%t2)
      else
         call recv_array(comm, xf%direct, size(amps%t2), amps%t2, 0, TAG_T2)
      end if
   end subroutine rccsd_xfer_amps_down

   subroutine rccsd_xfer_ladder_up(comm, xf, h, symmetric)
      !! A follower's partial ladder h(i,j,a,b) (on its device) to the first rank. Call `rccsd_xfer_ladder_in` there
      !! for the same follower.
      type(comm_t), intent(in) :: comm
      type(trc_rccsd_xfer_t), intent(inout) :: xf
      real(dp), intent(inout) :: h(:, :, :, :)
      logical, intent(in) :: symmetric

      if (symmetric) then
         call kernel_pack(xf%no, xf%nv, xf%nhalf, h, xf%stage)
         if (.not. xf%direct) then
            !$acc update self(xf%stage)
         end if
         call send_array(comm, xf%direct, xf%nhalf, xf%stage, 0, TAG_LADDER)
      else
         if (.not. xf%direct) then
            !$acc update self(h)
         end if
         call send_array(comm, xf%direct, size(h), h, 0, TAG_LADDER)
      end if
   end subroutine rccsd_xfer_ladder_up

   subroutine rccsd_xfer_ladder_in(comm, xf, h, source, symmetric, whole)
      !! The first rank: follower `source`'s partial ladder added to h on the device. `whole` is a device array of the
      !! size of h, used as the receive area of the general (non-symmetric) case only.
      type(comm_t), intent(in) :: comm
      type(trc_rccsd_xfer_t), intent(inout) :: xf
      real(dp), intent(inout) :: h(:, :, :, :)
      integer, intent(in) :: source
      logical, intent(in) :: symmetric
      real(dp), intent(inout) :: whole(:)

      if (symmetric) then
         call recv_array(comm, xf%direct, xf%nhalf, xf%stage, source, TAG_LADDER)
         call kernel_unpack_add(xf%no, xf%nv, xf%nhalf, xf%stage, h)
      else
         call recv_array(comm, xf%direct, size(h), whole, source, TAG_LADDER)
         call kernel_add(size(h), whole, h)
      end if
   end subroutine rccsd_xfer_ladder_in

   ! ---------------------------------------------------------------------
   ! the messages. Module procedures, no internal ones: this module makes collective-style calls under NVHPC's MPI.
   ! ---------------------------------------------------------------------

   subroutine send_host(comm, n, x, dest, tag)
      !! x(1:n), a host array, to `dest`
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, dest, tag
      real(dp), intent(inout), target :: x(n)

      type(request_t) :: req(1)

      call isend(comm, x, dest, tag, req(1))
      call waitall(req)
   end subroutine send_host

   subroutine recv_host(comm, n, x, source, tag)
      !! n doubles from `source` into the host array x
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, source, tag
      real(dp), intent(inout) :: x(n)

      call recv(comm, x, n, source, tag)
   end subroutine recv_host

   subroutine send_array(comm, direct, n, x, dest, tag)
      !! x(1:n), already in the form the path needs (the device array when `direct`, else its host copy brought up to
      !! date by the caller), to `dest` in MSG_CHUNK pieces. The sends are posted non-blocking and waited for before
      !! returning, as RI-MP2's ring does.
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: direct
      integer, intent(in) :: n, dest, tag
      real(dp), intent(inout), target :: x(n)

      type(request_t), allocatable :: reqs(:)
      integer :: nreq

      nreq = int((int(n, int64) + MSG_CHUNK - 1)/MSG_CHUNK)
      allocate (reqs(max(nreq, 1)))
      if (direct) then
         !$acc host_data use_device(x)
         call post_dev(comm, n, x, dest, tag, reqs)
         !$acc end host_data
      else
         call post_host(comm, n, x, dest, tag, reqs)
      end if
      if (nreq > 0) call waitall(reqs(1:nreq))
   end subroutine send_array

   subroutine recv_array(comm, direct, n, x, source, tag)
      !! n doubles from `source` into the device array x (`direct`) or into its host copy, then to the device
      type(comm_t), intent(in) :: comm
      logical, intent(in) :: direct
      integer, intent(in) :: n, source, tag
      real(dp), intent(inout), target :: x(n)

      if (direct) then
         !$acc host_data use_device(x)
         call take_dev(comm, n, x, source, tag)
         !$acc end host_data
      else
         call take_host(comm, n, x, source, tag)
         !$acc update device(x)
      end if
   end subroutine recv_array

   subroutine post_host(comm, n, x, dest, tag, reqs)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, dest, tag
      real(dp), intent(in), target :: x(n)
      type(request_t), intent(inout) :: reqs(:)

      integer(int64) :: k, m
      integer :: r

      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         call isend(comm, x(k:k + m - 1), dest, tag, reqs(r))
         k = k + m
      end do
   end subroutine post_host

   subroutine take_host(comm, n, x, source, tag)
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, source, tag
      real(dp), intent(inout) :: x(n)

      integer(int64) :: k, m

      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call recv(comm, x(k:k + m - 1), int(m), source, tag)
         k = k + m
      end do
   end subroutine take_host

   subroutine post_dev(comm, n, x, dest, tag, reqs)
      !! Called inside `host_data use_device`: pic-mpi's specific names, as RI-MP2's send_block_dev
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, dest, tag
      real(dp), intent(in), target :: x(n)
      type(request_t), intent(inout) :: reqs(:)

      integer(int64) :: k, m
      integer :: r

      k = 1
      r = 0
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         r = r + 1
         call comm_isend_real_dp_array_n(comm, x(k:k + m - 1), int(m), dest, tag, reqs(r))
         k = k + m
      end do
   end subroutine post_dev

   subroutine take_dev(comm, n, x, source, tag)
      !! Called inside `host_data use_device`, as RI-MP2's recv_block_dev
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n, source, tag
      real(dp), intent(inout), target :: x(n)

      integer(int64) :: k, m

      k = 1
      do while (k <= n)
         m = min(MSG_CHUNK, n - k + 1)
         call comm_recv_real_dp_array_n(comm, x(k:k + m - 1), int(m), source, tag)
         k = k + m
      end do
   end subroutine take_dev

   ! ---------------------------------------------------------------------
   ! the kernels: one flattened index each, items are pure `routine seq` with explicit-shape dummies
   ! ---------------------------------------------------------------------

   subroutine kernel_pack(no, nv, nhalf, x, packed)
      !! packed(p + q (q - 1) / 2) = x(i,j,a,b) for p = (i,a) <= q = (j,b)
      integer, intent(in) :: no, nv, nhalf
      real(dp), intent(in) :: x(no*no*nv*nv)
      real(dp), intent(out) :: packed(nhalf)
      integer :: t, k

      do concurrent(t=1:no*no*nv*nv) local(k)
         k = pair_first(t, no, nv)
         if (k > 0) packed(k) = x(t)
      end do
   end subroutine kernel_pack

   subroutine kernel_unpack(no, nv, nhalf, packed, x)
      !! x(i,j,a,b) = x(j,i,b,a) = the packed value of the pair
      integer, intent(in) :: no, nv, nhalf
      real(dp), intent(in) :: packed(nhalf)
      real(dp), intent(out) :: x(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         x(t) = packed(pair_any(t, no, nv))
      end do
   end subroutine kernel_unpack

   subroutine kernel_unpack_add(no, nv, nhalf, packed, x)
      !! x(i,j,a,b) += the packed value of the pair (i,a),(j,b), for every element
      integer, intent(in) :: no, nv, nhalf
      real(dp), intent(in) :: packed(nhalf)
      real(dp), intent(inout) :: x(no*no*nv*nv)
      integer :: t

      do concurrent(t=1:no*no*nv*nv)
         x(t) = x(t) + packed(pair_any(t, no, nv))
      end do
   end subroutine kernel_unpack_add

   subroutine kernel_add(n, x, h)
      !! h += x
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: h(n)
      integer :: t

      do concurrent(t=1:n)
         h(t) = h(t) + x(t)
      end do
   end subroutine kernel_add

   pure subroutine pair_of(t, no, nv, p, q)
      !! p = (i,a) and q = (j,b) of position t of (no, no, nv, nv), each as i + no (a - 1)
      !$acc routine seq
      integer, intent(in) :: t, no, nv
      integer, intent(out) :: p, q
      integer :: r, i, j, a, b

      r = t - 1
      i = mod(r, no) + 1
      r = r/no
      j = mod(r, no) + 1
      r = r/no
      a = mod(r, nv) + 1
      b = r/nv + 1
      p = i + no*(a - 1)
      q = j + no*(b - 1)
   end subroutine pair_of

   ! The loop index is a VALUE dummy in the item routines below: see the item routines of trc_rccsd_update (an ifx
   ! miscompilation of an index passed by reference).

   pure function pair_first(t, no, nv) result(k)
      !! The packed position of position t when p <= q, else 0
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      integer :: k
      integer :: p, q

      call pair_of(t, no, nv, p, q)
      k = 0
      if (p <= q) k = p + int((int(q, int64)*(q - 1))/2)
   end function pair_first

   pure function pair_any(t, no, nv) result(k)
      !! The packed position of the pair of position t, whichever of p, q is the larger
      !$acc routine seq
      integer, intent(in), value :: t
      integer, intent(in) :: no, nv
      integer :: k
      integer :: p, q

      call pair_of(t, no, nv, p, q)
      k = min(p, q) + int((int(max(p, q), int64)*(max(p, q) - 1))/2)
   end function pair_any

end module trc_rccsd_xfer
