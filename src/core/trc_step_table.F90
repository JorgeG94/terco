!
! A named set of timed steps, shared by every driver that reports a table
! like RI-MP2's: device-synced wall time per step, an optional flop count
! (FLOP/s), an optional byte count for a bandwidth-bound step (B/s
! instead, both SI-prefixed), an NVTX range opened and closed at the same device syncs as the
! timer, and a report routine that prints the table the way RI-MP2 always
! has. Lifted out of trc_rimp2_driver.F90 so trc_scf_df.F90 (DF-JK) reports
! in exactly the same shape rather than a second copy of the same bookkeeping.
!
! ONE OBJECT, N STEPS. `init` takes the step names and (for a bandwidth-bound
! step) which ones report B/s instead of FLOP/s; `step_on`/`step_off`
! bracket one call to a step, `add_flop`/`add_bytes` record its operation
! count once (added, so a step called many times across an iteration loop
! accumulates both time and operations the same way); `report` prints the
! table, an optional "all dense algebra"-style subtotal row over a caller-
! chosen subset of steps, and a total.
!
! OVER MPI RANKS each rank fills its own table; `reduce` (collective) turns
! them into the job's: a step's time is the slowest rank's, since the job
! waits for it, and its operations are summed over the ranks, except for a
! step every rank runs in full (`replicated_steps` in `init`), whose work
! counts once. A per-rank table instead shows one rank's share, and its
! rate stays flat however many ranks run.
!
module trc_step_table
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use trc_linalg, only: trc_linalg_t
   use pic_flop_recorder, only: flop_recorder_type
   use pic_profiler, only: profiler_start, profiler_stop
   use trc_si, only: si_string
   use pic_strings, only: to_string
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM, MPI_MAX
   implicit none
   private

   public :: step_table_t

   integer, parameter :: NAME_LEN = 24

   type :: step_table_t
      character(len=:), allocatable :: prefix
      !! The NVTX range / print-line prefix, e.g. "rimp2" or "dfjk".
      character(len=NAME_LEN), allocatable :: names(:)
      real(dp), allocatable :: t(:)
      type(flop_recorder_type), allocatable :: fl(:)
      integer(int64), allocatable :: bytes(:)
      logical, allocatable :: is_bytes(:)
      !! This step reports B/s (from `bytes`) rather than FLOP/s (from `fl`).
      logical, allocatable :: replicated(:)
      !! Every rank runs this step in full: `reduce` counts its work once.
      integer :: n = 0
      integer :: nranks_reduced = 0
      !! The rank count `reduce` summed over; 0 while the table is one rank's.
   contains
      procedure :: init => step_table_init
      procedure :: step_on => step_table_on
      procedure :: step_off => step_table_off
      procedure :: add_time => step_table_add_time
      procedure :: add_flop => step_table_add_flop
      procedure :: add_bytes => step_table_add_bytes
      procedure :: reduce => step_table_reduce
      procedure :: report => step_table_report
      procedure :: reset => step_table_reset
   end type step_table_t

contains

   !
   ! `names` fixes the step count and their print order; `bytes_steps`
   ! (optional, by 1-based index into `names`) marks which ones are
   ! bandwidth- rather than flop-bound, and `replicated_steps` which ones
   ! every rank runs in full. Safe to call again (e.g. a second build on the
   ! same object): resets every total to zero.
   !
   subroutine step_table_init(this, prefix, names, bytes_steps, replicated_steps)
      class(step_table_t), intent(inout) :: this
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: names(:)
      integer, intent(in), optional :: bytes_steps(:)
      integer, intent(in), optional :: replicated_steps(:)
      integer :: i
      ! Safe to call more than once on the same object (a df_jk_t released
      ! and set up again on a new geometry, say): deallocate whatever an
      ! earlier `init` left rather than erroring on a re-allocate.
      if (allocated(this%names)) deallocate (this%names, this%t, this%fl, this%bytes, this%is_bytes, this%replicated)
      this%prefix = trim(prefix)
      this%n = size(names)
      allocate (this%names(this%n))
      do i = 1, this%n
         this%names(i) = names(i)
      end do
      allocate (this%t(this%n), this%fl(this%n), this%bytes(this%n), this%is_bytes(this%n), &
                this%replicated(this%n))
      this%is_bytes = .false.
      if (present(bytes_steps)) then
         do i = 1, size(bytes_steps)
            this%is_bytes(bytes_steps(i)) = .true.
         end do
      end if
      this%replicated = .false.
      if (present(replicated_steps)) then
         do i = 1, size(replicated_steps)
            this%replicated(replicated_steps(i)) = .true.
         end do
      end if
      call this%reset()
   end subroutine step_table_init

   subroutine step_table_reset(this)
      !! Every step's time, flop and byte total back to zero -- a fresh
      !! accumulation window (a new SCF run on an object kept across runs).
      class(step_table_t), intent(inout) :: this
      integer :: i
      this%t = 0.0_dp
      this%bytes = 0_int64
      this%nranks_reduced = 0
      do i = 1, this%n
         call this%fl(i)%reset()
      end do
   end subroutine step_table_reset

   !
   ! Open step `k`: `tc` starts once everything queued on the device has
   ! finished -- cuBLAS/cuSOLVER return to the host before their kernels
   ! finish, so a host timer around them measures the launch, not the work,
   ! unless it syncs first -- and its NVTX range ("<prefix>: <name>", a no-op
   ! unless built with TERCO_NVTX) opens.
   !
   subroutine step_table_on(this, lp, k, tc)
      class(step_table_t), intent(in) :: this
      type(trc_linalg_t), intent(in) :: lp
      integer, intent(in) :: k
      type(timer_type), intent(out) :: tc
      call lp%sync()
      call tc%start()
      call profiler_start(this%prefix//": "//trim(this%names(k)), nvtx_only=.true.)
   end subroutine step_table_on

   !
   ! Close step `k`: after the same device sync, the time `tc` has run is
   ! added to the step's total, and its NVTX range closes.
   !
   subroutine step_table_off(this, lp, k, tc)
      class(step_table_t), intent(inout) :: this
      type(trc_linalg_t), intent(in) :: lp
      integer, intent(in) :: k
      type(timer_type), intent(in) :: tc
      call lp%sync()
      this%t(k) = this%t(k) + tc%get_elapsed_time()
      call profiler_stop(this%prefix//": "//trim(this%names(k)))
   end subroutine step_table_off

   !
   ! Add a time already measured elsewhere -- a caller that cannot bracket
   ! its own kernel with `step_on`/`step_off` (it only gets a finished
   ! elapsed time back, the way `trc_xc_rks`'s `t_points`/`t_pairs` or
   ! `trc_scf_context`'s grid-build timer already do) still gets its total
   ! into this table's `report`. No NVTX range: one cannot be drawn around
   ! work that has already finished by the time this is called.
   !
   subroutine step_table_add_time(this, k, t)
      class(step_table_t), intent(inout) :: this
      integer, intent(in) :: k
      real(dp), intent(in) :: t
      this%t(k) = this%t(k) + t
   end subroutine step_table_add_time

   subroutine step_table_add_flop(this, k, n)
      class(step_table_t), intent(inout) :: this
      integer, intent(in) :: k
      integer(int64), intent(in) :: n
      call this%fl(k)%add(n)
   end subroutine step_table_add_flop

   subroutine step_table_add_bytes(this, k, n)
      class(step_table_t), intent(inout) :: this
      integer, intent(in) :: k
      integer(int64), intent(in) :: n
      this%bytes(k) = this%bytes(k) + n
   end subroutine step_table_add_bytes

   !
   ! Collective over `comm`: every rank's table becomes the job's. Each
   ! step's time becomes the slowest rank's; its operations and bytes the
   ! sum over the ranks, or this rank's own for a replicated step, which
   ! every rank ran in full. Call it once, after the last step and before
   ! `report`; `report` then labels the total with the rank count.
   ! Operations travel as doubles (pic-mpi has no int64 reduction): exact
   ! below 2^53, about 9 PFLOP per step.
   !
   subroutine step_table_reduce(this, comm)
      class(step_table_t), intent(inout) :: this
      type(comm_t), intent(in) :: comm
      real(dp), allocatable :: w(:)
      integer :: k
      call allreduce(comm, this%t, op=MPI_MAX)
      allocate (w(2*this%n))
      do k = 1, this%n
         if (this%replicated(k) .and. comm%rank() /= 0) then
            w(k) = 0.0_dp
            w(this%n + k) = 0.0_dp
         else
            w(k) = real(this%fl(k)%get(), dp)
            w(this%n + k) = real(this%bytes(k), dp)
         end if
      end do
      call allreduce(comm, w, op=MPI_SUM)
      do k = 1, this%n
         call this%fl(k)%reset()
         call this%fl(k)%add(nint(w(k), int64))
         this%bytes(k) = nint(w(this%n + k), int64)
      end do
      this%nranks_reduced = comm%size()
   end subroutine step_table_reduce

   !
   ! One line per step: device time, operations (or bytes) executed and the
   ! rate, each count printed with the SI prefix that keeps it under 1000
   ! (trc_si); then, when `subtotal_mask` is given, one more line summing
   ! the time and flops of the masked steps under `subtotal_label` (RI-MP2's
   ! "all dense algebra" row); then a total, which says whose numbers these
   ! are: "all N ranks" after `reduce`, else "rank r of N", one rank's share.
   !
   subroutine step_table_report(this, t_total, rank, nranks, subtotal_mask, subtotal_label)
      class(step_table_t), intent(in) :: this
      real(dp), intent(in) :: t_total
      integer, intent(in) :: rank, nranks
      logical, intent(in), optional :: subtotal_mask(:)
      character(len=*), intent(in), optional :: subtotal_label
      integer, parameter :: TIME_W = 12, WORK_W = 15, RATE_W = 16
      integer :: k
      real(dp) :: tg, fg, fa
      character(len=:), allocatable :: whose
      if (this%nranks_reduced > 0) then
         whose = "   all "//to_string(this%nranks_reduced)//" ranks"
      else
         whose = "   rank "//to_string(rank)//" of "//to_string(nranks)
      end if
      call global_logger%performance("  "//this%prefix//": step"//repeat(" ", 3 + NAME_LEN - 4) &
                                     //right("time (s)", TIME_W)//right("work", WORK_W)//right("rate", RATE_W))
      do k = 1, this%n
         if (this%is_bytes(k)) then
            if (this%bytes(k) > 0 .and. this%t(k) > 0.0_dp) then
               call row(this%names(k), this%t(k), right("-", WORK_W), &
                        si_string(real(this%bytes(k), dp)/this%t(k), "B/s", RATE_W))
            else
               call row(this%names(k), this%t(k), right("-", WORK_W), right("-", RATE_W))
            end if
         else
            if (this%fl(k)%get() > 0 .and. this%t(k) > 0.0_dp) then
               call row(this%names(k), this%t(k), si_string(real(this%fl(k)%get(), dp), "FLOP", WORK_W), &
                        si_string(real(this%fl(k)%get(), dp)/this%t(k), "FLOP/s", RATE_W))
            else
               call row(this%names(k), this%t(k), right("-", WORK_W), right("-", RATE_W))
            end if
         end if
      end do
      if (present(subtotal_mask)) then
         tg = 0.0_dp
         fg = 0.0_dp
         do k = 1, this%n
            if (subtotal_mask(k)) then
               tg = tg + this%t(k)
               fg = fg + real(this%fl(k)%get(), dp)
            end if
         end do
         if (tg > 0.0_dp) call row(label(subtotal_label), tg, si_string(fg, "FLOP", WORK_W), &
                                   si_string(fg/tg, "FLOP/s", RATE_W))
      end if
      ! The whole table's work over its whole wall time (setup, transfers and
      ! untimed host work included): the rate the run actually achieved.
      fa = 0.0_dp
      do k = 1, this%n
         if (.not. this%is_bytes(k)) fa = fa + real(this%fl(k)%get(), dp)
      end do
      if (fa > 0.0_dp .and. t_total > 0.0_dp) then
         call row(label("total"), t_total, si_string(fa, "FLOP", WORK_W), si_string(fa/t_total, "FLOP/s", RATE_W), &
                  whose)
      else
         call row(label("total"), t_total, right("-", WORK_W), right("-", RATE_W), &
                  whose)
      end if

   contains

      subroutine row(name, t, work, rate, tail)
         character(len=*), intent(in) :: name, work, rate
         real(dp), intent(in) :: t
         character(len=*), intent(in), optional :: tail
         character(len=TIME_W) :: tbuf
         write (tbuf, "(f12.3)") t
         if (present(tail)) then
            call global_logger%performance("  "//this%prefix//":   "//label(name)//tbuf//work//rate//tail)
         else
            call global_logger%performance("  "//this%prefix//":   "//label(name)//tbuf//work//rate)
         end if
      end subroutine row

      pure function right(t, w) result(r)
         !! `t` right-justified in `w` characters.
         character(len=*), intent(in) :: t
         integer, intent(in) :: w
         character(len=:), allocatable :: r
         r = repeat(" ", max(0, w - len(t)))//t
      end function right

      pure function label(t) result(l)
         !! `t` left-justified in the table's column, as a24 would otherwise
         !! right-justify anything shorter.
         character(len=*), intent(in) :: t
         character(len=NAME_LEN) :: l
         l = t
      end function label

   end subroutine step_table_report

end module trc_step_table
