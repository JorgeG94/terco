!
! A named set of timed steps, shared by every driver that reports a table
! like RI-MP2's: device-synced wall time per step, an optional flop count
! (GFLOP/s), an optional byte count for a bandwidth-bound step (GB/s
! instead), an NVTX range opened and closed at the same device syncs as the
! timer, and a report routine that logs the table (global_logger, at
! PERFORMANCE) the way RI-MP2 always has. Lifted out of trc_rimp2_driver.F90
! so trc_scf_df.F90 (DF-JK) reports in exactly the same shape rather than a
! second copy of the same bookkeeping.
!
! ONE OBJECT, N STEPS. `init` takes the step names and (for a bandwidth-bound
! step) which ones report GB/s instead of GFLOP/s; `step_on`/`step_off`
! bracket one call to a step, `add_flop`/`add_bytes` record its operation
! count once (added, so a step called many times across an iteration loop
! accumulates both time and operations the same way); `report` prints the
! table, an optional "all dense algebra"-style subtotal row over a caller-
! chosen subset of steps, and a total.
!
module trc_step_table
   use pic_types, only: dp, int64
   use trc_linalg, only: trc_linalg_t
   use pic_flop_recorder, only: flop_recorder_type
   use pic_profiler, only: profiler_start, profiler_stop
   use pic_logger, only: global_logger
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
      !! This step reports GB/s (from `bytes`) rather than GFLOP/s (from `fl`).
      integer :: n = 0
   contains
      procedure :: init => step_table_init
      procedure :: step_on => step_table_on
      procedure :: step_off => step_table_off
      procedure :: add_time => step_table_add_time
      procedure :: add_flop => step_table_add_flop
      procedure :: add_bytes => step_table_add_bytes
      procedure :: report => step_table_report
      procedure :: reset => step_table_reset
   end type step_table_t

contains

   !
   ! `names` fixes the step count and their print order; `bytes_steps`
   ! (optional, by 1-based index into `names`) marks which ones are
   ! bandwidth- rather than flop-bound. Safe to call again (e.g. a second
   ! build on the same object): resets every total to zero.
   !
   subroutine step_table_init(this, prefix, names, bytes_steps)
      class(step_table_t), intent(inout) :: this
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: names(:)
      integer, intent(in), optional :: bytes_steps(:)
      integer :: i
      ! Safe to call more than once on the same object (a df_jk_t released
      ! and set up again on a new geometry, say): deallocate whatever an
      ! earlier `init` left rather than erroring on a re-allocate.
      if (allocated(this%names)) deallocate (this%names, this%t, this%fl, this%bytes, this%is_bytes)
      this%prefix = trim(prefix)
      this%n = size(names)
      allocate (this%names(this%n))
      do i = 1, this%n
         this%names(i) = names(i)
      end do
      allocate (this%t(this%n), this%fl(this%n), this%bytes(this%n), this%is_bytes(this%n))
      this%is_bytes = .false.
      if (present(bytes_steps)) then
         do i = 1, size(bytes_steps)
            this%is_bytes(bytes_steps(i)) = .true.
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
      do i = 1, this%n
         call this%fl(i)%reset()
      end do
   end subroutine step_table_reset

   !
   ! The clock after everything queued on the device has finished -- cuBLAS/
   ! cuSOLVER return to the host before their kernels finish, so a host timer
   ! around them measures the launch, not the work, unless it syncs first.
   !
   function tick(lp) result(t)
      type(trc_linalg_t), intent(in) :: lp
      real(dp) :: t
      integer(int64) :: cc, rate
      call lp%sync()
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function tick

   !
   ! Open step `k`: its start time in `tc`, and its NVTX range
   ! ("<prefix>: <name>", a no-op unless built with TERCO_NVTX).
   !
   subroutine step_table_on(this, lp, k, tc)
      class(step_table_t), intent(in) :: this
      type(trc_linalg_t), intent(in) :: lp
      integer, intent(in) :: k
      real(dp), intent(out) :: tc
      tc = tick(lp)
      call profiler_start(this%prefix//": "//trim(this%names(k)), nvtx_only=.true.)
   end subroutine step_table_on

   !
   ! Close step `k`: the time since `tc` added to its running total, its
   ! NVTX range closed.
   !
   subroutine step_table_off(this, lp, k, tc)
      class(step_table_t), intent(inout) :: this
      type(trc_linalg_t), intent(in) :: lp
      integer, intent(in) :: k
      real(dp), intent(in) :: tc
      this%t(k) = this%t(k) + (tick(lp) - tc)
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
   ! One line per step: device time, operations (or bytes) executed and the
   ! rate; then, when `subtotal_mask` is given, one more line summing the
   ! time and flops of the masked steps under `subtotal_label` (RI-MP2's
   ! "all dense algebra" row); then a total. `rank`/`nranks` are printed
   ! next to the total only -- the caller's own job to decide whether the
   ! per-step numbers are one rank's share or already summed.
   !
   subroutine step_table_report(this, t_total, rank, nranks, subtotal_mask, subtotal_label)
      class(step_table_t), intent(in) :: this
      real(dp), intent(in) :: t_total
      integer, intent(in) :: rank, nranks
      logical, intent(in), optional :: subtotal_mask(:)
      character(len=*), intent(in), optional :: subtotal_label
      integer :: k
      real(dp) :: gf, gb, tg, fg
      character(len=160) :: line
      write (line, "(a,a)") "  "//this%prefix//": step", "                   time (s)        GFLOP     GFLOP/s"
      call global_logger%performance(trim(line))
      do k = 1, this%n
         if (this%is_bytes(k)) then
            gb = real(this%bytes(k), dp)*1.0e-9_dp
            if (this%bytes(k) > 0 .and. this%t(k) > 0.0_dp) then
               write (line, "(a,a24,f12.3,a13,f12.1,a)") "  "//this%prefix//":   ", this%names(k), this%t(k), &
                  "-", gb/this%t(k), " GB/s"
               call global_logger%performance(trim(line))
            else
               write (line, "(a,a24,f12.3,a13,a12)") "  "//this%prefix//":   ", this%names(k), this%t(k), "-", "-"
               call global_logger%performance(trim(line))
            end if
         else
            gf = real(this%fl(k)%get(), dp)*1.0e-9_dp
            if (this%fl(k)%get() > 0 .and. this%t(k) > 0.0_dp) then
               write (line, "(a,a24,f12.3,f13.1,f12.1)") "  "//this%prefix//":   ", this%names(k), this%t(k), &
                  gf, gf/this%t(k)
               call global_logger%performance(trim(line))
            else
               write (line, "(a,a24,f12.3,a13,a12)") "  "//this%prefix//":   ", this%names(k), this%t(k), "-", "-"
               call global_logger%performance(trim(line))
            end if
         end if
      end do
      if (present(subtotal_mask)) then
         tg = 0.0_dp
         fg = 0.0_dp
         do k = 1, this%n
            if (subtotal_mask(k)) then
               tg = tg + this%t(k)
               fg = fg + real(this%fl(k)%get(), dp)*1.0e-9_dp
            end if
         end do
         if (tg > 0.0_dp) then
            write (line, "(a,a24,f12.3,f13.1,f12.1)") "  "//this%prefix//":   ", &
               label(subtotal_label), tg, fg, fg/tg
            call global_logger%performance(trim(line))
         end if
      end if
      write (line, "(a,a24,f12.3,a,i0,a,i0)") "  "//this%prefix//":   ", label("total"), t_total, &
         "   rank ", rank, " of ", nranks
      call global_logger%performance(trim(line))

   contains

      pure function label(t) result(l)
         !! `t` left-justified in the table's column, as a24 would otherwise
         !! right-justify anything shorter.
         character(len=*), intent(in) :: t
         character(len=NAME_LEN) :: l
         l = t
      end function label

   end subroutine step_table_report

end module trc_step_table
