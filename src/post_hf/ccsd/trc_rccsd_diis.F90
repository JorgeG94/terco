!
! Pulay's DIIS on the coupled-cluster amplitudes, and the size of an amplitude step,
! on the device.
!
! The history holds, for up to `size` iterations, the NEW amplitude vector that
! `update` produced and the ERROR vector, the step it took: e = new - old, over
! the concatenated (t1, t2). That is the choice of mqc (mqc_czt_rcc `pack_step`),
! of PySCF (lib.diis on the differences of successive vectors) and of the spin-orbital
! path of mqc: the step vanishes exactly at convergence. The B matrix is the table
! of dot products <e_i, e_j>, one new row and column per push, through the
! linear-algebra handle's dot (cuBLAS on the device) in a fixed order. The
! constrained system
!
!     [ B  -1 ] [c]   [ 0]
!     [-1^T 0 ] [l] = [-1]          sum_i c_i = 1,   x = sum_i c_i vec_i
!
! is solved on the host, assembled oldest entry first and scaled to O(1) as in mqc's
! `diis_coefficients` (an unscaled B of 1e-16 against the -1 border makes the weights
! noise), by Gaussian elimination with partial pivoting and the same pivot floor.
! trc_scf_diis has a `solve_diis` that does the same, but it is private to the
! SCF's own module and not worth publishing for ten lines.
!
! BREAKDOWN. A B matrix that is not finite, or whose elimination meets a pivot below
! the floor (two error vectors that are numerically parallel, which is what a
! converged or stagnated iteration produces), is never used: `next` RESTARTS the
! history (the entry just pushed is dropped too, the next push starts a new one),
! counts it in `restarts`, says so in `status`, and hands back the plain new
! amplitudes. A step that is itself not finite (a NaN or an infinity anywhere in the
! update) is not pushed: `push` sets `finite` to .false. and `status`, and the stored
! vectors, the overlap matrix and the count are unchanged. The bad step was measured in
! an error slot (the oldest entry's when the ring is full), and the next finite push
! overwrites that slot.
!
! THE STEP. Every push, with or without a history, reports in `step_max` the
! largest |e_p| and in `step_norm2` the Euclidean norm ||e||_2. PySCF's
! conv_tol_normt (default 1e-5) is on ||t_new - t_old||_2 of the same concatenated
! vector, before the extrapolation, so `step_norm2` is PySCF's number and `step_max`
! is the stricter-per-element quantity the driver tests (max |e_p| <= ||e||_2 <=
! sqrt(n) max |e_p|). `step_max` is the maximum of fixed chunks of CHUNK elements
! evaluated by one thread each and then of those partials on the host, so it does not
! depend on the thread count.
!
! MEMORY, o = active occupied, v = virtual, n = o v + o^2 v^2, 8 bytes a double,
! host and device each: 2 n size doubles of history (the vectors and the errors) and
! n/CHUNK partials. With size <= 1 there is no history of vectors and one error slot
! (n doubles). At no = 25, nv = 130, size 8: 2 * 8 * 10 565 750 doubles = 1.35 GB.
!
! RESIDENCY. The history exists on the host AND on the device (`enter data` in
! `init`, `exit data` in `release`, on the same `this%name` expression); the device
! copy is the one used, the host mirror is never read. A kernel gets an array only as
! an explicit-shape dummy of a `kernel_*` routine. The extrapolation is two GEMMs
! (one per amplitude block) writing straight into the device copy of the amplitudes
! the caller hands in; the coefficients go up through a small mapped array.
!
! A history vector is addressed with default integers: n above 2^31-1 is refused by
! `init` (status). Single rank only.
module trc_rccsd_diis
   use pic_types, only: dp, int64
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   integer, parameter :: CHUNK = 4096  !! elements per partial maximum of the step
   real(dp), parameter :: PIVOT_FLOOR = 1.0e-14_dp  !! mqc's: below this a pivot is singular

   public :: trc_rccsd_diis_t

   type :: trc_rccsd_diis_t
      !! The amplitude history and the last step
      type(error_t) :: status  !! set by `init`, and by `push` and `next` on a breakdown (the run goes on)
      integer :: no = 0
      integer :: nv = 0
      integer :: n = 0            !! length of one vector: o v + o^2 v^2
      integer :: size = 0         !! history length asked for; below 2 there is no extrapolation
      integer :: nstored = 0      !! entries in the history
      integer :: newest = 0       !! slot of the newest entry
      integer :: restarts = 0     !! breakdowns so far
      integer :: npart = 0        !! partial maxima
      logical :: finite = .true.  !! the last pushed step was finite
      real(dp) :: step_max = 0.0_dp     !! max |e| of the last pushed step
      real(dp) :: step_norm2 = 0.0_dp   !! ||e||_2 of the last pushed step
      real(dp), allocatable :: vec(:, :)      !! (n, size) new amplitudes, only when size >= 2
      real(dp), allocatable :: err(:, :)      !! (n, max(size, 1)) steps new - old
      real(dp), allocatable :: part(:)        !! partial maxima of the last step
      real(dp), allocatable :: coef(:)        !! (max(size, 1)) extrapolation weights, by slot
      real(dp), allocatable :: overlap(:, :)  !! (size, size) host only: <e_i, e_j> by slot
   contains
      procedure :: init => diis_init
      procedure :: release => diis_release
      procedure :: reset => diis_reset
      procedure :: push => diis_push
      procedure :: next => diis_next
   end type trc_rccsd_diis_t

contains

   subroutine diis_init(this, no, nv, size)
      !! An empty history for amplitudes of these dimensions. size < 2 keeps no vectors: only the step is measured.
      class(trc_rccsd_diis_t), intent(inout) :: this
      integer, intent(in) :: no, nv, size

      integer(int64) :: n64
      integer :: nslot

      call this%release()
      call this%status%clear()
      if (no < 1 .or. nv < 1 .or. size < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: need positive orbital counts and a history length >= 0")
         return
      end if
      n64 = int(no, int64)*nv + int(no, int64)*no*nv*nv
      if (n64 > int(huge(0), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: an amplitude vector has more than 2^31-1 elements")
         return
      end if
      this%no = no
      this%nv = nv
      this%n = int(n64)
      this%size = size
      this%npart = (this%n + CHUNK - 1)/CHUNK
      nslot = max(size, 1)
      allocate (this%err(this%n, nslot), this%part(this%npart), this%coef(nslot))
      !$acc enter data create(this%err, this%part, this%coef)
      if (size >= 2) then
         allocate (this%vec(this%n, size), this%overlap(size, size))
         !$acc enter data create(this%vec)
         this%overlap = 0.0_dp
      end if
      call this%reset()
   end subroutine diis_init

   subroutine diis_release(this)
      !! Frees host and device copies; safe on an empty type
      class(trc_rccsd_diis_t), intent(inout) :: this

      if (allocated(this%err)) then
         !$acc exit data delete(this%err)
         deallocate (this%err)
      end if
      if (allocated(this%part)) then
         !$acc exit data delete(this%part)
         deallocate (this%part)
      end if
      if (allocated(this%coef)) then
         !$acc exit data delete(this%coef)
         deallocate (this%coef)
      end if
      if (allocated(this%vec)) then
         !$acc exit data delete(this%vec)
         deallocate (this%vec)
      end if
      if (allocated(this%overlap)) deallocate (this%overlap)
      this%no = 0
      this%nv = 0
      this%n = 0
      this%size = 0
      this%nstored = 0
      this%newest = 0
      this%restarts = 0
      this%npart = 0
   end subroutine diis_release

   subroutine diis_reset(this)
      !! Forgets the history (the restart count and the last step stay)
      class(trc_rccsd_diis_t), intent(inout) :: this

      this%nstored = 0
      this%newest = 0
   end subroutine diis_reset

   function diis_ready(this, old, new) result(ok)
      !! True when the history and both amplitude sets have the same dimensions
      class(trc_rccsd_diis_t), intent(inout) :: this
      type(trc_rccsd_amps_t), intent(in) :: old, new
      logical :: ok

      ok = .false.
      if (.not. allocated(this%err)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: used before init")
         return
      end if
      if (.not. allocated(old%t1) .or. .not. allocated(new%t1) .or. .not. allocated(old%t2) .or. .not. allocated(new%t2)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: the amplitudes are not initialised")
         return
      end if
      if (old%no /= this%no .or. old%nv /= this%nv .or. new%no /= this%no .or. new%nv /= this%nv) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: (no, nv) of the amplitudes differ from the history's")
         return
      end if
      ok = .true.
   end function diis_ready

   subroutine diis_push(this, la, old, new, keep)
      !! Measures the step e = new - old (`step_max`, `step_norm2`, `finite`) and, when `keep` and the history exists,
      !! appends (new, e) to it. A step that is not finite is measured but never appended: the vectors, the
      !! overlaps and the count stay as they were, though the error slot the next push takes (the oldest entry's
      !! when the ring is full) may hold the bad step until that push overwrites it. Once `keep` has been true it
      !! must stay true: while it is false the step is measured in slot 1, which a non-empty history owns.
      class(trc_rccsd_diis_t), intent(inout) :: this
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_amps_t), intent(in) :: old, new
      logical, intent(in) :: keep

      integer :: slot, i, other
      real(dp) :: ee, dot

      call this%status%clear()
      this%finite = .false.
      this%step_max = huge(1.0_dp)
      this%step_norm2 = huge(1.0_dp)
      if (.not. diis_ready(this, old, new)) return
      slot = 1
      if (keep .and. this%size >= 2) slot = modulo(this%newest, this%size) + 1
      call kernel_step(this%no*this%nv, this%n - this%no*this%nv, old%t1, old%t2, new%t1, new%t2, this%err(1, slot))
      call kernel_chunk_max(this%n, CHUNK, this%npart, this%err(1, slot), this%part)
      ee = la%dot(this%n, this%err(1, slot), this%err(1, slot))
      !$acc update self(this%part)
      ! A NaN makes `ee` NaN and an infinity or an overflow makes it huge: neither compares <= huge
      this%finite = ee <= huge(1.0_dp) .and. ee >= 0.0_dp
      if (.not. this%finite) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: the amplitude step is not finite")
         return
      end if
      this%step_norm2 = sqrt(ee)
      this%step_max = maxval(this%part)
      if (.not. keep .or. this%size < 2) return

      call kernel_copy_amps(this%no*this%nv, this%n - this%no*this%nv, new%t1, new%t2, this%vec(1, slot))
      this%newest = slot
      if (this%nstored < this%size) this%nstored = this%nstored + 1
      ! the new row and column: <e_slot, e_j> for every entry held (the diagonal is `ee`)
      do i = 1, this%nstored
         other = slot_of_age(this, i)
         if (other == slot) then
            dot = ee
         else
            dot = la%dot(this%n, this%err(1, slot), this%err(1, other))
         end if
         this%overlap(slot, other) = dot
         this%overlap(other, slot) = dot
      end do
   end subroutine diis_push

   pure function slot_of_age(this, age) result(slot)
      !! Slot of the age-th oldest entry, age = 1 .. nstored (the history is a ring, nothing is shifted)
      class(trc_rccsd_diis_t), intent(in) :: this
      integer, intent(in) :: age
      integer :: slot

      slot = modulo(this%newest - this%nstored + age - 1, this%size) + 1
   end function slot_of_age

   subroutine diis_next(this, la, new, amps, extrapolated)
      !! The next iterate, into `amps`: the DIIS combination of the history when it holds two entries or more and
      !! its system is solvable, else a copy of `new`. A breakdown restarts the history (see the module header).
      !! `amps` may be the `old` of the last push (it is not read here).
      class(trc_rccsd_diis_t), intent(inout) :: this
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_amps_t), intent(in) :: new
      type(trc_rccsd_amps_t), intent(inout) :: amps
      logical, intent(out) :: extrapolated

      real(dp), allocatable :: weights(:)
      integer :: i, n1, n2
      logical :: ok

      extrapolated = .false.
      n1 = this%no*this%nv
      n2 = this%n - n1
      if (this%size >= 2 .and. this%nstored >= 2) then
         allocate (weights(this%nstored))
         call solve_weights(this, weights, ok)
         if (ok) then
            this%coef = 0.0_dp
            do i = 1, this%nstored
               this%coef(slot_of_age(this, i)) = weights(i)
            end do
            !$acc update device(this%coef)
            ! slots 1..nstored are exactly the ones in use (a ring fills them in order before it wraps)
            call la%gemm("N", "N", n1, 1, this%nstored, 1.0_dp, this%vec(1, 1), this%n, this%coef, this%nstored, &
                         0.0_dp, amps%t1, n1)
            call la%gemm("N", "N", n2, 1, this%nstored, 1.0_dp, this%vec(n1 + 1, 1), this%n, this%coef, this%nstored, &
                         0.0_dp, amps%t2, n2)
            extrapolated = .true.
         else
            this%restarts = this%restarts + 1
            call this%reset()
            call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: the DIIS system is singular or not finite; "// &
                                 "history restarted")
         end if
      end if
      if (.not. extrapolated) call kernel_copy_amps_out(n1, n2, new%t1, new%t2, amps%t1, amps%t2)
   end subroutine diis_next

   subroutine solve_weights(this, weights, ok)
      !! The DIIS weights, oldest entry first, from the cached overlaps; `ok` false on a non-finite B or a small pivot
      class(trc_rccsd_diis_t), intent(in) :: this
      real(dp), intent(out) :: weights(:)
      logical, intent(out) :: ok

      real(dp), allocatable :: aug(:, :), coeffs(:)
      real(dp) :: scale
      integer :: m, i, j

      m = this%nstored
      allocate (aug(m + 1, m + 2), coeffs(m + 1))
      aug = 0.0_dp
      do j = 1, m
         do i = 1, m
            aug(i, j) = this%overlap(slot_of_age(this, i), slot_of_age(this, j))
         end do
      end do
      weights = 0.0_dp
      ok = all(abs(aug(1:m, 1:m)) <= huge(1.0_dp))
      if (.not. ok) return
      ! scale the error block to O(1) against the -1 border; the weights are invariant, the multiplier alone rescales
      scale = maxval(abs(aug(1:m, 1:m)))
      if (scale > 0.0_dp) aug(1:m, 1:m) = aug(1:m, 1:m)/scale
      aug(1:m, m + 1) = -1.0_dp
      aug(m + 1, 1:m) = -1.0_dp
      aug(m + 1, m + 2) = -1.0_dp
      call eliminate(m + 1, aug, coeffs, ok)
      if (ok) ok = all(abs(coeffs) <= huge(1.0_dp))
      if (ok) weights = coeffs(1:m)
   end subroutine solve_weights

   subroutine eliminate(m, aug, x, ok)
      !! Gaussian elimination with partial pivoting of the m x m system in aug(:, 1:m) with right-hand side aug(:, m+1)
      integer, intent(in) :: m
      real(dp), intent(inout) :: aug(m, m + 1)
      real(dp), intent(out) :: x(m)
      logical, intent(out) :: ok

      integer :: i, j, p
      real(dp) :: pivot, factor

      ok = .true.
      x = 0.0_dp
      do i = 1, m
         p = i
         do j = i + 1, m
            if (abs(aug(j, i)) > abs(aug(p, i))) p = j
         end do
         if (p /= i) aug([i, p], :) = aug([p, i], :)
         pivot = aug(i, i)
         if (.not. (abs(pivot) >= PIVOT_FLOOR)) then
            ok = .false.
            return
         end if
         do j = i + 1, m
            factor = aug(j, i)/pivot
            aug(j, :) = aug(j, :) - factor*aug(i, :)
         end do
      end do
      do i = m, 1, -1
         x(i) = (aug(i, m + 1) - sum(aug(i, i + 1:m)*x(i + 1:m)))/aug(i, i)
      end do
   end subroutine eliminate

   ! ---------------------------------------------------------------------
   ! the kernels. Each takes its arrays as explicit-shape dummies and runs one
   ! `do concurrent` over ONE flattened index, whose body calls an item routine
   ! of this module. No `x%name` appears in any of them.
   ! ---------------------------------------------------------------------

   subroutine kernel_step(n1, n2, old1, old2, new1, new2, e)
      !! e = (new1 - old1, new2 - old2) as one vector
      integer, intent(in) :: n1, n2
      real(dp), intent(in) :: old1(n1), old2(n2), new1(n1), new2(n2)
      real(dp), intent(out) :: e(n1 + n2)
      integer :: t

      do concurrent(t=1:n1 + n2)
         e(t) = step_item(t, n1, n2, old1, old2, new1, new2)
      end do
   end subroutine kernel_step

   subroutine kernel_copy_amps(n1, n2, a1, a2, v)
      !! v = (a1, a2) as one vector
      integer, intent(in) :: n1, n2
      real(dp), intent(in) :: a1(n1), a2(n2)
      real(dp), intent(out) :: v(n1 + n2)
      integer :: t

      do concurrent(t=1:n1 + n2)
         v(t) = join_item(t, n1, n2, a1, a2)
      end do
   end subroutine kernel_copy_amps

   subroutine kernel_copy_amps_out(n1, n2, a1, a2, b1, b2)
      !! (b1, b2) = (a1, a2)
      integer, intent(in) :: n1, n2
      real(dp), intent(in) :: a1(n1), a2(n2)
      real(dp), intent(out) :: b1(n1), b2(n2)
      integer :: t

      do concurrent(t=1:n1)
         b1(t) = a1(t)
      end do
      do concurrent(t=1:n2)
         b2(t) = a2(t)
      end do
   end subroutine kernel_copy_amps_out

   subroutine kernel_chunk_max(n, chunk_size, npart, e, part)
      !! part(p) = max |e| over the p-th run of `chunk_size` elements
      integer, intent(in) :: n, chunk_size, npart
      real(dp), intent(in) :: e(n)
      real(dp), intent(out) :: part(npart)
      integer :: p

      do concurrent(p=1:npart)
         part(p) = chunk_max_item(p, n, chunk_size, e)
      end do
   end subroutine kernel_chunk_max

   pure function step_item(t, n1, n2, old1, old2, new1, new2) result(v)
      !! new - old at position t of the concatenated vector
      !$acc routine seq
      integer, intent(in) :: t, n1, n2
      real(dp), intent(in) :: old1(n1), old2(n2), new1(n1), new2(n2)
      real(dp) :: v

      if (t <= n1) then
         v = new1(t) - old1(t)
      else
         v = new2(t - n1) - old2(t - n1)
      end if
   end function step_item

   pure function join_item(t, n1, n2, a1, a2) result(v)
      !! element t of the concatenation (a1, a2)
      !$acc routine seq
      integer, intent(in) :: t, n1, n2
      real(dp), intent(in) :: a1(n1), a2(n2)
      real(dp) :: v

      if (t <= n1) then
         v = a1(t)
      else
         v = a2(t - n1)
      end if
   end function join_item

   pure function chunk_max_item(p, n, chunk_size, e) result(v)
      !! max |e| over elements (p-1) chunk_size + 1 .. min(p chunk_size, n), in ascending order
      !$acc routine seq
      integer, intent(in) :: p, n, chunk_size
      real(dp), intent(in) :: e(n)
      real(dp) :: v
      integer :: k

      v = 0.0_dp
      do k = (p - 1)*chunk_size + 1, min(p*chunk_size, n)
         v = max(v, abs(e(k)))
      end do
   end function chunk_max_item

end module trc_rccsd_diis
