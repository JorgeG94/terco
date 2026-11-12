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
! PACKED HISTORY (`init(..., packed=.true.)`). The solver's amplitudes have t2(i,j,a,b) = t2(j,i,b,a): the MP2
! start has it, the update keeps it to round-off and a DIIS combination is linear. With p = (i,a) and q = (j,b)
! (n_ov = o v values each) t2 is then a symmetric matrix T(p,q), and the history keeps only p >= q, column by
! column: a stored vector is the t1 block, then n_ov (n_ov + 1)/2 elements of t2, np = o v + n_ov (n_ov + 1)/2
! in all instead of n = o v + n_ov^2. The amplitude vector holds (T(p,q) + T(q,p))/2 off the diagonal and
! T(p,p) on it; the extrapolation combines the packed vectors (a GEMM for t1, one kernel for t2) and writes
! both (p,q) and (q,p) of the result, which is therefore exactly pair-symmetric (the iterates move at
! round-off only). The error vector holds the symmetrised step with its off-diagonal elements multiplied by
! sqrt(2), so that the plain dot of two packed errors is the full dot of the two symmetric steps
! (sum_pq e f = sum_p e_pp f_pp + 2 sum_{p>q} e_pq f_pq): B, hence the weights, are the full history's to
! round-off, through the same dot on vectors half as long. Packing is the caller's decision and is valid only
! for pair-symmetric amplitudes: the driver asks for it exactly when its ladder is told the amplitudes have
! the symmetry (`tau_symmetric`), so a run started from amplitudes of the caller's that lack it keeps the
! full history. The default is the full history, unchanged.
!
! THE STEP. Every push, with or without a history, reports in `step_max` the
! largest |e_p| and in `step_norm2` the Euclidean norm ||e||_2, both over the FULL vector also when the
! history is packed (the scaled packed dot is ||e||_2 of a symmetric step; the maximum is taken from the old
! and new amplitudes themselves, not from the averaged and scaled half). PySCF's
! conv_tol_normt (default 1e-5) is on the norm of t1 and the PACKED lower triangle of t2 (each off-diagonal pair
! once) of t_new - t_old, before the extrapolation, so `step_norm2`, which counts every pair of the full t2, is not
! PySCF's number: normt <= step_norm2 <= sqrt(2) normt (trc_rccsd_driver, CONVERGENCE). `step_max` is the
! stricter-per-element quantity the driver tests (max |e_p| <= ||e||_2 <= sqrt(n) max |e_p|).
! `step_max` is the maximum of npart = ceil(n_full/CHUNK) fixed partials, partial p over
! the elements p, p + npart, p + 2 npart, ..., one thread each, then of the partials on the host: a maximum is
! exact in any order, so it does not depend on the thread count. Interleaved, neighbouring threads read
! neighbouring elements: with a contiguous run of 4096 elements a thread, the partials of a push took 0.83 ms
! at o 20, v 265 (GH200), interleaved over 512 0.14 ms.
!
! MEMORY, o = active occupied, v = virtual, 8 bytes a double, host and device each: 2 np size doubles of
! history packed (2 n size full; the vectors and the errors) and npart partials. With size <= 1 there is no
! history of vectors and one error slot. At o 40, v 530, size 8 (10 waters, cc-pVTZ): packed 28.8 GB, full
! 57.5 GB.
!
! RESIDENCY. The history exists on the host AND on the device (`enter data` in
! `init`, `exit data` in `release`, on the same `this%name` expression); the device
! copy is the one used, the host mirror is never read.
!
! HOST HISTORY (`init(..., host=.true.)`). `vec` and `err` live on the host, and a column (one slot) is on the
! device only while it is used (`enter data` on the section), instead of 2 n size doubles resident: at o 40,
! v 530, size 8 packed, at most three vectors (5.4 GB) on the device instead of 28.8 GB. A push maps the next
! slot's two columns, computes the step and the packed new amplitudes into them with the same kernels, copies them
! down and unmaps them; the new row of B is computed on the device, each older error column copied up in turn and
! dotted with cuBLAS. The extrapolation is one host GEMV over the contiguous columns into `comb`, copied up and
! unpacked by the same kernel with one vector of weight 1 (exact). No host-side copy of a vector and no host dot:
! on GH200 those run on one thread (14-20 GB/s) while the link moves 125-300 GB/s, and a version that copied into
! the ring and dotted on the host cost 0.8 s an iteration at (H2O)9
! cc-pVTZ against 0.02 for the device history. `step_max` and `step_norm2` are computed on the device as before,
! so they are bitwise those of the device history; B and the weights differ from it at round-off (cuBLAS dots on
! the same vectors, but an extrapolation summed by the host BLAS). The host buffers are zeroed in `init`: a copy
! into host memory never written runs at a fraction of the link's speed on GH200. A kernel gets an array only as
! an explicit-shape dummy of a `kernel_*` routine. The extrapolation writes straight into
! the device copy of the amplitudes the caller hands in; the coefficients go up through a
! small mapped array.
!
! A stored vector is addressed with default integers: a stored length (np packed, n full) above 2^31-1 is
! refused by `init` (status); the full amplitudes are indexed with 64-bit integers.
! No MPI call: over several ranks the first rank alone holds the history (trc_rccsd_driver, RANKS).
#include "inc/trc_dc_index.inc"

module trc_rccsd_diis
   use pic_types, only: dp, int64
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use pic_strings, only: to_string
   use pic_blas_interfaces, only: pic_gemv
   implicit none
   private

   integer, parameter :: CHUNK = 512  !! elements per partial maximum of the step (every npart-th element)
   real(dp), parameter :: PIVOT_FLOOR = 1.0e-14_dp  !! mqc's: below this a pivot is singular
   real(dp), parameter :: SQRT2 = 1.4142135623730951_dp  !! off-diagonal weight of a packed error element

   public :: trc_rccsd_diis_t
   public :: rccsd_diis_history_bytes, rccsd_diis_free_device_bytes


   type :: trc_rccsd_diis_t
      !! The amplitude history and the last step
      type(error_t) :: status  !! set by `init`, and by `push` and `next` on a breakdown (the run goes on)
      integer :: no = 0
      integer :: nv = 0
      integer :: n = 0            !! length of one STORED vector: o v + n_ov (n_ov + 1)/2 packed, o v + o^2 v^2 full
      integer(int64) :: nfull = 0  !! length of the full vector, o v + o^2 v^2
      logical :: packed = .false.  !! t2 kept as the unique half (pair-symmetric amplitudes only, see the header)
      logical :: host = .false.   !! the history lives on the host, a slot on the device only while used (header)
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
      real(dp), allocatable :: comb(:)        !! (n) host history only: the extrapolation, host GEMV then up
      real(dp), allocatable :: unit(:)        !! (1) host history only: the weight 1 that writes `comb` out
   contains
      procedure :: init => diis_init
      procedure :: release => diis_release
      procedure :: reset => diis_reset
      procedure :: push => diis_push
      procedure :: next => diis_next
   end type trc_rccsd_diis_t

contains

   pure function rccsd_diis_history_bytes(no, nv, size, packed) result(bytes)
      !! Device bytes of the history kept on the device (`init` without `host`): 2 n size doubles, n the stored
      !! length (header, MEMORY); 0 below size 2, which keeps no vectors
      integer, intent(in) :: no, nv, size
      logical, intent(in) :: packed
      integer(int64) :: bytes
      integer(int64) :: nov, nstore

      bytes = 0
      if (size < 2) return
      nov = int(no, int64)*nv
      nstore = nov + nov*nov
      if (packed) nstore = nov + nov*(nov + 1)/2
      bytes = 16_int64*nstore*size
   end function rccsd_diis_history_bytes

   function rccsd_diis_free_device_bytes() result(bytes)
      !! The free device memory now; -1 on the host build or when the query gives nothing (no device to fit)
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, acc_property_free_memory
#endif
      integer(int64) :: bytes

      bytes = -1
#ifdef TRC_DEVICE
      bytes = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                   acc_property_free_memory), int64)
      if (bytes <= 0) bytes = -1
#endif
   end function rccsd_diis_free_device_bytes

   subroutine diis_init(this, no, nv, size, packed, host)
      !! An empty history for amplitudes of these dimensions. size < 2 keeps no vectors: only the step is measured.
      !! `packed` (default false) stores t2 as its unique half: only for amplitudes with t2(i,j,a,b) = t2(j,i,b,a).
      !! `host` (default false) keeps the history on the host (header, HOST HISTORY)
      class(trc_rccsd_diis_t), intent(inout) :: this
      integer, intent(in) :: no, nv, size
      logical, intent(in), optional :: packed
      logical, intent(in), optional :: host

      integer(int64) :: nstore, nov
      integer :: nslot

      call this%release()
      call this%status%clear()
      if (no < 1 .or. nv < 1 .or. size < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: need positive orbital counts and a history length >= 0")
         return
      end if
      this%packed = .false.
      if (present(packed)) this%packed = packed
      nov = int(no, int64)*nv
      this%nfull = nov + nov*nov
      nstore = this%nfull
      if (this%packed) nstore = nov + nov*(nov + 1)/2
      if (nstore > int(huge(0), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: a stored amplitude vector would have "// &
                              to_string(nstore)//" elements, above the limit of 2^31-1 (o v + o v (o v + 1)/2 packed, "// &
                              "o v + o^2 v^2 full)")
         return
      end if
      this%no = no
      this%nv = nv
      this%n = int(nstore)
      this%size = size
      this%npart = int((this%nfull + CHUNK - 1)/CHUNK)
      nslot = max(size, 1)
      this%host = .false.
      if (present(host)) this%host = host .and. size >= 2
      if (this%host) then
         ! the ring on the host, zeroed now (first touch); only the partials and the weights stay on the device
         allocate (this%err(this%n, nslot), this%vec(this%n, size), this%overlap(size, size))
         allocate (this%comb(this%n), this%part(this%npart), this%coef(nslot), this%unit(1))
         this%err = 0.0_dp
         this%vec = 0.0_dp
         this%comb = 0.0_dp
         this%overlap = 0.0_dp
         this%unit = 1.0_dp
         !$acc enter data create(this%part, this%coef)
         !$acc enter data copyin(this%unit)
         call this%reset()
         return
      end if
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

      if (this%host) then
         ! nothing of the ring is mapped between calls
         if (allocated(this%err)) deallocate (this%err)
         if (allocated(this%vec)) deallocate (this%vec)
         if (allocated(this%comb)) deallocate (this%comb)
         if (allocated(this%unit)) then
            !$acc exit data delete(this%unit)
            deallocate (this%unit)
         end if
      end if
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
      this%nfull = 0
      this%packed = .false.
      this%host = .false.
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

      integer :: slot, i, other, n1
      integer(int64) :: n2
      real(dp) :: ee, dot

      call this%status%clear()
      this%finite = .false.
      this%step_max = huge(1.0_dp)
      this%step_norm2 = huge(1.0_dp)
      if (.not. diis_ready(this, old, new)) return
      if (this%host) then
         call push_host(this, la, old, new, keep)
         return
      end if
      slot = 1
      if (keep .and. this%size >= 2) slot = modulo(this%newest, this%size) + 1
      n1 = this%no*this%nv
      n2 = this%nfull - n1
      if (this%packed) then
         call kernel_step_packed(this%no, this%nv, n2, this%n, old%t1, old%t2, new%t1, new%t2, this%err(1, slot))
         ! the maximum of the full step, not of the averaged and sqrt(2)-scaled half
         call kernel_chunk_max_diff(n1, n2, this%npart, old%t1, old%t2, new%t1, new%t2, this%part)
      else
         call kernel_step(n1, this%n - n1, old%t1, old%t2, new%t1, new%t2, this%err(1, slot))
         call kernel_chunk_max(this%n, this%npart, this%err(1, slot), this%part)
      end if
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

      if (this%packed) then
         call kernel_pack_amps(this%no, this%nv, n2, this%n, new%t1, new%t2, this%vec(1, slot))
      else
         call kernel_copy_amps(n1, this%n - n1, new%t1, new%t2, this%vec(1, slot))
      end if
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

   subroutine push_host(this, la, old, new, keep)
      !! `push` with the history on the host: the next slot's buffers are mapped, the step and the packed new
      !! amplitudes computed into them by the same kernels and measured there, copied down and unmapped; the new
      !! row of B is computed on the device, each older error mapped in turn
      class(trc_rccsd_diis_t), intent(inout) :: this
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_amps_t), intent(in) :: old, new
      logical, intent(in) :: keep

      integer :: slot, i, other, n1
      integer(int64) :: n2
      real(dp) :: ee, dot

      n1 = this%no*this%nv
      n2 = this%nfull - n1
      ! the slot the entry would take; before the first kept push it is free, after it `keep` stays true
      slot = modulo(this%newest, this%size) + 1
      !$acc enter data create(this%err(:, slot), this%vec(:, slot))
      if (this%packed) then
         call kernel_step_packed(this%no, this%nv, n2, this%n, old%t1, old%t2, new%t1, new%t2, this%err(1, slot))
         call kernel_chunk_max_diff(n1, n2, this%npart, old%t1, old%t2, new%t1, new%t2, this%part)
      else
         call kernel_step(n1, this%n - n1, old%t1, old%t2, new%t1, new%t2, this%err(1, slot))
         call kernel_chunk_max(this%n, this%npart, this%err(1, slot), this%part)
      end if
      ee = la%dot(this%n, this%err(1, slot), this%err(1, slot))
      !$acc update self(this%part)
      this%finite = ee <= huge(1.0_dp) .and. ee >= 0.0_dp
      if (.not. this%finite .or. .not. keep) then
         ! the host buffers were not written: a refused or unkept step leaves the ring as it was
         !$acc exit data delete(this%err(:, slot), this%vec(:, slot))
         if (.not. this%finite) call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: the amplitude step is not finite")
         if (this%finite) then
            this%step_norm2 = sqrt(ee)
            this%step_max = maxval(this%part)
         end if
         return
      end if
      this%step_norm2 = sqrt(ee)
      this%step_max = maxval(this%part)
      if (this%packed) then
         call kernel_pack_amps(this%no, this%nv, n2, this%n, new%t1, new%t2, this%vec(1, slot))
      else
         call kernel_copy_amps(n1, this%n - n1, new%t1, new%t2, this%vec(1, slot))
      end if
      !$acc exit data copyout(this%vec(:, slot))
      !$acc update self(this%err(:, slot))
      this%newest = slot
      if (this%nstored < this%size) this%nstored = this%nstored + 1
      ! the new row on the device: each older error up in turn, dotted with the new one where it already is
      do i = 1, this%nstored
         other = slot_of_age(this, i)
         if (other == slot) then
            dot = ee
         else
            !$acc enter data copyin(this%err(:, other))
            dot = la%dot(this%n, this%err(1, slot), this%err(1, other))
            !$acc exit data delete(this%err(:, other))
         end if
         this%overlap(slot, other) = dot
         this%overlap(other, slot) = dot
      end do
      !$acc exit data delete(this%err(:, slot))
   end subroutine push_host

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
            if (this%host) then
               call next_host(this, la, amps)
               extrapolated = .true.
               return
            end if
            !$acc update device(this%coef)
            ! slots 1..nstored are exactly the ones in use (a ring fills them in order before it wraps)
            call la%gemm("N", "N", n1, 1, this%nstored, 1.0_dp, this%vec(1, 1), this%n, this%coef, this%nstored, &
                         0.0_dp, amps%t1, n1)
            if (this%packed) then
               ! the half, combined and written to (i,j,a,b) and (j,i,b,a) in one pass: no packed scratch
               call kernel_extrapolate_unpack(this%no, this%nv, this%nfull - n1, this%n, this%nstored, this%vec, &
                                              this%coef, amps%t2)
            else
               call la%gemm("N", "N", n2, 1, this%nstored, 1.0_dp, this%vec(n1 + 1, 1), this%n, this%coef, &
                            this%nstored, 0.0_dp, amps%t2, n2)
            end if
            extrapolated = .true.
         else
            this%restarts = this%restarts + 1
            call this%reset()
            call this%status%set(ERROR_VALIDATION, "trc_rccsd_diis: the DIIS system is singular or not finite; "// &
                                 "history restarted")
         end if
      end if
      if (.not. extrapolated) call kernel_copy_amps_out(n1, this%nfull - n1, new%t1, new%t2, amps%t1, amps%t2)
   end subroutine diis_next

   subroutine next_host(this, la, amps)
      !! The extrapolation with the history on the host: comb = vec(:, 1:nstored) coef, one host GEMV over the
      !! contiguous columns (threaded where the host BLAS is), copied up and written to the amplitudes as one vector
      !! of weight 1, so the device side is a copy, exactly
      class(trc_rccsd_diis_t), intent(inout) :: this
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_amps_t), intent(inout) :: amps

      integer :: n1, n2

      n1 = this%no*this%nv
      n2 = this%n - n1
      ! slots 1..nstored are exactly the ones in use (a ring fills them in order before it wraps)
      call pic_gemv(this%vec(:, 1:this%nstored), this%coef(1:this%nstored), this%comb, alpha=1.0_dp, beta=0.0_dp)
      ! out through the unit weight: `coef` keeps the DIIS weights, by slot, as with the device history
      !$acc enter data copyin(this%comb)
      call la%gemm("N", "N", n1, 1, 1, 1.0_dp, this%comb, this%n, this%unit, 1, 0.0_dp, amps%t1, n1)
      if (this%packed) then
         call kernel_extrapolate_unpack(this%no, this%nv, this%nfull - n1, this%n, 1, this%comb, this%unit, amps%t2)
      else
         call la%gemm("N", "N", n2, 1, 1, 1.0_dp, this%comb(n1 + 1), this%n, this%unit, 1, 0.0_dp, amps%t2, n2)
      end if
      !$acc exit data delete(this%comb)
   end subroutine next_host

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
         e(t) = step_item(TRC_DC_INDEX(t), n1, n2, old1, old2, new1, new2)
      end do
   end subroutine kernel_step

   subroutine kernel_copy_amps(n1, n2, a1, a2, v)
      !! v = (a1, a2) as one vector
      integer, intent(in) :: n1, n2
      real(dp), intent(in) :: a1(n1), a2(n2)
      real(dp), intent(out) :: v(n1 + n2)
      integer :: t

      do concurrent(t=1:n1 + n2)
         v(t) = join_item(TRC_DC_INDEX(t), n1, n2, a1, a2)
      end do
   end subroutine kernel_copy_amps

   subroutine kernel_copy_amps_out(n1, n2, a1, a2, b1, b2)
      !! (b1, b2) = (a1, a2)
      integer, intent(in) :: n1
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: a1(n1), a2(n2)
      real(dp), intent(out) :: b1(n1), b2(n2)
      integer :: t
      integer(int64) :: s

      do concurrent(t=1:n1)
         b1(t) = a1(t)
      end do
      do concurrent(s=1:n2)
         b2(s) = a2(s)
      end do
   end subroutine kernel_copy_amps_out

   subroutine kernel_step_packed(no, nv, n2, nstore, old1, old2, new1, new2, e)
      !! e = (new1 - old1, the packed half of the symmetrised t2 step scaled by sqrt(2) off the diagonal)
      integer, intent(in) :: no, nv, nstore
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: old1(no*nv), old2(n2), new1(no*nv), new2(n2)
      real(dp), intent(out) :: e(nstore)
      integer :: t
      integer(int64) :: s, k, r

      do concurrent(t=1:no*nv)
         e(t) = new1(t) - old1(t)
      end do
      ! one thread per element of the FULL t2, in its storage order, the ones of the upper half idle: the partner
      ! (j,i,b,a) of a run over (i,j) lies in the one contiguous no x no block (b,a), so both reads stay local
      do concurrent(s=1:n2) local(k, r)
         call pair_indices(TRC_DC_INDEX(s), no, nv, k, r)
         if (k > 0) e(no*nv + k) = packed_step_item(TRC_DC_INDEX(s), TRC_DC_INDEX(r), n2, old2, new2)
      end do
   end subroutine kernel_step_packed

   subroutine kernel_pack_amps(no, nv, n2, nstore, a1, a2, v)
      !! v = (a1, the packed half of (a2(i,j,a,b) + a2(j,i,b,a))/2)
      integer, intent(in) :: no, nv, nstore
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: a1(no*nv), a2(n2)
      real(dp), intent(out) :: v(nstore)
      integer :: t
      integer(int64) :: s, k, r

      do concurrent(t=1:no*nv)
         v(t) = a1(t)
      end do
      do concurrent(s=1:n2) local(k, r)
         call pair_indices(TRC_DC_INDEX(s), no, nv, k, r)
         if (k > 0) v(no*nv + k) = 0.5_dp*(a2(s) + a2(r))
      end do
   end subroutine kernel_pack_amps

   subroutine kernel_extrapolate_unpack(no, nv, n2, nstore, nuse, vec, coef, t2)
      !! t2(i,j,a,b) = t2(j,i,b,a) = sum_m coef(m) vec(o v + k, m), m = 1 .. nuse, k the packed position of the pair
      integer, intent(in) :: no, nv, nstore, nuse
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: vec(nstore, nuse), coef(nuse)
      real(dp), intent(inout) :: t2(n2)
      integer(int64) :: s, k, r
      real(dp) :: x

      do concurrent(s=1:n2) local(k, r, x)
         call pair_indices(TRC_DC_INDEX(s), no, nv, k, r)
         if (k > 0) then
            x = combine_item(no*nv + k, nstore, nuse, vec, coef)
            t2(s) = x
            t2(r) = x
         end if
      end do
   end subroutine kernel_extrapolate_unpack

   subroutine kernel_chunk_max_diff(n1, n2, npart, old1, old2, new1, new2, part)
      !! part(p) = max |new - old| over every npart-th element of the full concatenated vector, from the p-th
      integer, intent(in) :: n1, npart
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: old1(n1), old2(n2), new1(n1), new2(n2)
      real(dp), intent(out) :: part(npart)
      integer :: p

      do concurrent(p=1:npart)
         part(p) = chunk_max_diff_item(TRC_DC_INDEX(p), n1, n2, npart, old1, old2, new1, new2)
      end do
   end subroutine kernel_chunk_max_diff

   subroutine kernel_chunk_max(n, npart, e, part)
      !! part(p) = max |e| over every npart-th element, from the p-th
      integer, intent(in) :: n, npart
      real(dp), intent(in) :: e(n)
      real(dp), intent(out) :: part(npart)
      integer :: p

      do concurrent(p=1:npart)
         part(p) = chunk_max_item(TRC_DC_INDEX(p), n, npart, e)
      end do
   end subroutine kernel_chunk_max

   ! A call in a `do concurrent` body passes the loop index through the TRC_DC_INDEX macro:
   ! see inc/trc_dc_index.inc (an ifx miscompilation of the plain variable).

   pure function step_item(t, n1, n2, old1, old2, new1, new2) result(v)
      !! new - old at position t of the concatenated vector
      !$acc routine seq
      integer, intent(in) :: t
      integer, intent(in) :: n1, n2
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
      integer, intent(in) :: t
      integer, intent(in) :: n1, n2
      real(dp), intent(in) :: a1(n1), a2(n2)
      real(dp) :: v

      if (t <= n1) then
         v = a1(t)
      else
         v = a2(t - n1)
      end if
   end function join_item

   pure function chunk_max_item(p, n, npart, e) result(v)
      !! max |e| over elements p, p + npart, p + 2 npart, ... <= n: neighbouring partials read neighbouring elements
      !$acc routine seq
      integer, intent(in) :: p
      integer, intent(in) :: n, npart
      real(dp), intent(in) :: e(n)
      real(dp) :: v
      integer(int64) :: k

      v = 0.0_dp
      do k = p, n, npart
         v = max(v, abs(e(k)))
      end do
   end function chunk_max_item

   pure subroutine pair_indices(s, no, nv, k, r)
      !! For position s of t2(i,j,a,b) with p = (i,a), q = (j,b): k, the position of (p,q) in the packed lower half
      !! (column-major, p >= q), 0 for p < q; r, the position of the partner t2(j,i,b,a). One 64-bit division:
      !! the rest of the decoding fits default integers (o^2 and v^2 do).
      !$acc routine seq
      integer(int64), intent(in) :: s
      integer, intent(in) :: no, nv
      integer(int64), intent(out) :: k, r
      integer(int64) :: ab, p, q, nov
      integer :: ij, i, j, a, b

      ab = (s - 1)/(no*no)
      ij = int(s - 1 - ab*no*no)
      i = mod(ij, no)
      j = ij/no
      a = mod(int(ab), nv)
      b = int(ab)/nv
      nov = int(no, int64)*nv
      p = i + 1 + int(no, int64)*a
      q = j + 1 + int(no, int64)*b
      k = 0
      if (p >= q) k = p + (q - 1)*(2*nov - q)/2
      r = 1 + j + no*i + int(no*no, int64)*(b + int(nv, int64)*a)
   end subroutine pair_indices

   pure function packed_step_item(s, r, n2, old2, new2) result(v)
      !! (e(s) + e(r))/2 for e = new2 - old2 and r the pair partner of s, times sqrt(2) off the diagonal (r /= s)
      !$acc routine seq
      integer(int64), intent(in) :: s, r
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: old2(n2), new2(n2)
      real(dp) :: v

      v = 0.5_dp*((new2(s) - old2(s)) + (new2(r) - old2(r)))
      if (r /= s) v = SQRT2*v
   end function packed_step_item

   pure function combine_item(row, nstore, nuse, vec, coef) result(v)
      !! sum_m coef(m) vec(row, m), m = 1 .. nuse in ascending order
      !$acc routine seq
      integer(int64), intent(in) :: row
      integer, intent(in) :: nstore, nuse
      real(dp), intent(in) :: vec(nstore, nuse), coef(nuse)
      real(dp) :: v
      integer :: m

      v = 0.0_dp
      do m = 1, nuse
         v = v + coef(m)*vec(row, m)
      end do
   end function combine_item

   pure function chunk_max_diff_item(p, n1, n2, npart, old1, old2, new1, new2) result(v)
      !! max |new - old| over elements p, p + npart, p + 2 npart, ... of the concatenation (t1, t2)
      !$acc routine seq
      integer, intent(in) :: p
      integer, intent(in) :: n1, npart
      integer(int64), intent(in) :: n2
      real(dp), intent(in) :: old1(n1), old2(n2), new1(n1), new2(n2)
      real(dp) :: v
      integer(int64) :: k

      v = 0.0_dp
      do k = p, n1 + n2, npart
         if (k <= n1) then
            v = max(v, abs(new1(k) - old1(k)))
         else
            v = max(v, abs(new2(k - n1) - old2(k - n1)))
         end if
      end do
   end function chunk_max_diff_item

end module trc_rccsd_diis
