!
! The perturbative triples (T) of closed-shell RI-CCSD(T) on the DEVICE: the same quantity as
! `rccsd_triples_reference` (trc_rccsd_triples), with the same algebra, restructured only in how many triples go
! through each gemm and each kernel launch.
!
! THE ALGEBRA (unchanged, see trc_rccsd_triples). For each virtual triple a >= b >= c and each of the six
! permutations p = (u,v,w) of (a,b,c) the six W of Rendell, Lee and Komornicki
!
!     W_p(i,j,k) = sum_f (iu|fv) t2(k,j,w,f) - sum_m (ij|mu)' t2(m,k,v,w),      (ij|mu)' = ooov(j,m,i,u)
!
! then V, Z and the thirty-six-term sum of the WMAP table. The a = b = c triples are kept (they contribute zero for a
! closed shell, but keeping them keeps every per-pair partial equal to the host routine's).
!
! TILES. A tile is a run of consecutive triples in the order a = 1..nv, b = 1..a, c = 1..b (c fastest), as many as
! the memory budget lets the W tile hold: 6 no^3 doubles per triple (default budget 1 GiB, `budget_bytes` to change
! it, `tile_triples` to force a size in tests). A tile is cut into SEGMENTS (a, b, c0..c1), the pieces of one
! (a,b) pair it contains: at no = 25, nv = 130 a 1 GiB tile holds 1431 triples, about 11 pairs' worth at the
! top of the range. Per segment, with nc = c1 - c0 + 1, o = no, v = nv, the W blocks W(i,j,k) of the nc triples
! are stored contiguously (i fastest) for each p, so the particle product writes them in place and the hole
! product accumulates into them with no scatter pass. Twelve library calls per segment:
!
!   particle, W = ovvv-block x tt2, tt2(f,(j,k,w)) = t2(k,j,w,f) (one repack of t2, built once per call):
!     p = 1, 3  (c is the third index)  one gemm (o x v)(v x o^2 nc): the integral block ovvv(:,u,:,v) is read in
!               place (lda = o v, no copy), tt2 supplies the nc columns blocks of c0..c1
!     p = 2, 4  (c in the pair, u fixed)  one gemm_batched over c: A_c = ovvv(:,u,:,c) in place (stride o v^2),
!               the tt2 block of w shared, m = o, n = o^2, k = v, nc items
!     p = 5, 6  the same with u = c varying: A_c = ovvv(:,c,:,v), stride o
!   hole, W -= oop x t2, oop(i,j,m,u) = ooov(j,m,i,u) (one repack, o^3 v) and t2r(m,k,c,v) = t2(m,k,v,c):
!     p = 1..4  one gemm (o^2 x o)(o x o nc): the oop block of u shared, the t2 or t2r columns (k,c) of c0..c1 in place
!     p = 5, 6  one gemm_batched over c: A_c = oop(:,:,:,c), the t2 block shared, o^2 x o x o, nc items
!
! The batched particle products are the shape of the host routine's stacked (nc o x v)(v x o^2) product without
! its repack (mmt), skinny in m (= o) and large in n (= o^2) and in the batch; stacking them into one tall gemm
! needs the integral columns (i,c) contiguous, which ovvv does not give, and a tall output (i,c,j,k) which is not the
! per-triple block the next stage reads.
!
! ENERGY. Two `do concurrent` per tile. The first runs over the flattened (i, j, triple) index; its item routine
! loops over k and, for each (i,j,k), forms the six X = W + V/2, the six Z, the 36-term sum (the WMAP table
! written out: six calls of `triples_pair_term`, each with its row of the table) and returns the sum over k. It
! reads W in place and builds nothing: V comes from ovov and t1, D3 from eps. One thread per (i,j,triple) rather
! than per (i,j,k,triple) so that no atomic and no per-point partial array is needed, and the row sums fix their
! own order. The second sums the no^2 row values of each triple in index order. The per-triple energies are copied
! to the host and added into the per-pair array in ascending triple order; the pairs are then summed in pair order,
! exactly as the host routine does, so the result does not depend on the tile size beyond round-off in the gemms,
! nor on run-to-run thread scheduling.
!
! DEVICE MEMORY, doubles unless noted: tt2 and t2r each no^2 nv^2, oop no^3 nv, the W tile 6 no^3 T, the row sums
! no^2 T, the triple energies T, and 3 T integers, T the tile. For no = 25, nv = 130 and the default budget
! (T = 1431): 84.5 MB + 84.5 MB + 16.9 MB + 1.07 GB + 7.2 MB + negligible.
!
! RESIDENCY. The blocks of `ints` and the amplitudes of `amps` must be on the device and current (they are after
! `trc_rccsd_solve` and after `ints%build`; a caller that edits an amplitude on the host `amps%push()`es first). The
! scratch type maps its own arrays (`enter data create` in `init`, `exit data delete` in `release`, on the same
! expression) and records `ints%stamp`: `rccsd_triples_dev_run_scratch` refuses a scratch whose stamp is not the
! integrals' (after `ints%build` again, `init` anew). Each kernel gets every array as an explicit-shape dummy and
! names no `x%component`. The library's own cuBLAS handle is made and destroyed inside `rccsd_triples_dev_run`.
!
! MEMORY NOT IN THE BUDGET. `budget_bytes` counts the W tile only. The three repacks (tt2, t2r: 2 no^2 nv^2 doubles,
! oop: no^3 nv doubles) and the row sums are on top of it. A device allocation failure in `enter data create` is not a
! refusal but an abort, as everywhere else in the code base.
!
! INDEX RANGES. The W tile is indexed with 64-bit integers throughout (6 no^3 T can pass 2^31 for a large budget);
! the gemm dimensions are 32-bit and checked in `init` (no^2 nv, no nv^2 and the sizes of the repacked arrays,
! no^2 nv^2 and no^3 nv, below 2^31 - 1).
module trc_rccsd_triples_dev
   use pic_types, only: dp, int32, int64
   use pic_timer, only: timer_type
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   use trc_rccsd_ints, only: trc_rccsd_ints_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t
   implicit none
   private

   public :: trc_rccsd_triples_dev_options_t, trc_rccsd_triples_scratch_t
   public :: rccsd_triples_dev_run, rccsd_triples_dev_run_scratch, rccsd_triples_dev_point_energy

   integer(int64), parameter :: DEFAULT_BUDGET = 1073741824_int64  !! bytes of the W tile when `budget_bytes` is 0
   integer, parameter :: SEG_A = 1, SEG_B = 2, SEG_C0 = 3, SEG_NC = 4, SEG_S0 = 5  !! rows of the segment list

   type :: trc_rccsd_triples_dev_options_t
      !! Plain values only
      integer(int64) :: budget_bytes = 0   !! bytes the W tile may take (6 no^3 doubles a triple); 0 means 1 GiB
      integer(int64) :: tile_triples = 0   !! triples per tile; above 0 it replaces the budget (tests)
      integer(int64) :: slice_bytes = 0    !! occupied-driven routine only: bytes of cached ovvv slices; 0 means 4 GiB
      logical :: diagonal_triples = .false.  !! occupied-driven routine only: also compute the i = j = k triples (tests)
      logical :: verbose = .false.         !! the summary: tiles, triples, times, flops, GFLOP/s
   end type trc_rccsd_triples_dev_options_t

   type :: trc_rccsd_triples_scratch_t
      !! The device arrays of one (T) and the plan they were sized for, see the module header
      type(error_t) :: status   !! set by `init`: a refusal leaves its message here and the type empty
      integer :: no = 0
      integer :: nv = 0
      integer :: stamp = 0               !! `ints%stamp` of the integrals `init` planned for
      integer(int64) :: ncap = 0         !! triples per tile
      integer(int64) :: ntotal = 0       !! nv (nv + 1) (nv + 2) / 6
      integer(int64) :: bytes = 0        !! device bytes of the arrays below
      integer(int64) :: ntiles = 0       !! tiles of the last run
      real(dp) :: seconds_particle = 0.0_dp, seconds_hole = 0.0_dp, seconds_energy = 0.0_dp  !! device seconds, last run
      real(dp) :: flops_particle = 0.0_dp, flops_hole = 0.0_dp                               !! counted flops, last run
      real(dp), allocatable :: tt2(:)    !! tt2(f,(j,k,w)) = t2(k,j,w,f), nv x no^2 nv
      real(dp), allocatable :: t2r(:)    !! t2r(m,k,c,v) = t2(m,k,v,c), no^2 nv^2
      real(dp), allocatable :: oop(:)    !! oop(i,j,m,u) = ooov(j,m,i,u), no^3 nv
      real(dp), allocatable :: w(:)      !! the W tile, (no^3, ncap, 6): block (s,p) at ((p-1) ncap + s - 1) no^3
      real(dp), allocatable :: part(:)   !! row sums of the energy, no^2 ncap
      real(dp), allocatable :: etrip(:)  !! energy of each triple of the tile, ncap
      integer, allocatable :: trip(:, :)  !! (a,b,c) of each triple of the tile, 3 x ncap
   contains
      procedure :: init => scratch_init
      procedure :: release => scratch_release
   end type trc_rccsd_triples_scratch_t

contains

   ! ---------------------------------------------------------------------
   ! the scratch
   ! ---------------------------------------------------------------------

   subroutine scratch_init(this, ints, opts)
      !! Plans the tile and maps the arrays; a refusal is in `this%status` and leaves the type empty
      class(trc_rccsd_triples_scratch_t), intent(inout) :: this
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_triples_dev_options_t), intent(in) :: opts

      integer(int64) :: no, nv, no2, no3, budget, per_triple
      integer :: ierr

      call this%release()
      call this%status%clear()
      if (ints%stamp == 0 .or. .not. allocated(ints%ovvv) .or. .not. allocated(ints%ooov) .or. &
          .not. allocated(ints%ovov) .or. .not. allocated(ints%eps_o) .or. .not. allocated(ints%eps_v)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: the integral blocks are empty (no successful build)")
         return
      end if
      if (ints%no < 1 .or. ints%nv < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: no occupied or no virtual orbitals: there are no triples")
         return
      end if
      if (opts%budget_bytes < 0 .or. opts%tile_triples < 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: a negative budget or tile size")
         return
      end if
      no = ints%no
      nv = ints%nv
      no2 = no*no
      no3 = no2*no
      if (no2*nv > int(huge(0_int32), int64) .or. no*nv*nv > int(huge(0_int32), int64) .or. &
          no2*nv*nv > int(huge(0_int32), int64) .or. no3*nv > int(huge(0_int32), int64)) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: a gemm dimension or a repacked array would pass "// &
                              "2^31-1 elements ("//to_string(ints%no)//" occupied, "//to_string(ints%nv)//" virtual)")
         return
      end if
      this%no = ints%no
      this%nv = ints%nv
      this%ntotal = nv*(nv + 1)*(nv + 2)/6
      per_triple = 48_int64*no3
      budget = merge(opts%budget_bytes, DEFAULT_BUDGET, opts%budget_bytes > 0)
      if (opts%tile_triples > 0) then
         this%ncap = min(opts%tile_triples, this%ntotal)
      else
         this%ncap = min(budget/per_triple, this%ntotal)
      end if
      if (this%ncap < 1) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: the budget of "//to_string(budget)// &
                              " bytes does not hold one triple ("//to_string(per_triple)//" bytes of W)")
         this%no = 0
         this%nv = 0
         return
      end if
      allocate (this%tt2(nv*no2*nv), stat=ierr)
      if (ierr == 0) allocate (this%t2r(no2*nv*nv), stat=ierr)
      if (ierr == 0) allocate (this%oop(no3*nv), stat=ierr)
      if (ierr == 0) allocate (this%w(6_int64*no3*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%part(no2*this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%etrip(this%ncap), stat=ierr)
      if (ierr == 0) allocate (this%trip(3, this%ncap), stat=ierr)
      if (ierr /= 0) then
         call this%status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: the host arrays of the plan could not be allocated")
         call scratch_free(this)
         return
      end if
      !$acc enter data create(this%tt2, this%t2r, this%oop, this%w, this%part, this%etrip, this%trip)
      this%stamp = ints%stamp
      this%bytes = 8_int64*(2_int64*no2*nv*nv + no3*nv + 6_int64*no3*this%ncap + no2*this%ncap + this%ncap) + &
                   12_int64*this%ncap
   end subroutine scratch_init

   subroutine scratch_free(this)
      !! Frees what is allocated; the device copies only for the arrays `init` mapped (all or none)
      class(trc_rccsd_triples_scratch_t), intent(inout) :: this

      if (this%stamp /= 0) then
         !$acc exit data delete(this%tt2, this%t2r, this%oop, this%w, this%part, this%etrip, this%trip)
      end if
      if (allocated(this%tt2)) deallocate (this%tt2)
      if (allocated(this%t2r)) deallocate (this%t2r)
      if (allocated(this%oop)) deallocate (this%oop)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%part)) deallocate (this%part)
      if (allocated(this%etrip)) deallocate (this%etrip)
      if (allocated(this%trip)) deallocate (this%trip)
      this%no = 0
      this%nv = 0
      this%stamp = 0
      this%ncap = 0
      this%ntotal = 0
      this%bytes = 0
   end subroutine scratch_free

   subroutine scratch_release(this)
      !! Frees both copies; safe on an empty type
      class(trc_rccsd_triples_scratch_t), intent(inout) :: this

      call scratch_free(this)
   end subroutine scratch_release

   ! ---------------------------------------------------------------------
   ! the entry points
   ! ---------------------------------------------------------------------

   subroutine rccsd_triples_dev_run(ints, amps, e_t, status, opts, seconds, triple_energies)
      !! E(T) in Hartree on the device from the DEVICE copies of `ints` and `amps` (module header)
      !!
      !! Refused, with `e_t` = 0 and the reason in `status`: empty integrals, amplitudes missing or of another size
      !! than the integrals, no occupied or no virtual orbitals, a budget that does not hold one triple, a gemm
      !! dimension above 2^31-1, and a result that is NaN or infinite. `triple_energies`, when present, is
      !! allocated to nv (nv + 1) (nv + 2) / 6 in the order of the virtual triples a >= b >= c (a = 1..nv, b = 1..a,
      !! c = 1..b, c fastest); each entry is the weighted energy of that one triple, and E(T) = 2 * their sum taken
      !! in that order.
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      real(dp), intent(out) :: e_t
      type(error_t), intent(out) :: status
      type(trc_rccsd_triples_dev_options_t), intent(in), optional :: opts
      real(dp), intent(out), optional :: seconds  !! device wall seconds of the whole call
      real(dp), intent(out), allocatable, optional :: triple_energies(:)

      type(trc_rccsd_triples_dev_options_t) :: o
      type(trc_rccsd_triples_scratch_t) :: ts
      type(timer_type) :: clock
      type(trc_linalg_t) :: la

      e_t = 0.0_dp
      if (present(seconds)) seconds = 0.0_dp
      if (present(opts)) o = opts
      call la%init(1)
      call clock%start()
      call ts%init(ints, o)
      if (ts%status%has_error()) then
         call status%set(ERROR_VALIDATION, ts%status%get_message())
      else
         call rccsd_triples_dev_run_scratch(ints, amps, ts, la, e_t, status, triple_energies)
         call la%sync()
         if (present(seconds)) seconds = clock%get_elapsed_time()
         if (o%verbose .and. .not. status%has_error()) call log_summary(ints, ts, e_t, clock%get_elapsed_time())
      end if
      call ts%release()
      call la%release()
   end subroutine rccsd_triples_dev_run

   subroutine rccsd_triples_dev_run_scratch(ints, amps, ts, la, e_t, status, triple_energies)
      !! One (T) with a scratch that `ts%init` planned for these integrals; see `rccsd_triples_dev_run`
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(out) :: e_t
      type(error_t), intent(out) :: status
      real(dp), intent(out), allocatable, optional :: triple_energies(:)

      real(dp), allocatable :: e_pair(:)
      integer, allocatable :: seg(:, :)
      character(len=:), allocatable :: why
      integer(int64) :: done, ipair, npair
      integer :: cursor(3)
      integer :: nt, nseg, s, a, b

      e_t = 0.0_dp
      why = inputs_refusal(ints, amps, ts)
      if (len(why) > 0) then
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: "//why)
         return
      end if
      call pack_inputs(ints, amps, ts)
      npair = int(ts%nv, int64)*(ts%nv + 1)/2
      allocate (e_pair(npair), seg(5, min(ts%ncap, npair)))
      e_pair = 0.0_dp
      if (present(triple_energies)) then
         if (allocated(triple_energies)) deallocate (triple_energies)
         allocate (triple_energies(ts%ntotal))
      end if
      ts%seconds_particle = 0.0_dp
      ts%seconds_hole = 0.0_dp
      ts%seconds_energy = 0.0_dp
      ts%flops_particle = 0.0_dp
      ts%flops_hole = 0.0_dp
      ts%ntiles = 0
      cursor = [1, 1, 1]
      done = 0
      do while (done < ts%ntotal)
         nt = int(min(ts%ncap, ts%ntotal - done))
         call tile_segments(ts%nv, nt, cursor, seg, nseg, ts%trip)
         !$acc update device(ts%trip(1:3, 1:nt))
         call tile_products(ints, amps, ts, la, seg, nseg, nt)
         call tile_energy(ints, amps, ts, la, nt)
         !$acc update self(ts%etrip(1:nt))
         do s = 1, nt
            a = ts%trip(1, s)
            b = ts%trip(2, s)
            ipair = int(a, int64)*(a - 1)/2 + b
            e_pair(ipair) = e_pair(ipair) + ts%etrip(s)
         end do
         if (present(triple_energies)) triple_energies(done + 1:done + nt) = ts%etrip(1:nt)
         done = done + nt
         ts%ntiles = ts%ntiles + 1
      end do
      ! Summed in pair order, as the host routine does, so the total is the same for every tile size up to the gemms.
      do ipair = 1, npair
         e_t = e_t + e_pair(ipair)
      end do
      e_t = 2.0_dp*e_t
      if (.not. abs(e_t) <= huge(1.0_dp)) then
         e_t = 0.0_dp
         call status%set(ERROR_VALIDATION, "trc_rccsd_triples_dev: E(T) is NaN or infinite (a zero denominator or a "// &
                         "non-finite amplitude or integral)")
      end if
      deallocate (e_pair, seg)
   end subroutine rccsd_triples_dev_run_scratch

   function inputs_refusal(ints, amps, ts) result(why)
      !! Empty when the inputs can be used with the scratch, else the reason
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_scratch_t), intent(in) :: ts
      character(len=:), allocatable :: why

      why = ""
      if (ts%stamp == 0 .or. .not. allocated(ts%w)) then
         why = "the scratch is not initialised"
      else if (ints%stamp == 0 .or. .not. allocated(ints%ovvv) .or. .not. allocated(ints%ooov) .or. &
               .not. allocated(ints%ovov) .or. .not. allocated(ints%eps_o) .or. .not. allocated(ints%eps_v)) then
         why = "the integral blocks are empty (no successful build)"
      else if (ts%stamp /= ints%stamp) then
         why = "the scratch was planned for an earlier build of the integrals (stale stamp): init it again"
      else if (ts%no /= ints%no .or. ts%nv /= ints%nv) then
         why = "the scratch has the sizes of other integrals"
      else if (.not. allocated(amps%t1) .or. .not. allocated(amps%t2)) then
         why = "the amplitudes are not allocated"
      else if (amps%no /= ints%no .or. amps%nv /= ints%nv) then
         why = "the amplitudes are of another size than the integrals"
      else if (any(shape(amps%t1) /= [ints%no, ints%nv]) .or. any(shape(amps%t2) /= [ints%no, ints%no, ints%nv, ints%nv])) then
         why = "the amplitude arrays do not have the shape (no, nv) and (no, no, nv, nv)"
      else if (any(shape(ints%ovvv) /= [ints%no, ints%nv, ints%nv, ints%nv]) .or. &
               any(shape(ints%ooov) /= [ints%no, ints%no, ints%no, ints%nv]) .or. &
               any(shape(ints%ovov) /= [ints%no, ints%nv, ints%no, ints%nv]) .or. &
               size(ints%eps_o) /= ints%no .or. size(ints%eps_v) /= ints%nv) then
         why = "the integral blocks do not have the shapes of (no, nv)"
      end if
      if (len(why) > 0) why = trim(why)//" ("//to_string(ints%no)//" occupied, "//to_string(ints%nv)//" virtual)"
   end function inputs_refusal

   ! ---------------------------------------------------------------------
   ! one call's steps
   ! ---------------------------------------------------------------------

   subroutine pack_inputs(ints, amps, ts)
      !! tt2, t2r and oop from the device amplitudes and integrals, once per call
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_scratch_t), intent(inout) :: ts

      call kernel_pack_tt2(ts%no, ts%nv, amps%t2, ts%tt2)
      call kernel_pack_t2r(ts%no, ts%nv, amps%t2, ts%t2r)
      call kernel_pack_oop(ts%no, ts%nv, ints%ooov, ts%oop)
   end subroutine pack_inputs

   subroutine tile_segments(nv, nt, cursor, seg, nseg, trip)
      !! The next `nt` triples from `cursor` = (a,b,c) as segments, and their (a,b,c) in `trip`
      integer, intent(in) :: nv, nt
      integer, intent(inout) :: cursor(3)
      integer, intent(out) :: seg(:, :)
      integer, intent(out) :: nseg
      integer, intent(out) :: trip(:, :)

      integer :: left, nc, s, c

      nseg = 0
      s = 0
      left = nt
      do while (left > 0)
         nc = min(left, cursor(2) - cursor(3) + 1)
         nseg = nseg + 1
         seg(SEG_A, nseg) = cursor(1)
         seg(SEG_B, nseg) = cursor(2)
         seg(SEG_C0, nseg) = cursor(3)
         seg(SEG_NC, nseg) = nc
         seg(SEG_S0, nseg) = s
         do c = cursor(3), cursor(3) + nc - 1
            s = s + 1
            trip(:, s) = [cursor(1), cursor(2), c]
         end do
         left = left - nc
         cursor(3) = cursor(3) + nc
         if (cursor(3) > cursor(2)) then
            cursor(3) = 1
            cursor(2) = cursor(2) + 1
            if (cursor(2) > cursor(1)) then
               cursor(2) = 1
               cursor(1) = min(cursor(1) + 1, nv + 1)
            end if
         end if
      end do
   end subroutine tile_segments

   subroutine tile_products(ints, amps, ts, la, seg, nseg, nt)
      !! The particle and the hole products of every segment of the tile, timed apart
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: seg(:, :)
      integer, intent(in) :: nseg, nt

      type(timer_type) :: clock
      real(dp) :: o, v

      o = real(ts%no, dp)
      v = real(ts%nv, dp)
      call la%sync()
      call clock%start()
      call particle_products(ints, ts, la, seg, nseg)
      call la%sync()
      ts%seconds_particle = ts%seconds_particle + clock%get_elapsed_time()
      call clock%start()
      call hole_products(amps, ts, la, seg, nseg)
      call la%sync()
      ts%seconds_hole = ts%seconds_hole + clock%get_elapsed_time()
      ts%flops_particle = ts%flops_particle + 12.0_dp*o*o*o*v*real(nt, dp)
      ts%flops_hole = ts%flops_hole + 12.0_dp*o*o*o*o*real(nt, dp)
   end subroutine tile_products

   subroutine tile_energy(ints, amps, ts, la, nt)
      !! The two energy kernels over the `nt` triples of the tile; the sums stay on the device
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: nt

      type(timer_type) :: clock

      call la%sync()
      call clock%start()
      call kernel_energy_rows(ts%no, ts%nv, nt, ts%ncap, ts%w, ints%ovov, amps%t1, ints%eps_o, ints%eps_v, ts%trip, ts%part)
      call kernel_energy_sum(ts%no, nt, ts%ncap, ts%part, ts%etrip)
      call la%sync()
      ts%seconds_energy = ts%seconds_energy + clock%get_elapsed_time()
   end subroutine tile_energy

   pure function wpos(ts_ncap, no, p, s) result(pos)
      !! Position (1-based) in the W tile of block (triple s, permutation p), s counted from 0
      integer(int64), intent(in) :: ts_ncap
      integer, intent(in) :: no, p, s
      integer(int64) :: pos

      pos = (int(p - 1, int64)*ts_ncap + s)*int(no, int64)**3 + 1_int64
   end function wpos

   subroutine particle_products(ints, ts, la, seg, nseg)
      !! W_p = ovvv block x tt2 of every segment: 2 gemms and 4 batched gemms each, beta = 0
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_triples_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: seg(:, :)
      integer, intent(in) :: nseg

      type(trc_gemm_batch_t) :: by_a, by_c
      integer :: g, a, b, c0, nc, s0, o, nv
      integer(int64) :: o2, o3, tt_c0, tt_a, tt_b

      o = ts%no
      nv = ts%nv
      o2 = int(o, int64)*o
      o3 = o2*o
      do g = 1, nseg
         a = seg(SEG_A, g)
         b = seg(SEG_B, g)
         c0 = seg(SEG_C0, g)
         nc = seg(SEG_NC, g)
         s0 = seg(SEG_S0, g)
         tt_c0 = 1_int64 + int(nv, int64)*o2*(c0 - 1)
         tt_a = 1_int64 + int(nv, int64)*o2*(a - 1)
         tt_b = 1_int64 + int(nv, int64)*o2*(b - 1)
         ! p = 1 (a,b,c) and 3 (b,a,c): the integral block is fixed, the c blocks of tt2 stack as columns
         call la%gemm("N", "N", o, int(o2)*nc, nv, 1.0_dp, ints%ovvv(1, a, 1, b), o*nv, ts%tt2(tt_c0), nv, 0.0_dp, &
                      ts%w(wpos(ts%ncap, o, 1, s0)), o)
         call la%gemm("N", "N", o, int(o2)*nc, nv, 1.0_dp, ints%ovvv(1, b, 1, a), o*nv, ts%tt2(tt_c0), nv, 0.0_dp, &
                      ts%w(wpos(ts%ncap, o, 3, s0)), o)
         ! p = 2 (a,c,b), 4 (b,c,a): u fixed, v = c varies (A stride o nv^2), the tt2 block of w shared
         by_a = trc_gemm_batch(m=o, n=int(o2), k=nv, count=nc, lda=o*nv, ldb=nv, ldc=o, &
                               sa=int(o, int64)*nv*nv, sb=0_int64, sc=o3)
         call la%gemm_batched(by_a, ints%ovvv(1, a, 1, c0), ts%tt2(tt_b), ts%w(wpos(ts%ncap, o, 2, s0)))
         call la%gemm_batched(by_a, ints%ovvv(1, b, 1, c0), ts%tt2(tt_a), ts%w(wpos(ts%ncap, o, 4, s0)))
         ! p = 5 (c,a,b), 6 (c,b,a): u = c varies (A stride o), v fixed
         by_c = trc_gemm_batch(m=o, n=int(o2), k=nv, count=nc, lda=o*nv, ldb=nv, ldc=o, &
                               sa=int(o, int64), sb=0_int64, sc=o3)
         call la%gemm_batched(by_c, ints%ovvv(1, c0, 1, a), ts%tt2(tt_b), ts%w(wpos(ts%ncap, o, 5, s0)))
         call la%gemm_batched(by_c, ints%ovvv(1, c0, 1, b), ts%tt2(tt_a), ts%w(wpos(ts%ncap, o, 6, s0)))
      end do
   end subroutine particle_products

   subroutine hole_products(amps, ts, la, seg, nseg)
      !! W_p -= oop x t2 of every segment: 4 gemms and 2 batched gemms each, alpha = -1, beta = 1
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_triples_scratch_t), intent(inout) :: ts
      type(trc_linalg_t), intent(in) :: la
      integer, intent(in) :: seg(:, :)
      integer, intent(in) :: nseg

      type(trc_gemm_batch_t) :: by_c
      integer :: g, a, b, c0, nc, s0, o, nv
      integer(int64) :: o2, o3, r_c0b, r_c0a, oop_a, oop_b, oop_c

      o = ts%no
      nv = ts%nv
      o2 = int(o, int64)*o
      o3 = o2*o
      by_c = trc_gemm_batch(m=int(o2), n=o, k=o, count=1, lda=int(o2), ldb=o, ldc=int(o2), sa=o3, sb=0_int64, sc=o3)
      do g = 1, nseg
         a = seg(SEG_A, g)
         b = seg(SEG_B, g)
         c0 = seg(SEG_C0, g)
         nc = seg(SEG_NC, g)
         s0 = seg(SEG_S0, g)
         oop_a = 1_int64 + o3*(a - 1)
         oop_b = 1_int64 + o3*(b - 1)
         oop_c = 1_int64 + o3*(c0 - 1)
         r_c0b = 1_int64 + o2*(c0 - 1) + o2*nv*(b - 1)
         r_c0a = 1_int64 + o2*(c0 - 1) + o2*nv*(a - 1)
         ! p = 1 (a,b,c): t2(m,k,b,c) = t2r(m,k,c,b); p = 3 (b,a,c): t2r(m,k,c,a); the oop block of u shared
         call la%gemm("N", "N", int(o2), o*nc, o, -1.0_dp, ts%oop(oop_a), int(o2), ts%t2r(r_c0b), o, 1.0_dp, &
                      ts%w(wpos(ts%ncap, o, 1, s0)), int(o2))
         call la%gemm("N", "N", int(o2), o*nc, o, -1.0_dp, ts%oop(oop_b), int(o2), ts%t2r(r_c0a), o, 1.0_dp, &
                      ts%w(wpos(ts%ncap, o, 3, s0)), int(o2))
         ! p = 2 (a,c,b): t2(m,k,c,b); p = 4 (b,c,a): t2(m,k,c,a), both in place in t2
         call la%gemm("N", "N", int(o2), o*nc, o, -1.0_dp, ts%oop(oop_a), int(o2), amps%t2(1, 1, c0, b), o, 1.0_dp, &
                      ts%w(wpos(ts%ncap, o, 2, s0)), int(o2))
         call la%gemm("N", "N", int(o2), o*nc, o, -1.0_dp, ts%oop(oop_b), int(o2), amps%t2(1, 1, c0, a), o, 1.0_dp, &
                      ts%w(wpos(ts%ncap, o, 4, s0)), int(o2))
         ! p = 5 (c,a,b), 6 (c,b,a): u = c varies, the t2 block (a,b) or (b,a) shared
         by_c%count = nc
         call la%gemm_batched(by_c, ts%oop(oop_c), amps%t2(1, 1, a, b), ts%w(wpos(ts%ncap, o, 5, s0)), &
                              alpha=-1.0_dp, beta=1.0_dp)
         call la%gemm_batched(by_c, ts%oop(oop_c), amps%t2(1, 1, b, a), ts%w(wpos(ts%ncap, o, 6, s0)), &
                              alpha=-1.0_dp, beta=1.0_dp)
      end do
   end subroutine hole_products

   subroutine log_summary(ints, ts, e_t, seconds)
      !! The verbose table: tiles, triples, device time by stage, counted flops and the rates
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_triples_scratch_t), intent(in) :: ts
      real(dp), intent(in) :: e_t, seconds

      character(len=240) :: line
      real(dp) :: gemm_seconds

      gemm_seconds = ts%seconds_particle + ts%seconds_hole
      write (line, "(a,i0,a,i0,a,f20.12,a,f9.2,a)") "  rccsd (T) device: no ", ints%no, " nv ", ints%nv, "   E(T) ", e_t, &
         "   ", seconds, " s"
      call global_logger%info(trim(line))
      write (line, "(a,i0,a,i0,a,i0,a,f9.1,a)") "    ", ts%ntiles, " tiles of up to ", ts%ncap, " triples (", ts%ntotal, &
         " in all), device scratch ", real(ts%bytes, dp)/1.0e6_dp, " MB"
      call global_logger%info(trim(line))
      write (line, "(a,f9.3,a,f9.3,a,f9.3,a)") "    seconds: particle gemms ", ts%seconds_particle, "  hole gemms ", &
         ts%seconds_hole, "  energy kernels ", ts%seconds_energy, ""
      call global_logger%info(trim(line))
      write (line, "(a,es10.3,a,f9.1,a,es10.3,a,f9.1,a,f9.1,a)") "    flops: particle ", ts%flops_particle, " (", &
         1.0e-9_dp*ts%flops_particle/max(ts%seconds_particle, tiny(1.0_dp)), " GFLOP/s)  hole ", ts%flops_hole, " (", &
         1.0e-9_dp*ts%flops_hole/max(ts%seconds_hole, tiny(1.0_dp)), " GFLOP/s)  both ", &
         1.0e-9_dp*(ts%flops_particle + ts%flops_hole)/max(gemm_seconds, tiny(1.0_dp)), " GFLOP/s"
      call global_logger%info(trim(line))
   end subroutine log_summary

   ! ---------------------------------------------------------------------
   ! the kernels. Each takes its arrays as explicit-shape dummies and runs one
   ! `do concurrent` over ONE flattened index, whose body calls an item routine
   ! of this module and assigns one array element. No `x%name` appears in any.
   ! ---------------------------------------------------------------------

   subroutine kernel_pack_tt2(no, nv, t2, tt2)
      !! tt2(f,(j,k,w)) = t2(k,j,w,f)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no*no*nv*nv)
      real(dp), intent(out) :: tt2(nv*no*no*nv)
      integer(int64) :: t, n

      n = int(nv, int64)*no*no*nv
      do concurrent(t=1:n)
         tt2(t) = tt2_item(t, no, nv, t2)
      end do
   end subroutine kernel_pack_tt2

   subroutine kernel_pack_t2r(no, nv, t2, t2r)
      !! t2r(m,k,c,v) = t2(m,k,v,c)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no*no*nv*nv)
      real(dp), intent(out) :: t2r(no*no*nv*nv)
      integer(int64) :: t, n

      n = int(no, int64)*no*nv*nv
      do concurrent(t=1:n)
         t2r(t) = t2r_item(t, no, nv, t2)
      end do
   end subroutine kernel_pack_t2r

   subroutine kernel_pack_oop(no, nv, ooov, oop)
      !! oop(i,j,m,u) = ooov(j,m,i,u)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: ooov(no*no*no*nv)
      real(dp), intent(out) :: oop(no*no*no*nv)
      integer(int64) :: t, n

      n = int(no, int64)*no*no*nv
      do concurrent(t=1:n)
         oop(t) = oop_item(t, no, nv, ooov)
      end do
   end subroutine kernel_pack_oop

   subroutine kernel_energy_rows(no, nv, nt, ncap, w, ovov, t1, eps_o, eps_v, trip, part)
      !! part(i,j,s) = sum_k of the energy of point (i,j,k) of triple s, s = 1..nt
      integer, intent(in) :: no, nv, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: w(6_int64*no*no*no*ncap)
      real(dp), intent(in) :: ovov(no, nv, no, nv), t1(no, nv), eps_o(no), eps_v(nv)
      integer, intent(in) :: trip(3, ncap)
      real(dp), intent(out) :: part(no*no*ncap)
      integer(int64) :: t, n

      n = int(no, int64)*no*nt
      do concurrent(t=1:n)
         part(t) = row_energy_item(t, no, nv, ncap, w, ovov, t1, eps_o, eps_v, trip)
      end do
   end subroutine kernel_energy_rows

   subroutine kernel_energy_sum(no, nt, ncap, part, etrip)
      !! etrip(s) = sum of the no^2 row energies of triple s, in index order
      integer, intent(in) :: no, nt
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: part(no*no*ncap)
      real(dp), intent(out) :: etrip(ncap)
      integer :: s

      do concurrent(s=1:nt)
         etrip(s) = sum_item(s, no, ncap, part)
      end do
   end subroutine kernel_energy_sum

   ! ---------------------------------------------------------------------
   ! the item routines: pure, `routine seq`, explicit-shape dummies, scalars
   ! only. A flat 1-based index t is the position in the array that is written,
   ! first index fastest.
   ! ---------------------------------------------------------------------

   pure subroutine split4(t, n1, n2, n3, i1, i2, i3, i4)
      !! t - 1 = (i1-1) + n1 ((i2-1) + n2 ((i3-1) + n3 (i4-1))), 64-bit t
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: n1, n2, n3
      integer, intent(out) :: i1, i2, i3, i4
      integer(int64) :: q

      q = t - 1
      i1 = int(mod(q, int(n1, int64))) + 1
      q = q/n1
      i2 = int(mod(q, int(n2, int64))) + 1
      q = q/n2
      i3 = int(mod(q, int(n3, int64))) + 1
      i4 = int(q/n3) + 1
   end subroutine split4

   pure function tt2_item(t, no, nv, t2) result(v)
      !! tt2(f,(j,k,w)) = t2(k,j,w,f), position t of (nv, no, no, nv)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no*no*nv*nv)
      real(dp) :: v
      integer :: f, j, k, w

      call split4(t, nv, no, no, f, j, k, w)
      v = t2(k + no*(j - 1) + int(no, int64)*no*(w - 1) + int(no, int64)*no*nv*(f - 1))
   end function tt2_item

   pure function t2r_item(t, no, nv, t2) result(v)
      !! t2r(m,k,c,vv) = t2(m,k,vv,c), position t of (no, no, nv, nv)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t2(no*no*nv*nv)
      real(dp) :: v
      integer :: m, k, c, vv

      call split4(t, no, no, nv, m, k, c, vv)
      v = t2(m + no*(k - 1) + int(no, int64)*no*(vv - 1) + int(no, int64)*no*nv*(c - 1))
   end function t2r_item

   pure function oop_item(t, no, nv, ooov) result(v)
      !! oop(i,j,m,u) = ooov(j,m,i,u), position t of (no, no, no, nv)
      !$acc routine seq
      integer(int64), intent(in) :: t
      integer, intent(in) :: no
      integer, intent(in) :: nv
      real(dp), intent(in) :: ooov(no*no*no*nv)
      real(dp) :: v
      integer :: i, j, m, u

      call split4(t, no, no, no, i, j, m, u)
      v = ooov(j + no*(m - 1) + int(no, int64)*no*(i - 1) + int(no, int64)*no*no*(u - 1))
   end function oop_item

   pure function sum_item(s, no, ncap, part) result(e)
      !! The no^2 row energies of triple s, added in index order
      !$acc routine seq
      integer, intent(in) :: s, no
      integer(int64), intent(in) :: ncap
      real(dp), intent(in) :: part(no*no*ncap)
      real(dp) :: e
      integer :: r

      e = 0.0_dp
      do r = 1, no*no
         e = e + part(int(s - 1, int64)*no*no + r)
      end do
   end function sum_item

   pure function row_energy_item(t, no, nv, ncap, w, ovov, t1, eps_o, eps_v, trip) result(e)
      !! The energy of the points (i, j, 1..no) of one triple: position t of (no[i], no[j], nt[s])
      !$acc routine seq
      integer(int64), intent(in) :: t, ncap
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: w(6_int64*no*no*no*ncap)
      real(dp), intent(in) :: ovov(no, nv, no, nv), t1(no, nv), eps_o(no), eps_v(nv)
      integer, intent(in) :: trip(3, ncap)
      real(dp) :: e
      integer :: i, j, k, s, unused
      integer(int64) :: no3

      call split4(t, no, no, 1, i, j, unused, s)
      no3 = int(no, int64)**3
      e = 0.0_dp
      do k = 1, no
         e = e + rccsd_triples_dev_point_energy(i, j, k, no, nv, no3*(s - 1), no3*ncap, trip(1, s), trip(2, s), trip(3, s), &
                                                w, 6_int64*no3*ncap, ovov, t1, eps_o, eps_v)
      end do
   end function row_energy_item

   pure function rccsd_triples_dev_point_energy(i, j, k, no, nv, sbase, pstride, a, b, c, w, wn, ovov, t1, eps_o, eps_v) &
      result(e)
      !! sum_p Z_p(i,j,k) sum_q W_wmap(p,q)(perm_q(i,j,k)) for the triple (a,b,c) whose W blocks start `sbase`
      !! elements into each permutation's plane of `w` (planes `pstride` apart): one point of the host routine's
      !! `triples_energy_block`, with the six X = W + V/2, the six Z and D3 formed on the fly
      !$acc routine seq
      integer, intent(in) :: i, j, k, no, nv, a, b, c
      integer(int64), intent(in) :: sbase, pstride, wn
      real(dp), intent(in) :: w(wn)
      real(dp), intent(in) :: ovov(no, nv, no, nv), t1(no, nv), eps_o(no), eps_v(nv)
      real(dp) :: e
      real(dp) :: d3, scale

      scale = 1.0_dp
      if (a == c) then
         scale = 6.0_dp          ! a == b == c
      else if (a == b .or. b == c) then
         scale = 2.0_dp
      end if
      d3 = scale*(eps_o(i) + eps_o(j) + eps_o(k) - eps_v(a) - eps_v(b) - eps_v(c))
      ! The WMAP table of trc_rccsd_triples written out: row p, the W plane paired with each occupied permutation
      ! (ijk, ikj, jik, jki, kij, kji); the virtual permutation p is (u,v,w) as the host routine's `perm`.
      e = triples_pair_term(1, 1, 2, 3, 4, 5, 6, a, b, c, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, t1)
      e = e + triples_pair_term(2, 2, 1, 5, 6, 3, 4, a, c, b, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, t1)
      e = e + triples_pair_term(3, 3, 4, 1, 2, 6, 5, b, a, c, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, t1)
      e = e + triples_pair_term(4, 4, 3, 6, 5, 1, 2, b, c, a, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, t1)
      e = e + triples_pair_term(5, 5, 6, 2, 1, 4, 3, c, a, b, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, t1)
      e = e + triples_pair_term(6, 6, 5, 4, 3, 2, 1, c, b, a, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, t1)
   end function rccsd_triples_dev_point_energy

   pure function triples_pair_term(p, r1, r2, r3, r4, r5, r6, u, v, ww, i, j, k, d3, no, nv, sbase, pstride, w, wn, ovov, &
                                   t1) result(term)
      !! Z_p(i,j,k) [W_r1(ijk) + W_r2(ikj) + W_r3(jik) + W_r4(jki) + W_r5(kij) + W_r6(kji)], (r1..r6) row p of WMAP
      !!
      !! Z_p = [4 X(ijk) + X(kij) + X(jki) - 2 X(kji) - 2 X(ikj) - 2 X(jik)] / D3 with X = W_p + V_p / 2 and
      !! V_p(x,y,z) = ovov(x,u,y,v) t1(z,ww): the host routine's `triples_v` and `triples_r3`, same order of operations
      !$acc routine seq
      integer, intent(in) :: p, r1, r2, r3, r4, r5, r6, u, v, ww, i, j, k, no, nv
      integer(int64), intent(in) :: sbase, pstride, wn
      real(dp), intent(in) :: d3
      real(dp), intent(in) :: w(wn)
      real(dp), intent(in) :: ovov(no, nv, no, nv), t1(no, nv)
      real(dp) :: term
      real(dp) :: x1, x2, x3, x4, x5, x6, z, s
      integer(int64) :: off

      off = int(p - 1, int64)*pstride + sbase
      x1 = triples_x_value(off, u, v, ww, i, j, k, no, nv, w, wn, ovov, t1)
      x2 = triples_x_value(off, u, v, ww, i, k, j, no, nv, w, wn, ovov, t1)
      x3 = triples_x_value(off, u, v, ww, j, i, k, no, nv, w, wn, ovov, t1)
      x4 = triples_x_value(off, u, v, ww, j, k, i, no, nv, w, wn, ovov, t1)
      x5 = triples_x_value(off, u, v, ww, k, i, j, no, nv, w, wn, ovov, t1)
      x6 = triples_x_value(off, u, v, ww, k, j, i, no, nv, w, wn, ovov, t1)
      z = (4.0_dp*x1 + x5 + x4 - 2.0_dp*x6 - 2.0_dp*x2 - 2.0_dp*x3)/d3
      s = triples_w_value(int(r1 - 1, int64)*pstride + sbase, i, j, k, no, w, wn) &
          + triples_w_value(int(r2 - 1, int64)*pstride + sbase, i, k, j, no, w, wn) &
          + triples_w_value(int(r3 - 1, int64)*pstride + sbase, j, i, k, no, w, wn) &
          + triples_w_value(int(r4 - 1, int64)*pstride + sbase, j, k, i, no, w, wn) &
          + triples_w_value(int(r5 - 1, int64)*pstride + sbase, k, i, j, no, w, wn) &
          + triples_w_value(int(r6 - 1, int64)*pstride + sbase, k, j, i, no, w, wn)
      term = z*s
   end function triples_pair_term

   pure function triples_w_value(off, x, y, z, no, w, wn) result(value)
      !! W(x,y,z) of the block that starts `off` elements into w
      !$acc routine seq
      integer(int64), intent(in) :: off, wn
      integer, intent(in) :: x, y, z, no
      real(dp), intent(in) :: w(wn)
      real(dp) :: value

      value = w(off + 1_int64 + (x - 1) + int(no, int64)*(y - 1) + int(no, int64)*no*(z - 1))
   end function triples_w_value

   pure function triples_x_value(off, u, v, ww, x, y, z, no, nv, w, wn, ovov, t1) result(value)
      !! X(x,y,z) = W(x,y,z) + 0.5 V(x,y,z), V = ovov(x,u,y,v) t1(z,ww)
      !$acc routine seq
      integer(int64), intent(in) :: off, wn
      integer, intent(in) :: u, v, ww, x, y, z, no, nv
      real(dp), intent(in) :: w(wn)
      real(dp), intent(in) :: ovov(no, nv, no, nv), t1(no, nv)
      real(dp) :: value

      value = triples_w_value(off, x, y, z, no, w, wn) + 0.5_dp*(ovov(x, u, y, v)*t1(z, ww))
   end function triples_x_value

end module trc_rccsd_triples_dev
