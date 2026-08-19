!
! Binned shell-pair container: quartets formed ON THE DEVICE.
!
! WHY THIS REPLACES THE WORK LIST
! -------------------------------
! At 75 waters / 6-31G the materialised quartet list cost 26 s of host time and
! 11.7 GB, against 15 s of GPU.  The host had become the bottleneck and the
! memory did not extrapolate: 150 waters would have wanted ~47 GB and ~105 s.
!
! The fix, from Galvez Vallejo et al. (EXESS), is to never enumerate quartets.
! Keep the O(N^2) shell-pair list, bin it, and let each kernel launch cover the
! cross product of two bins -- a thread derives its own (ab, cd) from its index.
! Quartet formation moves to the device and costs nothing to store.
!
! THE SCREENING TRICK
! -------------------
! Bin pairs by TYPE (la, lb) and by SIZE, the Schwarz bound bucketed into
! decades:
!
!     s_ab = int(-log10(Q_ab))          Q_ab = sqrt(max |(ab|ab)|)
!
! Then, since |(ab|cd)| <= Q_ab Q_cd,
!
!     keep the bin pair  <=>  s_ab + s_cd <= -log10(tau)
!
! is an INTEGER test on two bin indices, evaluated once per bin pair on the
! host.  A whole cross product is admitted or rejected together, so the kernel
! contains no screening conditional at all -- which is the point EXESS makes:
! conditionals in the kernel cost thread divergence, and here there are none.
!
! Truncation direction matters and is safe: int() rounds toward zero, so
! s_ab + s_cd <= -log10(Q_ab Q_cd).  A bin pair that should be kept therefore
! always passes.  The test is conservative -- it can admit slightly more than
! the threshold demands, never less.
!
! Binning by type as well as size is not free of charge either: it means every
! thread in a launch shares (la, lb, lc, ld), which is the class-sorting win
! from earlier arriving as a side effect of the data structure.
!
!
! [1] is the reference in README.md; see its `References` section.
!
module trc_bins
   use trc_boys, only: dp
   use trc_tables, only: LMAX
   use pic_types, only: int64, int_index
   use pic_sorting, only: sort_index
   implicit none
   private

   public :: SMAX, bin_key, build_binned_pairs, pair_bins_t, bin_dmax

   integer, parameter :: SMAX = 9    !! deepest size bucket, following [1]
   integer, parameter :: KCAP = 64   !! largest contraction degree binned exactly

   type :: pair_bins_t
      integer :: nbin = 0
      integer :: npair = 0
      integer, allocatable :: sp_i(:), sp_j(:)     !! shell indices, bin-sorted
      real(dp), allocatable :: sp_q(:)             !! Schwarz bound per pair
      integer, allocatable :: bin_off(:)           !! first pair of bin b
      integer, allocatable :: bin_cnt(:)           !! pairs in bin b
      integer, allocatable :: bin_s(:)             !! size index of bin b
      integer, allocatable :: bin_la(:), bin_lb(:) !! type of bin b
      integer, allocatable :: live(:)              !! indices of non-empty bins
      integer :: nlive = 0
      real(dp), allocatable :: bin_dm(:)           !! largest |D| block over a bin's pairs
   end type pair_bins_t

contains

   !
   ! Bin identity is (type, size), and TYPE CARRIES THE CONTRACTION DEGREE as
   ! well as the angular momenta, which is the shell-type assignment of [1].
   !
   ! That last part is not cosmetic.  Binning on (la, lb) alone lets a bin hold
   ! shell pairs with K = 36, 18, 6, 3 and 1 side by side, and 6-31G oxygen has
   ! all of those.  The kernel's primitive loops then run wildly different trip
   ! counts within one warp, which is the same divergence the class sort fixed
   ! at the quartet level, reappearing one loop deeper.
   !
   pure integer function bin_key(la, lb, kab, s)
      integer, intent(in) :: la, lb, kab, s
      bin_key = ((la*(LMAX + 1) + lb)*(KCAP + 1) + min(kab, KCAP))*(SMAX + 1) + s + 1
   end function bin_key

   !
   ! Build the shell-pair list and sort it into (type, size) bins.
   !
   ! O(N^2) throughout -- this is the whole point.  For 675 shells that is 228k
   ! pairs, against the 586M quartets the old path enumerated.
   !
   subroutine build_binned_pairs(nbas, sh_l, sh_np, sh_r, q, thresh, b)
      integer,  intent(in)  :: nbas
      integer,  intent(in)  :: sh_l(nbas), sh_np(nbas)
      real(dp), intent(in)  :: sh_r(3, nbas)      !! shell centres, for the intra-bin sort
      real(dp), intent(in)  :: q(:)               !! Schwarz bound, canonical pair index
      real(dp), intent(in)  :: thresh
      type(pair_bins_t), intent(out) :: b

      integer :: nb, i, j, s, k, p, n, key
      integer, allocatable :: cnt(:), pos(:), keyv(:), ti(:), tj(:)
      real(dp), allocatable :: tq(:)
      real(dp) :: qq, qmax

      nb = (LMAX + 1)*(LMAX + 1)*(KCAP + 1)*(SMAX + 1) + 1
      allocate (cnt(nb), pos(nb))
      cnt = 0

      !
      ! PRE-SCREEN before binning, guarded on significance as in [1]
      ! and which this originally lacked.  It matters because the size bucket
      ! is CLAMPED at SMAX: without a pre-screen, a pair whose true bucket is 15
      ! is stored as 9 and can then pair with a bucket-1 partner and pass a
      ! `sum <= 10` test it should have failed.  Clamping only over-admits for
      ! pairs the pre-screen has already removed.
      !
      qmax = 0.0_dp
      do i = 1, nbas
         do j = 1, i
            qmax = max(qmax, q(i*(i - 1)/2 + j))
         end do
      end do

      ! pass 1: count
      n = 0
      do i = 1, nbas
         do j = 1, i
            qq = q(i*(i - 1)/2 + j)
            if (qq*qmax <= thresh) cycle
            s = size_index(qq)
            key = bin_key(sh_l(i), sh_l(j), sh_np(i)*sh_np(j), s)
            cnt(key) = cnt(key) + 1
            n = n + 1
         end do
      end do

      b%nbin = nb
      b%npair = n
      allocate (b%sp_i(n), b%sp_j(n), b%sp_q(n))
      allocate (b%bin_off(nb), b%bin_cnt(nb), b%bin_s(nb), b%bin_la(nb), b%bin_lb(nb))
      b%bin_cnt = cnt

      b%bin_off(1) = 0
      do k = 2, nb
         b%bin_off(k) = b%bin_off(k - 1) + cnt(k - 1)
      end do
      pos = b%bin_off

      ! pass 2: scatter
      do i = 1, nbas
         do j = 1, i
            qq = q(i*(i - 1)/2 + j)
            if (qq*qmax <= thresh) cycle
            s = size_index(qq)
            key = bin_key(sh_l(i), sh_l(j), sh_np(i)*sh_np(j), s)
            pos(key) = pos(key) + 1
            b%sp_i(pos(key)) = i
            b%sp_j(pos(key)) = j
            b%sp_q(pos(key)) = qq
         end do
      end do

#ifdef TRC_SPATIAL_SORT
      call sort_bins_spatially(b, nbas, sh_r)
#endif

      ! bin metadata and the live list
      b%bin_s = -1; b%bin_la = -1; b%bin_lb = -1
      do i = 0, LMAX
         do j = 0, LMAX
            do k = 0, KCAP
               do s = 0, SMAX
                  key = bin_key(i, j, k, s)
                  b%bin_la(key) = i; b%bin_lb(key) = j; b%bin_s(key) = s
               end do
            end do
         end do
      end do

      allocate (b%live(nb))
      b%nlive = 0
      do k = 1, nb
         if (b%bin_cnt(k) > 0) then
            b%nlive = b%nlive + 1
            b%live(b%nlive) = k
         end if
      end do

      deallocate (cnt, pos)
   end subroutine build_binned_pairs

   !
   ! Order the pairs INSIDE each bin along a Morton curve through the pair
   ! centroid.
   !
   ! WHY
   ! ---
   ! Screening that removes work only pays if the removal is warp-uniform. A
   ! thread whose quartet fails the density test returns immediately, but it
   ! keeps its slot and its warp still runs until the slowest SURVIVING lane
   ! finishes. Scattered rejections therefore buy almost nothing: measured, 29%
   ! of quartets rejected at w128 bought about 4% of wall clock, and 7.5% at
   ! Gly30 bought 4%.
   !
   ! gpu4pyscf avoids this by compacting survivors into a task list before the
   ! integral kernel launches, so every thread does real work. Compaction here
   ! would cost two extra passes over ~2e9 work items, which is most of the
   ! saving.
   !
   ! This is the cheap way to the same place. Consecutive threads take
   ! consecutive ket indices within a bin, so if the ket pairs are ordered so
   ! that neighbours are spatially close, then for a fixed bra the density
   ! bound varies smoothly along the index and rejections come in runs rather
   ! than scattered. Whole warps then fail together and genuinely exit early.
   !
   ! This is a permutation within a bin. Bin identity (type, size class) is
   ! untouched, every pair still appears exactly once, and the diagonal
   ! lower-triangle indexing is order-independent -- so results are unchanged
   ! up to atomic accumulation order.
   !
   ! MEASURED, AND IT LOSES.  w128 went 5.03 -> 5.68 s (+13%), Gly30 unchanged,
   ! sum(G) bit-identical so the permutation itself is right.  The reason is
   ! that the ORIGINAL order was already the good one for memory: pairs are
   ! generated i-major, so consecutive ket indices share a shell and their
   ! primitive-pair data is largely reused across a warp.  Morton order
   ! scatters those loads, and the cache locality lost exceeds the warp
   ! clustering gained -- the survivor count is identical either way, so the
   ! clustering bought nothing measurable at all.
   !
   ! Kept behind -DTRC_SPATIAL_SORT rather than deleted: the idea is sound
   ! for a kernel whose per-thread data does not already come from a shared
   ! shell, and this records that it was tried.
   !
   subroutine sort_bins_spatially(b, nbas, sh_r)
      type(pair_bins_t), intent(inout) :: b
      integer,  intent(in) :: nbas
      real(dp), intent(in) :: sh_r(3, nbas)

      integer :: k, p, lo, n, m
      real(dp) :: lo3(3), hi3(3), span, c(3)
      integer(int64), allocatable :: code(:), kbuf(:)
      integer(int_index), allocatable :: idx(:)
      integer, allocatable :: ti(:), tj(:)
      real(dp), allocatable :: tq(:)

      if (b%npair <= 1) return

      ! bounding box over all shell centres, so the quantisation is global and
      ! bins remain comparable to one another
      lo3 = sh_r(:, 1); hi3 = sh_r(:, 1)
      do p = 2, nbas
         lo3 = min(lo3, sh_r(:, p))
         hi3 = max(hi3, sh_r(:, p))
      end do
      span = maxval(hi3 - lo3)
      if (span <= 0.0_dp) return

      allocate (code(b%npair))
      do p = 1, b%npair
         c = 0.5_dp*(sh_r(:, b%sp_i(p)) + sh_r(:, b%sp_j(p)))
         code(p) = morton3(int((c(1) - lo3(1))/span*1023.0_dp), &
                           int((c(2) - lo3(2))/span*1023.0_dp), &
                           int((c(3) - lo3(3))/span*1023.0_dp))
      end do

      allocate (ti(b%npair), tj(b%npair), tq(b%npair))
      ti = b%sp_i; tj = b%sp_j; tq = b%sp_q

      allocate (kbuf(maxval(b%bin_cnt)), idx(maxval(b%bin_cnt)))
      do k = 1, b%nbin
         n = b%bin_cnt(k)
         if (n <= 1) cycle
         lo = b%bin_off(k) + 1
         ! sort_index takes the keys 0-based and returns a 1-based permutation
         ! of them; it sorts `kbuf` in place, hence the copy.
         kbuf(1:n) = code(lo:lo + n - 1)
         call sort_index(kbuf(1:n), idx(1:n))
         do m = 1, n
            p = lo + int(idx(m)) - 1
            b%sp_i(lo + m - 1) = ti(p)
            b%sp_j(lo + m - 1) = tj(p)
            b%sp_q(lo + m - 1) = tq(p)
         end do
      end do

      deallocate (code, kbuf, idx, ti, tj, tq)
   end subroutine sort_bins_spatially

   !> Interleave the low ten bits of three coordinates.
   pure integer(int64) function morton3(ix, iy, iz)
      integer, intent(in) :: ix, iy, iz
      integer :: t, bit
      morton3 = 0_8
      do bit = 0, 9
         t = ibits(max(0, min(1023, ix)), bit, 1)
         if (t /= 0) morton3 = ibset(morton3, 3*bit)
         t = ibits(max(0, min(1023, iy)), bit, 1)
         if (t /= 0) morton3 = ibset(morton3, 3*bit + 1)
         t = ibits(max(0, min(1023, iz)), bit, 1)
         if (t /= 0) morton3 = ibset(morton3, 3*bit + 2)
      end do
   end function morton3

   !
   ! Per-bin density maximum, for the density screen.
   !
   ! The Fock contribution of a quartet is bounded by Q_ab Q_cd times the
   ! largest density block it can touch, so screening on the integral alone
   ! discards on the wrong quantity.  Taking the max over a whole bin keeps the
   ! test at bin-pair granularity, which is what keeps conditionals out of the
   ! kernel.
   !
   ! Recomputed every SCF iteration -- D changes, the bins do not.
   !
   subroutine bin_dmax(b, nbas, nao, sh_l, ao_off, dmat)
      type(pair_bins_t), intent(inout) :: b
      integer,  intent(in) :: nbas, nao
      integer,  intent(in) :: sh_l(nbas), ao_off(nbas)
      real(dp), intent(in) :: dmat(nao, nao)

      integer  :: k, p, i, j, a, c
      real(dp) :: m

      if (.not. allocated(b%bin_dm)) allocate (b%bin_dm(b%nbin))
      b%bin_dm = 0.0_dp

      do k = 1, b%nbin
         if (b%bin_cnt(k) == 0) cycle
         m = 0.0_dp
         do p = b%bin_off(k) + 1, b%bin_off(k) + b%bin_cnt(k)
            i = b%sp_i(p); j = b%sp_j(p)
            do a = 0, (sh_l(i) + 1)*(sh_l(i) + 2)/2 - 1
               do c = 0, (sh_l(j) + 1)*(sh_l(j) + 2)/2 - 1
                  m = max(m, abs(dmat(ao_off(i) + a, ao_off(j) + c)))
               end do
            end do
         end do
         b%bin_dm(k) = m
      end do
   end subroutine bin_dmax

   !
   ! Schwarz bound -> decade bucket.  int() toward zero keeps the later
   ! s_ab + s_cd test conservative.
   !
   !
   ! Schwarz bound -> decade bucket, ceiling and clamped to SMAX, matching
   ! [1].  Ceiling rather than truncation: truncation is conservative and
   ! over-admits by up to a decade in EACH factor, which measured as 37% more
   ! work than the exact canonical count at 75 waters.
   !
   pure integer function size_index(qq)
      real(dp), intent(in) :: qq
      real(dp) :: v
      if (qq >= 1.0_dp) then
         size_index = 0
      else
         v = -log10(qq)
         size_index = min(SMAX, int(ceiling(v)))
      end if
   end function size_index

end module trc_bins
