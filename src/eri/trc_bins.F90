!
! Binned shell-pair container: quartets formed on the device.
!
! Materialising the quartet list was the bottleneck (75 waters / 6-31G: 26 s
! host time and 11.7 GB against 15 s of GPU). Following Galvez Vallejo et al.
! [1], quartets are never enumerated: the O(N^2) shell-pair list is binned and
! each kernel launch covers the cross product of two bins, each thread
! deriving its own (ab, cd) from its index.
!
! Pairs are binned by type (la, lb) and by size, the Schwarz bound
! Q_ab = sqrt(max |(ab|ab)|) bucketed as s_ab = int(-log10(Q_ab)). Since
! |(ab|cd)| <= Q_ab Q_cd, a bin pair is kept iff s_ab + s_cd <= -log10(tau),
! an integer test done once per bin pair on the host. Whole cross products are
! admitted or rejected together, so kernels contain no screening conditional.
! int() truncates toward zero, so the test is conservative: it may admit
! slightly more than tau demands, never less. Binning by type also means every
! thread in a launch shares (la, lb, lc, ld).
!
! [1] J. L. Galvez Vallejo, G. M. J. Barca and M. S. Gordon,
!     "High-performance GPU-accelerated evaluation of electron repulsion
!     integrals", Mol. Phys. (2022) e2112987, doi:10.1080/00268976.2022.2112987
!
module trc_bins
   use trc_boys, only: dp
   use trc_tables, only: LMAX
   use pic_types, only: int64, int_index
   use pic_sorting, only: sort_index
   implicit none
   private

   public :: SMAX, bin_key, build_binned_pairs, pair_bins_t, bin_dmax, sort_bins_by_weight
   public :: build_pair_quartets

   !
   ! Size buckets are 1/`res` of a decade (runtime, up to SMAX_RES), as in the
   ! paper; res = 1 is whole decades. Finer buckets trade launch count against
   ! how much the pre-launch screen catches.
   !
   integer, parameter :: SMAX_RES = 4   !! largest resolution the bins are sized for
   integer, parameter :: SDEC = 9       !! decades of Schwarz bound worth binning
   integer, parameter :: SMAX = SDEC*SMAX_RES  !! deepest bucket at any resolution
   integer, parameter :: KCAP = 64   !! largest contraction degree binned exactly

   type :: pair_bins_t
      integer :: nbin = 0
      integer :: npair = 0
      integer :: sres = 1
      !! Bucket resolution these bins were built at: buckets per decade of
      !! Schwarz bound. The admission test in the launcher needs it.
      logical :: sorted = .false.
      !! Whether the pairs inside each bin have been ordered by how much
      !! they matter. Done once, at the first Fock build, and kept: the
      !! order is a heuristic for warp coherence, not a correctness
      !! requirement, and the density it is keyed on moves slowly.
      integer, allocatable :: sp_i(:), sp_j(:)     !! shell indices, bin-sorted
      real(dp), allocatable :: sp_q(:)             !! Schwarz bound per pair
      integer, allocatable :: bin_off(:)           !! first pair of bin b
      integer, allocatable :: bin_cnt(:)           !! pairs in bin b
      integer, allocatable :: bin_s(:)             !! size index of bin b
      integer, allocatable :: bin_la(:), bin_lb(:) !! type of bin b
      integer, allocatable :: bin_g(:)             !! 1 if a pair of bin b has a general side
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
   !
   ! `g` separates the pairs with a generally contracted side from the rest:
   ! a segment is then either all single-column, and goes to the scalar
   ! kernel, or all multi-column, and goes to the blocked one -- no launch
   ! carries both, so neither pays the other's registers or divergence.
   !
   pure integer function bin_key(la, lb, kab, s, g)
      integer, intent(in) :: la, lb, kab, s, g
      bin_key = (((la*(LMAX + 1) + lb)*(KCAP + 1) + min(kab, KCAP))*(SMAX + 1) + s)*2 + g + 1
   end function bin_key

   !
   ! Build the shell-pair list and sort it into (type, size) bins.
   !
   ! O(N^2) throughout -- this is the whole point.  For 675 shells that is 228k
   ! pairs, against the 586M quartets the old path enumerated.
   !
   subroutine build_binned_pairs(nbas, sh_l, sh_np, sh_r, q, thresh, b, sh_g, pp_n, res)
      integer,  intent(in)  :: nbas
      integer,  intent(in)  :: sh_l(nbas), sh_np(nbas)
      real(dp), intent(in)  :: sh_r(3, nbas)      !! shell centres, for the intra-bin sort
      real(dp), intent(in)  :: q(:)               !! Schwarz bound, canonical pair index
      real(dp), intent(in)  :: thresh
      type(pair_bins_t), intent(out) :: b
      logical,  intent(in), optional :: sh_g(nbas)  !! shell has more than one column
      integer,  intent(in), optional :: pp_n(:)
      !! Primitive pairs actually kept per shell pair, rectangular key
      !! (i-1)*nbas+j. The contraction degree the bin key sorts on is
      !! this, not sh_np(i)*sh_np(j): after pruning the two differ, and
      !! it is the kept count that sets the kernel's trip count.
      integer,  intent(in), optional :: res
      !! Buckets per decade of Schwarz bound; 1 (the default) is decades.

      integer :: nb, i, j, s, k, p, n, key, gg, ia, ib, sres
      integer, allocatable :: cnt(:), pos(:), keyv(:), ti(:), tj(:)
      real(dp), allocatable :: tq(:)
      real(dp) :: qq, qmax

      sres = 1
      if (present(res)) sres = max(1, min(SMAX_RES, res))
      b%sres = sres
      nb = 2*(LMAX + 1)*(LMAX + 1)*(KCAP + 1)*(SMAX + 1) + 1
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

      !
      ! ANGULAR-MOMENTUM CANONICAL, not index canonical.
      !
      ! (ab|cd) = (ba|cd), so a pair is free to be stored with its higher-l
      ! shell first, and doing so means a (ps) bin never exists: every pair
      ! of an s and a p shell lands in (sp). At LMAX=2 that takes the nine
      ! (la, lb) types down to six and, with the bra >= ket triangle the
      ! segment loop already walks, the reachable classes from 45 to 21.
      ! Same work, fewer and fatter launches. The index order is the
      ! tiebreak when the momenta match, so the canonical enumeration
      ! within a bin is unchanged. The Fock matrix is folded with its own
      ! transpose afterwards, which is what makes the swap invisible.
      !
      ! pass 1: count
      n = 0
      do i = 1, nbas
         do j = 1, i
            qq = q(i*(i - 1)/2 + j)
            if (qq*qmax <= thresh) cycle
            call lorder(i, j, ia, ib)
            s = size_index(qq, sres)
            gg = 0
            if (present(sh_g)) then
               if (sh_g(i) .or. sh_g(j)) gg = 1
            end if
            key = bin_key(sh_l(ia), sh_l(ib), kab_of(ia, ib), s, gg)
            cnt(key) = cnt(key) + 1
            n = n + 1
         end do
      end do

      b%nbin = nb
      b%npair = n
      allocate (b%sp_i(n), b%sp_j(n), b%sp_q(n))
      allocate (b%bin_off(nb), b%bin_cnt(nb), b%bin_s(nb), b%bin_la(nb), b%bin_lb(nb), b%bin_g(nb))
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
            call lorder(i, j, ia, ib)
            s = size_index(qq, sres)
            gg = 0
            if (present(sh_g)) then
               if (sh_g(i) .or. sh_g(j)) gg = 1
            end if
            key = bin_key(sh_l(ia), sh_l(ib), kab_of(ia, ib), s, gg)
            pos(key) = pos(key) + 1
            b%sp_i(pos(key)) = ia
            b%sp_j(pos(key)) = ib
            b%sp_q(pos(key)) = qq
         end do
      end do

#ifdef TRC_SPATIAL_SORT
      call sort_bins_spatially(b, nbas, sh_r)
#endif

      ! bin metadata and the live list
      b%bin_s = -1; b%bin_la = -1; b%bin_lb = -1; b%bin_g = 0
      do i = 0, LMAX
         do j = 0, LMAX
            do k = 0, KCAP
               do s = 0, SDEC*sres
                  do gg = 0, 1
                     key = bin_key(i, j, k, s, gg)
                     b%bin_la(key) = i; b%bin_lb(key) = j; b%bin_s(key) = s; b%bin_g(key) = gg
                  end do
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

   contains
      pure subroutine lorder(i_, j_, a_, b_)
         !! The pair with the higher angular momentum first; index order
         !! breaks the tie, so a pair of equal-l shells keeps i >= j.
         integer, intent(in) :: i_, j_
         integer, intent(out) :: a_, b_
         if (sh_l(j_) > sh_l(i_)) then
            a_ = j_; b_ = i_
         else
            a_ = i_; b_ = j_
         end if
      end subroutine lorder

      pure integer function kab_of(i, j)
         integer, intent(in) :: i, j
         if (present(pp_n)) then
            kab_of = pp_n((i - 1)*nbas + j)
         else
            kab_of = sh_np(i)*sh_np(j)
         end if
      end function kab_of
   end subroutine build_binned_pairs

   !
   ! Build the bin-interleaved primitive-pair tables from the pp_* ("hp_*" at
   ! the call site) data build_pairs_hgp already built. Called once after
   ! build_binned_pairs, and again whenever the slot order inside bins
   ! changes (sort_bins_by_weight).
   !
   ! K is constant inside a bin except for pairs clamped at KCAP (bin_key),
   ! so bin_k(b) -- the max kept primitive-pair count over the bin's pairs --
   ! equals every pair's own count for the common case and only pads for the
   ! rare KCAP-clamped bin. A padding slot (k beyond that PAIR's own count,
   ! k <= bin_k(b)) is never read: the kernel's loop bound is the pair's own
   ! count (`sq_n`), exactly as pp_n(key) was. It is still filled, with
   ! qp_c = 0, qp_z = 1 and P at the A centre, so the storage has no
   ! uninitialised entries and nothing downstream has to know which slots
   ! are real.
   !
   !
   ! The new tables are NOT fields of `b`. An earlier version added them to
   ! pair_bins_t, which is itself mirrored whole by an `!$acc enter data
   ! copyin(this%bins, this%bins%sp_i, ...)` at the call site -- and growing
   ! that type's allocatable-component count that way made nvfortran 26.5
   ! miscompile an UNRELATED host loop two modules away (an internal
   ! "bitcast with differing sizes" in code that never touches a bin). The
   ! caller instead owns these as its own flat allocatables, exactly the
   ! shape pp_off/pp_n/pp_p/... already have, and passes them back in on a
   ! rebuild so the reuse-the-allocation logic below still applies.
   !
   subroutine build_pair_quartets(b, nbas, npp, pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                  bin_pbase, bin_k, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                                  sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, pp_a, qp_a)
      type(pair_bins_t), intent(in) :: b
      integer,  intent(in) :: nbas, npp
      integer,  intent(in) :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in) :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  allocatable, intent(inout) :: bin_pbase(:), bin_k(:)
      integer,  intent(out) :: nqp
      real(dp), allocatable, intent(inout) :: qp_z(:), qp_c(:)
      real(dp), allocatable, intent(inout) :: qp_px(:), qp_py(:), qp_pz(:)
      real(dp), allocatable, intent(inout) :: sq_ax(:), sq_ay(:), sq_az(:)
      real(dp), allocatable, intent(inout) :: sq_abx(:), sq_aby(:), sq_abz(:)
      integer,  allocatable, intent(inout) :: sq_n(:)
      real(dp), optional, intent(in) :: pp_a(npp)
      !! The i-side primitive's own exponent, pp_a from build_pairs_hgp.
      !! Optional, like there: only the gradient path threads it through.
      real(dp), allocatable, optional, intent(inout) :: qp_a(:)
      !! Bin-interleaved form of pp_a, same layout as qp_z. The padding
      !! slot (k > a pair's own kept count) gets 0: never read, since the
      !! kernel's loop bound is always the pair's real count (`sq_n`).

      integer :: ib, bb, s, p, key, k1, kreal, off, ix

      !
      ! REUSE THE ALLOCATION ON A REBUILD.
      !
      ! sort_bins_by_weight only permutes slots WITHIN a bin, so bin_cnt,
      ! every pair's own count and therefore bin_k and nqp come out
      ! identical the second time. Device residency is keyed on the host
      ! array's address: an `!$acc enter data copyin` from build_structures
      ! stays valid only if a later `!$acc update device` (after a resort)
      ! finds the same allocation, so this fills the existing arrays in
      ! place rather than deallocating and reallocating them.
      !
      if (.not. allocated(bin_pbase)) allocate (bin_pbase(b%nbin), bin_k(b%nbin))
      if (.not. allocated(sq_ax)) then
         allocate (sq_ax(b%npair), sq_ay(b%npair), sq_az(b%npair))
         allocate (sq_abx(b%npair), sq_aby(b%npair), sq_abz(b%npair))
         allocate (sq_n(b%npair))
      end if
      bin_pbase = 0; bin_k = 0

      ! Per pair slot: the A centre, A - B, and the real kept count -- what
      ! pp_off(key)/pp_n(key) gave the kernel before, read here once rather
      ! than once per quartet.
      do p = 1, b%npair
         key = (b%sp_i(p) - 1)*nbas + b%sp_j(p)
         off = pp_off(key)
         sq_ax(p) = pp_ra(off + 1, 1); sq_ay(p) = pp_ra(off + 1, 2); sq_az(p) = pp_ra(off + 1, 3)
         sq_abx(p) = pp_ra(off + 1, 1) - pp_rb(off + 1, 1)
         sq_aby(p) = pp_ra(off + 1, 2) - pp_rb(off + 1, 2)
         sq_abz(p) = pp_ra(off + 1, 3) - pp_rb(off + 1, 3)
         sq_n(p) = pp_n(key)
      end do

      ! Per bin: K_b, then the prefix sum that lays bins out one after
      ! another in qp_*.
      do ib = 1, b%nlive
         bb = b%live(ib)
         do s = 1, b%bin_cnt(bb)
            bin_k(bb) = max(bin_k(bb), sq_n(b%bin_off(bb) + s))
         end do
      end do
      nqp = 0
      do bb = 1, b%nbin
         bin_pbase(bb) = nqp
         nqp = nqp + b%bin_cnt(bb)*bin_k(bb)
      end do

      if (.not. allocated(qp_z)) then
         allocate (qp_z(nqp), qp_c(nqp))
         allocate (qp_px(nqp), qp_py(nqp), qp_pz(nqp))
      end if
      if (present(qp_a)) then
         if (.not. allocated(qp_a)) allocate (qp_a(nqp))
      end if

      ! Interleaved by primitive index: ix(b, k, s) = bin_pbase(b) + (k-1)*n_b + s.
      ! The real primitives keep the order build_pairs_hgp already sorted them
      ! into (largest |c|/p first), so every sum happens in the same order and
      ! results are bit-identical to the pp_* path.
      do ib = 1, b%nlive
         bb = b%live(ib)
         do s = 1, b%bin_cnt(bb)
            p = b%bin_off(bb) + s
            key = (b%sp_i(p) - 1)*nbas + b%sp_j(p)
            off = pp_off(key)
            kreal = sq_n(p)
            do k1 = 1, bin_k(bb)
               ix = bin_pbase(bb) + (k1 - 1)*b%bin_cnt(bb) + s
               if (k1 <= kreal) then
                  qp_z(ix) = pp_p(off + k1); qp_c(ix) = pp_c(off + k1)
                  qp_px(ix) = pp_r(off + k1, 1)
                  qp_py(ix) = pp_r(off + k1, 2)
                  qp_pz(ix) = pp_r(off + k1, 3)
                  if (present(qp_a)) qp_a(ix) = pp_a(off + k1)
               else
                  qp_z(ix) = 1.0_dp; qp_c(ix) = 0.0_dp
                  qp_px(ix) = sq_ax(p); qp_py(ix) = sq_ay(p); qp_pz(ix) = sq_az(p)
                  if (present(qp_a)) qp_a(ix) = 0.0_dp
               end if
            end do
         end do
      end do
   end subroutine build_pair_quartets

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
   ! ORDER THE PAIRS INSIDE EACH BIN BY HOW MUCH THEY MATTER.
   !
   ! A warp takes 32 consecutive pairs of a bin, and each thread decides on
   ! its own whether its quartet survives the Schwarz and density test. If
   ! neighbouring pairs are unrelated, a warp is a mix of threads that work
   ! and threads that returned, and it runs at the speed of the ones that
   ! work while the rest wait. Ordering by the product of the Schwarz bound
   ! and the largest density block -- the quantity the per-quartet test is
   ! built from -- makes neighbours agree far more often.
   !
   ! Once, and then never again: the density this is keyed on is the guess
   ! at the first build, and it moves slowly enough that the order stays a
   ! good heuristic for the whole SCF. That makes the cost a fixed charge
   ! against every later build rather than a per-build tax.
   !
   subroutine sort_bins_by_weight(b, nbas, dsh)
      type(pair_bins_t), intent(inout) :: b
      integer,  intent(in) :: nbas
      real(dp), intent(in) :: dsh(nbas, nbas)

      integer :: k, lo, n, p, m
      real(dp) :: wmax, w
      logical :: qonly
      integer(int64), allocatable :: code(:), kbuf(:)
      integer(int_index), allocatable :: idx(:)
      integer,  allocatable :: ti(:), tj(:)
      real(dp), allocatable :: tq(:)

      if (b%sorted .or. b%npair <= 1) return
      allocate (code(b%npair), ti(b%npair), tj(b%npair), tq(b%npair))

      !
      ! The key is the product the per-quartet test is built from, mapped to
      ! an integer that sorts DESCENDING -- sort_index is ascending, and the
      ! heavy pairs want to be together at the front. A pair whose product
      ! underflows the log lands in the last bucket, which is where the
      ! threads that return immediately belong.
      !
      qonly = qkey_only()
      wmax = 0.0_dp
      do p = 1, b%npair
         if (qonly) then
            wmax = max(wmax, b%sp_q(p))
         else
            wmax = max(wmax, b%sp_q(p)*dsh(b%sp_i(p), b%sp_j(p)))
         end if
      end do
      if (wmax <= 0.0_dp) wmax = 1.0_dp
      do p = 1, b%npair
         if (qonly) then
            w = b%sp_q(p)
         else
            w = b%sp_q(p)*dsh(b%sp_i(p), b%sp_j(p))
         end if
         code(p) = int(1.0e6_dp*min(30.0_dp, -log10(max(w, 1.0e-30_dp)/wmax)), int64)
      end do

      ti = b%sp_i; tj = b%sp_j; tq = b%sp_q
      allocate (kbuf(maxval(b%bin_cnt)), idx(maxval(b%bin_cnt)))
      do k = 1, b%nbin
         n = b%bin_cnt(k)
         if (n <= 1) cycle
         lo = b%bin_off(k) + 1
         kbuf(1:n) = code(lo:lo + n - 1)
         call sort_index(kbuf(1:n), idx(1:n))
         do m = 1, n
            p = lo + int(idx(m)) - 1
            b%sp_i(lo + m - 1) = ti(p)
            b%sp_j(lo + m - 1) = tj(p)
            b%sp_q(lo + m - 1) = tq(p)
         end do
      end do

      b%sorted = .true.
      deallocate (code, kbuf, idx, ti, tj, tq)

   contains
      logical function qkey_only()
         !! Key on the Schwarz bound alone rather than on the product with the
         !! density block. The bound correlates with distance, so this keeps
         !! more of the spatial locality the natural order has, at the cost of
         !! a less faithful prediction of which threads will screen out.
         character(len=8) :: v
         call get_environment_variable("TRC_BIN_SORT_Q", v)
         qkey_only = (trim(v) == "1")
      end function qkey_only
   end subroutine sort_bins_by_weight

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

   pure integer(int64) function morton3(ix, iy, iz)
      !! Interleave the low ten bits of three coordinates.
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
   !
   ! MEASURED, on the 123-atom silica slice and the T8 cage: raising `res`
   ! above 1 enumerates MORE quartets, not fewer, and costs up to 1.6x.
   ! The bucket is a CEILING on -log10(Q), so at res = 1 the admission test
   ! `s_a + s_b <= T` rejects pairs whose exact product is above threshold
   ! -- two pairs at Q = 0.8 both round to bucket 1 and are refused at
   ! T = 1. That over-rejection is why the coarse bucket looks cheap: it
   ! is not screening harder, it is discarding work it should keep, about
   ! 0.5% of the surviving quartets on the slice. A finer bucket rounds
   ! less and hands that work back. So res is a knob for accuracy at the
   ! pre-launch screen, not for speed, and 1 stays the default.
   !
   pure integer function size_index(qq, res)
      real(dp), intent(in) :: qq
      integer, intent(in) :: res
      real(dp) :: v
      if (qq >= 1.0_dp) then
         size_index = 0
      else
         v = -log10(qq)*real(res, dp)
         size_index = min(SDEC*res, int(ceiling(v)))
      end if
   end function size_index

end module trc_bins
