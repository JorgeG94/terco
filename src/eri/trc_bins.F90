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
   use pic_types, only: dp
   use trc_tables, only: LMAX
   use pic_types, only: int64, int_index
   use pic_sorting, only: sort_index
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t
   implicit none
   private

   public :: SMAX, bins_build_binned_pairs, pair_bins_t, bins_sort_by_weight
   public :: bins_build_pair_quartets
   ! bin_key/bin_dmax stay private: nothing outside this module uses them.

   !
   ! Size buckets are 1/`res` of a decade (runtime, up to SMAX_RES), as in the
   ! paper; res = 1 is whole decades. Finer buckets trade launch count against
   ! how much the pre-launch screen catches.
   !
   integer, parameter :: SMAX_RES = 4   !! largest resolution the bins are sized for
   integer, parameter :: SDEC = 9       !! decades of Schwarz bound worth binning
   integer, parameter :: SMAX = SDEC*SMAX_RES  !! deepest bucket at any resolution
   integer, parameter :: KCAP = 64   !! largest contraction degree binned exactly
   integer, parameter :: BIN_CAP = 32768
   !! Most pairs one bin may hold. A segment is the cross product of two
   !! bins, and the kernels index within a segment by a default integer
   !! (each Fock kernel's within-segment `t`), so a segment must stay below 2**31 quartets:
   !! 32768**2 is 2**30. A bin past this is split by split_large_bins.
   !! TRC_BIN_CAP in the environment lowers it, for testing the split on a
   !! system small enough for CI (terco/check_binfock_split).

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
      integer, allocatable :: bin_la(:), bin_lb(:)  !! type of bin b
      integer, allocatable :: bin_g(:)             !! 1 if a pair of bin b has a general side
      integer, allocatable :: live(:)              !! indices of non-empty bins
      integer :: nlive = 0
      real(dp), allocatable :: bin_dm(:)           !! largest |D| block over a bin's pairs
   end type pair_bins_t

   !
   ! The launch plan: which bin pairs survive the admission test, how their
   ! quartets are numbered, and in what order the class-specialised kernels
   ! cover them. One segment is one admitted bin pair; a thread finds its own
   ! segment by searching the prefix sum `sOff`. binkernel_fock_bins and
   ! grad_bins build the same plan, so it is built here once for both.
   !
   ! The arrays keep the names the generated kernels give their dummies
   ! (sA, sNB, ...). Never mapped whole: a launcher enters the arrays it
   ! hands to its kernels one by one.
   !
   type, public :: trc_bin_plan_t
      integer :: nseg = 0
      !! Segments: admitted bin pairs.
      integer(int64) :: nwork = 0
      !! Work items over all segments: sOff(nseg + 1).
      integer, allocatable :: sA(:), sNB(:)
      !! (nseg) pairs in the bra bin and in the ket bin.
      integer, allocatable :: sOA(:), sOB(:)
      !! (nseg) first pair - 1 of the bra and ket bins in pair_bins_t's lists.
      integer, allocatable :: sPA(:), sPB(:)
      !! (nseg) first slot - 1 of the two bins in the bin-interleaved tables.
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      !! (nseg) the class (la, lb | lc, ld).
      logical, allocatable :: sD(:)
      !! (nseg) bra bin == ket bin: the segment is a triangle, not a rectangle.
      integer(int64), allocatable :: sOff(:)
      !! (nseg + 1) prefix sum of the work items; sOff(1) = 0.
   contains
      procedure :: build         => bins_plan_build
      procedure :: sort_by_class => bins_plan_sort_by_class
      procedure :: class_key     => bins_plan_class_key
      procedure :: release       => bins_plan_release
   end type trc_bin_plan_t

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
   pure integer function bin_key(la, lb, kab, s, g) result(res)
      integer, intent(in) :: la, lb, kab, s, g
      res = (((la*(LMAX + 1) + lb)*(KCAP + 1) + min(kab, KCAP))*(SMAX + 1) + s)*2 + g + 1
   end function bin_key

   !
   ! Build the shell-pair list and sort it into (type, size) bins.
   !
   ! O(N^2) throughout -- this is the whole point.  For 675 shells that is 228k
   ! pairs, against the 586M quartets the old path enumerated.
   !
   subroutine bins_build_binned_pairs(sh_l, sh_np, sh_r, q, thresh, b, pp_n, res)
      integer,  intent(in)  :: sh_l(:), sh_np(:)  !! (nbas) per shell; their size is the shell count
      real(dp), intent(in)  :: sh_r(:, :)         !! (3, nbas) shell centres, for the intra-bin sort
      real(dp), intent(in)  :: q(:)               !! Schwarz bound, canonical pair index
      real(dp), intent(in)  :: thresh
      type(pair_bins_t), intent(out) :: b
      integer,  intent(in), optional :: pp_n(:)
      !! Primitive pairs actually kept per shell pair, rectangular key
      !! (i-1)*nbas+j. The contraction degree the bin key sorts on is
      !! this, not sh_np(i)*sh_np(j): after pruning the two differ, and
      !! it is the kept count that sets the kernel's trip count.
      integer,  intent(in), optional :: res
      !! Buckets per decade of Schwarz bound; 1 (the default) is decades.

      integer :: nbas, nb, i, j, s, k, p, n, key, gg, ia, ib, sres
      integer, allocatable :: cnt(:), pos(:), keyv(:), ti(:), tj(:)
      real(dp), allocatable :: tq(:)
      real(dp) :: qq, qmax

      nbas = size(sh_l)
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
      b%bin_s = -1
      b%bin_la = -1
      b%bin_lb = -1
      b%bin_g = 0
      do i = 0, LMAX
         do j = 0, LMAX
            do k = 0, KCAP
               do s = 0, SDEC*sres
                  do gg = 0, 1
                     key = bin_key(i, j, k, s, gg)
                     b%bin_la(key) = i
                     b%bin_lb(key) = j
                     b%bin_s(key) = s
                     b%bin_g(key) = gg
                  end do
               end do
            end do
         end do
      end do

      call split_large_bins(b)

      allocate (b%live(b%nbin))
      b%nlive = 0
      do k = 1, b%nbin
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
            a_ = j_
            b_ = i_
         else
            a_ = i_
            b_ = j_
         end if
      end subroutine lorder

      pure integer function kab_of(i, j) result(res)
         integer, intent(in) :: i, j
         if (present(pp_n)) then
            res = pp_n((i - 1)*nbas + j)
         else
            res = sh_np(i)*sh_np(j)
         end if
      end function kab_of
   end subroutine bins_build_binned_pairs

   !
   ! Split every bin of more than BIN_CAP pairs into consecutive chunks of at
   ! most BIN_CAP, appended as new bins with the same (type, size, g).
   !
   ! Widening the kernels' within-segment index instead would cost registers
   ! in the hottest loop, which is the binding constraint there. Splitting
   ! costs nothing in the kernels: the quartets are the same ones. A split
   ! bin's own triangle becomes each chunk's triangle plus the rectangle
   ! between every two chunks, which the launcher's bra >= ket walk over
   ! bin pairs already produces, each quartet once. Everything per bin --
   ! primitive layout, K, density bound, class -- is derived later from
   ! bin_off/bin_cnt and the metadata copied here, never from the bin's
   ! index. Bins at or below the cap are untouched, so a system that never
   ! reaches it gets bit-identical results.
   !
   ! First reached around 14,000 basis functions (a 9514-shell system,
   ! where one bin outgrew 65,536 pairs and its triangle 2**31).
   !
   subroutine split_large_bins(b)
      type(pair_bins_t), intent(inout) :: b
      integer :: k, c, nextra, nb, ich, nch, cap
      integer, allocatable :: off(:), cnt(:), sz(:), la(:), lb(:), g(:)

      cap = BIN_CAP
      block
         character(len=16) :: env
         integer :: st, v, ios
         call get_environment_variable("TRC_BIN_CAP", env, status=st)
         if (st == 0) then
            read (env, *, iostat=ios) v
            if (ios == 0 .and. v >= 1) cap = min(cap, v)
         end if
      end block

      nextra = 0
      do k = 1, b%nbin
         if (b%bin_cnt(k) > cap) nextra = nextra + (b%bin_cnt(k) - 1)/cap
      end do
      if (nextra == 0) return

      nb = b%nbin + nextra
      allocate (off(nb), cnt(nb), sz(nb), la(nb), lb(nb), g(nb))
      off(1:b%nbin) = b%bin_off
      cnt(1:b%nbin) = b%bin_cnt
      sz(1:b%nbin) = b%bin_s
      la(1:b%nbin) = b%bin_la
      lb(1:b%nbin) = b%bin_lb
      g(1:b%nbin) = b%bin_g
      c = b%nbin
      do k = 1, b%nbin
         if (b%bin_cnt(k) <= cap) cycle
         nch = (b%bin_cnt(k) + cap - 1)/cap
         cnt(k) = cap
         do ich = 2, nch
            c = c + 1
            off(c) = b%bin_off(k) + (ich - 1)*cap
            cnt(c) = min(cap, b%bin_cnt(k) - (ich - 1)*cap)
            sz(c) = b%bin_s(k)
            la(c) = b%bin_la(k)
            lb(c) = b%bin_lb(k)
            g(c) = b%bin_g(k)
         end do
      end do
      b%nbin = nb
      call move_alloc(off, b%bin_off)
      call move_alloc(cnt, b%bin_cnt)
      call move_alloc(sz, b%bin_s)
      call move_alloc(la, b%bin_la)
      call move_alloc(lb, b%bin_lb)
      call move_alloc(g, b%bin_g)
   end subroutine split_large_bins

   !
   ! Build the bin-interleaved primitive-pair tables from the primitive
   ! pairs hgp_build_pairs already built. Called once after
   ! bins_build_binned_pairs, and again whenever the slot order inside bins
   ! changes (bins_sort_by_weight).
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
   ! "bitcast with differing sizes" in code that never touches a bin). They
   ! live in their own trc_pair_quartets_t instead, which is never mapped
   ! whole -- its arrays go to the device one by one -- and the caller
   ! passes it back in on a rebuild so the reuse-the-allocation logic below
   ! still applies.
   !
   subroutine bins_build_pair_quartets(b, pairs, quartets)
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      !! The primitive pairs from hgp_build_pairs.
      type(trc_pair_quartets_t), intent(inout) :: quartets
      !! Filled in place: `intent(inout)` because a rebuild after a resort
      !! reuses the allocation (see below). A padding slot (k past a pair's
      !! own kept count) is never read, since the kernel's loop bound is
      !! always the pair's real count (`n`); its `a` is 0.

      integer :: nbas, nqp, ib, bb, s, p, key, k1, kreal, off, ix

      !
      ! REUSE THE ALLOCATION ON A REBUILD.
      !
      ! bins_sort_by_weight only permutes slots WITHIN a bin, so bin_cnt,
      ! every pair's own count and therefore bin_k and nqp come out
      ! identical the second time. Device residency is keyed on the host
      ! array's address: an `!$acc enter data copyin` from build_structures
      ! stays valid only if a later `!$acc update device` (after a resort)
      ! finds the same allocation, so this fills the existing arrays in
      ! place rather than deallocating and reallocating them.
      !
      nbas = pairs%nbas
      if (.not. allocated(quartets%bin_pbase)) allocate (quartets%bin_pbase(b%nbin), quartets%bin_k(b%nbin))
      if (.not. allocated(quartets%ax)) then
         allocate (quartets%ax(b%npair), quartets%ay(b%npair), quartets%az(b%npair))
         allocate (quartets%abx(b%npair), quartets%aby(b%npair), quartets%abz(b%npair))
         allocate (quartets%n(b%npair))
      end if
      quartets%bin_pbase = 0
      quartets%bin_k = 0

      ! Per pair slot: the A centre, A - B, and the real kept count -- what
      ! pairs%off(key)/pairs%n(key) would give the kernel, read here once
      ! rather than once per quartet.
      do p = 1, b%npair
         key = (b%sp_i(p) - 1)*nbas + b%sp_j(p)
         off = pairs%off(key)
         quartets%ax(p) = pairs%ra(off + 1, 1)
         quartets%ay(p) = pairs%ra(off + 1, 2)
         quartets%az(p) = pairs%ra(off + 1, 3)
         quartets%abx(p) = pairs%ra(off + 1, 1) - pairs%rb(off + 1, 1)
         quartets%aby(p) = pairs%ra(off + 1, 2) - pairs%rb(off + 1, 2)
         quartets%abz(p) = pairs%ra(off + 1, 3) - pairs%rb(off + 1, 3)
         quartets%n(p) = pairs%n(key)
      end do

      ! Per bin: K_b, then the prefix sum that lays bins out one after
      ! another in the slot arrays.
      do ib = 1, b%nlive
         bb = b%live(ib)
         do s = 1, b%bin_cnt(bb)
            quartets%bin_k(bb) = max(quartets%bin_k(bb), quartets%n(b%bin_off(bb) + s))
         end do
      end do
      nqp = 0
      do bb = 1, b%nbin
         quartets%bin_pbase(bb) = nqp
         nqp = nqp + b%bin_cnt(bb)*quartets%bin_k(bb)
      end do

      if (.not. allocated(quartets%z)) then
         allocate (quartets%z(nqp), quartets%c(nqp), quartets%a(nqp))
         allocate (quartets%px(nqp), quartets%py(nqp), quartets%pz(nqp))
      end if
      quartets%nqp = nqp

      ! Interleaved by primitive index: ix(b, k, s) = bin_pbase(b) + (k-1)*n_b + s.
      ! The real primitives keep the order hgp_build_pairs already sorted them
      ! into (largest |c|/p first), so every sum happens in the same order and
      ! results are bit-identical to the pp_* path.
      do ib = 1, b%nlive
         bb = b%live(ib)
         do s = 1, b%bin_cnt(bb)
            p = b%bin_off(bb) + s
            key = (b%sp_i(p) - 1)*nbas + b%sp_j(p)
            off = pairs%off(key)
            kreal = quartets%n(p)
            do k1 = 1, quartets%bin_k(bb)
               ix = quartets%bin_pbase(bb) + (k1 - 1)*b%bin_cnt(bb) + s
               if (k1 <= kreal) then
                  quartets%z(ix) = pairs%p(off + k1)
                  quartets%c(ix) = pairs%c(off + k1)
                  quartets%px(ix) = pairs%r(off + k1, 1)
                  quartets%py(ix) = pairs%r(off + k1, 2)
                  quartets%pz(ix) = pairs%r(off + k1, 3)
                  quartets%a(ix) = pairs%a(off + k1)
               else
                  quartets%z(ix) = 1.0_dp
                  quartets%c(ix) = 0.0_dp
                  quartets%px(ix) = quartets%ax(p)
                  quartets%py(ix) = quartets%ay(p)
                  quartets%pz(ix) = quartets%az(p)
                  quartets%a(ix) = 0.0_dp
               end if
            end do
         end do
      end do
   end subroutine bins_build_pair_quartets

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
   subroutine bins_sort_by_weight(b, nbas, dsh)
      type(pair_bins_t), intent(inout) :: b
      integer,  intent(in) :: nbas
      real(dp), intent(in) :: dsh(nbas, nbas)

      integer :: k, lo, n, p, m
      real(dp) :: wmax, w
      logical :: qonly
      integer(int64), allocatable :: code(:), kbuf(:)
      integer(int_index), allocatable :: index(:)
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

      ti = b%sp_i
      tj = b%sp_j
      tq = b%sp_q
      allocate (kbuf(maxval(b%bin_cnt)), index(maxval(b%bin_cnt)))
      do k = 1, b%nbin
         n = b%bin_cnt(k)
         if (n <= 1) cycle
         lo = b%bin_off(k) + 1
         kbuf(1:n) = code(lo:lo + n - 1)
         call sort_index(kbuf(1:n), index(1:n))
         do m = 1, n
            p = lo + int(index(m)) - 1
            b%sp_i(lo + m - 1) = ti(p)
            b%sp_j(lo + m - 1) = tj(p)
            b%sp_q(lo + m - 1) = tq(p)
         end do
      end do

      b%sorted = .true.
      deallocate (code, kbuf, index, ti, tj, tq)

   contains
      logical function qkey_only() result(res)
         !! Key on the Schwarz bound alone rather than on the product with the
         !! density block. The bound correlates with distance, so this keeps
         !! more of the spatial locality the natural order has, at the cost of
         !! a less faithful prediction of which threads will screen out.
         character(len=8) :: v
         call get_environment_variable("TRC_BIN_SORT_Q", v)
         res = (trim(v) == "1")
      end function qkey_only
   end subroutine bins_sort_by_weight

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
      real(dp) :: span
      real(dp) :: lo3(3), hi3(3), c(3)
      integer(int64), allocatable :: code(:), kbuf(:)
      integer(int_index), allocatable :: index(:)
      integer, allocatable :: ti(:), tj(:)
      real(dp), allocatable :: tq(:)

      if (b%npair <= 1) return

      ! bounding box over all shell centres, so the quantisation is global and
      ! bins remain comparable to one another
      lo3 = sh_r(:, 1)
      hi3 = sh_r(:, 1)
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
      ti = b%sp_i
      tj = b%sp_j
      tq = b%sp_q

      allocate (kbuf(maxval(b%bin_cnt)), index(maxval(b%bin_cnt)))
      do k = 1, b%nbin
         n = b%bin_cnt(k)
         if (n <= 1) cycle
         lo = b%bin_off(k) + 1
         ! sort_index takes the keys 0-based and returns a 1-based permutation
         ! of them; it sorts `kbuf` in place, hence the copy.
         kbuf(1:n) = code(lo:lo + n - 1)
         call sort_index(kbuf(1:n), index(1:n))
         do m = 1, n
            p = lo + int(index(m)) - 1
            b%sp_i(lo + m - 1) = ti(p)
            b%sp_j(lo + m - 1) = tj(p)
            b%sp_q(lo + m - 1) = tq(p)
         end do
      end do

      deallocate (code, kbuf, index, ti, tj, tq)
   end subroutine sort_bins_spatially

   pure integer(int64) function morton3(ix, iy, iz) result(res)
      !! Interleave the low ten bits of three coordinates.
      integer, intent(in) :: ix, iy, iz
      integer :: t, bit
      res = 0_int64
      do bit = 0, 9
         t = ibits(max(0, min(1023, ix)), bit, 1)
         if (t /= 0) res = ibset(res, 3*bit)
         t = ibits(max(0, min(1023, iy)), bit, 1)
         if (t /= 0) res = ibset(res, 3*bit + 1)
         t = ibits(max(0, min(1023, iz)), bit, 1)
         if (t /= 0) res = ibset(res, 3*bit + 2)
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
            i = b%sp_i(p)
            j = b%sp_j(p)
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
   pure integer function size_index(qq, res) result(res_)
      real(dp), intent(in) :: qq
      integer, intent(in) :: res
      real(dp) :: v
      if (qq >= 1.0_dp) then
         res_ = 0
      else
         v = -log10(qq)*real(res, dp)
         res_ = min(SDEC*res, int(ceiling(v)))
      end if
   end function size_index

   !
   ! Admit bin pairs and number their quartets.
   !
   ! ONE launch over every admitted bin pair, not one launch each. Binning at
   ! (type, size) granularity fragments the work badly -- 484 launches
   ! averaging 2200 items apiece measured 12x slower than a single large
   ! launch. The screening granularity is worth keeping; the launch
   ! granularity is not. A prefix sum lets a thread find its own segment.
   !
   subroutine bins_plan_build(this, b, thresh, bin_pbase, use_density)
      class(trc_bin_plan_t), intent(inout) :: this
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer,  intent(in) :: bin_pbase(:)
      !! (nbin) first slot - 1 of each bin in the bin-interleaved tables
      !! (trc_pair_quartets_t's bin_pbase).
      logical,  intent(in) :: use_density
      !! Also apply the bin-level density screen (density_reject).

      integer :: ia, ib, ka, kb, smax_keep, nA, nB, nseg, is
      integer(int64) :: nt

      call this%release()

      ! In buckets, not decades: the bins carry the resolution they were
      ! built at, and the admission test has to speak the same units.
      smax_keep = int(-log10(thresh))*b%sres

      ! pass 1: count admitted bin pairs
      nseg = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            if (use_density) then
               if (density_reject(b, ka, kb, thresh)) cycle
            end if
            nseg = nseg + 1
         end do
      end do

      allocate (this%sA(nseg), this%sNB(nseg), this%sOA(nseg), this%sOB(nseg))
      allocate (this%sPA(nseg), this%sPB(nseg), this%sD(nseg))
      allocate (this%sLA(nseg), this%sLB(nseg), this%sLC(nseg), this%sLD(nseg))
      allocate (this%sOff(nseg + 1))

      ! pass 2: fill the descriptors and prefix-sum the work
      is = 0
      this%sOff(1) = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         nA = b%bin_cnt(ka)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            nB = b%bin_cnt(kb)
            if (nA == 0 .or. nB == 0) cycle
            if (use_density) then
               if (density_reject(b, ka, kb, thresh)) cycle
            end if
            is = is + 1
            this%sA(is) = nA
            this%sNB(is) = nB
            this%sOA(is) = b%bin_off(ka)
            this%sOB(is) = b%bin_off(kb)
            this%sPA(is) = bin_pbase(ka)
            this%sPB(is) = bin_pbase(kb)
            this%sD(is) = (ka == kb)
            this%sLA(is) = b%bin_la(ka)
            this%sLB(is) = b%bin_lb(ka)
            this%sLC(is) = b%bin_la(kb)
            this%sLD(is) = b%bin_lb(kb)
            if (ka == kb) then
               nt = int(nA, int64)*int(nA + 1, int64)/2
            else
               nt = int(nA, int64)*int(nB, int64)
            end if
            !
            ! The kernels take their within-segment index as a default integer,
            ! so a segment must fit in one. Widening it would cost registers
            ! in the hot kernel, which are the binding constraint; instead
            ! bins_build_binned_pairs caps a bin at BIN_CAP = 32768 pairs
            ! (split_large_bins), so no segment exceeds 2**30. This guard is
            ! the backstop if that cap is ever raised.
            !
            if (nt > 2147483000_int64) error stop &
               "terco: bin-pair segment exceeds a default integer; widen the kernels' t"
            this%sOff(is + 1) = this%sOff(is) + nt
         end do
      end do
      this%nseg = nseg
      this%nwork = this%sOff(nseg + 1)
   end subroutine bins_plan_build

   !
   ! Sort the segments by class, so each class occupies a contiguous run and
   ! one kernel per class covers it. Every thread in a launch then executes
   ! code compiled for exactly its (la,lb,lc,ld) -- no select case, and
   ! ptxas budgets registers per class instead of for the worst class.
   !
   ! Stable argsort (pic sort_index) on the class key, so the segments of one
   ! class keep their admission order. The prefix sum is rebuilt in the new
   ! order.
   !
   subroutine bins_plan_sort_by_class(this, radix)
      class(trc_bin_plan_t), intent(inout) :: this
      integer, intent(in) :: radix
      !! The class key is ((la*radix + lb)*radix + lc)*radix + ld: the
      !! radix of the dispatcher that will read the key.
      integer(int_index), allocatable :: ord(:)
      integer, allocatable :: ckey(:), t(:)
      logical, allocatable :: tl(:)
      integer :: a, k, n

      n = this%nseg
      allocate (ord(n), ckey(n), t(n), tl(n))
      do a = 1, n
         ckey(a) = this%class_key(a, radix)
      end do
      ! Stable like the insertion sort it replaces: equal keys keep segment
      ! order, so the plan is bit-identical.
      call sort_index(ckey, ord)

      t = this%sA(ord)
      this%sA = t
      t = this%sNB(ord)
      this%sNB = t
      t = this%sOA(ord)
      this%sOA = t
      t = this%sOB(ord)
      this%sOB = t
      t = this%sPA(ord)
      this%sPA = t
      t = this%sPB(ord)
      this%sPB = t
      t = this%sLA(ord)
      this%sLA = t
      t = this%sLB(ord)
      this%sLB = t
      t = this%sLC(ord)
      this%sLC = t
      t = this%sLD(ord)
      this%sLD = t
      tl = this%sD(ord)
      this%sD = tl
      this%sOff(1) = 0
      do k = 1, n
         if (this%sD(k)) then
            this%sOff(k + 1) = this%sOff(k) + int(this%sA(k), int64)*int(this%sA(k) + 1, int64)/2
         else
            this%sOff(k + 1) = this%sOff(k) + int(this%sA(k), int64)*int(this%sNB(k), int64)
         end if
      end do
   end subroutine bins_plan_sort_by_class

   pure integer function bins_plan_class_key(this, seg, radix) result(res)
      !! Segment `seg`'s class key at this radix (see sort_by_class).
      class(trc_bin_plan_t), intent(in) :: this
      integer, intent(in) :: seg, radix
      res = ((this%sLA(seg)*radix + this%sLB(seg))*radix + this%sLC(seg))*radix + this%sLD(seg)
   end function bins_plan_class_key

   subroutine bins_plan_release(this)
      class(trc_bin_plan_t), intent(inout) :: this
      if (allocated(this%sA)) deallocate (this%sA, this%sNB, this%sOA, this%sOB, this%sPA, this%sPB, &
                                          this%sD, this%sLA, this%sLB, this%sLC, this%sLD, this%sOff)
      this%nseg = 0
      this%nwork = 0
   end subroutine bins_plan_release

   !
   ! Density screen at bin-pair granularity.
   !
   ! A quartet's Fock contribution is bounded by Q_ab Q_cd times the largest
   ! density block it can reach.  Over a bin pair the reachable blocks are
   ! (ab), (cd) and the four cross combinations, and only the first two are
   ! known bin-locally -- the cross terms would need every (i,k) shell pair.
   ! Using the two that are known keeps this a bin-level test; it is therefore
   ! a partial screen, and deliberately conservative.
   !
   ! NOT ARMED, and it cannot be as written: bin_dm(k) bounds the density
   ! blocks of pairs WITHIN bin k, which bounds the Coulomb terms D_ab and
   ! D_cd but says nothing about D_ac, D_ad, D_bc, D_bd -- the exchange
   ! blocks, between a shell of the bra pair and a shell of the ket pair.
   ! Arming it on this bound drops real exchange contributions: measured on
   ! the 123-atom silica slice in cc-pVDZ, the third SCF energy moved by
   ! 0.5 Hartree. A sound bin-level test would need the largest density
   ! block between the two bins' shell sets, which is not a per-bin
   ! quantity. The per-quartet test in the kernel, which reads dsh for all
   ! six blocks, is the one that screens on the density.
   pure logical function density_reject(b, ka, kb, thresh) result(res)
      type(pair_bins_t), intent(in) :: b
      integer,  intent(in) :: ka, kb
      real(dp), intent(in) :: thresh
      real(dp) :: qq, dd
      res = .false.
      if (.not. allocated(b%bin_dm)) return
      ! bound on Q_ab Q_cd from the two size buckets, back off one decade each
      ! because the bucket is a ceiling
      qq = 10.0_dp**(-real(max(0, b%bin_s(ka) - 1), dp)) &
           *10.0_dp**(-real(max(0, b%bin_s(kb) - 1), dp))
      dd = 4.0_dp*max(b%bin_dm(ka), b%bin_dm(kb))
      res = (qq*dd <= thresh)
   end function density_reject

end module trc_bins
