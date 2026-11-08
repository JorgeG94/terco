!
! Exchange tiles: density-aware pruning of the bin launch plan, before launch.
!
! A segment of the launch plan is a whole bin pair, admitted on the Schwarz
! bound alone (trc_bins), and the density screen happens in the kernel, per
! thread. That is right for Coulomb, whose weight D_ab Q_ab Q_cd does not fall
! off with the distance between the two pairs, but exchange couples them
! through D_ac, D_ad, D_bc and D_bd, which do. When the MD J engine forms J
! (any basis with an f shell) the bins form exchange alone, and on an extended
! molecule most of the launched threads only find out that their quartet is
! negligible: measured on gly30 cc-pVTZ (thresh 3.16e-11, second SCF
! iteration), 27% of the 1.74e10 launched quartets survived the exchange
! screen.
!
! TILES. Each bin is cut into tiles: runs of at most TILE_PAIRS consecutive
! pairs whose leading shell stays within TILE_SPAN shell blocks (TILE_BLOCK
! consecutive shells each). A tile keeps its largest Schwarz bound and the
! hull [tb1, tb2] of the blocks its pairs touch. Per Fock build, the
! density bound dsh is reduced to its block maxima db(I, J), and for every
! level m (L_m = 10^-m) and block I the hull of the blocks J with
! db(I, J) > L_m; a tile's reach at level m is the union of those hulls over
! its own blocks.
!
! THE TEST. For a bra tile A and a ket tile B, every quartet's exchange
! weight is at most Q_A Q_B max db over A's blocks x B's blocks. Take the
! level m with L_m <= thresh/(Q_A Q_B): if B's block hull misses A's reach at
! m, every block pair has db <= L_m, so every quartet of A x B has
! Q_ab Q_cd D <= thresh and fails the kernel's own screen. Dropping the tile
! pair therefore drops nothing the kernel would have kept: the energies are
! the same, quartet for quartet. The kernel's per-quartet screen still runs on
! what survives.
!
! A bin pair whose surviving tile area is most of its full area is launched
! whole, as before: on a compact system, where little dies, tiling would only
! cut the launch into more and smaller segments.
!
! MEASURED (V100, RHF cc-pVTZ, thresh 3.16e-11, second SCF iteration, full
! J + K build; the J part is MD J's and unchanged): gly30 28.5 -> 26.6 s and
! gly50 59.7 -> 53.9 s against the same build untiled, the launch cut from
! 1.74e10 to 1.14e10 and from 4.88e10 to 2.10e10 quartets with the kernels'
! survivors unchanged. On the first iteration's SAD density, 6.3 -> 3.3 s
! and 14.6 -> 6.8 s. The rejected threads are cheap (about 0.1 ns each),
! which is why cutting 35-57% of the launch buys 7-10%, not more. A tile's
! shape is a trade between its hull and the segments it makes: 64 pairs
! with no span limit gave 28.3 / 57.3 s, span 4 33.8 / 67.4 s (tiny tiles,
! warps no longer share a ket pair), span 16 27.0 / 54.8 s, span 32 the
! times above; 32-pair tiles were no better than 64.
!
! LOCALITY comes from the shell order: tiles are runs of a bin's pairs in
! their canonical order, which follows the shell index, which follows the
! order of the atoms in the input. A chain listed residue by residue (the gly
! inputs) gives tiles that are compact in space; a scrambled input gives wide
! hulls and nothing pruned, never a wrong answer.
!
! The tiles live in trc_pair_quartets_t (trc_eri_data): like its slots they
! follow the pair order inside the bins and are rebuilt with them.
!
! Only exchange-only launches are tiled (jfac = 0, kfac /= 0): a Coulomb
! launch keeps nearly every bin pair, and TRC_EXCHANGE_TILES=0 turns the
! tiling off for measurement.
!
module trc_bin_tiles
   use pic_types, only: dp, int64
   use trc_bins, only: pair_bins_t, trc_bin_plan_t
   use trc_eri_data, only: trc_bin_tiles_t
   implicit none
   private

   public :: bin_tiles_requested, bin_tiles_set_enabled, bin_tiles_build, bin_tiles_reach, bin_tiles_plan

   integer, parameter :: TILE_PAIRS = 64
   !! Pairs per tile: two warps of bra pairs against one ket pair, so the
   !! ket-uniform decode of the kernels keeps its broadcast loads.
   integer, parameter :: TILE_BLOCK = 4
   !! Shells per density block.
   integer, parameter :: TILE_SPAN = 32
   !! Blocks the leading shell of a tile's pairs may cover (cut_bin).
   integer, parameter :: TILE_NLEV = 24
   !! Levels L_m = 10^-m, m = 0..TILE_NLEV. A tile pair whose cut lies below
   !! 10^-TILE_NLEV is not pruned.
   real(dp), parameter :: TILE_KEEP_MAX = 0.8_dp
   !! A bin pair is launched whole when its tiles keep more than this
   !! fraction of its quartets.

   logical :: tiles_enabled = .true.
   !! bin_tiles_set_enabled's switch; TRC_EXCHANGE_TILES=0 overrides it off.

   type :: segment_t
      !! One segment of the plan: pairs oa+1..oa+na of bin ka against
      !! ob+1..ob+nb of bin kb, a triangle when diag.
      integer :: ka = 0, kb = 0, oa = 0, na = 0, ob = 0, nb = 0
      logical :: diag = .false.
   end type segment_t

contains

   logical function bin_tiles_requested() result(res)
      !! On unless bin_tiles_set_enabled(.false.) or TRC_EXCHANGE_TILES=0
      !! turned it off.
      character(len=8) :: e
      e = " "
      call get_environment_variable("TRC_EXCHANGE_TILES", e)
      res = tiles_enabled .and. trim(e) /= "0"
   end function bin_tiles_requested

   subroutine bin_tiles_set_enabled(on)
      !! Turn the tiling on or off for the builds that follow, so a check can
      !! compare the two plans on one geometry.
      logical, intent(in) :: on
      tiles_enabled = on
   end subroutine bin_tiles_set_enabled

   subroutine bin_tiles_build(this, b, nbas)
      !! Cut every bin of `b` into tiles. Geometry only: once per bins, and
      !! again whenever the pair order inside the bins changes.
      type(trc_bin_tiles_t), intent(inout) :: this
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: nbas

      integer :: k, t, p, lo, hi, nt

      this = trc_bin_tiles_t()
      this%nbin = b%nbin
      this%nbas = nbas
      this%nblk = (nbas + TILE_BLOCK - 1)/TILE_BLOCK
      allocate (this%bt0(b%nbin), this%btn(b%nbin), this%qbin(b%nbin))
      nt = 0
      do k = 1, b%nbin
         this%bt0(k) = nt
         this%btn(k) = 0
         call cut_bin(b, k, .false., this, nt)
         this%btn(k) = nt - this%bt0(k)
      end do
      this%ntile = nt
      allocate (this%t0(nt), this%tn(nt), this%tb1(nt), this%tb2(nt), this%tq(nt))
      allocate (this%pmax2(nt), this%smin1(nt))
      allocate (this%rlo(0:TILE_NLEV, nt), this%rhi(0:TILE_NLEV, nt))
      allocate (this%db(this%nblk, this%nblk))
      this%qbin = 0.0_dp
      nt = 0
      do k = 1, b%nbin
         call cut_bin(b, k, .true., this, nt)
      end do
      do k = 1, b%nbin
         do t = this%bt0(k) + 1, this%bt0(k) + this%btn(k)
            lo = this%nblk
            hi = 1
            this%tq(t) = 0.0_dp
            do p = b%bin_off(k) + this%t0(t) + 1, b%bin_off(k) + this%t0(t) + this%tn(t)
               lo = min(lo, block_of(b%sp_i(p)), block_of(b%sp_j(p)))
               hi = max(hi, block_of(b%sp_i(p)), block_of(b%sp_j(p)))
               this%tq(t) = max(this%tq(t), b%sp_q(p))
            end do
            this%tb1(t) = lo
            this%tb2(t) = hi
            this%qbin(k) = max(this%qbin(k), this%tq(t))
         end do
         if (this%btn(k) == 0) cycle
         lo = this%bt0(k) + 1
         hi = this%bt0(k) + this%btn(k)
         this%pmax2(lo) = this%tb2(lo)
         do t = lo + 1, hi
            this%pmax2(t) = max(this%pmax2(t - 1), this%tb2(t))
         end do
         this%smin1(hi) = this%tb1(hi)
         do t = hi - 1, lo, -1
            this%smin1(t) = min(this%smin1(t + 1), this%tb1(t))
         end do
      end do
      this%built = .true.
   end subroutine bin_tiles_build

   subroutine cut_bin(b, k, fill, this, nt)
      !! The tiles of bin k: runs of consecutive pairs, each ended at
      !! TILE_PAIRS pairs or when the leading shell (the larger index, which
      !! the canonical pair order ascends) has moved TILE_SPAN blocks past the
      !! tile's first. Counted into nt, or (fill) also written from tile
      !! nt + 1. The span cut is what keeps a tile compact in a sparse bin:
      !! 64 pairs of a bin with a handful per atom would span the molecule.
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: k
      logical, intent(in) :: fill
      type(trc_bin_tiles_t), intent(inout) :: this
      integer, intent(inout) :: nt
      integer :: p, n, lead, first, start

      n = 0
      first = 0
      start = 0
      do p = 1, b%bin_cnt(k)
         lead = block_of(max(b%sp_i(b%bin_off(k) + p), b%sp_j(b%bin_off(k) + p)))
         if (n == 0 .or. n == TILE_PAIRS .or. lead >= first + TILE_SPAN) then
            if (n > 0 .and. fill) then
               this%t0(nt) = start
               this%tn(nt) = n
            end if
            nt = nt + 1
            n = 0
            first = lead
            start = p - 1
         end if
         n = n + 1
      end do
      if (n > 0 .and. fill) then
         this%t0(nt) = start
         this%tn(nt) = n
      end if
   end subroutine cut_bin

   pure integer function block_of(shell) result(res)
      integer, intent(in) :: shell
      res = (shell - 1)/TILE_BLOCK + 1
   end function block_of

   subroutine bin_tiles_reach(this, dsh)
      !! This build's reach of every tile, from the shell-pair density bound.
      !! `dsh` must be present on the device; only its block maxima come back.
      type(trc_bin_tiles_t), intent(inout) :: this
      real(dp), intent(in) :: dsh(this%nbas, this%nbas)

      integer, allocatable :: hlo(:, :), hhi(:, :)
      integer :: i, j, m, mv, t, nblk

      nblk = this%nblk
      !$acc enter data create(this%db)
      call block_maxima(this%nbas, nblk, dsh, this%db)
      !$acc update self(this%db)
      !$acc exit data delete(this%db)

      ! Per block row I and level m: the hull of the blocks J with
      ! db(I, J) > 10^-m. The sets grow with m, so a scan from each end
      ! stops at the first level already claimed.
      allocate (hlo(0:TILE_NLEV, nblk), hhi(0:TILE_NLEV, nblk))
      hlo = nblk + 1
      hhi = 0
      do i = 1, nblk
         do j = 1, nblk
            mv = lowest_level(this%db(i, j))
            do m = mv, TILE_NLEV
               if (hlo(m, i) <= nblk) exit
               hlo(m, i) = j
            end do
         end do
         do j = nblk, 1, -1
            mv = lowest_level(this%db(i, j))
            do m = mv, TILE_NLEV
               if (hhi(m, i) > 0) exit
               hhi(m, i) = j
            end do
         end do
      end do

      do t = 1, this%ntile
         do m = 0, TILE_NLEV
            this%rlo(m, t) = minval(hlo(m, this%tb1(t):this%tb2(t)))
            this%rhi(m, t) = maxval(hhi(m, this%tb1(t):this%tb2(t)))
         end do
      end do
      deallocate (hlo, hhi)
   end subroutine bin_tiles_reach

   pure integer function lowest_level(v) result(m)
      !! A level m <= the smallest with 10^-m < v, never above it: a block
      !! enters the hulls at m and every later level, so rounding in log10
      !! may only add it early. TILE_NLEV + 1 for a block that is zero.
      real(dp), intent(in) :: v
      if (v <= 0.0_dp) then
         m = TILE_NLEV + 1
      else if (v >= 1.0_dp) then
         m = 0
      else
         m = min(TILE_NLEV + 1, int(floor(-log10(v))))
      end if
   end function lowest_level

   pure integer function cut_level(cut) result(m)
      !! A level m with 10^-m <= cut: floor(-log10 cut) + 1, which rounding
      !! can only push one level deeper (more conservative), never shallower.
      !! TILE_NLEV + 1 when the cut lies below every level.
      real(dp), intent(in) :: cut
      if (cut >= 1.0_dp) then
         m = 0
      else if (cut <= 0.0_dp) then
         m = TILE_NLEV + 1
      else
         m = min(TILE_NLEV + 1, int(floor(-log10(cut))) + 1)
      end if
   end function cut_level

   subroutine block_maxima(nbas, nblk, dsh, db)
      !! db(I, J) = max dsh over shell block I x shell block J, on the device.
      integer,  intent(in)  :: nbas, nblk
      real(dp), intent(in)  :: dsh(nbas, nbas)
      real(dp), intent(out) :: db(nblk, nblk)
      integer :: ij
      do concurrent(ij=1:nblk*nblk)
         call block_max_item(ij, nbas, nblk, dsh, db)
      end do
   end subroutine block_maxima

   pure subroutine block_max_item(ij, nbas, nblk, dsh, db)
      !$acc routine seq
      integer,  intent(in)    :: ij, nbas, nblk
      real(dp), intent(in)    :: dsh(nbas, nbas)
      real(dp), intent(inout) :: db(nblk, nblk)
      integer :: bi, bj, i, j
      real(dp) :: m
      bi = mod(ij - 1, nblk) + 1
      bj = (ij - 1)/nblk + 1
      m = 0.0_dp
      do j = (bj - 1)*TILE_BLOCK + 1, min(nbas, bj*TILE_BLOCK)
         do i = (bi - 1)*TILE_BLOCK + 1, min(nbas, bi*TILE_BLOCK)
            m = max(m, dsh(i, j))
         end do
      end do
      db(bi, bj) = m
   end subroutine block_max_item

   subroutine bin_tiles_plan(this, plan, b, thresh, bin_pbase)
      !! The launch plan of trc_bin_plan_t%build, with every Schwarz-admitted
      !! bin pair either kept whole or cut into the tile pairs that can carry
      !! exchange at this density (see the module header). Needs this
      !! build's `reach`.
      type(trc_bin_tiles_t), intent(in) :: this
      type(trc_bin_plan_t), intent(inout) :: plan
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bin_pbase(:)
      !! (nbin) first slot - 1 of each bin in the bin-interleaved tables

      integer :: nseg

      call plan%release()
      nseg = 0
      call walk(this, plan, b, thresh, bin_pbase, .false., nseg)
      allocate (plan%sA(nseg), plan%sNB(nseg), plan%sOA(nseg), plan%sOB(nseg))
      allocate (plan%sPA(nseg), plan%sPB(nseg), plan%sSA(nseg), plan%sSB(nseg), plan%sD(nseg))
      allocate (plan%sLA(nseg), plan%sLB(nseg), plan%sLC(nseg), plan%sLD(nseg))
      allocate (plan%sOff(nseg + 1))
      plan%sOff(1) = 0
      nseg = 0
      call walk(this, plan, b, thresh, bin_pbase, .true., nseg)
      plan%nseg = nseg
      plan%nwork = plan%sOff(nseg + 1)
   end subroutine bin_tiles_plan

   subroutine walk(this, plan, b, thresh, bin_pbase, fill, nseg)
      !! One pass over the admitted bin pairs: count the segments, or (fill)
      !! write them. Both passes take the same decisions.
      type(trc_bin_tiles_t), intent(in) :: this
      type(trc_bin_plan_t), intent(inout) :: plan
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bin_pbase(:)
      logical, intent(in) :: fill
      integer, intent(inout) :: nseg

      integer :: ia, ib, ka, kb, smax_keep, ta, tb, t1, t2, tlast, run0, run1, nkept
      integer(int64) :: area, full
      logical :: same, keep

      smax_keep = int(-log10(thresh))*b%sres
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            same = ka == kb
            if (same) then
               full = int(b%bin_cnt(ka), int64)*int(b%bin_cnt(ka) + 1, int64)/2
            else
               full = int(b%bin_cnt(ka), int64)*int(b%bin_cnt(kb), int64)
            end if
            ! First pass over the tiles: the area they keep, to decide whether
            ! the bin pair is worth cutting at all.
            area = 0
            nkept = 0
            if (this%btn(ka) > 1 .or. this%btn(kb) > 1) then
               do ta = this%bt0(ka) + 1, this%bt0(ka) + this%btn(ka)
                  call candidates(this, ta, kb, thresh, t1, t2)
                  if (same) t2 = min(t2, ta)
                  do tb = t1, t2
                     if (.not. tile_pair_kept(this, ta, tb, thresh)) cycle
                     nkept = nkept + 1
                     if (same .and. tb == ta) then
                        area = area + int(this%tn(ta), int64)*int(this%tn(ta) + 1, int64)/2
                     else
                        area = area + int(this%tn(ta), int64)*int(this%tn(tb), int64)
                     end if
                  end do
               end do
            else
               area = full
            end if
            if (real(area, dp) > TILE_KEEP_MAX*real(full, dp)) then
               call emit(plan, b, bin_pbase, fill, nseg, segment_t(ka, kb, 0, b%bin_cnt(ka), 0, b%bin_cnt(kb), same))
               cycle
            end if
            if (nkept == 0) cycle
            ! Cut: per bra tile, each run of consecutive kept ket tiles is one
            ! rectangle; on a bin against itself the diagonal tile is its own
            ! triangle.
            do ta = this%bt0(ka) + 1, this%bt0(ka) + this%btn(ka)
               call candidates(this, ta, kb, thresh, t1, t2)
               if (same) t2 = min(t2, ta)
               run0 = 0
               run1 = -1
               tlast = t2
               if (same .and. t2 == ta) tlast = ta - 1
               do tb = t1, tlast + 1
                  keep = .false.
                  if (tb <= tlast) keep = tile_pair_kept(this, ta, tb, thresh)
                  if (keep) then
                     if (run1 < run0) run0 = tb
                     run1 = tb
                  else if (run1 >= run0) then
                     call emit(plan, b, bin_pbase, fill, nseg, &
                               segment_t(ka, kb, this%t0(ta), this%tn(ta), this%t0(run0), &
                                         this%t0(run1) + this%tn(run1) - this%t0(run0), .false.))
                     run0 = 0
                     run1 = -1
                  end if
               end do
               if (same .and. t2 == ta) then
                  if (tile_pair_kept(this, ta, ta, thresh)) then
                     call emit(plan, b, bin_pbase, fill, nseg, &
                               segment_t(ka, kb, this%t0(ta), this%tn(ta), this%t0(ta), this%tn(ta), .true.))
                  end if
               end if
            end do
         end do
      end do
   end subroutine walk

   subroutine candidates(this, ta, kb, thresh, t1, t2)
      !! The tiles t1..t2 of bin kb that can meet bra tile ta at all: those
      !! whose block hull overlaps ta's reach at the loosest cut the bin
      !! allows (its largest Schwarz bound). Empty when t1 > t2.
      type(trc_bin_tiles_t), intent(in) :: this
      integer, intent(in) :: ta, kb
      real(dp), intent(in) :: thresh
      integer, intent(out) :: t1, t2
      integer :: m, lo, hi, a, z, mid

      t1 = this%bt0(kb) + 1
      t2 = this%bt0(kb) + this%btn(kb)
      if (this%tq(ta)*this%qbin(kb) <= thresh) then
         t2 = t1 - 1
         return
      end if
      m = cut_level(thresh/(this%tq(ta)*this%qbin(kb)))
      if (m > TILE_NLEV) return
      lo = this%rlo(m, ta)
      hi = this%rhi(m, ta)
      if (lo > hi) then
         t2 = t1 - 1
         return
      end if
      ! first tile whose running max of tb2 reaches lo
      a = t1
      z = t2 + 1
      do while (a < z)
         mid = (a + z)/2
         if (this%pmax2(mid) >= lo) then
            z = mid
         else
            a = mid + 1
         end if
      end do
      t1 = a
      ! last tile whose running min of tb1 from the end is within hi
      a = this%bt0(kb)
      z = t2
      do while (a < z)
         mid = (a + z + 1)/2
         if (this%smin1(mid) <= hi) then
            a = mid
         else
            z = mid - 1
         end if
      end do
      t2 = a
   end subroutine candidates

   pure logical function tile_pair_kept(this, ta, tb, thresh) result(res)
      !! Whether bra tile ta and ket tile tb can hold a quartet that passes
      !! the kernel's exchange screen (see the module header).
      type(trc_bin_tiles_t), intent(in) :: this
      integer, intent(in) :: ta, tb
      real(dp), intent(in) :: thresh
      integer :: m
      res = .false.
      if (this%tq(ta)*this%tq(tb) <= thresh) return
      m = cut_level(thresh/(this%tq(ta)*this%tq(tb)))
      res = .true.
      if (m > TILE_NLEV) return
      res = this%rlo(m, ta) <= this%tb2(tb) .and. this%rhi(m, ta) >= this%tb1(tb)
   end function tile_pair_kept

   subroutine emit(plan, b, bin_pbase, fill, nseg, sg)
      !! Count one segment, or (fill) write it as segment nseg.
      type(trc_bin_plan_t), intent(inout) :: plan
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: bin_pbase(:)
      logical, intent(in) :: fill
      integer, intent(inout) :: nseg
      type(segment_t), intent(in) :: sg
      integer(int64) :: nt
      integer :: ka, kb, oa, na, ob, nb
      logical :: diag

      nseg = nseg + 1
      if (.not. fill) return
      ka = sg%ka
      kb = sg%kb
      oa = sg%oa
      na = sg%na
      ob = sg%ob
      nb = sg%nb
      diag = sg%diag
      plan%sA(nseg) = na
      plan%sNB(nseg) = nb
      plan%sOA(nseg) = b%bin_off(ka) + oa
      plan%sOB(nseg) = b%bin_off(kb) + ob
      plan%sPA(nseg) = bin_pbase(ka) + oa
      plan%sPB(nseg) = bin_pbase(kb) + ob
      plan%sSA(nseg) = b%bin_cnt(ka)
      plan%sSB(nseg) = b%bin_cnt(kb)
      plan%sD(nseg) = diag
      plan%sLA(nseg) = b%bin_la(ka)
      plan%sLB(nseg) = b%bin_lb(ka)
      plan%sLC(nseg) = b%bin_la(kb)
      plan%sLD(nseg) = b%bin_lb(kb)
      if (diag) then
         nt = int(na, int64)*int(na + 1, int64)/2
      else
         nt = int(na, int64)*int(nb, int64)
      end if
      ! See bins_plan_build: the kernels index within a segment by a
      ! default integer, and bins are capped so that no segment needs more.
      if (nt > 2147483000_int64) error stop &
         "terco: bin-pair segment exceeds a default integer; widen the kernels' t"
      plan%sOff(nseg + 1) = plan%sOff(nseg) + nt
   end subroutine emit

end module trc_bin_tiles
