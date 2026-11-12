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
   use trc_bins, only: pair_bins_t, trc_bin_plan_t, bins_skel_cut_first, bins_plan_emit, bins_plan_allocate, bins_any_masked, &
                       trc_bin_segment_t
   use trc_eri_data, only: trc_bin_tiles_t, trc_bin_orbits_t
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

   subroutine bin_tiles_build(this, b, nbas, orb)
      !! Cut every bin of `b` into tiles. Geometry only: once per bins, and
      !! again whenever the pair order inside the bins changes.
      !!
      !! With `orb` (bins closed under a point group, trc_bins) a tile never
      !! crosses a boundary between two orbit sizes or between the
      !! representatives and the other pairs, so that a representative tile
      !! has one weight; the representative tiles come first in a bin.
      type(trc_bin_tiles_t), intent(inout) :: this
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: nbas
      type(trc_bin_orbits_t), intent(in), optional :: orb
      logical :: skel

      integer :: k, t, p, lo, hi, nt

      this = trc_bin_tiles_t()
      skel = .false.
      if (present(orb)) skel = orb%built
      this%nbin = b%nbin
      this%nbas = nbas
      this%nblk = (nbas + TILE_BLOCK - 1)/TILE_BLOCK
      allocate (this%bt0(b%nbin), this%btn(b%nbin), this%qbin(b%nbin))
      nt = 0
      do k = 1, b%nbin
         this%bt0(k) = nt
         this%btn(k) = 0
         if (skel) then
            call cut_bin(b, k, .false., this, nt, orb)
         else
            call cut_bin(b, k, .false., this, nt)
         end if
         this%btn(k) = nt - this%bt0(k)
      end do
      this%ntile = nt
      allocate (this%t0(nt), this%tn(nt), this%tb1(nt), this%tb2(nt), this%tq(nt))
      allocate (this%pmax2(nt), this%smin1(nt))
      allocate (this%rlo(0:TILE_NLEV, nt), this%rhi(0:TILE_NLEV, nt))
      allocate (this%db(this%nblk, this%nblk))
      this%qbin = 0.0_dp
      if (skel) then
         allocate (this%brep(b%nbin), this%tw(nt))
         this%brep = 0
         this%tw = 0
      end if
      nt = 0
      do k = 1, b%nbin
         if (skel) then
            call cut_bin(b, k, .true., this, nt, orb)
         else
            call cut_bin(b, k, .true., this, nt)
         end if
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

   subroutine cut_bin(b, k, fill, this, nt, orb)
      !! The tiles of bin k: runs of consecutive pairs, each ended at
      !! TILE_PAIRS pairs or when the leading shell (the larger index, which
      !! the canonical pair order ascends) has moved TILE_SPAN blocks past the
      !! tile's first. Counted into nt, or (fill) also written from tile
      !! nt + 1. The span cut is what keeps a tile compact in a sparse bin:
      !! 64 pairs of a bin with a handful per atom would span the molecule.
      !! With `orb` the bin is cut one range at a time (each weight of the
      !! representatives, then the other pairs), and a representative tile
      !! records its weight and counts in `brep`.
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: k
      logical, intent(in) :: fill
      type(trc_bin_tiles_t), intent(inout) :: this
      integer, intent(inout) :: nt
      type(trc_bin_orbits_t), intent(in), optional :: orb
      integer :: r

      if (.not. present(orb)) then
         call cut_range(b, k, 0, b%bin_cnt(k), 0, fill, this, nt)
         return
      end if
      do r = orb%sec0(k) + 1, orb%sec0(k + 1)
         call cut_range(b, k, orb%sec_beg(r), orb%sec_cnt(r), orb%sec_w(r), fill, this, nt)
      end do
      if (fill) this%brep(k) = nt - this%bt0(k)
      if (orb%bin_nrep(k) < b%bin_cnt(k)) then
         call cut_range(b, k, orb%bin_nrep(k), b%bin_cnt(k) - orb%bin_nrep(k), 0, fill, this, nt)
      end if
   end subroutine cut_bin

   subroutine cut_range(b, k, beg, cnt, w, fill, this, nt)
      !! The tiles of the pairs beg+1..beg+cnt of bin k (cut_bin's rule); w is
      !! the weight recorded for them.
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: k, beg, cnt, w
      logical, intent(in) :: fill
      type(trc_bin_tiles_t), intent(inout) :: this
      integer, intent(inout) :: nt
      integer :: p, n, lead, first, start

      n = 0
      first = 0
      start = 0
      do p = beg + 1, beg + cnt
         lead = block_of(max(b%sp_i(b%bin_off(k) + p), b%sp_j(b%bin_off(k) + p)))
         if (n == 0 .or. n == TILE_PAIRS .or. lead >= first + TILE_SPAN) then
            if (n > 0 .and. fill) then
               this%t0(nt) = start
               this%tn(nt) = n
               if (allocated(this%tw)) this%tw(nt) = w
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
         if (allocated(this%tw)) this%tw(nt) = w
      end if
   end subroutine cut_range

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

   subroutine bin_tiles_plan(this, plan, b, thresh, bin_pbase, orb)
      !! The launch plan of trc_bin_plan_t%build, with every Schwarz-admitted
      !! bin pair either kept whole or cut into the tile pairs that can carry
      !! exchange at this density (see the module header). Needs this
      !! build's `reach`.
      !!
      !! With `orb` the plan is the skeleton plan of trc_bin_plan_t%build_skeleton
      !! (one side of each bin pair restricted to representatives, weights in
      !! sW), cut the same way: the tiles are those of the representative
      !! ranges against the whole other bin.
      !!
      !! A plan with orbit weights (skeleton) or parity masks is built here,
      !! by the host walk; any other goes to plan_device, which builds the
      !! same segments on the device.
      type(trc_bin_tiles_t), intent(in) :: this
      type(trc_bin_plan_t), intent(inout) :: plan
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bin_pbase(:)
      !! (nbin) first slot - 1 of each bin in the bin-interleaved tables
      type(trc_bin_orbits_t), intent(in), optional :: orb

      integer :: nseg
      logical :: skel

      skel = .false.
      if (present(orb)) skel = orb%built
      call plan%release()
      if (.not. (skel .or. bins_any_masked(b))) then
         ! No orbit weights and no parity masks: the device planner.
         call plan_device(this, plan, b, thresh, bin_pbase)
         return
      end if
      nseg = 0
      if (skel) then
         call walk_skeleton(this, plan, b, orb, thresh, bin_pbase, .false., nseg)
      else
         call walk(this, plan, b, thresh, bin_pbase, .false., nseg)
      end if
      call bins_plan_allocate(plan, nseg, skel)
      nseg = 0
      if (skel) then
         call walk_skeleton(this, plan, b, orb, thresh, bin_pbase, .true., nseg)
      else
         call walk(this, plan, b, thresh, bin_pbase, .true., nseg)
      end if
      plan%nseg = nseg
      plan%nwork = plan%sOff(nseg + 1)
   end subroutine bin_tiles_plan

   subroutine walk_skeleton(this, plan, b, orb, thresh, bin_pbase, fill, nseg)
      !! walk for the skeleton plan. For each admitted bin pair the bra side
      !! kx is the one bins_skel_cut_first picks, and only its representative
      !! tiles (the leading `brep` tiles) are bra tiles; the ket side is all
      !! of the other bin. A bin against itself is the same rectangle, with
      !! the half of the diagonal bin pair in the weight. Both passes take
      !! the same decisions.
      type(trc_bin_tiles_t), intent(in) :: this
      type(trc_bin_plan_t), intent(inout) :: plan
      type(pair_bins_t), intent(in) :: b
      type(trc_bin_orbits_t), intent(in) :: orb
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bin_pbase(:)
      logical, intent(in) :: fill
      integer, intent(inout) :: nseg

      integer :: ia, ib, ka, kb, kx, ky, smax_keep, ta, tb, t1, t2, run0, run1, nkept, r, ntx, tlast
      integer(int64) :: area, full
      logical :: keep
      real(dp) :: half, w

      smax_keep = int(-log10(thresh))*b%sres
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            kx = ka
            ky = kb
            if (.not. bins_skel_cut_first(b, orb, ka, kb)) then
               kx = kb
               ky = ka
            end if
            half = merge(0.5_dp, 1.0_dp, ka == kb)
            ntx = this%brep(kx)
            full = int(orb%bin_nrep(kx), int64)*int(b%bin_cnt(ky), int64)
            area = 0
            nkept = 0
            if (this%btn(kx) > 1 .or. this%btn(ky) > 1) then
               do ta = this%bt0(kx) + 1, this%bt0(kx) + ntx
                  call candidates(this, ta, ky, thresh, t1, t2)
                  do tb = t1, t2
                     if (.not. tile_pair_kept(this, ta, tb, thresh)) cycle
                     nkept = nkept + 1
                     area = area + int(this%tn(ta), int64)*int(this%tn(tb), int64)
                  end do
               end do
            else
               area = full
            end if
            if (real(area, dp) > TILE_KEEP_MAX*real(full, dp)) then
               do r = orb%sec0(kx) + 1, orb%sec0(kx + 1)
                  ! Not canonicalised (bins_plan_emit's `canon`): this is a
                  ! Fock exchange plan, whose kernels exist in both bra/ket
                  ! orderings, and reordering here would only move its
                  ! accumulation order.
                  call bins_plan_emit(plan, b, bin_pbase, fill, nseg, &
                                      trc_bin_segment_t(kx, ky, orb%sec_beg(r), orb%sec_cnt(r), 0, b%bin_cnt(ky), &
                                                        .false., half*real(orb%sec_w(r), dp)), .false.)
               end do
               cycle
            end if
            if (nkept == 0) cycle
            ! Cut: per representative bra tile, each run of consecutive kept
            ! ket tiles is one rectangle.
            do ta = this%bt0(kx) + 1, this%bt0(kx) + ntx
               call candidates(this, ta, ky, thresh, t1, t2)
               w = half*real(this%tw(ta), dp)
               run0 = 0
               run1 = -1
               tlast = t2
               do tb = t1, tlast + 1
                  keep = .false.
                  if (tb <= tlast) keep = tile_pair_kept(this, ta, tb, thresh)
                  if (keep) then
                     if (run1 < run0) run0 = tb
                     run1 = tb
                  else if (run1 >= run0) then
                     call bins_plan_emit(plan, b, bin_pbase, fill, nseg, &
                                         trc_bin_segment_t(kx, ky, this%t0(ta), this%tn(ta), this%t0(run0), &
                                                           this%t0(run1) + this%tn(run1) - this%t0(run0), .false., w), &
                                         .false.)
                     run0 = 0
                     run1 = -1
                  end if
               end do
            end do
         end do
      end do
   end subroutine walk_skeleton

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
               call emit(plan, b, bin_pbase, fill, nseg, trc_bin_segment_t(ka, kb, 0, b%bin_cnt(ka), 0, b%bin_cnt(kb), same))
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
                               trc_bin_segment_t(ka, kb, this%t0(ta), this%tn(ta), this%t0(run0), &
                                         this%t0(run1) + this%tn(run1) - this%t0(run0), .false.))
                     run0 = 0
                     run1 = -1
                  end if
               end do
               if (same .and. t2 == ta) then
                  if (tile_pair_kept(this, ta, ta, thresh)) then
                     call emit(plan, b, bin_pbase, fill, nseg, &
                               trc_bin_segment_t(ka, kb, this%t0(ta), this%tn(ta), this%t0(ta), this%tn(ta), .true.))
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
      type(trc_bin_segment_t), intent(in) :: sg
      call bins_plan_emit(plan, b, bin_pbase, fill, nseg, sg, .false.)
   end subroutine emit


   subroutine plan_device(this, plan, b, thresh, bin_pbase)
      !! bin_tiles_plan for a plan with neither orbit weights nor parity masks:
      !! every segment's mask (sM) is 0 and sW stays unallocated.
      !!
      !! The launch plan of trc_bin_plan_t%build, with every Schwarz-admitted
      !! bin pair either kept whole or cut into the tile pairs that can carry
      !! exchange at this density (see the module header). Needs this
      !! build's `reach`.
      !!
      !! ON THE DEVICE. The plan was a serial host walk over every admitted
      !! bin pair's tile pairs, twice (count, then fill), with the GPU idle:
      !! 0.14-0.40 s of a 4.4 s Fock build on gly50 cc-pVDZ (V100), growing
      !! as the density spreads. It is now two `do concurrent` passes over
      !! (bin pair, bra tile) slots -- count, then fill -- with only the
      !! per-pair decision and the prefix sums between them on the host,
      !! O(slots). The segments come out in the host walk's order: bin pairs
      !! as `walk` visited them, a cut pair's bra tiles ascending, each
      !! tile's runs ascending, its diagonal triangle last.
      type(trc_bin_tiles_t), intent(in) :: this
      type(trc_bin_plan_t), intent(inout) :: plan
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bin_pbase(:)
      !! (nbin) first slot - 1 of each bin in the bin-interleaved tables

      integer, allocatable :: pr_ka(:), pr_kb(:), pr_s0(:), pr_ns(:), pr_mode(:), pr_off(:)
      integer, allocatable :: sl_pr(:), sl_ta(:), sl_kept(:), sl_nrun(:), sl_off(:)
      integer(int64), allocatable :: sl_area(:)
      integer, allocatable :: qA(:), qNB(:), qOA(:), qOB(:), qPA(:), qPB(:), qSA(:), qSB(:)
      integer, allocatable :: qLA(:), qLB(:), qLC(:), qLD(:)
      logical, allocatable :: qD(:)
      ! Device copies of what the kernels read, as plain local arrays: the
      ! derived-type components themselves were not found by the do
      ! concurrent loops (nvfortran 26.5 read the host address under
      ! mem:separate), whether `copyin`-mapped or not.
      integer, allocatable :: xbt0(:), xbtn(:), xt0(:), xtn(:), xtb1(:), xtb2(:), xpmax2(:), xsmin1(:)
      integer, allocatable :: xrlo(:, :), xrhi(:, :), xcnt(:), xoff(:), xla(:), xlb(:), xpb(:)
      real(dp), allocatable :: xtq(:), xqbin(:)
      integer :: ia, ib, ka, kb, npr, nslot, p, s, t, nseg, nbin, ntile, smax_keep
      integer(int64) :: area, full, nt
      integer :: nkept

      call plan%release()
      nbin = this%nbin
      ntile = max(this%ntile, 1)

      ! The admitted bin pairs, in walk's order, and their slots: one per bra
      ! tile of a pair that can be cut, none for one that cannot (whole).
      smax_keep = int(-log10(thresh))*b%sres
      npr = 0
      do ia = 1, b%nlive
         do ib = 1, ia
            if (admitted(b, b%live(ia), b%live(ib), smax_keep)) npr = npr + 1
         end do
      end do
      allocate (pr_ka(max(npr, 1)), pr_kb(max(npr, 1)), pr_s0(max(npr, 1)), pr_ns(max(npr, 1)))
      allocate (pr_mode(max(npr, 1)), pr_off(max(npr, 1)))
      npr = 0
      nslot = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (.not. admitted(b, ka, kb, smax_keep)) cycle
            npr = npr + 1
            pr_ka(npr) = ka
            pr_kb(npr) = kb
            pr_s0(npr) = nslot
            pr_ns(npr) = 0
            if (this%btn(ka) > 1 .or. this%btn(kb) > 1) pr_ns(npr) = this%btn(ka)
            nslot = nslot + pr_ns(npr)
         end do
      end do
      allocate (sl_pr(max(nslot, 1)), sl_ta(max(nslot, 1)), sl_kept(max(nslot, 1)), sl_nrun(max(nslot, 1)))
      allocate (sl_off(max(nslot, 1)), sl_area(max(nslot, 1)))
      do p = 1, npr
         do t = 1, pr_ns(p)
            sl_pr(pr_s0(p) + t) = p
            sl_ta(pr_s0(p) + t) = this%bt0(pr_ka(p)) + t
         end do
      end do

      xbt0 = this%bt0
      xbtn = this%btn
      xqbin = this%qbin
      allocate (xt0(ntile), xtn(ntile), xtb1(ntile), xtb2(ntile), xpmax2(ntile), xsmin1(ntile), xtq(ntile))
      allocate (xrlo(0:TILE_NLEV, ntile), xrhi(0:TILE_NLEV, ntile))
      xt0 = 0
      xtn = 0
      xtb1 = 0
      xtb2 = 0
      xpmax2 = 0
      xsmin1 = 0
      xtq = 0.0_dp
      xrlo = 0
      xrhi = 0
      if (this%ntile > 0) then
         xt0 = this%t0
         xtn = this%tn
         xtb1 = this%tb1
         xtb2 = this%tb2
         xpmax2 = this%pmax2
         xsmin1 = this%smin1
         xtq = this%tq
         xrlo = this%rlo
         xrhi = this%rhi
      end if
      xcnt = b%bin_cnt(1:nbin)
      xoff = b%bin_off(1:nbin)
      xla = b%bin_la(1:nbin)
      xlb = b%bin_lb(1:nbin)
      xpb = bin_pbase(1:nbin)
      !$acc enter data copyin(xbt0, xbtn, xt0, xtn, xtb1, xtb2, xpmax2, xsmin1, xtq, xqbin, xrlo, xrhi, &
      !$acc                   pr_ka, pr_kb, sl_pr, sl_ta) &
      !$acc            create(sl_area, sl_kept, sl_nrun)
      ! The item routines take the loop index by value: passed by reference,
      ! optimised ifx 2026.0 ran every iteration as slot 1 (the counts came
      ! back 341 kept tile pairs instead of 662 in check_exchange_screen).
      do concurrent(s=1:nslot)
         call count_item(s, nbin, ntile, npr, nslot, xbt0, xbtn, xtn, &
                         xtb1, xtb2, xpmax2, xsmin1, xtq, xqbin, &
                         xrlo, xrhi, pr_ka, pr_kb, sl_pr, sl_ta, thresh, &
                         sl_area, sl_kept, sl_nrun)
      end do
      !$acc update self(sl_area, sl_kept, sl_nrun)

      ! Per pair, walk's decision; then each segment's place in walk's order.
      nseg = 0
      do p = 1, npr
         ka = pr_ka(p)
         kb = pr_kb(p)
         if (ka == kb) then
            full = int(b%bin_cnt(ka), int64)*int(b%bin_cnt(ka) + 1, int64)/2
         else
            full = int(b%bin_cnt(ka), int64)*int(b%bin_cnt(kb), int64)
         end if
         if (pr_ns(p) == 0) then
            area = full
            nkept = 0
         else
            area = sum(sl_area(pr_s0(p) + 1:pr_s0(p) + pr_ns(p)))
            nkept = sum(sl_kept(pr_s0(p) + 1:pr_s0(p) + pr_ns(p)))
         end if
         pr_off(p) = nseg
         if (real(area, dp) > TILE_KEEP_MAX*real(full, dp)) then
            pr_mode(p) = 1
            nseg = nseg + 1
         else if (nkept == 0) then
            pr_mode(p) = 0
         else
            pr_mode(p) = 2
            do s = pr_s0(p) + 1, pr_s0(p) + pr_ns(p)
               sl_off(s) = nseg
               nseg = nseg + sl_nrun(s)
            end do
         end if
      end do

      ! Filled as plain local arrays and moved into `plan` afterwards: the
      ! plan's own allocatable components, mapped with `create`, were not
      ! found by the do concurrent below (nvfortran 26.5 wrote to the host
      ! address: an illegal address under mem:separate).
      allocate (qA(max(nseg, 1)), qNB(max(nseg, 1)), qOA(max(nseg, 1)), qOB(max(nseg, 1)))
      allocate (qPA(max(nseg, 1)), qPB(max(nseg, 1)), qSA(max(nseg, 1)), qSB(max(nseg, 1)), qD(max(nseg, 1)))
      allocate (qLA(max(nseg, 1)), qLB(max(nseg, 1)), qLC(max(nseg, 1)), qLD(max(nseg, 1)))
      allocate (plan%sOff(nseg + 1))
      !$acc enter data copyin(pr_mode, pr_off, sl_off, xcnt, xoff, xla, xlb, xpb) &
      !$acc            create(qA, qNB, qOA, qOB, qPA, qPB, qSA, qSB, qD, qLA, qLB, qLC, qLD)
      do concurrent(p=1:npr)
         call whole_item(p, nbin, npr, max(nseg, 1), pr_ka, pr_kb, pr_mode, pr_off, &
                         xcnt, xoff, xla, xlb, xpb, &
                         qA, qNB, qOA, qOB, qPA, qPB, qSA, qSB, &
                         qD, qLA, qLB, qLC, qLD)
      end do
      do concurrent(s=1:nslot)
         call fill_item(s, nbin, ntile, npr, nslot, max(nseg, 1), xbt0, xbtn, xt0, xtn, &
                        xtb1, xtb2, xpmax2, xsmin1, xtq, xqbin, &
                        xrlo, xrhi, pr_ka, pr_kb, pr_mode, sl_pr, sl_ta, sl_off, thresh, &
                        xcnt, xoff, xla, xlb, xpb, &
                        qA, qNB, qOA, qOB, qPA, qPB, qSA, qSB, &
                        qD, qLA, qLB, qLC, qLD)
      end do
      !$acc update self(qA, qNB, qOA, qOB, qPA, qPB, qSA, qSB, &
      !$acc             qD, qLA, qLB, qLC, qLD)
      !$acc exit data delete(qA, qNB, qOA, qOB, qPA, qPB, qSA, qSB, &
      !$acc                  qD, qLA, qLB, qLC, qLD, &
      !$acc                  pr_mode, pr_off, sl_off, xcnt, xoff, xla, xlb, xpb, &
      !$acc                  xbt0, xbtn, xt0, xtn, xtb1, xtb2, &
      !$acc                  xpmax2, xsmin1, xtq, xqbin, xrlo, xrhi, &
      !$acc                  pr_ka, pr_kb, sl_pr, sl_ta, sl_area, sl_kept, sl_nrun)

      call move_alloc(qA, plan%sA)
      call move_alloc(qNB, plan%sNB)
      call move_alloc(qOA, plan%sOA)
      call move_alloc(qOB, plan%sOB)
      call move_alloc(qPA, plan%sPA)
      call move_alloc(qPB, plan%sPB)
      call move_alloc(qSA, plan%sSA)
      call move_alloc(qSB, plan%sSB)
      call move_alloc(qD, plan%sD)
      call move_alloc(qLA, plan%sLA)
      call move_alloc(qLB, plan%sLB)
      call move_alloc(qLC, plan%sLC)
      call move_alloc(qLD, plan%sLD)
      plan%sOff(1) = 0
      do s = 1, nseg
         if (plan%sD(s)) then
            nt = int(plan%sA(s), int64)*int(plan%sA(s) + 1, int64)/2
         else
            nt = int(plan%sA(s), int64)*int(plan%sNB(s), int64)
         end if
         ! See bins_plan_build: the kernels index within a segment by a
         ! default integer, and bins are capped so that no segment needs more.
         if (nt > 2147483000_int64) error stop &
            "terco: bin-pair segment exceeds a default integer; widen the kernels' t"
         plan%sOff(s + 1) = plan%sOff(s) + nt
      end do
      ! The fork's parity masks: none on a plan this routine builds.
      allocate (plan%sM(max(nseg, 1)))
      plan%sM = 0
      plan%nseg = nseg
      plan%nwork = plan%sOff(nseg + 1)
      deallocate (pr_ka, pr_kb, pr_s0, pr_ns, pr_mode, pr_off, sl_pr, sl_ta, sl_kept, sl_nrun, sl_off, sl_area)
      deallocate (xbt0, xbtn, xt0, xtn, xtb1, xtb2, xpmax2, xsmin1, xtq, xqbin, xrlo, xrhi, xcnt, xoff, xla, xlb, xpb)
   end subroutine plan_device

   pure logical function admitted(b, ka, kb, smax_keep) result(res)
      !! Whether bin pair (ka, kb) passes the Schwarz cut and has pairs at all.
      type(pair_bins_t), intent(in) :: b
      integer, intent(in) :: ka, kb, smax_keep
      res = b%bin_s(ka) + b%bin_s(kb) <= smax_keep .and. b%bin_cnt(ka) > 0 .and. b%bin_cnt(kb) > 0
   end function admitted

   pure subroutine count_item(s, nbin, ntile, npr, nslot, bt0, btn, tn, tb1, tb2, pmax2, smin1, &
                              tq, qbin, rlo, rhi, pr_ka, pr_kb, sl_pr, sl_ta, thresh, sl_area, sl_kept, sl_nrun)
      !! Slot s, bra tile ta of bin pair p: the area and the kept tile pairs it
      !! adds to p's whole-or-cut test, and the segments it emits if p is cut.
      !$acc routine seq
      integer, value :: s
      integer, intent(in) :: nbin, ntile, npr, nslot
      integer, intent(in) :: bt0(nbin), btn(nbin), tn(ntile), tb1(ntile), tb2(ntile)
      integer, intent(in) :: pmax2(ntile), smin1(ntile), rlo(0:TILE_NLEV, ntile), rhi(0:TILE_NLEV, ntile)
      real(dp), intent(in) :: tq(ntile), qbin(nbin)
      real(dp), intent(in) :: thresh
      integer, intent(in) :: pr_ka(npr), pr_kb(npr), sl_pr(nslot), sl_ta(nslot)
      integer(int64), intent(inout) :: sl_area(nslot)
      integer, intent(inout) :: sl_kept(nslot), sl_nrun(nslot)
      integer :: p, ka, kb, ta, tb, t1, t2, tlast, nrun
      integer(int64) :: area
      logical :: same, keep, inrun

      p = sl_pr(s)
      ka = pr_ka(p)
      kb = pr_kb(p)
      ta = sl_ta(s)
      same = ka == kb
      call cand_d(nbin, ntile, bt0, btn, pmax2, smin1, tq, qbin, rlo, rhi, ta, kb, thresh, t1, t2)
      if (same) t2 = min(t2, ta)
      area = 0
      sl_kept(s) = 0
      do tb = t1, t2
         if (.not. kept_d(ntile, tq, tb1, tb2, rlo, rhi, ta, tb, thresh)) cycle
         sl_kept(s) = sl_kept(s) + 1
         if (same .and. tb == ta) then
            area = area + int(tn(ta), int64)*int(tn(ta) + 1, int64)/2
         else
            area = area + int(tn(ta), int64)*int(tn(tb), int64)
         end if
      end do
      sl_area(s) = area
      ! The cut's segments: one per run of consecutive kept ket tiles, and
      ! the diagonal triangle of a bin against itself.
      nrun = 0
      tlast = t2
      if (same .and. t2 == ta) tlast = ta - 1
      inrun = .false.
      do tb = t1, tlast
         keep = kept_d(ntile, tq, tb1, tb2, rlo, rhi, ta, tb, thresh)
         if (keep .and. .not. inrun) nrun = nrun + 1
         inrun = keep
      end do
      if (same .and. t2 == ta) then
         if (kept_d(ntile, tq, tb1, tb2, rlo, rhi, ta, ta, thresh)) nrun = nrun + 1
      end if
      sl_nrun(s) = nrun
   end subroutine count_item

   pure subroutine fill_item(s, nbin, ntile, npr, nslot, nseg, bt0, btn, t0, tn, tb1, tb2, pmax2, smin1, &
                             tq, qbin, rlo, rhi, pr_ka, pr_kb, pr_mode, sl_pr, sl_ta, sl_off, thresh, &
                             bin_cnt, bin_off, bin_la, bin_lb, bin_pbase, &
                             sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
      !! Slot s of a cut bin pair: its runs, then its diagonal triangle, from
      !! segment sl_off(s) + 1 on -- count_item's segments, in walk's order.
      !$acc routine seq
      integer, value :: s
      integer, intent(in) :: nbin, ntile, npr, nslot, nseg
      integer, intent(in) :: bt0(nbin), btn(nbin), t0(ntile), tn(ntile), tb1(ntile), tb2(ntile)
      integer, intent(in) :: pmax2(ntile), smin1(ntile), rlo(0:TILE_NLEV, ntile), rhi(0:TILE_NLEV, ntile)
      real(dp), intent(in) :: tq(ntile), qbin(nbin)
      real(dp), intent(in) :: thresh
      integer, intent(in) :: pr_ka(npr), pr_kb(npr), pr_mode(npr), sl_pr(nslot), sl_ta(nslot), sl_off(nslot)
      integer, intent(in) :: bin_cnt(nbin), bin_off(nbin), bin_la(nbin), bin_lb(nbin), bin_pbase(nbin)
      integer, intent(inout) :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg), sPA(nseg), sPB(nseg)
      integer, intent(inout) :: sSA(nseg), sSB(nseg), sLA(nseg), sLB(nseg), sLC(nseg), sLD(nseg)
      logical, intent(inout) :: sD(nseg)
      integer :: p, ka, kb, ta, tb, t1, t2, tlast, run0, run1, k
      logical :: same, keep, inrun

      p = sl_pr(s)
      if (pr_mode(p) /= 2) return
      ka = pr_ka(p)
      kb = pr_kb(p)
      ta = sl_ta(s)
      same = ka == kb
      k = sl_off(s)
      call cand_d(nbin, ntile, bt0, btn, pmax2, smin1, tq, qbin, rlo, rhi, ta, kb, thresh, t1, t2)
      if (same) t2 = min(t2, ta)
      ! Runs of consecutive kept ket tiles, tracked with an explicit in-run
      ! flag as count_item does. The host walk closed a run on an extra
      ! iteration past the end with keep forced false; on the device (nvfortran
      ! 26.5) that lost runs, leaving their segments unwritten.
      tlast = t2
      if (same .and. t2 == ta) tlast = ta - 1
      inrun = .false.
      run0 = 0
      run1 = 0
      do tb = t1, tlast
         keep = kept_d(ntile, tq, tb1, tb2, rlo, rhi, ta, tb, thresh)
         if (keep) then
            if (.not. inrun) run0 = tb
            run1 = tb
         else if (inrun) then
            k = k + 1
            call put_seg(k, nseg, nbin, ka, kb, t0(ta), tn(ta), t0(run0), t0(run1) + tn(run1) - t0(run0), &
                         .false., bin_cnt, bin_off, bin_la, bin_lb, bin_pbase, &
                         sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
         end if
         inrun = keep
      end do
      if (inrun) then
         k = k + 1
         call put_seg(k, nseg, nbin, ka, kb, t0(ta), tn(ta), t0(run0), t0(run1) + tn(run1) - t0(run0), &
                      .false., bin_cnt, bin_off, bin_la, bin_lb, bin_pbase, &
                      sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
      end if
      if (same .and. t2 == ta) then
         if (kept_d(ntile, tq, tb1, tb2, rlo, rhi, ta, ta, thresh)) then
            k = k + 1
            call put_seg(k, nseg, nbin, ka, kb, t0(ta), tn(ta), t0(ta), tn(ta), .true., &
                         bin_cnt, bin_off, bin_la, bin_lb, bin_pbase, &
                         sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
         end if
      end if
   end subroutine fill_item

   pure subroutine whole_item(p, nbin, npr, nseg, pr_ka, pr_kb, pr_mode, pr_off, &
                              bin_cnt, bin_off, bin_la, bin_lb, bin_pbase, &
                              sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
      !! A bin pair launched whole: its one segment.
      !$acc routine seq
      integer, value :: p
      integer, intent(in) :: nbin, npr, nseg
      integer, intent(in) :: pr_ka(npr), pr_kb(npr), pr_mode(npr), pr_off(npr)
      integer, intent(in) :: bin_cnt(nbin), bin_off(nbin), bin_la(nbin), bin_lb(nbin), bin_pbase(nbin)
      integer, intent(inout) :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg), sPA(nseg), sPB(nseg)
      integer, intent(inout) :: sSA(nseg), sSB(nseg), sLA(nseg), sLB(nseg), sLC(nseg), sLD(nseg)
      logical, intent(inout) :: sD(nseg)
      integer :: ka, kb
      if (pr_mode(p) /= 1) return
      ka = pr_ka(p)
      kb = pr_kb(p)
      call put_seg(pr_off(p) + 1, nseg, nbin, ka, kb, 0, bin_cnt(ka), 0, bin_cnt(kb), ka == kb, &
                   bin_cnt, bin_off, bin_la, bin_lb, bin_pbase, &
                   sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
   end subroutine whole_item

   pure subroutine put_seg(k, nseg, nbin, ka, kb, oa, na, ob, nb, diag, bin_cnt, bin_off, bin_la, bin_lb, &
                           bin_pbase, sA, sNB, sOA, sOB, sPA, sPB, sSA, sSB, sD, sLA, sLB, sLC, sLD)
      !! Segment k: pairs oa+1..oa+na of bin ka against ob+1..ob+nb of bin kb
      !! (walk's emit, less the work offset, which the host adds after).
      !$acc routine seq
      integer, intent(in) :: k, nseg, nbin, ka, kb, oa, na, ob, nb
      logical, intent(in) :: diag
      integer, intent(in) :: bin_cnt(nbin), bin_off(nbin), bin_la(nbin), bin_lb(nbin), bin_pbase(nbin)
      integer, intent(inout) :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg), sPA(nseg), sPB(nseg)
      integer, intent(inout) :: sSA(nseg), sSB(nseg), sLA(nseg), sLB(nseg), sLC(nseg), sLD(nseg)
      logical, intent(inout) :: sD(nseg)
      sA(k) = na
      sNB(k) = nb
      sOA(k) = bin_off(ka) + oa
      sOB(k) = bin_off(kb) + ob
      sPA(k) = bin_pbase(ka) + oa
      sPB(k) = bin_pbase(kb) + ob
      sSA(k) = bin_cnt(ka)
      sSB(k) = bin_cnt(kb)
      sD(k) = diag
      sLA(k) = bin_la(ka)
      sLB(k) = bin_lb(ka)
      sLC(k) = bin_la(kb)
      sLD(k) = bin_lb(kb)
   end subroutine put_seg

   pure subroutine cand_d(nbin, ntile, bt0, btn, pmax2, smin1, tq, qbin, rlo, rhi, ta, kb, thresh, t1, t2)
      !! `candidates`, on explicit arrays for the device.
      !$acc routine seq
      integer, intent(in) :: nbin, ntile, ta, kb
      integer, intent(in) :: bt0(nbin), btn(nbin), pmax2(ntile), smin1(ntile)
      integer, intent(in) :: rlo(0:TILE_NLEV, ntile), rhi(0:TILE_NLEV, ntile)
      real(dp), intent(in) :: tq(ntile), qbin(nbin)
      real(dp), intent(in) :: thresh
      integer, intent(out) :: t1, t2
      integer :: m, lo, hi, a, z, mid

      t1 = bt0(kb) + 1
      t2 = bt0(kb) + btn(kb)
      if (tq(ta)*qbin(kb) <= thresh) then
         t2 = t1 - 1
         return
      end if
      m = cut_level_d(thresh/(tq(ta)*qbin(kb)))
      if (m > TILE_NLEV) return
      lo = rlo(m, ta)
      hi = rhi(m, ta)
      if (lo > hi) then
         t2 = t1 - 1
         return
      end if
      a = t1
      z = t2 + 1
      do while (a < z)
         mid = (a + z)/2
         if (pmax2(mid) >= lo) then
            z = mid
         else
            a = mid + 1
         end if
      end do
      t1 = a
      a = bt0(kb)
      z = t2
      do while (a < z)
         mid = (a + z + 1)/2
         if (smin1(mid) <= hi) then
            a = mid
         else
            z = mid - 1
         end if
      end do
      t2 = a
   end subroutine cand_d

   pure logical function kept_d(ntile, tq, tb1, tb2, rlo, rhi, ta, tb, thresh) result(res)
      !! `tile_pair_kept`, on explicit arrays for the device.
      !$acc routine seq
      integer, intent(in) :: ntile, ta, tb
      integer, intent(in) :: tb1(ntile), tb2(ntile), rlo(0:TILE_NLEV, ntile), rhi(0:TILE_NLEV, ntile)
      real(dp), intent(in) :: tq(ntile)
      real(dp), intent(in) :: thresh
      integer :: m
      res = .false.
      if (tq(ta)*tq(tb) <= thresh) return
      m = cut_level_d(thresh/(tq(ta)*tq(tb)))
      res = .true.
      if (m > TILE_NLEV) return
      res = rlo(m, ta) <= tb2(tb) .and. rhi(m, ta) >= tb1(tb)
   end function kept_d

   pure integer function cut_level_d(cut) result(m)
      !! `cut_level`, as a device routine.
      !$acc routine seq
      real(dp), intent(in) :: cut
      if (cut >= 1.0_dp) then
         m = 0
      else if (cut <= 0.0_dp) then
         m = TILE_NLEV + 1
      else
         m = min(TILE_NLEV + 1, int(floor(-log10(cut))) + 1)
      end if
   end function cut_level_d

end module trc_bin_tiles
