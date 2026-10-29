!
! McMurchie-Davidson exchange (and Coulomb), per class, in do concurrent.
!
! For the g classes above TERCO_RYS_G_LTOT (the only ones scripts/gen_md_k.py
! generates), fock_bins skips the class and this pass adds it into the same
! folded matrix. The per-window algebra is the generator's: E matrices, R per
! primitive quartet, S = PQ Hk^T without storing PQ, the integral rows and
! their digestion. What lives here is the ORDER OF THE WORK.
!
! GROUPS. The shell pairs of the MD J object are grouped by (l_a, l_b) type and
! primitive count, so every window of a group pair has one shape, and each
! group is sorted by its Schwarz bound, largest first.
!
! STAIRCASE. A bra block of up to MDG_BLOCK sorted pairs, with largest bound
! Q_1, can only meet ket pairs with 4 Q_1 Q_cd max|D| > thresh: a PREFIX of
! the sorted ket list, found by bisection. A block with an empty prefix ends
! the bra loop, since every later block is smaller still: the early exit of
! Ufimtsev and Martinez's exchange, here at the granularity of a window. Inside
! a window the HGP kernels' per-quartet screen still zeroes what the coarse
! cut lets through.
!
! The fock_bins segments this replaces are split by Schwarz decade and by
! primitive count, and the first prototype ran one window per segment: 1462
! windows per Fock build on w20 cc-pVTZ, each far too small to fill the GPU
! (nsys, V100).
!
! REFERENCES (inspiration; the formulation and code here are terco's own)
!   L. E. McMurchie, E. R. Davidson, J. Comput. Phys. 26, 218 (1978).
!   I. S. Ufimtsev, T. J. Martinez, J. Chem. Theory Comput. 4, 222 (2008);
!     5, 1004 (2009).
!   E. Palethorpe, G. M. J. Barca, J. Chem. Theory Comput. 21, 9388 (2025).
!   E. Palethorpe, R. Stocks, G. M. J. Barca, ChemRxiv (2026),
!     doi:10.26434/chemrxiv.15004739/v2.
!   A. Asadchev, E. F. Valeev, J. Phys. Chem. A 127, 10889 (2023);
!     J. Chem. Phys. 160, 244109 (2024).
!
module trc_mdg
   use pic_types, only: dp, int64, int_index
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use pic_profiler, only: profiler_start, profiler_stop
   use pic_sorting, only: sort_index
   use trc_mdj, only: trc_mdj_t
   use trc_mdg_kernels, only: mdg_class_count, mdg_class_key, mdg_sizes, mdg_fused
   use trc_mdg_h, only: mdg_h, mdg_hsph, mdg_dsph, mdg_kcart
   use trc_api, only: trc_basis_t
   implicit none
   private
   public :: trc_mdg_t, trc_mdg_scales_t, mdg_routes

   integer, parameter :: CLASS_RADIX = 11
   !! trc_binkernel's class key radix.
   integer, parameter :: MDG_BLOCK = 256
   integer, parameter :: N_SIZES = 8
   !! entries of mdg_sizes
   integer, parameter :: N_WINDOW = 6
   !! entries of a window descriptor
   !! bra pairs per staircase block
   integer(int64), parameter :: MDG_BUDGET = 2_int64**26
   !! doubles in one window's buffers together (512 MB)

   type :: trc_mdg_scales_t
      real(dp) :: dmax = 0.0_dp
      !! largest shell-block density element, set by mdg_fock
      real(dp) :: thresh = 0.0_dp
      real(dp) :: pcut = 0.0_dp
      real(dp) :: jfac = 0.0_dp
      real(dp) :: kfac = 0.0_dp
      integer :: rank = 0
      integer :: nranks = 1
      !! MPI ranks: each takes every nranks-th window, and the caller sums
      !! the ranks' matrices, as for the bins
   end type trc_mdg_scales_t

   type :: trc_mdg_t
      logical :: built = .false.
      integer :: ngrp = 0
      integer, allocatable :: gr_ty(:), gr_k(:), gr_off(:), gr_n(:)
      !! per group: type index in the MD J object, primitive count, offset into gl, size
      integer, allocatable :: gl(:)
      !! shell-pair indices, group by group, each group by Schwarz bound descending
      real(dp), allocatable :: gq(:)
      !! their bounds, host only
      integer(int64), allocatable :: gh_off(:)
      !! per group: first element of its E matrices in gh, less one (-1: not formed)
      real(dp), allocatable :: gh(:)
      !! E matrices of every group a routed class reads, h(nab, np*k, pair), in
      !! each group's sorted order: geometry only, formed once. Per window they
      !! were 960 of 2010 ms per Fock build on w5 cc-pVQZ (V100).
      logical :: sph = .false.
      !! Exchange in spherical components: the caller's basis is spherical.
      !! The density goes in as D_sph = S^T D_cart S per shell block and the
      !! result comes back as S K_sph S^T, which the caller's fold T^T . T
      !! returns to exactly K_sph (S = T (T^T T)^-1). At (gg) that is 81
      !! component pairs where Cartesian is 225, in kernels B and C and the digestion.
      integer :: nsh = 0, naos = 0
      integer, allocatable :: sh_l(:), sh_ao(:), sh_ao_sph(:)
      !! the working basis's shells, for the block transforms
      real(dp), allocatable :: wrb(:), wsm(:)
      !! Window buffers, kept on the device across windows and builds and grown
      !! only when a group pair needs more: made per group pair, the device
      !! allocations were ~5 of 20.5 s per Fock build on w20 cc-pVQZ.
      logical :: stats = .false.
      !! TRC_MDG_STATS=1: count windows and launched and surviving quartets per
      !! build, and log them with per-class times. The survivor count is a host loop.
      integer(int64) :: st_win = 0, st_quart = 0, st_pass = 0
      integer(int64) :: nwin = 0
      !! windows met so far this build, on every rank alike: the rank split
      real(dp) :: st_cpass = 0.0_dp, st_cquart = 0.0_dp
      !! the current class's share, for the per-class line
   contains
      procedure :: build => mdg_build
      procedure :: fock => mdg_fock
      procedure :: release => mdg_release
   end type trc_mdg_t

contains

   logical function mdg_routes(key) result(res)
      !! Whether this class has an MD exchange kernel: the build generates them
      !! for exactly the g classes above TERCO_RYS_G_LTOT, which have no Rys one.
      integer, intent(in) :: key
      integer :: sz(N_SIZES)
      sz = mdg_sizes(key)
      res = sz(1) /= 0
   end function mdg_routes

   subroutine mdg_build(this, md, wb, sph)
      !! Groups of the MD J object's shell pairs, sorted, on the device, and
      !! their E matrices.
      class(trc_mdg_t), intent(inout) :: this
      type(trc_mdj_t), intent(in) :: md
      type(trc_basis_t), intent(in) :: wb
      !! the basis md was built over (the split one on a general contraction)
      logical, intent(in) :: sph
      !! the CALLER's basis is spherical
      integer :: t, s, n, k, g, m, kmax
      integer(int64), allocatable :: key(:)
      integer(int_index), allocatable :: sort_order(:)
      integer, allocatable :: bucket_pair(:)

      call this%release()
      kmax = 0
      if (md%nsp > 0) kmax = maxval(md%sp_nk(1:md%nsp))
      ! Count the groups: one per (type, primitive-count bucket) that has pairs.
      ! Buckets 1, 2, 3-4, 5-8, ...: a group pads its pairs to its largest
      ! count (at most 2x), and gets far fewer group pairs, each one window or
      ! more and three launches a window. Split by exact count, w20 cc-pVQZ
      ! had ~2000 windows per build whatever the block size, ~5 s of launches.
      this%ngrp = 0
      do t = 1, md%ntype
         do k = 1, kbucket(kmax)
            if (any(kbucket_v(md%sp_nk(md%ty_s0(t) + 1:md%ty_s0(t) + md%ty_ns(t))) == k)) then
              this%ngrp = this%ngrp + 1
            end if
         end do
      end do
      allocate (this%gr_ty(this%ngrp), this%gr_k(this%ngrp), this%gr_off(this%ngrp), this%gr_n(this%ngrp))
      allocate (this%gl(max(md%nsp, 1)), this%gq(max(md%nsp, 1)), bucket_pair(max(md%nsp, 1)))
      allocate (key(max(md%nsp, 1)), sort_order(max(md%nsp, 1)))
      g = 0
      m = 0
      do t = 1, md%ntype
         do k = 1, kbucket(kmax)
            n = 0
            do s = md%ty_s0(t) + 1, md%ty_s0(t) + md%ty_ns(t)
               if (kbucket(md%sp_nk(s)) /= k) cycle
               n = n + 1
               bucket_pair(n) = s
               ! Ascending integer key on -log10(Q): the largest bound first.
               key(n) = int(1.0e6_dp*min(60.0_dp, -log10(max(md%sp_q(s), 1.0e-60_dp))), int64)
            end do
            if (n == 0) cycle
            call sort_index(key(1:n), sort_order(1:n))
            g = g + 1
            this%gr_ty(g) = t
            this%gr_k(g) = maxval(md%sp_nk(bucket_pair(1:n)))
            this%gr_off(g) = m
            this%gr_n(g) = n
            do s = 1, n
               this%gl(m + s) = bucket_pair(sort_order(s))
               this%gq(m + s) = md%sp_q(bucket_pair(sort_order(s)))
            end do
            m = m + n
         end do
      end do
      deallocate (key, sort_order, bucket_pair)
      this%sph = sph
      this%nsh = wb%nshell
      this%naos = wb%nao_sph
      allocate (this%sh_l(this%nsh), this%sh_ao(this%nsh), this%sh_ao_sph(this%nsh))
      this%sh_l = wb%sh_l
      this%sh_ao = wb%sh_ao
      this%sh_ao_sph = wb%sh_ao_sph
      !$acc enter data copyin(this%gl, this%sh_l, this%sh_ao, this%sh_ao_sph)
      call form_e_matrices(this, md)
      this%built = .true.
   end subroutine mdg_build

   subroutine form_e_matrices(this, md)
      !! E matrices of the groups whose type a generated, routed class reads.
      type(trc_mdg_t), intent(inout) :: this
      type(trc_mdj_t), intent(in) :: md
      integer :: g, c, key, t1, t2, la, lb, ncol
      integer(int64) :: total
      logical, allocatable :: need(:)
      real(dp), allocatable :: hc(:)

      allocate (need(md%ntype), this%gh_off(this%ngrp))
      need = .false.
      do c = 1, mdg_class_count()
         key = mdg_class_key(c)
         if (.not. mdg_routes(key)) cycle
         call class_types(md, key, t1, t2)
         if (t1 > 0) need(t1) = .true.
         if (t2 > 0) need(t2) = .true.
      end do
      total = 0
      do g = 1, this%ngrp
         this%gh_off(g) = -1
         if (.not. need(this%gr_ty(g))) cycle
         la = md%ty_la(this%gr_ty(g))
         lb = md%ty_lb(this%gr_ty(g))
         this%gh_off(g) = total
         total = total + int(ncomp(this, la)*ncomp(this, lb), int64)*nherm(la + lb)*this%gr_k(g)*this%gr_n(g)
      end do
      allocate (this%gh(max(total, 1_int64)))
      !$acc enter data create(this%gh)
      do g = 1, this%ngrp
         if (this%gh_off(g) < 0) cycle
         la = md%ty_la(this%gr_ty(g))
         lb = md%ty_lb(this%gr_ty(g))
         ncol = nherm(la + lb)*this%gr_k(g)*this%gr_n(g)
         if (this%sph) then
            ! Cartesian first, into scratch, then the rows to spherical.
            allocate (hc(int(ncart(la)*ncart(lb), int64)*ncol))
            !$acc enter data create(hc)
            call mdg_h(la, lb, md%nhp, md%mp_z, md%mp_px, md%mp_py, md%mp_pz, md%mp_ax, md%mp_ay, md%mp_az, &
                       md%mp_bx, md%mp_by, md%mp_bz, md%mp_c, md%nsp, md%sp_ao1, md%sp_ao2, md%sp_as1, &
                       md%sp_as2, md%sp_k0, md%sp_nk, md%sp_sh1, md%sp_sh2, md%sp_q, size(this%gl), this%gl, &
                       this%gr_off(g), this%gr_n(g), this%gr_k(g), hc)
            call mdg_hsph(la, lb, ncol, hc, this%gh(this%gh_off(g) + 1))
            !$acc exit data delete(hc)
            deallocate (hc)
         else
            call mdg_h(la, lb, md%nhp, md%mp_z, md%mp_px, md%mp_py, md%mp_pz, md%mp_ax, md%mp_ay, md%mp_az, &
                       md%mp_bx, md%mp_by, md%mp_bz, md%mp_c, md%nsp, md%sp_ao1, md%sp_ao2, md%sp_as1, &
                       md%sp_as2, md%sp_k0, md%sp_nk, md%sp_sh1, md%sp_sh2, md%sp_q, size(this%gl), this%gl, &
                       this%gr_off(g), this%gr_n(g), this%gr_k(g), this%gh(this%gh_off(g) + 1))
         end if
      end do
      deallocate (need)
   end subroutine form_e_matrices

   pure integer function ncart(l) result(n)
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)/2
   end function ncart

   pure integer function ncomp(this, l) result(n)
      !! Components of a shell in the representation the exchange runs in.
      type(trc_mdg_t), intent(in) :: this
      integer, intent(in) :: l
      if (this%sph) then
         n = 2*l + 1
      else
         n = ncart(l)
      end if
   end function ncomp

   pure integer function nherm(l) result(n)
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)*(l + 3)/6
   end function nherm

   subroutine mdg_fock(this, md, dsh, dmat, jmat, sc)
      !! Every routed class's contribution into jmat (fock_bins' folded
      !! layout, ndens densities). dsh, dmat and jmat are on the device.
      class(trc_mdg_t), intent(inout) :: this
      type(trc_mdj_t), intent(in) :: md
      real(dp), intent(inout) :: dsh(:, :)
      real(dp), intent(in)    :: dmat(:, :, :)
      real(dp), intent(inout) :: jmat(:, :, :)
      type(trc_mdg_scales_t), intent(in) :: sc
      type(trc_mdg_scales_t) :: scd
      type(timer_type) :: tm
      integer :: c, key, t1, t2, g1, g2

      !$acc update self(dsh)
      scd = sc
      scd%dmax = maxval(dsh)
      this%stats = stats_on()
      this%st_win = 0
      this%nwin = 0
      this%st_quart = 0
      this%st_pass = 0
      if (this%sph) then
         call fock_spherical(this, md, dsh, dmat, jmat, scd)
         call report_stats(this)
         return
      end if
      do c = 1, mdg_class_count()
         key = mdg_class_key(c)
         if (.not. mdg_routes(key)) cycle
         call class_types(md, key, t1, t2)
         if (t1 == 0 .or. t2 == 0) cycle
         call profiler_start("mdg: "//class_label(key), nvtx_only=.true.)
         if (this%stats) then
            !$acc wait
            call tm%start()
            this%st_cpass = real(this%st_pass, dp)
            this%st_cquart = real(this%st_quart, dp)
         end if
         do g1 = 1, this%ngrp
            if (this%gr_ty(g1) /= t1) cycle
            do g2 = 1, this%ngrp
               if (this%gr_ty(g2) /= t2) cycle
               ! One type against itself: each unordered pair of groups once,
               ! a group against itself as a full square at half weight.
               if (t1 == t2 .and. g2 > g1) cycle
               call group_pair(this, md, [key, g1, g2], dsh, dmat, jmat, scd)
            end do
         end do
         if (this%stats) call class_line(this, key, tm)
         !$acc wait
         call profiler_stop("mdg: "//class_label(key))
      end do
      call report_stats(this)
   end subroutine mdg_fock

   function class_label(key) result(lab)
      !! "(gs|ps)" for a class key.
      integer, intent(in) :: key
      character(len=7) :: lab
      character(len=*), parameter :: L = "spdfghi"
      lab = "(" // L(key/(CLASS_RADIX**3) + 1:key/(CLASS_RADIX**3) + 1) &
            // L(mod(key/(CLASS_RADIX**2), CLASS_RADIX) + 1:mod(key/(CLASS_RADIX**2), CLASS_RADIX) + 1) // "|" &
            // L(mod(key/CLASS_RADIX, CLASS_RADIX) + 1:mod(key/CLASS_RADIX, CLASS_RADIX) + 1) &
            // L(mod(key, CLASS_RADIX) + 1:mod(key, CLASS_RADIX) + 1) // ")"
   end function class_label

   pure integer function kbucket(k) result(b)
      !! 1 -> 1, 2 -> 2, 3-4 -> 3, 5-8 -> 4, ...
      integer, intent(in) :: k
      b = 1
      do while (2**(b - 1) < k)
         b = b + 1
      end do
   end function kbucket

   pure function kbucket_v(k) result(b)
      integer, intent(in) :: k(:)
      integer :: b(size(k))
      integer :: i
      do i = 1, size(k)
         b(i) = kbucket(k(i))
      end do
   end function kbucket_v

   integer function block_size() result(res)
      !! Bra pairs per staircase block: MDG_BLOCK unless TRC_MDG_BLOCK says.
      character(len=16) :: e
      integer :: ios
      e = " "
      res = MDG_BLOCK
      call get_environment_variable("TRC_MDG_BLOCK", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) res
         if (ios /= 0 .or. res < 1) res = MDG_BLOCK
      end if
   end function block_size

   logical function stats_on() result(res)
      character(len=8) :: e
      e = " "
      call get_environment_variable("TRC_MDG_STATS", e)
      res = trim(e) == "1"
   end function stats_on

   subroutine class_line(this, key, tm)
      !! Stats: one class's wall time (device synchronised) and screen pass rate.
      type(trc_mdg_t), intent(in) :: this
      integer, intent(in) :: key
      type(timer_type), intent(inout) :: tm
      character(len=256) :: line
      real(dp) :: t, q, p
      !$acc wait
      t = tm%get_elapsed_time()
      q = real(this%st_quart, dp) - this%st_cquart
      p = real(this%st_pass, dp) - this%st_cpass
      write (line, "(a,4i1,a,f8.3,a,es9.2,a,f6.1,a)") "  [mdg] class ", key/(CLASS_RADIX**3), &
         mod(key/(CLASS_RADIX**2), CLASS_RADIX), mod(key/CLASS_RADIX, CLASS_RADIX), mod(key, CLASS_RADIX), &
         " s ", t, "  quartets ", q, "  pass ", 100.0_dp*p/max(q, 1.0_dp), "%"
      call global_logger%info(trim(line))
   end subroutine class_line

   subroutine report_stats(this)
      type(trc_mdg_t), intent(in) :: this
      character(len=256) :: line
      if (.not. this%stats) return
      write (line, "(a,i0,a,i0,a,i0,a,f6.2,a)") "  [mdg] windows ", this%st_win, "  quartets ", &
         this%st_quart, "  pass screen ", this%st_pass, "  (", &
         100.0_dp*real(this%st_pass, dp)/max(real(this%st_quart, dp), 1.0_dp), "%)"
      call global_logger%info(trim(line))
   end subroutine report_stats

   subroutine fock_spherical(this, md, dsh, dmat, jmat, sc)
      !! The classes in spherical components: D_sph = S^T D S in, the windows on
      !! (ndens, nao_sph, nao_sph) arrays, S K_sph S^T added back to jmat.
      type(trc_mdg_t), intent(inout) :: this
      type(trc_mdj_t), intent(in) :: md
      real(dp), intent(inout) :: dsh(:, :)
      real(dp), intent(in)    :: dmat(:, :, :)
      real(dp), intent(inout) :: jmat(:, :, :)
      type(trc_mdg_scales_t), intent(in) :: sc
      real(dp), allocatable :: ds(:, :, :), ks(:, :, :)
      type(timer_type) :: tm
      integer :: nd, c, key, t1, t2, g1, g2
      nd = size(dmat, 1)
      allocate (ds(nd, this%naos, this%naos), ks(nd, this%naos, this%naos))
      !$acc enter data create(ds, ks)
      call mdg_dsph(this%nsh, this%sh_l, this%sh_ao, this%sh_ao_sph, nd, size(dmat, 2), this%naos, dmat, ds)
      call zero3(ks)
      do c = 1, mdg_class_count()
         key = mdg_class_key(c)
         if (.not. mdg_routes(key)) cycle
         call class_types(md, key, t1, t2)
         if (t1 == 0 .or. t2 == 0) cycle
         call profiler_start("mdg: "//class_label(key), nvtx_only=.true.)
         if (this%stats) then
            !$acc wait
            call tm%start()
            this%st_cpass = real(this%st_pass, dp)
            this%st_cquart = real(this%st_quart, dp)
         end if
         do g1 = 1, this%ngrp
            if (this%gr_ty(g1) /= t1) cycle
            do g2 = 1, this%ngrp
               if (this%gr_ty(g2) /= t2) cycle
               if (t1 == t2 .and. g2 > g1) cycle
               call group_pair(this, md, [key, g1, g2], dsh, ds, ks, sc)
            end do
         end do
         if (this%stats) call class_line(this, key, tm)
         !$acc wait
         call profiler_stop("mdg: "//class_label(key))
      end do
      call mdg_kcart(this%nsh, this%sh_l, this%sh_ao, this%sh_ao_sph, nd, size(jmat, 2), this%naos, ks, jmat)
      !$acc exit data delete(ds, ks)
      deallocate (ds, ks)
   end subroutine fock_spherical

   subroutine zero3(a)
      real(dp), intent(inout) :: a(:, :, :)
      call zero_flat(size(a), a)
   end subroutine zero3

   subroutine zero_flat(n, a)
      integer,  intent(in)    :: n
      real(dp), intent(inout) :: a(n)
      integer :: i
      do concurrent(i = 1:n)
         a(i) = 0.0_dp
      end do
   end subroutine zero_flat

   subroutine group_pair(this, md, job, dsh, dmat, jmat, sc)
      !! The staircase over one bra group and one ket group.
      type(trc_mdg_t), intent(inout) :: this
      type(trc_mdj_t), intent(in) :: md
      integer, intent(in) :: job(3)
      !! class key, bra group, ket group
      real(dp), intent(inout) :: dsh(:, :)
      real(dp), intent(in)    :: dmat(:, :, :)
      real(dp), intent(inout) :: jmat(:, :, :)
      type(trc_mdg_scales_t), intent(in) :: sc
      integer :: sz(N_SIZES)
      integer :: ka, kb, ib, ni, nimax, jmax, j0, nj, njmax, key, g1, g2
      integer(int64) :: per
      real(dp) :: wdiag
      integer(int64) :: hb0, hk0, perb, perk

      key = job(1)
      g1 = job(2)
      g2 = job(3)
      sz = mdg_sizes(key)
      ! Component-pair counts in the representation the windows run in.
      if (this%sph) then
         sz(4) = sz(7)
         sz(5) = sz(8)
      end if
      ka = this%gr_k(g1)
      kb = this%gr_k(g2)
      wdiag = 1.0_dp
      if (g1 == g2) wdiag = 0.5_dp
      nimax = min(block_size(), this%gr_n(g1))
      ! Per ket pair of a window: its share of R and S for a full bra block.
      per = int(nimax, int64)*ka*kb*sz(3) + int(nimax, int64)*sz(1)*ka*sz(5)
      njmax = int(max(1_int64, min(int(this%gr_n(g2), int64), MDG_BUDGET/per)))
      call grow(this%wrb, int(sz(3), int64)*nimax*ka*njmax*kb)
      call grow(this%wsm, int(sz(1), int64)*ka*nimax*sz(5)*njmax)
      hb0 = this%gh_off(g1)
      hk0 = this%gh_off(g2)
      perb = int(sz(4), int64)*sz(1)*ka
      perk = int(sz(5), int64)*sz(2)*kb
      do ib = 0, this%gr_n(g1) - 1, nimax
         ni = min(nimax, this%gr_n(g1) - ib)
         jmax = ket_prefix(this, g2, 4.0_dp*this%gq(this%gr_off(g1) + ib + 1)*sc%dmax, sc%thresh)
         ! Sorted: no later bra block can do better than this one.
         if (jmax == 0) exit
         do j0 = 0, jmax - 1, njmax
            nj = min(njmax, jmax - j0)
            ! Every rank walks the same windows in the same order (same
            ! density, same sorted groups), so the count picks a disjoint share.
            this%nwin = this%nwin + 1
            if (mod(this%nwin - 1, int(sc%nranks, int64)) /= sc%rank) cycle
            if (this%stats) call count_window(this, md, [this%gr_off(g1) + ib, ni, ka, this%gr_off(g2) + j0, nj, kb], &
                                              sz, dsh, sc%thresh)
            call mdg_fused(key, this%sph, md%nhp, md%mp_z, md%mp_px, md%mp_py, md%mp_pz, md%mp_ax, md%mp_ay, &
                        md%mp_az, md%mp_bx, md%mp_by, md%mp_bz, md%mp_c, md%nsp, md%sp_ao1, md%sp_ao2, &
                        md%sp_as1, md%sp_as2, md%sp_k0, md%sp_nk, md%sp_sh1, md%sp_sh2, md%sp_q, &
                        size(this%gl), this%gl, this%gr_off(g1) + ib, ni, ka, this%gr_off(g2) + j0, nj, kb, &
                        sc%thresh, sc%pcut, sc%jfac, sc%kfac, wdiag, size(dsh, 1), dsh, &
                        size(dmat, 2), size(dmat, 1), dmat, jmat, this%gh(hb0 + ib*perb + 1), &
                        this%gh(hk0 + j0*perk + 1), this%wrb, this%wsm)
         end do
      end do
   end subroutine group_pair

   subroutine grow(a, n)
      !! a, on the device, holds at least n elements; contents not kept.
      real(dp), allocatable, intent(inout) :: a(:)
      integer(int64), intent(in) :: n
      if (allocated(a)) then
         if (size(a, kind=int64) >= n) return
         !$acc exit data delete(a)
         deallocate (a)
      end if
      allocate (a(max(n, 1_int64)))
      !$acc enter data create(a)
   end subroutine grow

   subroutine count_window(this, md, w, sz, dsh, thresh)
      !! Stats: a window's quartets, how many pass the per-quartet screen
      !! (same test as the R kernel).
      type(trc_mdg_t), intent(inout) :: this
      type(trc_mdj_t), intent(in) :: md
      integer, intent(in) :: w(N_WINDOW), sz(N_SIZES)
      !! ib0, ni, ka, jb0, nj, kb
      real(dp), intent(in) :: dsh(:, :)
      real(dp), intent(in) :: thresh
      integer :: i, j, sa, sc2, si, sj, sk, sl
      real(dp) :: q
      this%st_win = this%st_win + 1
      this%st_quart = this%st_quart + int(w(2), int64)*w(5)
      do j = 1, w(5)
         sc2 = this%gl(w(4) + j)
         sk = md%sp_sh1(sc2)
         sl = md%sp_sh2(sc2)
         do i = 1, w(2)
            sa = this%gl(w(1) + i)
            si = md%sp_sh1(sa)
            sj = md%sp_sh2(sa)
            q = md%sp_q(sa)*md%sp_q(sc2)
            if (q <= thresh) cycle
            if (q*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), dsh(si, sk), dsh(si, sl), &
                      dsh(sj, sk), dsh(sj, sl)) <= thresh) cycle
            this%st_pass = this%st_pass + 1
         end do
      end do
   end subroutine count_window

   integer function ket_prefix(this, g, qd, thresh) result(n)
      !! How many of group g's sorted pairs have qd * Q_cd > thresh: bisection.
      type(trc_mdg_t), intent(in) :: this
      integer, intent(in) :: g
      real(dp), intent(in) :: qd, thresh
      integer :: lo, hi, mid
      lo = 0
      hi = this%gr_n(g)
      do while (lo < hi)
         mid = (lo + hi + 1)/2
         if (qd*this%gq(this%gr_off(g) + mid) > thresh) then
            lo = mid
         else
            hi = mid - 1
         end if
      end do
      n = lo
   end function ket_prefix

   subroutine class_types(md, key, t1, t2)
      !! The MD J object's type indices of the class's bra and ket pairs.
      type(trc_mdj_t), intent(in) :: md
      integer, intent(in) :: key
      integer, intent(out) :: t1, t2
      integer :: la, lb, lc, ld, t
      la = key/(CLASS_RADIX**3)
      lb = mod(key/(CLASS_RADIX**2), CLASS_RADIX)
      lc = mod(key/CLASS_RADIX, CLASS_RADIX)
      ld = mod(key, CLASS_RADIX)
      t1 = 0
      t2 = 0
      do t = 1, md%ntype
         if (md%ty_la(t) == la .and. md%ty_lb(t) == lb) t1 = t
         if (md%ty_la(t) == lc .and. md%ty_lb(t) == ld) t2 = t
      end do
   end subroutine class_types

   subroutine mdg_release(this)
      class(trc_mdg_t), intent(inout) :: this
      if (this%built) then
         !$acc exit data delete(this%gl, this%gh, this%sh_l, this%sh_ao, this%sh_ao_sph)
      end if
      if (allocated(this%gr_ty)) deallocate (this%gr_ty, this%gr_k, this%gr_off, this%gr_n)
      if (allocated(this%gl)) deallocate (this%gl, this%gq)
      if (allocated(this%gh)) deallocate (this%gh, this%gh_off)
      if (allocated(this%sh_l)) deallocate (this%sh_l, this%sh_ao, this%sh_ao_sph)
      if (allocated(this%wrb)) then
         !$acc exit data delete(this%wrb, this%wsm)
         deallocate (this%wrb, this%wsm)
      end if
      this%ngrp = 0
      this%built = .false.
   end subroutine mdg_release

end module trc_mdg
