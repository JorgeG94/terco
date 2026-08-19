!
! HGP over binned shell pairs: the quartet is formed on the device.
!
! One launch covers the cross product of two bins.  A thread turns its own
! linear index into (ab, cd) pair indices, reads two shell pairs, and evaluates
! the quartet.  Nothing enumerates quartets on the host and nothing stores
! them.
!
! Screening happened when the bin pair was admitted (trc_bins), so there is
! no screening conditional here and every thread has real work.
!
! The digestion is fused and symmetry-aware, so the integral tensor is never
! materialised either -- which is what removed the last large host allocation.
!
!
! [1] J. L. Galvez Vallejo, G. M. J. Barca and M. S. Gordon,
!     "High-performance GPU-accelerated evaluation of electron repulsion
!     integrals", Mol. Phys. (2022) e2112987, doi:10.1080/00268976.2022.2112987
!
module trc_binkernel
   use trc_boys, only: dp, boys_eval, BOYS_MMAX
   use trc_tables, only: LMAX
   use trc_cart, only: NCUM, cidx, cnx, cny, cnz, cll, cdir, cdn1, cdn2, cf2, &
                         ncum_of, ncart_of
   use trc_bins, only: pair_bins_t, SMAX
#ifdef TRC_PERCLASS
   use trc_pc_kernels, only: pc_dispatch
#endif
   implicit none
   private

   public :: fock_bins

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp


   !
   ! The launch plan: which bin pairs survive, how their quartets are numbered,
   ! and in what order the class-specialised kernels cover them.
   !
   ! This used to be rebuilt inside every `fock_bins` call -- a double loop
   ! over live bin pairs, ten allocations, a sort by class, a host-to-device
   ! copy, and the matching frees -- once per SCF ITERATION, although it
   ! depends only on the bins and the threshold. The density never enters it
   ! unless bin-level density screening is on, and that is the one case where
   ! a plan has to be rebuilt when D changes.
   !
   ! Promoting it to an object is the same move cuEST makes with
   ! `cuestOEIntPlanCreate`: setup that depends on the basis and not on the
   ! data is done once and handed back.
   !
   type, public :: trc_plan_t
      integer :: nseg = 0, nlaunch = 0
      integer(kind=8) :: nwork = 0
      integer, allocatable :: sA(:), sB(:), sOA(:), sOB(:), sNB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(kind=8), allocatable :: sOff(:)
      logical :: on_device = .false.
   contains
      procedure :: build   => plan_build
      procedure :: release => plan_release
   end type trc_plan_t


contains

   !
   ! LOCAL COPIES of trc_hgp's comp/binom/powi.
   !
   ! NVHPC does not inline an `!$acc routine seq` helper across a module
   ! boundary, and these are called from the innermost loop of the HRR --
   ! `binom` and `powi` six times per (b', d') term, which is up to sixteen
   ! terms per component.  Cross-module they were real calls in the hottest
   ! loop in the digestion.
   !
   ! Duplication, and only justified because the ablation put the digestion at
   ! 31% of runtime.  trc_hgp keeps the originals for the standalone HGP path
   ! and check_hgp still validates those.
   !
   pure subroutine comp(l, ic, px, py, pz)
      !$acc routine seq
      integer, intent(in)  :: l, ic
      integer, intent(out) :: px, py, pz
      integer :: n, lx, ly
      n = 0
      do lx = l, 0, -1
         do ly = l - lx, 0, -1
            if (n == ic) then
               px = lx; py = ly; pz = l - lx - ly
               return
            end if
            n = n + 1
         end do
      end do
      px = 0; py = 0; pz = 0
   end subroutine comp

   pure real(dp) function binom(n, k)
      !$acc routine seq
      integer, intent(in) :: n, k
      integer :: i
      real(dp) :: r
      r = 1.0_dp
      do i = 1, k
         r = r*real(n - k + i, dp)/real(i, dp)
      end do
      binom = r
   end function binom

   pure real(dp) function powi(x, n)
      !$acc routine seq
      real(dp), intent(in) :: x
      integer,  intent(in) :: n
      integer :: i
      powi = 1.0_dp
      do i = 1, n
         powi = powi*x
      end do
   end function powi

   !
   ! Whole Fock build: walk admitted bin pairs on the host, launch one kernel
   ! each.  The host loop is over BINS -- tens of thousands of iterations at
   ! most -- not over quartets.
   !
   subroutine fock_bins(b, nbas, npp, nao, sh_l, ao_off, thresh, use_dens, &
                        jfac, kfac, nosym, dsh, &
                        pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                        ndens, dmat, jmat, kmat, nlaunch, nwork, nkept)
      type(pair_bins_t), intent(in) :: b
      integer,  intent(in)    :: nbas, npp, nao
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      real(dp), intent(in)    :: thresh
      !! Coulomb and exchange scalings; 1 and 1 is Hartree-Fock.
      real(dp), intent(in)    :: jfac, kfac
      !! Contract a density that is NOT symmetric, through the enumerated
      !! kernel. Costs what folding bought; needed by the response Hessian.
      logical,  intent(in)    :: nosym
      !! max |D| over each shell-pair block; the per-quartet density screen.
      !! Pass it filled with huge() to disable density screening entirely.
      real(dp), intent(in)    :: dsh(nbas, nbas)
      logical,  intent(in)    :: use_dens
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)
      integer,  intent(out)   :: nlaunch
      integer(kind=8), intent(out) :: nwork
      !! Diagnostic: how many of `nwork` actually survive the in-kernel
      !! Schwarz and density tests.  Costs a separate cheap pass, so it is
      !! optional and off in the timed path.
      integer(kind=8), intent(out), optional :: nkept

      integer :: ia, ib, ka, kb, smax_keep, nA, nB, nseg, is
      integer(kind=8) :: nt
      integer, allocatable :: sA(:), sB(:), sOA(:), sOB(:), sNB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(kind=8), allocatable :: sOff(:)

      smax_keep = int(-log10(thresh))

      ! --- pass 1: count admitted bin pairs ---
      nseg = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            if (use_dens) then
               if (dens_reject(b, ka, kb, thresh)) cycle
            end if
            nseg = nseg + 1
         end do
      end do

      allocate (sA(nseg), sB(nseg), sOA(nseg), sOB(nseg), sNB(nseg), sD(nseg))
      allocate (sLA(nseg), sLB(nseg), sLC(nseg), sLD(nseg))
      allocate (sOff(nseg + 1))

      ! --- pass 2: fill descriptors and prefix-sum the work ---
      !
      ! ONE launch over every admitted bin pair, not one launch each.  Binning
      ! at (type, size) granularity fragments the work badly -- 484 launches
      ! averaging 2200 items apiece measured 12x slower than a single large
      ! launch.  The screening granularity is worth keeping; the launch
      ! granularity is not.  A prefix sum lets a thread find its own segment.
      !
      is = 0; sOff(1) = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         nA = b%bin_cnt(ka)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            nB = b%bin_cnt(kb)
            if (nA == 0 .or. nB == 0) cycle
            if (use_dens) then
               if (dens_reject(b, ka, kb, thresh)) cycle
            end if
            is = is + 1
            sA(is) = nA; sNB(is) = nB
            sOA(is) = b%bin_off(ka); sOB(is) = b%bin_off(kb)
            sD(is) = (ka == kb)
            sLA(is) = b%bin_la(ka); sLB(is) = b%bin_lb(ka)
            sLC(is) = b%bin_la(kb); sLD(is) = b%bin_lb(kb)
            if (ka == kb) then
               nt = int(nA, 8)*int(nA + 1, 8)/2
            else
               nt = int(nA, 8)*int(nB, 8)
            end if
            !
            ! one_bin_item takes its within-segment index as a default integer,
            ! so a segment must fit in one.  At 75 waters the largest bin holds
            ! 3617 pairs and the largest segment is 6.5e6, against a 2.1e9
            ! limit -- about 65,000 pairs in a single bin before it bites,
            ! which is a system some 20x larger.  Guarded rather than widened
            ! because widening the index costs registers in the hot kernel and
            ! the register budget is already the binding constraint.
            !
            if (nt > 2147483000_8) error stop &
               "terco: bin-pair segment exceeds a default integer; widen t in one_bin_item"
            sOff(is + 1) = sOff(is) + nt
         end do
      end do

      nwork = sOff(nseg + 1)
      nlaunch = 1
      if (nwork == 0) return

      if (nosym) then
         ! Every per-class kernel is folded, so a general density goes through
         ! the enumerated shared kernel instead.
         !$acc enter data copyin(sA, sNB, sOA, sOB, sD, sOff)
         call fock_all_nosym(nseg, nwork, sA, sNB, sOA, sOB, sD, sOff, &
                             b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, jfac, kfac, dsh, &
                             nbas, npp, nao, sh_l, ao_off, &
                             pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                             ndens, dmat, jmat, kmat)
         nlaunch = 1
         !$acc exit data delete(sA, sNB, sOA, sOB, sD, sOff)
         deallocate (sA, sB, sOA, sOB, sNB, sD, sOff, sLA, sLB, sLC, sLD)
         return
      end if

#ifdef TRC_PERCLASS
      !
      ! Per-class launches.  The descriptors are sorted by class so each class
      ! occupies a contiguous run of segments, then one kernel per class covers
      ! its own run.  Every thread in a launch then executes code compiled for
      ! exactly its (la,lb,lc,ld) -- no select case, and ptxas budgets registers
      ! per class instead of for the worst class across all of them.
      !
      block
         integer, allocatable :: ord(:), ckey(:)
         integer :: a2, b2, t2, c0, c1, nl
         logical :: l2
         integer(kind=8) :: o2
         allocate (ord(nseg), ckey(nseg))
         do a2 = 1, nseg
            ckey(a2) = ((sLA(a2)*(LMAX + 1) + sLB(a2))*(LMAX + 1) &
                        + sLC(a2))*(LMAX + 1) + sLD(a2)
            ord(a2) = a2
         end do
         ! insertion sort on the class key; nseg is small (hundreds)
         do a2 = 2, nseg
            t2 = ord(a2)
            b2 = a2 - 1
            do while (b2 >= 1)
               if (ckey(ord(b2)) <= ckey(t2)) exit
               ord(b2 + 1) = ord(b2); b2 = b2 - 1
            end do
            ord(b2 + 1) = t2
         end do
         call permute_segments(nseg, ord, sA, sNB, sOA, sOB, sD, sOff, &
                               sLA, sLB, sLC, sLD)

         !$acc enter data copyin(sA, sNB, sOA, sOB, sD, sOff)
         if (present(nkept)) call count_kept(nseg, nwork, sNB, sOA, sOB, sD, &
                                             sOff, b%npair, b%sp_i, b%sp_j, &
                                             b%sp_q, thresh, dsh, nbas, nkept)
         nl = 0
         c0 = 1
         do while (c0 <= nseg)
            c1 = c0
            do while (c1 < nseg)
               if (((sLA(c1 + 1)*(LMAX + 1) + sLB(c1 + 1))*(LMAX + 1) &
                    + sLC(c1 + 1))*(LMAX + 1) + sLD(c1 + 1) /= &
                   ((sLA(c0)*(LMAX + 1) + sLB(c0))*(LMAX + 1) &
                    + sLC(c0))*(LMAX + 1) + sLD(c0)) exit
               c1 = c1 + 1
            end do
            nl = nl + 1
            call pc_dispatch(((sLA(c0)*(LMAX + 1) + sLB(c0))*(LMAX + 1) &
                              + sLC(c0))*(LMAX + 1) + sLD(c0), &
                             c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sD, &
                             b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                             pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, ndens, dmat, jmat)
            c0 = c1 + 1
         end do
         nlaunch = nl
         !$acc exit data delete(sA, sNB, sOA, sOB, sD, sOff)
         deallocate (ord, ckey)
      end block
      deallocate (sA, sB, sOA, sOB, sNB, sD, sOff, sLA, sLB, sLC, sLD)
      return
#endif

      !$acc enter data copyin(sA, sNB, sOA, sOB, sD, sOff)
      call fock_all(nseg, nwork, sA, sNB, sOA, sOB, sD, sOff, &
                    b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                    pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, ndens, dmat, jmat, kmat)
      !$acc exit data delete(sA, sNB, sOA, sOB, sD, sOff)

      deallocate (sA, sB, sOA, sOB, sNB, sD, sOff, sLA, sLB, sLC, sLD)
   end subroutine fock_bins

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
   pure logical function dens_reject(b, ka, kb, thresh)
      type(pair_bins_t), intent(in) :: b
      integer,  intent(in) :: ka, kb
      real(dp), intent(in) :: thresh
      real(dp) :: qq, dd
      dens_reject = .false.
      if (.not. allocated(b%bin_dm)) return
      ! bound on Q_ab Q_cd from the two size buckets, back off one decade each
      ! because the bucket is a ceiling
      qq = 10.0_dp**(-real(max(0, b%bin_s(ka) - 1), dp)) &
           *10.0_dp**(-real(max(0, b%bin_s(kb) - 1), dp))
      dd = 4.0_dp*max(b%bin_dm(ka), b%bin_dm(kb))
      dens_reject = (qq*dd <= thresh)
   end function dens_reject

   !
   ! Count the quartets that survive screening, using the same index
   ! arithmetic and the same two tests as the kernel.  Diagnostic only: it
   ! exists so the effect of density screening is a measured number rather
   ! than an inference from wall-clock.
   !
   subroutine count_kept(nseg, nwork, sNB, sOA, sOB, sD, sOff, &
                         npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, nkept)
      integer,  intent(in) :: nseg, npair, nbas
      integer(kind=8), intent(in) :: nwork
      integer,  intent(in) :: sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in) :: sD(nseg)
      integer(kind=8), intent(in) :: sOff(nseg + 1)
      integer,  intent(in) :: sp_i(npair), sp_j(npair)
      real(dp), intent(in) :: sp_q(npair), thresh, dsh(nbas, nbas)
      integer(kind=8), intent(out) :: nkept

      integer(kind=8) :: gt, nk
      integer :: lo, hi, mid, seg, t, iab, icd, si, sj, sk, sl
      real(dp) :: qcut

      nk = 0
      !
      ! THE ONE PLACE A DIRECTIVE CARRIES A LOOP.
      !
      ! `do concurrent` gained `reduce()` in Fortran 2023, which is exactly
      ! what this wants -- but gfortran did not take locality specifiers until
      ! 15 and `reduce` later still, and terco is meant to build on the host
      ! with whatever a CI runner has. A reduction is the one thing the
      ! portable subset of the language cannot say, so it is said here instead
      ! and nowhere else.
      !
      ! Worth keeping in proportion: this counts kept quartets for a work
      ! estimate. It is not on the path any integral takes.
      !
      !$acc parallel loop reduction(+:nk) default(present) &
      !$acc   private(lo, hi, mid, seg, t, iab, icd, si, sj, sk, sl, qcut)
      do gt = 1, nwork
         lo = 1; hi = nseg
         do while (lo < hi)
            mid = (lo + hi + 1)/2
            if (sOff(mid) < gt) then
               lo = mid
            else
               hi = mid - 1
            end if
         end do
         seg = lo
         t = int(gt - sOff(seg))
         if (sD(seg)) then
            iab = int((1.0_dp + sqrt(1.0_dp + 8.0_dp*real(t - 1, dp)))/2.0_dp)
            if (iab*(iab + 1)/2 > t - 1) iab = iab - 1
            icd = (t - 1) - iab*(iab + 1)/2 + 1
            iab = iab + 1
         else
            iab = (t - 1)/sNB(seg) + 1
            icd = t - (iab - 1)*sNB(seg)
         end if
         qcut = sp_q(sOA(seg) + iab)*sp_q(sOB(seg) + icd)
         if (qcut > thresh) then
            si = sp_i(sOA(seg) + iab); sj = sp_j(sOA(seg) + iab)
            sk = sp_i(sOB(seg) + icd); sl = sp_j(sOB(seg) + icd)
            if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                         dsh(si, sk), dsh(si, sl), &
                         dsh(sj, sk), dsh(sj, sl)) > thresh) nk = nk + 1
         end if
      end do
      nkept = nk
   end subroutine count_kept

   !> Apply a permutation to the segment descriptor arrays, in place.
   subroutine permute_segments(n, ord, sA, sNB, sOA, sOB, sD, sOff, sLA, sLB, sLC, sLD)
      integer, intent(in) :: n, ord(n)
      integer, intent(inout) :: sA(n), sNB(n), sOA(n), sOB(n)
      integer, intent(inout) :: sLA(n), sLB(n), sLC(n), sLD(n)
      logical, intent(inout) :: sD(n)
      integer(kind=8), intent(inout) :: sOff(n + 1)
      integer :: k
      integer, allocatable :: t(:)
      logical, allocatable :: tl(:)
      integer(kind=8), allocatable :: to(:)
      allocate (t(n), tl(n), to(n + 1))
      t = sA(ord);  sA = t
      t = sNB(ord); sNB = t
      t = sOA(ord); sOA = t
      t = sOB(ord); sOB = t
      t = sLA(ord); sLA = t
      t = sLB(ord); sLB = t
      t = sLC(ord); sLC = t
      t = sLD(ord); sLD = t
      tl = sD(ord); sD = tl
      ! rebuild the prefix sum in the new order
      to(1) = 0
      do k = 1, n
         if (sD(k)) then
            to(k + 1) = to(k) + int(sA(k), 8)*int(sA(k) + 1, 8)/2
         else
            to(k + 1) = to(k) + int(sA(k), 8)*int(sNB(k), 8)
         end if
      end do
      sOff = to
      deallocate (t, tl, to)
   end subroutine permute_segments

   !
   ! The single launch.  Each thread locates its segment by binary search over
   ! the prefix sums -- about log2(nseg) steps, and uniform within a warp since
   ! neighbouring threads almost always share a segment.
   !
   subroutine fock_all(nseg, nwork, sA, sNB, sOA, sOB, sD, sOff, &
                       npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                       pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, ndens, dmat, jmat, kmat)
      integer,  intent(in)    :: nseg, npair, nbas, npp, nao
      integer(kind=8), intent(in) :: nwork
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer(kind=8), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

      integer(kind=8) :: gt
      integer :: lo, hi, mid, seg, t

      !
      ! `do concurrent`, and only that.
      !
      ! This was written twice for a while -- an OpenACC `parallel loop`
      ! behind -DTRC_ACC_LOOP -- to compare the two offload models on
      ! identical work.  The comparison has been made; carrying a directive
      ! version afterwards is a second thing to keep correct for no result.
      !
      ! OpenACC remains for the two jobs `do concurrent` genuinely cannot do:
      ! `routine seq` on the device helpers, and the atomic updates in the
      ! digestion, where a quartet scatters into six overlapping blocks and
      ! the standard has no atomic to offer.
      !
      ! `local()` names every per-iteration scalar explicitly. A scalar
      ! assigned before it is read is privatised without one, so this is not
      ! load-bearing -- it is the difference between the compiler inferring
      ! the intent and the source stating it, and the second is what someone
      ! reading this in a year needs.
      !
      ! Fortran 2018, and gfortran took it in 15, which is this project's
      ! minimum for that reason among others.
      do concurrent(gt=1:nwork) local(lo, hi, mid, seg, t)
         lo = 1; hi = nseg
         do while (lo < hi)
            mid = (lo + hi + 1)/2
            if (sOff(mid) < gt) then
               lo = mid
            else
               hi = mid - 1
            end if
         end do
         seg = lo
         t = int(gt - sOff(seg))
         call one_bin_item(t, sD(seg), sA(seg), sNB(seg), sOA(seg), sOB(seg), &
                           npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                           pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                           ndens, dmat, jmat, kmat)
      end do
   end subroutine fock_all

   !> Whole-worklist driver for the general-density kernel.
   subroutine fock_all_nosym(nseg, nwork, sA, sNB, sOA, sOB, sD, sOff, &
                             npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, &
                             nbas, npp, nao, sh_l, ao_off, &
                             pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                             ndens, dmat, jmat, kmat)
      integer,  intent(in)    :: nseg, npair, nbas, npp, nao, ndens
      integer(kind=8), intent(in) :: nwork
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer(kind=8), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh, jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

      integer(kind=8) :: gt
      integer :: lo, hi, mid, seg, t

      do concurrent(gt=1:nwork)
         lo = 1; hi = nseg
         do while (lo < hi)
            mid = (lo + hi + 1)/2
            if (sOff(mid) < gt) then
               lo = mid
            else
               hi = mid - 1
            end if
         end do
         seg = lo
         t = int(gt - sOff(seg))
         call one_bin_item_nosym(t, sD(seg), sA(seg), sNB(seg), sOA(seg), sOB(seg), &
                                 npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, dsh, &
                                 nbas, npp, nao, sh_l, ao_off, &
                                 pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                 ndens, dmat, jmat, kmat)
      end do
   end subroutine fock_all_nosym

   subroutine fock_one_bin(nt, diag, nA, nB, offA, offB, npair, sp_i, sp_j, &
                           sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                           pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                           ndens, dmat, jmat, kmat)
      integer,  intent(in)    :: nt, nA, nB, offA, offB, npair, nbas, npp, nao
      logical,  intent(in)    :: diag
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

      integer :: t

      do concurrent(t=1:nt)
         call one_bin_item(t, diag, nA, nB, offA, offB, npair, sp_i, sp_j, &
                           sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                           pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                           ndens, dmat, jmat, kmat)
      end do
   end subroutine fock_one_bin

   !
   ! General-density item routine: no assumption that D is symmetric.
   !
   ! The folded six-atomic digestion collapses the eight permutations using
   ! D_mn = D_nm and the caller then symmetrises the result. For the coupled-
   ! perturbed Hessian's response density both steps are wrong: J survives,
   ! being contracted against an index pair the integral is symmetric in, but
   ! K loses its antisymmetric part.
   !
   ! This is the ENUMERATED digestion, which was always general -- it reads
   ! D_mn and D_nm as distinct and writes J and K into separate matrices. It
   ! lives in the source behind the TRC_FOCK6 `#else` and stopped being
   ! compiled when folding became the default; this is that branch, forced,
   ! as its own routine.
   !
   ! SINGLE DENSITY by construction: mqc's `build_fock_direct_nosym` takes
   ! one, and the enumerated digestion nests deeply enough inside its
   ! permutation guards that threading a batch loop through it is more risk
   ! than batching would be worth. Ranks match the batched kernels so one
   ! driver serves both; the density index is a literal 1.
   !
   ! Slower than the folded path by about what folding bought. The Hessian is
   ! built once per makeEFP; the SCF Fock a hundred times.
   !
   pure subroutine one_bin_item_nosym(t, diag, nA, nB, offA, offB, npair, sp_i, sp_j, &
                                sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                                pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                ndens, dmat, jmat, kmat)
      !$acc routine seq
      integer,  intent(in)    :: t, nA, nB, offA, offB, npair, nbas, npp, nao
      logical,  intent(in)    :: diag
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

      integer  :: iab, icd, si, sj, sk, sl
      real(dp) :: qcut
      integer  :: la, lb, lc, ld, lab, lcd, lt, nca, ncc
      integer  :: keyab, keycd, offab, offcd, nab, ncd
      integer  :: kp, kq, m, d, x, cur, nxt, lim, dn
      integer  :: ia, ic, xd1, xd2, xc1, xc2, xac, ai, dd
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc, acc
      real(dp) :: pa(3), qc(3), wp(3), wq(3), pq(3)
      real(dp) :: pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: f(0:BOYS_MMAX)
      real(dp) :: v(NCUM*NCUM, 0:1), g(NCUM*NCUM)
      integer  :: na_, nb_, nc_, nd_, ib_, id_, mu, nu, lam, sig, idx
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: px, py, pz, qx, qy, qz, ja, jc
      real(dp) :: cb, cdc, vv, jsink, sc
      integer  :: idens
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      integer  :: mui, nuj, lamk, sigl
      ! largest Cartesian block this LMAX allows: ((L+1)(L+2)/2)^4
      real(dp) :: vbuf((((LMAX + 1)*(LMAX + 2)/2)**4))
      integer  :: pij, pkl
      logical  :: dij, dkl, dpq

      ! --- linear index -> (ab, cd) within the bin cross product ---
      if (diag) then
         ! lower triangle including the diagonal, row-major by ab
         iab = int((1.0_dp + sqrt(1.0_dp + 8.0_dp*real(t - 1, dp)))/2.0_dp)
         if (iab*(iab + 1)/2 > t - 1) iab = iab - 1
         icd = (t - 1) - iab*(iab + 1)/2 + 1
         iab = iab + 1
      else
         iab = (t - 1)/nB + 1
         icd = t - (iab - 1)*nB
      end if

      ! exact Schwarz, per quartet.  The host bin-pair test is only
      ! decade-granular (`int(-log10 Q)` clamped to 9) and so admits up to two
      ! decades too much; at RNA3 that was 46% of all quartets.
      qcut = sp_q(offA + iab)*sp_q(offB + icd)
      if (qcut <= thresh) return

      si = sp_i(offA + iab); sj = sp_j(offA + iab)
      sk = sp_i(offB + icd); sl = sp_j(offB + icd)

      ! Density-weighted screen -- see the note in scripts/gen_perclass.py.
      if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                   dsh(si, sk), dsh(si, sl), &
                   dsh(sj, sk), dsh(sj, sl)) <= thresh) return

      !
      ! NO canonical swap between bra and ket.
      !
      ! An earlier version enforced pair_ab >= pair_cd by swapping the two
      ! pairs per thread.  That silently destroyed the property the bins exist
      ! to provide: bin A and bin B have different (l, K) types, so after a
      ! per-thread swap some threads in a warp hold (A|B) and others (B|A),
      ! with different loop bounds and different primitive counts.  Nsight
      ! measured 21.3 of 32 threads active, against 30.1 before the binning.
      !
      ! The swap is not needed.  Symmetry only requires each unordered pair of
      ! shell pairs once, and the bin loop already guarantees that: distinct
      ! bins contribute every (A_i, B_j) exactly once, and the diagonal bin
      ! walks its own triangle.  So "the two pairs differ" -- which is all the
      ! digestion needs -- is a property of the bin pair, not of the shell
      ! indices.
      !
      pij = si*(si - 1)/2 + sj
      pkl = sk*(sk - 1)/2 + sl

      la = sh_l(si); lb = sh_l(sj); lc = sh_l(sk); ld = sh_l(sl)
      lab = la + lb; lcd = lc + ld; lt = lab + lcd
      nca = ncum_of(lab); ncc = ncum_of(lcd)

      keyab = (si - 1)*nbas + sj
      keycd = (sk - 1)*nbas + sl
      offab = pp_off(keyab); nab = pp_n(keyab)
      offcd = pp_off(keycd); ncd = pp_n(keycd)

      jsink = 0.0_dp
      do x = 1, nca*ncc
         g(x) = 0.0_dp
      end do

      do kp = offab + 1, offab + nab
         zeta = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            zpe = zeta + eta
            rho = zeta*eta/zpe
            do d = 1, 3
               pq(d) = pp_r(kp, d) - pp_r(kq, d)
               wc = (zeta*pp_r(kp, d) + eta*pp_r(kq, d))/zpe
               pa(d) = pp_r(kp, d) - pp_ra(kp, d)
               qc(d) = pp_r(kq, d) - pp_ra(kq, d)
               wp(d) = wc - pp_r(kp, d)
               wq(d) = wc - pp_r(kq, d)
            end do
            tval = rho*(pq(1)*pq(1) + pq(2)*pq(2) + pq(3)*pq(3))
            call boys_eval(lt, tval, f)

            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
            rz = rho/zeta; re = rho/eta
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

#ifdef TRC_UNROLL_VRR
            !
            ! Unrolled VRR: every index a literal, no table loads, no branches.
            ! See scripts/gen_vrr.py.
            !
            ! Straight-line VRR, textually inside this procedure so `v` never
            ! crosses a call boundary.  See scripts/gen_vrr.py --mode include.
            pax = pa(1); pay = pa(2); paz = pa(3)
            qcx = qc(1); qcy = qc(2); qcz = qc(3)
            wpx = wp(1); wpy = wp(2); wpz = wp(3)
            wqx = wq(1); wqy = wq(2); wqz = wq(3)
#include "trc_vrr_cases.inc"
#else
            cur = 0
            v(1, cur) = pref*f(lt)
            do m = lt - 1, 0, -1
               nxt = 1 - cur
               lim = lt - m
               v(1, nxt) = pref*f(m)
               do ia = 2, nca
                  if (cll(ia) > lim) exit
                  dd = cdir(ia); xd1 = cdn1(ia); xd2 = cdn2(ia)
                  acc = pa(dd)*v(xd1, nxt) + wp(dd)*v(xd1, cur)
                  if (xd2 > 0) acc = acc + cf2(ia)*oo2z*(v(xd2, nxt) - rz*v(xd2, cur))
                  v(ia, nxt) = acc
               end do
               do ic = 2, ncc
                  if (cll(ic) > lim) exit
                  dd = cdir(ic); xc1 = cdn1(ic); xc2 = cdn2(ic)
                  do ia = 1, nca
                     if (cll(ia) + cll(ic) > lim) exit
                     x = ia + (ic - 1)*nca
                     acc = qc(dd)*v(ia + (xc1 - 1)*nca, nxt) &
                           + wq(dd)*v(ia + (xc1 - 1)*nca, cur)
                     if (xc2 > 0) acc = acc + cf2(ic)*oo2e* &
                        (v(ia + (xc2 - 1)*nca, nxt) - re*v(ia + (xc2 - 1)*nca, cur))
                     ai = 0; xac = 0
                     if (dd == 1) then
                        ai = cnx(ia)
                        if (ai > 0) xac = cidx(cnx(ia) - 1, cny(ia), cnz(ia))
                     else if (dd == 2) then
                        ai = cny(ia)
                        if (ai > 0) xac = cidx(cnx(ia), cny(ia) - 1, cnz(ia))
                     else
                        ai = cnz(ia)
                        if (ai > 0) xac = cidx(cnx(ia), cny(ia), cnz(ia) - 1)
                     end if
                     if (xac > 0) acc = acc + real(ai, dp)*oo2ze*v(xac + (xc1 - 1)*nca, cur)
                     v(x, nxt) = acc
                  end do
               end do
               cur = nxt
            end do

#endif
            do x = 1, nca*ncc
               g(x) = g(x) + v(x, cur)
            end do
         end do
      end do

      ! --- HRR and fused digestion, straight from g ---
      na_ = ncart_of(la); nb_ = ncart_of(lb)
      nc_ = ncart_of(lc); nd_ = ncart_of(ld)
      dij = (si /= sj); dkl = (sk /= sl)
      ! distinct bins are always distinct pairs; the diagonal bin's
      ! triangle makes them equal only on iab == icd
      dpq = (.not. diag) .or. (iab /= icd)

#ifdef TRC_UNROLL_HRR
      ! Straight-line HRR + digestion, textually here so `g` never crosses a
      ! call boundary.  See scripts/gen_hrr.py.
      mui = ao_off(si); nuj = ao_off(sj)
      lamk = ao_off(sk); sigl = ao_off(sl)
      abx = pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1)
      aby = pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2)
      abz = pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3)
      cdx = pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1)
      cdy = pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2)
      cdz = pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3)
#include "trc_hrr_cases.inc"
      !
      ! ONE copy of the scatter, looping over the unrolled dot products.
      !
      ! Emitting the digestion per component beside each dot product -- 81
      ! copies for (pp|pp) -- measured 17% SLOWER than the generic HRR (5.76 s
      ! against 4.91 s at 75 waters).  The atomics are free per the ablation,
      ! so replicating them is instruction-cache cost with no work removed.
      ! Unroll the arithmetic; loop the scatter.
      !
      idx = 0
      do id_ = 0, nd_ - 1
         sig = sigl + id_
      do ic = 0, nc_ - 1
         lam = lamk + ic
      do ib_ = 0, nb_ - 1
         nu = nuj + ib_
      do ia = 0, na_ - 1
         mu = mui + ia
         idx = idx + 1
         vv = vbuf(idx)
#if defined(TRC_NO_DIGEST)
         jsink = jsink + vv
#else
#if defined(TRC_NO_ATOMIC)
#define FG_ATOMIC
#else
#define FG_ATOMIC !$acc atomic update
#endif

         FG_ATOMIC
         jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, lam, sig)
         FG_ATOMIC
         kmat(1, mu, lam) = kmat(1, mu, lam) + vv*dmat(1, nu, sig)
         if (dij) then
            FG_ATOMIC
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, lam, sig)
            FG_ATOMIC
            kmat(1, nu, lam) = kmat(1, nu, lam) + vv*dmat(1, mu, sig)
         end if
         if (dkl) then
            FG_ATOMIC
            jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, sig, lam)
            FG_ATOMIC
            kmat(1, mu, sig) = kmat(1, mu, sig) + vv*dmat(1, nu, lam)
         end if
         if (dij .and. dkl) then
            FG_ATOMIC
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, sig, lam)
            FG_ATOMIC
            kmat(1, nu, sig) = kmat(1, nu, sig) + vv*dmat(1, mu, lam)
         end if
         if (dpq) then
            FG_ATOMIC
            jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, mu, nu)
            FG_ATOMIC
            kmat(1, lam, mu) = kmat(1, lam, mu) + vv*dmat(1, sig, nu)
            if (dkl) then
               FG_ATOMIC
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, mu, nu)
               FG_ATOMIC
               kmat(1, sig, mu) = kmat(1, sig, mu) + vv*dmat(1, lam, nu)
            end if
            if (dij) then
               FG_ATOMIC
               jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, nu, mu)
               FG_ATOMIC
               kmat(1, lam, nu) = kmat(1, lam, nu) + vv*dmat(1, sig, mu)
            end if
            if (dij .and. dkl) then
               FG_ATOMIC
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, nu, mu)
               FG_ATOMIC
               kmat(1, sig, nu) = kmat(1, sig, nu) + vv*dmat(1, lam, mu)
            end if
         end if
#endif
      end do
      end do
      end do
      end do
#else
      do id_ = 0, nd_ - 1
         call comp(ld, id_, dx, dy, dz)
         sig = ao_off(sl) + id_
      do ic = 0, nc_ - 1
         call comp(lc, ic, cx, cy, cz)
         lam = ao_off(sk) + ic
      do ib_ = 0, nb_ - 1
         call comp(lb, ib_, bx, by, bz)
         nu = ao_off(sj) + ib_
      do ia = 0, na_ - 1
         call comp(la, ia, ax, ay, az)
         mu = ao_off(si) + ia

         acc = 0.0_dp
         do pz = 0, bz
         do py = 0, by
         do px = 0, bx
            cb = binom(bx, px)*binom(by, py)*binom(bz, pz) &
                 *powi(pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1), bx - px) &
                 *powi(pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2), by - py) &
                 *powi(pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3), bz - pz)
            if (cb == 0.0_dp) cycle
            ja = cidx(ax + px, ay + py, az + pz)
            do qz = 0, dz
            do qy = 0, dy
            do qx = 0, dx
               cdc = binom(dx, qx)*binom(dy, qy)*binom(dz, qz) &
                     *powi(pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1), dx - qx) &
                     *powi(pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2), dy - qy) &
                     *powi(pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3), dz - qz)
               if (cdc == 0.0_dp) cycle
               jc = cidx(cx + qx, cy + qy, cz + qz)
               acc = acc + cb*cdc*g(ja + (jc - 1)*nca)
            end do
            end do
            end do
         end do
         end do
         end do
         vv = acc

!
! DIAGNOSTIC ABLATIONS -- both produce WRONG RESULTS on purpose.
!
!   -DTRC_NO_ATOMIC   same stores, no atomic directives.  Isolates the cost
!                       of the atomic instruction itself from the address
!                       arithmetic around it.  Races, so J and K are garbage.
!   -DTRC_NO_DIGEST   no scatter at all, just a scalar accumulation.
!                       Isolates the whole digestion, addresses included.
!
! Only ever for reading `ptxas` register counts and kernel timings.  Neither is
! a correctness path and check_binfock will fail on both, which is the point.
!
#if defined(TRC_NO_DIGEST)
         jsink = jsink + vv
#else
#if defined(TRC_NO_ATOMIC)
#define FG_ATOMIC
#else
#define FG_ATOMIC !$acc atomic update
#endif

         FG_ATOMIC
         jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, lam, sig)
         FG_ATOMIC
         kmat(1, mu, lam) = kmat(1, mu, lam) + vv*dmat(1, nu, sig)
         if (dij) then
            FG_ATOMIC
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, lam, sig)
            FG_ATOMIC
            kmat(1, nu, lam) = kmat(1, nu, lam) + vv*dmat(1, mu, sig)
         end if
         if (dkl) then
            FG_ATOMIC
            jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, sig, lam)
            FG_ATOMIC
            kmat(1, mu, sig) = kmat(1, mu, sig) + vv*dmat(1, nu, lam)
         end if
         if (dij .and. dkl) then
            FG_ATOMIC
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, sig, lam)
            FG_ATOMIC
            kmat(1, nu, sig) = kmat(1, nu, sig) + vv*dmat(1, mu, lam)
         end if
         if (dpq) then
            FG_ATOMIC
            jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, mu, nu)
            FG_ATOMIC
            kmat(1, lam, mu) = kmat(1, lam, mu) + vv*dmat(1, sig, nu)
            if (dkl) then
               FG_ATOMIC
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, mu, nu)
               FG_ATOMIC
               kmat(1, sig, mu) = kmat(1, sig, mu) + vv*dmat(1, lam, nu)
            end if
            if (dij) then
               FG_ATOMIC
               jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, nu, mu)
               FG_ATOMIC
               kmat(1, lam, nu) = kmat(1, lam, nu) + vv*dmat(1, sig, mu)
            end if
            if (dij .and. dkl) then
               FG_ATOMIC
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, nu, mu)
               FG_ATOMIC
               kmat(1, sig, nu) = kmat(1, sig, nu) + vv*dmat(1, lam, mu)
            end if
         end if
#endif
      end do
      end do
      end do
      end do
#endif
#if defined(TRC_NO_DIGEST)
      ! keep jsink live so the whole evaluation is not dead-code eliminated
      if (jsink /= 0.0_dp) then
         !$acc atomic update
         jmat(1, 1) = jmat(1, 1) + jsink*1.0e-300_dp
      end if
#endif
   end subroutine one_bin_item_nosym

   pure subroutine one_bin_item(t, diag, nA, nB, offA, offB, npair, sp_i, sp_j, &
                                sp_q, thresh, jfac, kfac, dsh, nbas, npp, nao, sh_l, ao_off, &
                                pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                ndens, dmat, jmat, kmat)
      !$acc routine seq
      integer,  intent(in)    :: t, nA, nB, offA, offB, npair, nbas, npp, nao
      logical,  intent(in)    :: diag
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

      integer  :: iab, icd, si, sj, sk, sl
      real(dp) :: qcut
      integer  :: la, lb, lc, ld, lab, lcd, lt, nca, ncc
      integer  :: keyab, keycd, offab, offcd, nab, ncd
      integer  :: kp, kq, m, d, x, cur, nxt, lim, dn
      integer  :: ia, ic, xd1, xd2, xc1, xc2, xac, ai, dd
      real(dp) :: zeta, eta, zpe, rho, tval, pref, wc, acc
      real(dp) :: pa(3), qc(3), wp(3), wq(3), pq(3)
      real(dp) :: pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: f(0:BOYS_MMAX)
      real(dp) :: v(NCUM*NCUM, 0:1), g(NCUM*NCUM)
      integer  :: na_, nb_, nc_, nd_, ib_, id_, mu, nu, lam, sig, idx
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: px, py, pz, qx, qy, qz, ja, jc
      real(dp) :: cb, cdc, vv, jsink, sc
      integer  :: idens
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      integer  :: mui, nuj, lamk, sigl
      ! largest Cartesian block this LMAX allows: ((L+1)(L+2)/2)^4
      real(dp) :: vbuf((((LMAX + 1)*(LMAX + 2)/2)**4))
      integer  :: pij, pkl
      logical  :: dij, dkl, dpq

      ! --- linear index -> (ab, cd) within the bin cross product ---
      if (diag) then
         ! lower triangle including the diagonal, row-major by ab
         iab = int((1.0_dp + sqrt(1.0_dp + 8.0_dp*real(t - 1, dp)))/2.0_dp)
         if (iab*(iab + 1)/2 > t - 1) iab = iab - 1
         icd = (t - 1) - iab*(iab + 1)/2 + 1
         iab = iab + 1
      else
         iab = (t - 1)/nB + 1
         icd = t - (iab - 1)*nB
      end if

      ! exact Schwarz, per quartet.  The host bin-pair test is only
      ! decade-granular (`int(-log10 Q)` clamped to 9) and so admits up to two
      ! decades too much; at RNA3 that was 46% of all quartets.
      qcut = sp_q(offA + iab)*sp_q(offB + icd)
      if (qcut <= thresh) return

      si = sp_i(offA + iab); sj = sp_j(offA + iab)
      sk = sp_i(offB + icd); sl = sp_j(offB + icd)

      ! Density-weighted screen -- see the note in scripts/gen_perclass.py.
      if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                   dsh(si, sk), dsh(si, sl), &
                   dsh(sj, sk), dsh(sj, sl)) <= thresh) return

      !
      ! NO canonical swap between bra and ket.
      !
      ! An earlier version enforced pair_ab >= pair_cd by swapping the two
      ! pairs per thread.  That silently destroyed the property the bins exist
      ! to provide: bin A and bin B have different (l, K) types, so after a
      ! per-thread swap some threads in a warp hold (A|B) and others (B|A),
      ! with different loop bounds and different primitive counts.  Nsight
      ! measured 21.3 of 32 threads active, against 30.1 before the binning.
      !
      ! The swap is not needed.  Symmetry only requires each unordered pair of
      ! shell pairs once, and the bin loop already guarantees that: distinct
      ! bins contribute every (A_i, B_j) exactly once, and the diagonal bin
      ! walks its own triangle.  So "the two pairs differ" -- which is all the
      ! digestion needs -- is a property of the bin pair, not of the shell
      ! indices.
      !
      pij = si*(si - 1)/2 + sj
      pkl = sk*(sk - 1)/2 + sl

      la = sh_l(si); lb = sh_l(sj); lc = sh_l(sk); ld = sh_l(sl)
      lab = la + lb; lcd = lc + ld; lt = lab + lcd
      nca = ncum_of(lab); ncc = ncum_of(lcd)

      keyab = (si - 1)*nbas + sj
      keycd = (sk - 1)*nbas + sl
      offab = pp_off(keyab); nab = pp_n(keyab)
      offcd = pp_off(keycd); ncd = pp_n(keycd)

      jsink = 0.0_dp
      do x = 1, nca*ncc
         g(x) = 0.0_dp
      end do

      do kp = offab + 1, offab + nab
         zeta = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            zpe = zeta + eta
            rho = zeta*eta/zpe
            do d = 1, 3
               pq(d) = pp_r(kp, d) - pp_r(kq, d)
               wc = (zeta*pp_r(kp, d) + eta*pp_r(kq, d))/zpe
               pa(d) = pp_r(kp, d) - pp_ra(kp, d)
               qc(d) = pp_r(kq, d) - pp_ra(kq, d)
               wp(d) = wc - pp_r(kp, d)
               wq(d) = wc - pp_r(kq, d)
            end do
            tval = rho*(pq(1)*pq(1) + pq(2)*pq(2) + pq(3)*pq(3))
            call boys_eval(lt, tval, f)

            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = 0.5_dp/zpe
            rz = rho/zeta; re = rho/eta
            pref = TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

#ifdef TRC_UNROLL_VRR
            !
            ! Unrolled VRR: every index a literal, no table loads, no branches.
            ! See scripts/gen_vrr.py.
            !
            ! Straight-line VRR, textually inside this procedure so `v` never
            ! crosses a call boundary.  See scripts/gen_vrr.py --mode include.
            pax = pa(1); pay = pa(2); paz = pa(3)
            qcx = qc(1); qcy = qc(2); qcz = qc(3)
            wpx = wp(1); wpy = wp(2); wpz = wp(3)
            wqx = wq(1); wqy = wq(2); wqz = wq(3)
#include "trc_vrr_cases.inc"
#else
            cur = 0
            v(1, cur) = pref*f(lt)
            do m = lt - 1, 0, -1
               nxt = 1 - cur
               lim = lt - m
               v(1, nxt) = pref*f(m)
               do ia = 2, nca
                  if (cll(ia) > lim) exit
                  dd = cdir(ia); xd1 = cdn1(ia); xd2 = cdn2(ia)
                  acc = pa(dd)*v(xd1, nxt) + wp(dd)*v(xd1, cur)
                  if (xd2 > 0) acc = acc + cf2(ia)*oo2z*(v(xd2, nxt) - rz*v(xd2, cur))
                  v(ia, nxt) = acc
               end do
               do ic = 2, ncc
                  if (cll(ic) > lim) exit
                  dd = cdir(ic); xc1 = cdn1(ic); xc2 = cdn2(ic)
                  do ia = 1, nca
                     if (cll(ia) + cll(ic) > lim) exit
                     x = ia + (ic - 1)*nca
                     acc = qc(dd)*v(ia + (xc1 - 1)*nca, nxt) &
                           + wq(dd)*v(ia + (xc1 - 1)*nca, cur)
                     if (xc2 > 0) acc = acc + cf2(ic)*oo2e* &
                        (v(ia + (xc2 - 1)*nca, nxt) - re*v(ia + (xc2 - 1)*nca, cur))
                     ai = 0; xac = 0
                     if (dd == 1) then
                        ai = cnx(ia)
                        if (ai > 0) xac = cidx(cnx(ia) - 1, cny(ia), cnz(ia))
                     else if (dd == 2) then
                        ai = cny(ia)
                        if (ai > 0) xac = cidx(cnx(ia), cny(ia) - 1, cnz(ia))
                     else
                        ai = cnz(ia)
                        if (ai > 0) xac = cidx(cnx(ia), cny(ia), cnz(ia) - 1)
                     end if
                     if (xac > 0) acc = acc + real(ai, dp)*oo2ze*v(xac + (xc1 - 1)*nca, cur)
                     v(x, nxt) = acc
                  end do
               end do
               cur = nxt
            end do

#endif
            do x = 1, nca*ncc
               g(x) = g(x) + v(x, cur)
            end do
         end do
      end do

      ! --- HRR and fused digestion, straight from g ---
      na_ = ncart_of(la); nb_ = ncart_of(lb)
      nc_ = ncart_of(lc); nd_ = ncart_of(ld)
      dij = (si /= sj); dkl = (sk /= sl)
      ! distinct bins are always distinct pairs; the diagonal bin's
      ! triangle makes them equal only on iab == icd
      dpq = (.not. diag) .or. (iab /= icd)

#ifdef TRC_UNROLL_HRR
      ! Straight-line HRR + digestion, textually here so `g` never crosses a
      ! call boundary.  See scripts/gen_hrr.py.
      mui = ao_off(si); nuj = ao_off(sj)
      lamk = ao_off(sk); sigl = ao_off(sl)
      abx = pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1)
      aby = pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2)
      abz = pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3)
      cdx = pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1)
      cdy = pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2)
      cdz = pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3)
#include "trc_hrr_cases.inc"
      !
      ! ONE copy of the scatter, looping over the unrolled dot products.
      !
      ! Emitting the digestion per component beside each dot product -- 81
      ! copies for (pp|pp) -- measured 17% SLOWER than the generic HRR (5.76 s
      ! against 4.91 s at 75 waters).  The atomics are free per the ablation,
      ! so replicating them is instruction-cache cost with no work removed.
      ! Unroll the arithmetic; loop the scatter.
      !
      idx = 0
      do id_ = 0, nd_ - 1
         sig = sigl + id_
      do ic = 0, nc_ - 1
         lam = lamk + ic
      do ib_ = 0, nb_ - 1
         nu = nuj + ib_
      do ia = 0, na_ - 1
         mu = mui + ia
         idx = idx + 1
         vv = vbuf(idx)
#if defined(TRC_FOCK6)
         !
         ! SIX atomic updates into a single G = 2J - K, the form used in [1].
         !
         ! The eight-permutation version below writes sixteen times per
         ! component.  Ten of those are transposes of the other six, so
         ! accumulating into symmetric storage and symmetrising once at the end
         ! cuts the digestion's global read-modify-writes by 2.7x.  The earlier
         ! ablation said atomics cost nothing, but it removed only the atomic
         ! DIRECTIVE and kept the stores -- it measured the instruction, not the
         ! traffic.  The traffic is the part that matters.
         !
         ! The 0.5 factors are the standard degeneracy scales: a quartet with
         ! i==j, k==l or (ij)==(kl) stands for fewer than eight ordered tuples,
         ! and the halvings compensate exactly.  This is where prefactor-based
         ! digestions usually go wrong, which is why it is validated against the
         ! same libfint reference as the enumerated form.
         !
         sc = vv
         if (.not. dij) sc = sc*0.5_dp
         if (.not. dkl) sc = sc*0.5_dp
         if (.not. dpq) sc = sc*0.5_dp
         do idens = 1, ndens
         !$acc atomic update
         jmat(idens, mu, nu) = jmat(idens, mu, nu) &
                               + 4.0_dp*jfac*sc*dmat(idens, lam, sig)
         !$acc atomic update
         jmat(idens, lam, sig) = jmat(idens, lam, sig) + 4.0_dp*jfac*sc*dmat(idens, mu, nu)
         !$acc atomic update
         jmat(idens, mu, lam) = jmat(idens, mu, lam) - kfac*sc*dmat(idens, nu, sig)
         !$acc atomic update
         jmat(idens, mu, sig) = jmat(idens, mu, sig) - kfac*sc*dmat(idens, nu, lam)
         !$acc atomic update
         jmat(idens, nu, lam) = jmat(idens, nu, lam) - kfac*sc*dmat(idens, mu, sig)
         !$acc atomic update
         jmat(idens, nu, sig) = jmat(idens, nu, sig) - kfac*sc*dmat(idens, mu, lam)
      end do
#else
#if defined(TRC_NO_DIGEST)
         jsink = jsink + vv
#else
#if defined(TRC_NO_ATOMIC)
#define FG_ATOMIC
#else
#define FG_ATOMIC !$acc atomic update
#endif

         FG_ATOMIC
         jmat(mu, nu) = jmat(mu, nu) + vv*dmat(lam, sig)
         FG_ATOMIC
         kmat(mu, lam) = kmat(mu, lam) + vv*dmat(nu, sig)
         if (dij) then
            FG_ATOMIC
            jmat(nu, mu) = jmat(nu, mu) + vv*dmat(lam, sig)
            FG_ATOMIC
            kmat(nu, lam) = kmat(nu, lam) + vv*dmat(mu, sig)
         end if
         if (dkl) then
            FG_ATOMIC
            jmat(mu, nu) = jmat(mu, nu) + vv*dmat(sig, lam)
            FG_ATOMIC
            kmat(mu, sig) = kmat(mu, sig) + vv*dmat(nu, lam)
         end if
         if (dij .and. dkl) then
            FG_ATOMIC
            jmat(nu, mu) = jmat(nu, mu) + vv*dmat(sig, lam)
            FG_ATOMIC
            kmat(nu, sig) = kmat(nu, sig) + vv*dmat(mu, lam)
         end if
         if (dpq) then
            FG_ATOMIC
            jmat(lam, sig) = jmat(lam, sig) + vv*dmat(mu, nu)
            FG_ATOMIC
            kmat(lam, mu) = kmat(lam, mu) + vv*dmat(sig, nu)
            if (dkl) then
               FG_ATOMIC
               jmat(sig, lam) = jmat(sig, lam) + vv*dmat(mu, nu)
               FG_ATOMIC
               kmat(sig, mu) = kmat(sig, mu) + vv*dmat(lam, nu)
            end if
            if (dij) then
               FG_ATOMIC
               jmat(lam, sig) = jmat(lam, sig) + vv*dmat(nu, mu)
               FG_ATOMIC
               kmat(lam, nu) = kmat(lam, nu) + vv*dmat(sig, mu)
            end if
            if (dij .and. dkl) then
               FG_ATOMIC
               jmat(sig, lam) = jmat(sig, lam) + vv*dmat(nu, mu)
               FG_ATOMIC
               kmat(sig, nu) = kmat(sig, nu) + vv*dmat(lam, mu)
            end if
         end if
#endif
#endif
      end do
      end do
      end do
      end do
#else
      do id_ = 0, nd_ - 1
         call comp(ld, id_, dx, dy, dz)
         sig = ao_off(sl) + id_
      do ic = 0, nc_ - 1
         call comp(lc, ic, cx, cy, cz)
         lam = ao_off(sk) + ic
      do ib_ = 0, nb_ - 1
         call comp(lb, ib_, bx, by, bz)
         nu = ao_off(sj) + ib_
      do ia = 0, na_ - 1
         call comp(la, ia, ax, ay, az)
         mu = ao_off(si) + ia

         acc = 0.0_dp
         do pz = 0, bz
         do py = 0, by
         do px = 0, bx
            cb = binom(bx, px)*binom(by, py)*binom(bz, pz) &
                 *powi(pp_ra(offab + 1, 1) - pp_rb(offab + 1, 1), bx - px) &
                 *powi(pp_ra(offab + 1, 2) - pp_rb(offab + 1, 2), by - py) &
                 *powi(pp_ra(offab + 1, 3) - pp_rb(offab + 1, 3), bz - pz)
            if (cb == 0.0_dp) cycle
            ja = cidx(ax + px, ay + py, az + pz)
            do qz = 0, dz
            do qy = 0, dy
            do qx = 0, dx
               cdc = binom(dx, qx)*binom(dy, qy)*binom(dz, qz) &
                     *powi(pp_ra(offcd + 1, 1) - pp_rb(offcd + 1, 1), dx - qx) &
                     *powi(pp_ra(offcd + 1, 2) - pp_rb(offcd + 1, 2), dy - qy) &
                     *powi(pp_ra(offcd + 1, 3) - pp_rb(offcd + 1, 3), dz - qz)
               if (cdc == 0.0_dp) cycle
               jc = cidx(cx + qx, cy + qy, cz + qz)
               acc = acc + cb*cdc*g(ja + (jc - 1)*nca)
            end do
            end do
            end do
         end do
         end do
         end do
         vv = acc

!
! DIAGNOSTIC ABLATIONS -- both produce WRONG RESULTS on purpose.
!
!   -DTRC_NO_ATOMIC   same stores, no atomic directives.  Isolates the cost
!                       of the atomic instruction itself from the address
!                       arithmetic around it.  Races, so J and K are garbage.
!   -DTRC_NO_DIGEST   no scatter at all, just a scalar accumulation.
!                       Isolates the whole digestion, addresses included.
!
! Only ever for reading `ptxas` register counts and kernel timings.  Neither is
! a correctness path and check_binfock will fail on both, which is the point.
!
#if defined(TRC_NO_DIGEST)
         jsink = jsink + vv
#else
#if defined(TRC_NO_ATOMIC)
#define FG_ATOMIC
#else
#define FG_ATOMIC !$acc atomic update
#endif

         !
         ! ONE UPDATE PER DENSITY. These references carried two subscripts
         ! against rank-three arrays -- `jmat(mu, nu)` on
         ! `jmat(ndens, nao, nao)` -- left behind when batching added the
         ! density index. nvfortran accepts that and only WARNS
         ! (NVFORTRAN-W-0155), so it compiled, and would have written to
         ! whatever those two subscripts happened to address.
         !
         ! It stayed hidden because this is the GENERIC digestion, compiled
         ! only when the unrolled HRR is off -- which no default build has
         ! been. Forty-eight warnings appeared the moment LMAX=2 stopped
         ! unrolling the shared kernel, and not before.
         !
         ! The same trap cost a day in Phase 3 on a different array. The
         ! lesson has not changed: a rank mismatch here is a warning, not
         ! an error, and warnings scroll past.
         !
         do dn = 1, ndens
            FG_ATOMIC
            jmat(dn, mu, nu) = jmat(dn, mu, nu) + vv*dmat(dn, lam, sig)
            FG_ATOMIC
            kmat(dn, mu, lam) = kmat(dn, mu, lam) + vv*dmat(dn, nu, sig)
            if (dij) then
               FG_ATOMIC
               jmat(dn, nu, mu) = jmat(dn, nu, mu) + vv*dmat(dn, lam, sig)
               FG_ATOMIC
               kmat(dn, nu, lam) = kmat(dn, nu, lam) + vv*dmat(dn, mu, sig)
            end if
            if (dkl) then
               FG_ATOMIC
               jmat(dn, mu, nu) = jmat(dn, mu, nu) + vv*dmat(dn, sig, lam)
               FG_ATOMIC
               kmat(dn, mu, sig) = kmat(dn, mu, sig) + vv*dmat(dn, nu, lam)
            end if
            if (dij .and. dkl) then
               FG_ATOMIC
               jmat(dn, nu, mu) = jmat(dn, nu, mu) + vv*dmat(dn, sig, lam)
               FG_ATOMIC
               kmat(dn, nu, sig) = kmat(dn, nu, sig) + vv*dmat(dn, mu, lam)
            end if
            if (dpq) then
               FG_ATOMIC
               jmat(dn, lam, sig) = jmat(dn, lam, sig) + vv*dmat(dn, mu, nu)
               FG_ATOMIC
               kmat(dn, lam, mu) = kmat(dn, lam, mu) + vv*dmat(dn, sig, nu)
               if (dkl) then
                  FG_ATOMIC
                  jmat(dn, sig, lam) = jmat(dn, sig, lam) + vv*dmat(dn, mu, nu)
                  FG_ATOMIC
                  kmat(dn, sig, mu) = kmat(dn, sig, mu) + vv*dmat(dn, lam, nu)
               end if
               if (dij) then
                  FG_ATOMIC
                  jmat(dn, lam, sig) = jmat(dn, lam, sig) + vv*dmat(dn, nu, mu)
                  FG_ATOMIC
                  kmat(dn, lam, nu) = kmat(dn, lam, nu) + vv*dmat(dn, sig, mu)
               end if
               if (dij .and. dkl) then
                  FG_ATOMIC
                  jmat(dn, sig, lam) = jmat(dn, sig, lam) + vv*dmat(dn, nu, mu)
                  FG_ATOMIC
                  kmat(dn, sig, nu) = kmat(dn, sig, nu) + vv*dmat(dn, lam, mu)
               end if
            end if
         end do
#endif
      end do
      end do
      end do
      end do
#endif
#if defined(TRC_NO_DIGEST)
      ! keep jsink live so the whole evaluation is not dead-code eliminated
      if (jsink /= 0.0_dp) then
         !$acc atomic update
         jmat(1, 1) = jmat(1, 1) + jsink*1.0e-300_dp
      end if
#endif
   end subroutine one_bin_item


   subroutine plan_build(this, b, thresh, use_dens)
      class(trc_plan_t), intent(inout) :: this
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      logical, intent(in)  :: use_dens

      integer :: ia, ib, ka, kb, smax_keep, nA, nB, nseg, is
      integer(kind=8) :: nt
      integer, allocatable :: sA(:), sB(:), sOA(:), sOB(:), sNB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(kind=8), allocatable :: sOff(:)

      call this%release()

      smax_keep = int(-log10(thresh))

      ! --- pass 1: count admitted bin pairs ---
      nseg = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            if (use_dens) then
               if (dens_reject(b, ka, kb, thresh)) cycle
            end if
            nseg = nseg + 1
         end do
      end do

      allocate (sA(nseg), sB(nseg), sOA(nseg), sOB(nseg), sNB(nseg), sD(nseg))
      allocate (sLA(nseg), sLB(nseg), sLC(nseg), sLD(nseg))
      allocate (sOff(nseg + 1))

      ! --- pass 2: fill descriptors and prefix-sum the work ---
      !
      ! ONE launch over every admitted bin pair, not one launch each.  Binning
      ! at (type, size) granularity fragments the work badly -- 484 launches
      ! averaging 2200 items apiece measured 12x slower than a single large
      ! launch.  The screening granularity is worth keeping; the launch
      ! granularity is not.  A prefix sum lets a thread find its own segment.
      !
      is = 0; sOff(1) = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         nA = b%bin_cnt(ka)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            nB = b%bin_cnt(kb)
            if (nA == 0 .or. nB == 0) cycle
            if (use_dens) then
               if (dens_reject(b, ka, kb, thresh)) cycle
            end if
            is = is + 1
            sA(is) = nA; sNB(is) = nB
            sOA(is) = b%bin_off(ka); sOB(is) = b%bin_off(kb)
            sD(is) = (ka == kb)
            sLA(is) = b%bin_la(ka); sLB(is) = b%bin_lb(ka)
            sLC(is) = b%bin_la(kb); sLD(is) = b%bin_lb(kb)
            if (ka == kb) then
               nt = int(nA, 8)*int(nA + 1, 8)/2
            else
               nt = int(nA, 8)*int(nB, 8)
            end if
            !
            ! one_bin_item takes its within-segment index as a default integer,
            ! so a segment must fit in one.  At 75 waters the largest bin holds
            ! 3617 pairs and the largest segment is 6.5e6, against a 2.1e9
            ! limit -- about 65,000 pairs in a single bin before it bites,
            ! which is a system some 20x larger.  Guarded rather than widened
            ! because widening the index costs registers in the hot kernel and
            ! the register budget is already the binding constraint.
            !
            if (nt > 2147483000_8) error stop &
               "terco: bin-pair segment exceeds a default integer; widen t in one_bin_item"
            sOff(is + 1) = sOff(is) + nt
         end do
      end do


      this%nseg  = nseg
      this%nwork = sOff(nseg + 1)
      call move_alloc(sA,  this%sA)
      call move_alloc(sB,  this%sB)
      call move_alloc(sOA, this%sOA)
      call move_alloc(sOB, this%sOB)
      call move_alloc(sNB, this%sNB)
      call move_alloc(sD,  this%sD)
      call move_alloc(sOff, this%sOff)
      call move_alloc(sLA, this%sLA)
      call move_alloc(sLB, this%sLB)
      call move_alloc(sLC, this%sLC)
      call move_alloc(sLD, this%sLD)
   end subroutine plan_build

   subroutine plan_release(this)
      class(trc_plan_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%sA, this%sNB, this%sOA, this%sOB, &
         !$acc                  this%sD, this%sOff)
         this%on_device = .false.
      end if
      if (allocated(this%sA))  deallocate (this%sA)
      if (allocated(this%sB))  deallocate (this%sB)
      if (allocated(this%sOA)) deallocate (this%sOA)
      if (allocated(this%sOB)) deallocate (this%sOB)
      if (allocated(this%sNB)) deallocate (this%sNB)
      if (allocated(this%sD))  deallocate (this%sD)
      if (allocated(this%sOff)) deallocate (this%sOff)
      if (allocated(this%sLA)) deallocate (this%sLA)
      if (allocated(this%sLB)) deallocate (this%sLB)
      if (allocated(this%sLC)) deallocate (this%sLC)
      if (allocated(this%sLD)) deallocate (this%sLD)
      this%nseg = 0; this%nwork = 0; this%nlaunch = 0
   end subroutine plan_release

end module trc_binkernel
