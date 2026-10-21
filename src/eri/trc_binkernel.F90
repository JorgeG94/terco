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
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use trc_boys, only: BOYS_MMAX, boys_table, BOYS_NCHEB, BOYS_NGRID, BOYS_TMAX, BOYS_DT, BOYS_DTINV
   use trc_tables, only: LMAX
   use trc_cart, only: NCUM, cidx, cnx, cny, cnz, cll, cdir, cdn1, cdn2, cf2, &
                         ncum_of, ncart_of
   use trc_bins, only: pair_bins_t, SMAX, sort_bins_by_weight, build_pair_quartets
   use trc_pc_kernels, only: pc_dispatch, CLASS_RADIX
   use trc_rys_select, only: rys_select_init, rys_route
   use trc_rys_kernels, only: rys_dispatch
   implicit none
   private

   real(dp), protected, public :: TRC_PRIM_MARGIN = 1.0e-3_dp
   !! Fraction of the integral threshold at which a PRIMITIVE quartet is
   !! abandoned before the Boys function. The prefactor bounds the (ss|ss)
   !! primitive integral exactly and the higher classes only up to the
   !! polynomial factors of the recurrences, so the margin is empirical: a
   !! thousandfold was the original guess and nothing had measured what
   !! tightening it costs. This makes that measurable.

   public :: trc_set_prim_margin
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
      integer(int64) :: nwork = 0
      integer, allocatable :: sA(:), sB(:), sOA(:), sOB(:), sNB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(int64), allocatable :: sOff(:)
      logical :: on_device = .false.
   contains
      procedure :: build   => plan_build
      procedure :: release => plan_release
   end type trc_plan_t


contains

   integer(int64) function launch_cap() result(res)
      !! The one-time pair ordering, off unless asked for. See the comment at
      !! its call site for what it costs and buys.
      !! TRC_LAUNCH_CAP caps the quartets in one launch; 0 (default) is one
      !! launch per angular-momentum class.
      character(len=32) :: e
      integer :: ios
      e = " "
      call get_environment_variable("TRC_LAUNCH_CAP", e)
      res = 0_int64
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) res
         if (ios /= 0) res = 0_int64
      end if
   end function launch_cap

   logical function seg_stats_on() result(res)
      !! TRC_SEG_STATS=1 prints the segment shape of each Fock build.
      character(len=8) :: e
      e = " "
      call get_environment_variable("TRC_SEG_STATS", e)
      res = len_trim(e) > 0
   end function seg_stats_on

   logical function bin_sort_on() result(res)
      character(len=8) :: v
      call get_environment_variable("TRC_BIN_SORT", v)
      res = (trim(v) == "1")
   end function bin_sort_on

   subroutine trc_set_prim_margin(f)
      real(dp), intent(in) :: f
      if (f > 0.0_dp) TRC_PRIM_MARGIN = f
   end subroutine trc_set_prim_margin

#include "inc/trc_boys_eval.inc"

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
   !
   ! (l, ic) -> (px,py,pz), the component AT AN AO OFFSET (`ic` becomes
   ! `mu = ao_off(si) + ic` etc below) -- OUTPUT layout, follows the build's
   ! Cartesian order.  Generated so this deliberate second copy of trc_hgp's
   ! `comp` cannot drift from it; see scripts/gen_cartcomp.py.
   !
   pure subroutine comp(l, ic, px, py, pz)
      !$acc routine seq
      integer, intent(in)  :: l, ic
      integer, intent(out) :: px, py, pz
#include "trc_cart_comp.inc"
   end subroutine comp

   pure real(dp) function binom(n, k) result(res)
      !$acc routine seq
      integer, intent(in) :: n, k
      integer :: i
      real(dp) :: r
      r = 1.0_dp
      do i = 1, k
         r = r*real(n - k + i, dp)/real(i, dp)
      end do
      res = r
   end function binom

   pure real(dp) function powi(x, n) result(res)
      !$acc routine seq
      real(dp), intent(in) :: x
      integer,  intent(in) :: n
      integer :: i
      res = 1.0_dp
      do i = 1, n
         res = res*x
      end do
   end function powi

   !
   ! Whole Fock build: walk admitted bin pairs on the host, launch one kernel
   ! each.  The host loop is over BINS -- tens of thousands of iterations at
   ! most -- not over quartets.
   !
   subroutine fock_bins(b, nbas, npp, nao, sh_l, ao_off, thresh, use_dens, &
                        jfac, kfac, omega, nosym, dsh, &
                        pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                        qp_bin_pbase, qp_bin_k, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                        sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                        ndens, dmat, jmat, kmat, rank, nranks, nlaunch, nwork, nkept)
      type(pair_bins_t), intent(inout) :: b
      integer,  intent(in)    :: nbas, npp, nao
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      real(dp), intent(in)    :: thresh
      real(dp), intent(in)    :: jfac, kfac
      !! Coulomb and exchange scalings; 1 and 1 is Hartree-Fock.
      real(dp), intent(in)    :: omega
      !! Range-separation parameter, bohr^-1; zero is the Coulomb operator.
      !! See scripts/gen_perclass.py's primitive loop for the substitution,
      !! shared verbatim with the enumerated kernels below.
      logical,  intent(in)    :: nosym
      !! Contract a density that is NOT symmetric, through the enumerated
      !! kernel. Costs what folding bought; needed by the response Hessian.
      real(dp), intent(in)    :: dsh(nbas, nbas)
      !! max |D| over each shell-pair block; the per-quartet density screen.
      !! Pass it filled with huge() to disable density screening entirely.
      logical,  intent(in)    :: use_dens
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      integer,  allocatable, intent(inout) :: qp_bin_pbase(:), qp_bin_k(:)
      !! Bin-interleaved second copy of the above (trc_eri_t's own flat
      !! fields -- see build_pair_quartets for why they are not part of `b`).
      !! `intent(inout)`, not `in`: a resort rebuilds them in place.
      integer,  intent(inout) :: nqp
      real(dp), allocatable, intent(inout) :: qp_z(:), qp_c(:), qp_px(:), qp_py(:), qp_pz(:)
      real(dp), allocatable, intent(inout) :: sq_ax(:), sq_ay(:), sq_az(:)
      real(dp), allocatable, intent(inout) :: sq_abx(:), sq_aby(:), sq_abz(:)
      integer,  allocatable, intent(inout) :: sq_n(:)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)
      integer,  intent(in)    :: rank, nranks
      !! This rank's share of the sorted work list: every nranks-th item from
      !! rank. The caller reduces the result across ranks.
      integer,  intent(out)   :: nlaunch
      integer(int64), intent(out) :: nwork
      integer(int64), intent(out), optional :: nkept
      !! Diagnostic: how many of `nwork` actually survive the in-kernel
      !! Schwarz and density tests.  Costs a separate cheap pass, so it is
      !! optional and off in the timed path.
      integer :: ia, ib, ka, kb, smax_keep, nA, nB, nseg, is
      integer(int64) :: nt
      integer, allocatable :: sA(:), sB(:), sOA(:), sOB(:), sNB(:), sPA(:), sPB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(int64), allocatable :: sOff(:)

      ! Order the pairs inside each bin by how much they matter, once, so
      ! that neighbouring threads of a warp agree about screening. The
      ! device copy has to be refreshed with the new order.
      !
      ! OFF by default, on measurement. Ordering the pairs inside a bin by
      ! how much they matter should make a warp's threads agree about
      ! screening, and on the 28-atom cage it does: 3.7% off the cc-pVDZ
      ! Fock time. On the 123-atom slice the same change COSTS 2.8%, with
      ! either key -- the product with the density block, or the Schwarz
      ! bound alone. The natural order follows the shell index, which
      ! carries spatial locality into the density and Fock accesses, and on
      ! a large system losing that costs more than the divergence saves.
      ! Set TRC_BIN_SORT=1 to turn it on; TRC_BIN_SORT_Q=1 keys on the
      ! bound alone.
      !
      if (.not. b%sorted .and. bin_sort_on()) then
         call sort_bins_by_weight(b, nbas, dsh)
         !$acc update device(b%sp_i, b%sp_j, b%sp_q)
         ! The slot order inside a bin just changed, and the bin-interleaved
         ! tables are keyed on it -- rebuild them from the same pp_* data
         ! build_pairs_hgp filled, and refresh the device copy.
         call build_pair_quartets(b, nbas, npp, pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                  qp_bin_pbase, qp_bin_k, nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                                  sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n)
         !$acc update device(qp_z, qp_c, qp_px, qp_py, qp_pz, &
         !$acc                sq_ax, sq_ay, sq_az, &
         !$acc                sq_abx, sq_aby, sq_abz, sq_n)
      end if

      ! In buckets, not decades: the bins carry the resolution they were
      ! built at, and the admission test has to speak the same units.
      smax_keep = int(-log10(thresh))*b%sres

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
      allocate (sPA(nseg), sPB(nseg))
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
            sPA(is) = qp_bin_pbase(ka); sPB(is) = qp_bin_pbase(kb)
            sD(is) = (ka == kb)
            sLA(is) = b%bin_la(ka); sLB(is) = b%bin_lb(ka)
            sLC(is) = b%bin_la(kb); sLD(is) = b%bin_lb(kb)
            if (ka == kb) then
               nt = int(nA, int64)*int(nA + 1, int64)/2
            else
               nt = int(nA, int64)*int(nB, int64)
            end if
            !
            ! The kernels take their within-segment index as a default integer,
            ! so a segment must fit in one. Widening it would cost registers
            ! in the hot kernel, which are the binding constraint; instead
            ! build_binned_pairs caps a bin at BIN_CAP = 32768 pairs
            ! (split_large_bins, trc_bins.F90), so no segment exceeds 2**30.
            ! This guard is the backstop if that cap is ever raised.
            !
            if (nt > 2147483000_int64) error stop &
               "terco: bin-pair segment exceeds a default integer; widen the kernels' t"
            sOff(is + 1) = sOff(is) + nt
         end do
      end do

      nwork = sOff(nseg + 1)
      nlaunch = 1
      ! Segment shape, for comparing this against a batching scheme that caps
      ! its batches (gmshpc: 2560 shell pairs per bra batch). A segment here
      ! is one admitted bin pair and is homogeneous by construction -- same
      ! four angular momenta, same contraction degree, same magnitude bucket
      ! -- so what matters is how big they are and how much of a warp
      ! straddles a boundary.
      if (seg_stats_on()) then
         block
            integer :: q, nsmall
            integer(int64) :: mn, mx, sz
            mn = huge(0_int64); mx = 0_int64; nsmall = 0
            do q = 1, nseg
               sz = sOff(q + 1) - sOff(q)
               mn = min(mn, sz); mx = max(mx, sz)
               if (sz < 32) nsmall = nsmall + 1
            end do
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,i0,a,i0,a,i0,a,f8.1,a,i0,a,f5.2,a)") &
                  "  [seg] nseg ", nseg, "  work ", nwork, "  min ", mn, "  max ", mx, &
                  "  mean ", real(nwork, dp)/max(nseg, 1), "  <32 ", nsmall, &
                  "  boundary warps ", 100.0_dp*nseg*32.0_dp/max(real(nwork, dp), 1.0_dp), "% of threads"
               call global_logger%info(trim(logline))
            end block
         end block
      end if
      if (nwork == 0) return

      if (nosym) then
         ! Every per-class kernel is folded, so a general density goes through
         ! the enumerated non-symmetric kernel, one_bin_item_nosym, instead.
         !$acc enter data copyin(sA, sNB, sOA, sOB, sD, sOff)
         call fock_all_nosym(nseg, nwork, sA, sNB, sOA, sOB, sD, sOff, &
                             b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, jfac, kfac, omega, dsh, &
                             nbas, npp, nao, sh_l, ao_off, &
                             pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                             ndens, dmat, jmat, kmat, rank, nranks)
         nlaunch = 1
         !$acc exit data delete(sA, sNB, sOA, sOB, sD, sOff)
         deallocate (sA, sB, sOA, sOB, sNB, sPA, sPB, sD, sOff, sLA, sLB, sLC, sLD)
         return
      end if

      !
      ! Per-class launches.  The descriptors are sorted by class so each class
      ! occupies a contiguous run of segments, then one kernel per class covers
      ! its own run.  Every thread in a launch then executes code compiled for
      ! exactly its (la,lb,lc,ld) -- no select case, and ptxas budgets registers
      ! per class instead of for the worst class across all of them.
      !
      ! TRC_RYS_CLASSES, read once per Fock build -- `fock_bins` IS one Fock
      ! build, called once per SCF iteration (or once, standalone, by a
      ! correctness check) -- not once per launch below. See trc_rys_select.
      call rys_select_init()
      block
         integer, allocatable :: ord(:), ckey(:)
         integer :: a2, b2, t2, c0, c1, nl, key2
         integer(int64) :: lcap
         logical :: l2
         integer(int64) :: o2
         ! The trivial primitive-shell view used to be built here -- eight
         ! arrays, allocated, filled, and pushed to the device on every Fock
         ! build -- purely to fill dummy arguments the scalar kernel never
         ! read. They existed for the blocked kernel, which is gone with the
         ! general contraction it was written for; `trc_decontract` splits
         ! the basis instead. Two of them were once added to `enter data`
         ! and not to `exit data`, which is a stale mapping and an illegal
         ! address rather than a wrong number.
         allocate (ord(nseg), ckey(nseg))
         do a2 = 1, nseg
            ckey(a2) = ((sLA(a2)*CLASS_RADIX + sLB(a2))*CLASS_RADIX &
                        + sLC(a2))*CLASS_RADIX + sLD(a2)
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
         call permute_segments(nseg, ord, sA, sNB, sOA, sOB, sPA, sPB, sD, sOff, &
                               sLA, sLB, sLC, sLD)

         !$acc enter data copyin(sA, sNB, sOA, sOB, sPA, sPB, sD, sOff)
         if (present(nkept)) call count_kept(nseg, nwork, sNB, sOA, sOB, sD, &
                                             sOff, b%npair, b%sp_i, b%sp_j, &
                                             b%sp_q, thresh, dsh, nbas, nkept)
         nl = 0
         c0 = 1

         lcap = launch_cap()
         do while (c0 <= nseg)
            c1 = c0
            do while (c1 < nseg)
               if (((sLA(c1 + 1)*CLASS_RADIX + sLB(c1 + 1))*CLASS_RADIX &
                    + sLC(c1 + 1))*CLASS_RADIX + sLD(c1 + 1) /= &
                   ((sLA(c0)*CLASS_RADIX + sLB(c0))*CLASS_RADIX &
                    + sLC(c0))*CLASS_RADIX + sLD(c0)) exit
               if (lcap > 0 .and. sOff(c1 + 2) - sOff(c0) > lcap) exit
               c1 = c1 + 1
            end do
            nl = nl + 1
            key2 = ((sLA(c0)*CLASS_RADIX + sLB(c0))*CLASS_RADIX &
                    + sLC(c0))*CLASS_RADIX + sLD(c0)
            ! TRC_RYS_CLASSES routes this launch's whole class run to the
            ! Rys kernel instead of HGP's, for classes that have one --
            ! `rys_route` already checked that. Out of scope otherwise: Rys
            ! never changes HGP's own dispatch, decode or launch geometry.
            if (rys_route(key2)) then
               call rys_dispatch(key2, c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                 b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*TRC_PRIM_MARGIN, jfac, kfac, &
                                 omega, dsh, nbas, nao, sh_l, ao_off, &
                                 nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                                 sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                 ndens, dmat, jmat, rank, nranks)
            else
               call pc_dispatch(key2, c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*TRC_PRIM_MARGIN, jfac, kfac, &
                                omega, dsh, nbas, nao, sh_l, ao_off, &
                                nqp, qp_z, qp_c, qp_px, qp_py, qp_pz, &
                                sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                ndens, dmat, jmat, rank, nranks)
            end if
            c0 = c1 + 1
         end do
         nlaunch = nl
         !$acc exit data delete(sA, sNB, sOA, sOB, sPA, sPB, sD, sOff)
         deallocate (ord, ckey)
      end block
      deallocate (sA, sB, sOA, sOB, sNB, sPA, sPB, sD, sOff, sLA, sLB, sLC, sLD)
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
   pure logical function dens_reject(b, ka, kb, thresh) result(res)
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
      integer(int64), intent(in) :: nwork
      integer,  intent(in) :: sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in) :: sD(nseg)
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in) :: sp_i(npair), sp_j(npair)
      real(dp), intent(in) :: thresh
      real(dp), intent(in) :: sp_q(npair), dsh(nbas, nbas)
      integer(int64), intent(out) :: nkept

      integer(int64) :: gt, nk
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

   subroutine permute_segments(n, ord, sA, sNB, sOA, sOB, sPA, sPB, sD, sOff, sLA, sLB, sLC, sLD)
      !! Apply a permutation to the segment descriptor arrays, in place.
      integer, intent(in) :: n
      integer, intent(in) :: ord(n)
      integer, intent(inout) :: sA(n), sNB(n), sOA(n), sOB(n), sPA(n), sPB(n)
      integer, intent(inout) :: sLA(n), sLB(n), sLC(n), sLD(n)
      logical, intent(inout) :: sD(n)
      integer(int64), intent(inout) :: sOff(n + 1)
      integer :: k
      integer, allocatable :: t(:)
      logical, allocatable :: tl(:)
      integer(int64), allocatable :: to(:)
      allocate (t(n), tl(n), to(n + 1))
      t = sA(ord);  sA = t
      t = sNB(ord); sNB = t
      t = sOA(ord); sOA = t
      t = sOB(ord); sOB = t
      t = sPA(ord); sPA = t
      t = sPB(ord); sPB = t
      t = sLA(ord); sLA = t
      t = sLB(ord); sLB = t
      t = sLC(ord); sLC = t
      t = sLD(ord); sLD = t
      tl = sD(ord); sD = tl
      ! rebuild the prefix sum in the new order
      to(1) = 0
      do k = 1, n
         if (sD(k)) then
            to(k + 1) = to(k) + int(sA(k), int64)*int(sA(k) + 1, int64)/2
         else
            to(k + 1) = to(k) + int(sA(k), int64)*int(sNB(k), int64)
         end if
      end do
      sOff = to
      deallocate (t, tl, to)
   end subroutine permute_segments

   subroutine fock_all_nosym(nseg, nwork, sA, sNB, sOA, sOB, sD, sOff, &
                             npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, omega, dsh, &
                             nbas, npp, nao, sh_l, ao_off, &
                             pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                             ndens, dmat, jmat, kmat, rank, nranks)
      !! Whole-worklist driver for the general-density kernel.
      integer,  intent(in)    :: nseg, npair, nbas, npp, nao, ndens, rank, nranks
      integer(int64), intent(in) :: nwork
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in) :: thresh, jfac, kfac, omega
      real(dp), intent(in) :: sp_q(npair)
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas)
      integer,  intent(in)    :: pp_off(nbas*nbas), pp_n(nbas*nbas)
      real(dp), intent(in)    :: pp_p(npp), pp_r(npp, 3), pp_ra(npp, 3), pp_rb(npp, 3), pp_c(npp)
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

      integer(int64) :: gt, nr, i
      integer :: lo, hi, mid, seg, t

      ! `local()` is NOT optional here, and the reason is worth stating.
      !
      ! `do concurrent` privatises a scalar assigned before it is read, so
      ! this loop is correct without it. But the OpenMP port converts the
      ! construct to `!$omp parallel do`, where the default is SHARED -- and
      ! then `lo`, `hi`, `mid`, `seg` and `t` are one set of variables torn
      ! between every thread. That is a race, and it showed up as this
      ! routine disagreeing with the folded digestion on a symmetric density
      ! after conversion, which is exactly the failure `check_nosym` exists
      ! to catch.
      !
      ! `tools/dc_locality_lint.py` now refuses a `do concurrent` that
      ! assigns a scalar it has not named.
      nr = 0
      if (nwork >= rank + 1) nr = (nwork - rank - 1)/nranks + 1
      do concurrent(i=1:nr) local(gt, lo, hi, mid, seg, t)
         gt = rank + 1 + (i - 1)*nranks
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
                                 npair, sp_i, sp_j, sp_q, thresh, jfac, kfac, omega, dsh, &
                                 nbas, npp, nao, sh_l, ao_off, &
                                 pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                 ndens, dmat, jmat, kmat)
      end do
   end subroutine fock_all_nosym

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
                                sp_q, thresh, jfac, kfac, omega, dsh, nbas, npp, nao, sh_l, ao_off, &
                                pp_off, pp_n, pp_p, pp_r, pp_ra, pp_rb, pp_c, &
                                ndens, dmat, jmat, kmat)
      !$acc routine seq
      integer,  intent(in)    :: t, nA, nB, offA, offB, npair, nbas, npp, nao
      logical,  intent(in)    :: diag
      integer,  intent(in)    :: sp_i(npair), sp_j(npair)
      real(dp), intent(in) :: thresh
      real(dp), intent(in) :: sp_q(npair)
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: omega
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
      real(dp) :: theta, rho_theta
      real(dp) :: pa(3), qc(3), wp(3), wq(3), pq(3)
      real(dp) :: pax, pay, paz, qcx, qcy, qcz
      real(dp) :: wpx, wpy, wpz, wqx, wqy, wqz
      real(dp) :: oo2z, oo2e, oo2ze, rz, re
      real(dp) :: f(0:BOYS_MMAX)
      real(dp) :: v(NCUM*NCUM, 0:1), g(NCUM*NCUM)
      integer  :: na_, nb_, nc_, nd_, ib_, id_, mu, nu, lam, sig, idx
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: px, py, pz, qx, qy, qz, ja, jc
      real(dp) :: cb, cdc, vv
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

      do x = 1, nca*ncc
         g(x) = 0.0_dp
      end do

      do kp = offab + 1, offab + nab
         zeta = pp_p(kp)
         do kq = offcd + 1, offcd + ncd
            eta = pp_p(kq)
            zpe = zeta + eta
            rho = zeta*eta/zpe
            ! RANGE SEPARATION -- same substitution as the generated kernels
            ! (scripts/gen_perclass.py's primitive loop): erf(omega*r12)/r12
            ! is the Coulomb integral with rho replaced everywhere by
            ! rho_theta = theta*rho, theta = omega^2/(rho + omega^2), and the
            ! prefactor scaled by sqrt(theta). omega == 0 sets theta = 1.0_dp
            ! exactly, so every "theta*"/"rho_theta" below reproduces the
            ! Coulomb value bit-for-bit -- there is no screen here to protect
            ! (this whole bin pair already survived one at admission), so
            ! theta is computed unconditionally rather than branched.
            if (omega > 0.0_dp) then
               theta = omega*omega/(rho + omega*omega)
            else
               theta = 1.0_dp
            end if
            rho_theta = theta*rho
            do d = 1, 3
               pq(d) = pp_r(kp, d) - pp_r(kq, d)
               wc = (zeta*pp_r(kp, d) + eta*pp_r(kq, d))/zpe
               pa(d) = pp_r(kp, d) - pp_ra(kp, d)
               qc(d) = pp_r(kq, d) - pp_ra(kq, d)
               wp(d) = theta*(wc - pp_r(kp, d))
               wq(d) = theta*(wc - pp_r(kq, d))
            end do
            tval = rho_theta*(pq(1)*pq(1) + pq(2)*pq(2) + pq(3)*pq(3))
            call boys_eval(lt, tval, f)

            oo2z = 0.5_dp/zeta; oo2e = 0.5_dp/eta; oo2ze = theta*(0.5_dp/zpe)
            rz = rho_theta/zeta; re = rho_theta/eta
            pref = sqrt(theta)*TWO_PI_2_5/(zeta*eta*sqrt(zpe))*pp_c(kp)*pp_c(kq)

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

         !$acc atomic update
         jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, lam, sig)
         !$acc atomic update
         kmat(1, mu, lam) = kmat(1, mu, lam) + vv*dmat(1, nu, sig)
         if (dij) then
            !$acc atomic update
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, lam, sig)
            !$acc atomic update
            kmat(1, nu, lam) = kmat(1, nu, lam) + vv*dmat(1, mu, sig)
         end if
         if (dkl) then
            !$acc atomic update
            jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, sig, lam)
            !$acc atomic update
            kmat(1, mu, sig) = kmat(1, mu, sig) + vv*dmat(1, nu, lam)
         end if
         if (dij .and. dkl) then
            !$acc atomic update
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, sig, lam)
            !$acc atomic update
            kmat(1, nu, sig) = kmat(1, nu, sig) + vv*dmat(1, mu, lam)
         end if
         if (dpq) then
            !$acc atomic update
            jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, mu, nu)
            !$acc atomic update
            kmat(1, lam, mu) = kmat(1, lam, mu) + vv*dmat(1, sig, nu)
            if (dkl) then
               !$acc atomic update
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, mu, nu)
               !$acc atomic update
               kmat(1, sig, mu) = kmat(1, sig, mu) + vv*dmat(1, lam, nu)
            end if
            if (dij) then
               !$acc atomic update
               jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, nu, mu)
               !$acc atomic update
               kmat(1, lam, nu) = kmat(1, lam, nu) + vv*dmat(1, sig, mu)
            end if
            if (dij .and. dkl) then
               !$acc atomic update
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, nu, mu)
               !$acc atomic update
               kmat(1, sig, nu) = kmat(1, sig, nu) + vv*dmat(1, lam, mu)
            end if
         end if
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


         !$acc atomic update
         jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, lam, sig)
         !$acc atomic update
         kmat(1, mu, lam) = kmat(1, mu, lam) + vv*dmat(1, nu, sig)
         if (dij) then
            !$acc atomic update
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, lam, sig)
            !$acc atomic update
            kmat(1, nu, lam) = kmat(1, nu, lam) + vv*dmat(1, mu, sig)
         end if
         if (dkl) then
            !$acc atomic update
            jmat(1, mu, nu) = jmat(1, mu, nu) + vv*dmat(1, sig, lam)
            !$acc atomic update
            kmat(1, mu, sig) = kmat(1, mu, sig) + vv*dmat(1, nu, lam)
         end if
         if (dij .and. dkl) then
            !$acc atomic update
            jmat(1, nu, mu) = jmat(1, nu, mu) + vv*dmat(1, sig, lam)
            !$acc atomic update
            kmat(1, nu, sig) = kmat(1, nu, sig) + vv*dmat(1, mu, lam)
         end if
         if (dpq) then
            !$acc atomic update
            jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, mu, nu)
            !$acc atomic update
            kmat(1, lam, mu) = kmat(1, lam, mu) + vv*dmat(1, sig, nu)
            if (dkl) then
               !$acc atomic update
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, mu, nu)
               !$acc atomic update
               kmat(1, sig, mu) = kmat(1, sig, mu) + vv*dmat(1, lam, nu)
            end if
            if (dij) then
               !$acc atomic update
               jmat(1, lam, sig) = jmat(1, lam, sig) + vv*dmat(1, nu, mu)
               !$acc atomic update
               kmat(1, lam, nu) = kmat(1, lam, nu) + vv*dmat(1, sig, mu)
            end if
            if (dij .and. dkl) then
               !$acc atomic update
               jmat(1, sig, lam) = jmat(1, sig, lam) + vv*dmat(1, nu, mu)
               !$acc atomic update
               kmat(1, sig, nu) = kmat(1, sig, nu) + vv*dmat(1, lam, mu)
            end if
         end if
      end do
      end do
      end do
      end do
#endif
   end subroutine one_bin_item_nosym



   subroutine plan_build(this, b, thresh, use_dens)
      class(trc_plan_t), intent(inout) :: this
      type(pair_bins_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      logical, intent(in)  :: use_dens

      integer :: ia, ib, ka, kb, smax_keep, nA, nB, nseg, is
      integer(int64) :: nt
      integer, allocatable :: sA(:), sB(:), sOA(:), sOB(:), sNB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(int64), allocatable :: sOff(:)

      call this%release()

      ! In buckets, not decades: the bins carry the resolution they were
      ! built at, and the admission test has to speak the same units.
      smax_keep = int(-log10(thresh))*b%sres

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
               nt = int(nA, int64)*int(nA + 1, int64)/2
            else
               nt = int(nA, int64)*int(nB, int64)
            end if
            !
            ! The kernels take their within-segment index as a default integer,
            ! so a segment must fit in one. Widening it would cost registers
            ! in the hot kernel, which are the binding constraint; instead
            ! build_binned_pairs caps a bin at BIN_CAP = 32768 pairs
            ! (split_large_bins, trc_bins.F90), so no segment exceeds 2**30.
            ! This guard is the backstop if that cap is ever raised.
            !
            if (nt > 2147483000_int64) error stop &
               "terco: bin-pair segment exceeds a default integer; widen the kernels' t"
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
