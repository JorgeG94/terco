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
                         ncum_of, cart_ncart_of
   use trc_bins, only: pair_bins_t, trc_bin_plan_t, bins_sort_by_weight, bins_build_pair_quartets
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t, trc_bin_launch_t
   use trc_pc_kernels, only: pc_dispatch, CLASS_RADIX
   use trc_rys_select, only: rys_select_init, rys_select_route
   use trc_rys_kernels, only: rys_dispatch
   use trc_mdg, only: mdg_routes
   implicit none
   private

   real(dp), protected, public :: TRC_PRIM_MARGIN = 1.0e-3_dp

   integer, parameter :: LNOSYM = min(LMAX, 3)
   !! The enumerated non-symmetric kernel stops at f: a basis with g refuses
   !! eri%fock_nosym before reaching it (trc_eri). Sized by LMAX at LMAX = 4
   !! its per-thread scratch was 654 KB, past the 512 KB GPU stack limit.
   integer, parameter :: NCUM_NOSYM = (2*LNOSYM + 1)*(2*LNOSYM + 2)*(2*LNOSYM + 3)/6
   !! trc_cart's NCUM at the capped ceiling.
   !! Fraction of the integral threshold at which a PRIMITIVE quartet is
   !! abandoned before the Boys function. The prefactor bounds the (ss|ss)
   !! primitive integral exactly and the higher classes only up to the
   !! polynomial factors of the recurrences, so the margin is empirical: a
   !! thousandfold was the original guess and nothing had measured what
   !! tightening it costs. This makes that measurable.

   public :: trc_set_prim_margin
   public :: binkernel_fock_bins

   real(dp), parameter :: TWO_PI_2_5 = 34.986836655249725_dp


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
   subroutine binkernel_fock_bins(b, pairs, quartets, launch, dsh, dmat, jmat, kmat)
      type(pair_bins_t), intent(inout) :: b
      !! `intent(inout)`: the optional one-time resort below reorders it.
      type(trc_shell_pairs_t), intent(in) :: pairs
      !! The shells and their primitive pairs, device-resident.
      type(trc_pair_quartets_t), intent(inout) :: quartets
      !! The bin-interleaved pairs, device-resident. `intent(inout)`, not
      !! `in`: a resort rebuilds them in place.
      type(trc_bin_launch_t), intent(inout) :: launch
      !! What to build (thresh, jfac/kfac, omega, nosym, ranks, operand
      !! shape); nlaunch, nwork and, if launch%count_kept, nkept come back.
      !! jfac/kfac 1 and 1 is Hartree-Fock. omega is the range-separation
      !! parameter, bohr^-1, zero the Coulomb operator: see
      !! scripts/gen_perclass.py's primitive loop for the substitution,
      !! shared verbatim with the enumerated kernel below. nosym contracts a
      !! density that is NOT symmetric, through the enumerated kernel; it
      !! costs what folding bought and is needed by the response Hessian.
      !! The rank takes every nranks-th item of the sorted work list; the
      !! caller reduces the result across ranks.
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      !! max |D| over each shell-pair block; the per-quartet density screen.
      !! Pass it filled with huge() to disable density screening entirely.
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(inout) :: jmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(inout) :: kmat(launch%ndens, launch%nao, launch%nao)

      type(trc_bin_plan_t) :: plan
      real(dp) :: thresh
      integer :: nseg

      thresh = launch%thresh
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
         call bins_sort_by_weight(b, pairs%nbas, dsh)
         !$acc update device(b%sp_i, b%sp_j, b%sp_q)
         ! The slot order inside a bin just changed, and the bin-interleaved
         ! tables are keyed on it -- rebuild them from the same pairs
         ! hgp_build_pairs filled, and refresh the device copy.
         call bins_build_pair_quartets(b, pairs, quartets)
         !$acc update device(quartets%z, quartets%c, quartets%a, &
         !$acc                quartets%px, quartets%py, quartets%pz, &
         !$acc                quartets%ax, quartets%ay, quartets%az, &
         !$acc                quartets%abx, quartets%aby, quartets%abz, quartets%n)
      end if

      call plan%build(b, thresh, quartets%bin_pbase, launch%use_density)
      nseg = plan%nseg
      launch%nwork = plan%nwork
      launch%nlaunch = 1
      if (seg_stats_on()) call log_segment_shape(plan)
      if (plan%nwork == 0) return

      if (launch%nosym) then
         ! Every per-class kernel is folded, so a general density goes through
         ! the enumerated non-symmetric kernel, one_bin_item_nosym, instead.
         !$acc enter data copyin(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sD, plan%sOff)
         call fock_all_nosym(nseg, plan%nwork, plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sD, plan%sOff, &
                             b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, launch%jfac, launch%kfac, &
                             launch%omega, dsh, pairs%nbas, pairs%npp, launch%nao, pairs%sh_l, pairs%ao_off, &
                             pairs%off, pairs%n, pairs%p, pairs%r, pairs%ra, pairs%rb, pairs%c, &
                             launch%ndens, dmat, jmat, kmat, launch%rank, launch%nranks)
         !$acc exit data delete(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sD, plan%sOff)
         call plan%release()
         return
      end if

      ! TRC_RYS_CLASSES, read once per Fock build -- `binkernel_fock_bins` IS one Fock
      ! build, called once per SCF iteration (or once, standalone, by a
      ! correctness check) -- not once per launch below. See trc_rys_select.
      call rys_select_init()
      ! The trivial primitive-shell view used to be built here -- eight
      ! arrays, allocated, filled, and pushed to the device on every Fock
      ! build -- purely to fill dummy arguments the scalar kernel never
      ! read. They existed for the blocked kernel, which is gone with the
      ! general contraction it was written for; `trc_decontract` splits
      ! the basis instead. Two of them were once added to `enter data`
      ! and not to `exit data`, which is a stale mapping and an illegal
      ! address rather than a wrong number.
      call plan%sort_by_class(CLASS_RADIX)

      !$acc enter data copyin(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, plan%sOff)
      if (launch%count_kept) call count_kept(nseg, plan%nwork, plan%sNB, plan%sOA, plan%sOB, plan%sD, &
                                             plan%sOff, b%npair, b%sp_i, b%sp_j, &
                                             b%sp_q, thresh, dsh, pairs%nbas, launch%nkept)
      call launch_classes(b, pairs, quartets, launch, plan, dsh, dmat, jmat)
      !$acc exit data delete(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, plan%sOff)
      call plan%release()
   end subroutine binkernel_fock_bins

   !
   ! One launch per contiguous class run of a class-sorted plan (or per
   ! TRC_LAUNCH_CAP quartets of one), each through the class's own HGP or
   ! Rys kernel. The plan's arrays are already on the device.
   !
   subroutine launch_classes(b, pairs, quartets, launch, plan, dsh, dmat, jmat)
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      type(trc_pair_quartets_t), intent(in) :: quartets
      type(trc_bin_launch_t), intent(inout) :: launch
      !! launch%nlaunch comes back.
      type(trc_bin_plan_t), intent(in) :: plan
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(inout) :: jmat(launch%ndens, launch%nao, launch%nao)
      !! The folded accumulator: the per-class kernels add J and K into it.

      integer :: c0, c1, nl, key, nseg
      integer(int64) :: lcap
      real(dp) :: thresh

      nseg = plan%nseg
      thresh = launch%thresh
      nl = 0
      c0 = 1
      lcap = launch_cap()
      do while (c0 <= nseg)
         key = plan%class_key(c0, CLASS_RADIX)
         c1 = c0
         do while (c1 < nseg)
            if (plan%class_key(c1 + 1, CLASS_RADIX) /= key) exit
            if (lcap > 0 .and. plan%sOff(c1 + 2) - plan%sOff(c0) > lcap) exit
            c1 = c1 + 1
         end do
         nl = nl + 1
         ! TRC_RYS_CLASSES routes this launch's whole class run to the
         ! Rys kernel instead of HGP's, for classes that have one --
         ! `rys_select_route` already checked that. Out of scope otherwise: Rys
         ! never changes HGP's own dispatch, decode or launch geometry.
         if (mdg_routes(key)) then
            ! A g class above TERCO_RYS_G_LTOT: done after the bins by the
            ! MD exchange pass (trc_mdg), which orders the work its own way.
            continue
         else if (rys_select_route(key)) then
            call rys_dispatch(key, c0, c1, nseg, plan%sOff, plan%sA, plan%sNB, plan%sOA, plan%sOB, &
                              plan%sPA, plan%sPB, plan%sD, &
                              b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*TRC_PRIM_MARGIN, &
                              launch%jfac, launch%kfac, launch%omega, dsh, pairs%nbas, launch%nao, &
                              pairs%sh_l, pairs%ao_off, &
                              quartets%nqp, quartets%z, quartets%c, quartets%px, quartets%py, quartets%pz, &
                              quartets%ax, quartets%ay, quartets%az, &
                              quartets%abx, quartets%aby, quartets%abz, quartets%n, &
                              launch%ndens, dmat, jmat, launch%rank, launch%nranks)
         else
            call pc_dispatch(key, c0, c1, nseg, plan%sOff, plan%sA, plan%sNB, plan%sOA, plan%sOB, &
                             plan%sPA, plan%sPB, plan%sD, &
                             b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*TRC_PRIM_MARGIN, &
                             launch%jfac, launch%kfac, launch%omega, dsh, pairs%nbas, launch%nao, &
                             pairs%sh_l, pairs%ao_off, &
                             quartets%nqp, quartets%z, quartets%c, quartets%px, quartets%py, quartets%pz, &
                             quartets%ax, quartets%ay, quartets%az, &
                             quartets%abx, quartets%aby, quartets%abz, quartets%n, &
                             launch%ndens, dmat, jmat, launch%rank, launch%nranks)
         end if
         c0 = c1 + 1
      end do
      launch%nlaunch = nl
   end subroutine launch_classes

   !
   ! Segment shape, for comparing this against a batching scheme that caps
   ! its batches (gmshpc: 2560 shell pairs per bra batch). A segment here is
   ! one admitted bin pair and is homogeneous by construction -- same four
   ! angular momenta, same contraction degree, same magnitude bucket -- so
   ! what matters is how big they are and how much of a warp straddles a
   ! boundary.
   !
   subroutine log_segment_shape(plan)
      type(trc_bin_plan_t), intent(in) :: plan
      integer :: q, nsmall
      integer(int64) :: mn, mx, sz
      character(len=1024) :: logline
      mn = huge(0_int64)
      mx = 0_int64
      nsmall = 0
      do q = 1, plan%nseg
         sz = plan%sOff(q + 1) - plan%sOff(q)
         mn = min(mn, sz)
         mx = max(mx, sz)
         if (sz < 32) nsmall = nsmall + 1
      end do
      write (logline, "(a,i0,a,i0,a,i0,a,i0,a,f8.1,a,i0,a,f5.2,a)") &
         "  [seg] nseg ", plan%nseg, "  work ", plan%nwork, "  min ", mn, "  max ", mx, &
         "  mean ", real(plan%nwork, dp)/max(plan%nseg, 1), "  <32 ", nsmall, &
         "  boundary warps ", 100.0_dp*plan%nseg*32.0_dp/max(real(plan%nwork, dp), 1.0_dp), "% of threads"
      call global_logger%info(trim(logline))
   end subroutine log_segment_shape


   !
   ! Count the quartets that survive screening, using the same index
   ! arithmetic and the same two tests as the kernel.  Diagnostic only: it
   ! exists so the effect of density screening is a measured number rather
   ! than an inference from wall-clock.
   !
   ! A REDUCTION WITHOUT A REDUCTION CLAUSE. `do concurrent` gained
   ! `reduce()` in Fortran 2023, but gfortran does not take it, and this used
   ! to be the one OpenACC compute construct in terco (a `parallel loop
   ! reduction`). Instead each thread counts one chunk of KEPT_CHUNK work
   ! items into its own slot and the host adds the slots up: the same count,
   ! exactly, since it is an integer. 4096 items a thread keeps the partial
   ! array to 2 MB at 2e9 work items and the launch at half a million
   ! threads, plenty to fill the GPU.
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

      integer(int64), parameter :: KEPT_CHUNK = 4096
      integer(int64), allocatable :: part(:)
      integer(int64) :: nchunk, ic

      nchunk = (nwork + KEPT_CHUNK - 1)/KEPT_CHUNK
      allocate (part(nchunk))
      !$acc enter data create(part)
      do concurrent(ic=1:nchunk)
         call kept_in_chunk(ic, KEPT_CHUNK, nchunk, nwork, nseg, sNB, sOA, sOB, sD, sOff, &
                            npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, part)
      end do
      !$acc update self(part)
      !$acc exit data delete(part)
      nkept = sum(part)
      deallocate (part)
   end subroutine count_kept

   pure subroutine kept_in_chunk(ic, chunk, nchunk, nwork, nseg, sNB, sOA, sOB, sD, sOff, &
                                 npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, part)
      !! Work items (ic-1)*chunk+1 .. ic*chunk of count_kept, into part(ic).
      !$acc routine seq
      integer(int64), intent(in) :: ic, chunk, nchunk, nwork
      integer,  intent(in) :: nseg, npair, nbas
      integer,  intent(in) :: sNB(nseg), sOA(nseg), sOB(nseg)
      logical,  intent(in) :: sD(nseg)
      integer(int64), intent(in) :: sOff(nseg + 1)
      integer,  intent(in) :: sp_i(npair), sp_j(npair)
      real(dp), intent(in) :: thresh
      real(dp), intent(in) :: sp_q(npair), dsh(nbas, nbas)
      integer(int64), intent(inout) :: part(nchunk)

      integer(int64) :: gt, nk
      integer :: lo, hi, mid, seg, t, iab, icd, si, sj, sk, sl
      real(dp) :: qcut

      nk = 0
      do gt = (ic - 1)*chunk + 1, min(ic*chunk, nwork)
         lo = 1
         hi = nseg
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
            si = sp_i(sOA(seg) + iab)
            sj = sp_j(sOA(seg) + iab)
            sk = sp_i(sOB(seg) + icd)
            sl = sp_j(sOB(seg) + icd)
            if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                         dsh(si, sk), dsh(si, sl), &
                         dsh(sj, sk), dsh(sj, sl)) > thresh) nk = nk + 1
         end if
      end do
      part(ic) = nk
   end subroutine kept_in_chunk

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
      real(dp) :: v(NCUM_NOSYM*NCUM_NOSYM, 0:1), g(NCUM_NOSYM*NCUM_NOSYM)
      integer  :: na_, nb_, nc_, nd_, ib_, id_, mu, nu, lam, sig, index
      integer  :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz
      integer  :: px, py, pz, qx, qy, qz, ja, jc
      real(dp) :: cb, cdc, vv
      integer  :: idens
      real(dp) :: abx, aby, abz, cdx, cdy, cdz
      integer  :: mui, nuj, lamk, sigl
      ! largest Cartesian block this LMAX allows: ((L+1)(L+2)/2)^4
      real(dp) :: vbuf((((LNOSYM + 1)*(LNOSYM + 2)/2)**4))
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
      na_ = cart_ncart_of(la); nb_ = cart_ncart_of(lb)
      nc_ = cart_ncart_of(lc); nd_ = cart_ncart_of(ld)
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
      index = 0
      do id_ = 0, nd_ - 1
         sig = sigl + id_
      do ic = 0, nc_ - 1
         lam = lamk + ic
      do ib_ = 0, nb_ - 1
         nu = nuj + ib_
      do ia = 0, na_ - 1
         mu = mui + ia
         index = index + 1
         vv = vbuf(index)

         call scatter_nosym(mu, nu, lam, sig, vv, dij, dkl, dpq, ndens, nao, dmat, jmat, kmat)
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
         call scatter_nosym(mu, nu, lam, sig, vv, dij, dkl, dpq, ndens, nao, dmat, jmat, kmat)
      end do
      end do
      end do
      end do
#endif
   end subroutine one_bin_item_nosym

   !
   ! One Cartesian integral (mu nu|lam sig) into J and K, every permutation
   ! the bra/ket pair symmetry allows, with no assumption that D is
   ! symmetric: D_mn and D_nm are read separately. `dij`, `dkl` and `dpq`
   ! say which shells and pairs differ, so a permutation that maps the
   ! quartet onto itself is not counted twice.
   !
   pure subroutine scatter_nosym(mu, nu, lam, sig, vv, dij, dkl, dpq, ndens, nao, dmat, jmat, kmat)
      !$acc routine seq
      integer,  intent(in)    :: mu, nu, lam, sig, ndens, nao
      real(dp), intent(in)    :: vv
      logical,  intent(in)    :: dij, dkl, dpq
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: jmat(ndens, nao, nao), kmat(ndens, nao, nao)

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
   end subroutine scatter_nosym


end module trc_binkernel
