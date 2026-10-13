!
! Two-electron gradient, over the same binned shell-pair plan the Fock
! build uses: one launch per angular-momentum class, Rys only -- there is
! no HGP gradient (out of scope).
!
! SCREENING: the same bound as the energy/Fock build (`thresh`), not a
! separate gradient-specific tolerance. Gamma already carries the density
! weighting per quartet exactly as the Fock digestion's jab/jcd/kac/kad/kbc/
! kbd blocks do (see scripts/gen_rys_grad.py), so a quartet the Schwarz
! bound already lets through at `thresh` is precisely the set whose energy
! contribution the SCF trusted; there is no argument for trusting its
! derivative less. A modest tightening factor was considered and rejected
! for lack of any measurement showing the plain bound insufficient for the
! ~1e-7 relative finite-difference tolerance -- see
! test/check_rys_grad_fd.F90.
!
module trc_grad
   use trc_boys, only: dp
   use trc_bins, only: pair_bins_t
#ifdef TERCO_RYS
   use trc_rys_grad_kernels, only: rysg_dispatch, rysg_has_class, &
                                   rysg_screen_dispatch, rysg_compact_dispatch
#endif
   implicit none
   private
   public :: grad_bins

   !: Same fixed radix as trc_rys_select's RYS_CLASS_RADIX -- duplicated as a
   !: literal rather than imported, for the same reason: importing would
   !: pull in a module this one has no other use for, for one integer.
   integer, parameter :: GCLASS_RADIX = 11

   !: Same prefactor as trc_pc_rysg0000's TWO_PI_2_5 -- duplicated for the
   !: same reason GCLASS_RADIX is: this diagnostic has no other use for a
   !: module that carries the whole class kernel.
   real(dp), parameter :: GDIAG_TWO_PI_2_5 = 34.986836655249725_dp

   !: Memory bound for ORDERED SURVIVOR COMPACTION's per-batch survivor-list
   !: scratch (out_seg/out_iab/out_icd: three 4-byte integer arrays, plus an
   !: 8-byte atomic counter). A class run is walked in segment-aligned
   !: batches of at most this many work items (`grad_class_compact` below),
   !: so that scratch never exceeds 3*4*GRAD_COMPACT_CAP bytes (~384 MiB at
   !: the default) regardless of how large a single class's own segment
   !: run is. 2**25 is a round number comfortably inside a GPU's memory for
   !: three such arrays; no class run on w150 6-31G* is anywhere near large
   !: enough to need more than one batch -- this exists so a larger system
   !: does not exhaust device memory, not because it has been exercised.
   integer(kind=8), parameter :: GRAD_COMPACT_CAP = 2_8**25

   !: LANE-PRIVATE GRADIENT COPIES: how many copies gp(3, natm, GRAD_NPRIV)
   !: of the gradient the compact pass adds into (survivor i into copy
   !: mod(i-1, GRAD_NPRIV)+1; see scripts/gen_rys_grad.py's compact
   !: driver). One per lane of a warp. Measured on (ss|ss)+(ps|ss), w150
   !: 6-31G* 1e-13, V100, compact pass: 3.11 s with one grad, 1.93 s with
   !: 8 copies, 1.91 s with 32, 1.91 s with 64. 32 copies are 345 KB at
   !: 450 atoms: L2-resident, zeroed and summed once per eri%gradient
   !: (the two tiny kernels are microseconds). Whole eri%gradient on w150:
   !: 17.9 -> 16.4 s at 1e-13, 9.29 -> 8.80 s at 1e-10, with classes of
   !: L > 3 left on copy 1 (gen_rys_grad.py).
   integer, parameter :: GRAD_NPRIV = 32

contains

   !
   ! Walk admitted bin pairs on the host (same plan fock_bins builds), sort
   ! by class, launch `rysg_dispatch` once per contiguous class run. No
   ! density-weighted bin screen: Gamma is not guaranteed to admit the same
   ! bin-level bound a symmetric Fock density does (same reasoning as
   ! `eri_fock_nosym`), so this always walks the plain Schwarz-screened plan.
   !
   subroutine grad_bins(b, nbas, nao, sh_l, ao_off, sh_at, natm, thresh, &
                        jfac, kfac, dsh, &
                        qp_bin_pbase, qp_bin_k, nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                        sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                        ndens, dmat, grad, rank, nranks)
      type(pair_bins_t), intent(in) :: b
      integer,  intent(in)    :: nbas, nao, natm
      integer,  intent(in)    :: sh_l(nbas), ao_off(nbas), sh_at(nbas)
      real(dp), intent(in)    :: thresh, jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: qp_bin_pbase(:), qp_bin_k(:)
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(b%npair), sq_ay(b%npair), sq_az(b%npair)
      real(dp), intent(in)    :: sq_abx(b%npair), sq_aby(b%npair), sq_abz(b%npair)
      integer,  intent(in)    :: sq_n(b%npair)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: grad(3, natm)
      integer,  intent(in)    :: rank, nranks

      integer :: ia, ib, ka, kb, smax_keep, nA, nB, nseg, is
      integer(kind=8) :: nt
      integer, allocatable :: sA(:), sOA(:), sOB(:), sNB(:), sPA(:), sPB(:)
      integer, allocatable :: sLA(:), sLB(:), sLC(:), sLD(:)
      logical, allocatable :: sD(:)
      integer(kind=8), allocatable :: sOff(:)
      integer(kind=8) :: nwork

      smax_keep = int(-log10(thresh))*b%sres

      nseg = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            nseg = nseg + 1
         end do
      end do

      allocate (sA(nseg), sOA(nseg), sOB(nseg), sNB(nseg), sD(nseg))
      allocate (sLA(nseg), sLB(nseg), sLC(nseg), sLD(nseg))
      allocate (sPA(nseg), sPB(nseg), sOff(nseg + 1))

      is = 0; sOff(1) = 0
      do ia = 1, b%nlive
         ka = b%live(ia)
         nA = b%bin_cnt(ka)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            nB = b%bin_cnt(kb)
            if (nA == 0 .or. nB == 0) cycle
            is = is + 1
            sA(is) = nA; sNB(is) = nB
            sOA(is) = b%bin_off(ka); sOB(is) = b%bin_off(kb)
            sPA(is) = qp_bin_pbase(ka); sPB(is) = qp_bin_pbase(kb)
            sD(is) = (ka == kb)
            sLA(is) = b%bin_la(ka); sLB(is) = b%bin_lb(ka)
            sLC(is) = b%bin_la(kb); sLD(is) = b%bin_lb(kb)
            if (ka == kb) then
               nt = int(nA, 8)*int(nA + 1, 8)/2
            else
               nt = int(nA, 8)*int(nB, 8)
            end if
            if (nt > 2147483000_8) error stop &
               "terco: bin-pair segment exceeds a default integer in grad_bins"
            sOff(is + 1) = sOff(is) + nt
         end do
      end do

      nwork = sOff(nseg + 1)
      if (nwork == 0) then
         deallocate (sA, sOA, sOB, sNB, sPA, sPB, sD, sOff, sLA, sLB, sLC, sLD)
         return
      end if

#ifdef TERCO_RYS
      block
         integer, allocatable :: ord(:), ckey(:)
         integer :: a2, b2, t2, c0, c1, key2, ga, gc, gs
         real(dp) :: gsum
         real(dp), allocatable :: gp(:, :, :)
         !! The compact pass's lane-private gradient copies (GRAD_NPRIV).
         !! A local, not a component of trc_eri_t or an argument of
         !! grad_bins: nvfortran 26.5 has hung or ICE'd on both.
         allocate (ord(nseg), ckey(nseg))
         do a2 = 1, nseg
            ckey(a2) = ((sLA(a2)*GCLASS_RADIX + sLB(a2))*GCLASS_RADIX &
                        + sLC(a2))*GCLASS_RADIX + sLD(a2)
            ord(a2) = a2
         end do
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
         allocate (gp(3, natm, GRAD_NPRIV))
         !$acc enter data create(gp)
         do concurrent(gs=1:GRAD_NPRIV, ga=1:natm, gc=1:3)
            gp(gc, ga, gs) = 0.0_dp
         end do
         c0 = 1
         do while (c0 <= nseg)
            c1 = c0
            do while (c1 < nseg)
               if (((sLA(c1 + 1)*GCLASS_RADIX + sLB(c1 + 1))*GCLASS_RADIX &
                    + sLC(c1 + 1))*GCLASS_RADIX + sLD(c1 + 1) /= &
                   ((sLA(c0)*GCLASS_RADIX + sLB(c0))*GCLASS_RADIX &
                    + sLC(c0))*GCLASS_RADIX + sLD(c0)) exit
               c1 = c1 + 1
            end do
            key2 = ((sLA(c0)*GCLASS_RADIX + sLB(c0))*GCLASS_RADIX &
                    + sLC(c0))*GCLASS_RADIX + sLD(c0)
            if (.not. rysg_has_class(key2)) error stop &
               "terco: grad_bins met a class TERCO_RYS did not generate a gradient kernel for"
            if (key2 == 0 .and. grad_diag_on()) &
               call grad_diag_counts(c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                     b%npair, b%sp_i, b%sp_j, b%sp_q, sq_n, thresh, thresh*1.0e-3_dp, &
                                     dsh, nbas, sh_at, nqp, qp_z, qp_c, qp_a)
            if (grad_fused_on()) then
               call rysg_dispatch(key2, c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                 b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*1.0e-3_dp, &
                                 jfac, kfac, dsh, nbas, nao, sh_l, ao_off, sh_at, natm, &
                                 nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                                 sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                 ndens, dmat, grad, rank, nranks)
            else
               call grad_class_compact(key2, c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                 b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*1.0e-3_dp, &
                                 jfac, kfac, dsh, nbas, nao, ao_off, sh_at, natm, &
                                 nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                                 sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                 ndens, dmat, gp, rank, nranks)
            end if
            c0 = c1 + 1
         end do
         ! Sum the copies into grad, on the device (grad is already there:
         ! eri_gradient's enter data). The fused path never touches gp, so
         ! this adds zeros under TRC_GRAD_FUSED=1.
         do concurrent(ga=1:natm, gc=1:3) local(gsum, gs)
            gsum = 0.0_dp
            do gs = 1, GRAD_NPRIV
               gsum = gsum + gp(gc, ga, gs)
            end do
            grad(gc, ga) = grad(gc, ga) + gsum
         end do
         !$acc exit data delete(gp)
         deallocate (gp)
         !$acc exit data delete(sA, sNB, sOA, sOB, sPA, sPB, sD, sOff)
         deallocate (ord, ckey)
      end block
#else
      error stop "terco: grad_bins called without TERCO_RYS -- there is no HGP gradient"
#endif
      deallocate (sA, sOA, sOB, sNB, sPA, sPB, sD, sOff, sLA, sLB, sLC, sLD)
   end subroutine grad_bins

   subroutine permute_segments(n, ord, sA, sNB, sOA, sOB, sPA, sPB, sD, sOff, sLA, sLB, sLC, sLD)
      !! Same permutation fock_bins applies to its own segment descriptors.
      integer, intent(in) :: n, ord(n)
      integer, intent(inout) :: sA(n), sNB(n), sOA(n), sOB(n), sPA(n), sPB(n)
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
      t = sPA(ord); sPA = t
      t = sPB(ord); sPB = t
      t = sLA(ord); sLA = t
      t = sLB(ord); sLB = t
      t = sLC(ord); sLC = t
      t = sLD(ord); sLD = t
      tl = sD(ord); sD = tl
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

   logical function grad_fused_on()
      !! TRC_GRAD_FUSED=1 runs the fused rysg{tag} kernel (decode, the
      !! three quartet-level screens and the integral build all in one
      !! `do concurrent`, as the generator always emitted before ORDERED
      !! SURVIVOR COMPACTION) instead of the default screen/compact split
      !! below -- kept for A/B timing, never for correctness: both paths
      !! run the same generated integral body over the same survivors.
      character(len=8) :: v
      call get_environment_variable("TRC_GRAD_FUSED", v)
      grad_fused_on = (trim(v) == "1")
   end function grad_fused_on

   !
   ! ORDERED SURVIVOR COMPACTION glue (the default path for one class run;
   ! TRC_GRAD_FUSED=1 takes grad_bins back to rysg_dispatch above). Screens
   ! class run [c0, c1] into a dense (seg, iab, icd) survivor list with
   ! `rysg_screen_dispatch`, then runs `rysg_compact_dispatch` over exactly
   ! that many survivors -- see scripts/gen_rys_grad.py's
   ! screen_and_compact_text for the technique the per-class kernels
   ! implement and scripts/gen_rys_grad_dispatch.py for the two
   ! dispatchers called here.
   !
   ! RANKS: an MPI rank takes a CONTIGUOUS slice of this class's own
   ! segment range [c0, c1] -- not the fused kernel's interleaved
   ! work-item stride (rysgi{tag}'s `gt = g0 + (i-1)*nranks`). Every
   ! segment is assigned to exactly one rank, so every work item still
   ! runs exactly once across all ranks and the summed gradient is
   ! unchanged; load is balanced by segment count rather than by work
   ! item, and a rank's own slice stays segment-contiguous so the screen
   ! pass's RUNLEN=32 chunks keep sharing a ket within it. This has not
   ! been measured against the fused stride's balance -- MPI is outside
   ! today's validation (TERCO_ENABLE_MPI=OFF); nranks=1 reduces this to
   ! the whole [c0, c1] range, unchanged from before.
   !
   ! BATCHING: a rank's own slice is then walked forward in segment-
   ! aligned batches of at most GRAD_COMPACT_CAP work items, bounding the
   ! survivor-list scratch's memory regardless of how large a single
   ! class's segment run is -- see GRAD_COMPACT_CAP's own comment.
   !
#ifdef TERCO_RYS
   subroutine grad_class_compact(key, c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                 npair, sp_i, sp_j, sp_q, thresh, pcut, jfac, kfac, &
                                 dsh, nbas, nao, ao_off, sh_at, natm, &
                                 nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                                 sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                 ndens, dmat, gp, rank, nranks)
      integer,  intent(in)    :: key, c0, c1, nseg, npair, nbas, nao, natm
      integer(kind=8), intent(in) :: sOff(nseg + 1)
      integer,  intent(in)    :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in)    :: sPA(nseg), sPB(nseg)
      logical,  intent(in)    :: sD(nseg)
      integer,  intent(in)    :: sp_i(npair), sp_j(npair), sq_n(npair)
      real(dp), intent(in)    :: sp_q(npair), thresh, pcut, jfac, kfac
      real(dp), intent(in)    :: dsh(nbas, nbas)
      integer,  intent(in)    :: ao_off(nbas), sh_at(nbas)
      integer,  intent(in)    :: nqp
      real(dp), intent(in)    :: qp_z(nqp), qp_c(nqp), qp_a(nqp)
      real(dp), intent(in)    :: qp_px(nqp), qp_py(nqp), qp_pz(nqp)
      real(dp), intent(in)    :: sq_ax(npair), sq_ay(npair), sq_az(npair)
      real(dp), intent(in)    :: sq_abx(npair), sq_aby(npair), sq_abz(npair)
      integer,  intent(in)    :: ndens
      real(dp), intent(in)    :: dmat(ndens, nao, nao)
      real(dp), intent(inout) :: gp(3, natm, GRAD_NPRIV)
      !! The lane-private gradient copies grad_bins sums into grad.
      integer,  intent(in)    :: rank, nranks

      integer :: nseg_run, seg_lo, seg_hi, b0, b1
      integer, allocatable :: out_seg(:), out_iab(:), out_icd(:)
      integer :: nmax
      integer(kind=8) :: g0, g1, nout
      integer(kind=8) :: ctr(1)
      !! A 1-element array, not a bare scalar: nvfortran's mem:separate
      !! present-table lookup for an atomic-capture target threaded through
      !! two subroutine layers hangs on a bare scalar here (same finding as
      !! the (ps|ss) prototype's own ctr).

      nseg_run = c1 - c0 + 1
      seg_lo = c0 + int((int(rank, 8)*int(nseg_run, 8))/int(nranks, 8))
      seg_hi = c0 + int((int(rank + 1, 8)*int(nseg_run, 8))/int(nranks, 8)) - 1
      if (seg_hi < seg_lo) return

      b0 = seg_lo
      do while (b0 <= seg_hi)
         b1 = b0
         do while (b1 < seg_hi)
            if (sOff(b1 + 2) - sOff(b0) > GRAD_COMPACT_CAP) exit
            b1 = b1 + 1
         end do

         g0 = sOff(b0) + 1; g1 = sOff(b1 + 1)
         nmax = 0
         if (g1 >= g0) nmax = int(g1 - g0 + 1)
         if (nmax > 0) then
            allocate (out_seg(nmax), out_iab(nmax), out_icd(nmax))
            ctr(1) = 0_8
            !$acc enter data create(out_seg, out_iab, out_icd) copyin(ctr)
            call rysg_screen_dispatch(key, b0, b1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                      npair, sp_i, sp_j, sp_q, thresh, dsh, nbas, sh_at, &
                                      out_seg, out_iab, out_icd, nmax, ctr)
            !$acc update host(ctr)
            nout = min(ctr(1), int(nmax, 8))
            call rysg_compact_dispatch(key, nout, out_seg, out_iab, out_icd, nseg, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                                      npair, sp_i, sp_j, sp_q, pcut, jfac, kfac, nbas, nao, ao_off, sh_at, natm, &
                                      nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                                      sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                                      ndens, dmat, GRAD_NPRIV, gp)
            !$acc exit data delete(out_seg, out_iab, out_icd, ctr)
            deallocate (out_seg, out_iab, out_icd)
         end if
         b0 = b1 + 1
      end do
   end subroutine grad_class_compact
#endif

   logical function grad_diag_on()
      !! TRC_GRAD_DIAG=1 prints, once, how many (ss|ss) work items `grad_bins`
      !! launches and how many pass each screen: exact Schwarz,
      !! density-weighted, the same-atom skip, and the primitive-quartet
      !! prescreen inside the bra/ket loop. Diagnostic only -- never on the
      !! path an integral takes.
      character(len=8) :: v
      call get_environment_variable("TRC_GRAD_DIAG", v)
      grad_diag_on = (trim(v) == "1")
   end function grad_diag_on

   !
   ! Mirrors rysgi0000's own screens, stage by stage, over one class's
   ! segment run [c0, c1] -- same decode, same tests, same constant -- so the
   ! counts describe exactly what that kernel does rather than an
   ! approximation of it. Diagnostic only: this is not part of grad_bins'
   ! normal control flow and changes no result; it is called once, behind
   ! TRC_GRAD_DIAG, right before the (ss|ss) dispatch it is measuring.
   !
   subroutine grad_diag_counts(c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                               npair, sp_i, sp_j, sp_q, sq_n, thresh, pcut, dsh, nbas, sh_at, &
                               nqp, qp_z, qp_c, qp_a)
      integer,  intent(in) :: c0, c1, nseg, npair, nbas
      integer(kind=8), intent(in) :: sOff(nseg + 1)
      integer,  intent(in) :: sA(nseg), sNB(nseg), sOA(nseg), sOB(nseg)
      integer,  intent(in) :: sPA(nseg), sPB(nseg)
      logical,  intent(in) :: sD(nseg)
      integer,  intent(in) :: sp_i(npair), sp_j(npair), sq_n(npair)
      real(dp), intent(in) :: sp_q(npair), thresh, pcut
      real(dp), intent(in) :: dsh(nbas, nbas)
      integer,  intent(in) :: sh_at(nbas)
      integer,  intent(in) :: nqp
      real(dp), intent(in) :: qp_z(nqp), qp_c(nqp), qp_a(nqp)

      integer(kind=8) :: gt, g0, g1, nw
      integer(kind=8) :: n_total, n_schwarz, n_density, n_atom
      integer(kind=8) :: n_bra_tot, n_bra_pass, n_pair_tot, n_pair_pass
      integer :: lo, hi, mid, seg, t, iab, icd, si, sj, sk, sl
      integer :: slotab, slotcd, nab, ncd, k1, k2, kp, kq
      integer(kind=8) :: nsa, u, kx
      real(dp) :: qcut, zeta, eta, zpe, bnd, pref
      integer :: which

      g0 = sOff(c0) + 1; g1 = sOff(c1 + 1)
      nw = 0
      if (g1 >= g0) nw = g1 - g0 + 1
      n_total = nw; n_schwarz = 0; n_density = 0; n_atom = 0
      n_bra_tot = 0; n_bra_pass = 0; n_pair_tot = 0; n_pair_pass = 0

      ! Seven passes, each with exactly ONE OpenACC reduction variable --
      ! the shape count_kept (trc_binkernel.F90) already trusts in
      ! production. A single pass with several reduction variables at once
      ! gave inconsistent, run-to-run-varying counts here (density came back
      ! LARGER than Schwarz, which the control flow makes impossible), so
      ! this stays one reduction per launch rather than chase that further --
      ! diagnostic-only code, correctness over launch count.
      do which = 1, 7
         block
            integer(kind=8) :: acc
            acc = 0
            !$acc parallel loop default(present) reduction(+:acc) &
            !$acc   private(lo, hi, mid, seg, t, iab, icd, si, sj, sk, sl, qcut, &
            !$acc           slotab, slotcd, nab, ncd, k1, k2, kp, kq, nsa, u, kx, &
            !$acc           zeta, eta, zpe, bnd, pref)
            do gt = g0, g1
               lo = c0; hi = c1
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
                  nsa = int(sA(seg), kind(nsa))
                  u = int(t - 1, kind(u))
                  kx = int((real(2*nsa + 1, dp) - sqrt(real(2*nsa + 1, dp)**2 - 8.0_dp*real(u, dp)))/2.0_dp, kind(kx))
                  if (kx < 0) kx = 0
                  do while (kx > 0)
                     if ((kx*(2*nsa + 1) - kx*kx)/2 <= u) exit
                     kx = kx - 1
                  end do
                  do while (((kx + 1)*(2*nsa + 1) - (kx + 1)*(kx + 1))/2 <= u)
                     kx = kx + 1
                  end do
                  icd = int(kx) + 1
                  iab = icd + int(u - (kx*(2*nsa + 1) - kx*kx)/2)
               else
                  icd = (t - 1)/sA(seg) + 1
                  iab = t - (icd - 1)*sA(seg)
               end if

               qcut = sp_q(sOA(seg) + iab)*sp_q(sOB(seg) + icd)
               if (qcut <= thresh) cycle
               if (which == 1) then
                  acc = acc + 1; cycle
               end if

               si = sp_i(sOA(seg) + iab); sj = sp_j(sOA(seg) + iab)
               sk = sp_i(sOB(seg) + icd); sl = sp_j(sOB(seg) + icd)

               if (qcut*max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl), &
                            dsh(si, sk), dsh(si, sl), &
                            dsh(sj, sk), dsh(sj, sl)) <= thresh) cycle
               if (which == 2) then
                  acc = acc + 1; cycle
               end if

               if (sh_at(si) == sh_at(sj) .and. sh_at(si) == sh_at(sk) &
                   .and. sh_at(si) == sh_at(sl)) cycle
               if (which == 3) then
                  acc = acc + 1; cycle
               end if

               slotab = sOA(seg) + iab; slotcd = sOB(seg) + icd
               nab = sq_n(slotab); ncd = sq_n(slotcd)
               if (which == 4) then
                  acc = acc + nab; cycle
               end if
               if (which == 6) then
                  acc = acc + int(nab, kind(acc))*int(ncd, kind(acc)); cycle
               end if

               do k1 = 1, nab
                  kp = sPA(seg) + (k1 - 1)*sA(seg) + iab
                  zeta = qp_z(kp)
                  bnd = GDIAG_TWO_PI_2_5*abs(qp_c(kp))/(zeta*sqrt(zeta))
                  if (bnd*abs(qp_c(sPB(seg) + icd))/qp_z(sPB(seg) + icd) <= pcut) cycle
                  if (which == 5) then
                     acc = acc + 1; cycle
                  end if
                  do k2 = 1, ncd
                     kq = sPB(seg) + (k2 - 1)*sNB(seg) + icd
                     eta = qp_z(kq)
                     zpe = zeta + eta
                     pref = GDIAG_TWO_PI_2_5/(zeta*eta*sqrt(zpe))*qp_c(kp)*qp_c(kq)
                     if (bnd*abs(qp_c(kq))/eta <= pcut) exit
                     if (abs(pref) <= pcut) cycle
                     if (which == 7) acc = acc + 1
                  end do
               end do
            end do
            select case (which)
            case (1); n_schwarz = acc
            case (2); n_density = acc
            case (3); n_atom = acc
            case (4); n_bra_tot = acc
            case (5); n_bra_pass = acc
            case (6); n_pair_tot = acc
            case (7); n_pair_pass = acc
            end select
         end block
      end do

      print '(A,I0,A,I0,A,I0,A,I0)', &
         "TRC_GRAD_DIAG (ss|ss): total=", n_total, " schwarz=", n_schwarz, &
         " density=", n_density, " not_same_atom=", n_atom
      print '(A,I0,A,I0,A,I0,A,I0)', &
         "TRC_GRAD_DIAG (ss|ss) primitive pairs: bra_tot=", n_bra_tot, &
         " bra_pass=", n_bra_pass, " pair_tot=", n_pair_tot, " pair_pass=", n_pair_pass
   end subroutine grad_diag_counts

end module trc_grad
