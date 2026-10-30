!
! Two-electron gradient, over the same binned shell-pair plan the Fock
! build uses: one launch per angular-momentum class, Rys only -- there is
! no HGP gradient (out of scope).
!
! RANGE SEPARATION: `omega` switches the operator to erf(omega*r12)/r12,
! same meaning and same rho -> rho_theta substitution as eri_fock's own
! (scripts/gen_perclass.py) -- threaded only into the COMPACT pass
! (grad_classes/rysg_compact_dispatch), never the screen pass: the
! Schwarz bound stays the plain Coulomb one, which bounds the attenuated
! integral too (eri_fock's own comment). omega == 0 is bit-identical to
! the unattenuated path -- see scripts/gen_rys_grad.py's `_RANGE_SEP_BLOCK`.
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
! test/check_rys_grad_fd.F90. Measured on the 450-atom water cluster,
! 6-31G*, same converged 1e-13 density: the gradient at thresh = 1e-10
! differs from the one at 1e-13 by at most 3.9e-7 Eh/Bohr (rms 8.5e-8,
! against components up to 92), about 40x below a tight geometry
! optimisation's max-force criterion (1.5e-5).
!
module trc_grad
   use pic_types, only: dp, int64
   use trc_bins, only: pair_bins_t, trc_bin_plan_t
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t, trc_bin_launch_t
   use trc_rys_grad_kernels, only: rysg_has_class, rysg_screen_dispatch, &
                                   rysg_compact_dispatch
   implicit none
   private
   public :: grad_bins

   integer, parameter :: GCLASS_RADIX = 11
   !! Same fixed radix as trc_rys_select's RYS_CLASS_RADIX -- duplicated as a
   !! literal rather than imported, for the same reason: importing would
   !! pull in a module this one has no other use for, for one integer.

   integer(int64), parameter :: GRAD_COMPACT_CAP = 2_int64**25
   !! Memory bound for ORDERED SURVIVOR COMPACTION's per-batch survivor-list
   !! scratch (out_seg/out_iab/out_icd: three 4-byte integer arrays, plus an
   !! 8-byte atomic counter). A class run is walked in segment-aligned
   !! batches of at most this many work items (`grad_classes` below),
   !! and a single segment larger than this is itself walked in windows of
   !! this many, so that scratch never exceeds 3*4*GRAD_COMPACT_CAP bytes
   !! (~384 MiB at the default) however large a class run or a segment is.
   !! 2**25 is a round number comfortably inside a GPU's memory for three
   !! such arrays; no class run on w150 6-31G* is anywhere near large enough
   !! to need more than one batch -- this exists so a larger system does
   !! not exhaust device memory.

   integer, parameter :: GRAD_NPRIV = 32
   !! LANE-PRIVATE GRADIENT COPIES: how many copies gp(3, natm, GRAD_NPRIV)
   !! of the gradient the compact pass adds into (survivor i into copy
   !! mod(i-1, GRAD_NPRIV)+1; see scripts/gen_rys_grad.py's compact
   !! driver). One per lane of a warp. Measured on (ss|ss)+(ps|ss), w150
   !! 6-31G* 1e-13, V100, compact pass: 3.11 s with one grad, 1.93 s with
   !! 8 copies, 1.91 s with 32, 1.91 s with 64. 32 copies are 345 KB at
   !! 450 atoms: L2-resident, zeroed and summed once per eri%gradient
   !! (the two tiny kernels are microseconds). Whole eri%gradient on w150:
   !! 17.9 -> 16.4 s at 1e-13, 9.29 -> 8.80 s at 1e-10, with classes of
   !! L > 3 left on copy 1 (gen_rys_grad.py).

contains

   !
   ! Walk admitted bin pairs on the host (same plan binkernel_fock_bins builds), sort
   ! by class, run `grad_classes` over the contiguous class runs. No
   ! density-weighted bin screen: Gamma is not guaranteed to admit the same
   ! bin-level bound a symmetric Fock density does (same reasoning as
   ! `eri_fock_nosym`), so this always walks the plain Schwarz-screened plan.
   !
   subroutine grad_bins(b, pairs, quartets, launch, dsh, dmat, grad)
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      !! The shells (ao_off, sh_at, natm) and their pairs, device-resident.
      type(trc_pair_quartets_t), intent(in) :: quartets
      !! The bin-interleaved pairs, `a` included, device-resident.
      type(trc_bin_launch_t), intent(in) :: launch
      !! thresh, jfac, kfac, omega, ranks and the density's shape. omega is
      !! the range-separation parameter, bohr^-1, as eri_fock's; zero is the
      !! ordinary Coulomb operator the screen is already built for -- see
      !! grad_classes for why the SCREEN pass takes no omega of its
      !! own.
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(inout) :: grad(3, pairs%natm)

      type(trc_bin_plan_t) :: plan
      integer :: ga, gc, gs, natm
      real(dp) :: gsum
      real(dp), allocatable :: gp(:, :, :)
      !! The compact pass's lane-private gradient copies (GRAD_NPRIV).
      !! A local, not a component of trc_eri_t or an argument of
      !! grad_bins: nvfortran 26.5 has hung or ICE'd on both.

      natm = pairs%natm
      call plan%build(b, launch%thresh, quartets%bin_pbase, .false.)
      if (plan%nwork == 0) then
         call plan%release()
         return
      end if
      call plan%sort_by_class(GCLASS_RADIX)

      !$acc enter data copyin(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, plan%sOff)
      allocate (gp(3, natm, GRAD_NPRIV))
      !$acc enter data create(gp)
      do concurrent(gs=1:GRAD_NPRIV, ga=1:natm, gc=1:3)
         gp(gc, ga, gs) = 0.0_dp
      end do
      call grad_classes(plan, b, pairs, quartets, launch, dsh, dmat, gp)
      ! Sum the copies into grad, on the device (grad is already there:
      ! eri_gradient's enter data).
      do concurrent(ga=1:natm, gc=1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, GRAD_NPRIV
            gsum = gsum + gp(gc, ga, gs)
         end do
         grad(gc, ga) = grad(gc, ga) + gsum
      end do
      !$acc exit data delete(gp)
      deallocate (gp)
      !$acc exit data delete(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, plan%sOff)
      call plan%release()
   end subroutine grad_bins

   !
   ! ORDERED SURVIVOR COMPACTION glue, one class run [c0, c1] of the
   ! class-sorted plan at a time. Screens the run into a dense (seg, iab, icd) survivor list with
   ! `rysg_screen_dispatch`, then runs `rysg_compact_dispatch` over exactly
   ! that many survivors -- see scripts/gen_rys_grad.py's
   ! screen_and_compact_text for the technique the per-class kernels
   ! implement and scripts/gen_rys_grad_dispatch.py for the two
   ! dispatchers called here.
   !
   ! RANKS: an MPI rank takes a CONTIGUOUS slice of this class's own
   ! segment range [c0, c1] rather than an interleaved work-item stride.
   ! Every
   ! segment is assigned to exactly one rank, so every work item still
   ! runs exactly once across all ranks and the summed gradient is
   ! unchanged; load is balanced by segment count rather than by work
   ! item, and a rank's own slice stays segment-contiguous so the screen
   ! pass's RUNLEN=32 chunks keep sharing a ket within it. Its balance has
   ! not been measured -- MPI is outside today's validation
   ! (TERCO_ENABLE_MPI=OFF); nranks=1 reduces this to the whole [c0, c1]
   ! range.
   !
   ! BATCHING: a rank's own slice is then walked forward in segment-
   ! aligned batches of at most GRAD_COMPACT_CAP work items; a segment
   ! larger than that on its own is walked in windows of GRAD_COMPACT_CAP
   ! work items, which the screen pass takes as an explicit [w0, w1]. This
   ! bounds the survivor-list scratch's memory however large a class run
   ! or a single segment is -- see GRAD_COMPACT_CAP's own comment.
   !
   ! OMEGA reaches only `rysg_compact_dispatch` below, never
   ! `rysg_screen_dispatch`: the screen's Schwarz bound is the plain Coulomb
   ! one and stays valid at omega > 0 (pointwise larger -- see eri_fock's
   ! own comment), so the screen pass takes no omega of its own.
   !
   subroutine grad_classes(plan, b, pairs, quartets, launch, dsh, dmat, gp)
      type(trc_bin_plan_t), intent(in) :: plan
      !! Sorted by class, its arrays already on the device.
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      type(trc_pair_quartets_t), intent(in) :: quartets
      type(trc_bin_launch_t), intent(in) :: launch
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(inout) :: gp(3, pairs%natm, GRAD_NPRIV)
      !! The lane-private gradient copies grad_bins sums into grad.

      integer :: c0, c1, key, nseg, nseg_run, seg_lo, seg_hi, b0, b1
      integer, allocatable :: out_seg(:), out_iab(:), out_icd(:)
      integer :: nmax
      integer(int64) :: w0, w1, gend, nout
      integer(int64) :: ctr(1)
      !! A 1-element array, not a bare scalar: nvfortran's mem:separate
      !! present-table lookup for an atomic-capture target threaded through
      !! two subroutine layers hangs on a bare scalar here (same finding as
      !! the (ps|ss) prototype's own ctr).

      nseg = plan%nseg
      c0 = 1
      do while (c0 <= nseg)
         key = plan%class_key(c0, GCLASS_RADIX)
         c1 = c0
         do while (c1 < nseg)
            if (plan%class_key(c1 + 1, GCLASS_RADIX) /= key) exit
            c1 = c1 + 1
         end do
         if (.not. rysg_has_class(key)) error stop &
            "terco: grad_bins met a class the build did not generate a gradient kernel for"
         nseg_run = c1 - c0 + 1
         seg_lo = c0 + int((int(launch%rank, int64)*int(nseg_run, int64))/int(launch%nranks, int64))
         seg_hi = c0 + int((int(launch%rank + 1, int64)*int(nseg_run, int64))/int(launch%nranks, int64)) - 1

         b0 = seg_lo
         do while (b0 <= seg_hi)
            b1 = b0
            do while (b1 < seg_hi)
               if (plan%sOff(b1 + 2) - plan%sOff(b0) > GRAD_COMPACT_CAP) exit
               b1 = b1 + 1
            end do

            ! Batches of more than one segment are at most GRAD_COMPACT_CAP
            ! work items by construction, so this inner loop runs more than
            ! once only for a single segment larger than that.
            w0 = plan%sOff(b0) + 1
            gend = plan%sOff(b1 + 1)
            do while (w0 <= gend)
               w1 = min(gend, w0 + GRAD_COMPACT_CAP - 1)
               nmax = int(w1 - w0 + 1)
               allocate (out_seg(nmax), out_iab(nmax), out_icd(nmax))
               ctr(1) = 0_int64
               !$acc enter data create(out_seg, out_iab, out_icd) copyin(ctr)
               call rysg_screen_dispatch(key, b0, b1, w0, w1, nseg, plan%sOff, plan%sA, plan%sNB, &
                                         plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, &
                                         b%npair, b%sp_i, b%sp_j, b%sp_q, launch%thresh, dsh, pairs%nbas, &
                                         pairs%sh_at, out_seg, out_iab, out_icd, nmax, ctr)
               !$acc update host(ctr)
               nout = min(ctr(1), int(nmax, int64))
               call rysg_compact_dispatch(key, nout, out_seg, out_iab, out_icd, nseg, plan%sA, plan%sNB, &
                                          plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, &
                                          b%npair, b%sp_i, b%sp_j, b%sp_q, launch%thresh*1.0e-3_dp, &
                                          launch%jfac, launch%kfac, launch%omega, pairs%nbas, launch%nao, &
                                          pairs%ao_off, pairs%sh_at, pairs%natm, &
                                          quartets%nqp, quartets%z, quartets%c, quartets%a, &
                                          quartets%px, quartets%py, quartets%pz, &
                                          quartets%ax, quartets%ay, quartets%az, &
                                          quartets%abx, quartets%aby, quartets%abz, quartets%n, &
                                          launch%ndens, dmat, GRAD_NPRIV, gp)
               !$acc exit data delete(out_seg, out_iab, out_icd, ctr)
               deallocate (out_seg, out_iab, out_icd)
               w0 = w1 + 1
            end do
            b0 = b1 + 1
         end do
         c0 = c1 + 1
      end do
   end subroutine grad_classes

end module trc_grad
