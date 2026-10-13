!
! Two-electron gradient, over the same binned shell-pair plan the Fock
! build uses: one launch per angular-momentum class, Rys only -- there is
! no HGP gradient (out of scope; see docs/RYS_PLAN.md and the task plan
! this follows).
!
! SCREENING: the same bound as the energy/Fock build (`thresh`), not a
! separate gradient-specific tolerance. Gamma already carries the density
! weighting per quartet exactly as the Fock digestion's jab/jcd/kac/kad/kbc/
! kbd blocks do (see scripts/gen_rys_grad.py), so a quartet the Schwarz
! bound already lets through at `thresh` is precisely the set whose energy
! contribution the SCF trusted; there is no argument for trusting its
! derivative less. A modest tightening factor was considered (the task plan
! asks this to be stated) and rejected for lack of any measurement showing
! the plain bound insufficient for the ~1e-7 relative finite-difference
! gate -- see test/check_rys_gradient_fd.F90.
!
module trc_grad
   use trc_boys, only: dp
   use trc_bins, only: pair_bins_t
#ifdef TERCO_RYS
   use trc_rys_grad_kernels, only: rysg_dispatch, rysg_has_class
#endif
   implicit none
   private
   public :: grad_bins

   !: Same fixed radix as trc_rys_select's RYS_CLASS_RADIX -- duplicated as a
   !: literal rather than imported, for the same reason: importing would
   !: pull in a module this one has no other use for, for one integer.
   integer, parameter :: GCLASS_RADIX = 11

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
         integer :: a2, b2, t2, c0, c1, key2
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
            call rysg_dispatch(key2, c0, c1, nseg, sOff, sA, sNB, sOA, sOB, sPA, sPB, sD, &
                              b%npair, b%sp_i, b%sp_j, b%sp_q, thresh, thresh*1.0e-3_dp, &
                              jfac, kfac, dsh, nbas, nao, sh_l, ao_off, sh_at, natm, &
                              nqp, qp_z, qp_c, qp_a, qp_px, qp_py, qp_pz, &
                              sq_ax, sq_ay, sq_az, sq_abx, sq_aby, sq_abz, sq_n, &
                              ndens, dmat, grad, rank, nranks)
            c0 = c1 + 1
         end do
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

end module trc_grad
