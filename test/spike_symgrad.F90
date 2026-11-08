!
! SKELETON GRADIENT: design prototype (docs/SYMMETRY_GRADIENT.md, S16).
! Registered DISABLED; not part of the suite. CPU, plain Fortran.
!
! What it answers, with numbers, for one molecule:
!
!  1. COUNTS. The work items (quartets) the gradient's own plan launches
!     today -- plan%build(..., use_density=.false.), what grad_bins calls --
!     against the skeleton plan (plan%build_skeleton), and the class runs and
!     kernel launches each gives. Also the survivors of the gradient's
!     in-kernel screen, which no counter in src reports (there is no
!     TRC_KEPT_STATS for the gradient).
!  2. THE CLASS FILTER. How many skeleton segments carry (lc, ld) > (la, lb)
!     (the bra/ket swap bins_skel_cut_first makes), split by whether the
!     class is one src/CMakeLists.txt builds a gradient kernel for.
!  3. EXACTNESS. The symmetrized skeleton gradient against the full one.
!     There is no weighted gradient kernel, so the weight is applied by
!     REPLICATION: a segment of weight w is emitted nint(2*w) times (every
!     weight is a multiple of one half, so 2*w is an integer) and the sum is
!     halved. That costs about twice the full build's work and is useless as
!     an implementation, but it is an exact stand-in for the per-survivor
!     weight the design proposes, with no kernel change.
!
! usage: spike_symgrad <xyz> <basis.json> <sph|cart> <mode> [thresh]
!
!   mode = plan     plan counts only (no integrals, no screen pass): work
!                   items and segments of the full plan, of the skeleton plan,
!                   and of a skeleton plan with the bra/ket swap forbidden
!   mode = count    plan, plus the screen pass's survivor counts and launches
!   mode = verify   count, plus the skeleton gradient run through the
!                   PRODUCTION whole-bin stride -- i.e. the wrong answer, kept
!                   because it is the evidence that the gradient kernels need
!                   sSA/sSB (section 2 of the design)
!   mode = verify2  the exact oracle: one segment per quartet over a
!                   stride-one primitive layout, the weight by repetition,
!                   symmetrized against the full gradient
!
! The group projection of the density here is hand-rolled from ao_map/ao_sign
! and is NOT invariant in the spherical AO space (it is, to 6e-17, in the
! Cartesian one), so run verify/verify2 with `cart`. The production
! symmetrizer (trc_sym_scf_t%symmetrize, which check_sym_fock uses) has no
! such limitation.
!
module spike_symgrad_mod
   use pic_types, only: dp, int64
   use trc_bins, only: pair_bins_t, trc_bin_plan_t
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t, trc_bin_launch_t
   use trc_rys_grad_kernels, only: rysg_has_class, rysg_screen_dispatch, rysg_compact_dispatch
   implicit none

   integer, parameter :: GCLASS_RADIX = 11
   !! trc_grad's own.
   integer(int64), parameter :: GRAD_COMPACT_CAP = 2_int64**25
   !! trc_grad's own: the survivor-list batch bound.
   integer, parameter :: GRAD_NPRIV = 32
   !! trc_grad's own: the lane-private gradient copies.

contains

   !
   ! trc_grad's grad_classes, with the counters it does not keep: screen and
   ! compact launches, and survivors. `run` false screens only (no integrals),
   ! which is what the count mode wants.
   !
   subroutine spike_grad_classes(plan, b, pairs, quartets, launch, dsh, dmat, gp, run, nl, nsurv)
      type(trc_bin_plan_t), intent(in) :: plan
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      type(trc_pair_quartets_t), intent(in) :: quartets
      type(trc_bin_launch_t), intent(in) :: launch
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(inout) :: gp(3, pairs%natm, GRAD_NPRIV)
      logical, intent(in) :: run
      integer, intent(out) :: nl
      integer(int64), intent(out) :: nsurv

      integer :: c0, c1, key, nseg, b0, b1
      integer, allocatable :: out_seg(:), out_iab(:), out_icd(:)
      integer :: nmax
      integer(int64) :: w0, w1, gend, nout
      integer(int64) :: ctr(1)

      nseg = plan%nseg
      nl = 0
      nsurv = 0_int64
      c0 = 1
      do while (c0 <= nseg)
         key = plan%class_key(c0, GCLASS_RADIX)
         c1 = c0
         do while (c1 < nseg)
            if (plan%class_key(c1 + 1, GCLASS_RADIX) /= key) exit
            c1 = c1 + 1
         end do
         if (.not. rysg_has_class(key)) then
            print '(a,i0)', '  *** no gradient kernel for class key ', key
            c0 = c1 + 1
            cycle
         end if
         b0 = c0
         do while (b0 <= c1)
            b1 = b0
            do while (b1 < c1)
               if (plan%sOff(b1 + 2) - plan%sOff(b0) > GRAD_COMPACT_CAP) exit
               b1 = b1 + 1
            end do
            w0 = plan%sOff(b0) + 1
            gend = plan%sOff(b1 + 1)
            do while (w0 <= gend)
               w1 = min(gend, w0 + GRAD_COMPACT_CAP - 1)
               nmax = int(w1 - w0 + 1)
               allocate (out_seg(nmax), out_iab(nmax), out_icd(nmax))
               ctr(1) = 0_int64
               call rysg_screen_dispatch(key, b0, b1, w0, w1, nseg, plan%sOff, plan%sA, plan%sNB, &
                                         plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, &
                                         b%npair, b%sp_i, b%sp_j, b%sp_q, launch%thresh, dsh, pairs%nbas, &
                                         pairs%sh_at, out_seg, out_iab, out_icd, nmax, ctr)
               nout = min(ctr(1), int(nmax, int64))
               nsurv = nsurv + nout
               nl = nl + 1
               if (run) then
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
                  nl = nl + 1
               end if
               deallocate (out_seg, out_iab, out_icd)
               w0 = w1 + 1
            end do
            b0 = b1 + 1
         end do
         c0 = c1 + 1
      end do
   end subroutine spike_grad_classes

   subroutine spike_grad_run(plan, b, pairs, quartets, launch, dsh, dmat, grad, run, nl, nsurv)
      !! trc_grad's grad_bins over a plan the caller built.
      type(trc_bin_plan_t), intent(in) :: plan
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      type(trc_pair_quartets_t), intent(in) :: quartets
      type(trc_bin_launch_t), intent(in) :: launch
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      real(dp), intent(out)   :: grad(3, pairs%natm)
      logical, intent(in) :: run
      integer, intent(out) :: nl
      integer(int64), intent(out) :: nsurv
      real(dp), allocatable :: gp(:, :, :)
      integer :: ga, gc, gs

      grad = 0.0_dp
      allocate (gp(3, pairs%natm, GRAD_NPRIV))
      gp = 0.0_dp
      call spike_grad_classes(plan, b, pairs, quartets, launch, dsh, dmat, gp, run, nl, nsurv)
      do ga = 1, pairs%natm
         do gc = 1, 3
            do gs = 1, GRAD_NPRIV
               grad(gc, ga) = grad(gc, ga) + gp(gc, ga, gs)
            end do
         end do
      end do
      deallocate (gp)
   end subroutine spike_grad_run

   subroutine replicate_plan(src, dst)
      !! Every segment of `src` repeated nint(2*sW) times, with no weights:
      !! the unweighted kernels then accumulate 2*w times the segment's
      !! contribution, so the sum is twice the weighted (skeleton) one.
      type(trc_bin_plan_t), intent(in) :: src
      type(trc_bin_plan_t), intent(out) :: dst
      integer :: s, r, m, n, k
      integer(int64) :: nt

      n = 0
      do s = 1, src%nseg
         n = n + max(1, nint(2.0_dp*src%sW(s)))
      end do
      allocate (dst%sA(n), dst%sNB(n), dst%sOA(n), dst%sOB(n))
      allocate (dst%sPA(n), dst%sPB(n), dst%sSA(n), dst%sSB(n), dst%sD(n))
      allocate (dst%sLA(n), dst%sLB(n), dst%sLC(n), dst%sLD(n), dst%sM(n))
      allocate (dst%sOff(n + 1))
      dst%sOff(1) = 0
      k = 0
      do s = 1, src%nseg
         m = max(1, nint(2.0_dp*src%sW(s)))
         do r = 1, m
            k = k + 1
            dst%sA(k) = src%sA(s); dst%sNB(k) = src%sNB(s)
            dst%sOA(k) = src%sOA(s); dst%sOB(k) = src%sOB(s)
            dst%sPA(k) = src%sPA(s); dst%sPB(k) = src%sPB(s)
            dst%sSA(k) = src%sSA(s); dst%sSB(k) = src%sSB(s)
            dst%sD(k) = src%sD(s)
            dst%sLA(k) = src%sLA(s); dst%sLB(k) = src%sLB(s)
            dst%sLC(k) = src%sLC(s); dst%sLD(k) = src%sLD(s)
            dst%sM(k) = src%sM(s)
            nt = src%sOff(s + 1) - src%sOff(s)
            dst%sOff(k + 1) = dst%sOff(k) + nt
         end do
      end do
      dst%nseg = n
      dst%nwork = dst%sOff(n + 1)
   end subroutine replicate_plan


   !
   ! ONE SEGMENT PER QUARTET, with the primitive data re-laid out so that each
   ! pair's primitives are contiguous (stride one). This exists because the
   ! GRADIENT KERNELS TAKE NO sSA/sSB: they index the bin-interleaved tables
   ! as sPA(seg) + (k-1)*sA(seg) + iab, so a segment that covers only part of
   ! a bin -- which every skeleton segment does -- reads the wrong primitives.
   ! A segment of one pair against one pair has sA = sNB = 1, which is the
   ! stride of the re-laid-out table, so the kernels are exact on it. Only a
   ! prototype: it is one segment per quartet.
   !
   subroutine relayout_pairs(b, q, base, nqp, z, c, a, px, py, pz)
      type(pair_bins_t), intent(in) :: b
      type(trc_pair_quartets_t), intent(in) :: q
      integer, allocatable, intent(out) :: base(:)
      integer, intent(out) :: nqp
      real(dp), allocatable, intent(out) :: z(:), c(:), a(:), px(:), py(:), pz(:)
      integer :: k, p, kk, n, src, dst
      allocate (base(b%npair))
      nqp = 0
      do k = 1, b%nbin
         do p = b%bin_off(k) + 1, b%bin_off(k) + b%bin_cnt(k)
            base(p) = nqp
            nqp = nqp + q%n(p)
         end do
      end do
      allocate (z(nqp), c(nqp), a(nqp), px(nqp), py(nqp), pz(nqp))
      do k = 1, b%nbin
         n = b%bin_cnt(k)
         do p = b%bin_off(k) + 1, b%bin_off(k) + n
            do kk = 1, q%n(p)
               src = q%bin_pbase(k) + (kk - 1)*n + (p - b%bin_off(k))
               dst = base(p) + kk
               z(dst) = q%z(src); c(dst) = q%c(src); a(dst) = q%a(src)
               px(dst) = q%px(src); py(dst) = q%py(src); pz(dst) = q%pz(src)
            end do
         end do
      end do
   end subroutine relayout_pairs

   subroutine expand_plan(src, base, dst, rep)
      !! Every work item of `src` as its own segment of one bra pair against
      !! one ket pair, with stride-one primitive offsets from `base`. `rep`
      !! repeats each quartet nint(2*sW) times (the weight, so the caller
      !! halves the sum); without it every quartet appears once.
      type(trc_bin_plan_t), intent(in) :: src
      integer, intent(in) :: base(:)
      type(trc_bin_plan_t), intent(out) :: dst
      logical, intent(in) :: rep
      integer :: s, iab, icd, m, r, k, icd_hi
      integer(int64) :: n

      n = 0_int64
      do s = 1, src%nseg
         m = 1
         if (rep .and. allocated(src%sW)) m = max(1, nint(2.0_dp*src%sW(s)))
         n = n + (src%sOff(s + 1) - src%sOff(s))*int(m, int64)
      end do
      if (n > 2000000000_int64) then
         print '(a,i0)', 'expand_plan: too many segments ', n
         stop 1
      end if
      k = int(n)
      allocate (dst%sA(k), dst%sNB(k), dst%sOA(k), dst%sOB(k))
      allocate (dst%sPA(k), dst%sPB(k), dst%sSA(k), dst%sSB(k), dst%sD(k))
      allocate (dst%sLA(k), dst%sLB(k), dst%sLC(k), dst%sLD(k), dst%sM(k))
      allocate (dst%sOff(k + 1))
      dst%sOff(1) = 0
      k = 0
      do s = 1, src%nseg
         m = 1
         if (rep .and. allocated(src%sW)) m = max(1, nint(2.0_dp*src%sW(s)))
         do iab = 1, src%sA(s)
            icd_hi = src%sNB(s)
            if (src%sD(s)) icd_hi = iab
            do icd = 1, icd_hi
               do r = 1, m
                  k = k + 1
                  dst%sA(k) = 1; dst%sNB(k) = 1
                  dst%sOA(k) = src%sOA(s) + iab - 1
                  dst%sOB(k) = src%sOB(s) + icd - 1
                  dst%sPA(k) = base(src%sOA(s) + iab)
                  dst%sPB(k) = base(src%sOB(s) + icd)
                  dst%sSA(k) = 1; dst%sSB(k) = 1
                  dst%sD(k) = src%sD(s) .and. (iab == icd)
                  dst%sLA(k) = src%sLA(s); dst%sLB(k) = src%sLB(s)
                  dst%sLC(k) = src%sLC(s); dst%sLD(k) = src%sLD(s)
                  dst%sM(k) = src%sM(s)
                  dst%sOff(k + 1) = dst%sOff(k) + 1
               end do
            end do
         end do
      end do
      dst%nseg = k
      dst%nwork = dst%sOff(k + 1)
   end subroutine expand_plan

end module spike_symgrad_mod

program spike_symgrad
   use pic_types, only: dp, int64
   use spike_symgrad_mod
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_bins, only: trc_bin_plan_t
   use trc_eri, only: trc_eri_t
   use trc_eri_data, only: trc_bin_launch_t, trc_pair_quartets_t
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_sym_pairs, only: sym_pairs_symmetrize_dsh
   use trc_sym_scf, only: sym_grad_symmetrize
   use trc_test_basis, only: read_xyz
   use trc_decontract, only: decontract_expand_host
   implicit none

   character(len=256) :: xyz, bfile, aospace, mode, arg
   real(dp) :: thresh
   integer :: nat, i, j, g, dev, nw_sw_ok, nw_sw_bad, s
   integer :: nl_full, nl_skel, nseg_sw, nmaxl_sw
   integer, allocatable :: elem(:)
   real(dp), allocatable :: r(:, :), dens(:, :), d0(:, :)
   real(dp), allocatable :: dwork(:, :, :), gfull(:, :), gskel(:, :), graw(:, :), gref(:, :)
   real(dp), allocatable :: dshsym(:, :)
   type(trc_basis_t) :: bas
   type(trc_sym_geom_t) :: geom
   type(trc_sym_ao_t) :: ao
   type(trc_eri_t) :: eri
   type(trc_bin_plan_t) :: pfull, pskel, prep
   type(trc_bin_launch_t) :: launch
   type(error_t) :: err
   logical :: sph, verify, verify2
   type(trc_pair_quartets_t) :: qx
   type(trc_bin_plan_t) :: pfx, psx
   integer, allocatable :: pbase(:)
   integer :: nqpx
   real(dp), allocatable :: xz(:), xc(:), xa(:), xpx(:), xpy(:), xpz(:)
   integer(int64) :: ns_full, ns_skel, ncwork
   integer :: ncseg
   real(dp) :: emax, gmax, e2

   dev = trc_bind_device(0)
   call get_command_argument(1, xyz)
   call get_command_argument(2, bfile)
   call get_command_argument(3, aospace)
   call get_command_argument(4, mode)
   thresh = 1.0e-10_dp
   call get_command_argument(5, arg)
   if (len_trim(arg) > 0) read (arg, *) thresh
   sph = index(aospace, 'sph') > 0
   verify = index(mode, 'verify') > 0
   verify2 = index(mode, 'verify2') > 0

   call read_xyz(trim(xyz), nat, elem, r)
   geom = sym_detect(r, elem, err, basis_id=elem)
   if (err%has_error()) stop 'sym_detect failed'
   call trc_basis_from_json(trim(bfile), nat, elem, geom%r_std, bas, err, force_cart=.not. sph)
   if (err%has_error()) stop 'basis failed'
   bas%spherical = sph
   call bas%to_device()
   ao = sym_ao_build(bas, geom, err)
   if (err%has_error()) stop 'sym_ao_build failed'
   print '(a,a,a,i0,a,a,a,i0,a,i0,a,a)', 'molecule ', trim(xyz), '  natm ', nat, '  group ', trim(geom%group%name), &
      '  |G| ', geom%group%order, '  nao ', bas%nao, '  ', trim(aospace)

   call eri%build(bas, thresh, sym=geom)
   if (.not. eri%orb%built) stop 'no orbit tables (C1?)'

   ! An invariant density: a deterministic pseudo-random symmetric matrix,
   ! projected onto the group with the AO maps and signs. The gradient's
   ! in-kernel screen is density weighted, so the counts below are the
   ! counts for a density of this magnitude; the plan itself
   ! (plan%build(..., .false.)) does not see it.
   allocate (dens(bas%nao, bas%nao), d0(bas%nao, bas%nao))
   do j = 1, bas%nao
      do i = 1, bas%nao
         d0(i, j) = 0.5_dp*sin(1.7_dp*real(i*7 + j*13, dp))
      end do
   end do
   d0 = 0.5_dp*(d0 + transpose(d0))
   dens = d0
   do g = 2, ao%nop
      do j = 1, bas%nao
         do i = 1, bas%nao
            dens(i, j) = dens(i, j) + real(ao%ao_sign(i, g)*ao%ao_sign(j, g), dp)* &
                         d0(ao%ao_map(i, g), ao%ao_map(j, g))
         end do
      end do
   end do
   dens = dens/real(ao%nop, dp)
   ! invariance diagnostic
   emax = 0.0_dp
   do g = 2, ao%nop
      do j = 1, bas%nao
         do i = 1, bas%nao
            emax = max(emax, abs(dens(i, j) - real(ao%ao_sign(i, g)*ao%ao_sign(j, g), dp)* &
                                 dens(ao%ao_map(i, g), ao%ao_map(j, g))))
         end do
      end do
   end do
   print '(a,es11.3,a,es11.3)', 'density: max|D| ', maxval(abs(dens)), '  max|D - R D R^T| ', emax

   ! eri%gradient completes the pair data (the unique-pair build forms the
   ! representatives' only) and leaves eri%dsh holding this density's shell
   ! blocks, which is what the screen pass reads. Its own result is the
   ! reference.
   allocate (gfull(3, nat), gskel(3, nat), graw(3, nat), gref(3, nat))
   allocate (dwork(1, eri%nao_w, eri%nao_w))
   call eri%gradient(bas, reshape(dens, [bas%nao, bas%nao, 1]), gref)
   if (eri%dec%active) then
      call decontract_expand_host(eri%dec, dens, dwork(1, :, :))
   else
      dwork(1, :, :) = dens
   end if
   ! The skeleton build needs an invariant screen bound, so that the
   ! surviving quartet set is closed under the group (trc_eri's skeleton_dsh).
   ! Both runs below use the same symmetrized bound, so the comparison is of
   ! the schemes and not of their screens.
   graw = gref
   call sym_grad_symmetrize(geom%group, geom%perm, graw)
   print '(a,es11.3)', 'full gradient already group symmetric? max|g - sym(g)| ', maxval(abs(graw - gref))
   dshsym = eri%dsh
   call sym_pairs_symmetrize_dsh(eri%orb%shmap, dshsym)

   launch%thresh = thresh
   launch%jfac = 1.0_dp
   launch%kfac = 1.0_dp
   launch%omega = 0.0_dp
   launch%ndens = 1
   launch%nao = eri%nao_w
   launch%rank = 0
   launch%nranks = 1

   ! --- the two plans -------------------------------------------------
   call pfull%build(eri%bins, thresh, eri%quartets%bin_pbase, .false.)
   call pfull%sort_by_class(GCLASS_RADIX)
   call pskel%build_skeleton(eri%bins, eri%orb, thresh, eri%quartets%bin_pbase)
   call pskel%sort_by_class(GCLASS_RADIX)

   print '(a,i0,a,i0,a,f6.3,a)', 'plan full     segments ', pfull%nseg, '  work items ', pfull%nwork, &
      '  (1.000x)'
   print '(a,i0,a,i0,a,f6.3,a)', 'plan skeleton segments ', pskel%nseg, '  work items ', pskel%nwork, &
      '  reduction ', real(pfull%nwork, dp)/max(real(pskel%nwork, dp), 1.0_dp), 'x'

   ! --- the same plan with the bra/ket swap forbidden -------------------
   ! bins_skel_cut_first picks the side with fewer representatives x pairs,
   ! which can make (lc, ld) > (la, lb) -- a class the gradient's filter
   ! drops. Cutting always on the FULL PLAN'S OWN BRA (the first bin of the
   ! admitted pair) keeps the full plan's class set exactly, so it needs no
   ! new kernel and no new class run. This is what that costs in work items.
   call canon_counts(eri%bins, eri%orb, thresh, ncseg, ncwork)
   print '(a,i0,a,i0,a,f6.3,a)', 'plan skeleton, swap forbidden: segments ', ncseg, '  work items ', ncwork, &
      '  reduction ', real(pfull%nwork, dp)/max(real(ncwork, dp), 1.0_dp), 'x'

   ! --- the bra/ket swap against the gradient's class filter -----------
   nseg_sw = 0
   nw_sw_ok = 0
   nw_sw_bad = 0
   nmaxl_sw = 0
   do s = 1, pskel%nseg
      if (pskel%sLC(s)*16 + pskel%sLD(s) > pskel%sLA(s)*16 + pskel%sLB(s)) then
         nseg_sw = nseg_sw + 1
         nmaxl_sw = max(nmaxl_sw, max(pskel%sLA(s), pskel%sLB(s), pskel%sLC(s), pskel%sLD(s)))
         if (rysg_has_class(pskel%class_key(s, GCLASS_RADIX))) then
            nw_sw_ok = nw_sw_ok + 1
         else
            nw_sw_bad = nw_sw_bad + 1
         end if
      end if
   end do
   print '(a,i0,a,i0,a,i0,a,i0)', 'skeleton segments with (lc,ld) > (la,lb): ', nseg_sw, &
      '  of ', pskel%nseg, ';  with a gradient kernel ', nw_sw_ok, ', without ', nw_sw_bad
   print '(a,i0)', '  highest l in a swapped segment: ', nmaxl_sw

   ! --- weights present ------------------------------------------------
   call report_weights(pskel)

   ! --- survivors and launches -----------------------------------------
   call spike_grad_run(pfull, eri%bins, eri%pairs, eri%quartets, launch, dshsym, dwork, graw, .false., &
                       nl_full, ns_full)
   call spike_grad_run(pskel, eri%bins, eri%pairs, eri%quartets, launch, dshsym, dwork, graw, .false., &
                       nl_skel, ns_skel)
   print '(a,i0,a,i0,a,f6.3,a)', 'screen survivors full ', ns_full, '  skeleton ', ns_skel, &
      '  reduction ', real(ns_full, dp)/max(real(ns_skel, dp), 1.0_dp), 'x'
   print '(a,i0,a,i0)', 'screen launches (one per batch per class) full ', nl_full, '  skeleton ', nl_skel

   if (index(mode, 'plan') > 0) then
      call pfull%release()
      call pskel%release()
      call eri%release()
      call bas%release()
      stop 0
   end if

   if (verify) then
      ! The weighted skeleton gradient, by replication (see the header).
      call replicate_plan(pskel, prep)
      call prep%sort_by_class(GCLASS_RADIX)
      print '(a,i0,a,i0)', 'replicated plan segments ', prep%nseg, '  work items ', prep%nwork
      call spike_grad_run(prep, eri%bins, eri%pairs, eri%quartets, launch, dshsym, dwork, gskel, .true., &
                          nl_skel, ns_skel)
      gskel = 0.5_dp*gskel
      ! The same screen bound for the reference, so only the scheme differs.
      call spike_grad_run(pfull, eri%bins, eri%pairs, eri%quartets, launch, dshsym, dwork, gfull, .true., &
                          nl_full, ns_full)
      print '(a,es11.3)', 'harness check: my full plan against eri%gradient ', maxval(abs(gfull - gref))
      graw = gskel
      call sym_grad_symmetrize(geom%group, geom%perm, gskel)
      gmax = maxval(abs(gfull))
      emax = maxval(abs(gskel - gfull))
      e2 = maxval(abs(graw - gfull))
      print '(a,es11.3,a,es11.3,a,es11.3)', 'max|grad_full| ', gmax, &
         '   raw skeleton err ', e2, '   symmetrized err ', emax
      print '(a,es11.3)', 'relative to max|grad_full|: ', emax/max(gmax, 1.0_dp)
      do i = 1, nat
         print '(a,i0,a,3f14.8,a,3f14.8,a,3f14.8)', '  atom ', i, ' full ', gfull(:, i), ' | skel raw ', graw(:, i), &
            ' | skel sym ', gskel(:, i)
      end do
      print '(a,es11.3,a,es11.3)', 'translational invariance: full ', maxval(abs(sum(gfull, dim=2))), &
         '  skeleton ', maxval(abs(sum(graw, dim=2)))
      call prep%release()
   end if

   if (verify2) then
      ! EXACT ORACLE. One segment per quartet over a stride-one primitive
      ! layout, so the gradient kernels (which take no sSA/sSB) are correct
      ! on a skeleton segment too; the orbit weight applied by repetition.
      call relayout_pairs(eri%bins, eri%quartets, pbase, nqpx, xz, xc, xa, xpx, xpy, xpz)
      qx = eri%quartets
      qx%nqp = nqpx
      qx%z = xz; qx%c = xc; qx%a = xa
      qx%px = xpx; qx%py = xpy; qx%pz = xpz
      call expand_plan(pfull, pbase, pfx, .false.)
      call pfx%sort_by_class(GCLASS_RADIX)
      call expand_plan(pskel, pbase, psx, .true.)
      call psx%sort_by_class(GCLASS_RADIX)
      print '(a,i0,a,i0)', 'expanded plans: full segments ', pfx%nseg, '  skeleton (2x weighted) segments ', psx%nseg
      launch%nao = eri%nao_w
      call spike_grad_run(pfx, eri%bins, eri%pairs, qx, launch, dshsym, dwork, gfull, .true., nl_full, ns_full)
      print '(a,es11.3)', 'expanded full against eri%gradient: ', maxval(abs(gfull - gref))
      call spike_grad_run(psx, eri%bins, eri%pairs, qx, launch, dshsym, dwork, gskel, .true., nl_skel, ns_skel)
      gskel = 0.5_dp*gskel
      graw = gskel
      call sym_grad_symmetrize(geom%group, geom%perm, gskel)
      gmax = maxval(abs(gref))
      print '(a,es11.3,a,es11.3,a,es11.3)', 'max|grad_full| ', gmax, &
         '   raw skeleton err ', maxval(abs(graw - gref)), '   symmetrized err ', maxval(abs(gskel - gref))
      print '(a,es11.3)', 'symmetrized err relative to max|grad_full|: ', maxval(abs(gskel - gref))/max(gmax, 1.0_dp)
      print '(a,es11.3,a,es11.3)', 'translational invariance: full ', maxval(abs(sum(gref, dim=2))), &
         '  raw skeleton ', maxval(abs(sum(graw, dim=2)))
      call pfx%release()
      call psx%release()
   end if

   call pfull%release()
   call pskel%release()
   call eri%release()
   call bas%release()

contains

   subroutine canon_counts(b, orb, thresh, nseg, nwork)
      !! Segments and work items of a skeleton plan that always cuts on the
      !! first bin of an admitted bin pair (no bra/ket swap).
      use trc_bins, only: pair_bins_t
      use trc_eri_data, only: trc_bin_orbits_t
      type(pair_bins_t), intent(in) :: b
      type(trc_bin_orbits_t), intent(in) :: orb
      real(dp), intent(in) :: thresh
      integer, intent(out) :: nseg
      integer(int64), intent(out) :: nwork
      integer :: ia, ib, ka, kb, smax_keep, rr
      smax_keep = int(-log10(thresh))*b%sres
      nseg = 0
      nwork = 0_int64
      do ia = 1, b%nlive
         ka = b%live(ia)
         do ib = 1, ia
            kb = b%live(ib)
            if (b%bin_s(ka) + b%bin_s(kb) > smax_keep) cycle
            if (b%bin_cnt(ka) == 0 .or. b%bin_cnt(kb) == 0) cycle
            do rr = orb%sec0(ka) + 1, orb%sec0(ka + 1)
               nseg = nseg + 1
               nwork = nwork + int(orb%sec_cnt(rr), int64)*int(b%bin_cnt(kb), int64)
            end do
         end do
      end do
   end subroutine canon_counts

   subroutine report_weights(plan)
      type(trc_bin_plan_t), intent(in) :: plan
      integer :: k, nk
      real(dp) :: wv(32)
      integer(int64) :: wn(32)
      integer :: t
      character(len=512) :: line
      nk = 0
      wn = 0_int64
      do t = 1, plan%nseg
         do k = 1, nk
            if (abs(wv(k) - plan%sW(t)) < 1.0e-12_dp) exit
         end do
         if (k > nk) then
            nk = min(nk + 1, 32)
            wv(nk) = plan%sW(t)
            k = nk
         end if
         wn(k) = wn(k) + (plan%sOff(t + 1) - plan%sOff(t))
      end do
      line = 'skeleton weights (work items each):'
      do k = 1, nk
         write (line, '(a,1x,f5.2,a,i0,a)') trim(line), wv(k), ':', wn(k), ''
      end do
      print '(a)', trim(line)
   end subroutine report_weights

end program spike_symgrad
