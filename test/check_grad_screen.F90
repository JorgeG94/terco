module check_grad_screen_suite
   !!
   !! The density screen of the two-electron gradient (trc_grad's SCREENING
   !! paragraph): a quartet is dropped when
   !!     Q * max(blocks of the terms the launch has) <= thresh,
   !! Q the energy Schwarz product and the blocks 4 D_ab, 4 D_cd (Coulomb,
   !! when jfac /= 0) and D_ac, D_ad, D_bc, D_bd (exchange, when kfac /= 0),
   !! D_xy the shell-block maximum `eri%dsh`. A pair table has a term when
   !! some pair of the sweep weights it.
   !!
   !! What each test pins, and the change that makes it fail:
   !!   - FUSED launches (both terms) are what the six-block screen gave:
   !!     `eri%grad_screen_all_blocks` forces the old test, so the comparison
   !!     cannot pass vacuously (the same run is also checked to be screening
   !!     at all, fewer survivors than with no density screen). The survivor
   !!     counts are equal exactly; the gradient is equal bit for bit on the
   !!     serial host build, and to round-off where atomics order the sums
   !!     (`same_gradient`). Counting too few blocks in a fused launch changes
   !!     the count.
   !!   - EXCHANGE-ONLY (jfac = 0, what the MD Coulomb gradient leaves to the
   !!     Rys kernels and the long-range pass of a range-separated hybrid
   !!     is) and COULOMB-ONLY (kfac = 0) launches keep strictly fewer
   !!     quartets than the six-block screen, and their gradient is within
   !!     the dropped quartets' own bound of the one with no density screen
   !!     (`eri%grad_density_screen`). Counting the wrong terms' blocks, or
   !!     a term as absent while it is present, drops quartets that matter
   !!     and breaks that bound by orders of magnitude. The kernel's survivor
   !!     count is also compared, EXACTLY, with a host loop over the work
   !!     items of the very plan the launch walks, applying the same
   !!     inequality: a screen that drops a few percent too many quartets
   !!     passes the (loose) accuracy bound but not the count.
   !!   - The PAIR TABLE: a sweep whose pairs all have pj = 0 (or all pk = 0)
   !!     as the above; a mixed sweep, one nonzero weight among zeros
   !!     included, is the six-block screen (counts exact, gradient as above).
   !!   - The MD Coulomb switch (`grad_j_engine`): full HF, global-hybrid and
   !!     range-separated-hybrid gradients equal the fused ones to check_mdj_grad's
   !!     tolerance, RHF and UHF, with the exchange-only Rys pass keeping
   !!     fewer quartets than the fused pass.
   !!
   !! The bound of the accuracy checks. A quartet dropped at Q*dmax <= thresh
   !! contributes 8 wq Gamma d(ab|cd)/dR with |Gamma| <= dmax (the blocks of
   !! the terms present) and, as the screen assumes throughout, a derivative
   !! integral of the size of the energy bound Q, so at most 8 wq Q dmax <=
   !! 8 wq thresh. The host loop sums 8 wq Q dmax over the dropped quartets
   !! and that sum, plus the round-off of a different summation order, is the
   !! tolerance: no number of its own.
   !!
   !! The densities are deterministic synthetic matrices that decay with the
   !! distance between the atoms of their two functions, so the screen bites
   !! on an extended molecule (glycine tripeptide); this checks the screen,
   !! not an SCF.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t
   use trc_tables, only: LMAX
   use trc_cart, only: cart_ncart_of
   use trc_bins, only: trc_bin_plan_t
   use trc_eri, only: trc_eri_t, trc_grad_pair_t, TRC_J_ENGINE_MDJ, TRC_J_ENGINE_BINS
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_grad_screen_tests, release_contexts

   character(len=*), parameter :: XYZ_SCREEN = "gly3.xyz", XYZ_TIGHT = "benzene.xyz"
   real(dp), parameter :: THRESH_SCREEN = 1.0e-4_dp
   !! Where the screen bites on the tripeptide.
   real(dp), parameter :: THRESH_TIGHT = 1.0e-14_dp
   !! check_mdj_grad's threshold, for the MD-against-fused comparison.
   real(dp), parameter :: ROUNDOFF = 1.0e-11_dp
   !! Two runs of the same derivative integrals over different survivor
   !! lists, atomics summed in a different order, relative to max(1, max|g|).
   real(dp), parameter :: AGREE = 3.0e-13_dp
   !! check_mdj_grad's MD-against-Rys tolerance, relative to max(1, max|g|).

   type :: ctx_t
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      integer :: nat = 0
      real(dp), allocatable :: r(:, :)
      type(trc_bin_plan_t) :: plan
      !! The plan the gradient launches walk, rebuilt on the host for the exact count.
      logical :: plan_built = .false.
      logical :: built = .false.
   end type ctx_t

   type(ctx_t), save, target :: cache_screen, cache_tight

   type :: count_t
      !! What one host loop over the bin pairs finds for one set of terms.
      integer(int64) :: n_items = 0, n_schwarz = 0, n_kept = 0
      !! Work items of the plan; those past the Schwarz and same-atom
      !! screens; those past the density screen too.
      real(dp) :: dropped = 0.0_dp
      !! Sum over the quartets Schwarz keeps and the screen drops of
      !! wq*Q*dmax, dmax the blocks of the terms present.
   end type count_t

contains

   subroutine collect_grad_screen_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fused launches (J and K) keep the six-block screen", test_fused), &
                  new_unittest("exchange-only launch screens on the exchange blocks", test_exchange_only), &
                  new_unittest("Coulomb-only launch screens on the Coulomb blocks", test_coulomb_only), &
                  new_unittest("pair table: J-free, K-free and mixed sweeps", test_pairs), &
                  new_unittest("MD Coulomb switch: full gradients equal the fused ones", test_mdj_switch) &
                  ]
   end subroutine collect_grad_screen_tests

   subroutine release_contexts()
      if (cache_screen%plan_built) call cache_screen%plan%release()
      if (cache_screen%built) then
         call cache_screen%eri%release()
         call cache_screen%bas%release()
      end if
      if (cache_tight%built) then
         call cache_tight%eri%release()
         call cache_tight%bas%release()
      end if
   end subroutine release_contexts

   subroutine get_ctx(ctx, xyz, thresh, error)
      !! The tripeptide in the command-line basis, built once per threshold
      !! and reused by every test.
      type(ctx_t), intent(inout) :: ctx
      character(*), intent(in) :: xyz
      real(dp), intent(in) :: thresh
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      integer, allocatable :: at_z(:)
      type(error_t) :: err
      if (ctx%built) return
      basis_json = ''
      if (command_argument_count() >= 1) call get_command_argument(1, basis_json)
      call check(error, len_trim(basis_json) > 0, "usage: check_grad_screen <basis.json>")
      if (allocated(error)) return
      call read_xyz(xyz, ctx%nat, at_z, ctx%r)
      call trc_basis_from_json(trim(basis_json), ctx%nat, at_z, ctx%r, ctx%bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "check_grad_screen: "//err%get_message())
      if (allocated(error)) return
      if (maxval(ctx%bas%sh_l) > LMAX) then
         call skip_test(error, "SKIP: basis needs l = "//to_string(maxval(ctx%bas%sh_l)))
         call ctx%bas%release()
         return
      end if
      call ctx%bas%to_device()
      call ctx%eri%build(ctx%bas, thresh)
      call check(error, .not. ctx%eri%dec%active, "check_grad_screen: expects a segmented basis")
      ctx%built = .true.
   end subroutine get_ctx

   subroutine synth(ctx, centre, width, rate, phase, d)
      !! A symmetric density that decays with the distance between the two
      !! functions' atoms (exp(-rate*dist)) and is concentrated around the
      !! fraction `centre` of the molecule's x extent (Gaussian of relative
      !! `width`), the same shape on both sides of any shell pair, so the
      !! shell-block maxima follow the decay: J blocks large where K blocks
      !! are small and the other way round.
      type(ctx_t), intent(in) :: ctx
      real(dp), intent(in) :: centre, width, rate
      integer, intent(in) :: phase
      real(dp), intent(out) :: d(:, :)
      integer :: i, j, ai, aj, si
      real(dp), allocatable :: wa(:)
      real(dp) :: x0, xlo, xhi, dist
      integer, allocatable :: atom_of(:)
      allocate (atom_of(size(d, 1)), wa(ctx%nat))
      do si = 1, ctx%eri%nbas
         atom_of(ctx%eri%pairs%ao_off(si):ctx%eri%pairs%ao_off(si) + cart_ncart_of(ctx%eri%pairs%sh_l(si)) - 1) = &
            ctx%eri%pairs%sh_at(si)
      end do
      xlo = minval(ctx%r(1, :))
      xhi = maxval(ctx%r(1, :))
      x0 = xlo + centre*(xhi - xlo)
      do ai = 1, ctx%nat
         wa(ai) = exp(-((ctx%r(1, ai) - x0)/(width*(xhi - xlo)))**2)
      end do
      do j = 1, size(d, 2)
         aj = atom_of(j)
         do i = 1, size(d, 1)
            ai = atom_of(i)
            dist = norm2(ctx%r(:, ai) - ctx%r(:, aj))
            d(i, j) = wa(ai)*wa(aj)*exp(-rate*dist)*(0.7_dp + 0.2_dp*cos(real(3*i + 5*j + 7*phase, dp)))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth

   subroutine host_count(ctx, dsh, cj, ck, cnt)
      !! The survivors of the Schwarz, same-atom and density screens, by a host
      !! loop over EXACTLY the work items the gradient's plan admits (the plan
      !! is built here as grad_bins builds it, so its bin-level rounding is
      !! the kernel's own), with the library's own dsh from the run being
      !! judged. The kernel's count must equal n_kept with no slack.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dsh(:, :)
      logical, intent(in) :: cj, ck
      type(count_t), intent(out) :: cnt
      integer :: seg, iab, icd, ilo, si, sj, sk, sl, nab, ncd
      real(dp) :: qcut, dmax, wq, th
      th = ctx%eri%thresh
      if (.not. ctx%plan_built) then
         call ctx%plan%build(ctx%eri%bins, th, ctx%eri%quartets%bin_pbase, .false., canonical=.true.)
         ctx%plan_built = .true.
      end if
      do seg = 1, ctx%plan%nseg
         nab = ctx%plan%sA(seg)
         ncd = ctx%plan%sNB(seg)
         do icd = 1, ncd
            ilo = 1
            if (ctx%plan%sD(seg)) ilo = icd
            do iab = ilo, nab
               cnt%n_items = cnt%n_items + 1
               qcut = ctx%eri%bins%sp_q(ctx%plan%sOA(seg) + iab)*ctx%eri%bins%sp_q(ctx%plan%sOB(seg) + icd)
               if (qcut <= th) cycle
               si = ctx%eri%bins%sp_i(ctx%plan%sOA(seg) + iab)
               sj = ctx%eri%bins%sp_j(ctx%plan%sOA(seg) + iab)
               sk = ctx%eri%bins%sp_i(ctx%plan%sOB(seg) + icd)
               sl = ctx%eri%bins%sp_j(ctx%plan%sOB(seg) + icd)
               dmax = 0.0_dp
               if (cj) dmax = max(4.0_dp*dsh(si, sj), 4.0_dp*dsh(sk, sl))
               if (ck) dmax = max(dmax, dsh(si, sk), dsh(si, sl), dsh(sj, sk), dsh(sj, sl))
               if (qcut*dmax <= th) then
                  ! Dropped by the density screen (a same-atom quartet is
                  ! dropped after it, so it adds to neither count).
                  if (.not. same_atom(ctx, si, sj, sk, sl)) then
                     cnt%n_schwarz = cnt%n_schwarz + 1
                     wq = 1.0_dp
                     if (si == sj) wq = 0.5_dp*wq
                     if (sk == sl) wq = 0.5_dp*wq
                     if (ctx%plan%sD(seg) .and. iab == icd) wq = 0.5_dp*wq
                     cnt%dropped = cnt%dropped + 8.0_dp*wq*qcut*dmax
                  end if
                  cycle
               end if
               if (same_atom(ctx, si, sj, sk, sl)) cycle
               cnt%n_schwarz = cnt%n_schwarz + 1
               cnt%n_kept = cnt%n_kept + 1
            end do
         end do
      end do
   end subroutine host_count

   logical function same_atom(ctx, si, sj, sk, sl)
      type(ctx_t), intent(in) :: ctx
      integer, intent(in) :: si, sj, sk, sl
      same_atom = ctx%eri%pairs%sh_at(si) == ctx%eri%pairs%sh_at(sj) .and. &
                  ctx%eri%pairs%sh_at(si) == ctx%eri%pairs%sh_at(sk) .and. &
                  ctx%eri%pairs%sh_at(si) == ctx%eri%pairs%sh_at(sl)
   end function same_atom

   subroutine run_gradient(ctx, dens, jfac, kfac, omega, mode, g, nkept)
      !! eri%gradient under one of three screens: 0 the new one, 1 the
      !! six-block one (the old behaviour), 2 no density screen (Schwarz
      !! only). Returns the gradient, the survivor count and, as `dsh`
      !! after the call, the block maxima of mode 0 and 1 (kept apart in
      !! ctx%eri%dsh by the caller before another mode overwrites it).
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dens(:, :, :), jfac, kfac, omega
      integer, intent(in) :: mode
      real(dp), intent(out) :: g(:, :)
      integer(int64), intent(out) :: nkept
      ! The fused Rys gradient named explicitly, so a run with TRC_MDJ_GRAD=1
      ! in the environment still tests the screen of the launches it names.
      ctx%eri%grad_j_engine = TRC_J_ENGINE_BINS
      ctx%eri%grad_screen_all_blocks = mode == 1
      ctx%eri%grad_density_screen = mode /= 2
      call ctx%eri%gradient(ctx%bas, dens, g, k_scale=kfac, j_scale=jfac, omega=omega)
      ctx%eri%grad_screen_all_blocks = .false.
      ctx%eri%grad_density_screen = .true.
      nkept = ctx%eri%nkept
   end subroutine run_gradient

   logical function same_gradient(a, b)
      !! Two gradients over the same survivors. Exactly equal on the serial
      !! host build; where the per-atom sums are atomics (the device build,
      !! the OpenMP port) their thread order differs between two runs, so
      !! relative to max(1, max|g|) at round-off, the tolerance
      !! check_fock_split uses for two orders of the same sum.
      real(dp), intent(in) :: a(:, :), b(:, :)
#if defined(_OPENACC) || defined(_OPENMP)
      same_gradient = maxval(abs(a - b)) <= 1.0e-13_dp*max(1.0_dp, maxval(abs(b)))
#else
      same_gradient = all(a == b)
#endif
   end function same_gradient

   subroutine densities(ctx, ndens, dens, rate)
      !! RHF: one density. UHF: alpha concentrated left, beta right, each
      !! dominant where the other is small.
      type(ctx_t), intent(in) :: ctx
      integer, intent(in) :: ndens
      real(dp), allocatable, intent(out) :: dens(:, :, :)
      real(dp), intent(in) :: rate
      !! Decay rate of the density with the distance between atoms.
      allocate (dens(ctx%eri%nao_w, ctx%eri%nao_w, ndens))
      if (ndens == 1) then
         call synth(ctx, 0.5_dp, 0.6_dp, rate, 1, dens(:, :, 1))
      else
         call synth(ctx, 0.25_dp, 0.5_dp, rate, 1, dens(:, :, 1))
         call synth(ctx, 0.75_dp, 0.5_dp, rate, 2, dens(:, :, 2))
      end if
   end subroutine densities

   subroutine exact_counts(label, ctx, dsh, cj, ck, nn, no, np, cnt, error)
      !! The kernel's three survivor counts against the host loop over the
      !! plan's own work items, EXACTLY: the new screen (terms cj, ck), the
      !! six-block one and Schwarz only. cnt is the host result for (cj, ck).
      character(*), intent(in) :: label
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dsh(:, :)
      logical, intent(in) :: cj, ck
      integer(int64), intent(in) :: nn, no, np
      type(count_t), intent(out) :: cnt
      type(error_type), allocatable, intent(out) :: error
      type(count_t) :: cnt6
      call host_count(ctx, dsh, cj, ck, cnt)
      call host_count(ctx, dsh, .true., .true., cnt6)
      call check(error, cnt%n_items == ctx%eri%nwork, label//": the host plan has "//to_string(cnt%n_items)// &
                 " work items, the launch "//to_string(ctx%eri%nwork))
      if (allocated(error)) return
      call check(error, nn == cnt%n_kept, label//": the kernel kept "//to_string(nn)//" quartets, the host loop over the "// &
                 "plan's work items with the same inequality "//to_string(cnt%n_kept))
      if (allocated(error)) return
      call check(error, no == cnt6%n_kept, label//": the six-block run kept "//to_string(no)//", the host loop "// &
                 to_string(cnt6%n_kept))
      if (allocated(error)) return
      call check(error, np == cnt%n_schwarz, label//": the Schwarz-only run kept "//to_string(np)//", the host loop "// &
                 to_string(cnt%n_schwarz))
   end subroutine exact_counts

   subroutine test_fused(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :, :), gn(:, :), go(:, :), gp(:, :), dsh(:, :)
      integer(int64) :: nn, no, np
      type(count_t) :: cnt
      integer :: ndens, ic, ncase
      ! (jfac, kfac): Hartree-Fock, then a Coulomb term present only by a
      ! hair, which must still count its blocks (RHF; UHF takes Hartree-Fock).
      real(dp), parameter :: JF(2) = [1.0_dp, 1.0e-9_dp]
      real(dp), parameter :: KF(2) = [1.0_dp, 0.25_dp]

      call get_ctx(cache_screen, XYZ_SCREEN, THRESH_SCREEN, error)
      if (allocated(error)) return
      allocate (gn(3, cache_screen%nat), go(3, cache_screen%nat), gp(3, cache_screen%nat))
      do ndens = 1, 2
         call densities(cache_screen, ndens, d, 1.2_dp)
         ncase = merge(2, 1, ndens == 1)
         do ic = 1, ncase
            call run_gradient(cache_screen, d, JF(ic), KF(ic), 0.0_dp, 0, gn, nn)
            dsh = cache_screen%eri%dsh
            call run_gradient(cache_screen, d, JF(ic), KF(ic), 0.0_dp, 1, go, no)
            call run_gradient(cache_screen, d, JF(ic), KF(ic), 0.0_dp, 2, gp, np)
            call exact_counts("fused ndens "//to_string(ndens)//" case "//to_string(ic), cache_screen, dsh, &
                              .true., .true., nn, no, np, cnt, error)
            if (allocated(error)) return
            print '(a,i0,a,i0,a,i0,a,i0,a,i0,a,l1)', '  fused ndens ', ndens, ' case ', ic, ': survivors new ', nn, &
               ', six-block ', no, ', Schwarz only ', np, '  same gradient ', same_gradient(gn, go)
            call check(error, same_gradient(gn, go), "fused launch, ndens "//to_string(ndens)//" case "//to_string(ic)// &
                       ": the gradient differs from the six-block screen's by "//to_string(maxval(abs(gn - go)), "es10.2"))
            if (allocated(error)) return
            call check(error, nn == no, "fused launch: "//to_string(nn)//" survivors against the six-block "//to_string(no))
            if (allocated(error)) return
            ! The six-block run really is screening (else the equality would be vacuous).
            call check(error, 10*no < 9*np, "the density screen removes under 10% of the Schwarz survivors ("// &
                       to_string(no)//" of "//to_string(np)//"): the system does not exercise it")
            if (allocated(error)) return
         end do
         deallocate (d)
      end do
   end subroutine test_fused

   subroutine one_term_launch(label, ndens, jfac, kfac, cj, ck, error)
      !! A launch with one term: survivors against the six-block screen and,
      !! exactly, against the host loop, gradient against the no-density-screen run
      !! within the dropped quartets' bound.
      character(*), intent(in) :: label
      integer, intent(in) :: ndens
      real(dp), intent(in) :: jfac, kfac
      logical, intent(in) :: cj, ck
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :, :), gn(:, :), go(:, :), gp(:, :)
      real(dp), allocatable :: dsh(:, :)
      integer(int64) :: nn, no, np
      type(count_t) :: cnt
      real(dp) :: diff, tol, scale

      call densities(cache_screen, ndens, d, 1.2_dp)
      allocate (gn(3, cache_screen%nat), go(3, cache_screen%nat), gp(3, cache_screen%nat))
      call run_gradient(cache_screen, d, jfac, kfac, 0.0_dp, 0, gn, nn)
      dsh = cache_screen%eri%dsh
      call run_gradient(cache_screen, d, jfac, kfac, 0.0_dp, 1, go, no)
      call run_gradient(cache_screen, d, jfac, kfac, 0.0_dp, 2, gp, np)
      call exact_counts(label, cache_screen, dsh, cj, ck, nn, no, np, cnt, error)
      if (allocated(error)) return
      diff = maxval(abs(gn - gp))
      scale = max(1.0_dp, maxval(abs(gp)))
      tol = cnt%dropped + ROUNDOFF*scale
      print '(a,a,a,i0,a,i0,a,i0,a,i0)', '  ', label, ': survivors new ', nn, ', six-block ', no, ', Schwarz only ', np, &
         ', host loop ', cnt%n_kept
      print '(a,es10.2,a,es10.2,a,es10.2)', '    max |g_screen - g_noscreen| ', diff, ' against bound ', tol, &
         '  (max |g| ', scale
      call check(error, nn < no, label//": the screen keeps "//to_string(nn)//" quartets, the six-block screen "// &
                 to_string(no)//": no tighter")
      if (allocated(error)) return
      call check(error, diff <= tol, label//": the screened gradient is "//to_string(diff, "es10.2")// &
                 " off the unscreened one, bound "//to_string(tol, "es10.2"))
   end subroutine one_term_launch

   subroutine test_exchange_only(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :, :), g(:, :)
      integer(int64) :: nk

      call get_ctx(cache_screen, XYZ_SCREEN, THRESH_SCREEN, error)
      if (allocated(error)) return
      call one_term_launch("exchange only RHF", 1, 0.0_dp, 0.8_dp, .false., .true., error)
      if (allocated(error)) return
      call one_term_launch("exchange only UHF", 2, 0.0_dp, 0.8_dp, .false., .true., error)
      if (allocated(error)) return
      ! Neither term: nothing to contract, nothing launched.
      call densities(cache_screen, 1, d, 1.2_dp)
      allocate (g(3, cache_screen%nat))
      call run_gradient(cache_screen, d, 0.0_dp, 0.0_dp, 0.0_dp, 0, g, nk)
      call check(error, nk == 0 .and. cache_screen%eri%nwork == 0 .and. all(g == 0.0_dp), &
                 "jfac = kfac = 0 launched "//to_string(nk)//" quartets, nwork "//to_string(cache_screen%eri%nwork))
   end subroutine test_exchange_only

   subroutine test_coulomb_only(error)
      type(error_type), allocatable, intent(out) :: error
      call get_ctx(cache_screen, XYZ_SCREEN, THRESH_SCREEN, error)
      if (allocated(error)) return
      call one_term_launch("Coulomb only UHF", 2, 1.0_dp, 0.0_dp, .true., .false., error)
   end subroutine test_coulomb_only

   subroutine pair_sweep(ctx, pool, ps, mode, g, nkept, dsh)
      !! eri%gradient_pairs under a screen mode (as run_gradient's).
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: pool(:, :, :)
      type(trc_grad_pair_t), intent(in) :: ps(:)
      integer, intent(in) :: mode
      real(dp), intent(out) :: g(:, :)
      integer(int64), intent(out) :: nkept
      real(dp), allocatable, intent(out), optional :: dsh(:, :)
      real(dp) :: gg(3, ctx%nat, 1)
      type(error_t) :: err
      ctx%eri%grad_screen_all_blocks = mode == 1
      ctx%eri%grad_density_screen = mode /= 2
      call ctx%eri%gradient_pairs(ctx%bas, pool, ps, 1, gg, err)
      if (err%has_error()) error stop "check_grad_screen: gradient_pairs: "//err%get_message()
      if (present(dsh)) dsh = ctx%eri%dsh
      ctx%eri%grad_screen_all_blocks = .false.
      ctx%eri%grad_density_screen = .true.
      nkept = ctx%eri%nkept
      g = gg(:, :, 1)
   end subroutine pair_sweep

   subroutine pair_launch(label, pool, ps, cj, ck, expect_tighter, error)
      !! One sweep: against the six-block screen, the Schwarz-only run and
      !! the host loop, as one_term_launch does for eri%gradient.
      character(*), intent(in) :: label
      real(dp), intent(in) :: pool(:, :, :)
      type(trc_grad_pair_t), intent(in) :: ps(:)
      logical, intent(in) :: cj, ck, expect_tighter
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gn(:, :), go(:, :), gp(:, :), dsh(:, :)
      integer(int64) :: nn, no, np
      type(count_t) :: cnt
      real(dp) :: diff, tol, scale

      allocate (gn(3, cache_screen%nat), go(3, cache_screen%nat), gp(3, cache_screen%nat))
      call pair_sweep(cache_screen, pool, ps, 0, gn, nn, dsh)
      call pair_sweep(cache_screen, pool, ps, 1, go, no)
      call pair_sweep(cache_screen, pool, ps, 2, gp, np)
      call exact_counts("pairs "//label, cache_screen, dsh, cj, ck, nn, no, np, cnt, error)
      if (allocated(error)) return
      diff = maxval(abs(gn - gp))
      scale = max(1.0_dp, maxval(abs(gp)))
      tol = cnt%dropped + ROUNDOFF*scale
      print '(a,a,a,i0,a,i0,a,i0,a,i0,a,l1)', '  pairs ', label, ': survivors new ', nn, ', six-block ', no, &
         ', Schwarz only ', np, ', host loop ', cnt%n_kept, '  same gradient as six-block ', same_gradient(gn, go)
      print '(a,es10.2,a,es10.2)', '    max |g_screen - g_noscreen| ', diff, ' against bound ', tol
      if (expect_tighter) then
         call check(error, nn < no, "pairs "//label//": "//to_string(nn)//" survivors, the six-block screen "//to_string(no))
      else
         call check(error, nn == no .and. same_gradient(gn, go), "pairs "//label//": a mixed sweep differs from the six-block "// &
                    "screen ("//to_string(nn)//" against "//to_string(no)//" survivors, gradient by "// &
                    to_string(maxval(abs(gn - go)), "es10.2")//")")
      end if
      if (allocated(error)) return
      call check(error, 10*no < 9*np, "pairs "//label//": the density screen removes under 10% of the Schwarz survivors")
      if (allocated(error)) return
      call check(error, diff <= tol, "pairs "//label//": the screened gradient is "//to_string(diff, "es10.2")// &
                 " off the unscreened one, bound "//to_string(tol, "es10.2"))
   end subroutine pair_launch

   subroutine test_pairs(error)
      !! A pool of two densities (alpha-like and beta-like), one column.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: pool(:, :, :)
      type(trc_grad_pair_t) :: ps(3)

      call get_ctx(cache_screen, XYZ_SCREEN, THRESH_SCREEN, error)
      if (allocated(error)) return
      call densities(cache_screen, 2, pool, 1.2_dp)
      ! Exchange only: every Coulomb weight zero.
      ps(1) = trc_grad_pair_t(left=1, right=1, column=1, coulomb=0.0_dp, exchange=0.8_dp)
      ps(2) = trc_grad_pair_t(left=2, right=2, column=1, coulomb=0.0_dp, exchange=0.8_dp)
      ps(3) = trc_grad_pair_t(left=1, right=2, column=1, coulomb=0.0_dp, exchange=0.4_dp)
      call pair_launch("exchange only", pool, ps, .false., .true., .true., error)
      if (allocated(error)) return
      ! Coulomb only: every exchange weight zero.
      ps(1)%coulomb = 0.5_dp; ps(1)%exchange = 0.0_dp
      ps(2)%coulomb = 0.5_dp; ps(2)%exchange = 0.0_dp
      ps(3)%coulomb = 0.25_dp; ps(3)%exchange = 0.0_dp
      call pair_launch("Coulomb only", pool, ps, .true., .false., .true., error)
      if (allocated(error)) return
      ! One nonzero Coulomb weight among zeros: the sweep has a Coulomb term,
      ! so it screens as the six-block test.
      ps(1)%coulomb = 0.0_dp; ps(1)%exchange = 0.8_dp
      ps(2)%coulomb = 0.0_dp; ps(2)%exchange = 0.8_dp
      ps(3)%coulomb = 1.0e-9_dp; ps(3)%exchange = 0.4_dp
      call pair_launch("one nonzero Coulomb weight", pool, ps, .true., .true., .false., error)
      if (allocated(error)) return
      ! The same for one nonzero exchange weight among zeros.
      ps(1)%coulomb = 0.5_dp; ps(1)%exchange = 0.0_dp
      ps(2)%coulomb = 0.5_dp; ps(2)%exchange = 0.0_dp
      ps(3)%coulomb = 0.25_dp; ps(3)%exchange = 1.0e-9_dp
      call pair_launch("one nonzero exchange weight", pool, ps, .true., .true., .false., error)
   end subroutine test_pairs

   subroutine test_mdj_switch(error)
      !! The MD Coulomb gradient plus the exchange-only Rys pass against the
      !! fused kernels (check_mdj_grad's cases and tolerance), RHF (all three
      !! cases) and UHF (the range-separated one, which has a Coulomb-plus-exchange
      !! call and an exchange-only one),
      !! at check_mdj_grad's threshold, and the exchange-only pass keeps
      !! fewer quartets than the fused one.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :, :), gmd(:, :), gfu(:, :), gm2(:, :), gf2(:, :)
      integer(int64) :: nmd, nfu, nmd2, nfu2
      integer :: ndens, ic
      real(dp) :: dmax
      character(len=24), parameter :: NAMES(3) = [character(len=24) :: "HF", "global hybrid", "RSH (SR + LR calls)"]

      call get_ctx(cache_tight, XYZ_TIGHT, THRESH_TIGHT, error)
      if (allocated(error)) return
      allocate (gmd(3, cache_tight%nat), gfu(3, cache_tight%nat), gm2(3, cache_tight%nat), gf2(3, cache_tight%nat))
      do ndens = 1, 2
         call densities(cache_tight, ndens, d, 3.0_dp)
         do ic = 1, 3
            if (ndens == 2 .and. ic /= 3) cycle
            call md_and_fused(cache_tight, d, merge(1.0_dp, merge(0.25_dp, 0.2_dp, ic == 2), ic == 1), 1.0_dp, 0.0_dp, &
                              gmd, gfu, nmd, nfu)
            if (ic == 3) then
               ! The long-range exchange alone (jfac = 0, omega > 0), as the SCF makes it.
               call md_and_fused(cache_tight, d, 0.6_dp, 0.0_dp, 0.33_dp, gm2, gf2, nmd2, nfu2)
               gmd = gmd + gm2
               gfu = gfu + gf2
            end if
            dmax = maxval(abs(gmd - gfu))
            print '(a,i0,a,a,a,es10.2,a,i0,a,i0)', '  MD switch ndens ', ndens, ' ', trim(NAMES(ic)), &
               ': max |MD - fused| ', dmax, '  Rys survivors MD ', nmd, ', fused ', nfu
            call check(error, dmax <= AGREE*max(1.0_dp, maxval(abs(gfu))), "ndens "//to_string(ndens)//" "//trim(NAMES(ic))// &
                       ": max |MD - fused| "//to_string(dmax, "es10.2"))
            if (allocated(error)) return
            call check(error, nmd < nfu, "ndens "//to_string(ndens)//" "//trim(NAMES(ic))//": the exchange-only Rys pass kept "// &
                       to_string(nmd)//" quartets, the fused pass "//to_string(nfu))
            if (allocated(error)) return
         end do
         deallocate (d)
      end do
   end subroutine test_mdj_switch

   subroutine md_and_fused(ctx, d, kfac, jfac, omega, gmd, gfu, nmd, nfu)
      !! One call per engine; nmd and nfu the Rys survivors of each.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: d(:, :, :), kfac, jfac, omega
      real(dp), intent(out) :: gmd(:, :), gfu(:, :)
      integer(int64), intent(out) :: nmd, nfu
      ctx%eri%grad_j_engine = TRC_J_ENGINE_MDJ
      call ctx%eri%gradient(ctx%bas, d, gmd, k_scale=kfac, j_scale=jfac, omega=omega)
      nmd = ctx%eri%nkept
      ctx%eri%grad_j_engine = TRC_J_ENGINE_BINS
      call ctx%eri%gradient(ctx%bas, d, gfu, k_scale=kfac, j_scale=jfac, omega=omega)
      nfu = ctx%eri%nkept
   end subroutine md_and_fused

end module check_grad_screen_suite

program check_grad_screen
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_grad_screen_suite, only: collect_grad_screen_tests, release_contexts
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("grad_screen", collect_grad_screen_tests)]
   do is = 1, size(testsuites)
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   call release_contexts()
   if (stat > 0) error stop 1
end program check_grad_screen
