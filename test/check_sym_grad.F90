!
! The skeleton nuclear gradient: the S16 oracle.
!
! `eri%gradient(..., skeleton=.true.)` launches only the orbit
! representatives of the cheaper side of each admitted bin pair, weighted by
! their orbit size, so what it returns is grad_skel and NOT a gradient. The
! full gradient is the projection
!
!   grad = (1/|G|) sum_g P_g grad_skel,   (P_g X)(:, a) = s_g X(:, perm(a, g))
!
! which is `sym_scf_grad_symmetrize` of trc_sym_scf, the same routine
! trc_scf_gradient already applied to the finished gradient. This check
! builds both on the SAME engine (bins closed under the group,
! `eri%build(..., sym=geom)`; the full build walks the same quartets in
! another order) and compares them:
!
!   |sym(grad_skel) - grad_full|_max  <=  TOL_REL * max(1, |grad_full|_max)
!
! Unlike the Fock skeleton, the projection here is not a tidy-up: the raw
! skeleton gradient is wrong by order one (0.58 to 1.84 against components of
! 0.27 to 2.14), so `raw_case` asserts
! that it IS wrong -- a test that accidentally compared two full builds would
! pass everything else.
!
! The axes, per molecule:
!   - Cartesian and spherical caller AO spaces. The density is made invariant
!     by the production symmetrizer (`trc_sym_scf_t%symmetrize`) in the
!     caller's space and folded to the Cartesian working basis the integrals
!     are formed in -- section 10 of the design records that the spike's own
!     hand-rolled prototype projection was NOT invariant in the spherical space, which
!     is why both spaces are swept here;
!   - J+K, J only (k_scale = 0), K only (j_scale = 0) and long range
!     (j_scale = 0, k_scale, omega > 0);
!   - RHF (ndens = 1, the total density) and UHF (ndens = 2, two spin
!     densities), which must agree with each other as well as with their own
!     full build;
!   - a random symmetrized density and the symmetric SCF's converged one;
!   - the range-separated hybrid's TWO eri%gradient passes, summed and
!     projected ONCE;
!   - the work and survivor counts (`eri%nwork`, `eri%nkept`): the skeleton
!     plan must be SMALLER. A silent fall back to the full plan is otherwise
!     invisible, and it is the failure mode that made parity skipping look
!     active when it was not;
!   - translational invariance of the RAW skeleton gradient, which holds on
!     its own (every quartet's contribution is translation invariant) and so
!     does not depend on the projection;
!   - the production path end to end: trc_scf_gradient with `sym_geom`, whose
!     four-centre term is a skeleton build projected inside
!     two_electron_4c, against the same molecule's C1 gradient, with the
!     components the group forbids and the equivariance of the result;
!   - and the same path with the skeleton DECLINED by the size rule
!     (`sym_skeleton`), where the four-centre term is the full C1 build and
!     the end-of-routine projection is the only one left. That case pins the
!     ROUTING, not the numbers: the engine the SCF already built is reused
!     instead of a second one being built with the group;
!   - a central finite difference of terco's own reconverged SCF energy
!     against that production gradient -- the only check here that would
!     catch an error consistent between the skeleton and the full path.
!
! Registered as check_sym_grad (mode "main"), check_sym_grad_cap (mode
! "quick" with TRC_BIN_CAP lowered, the orbit-closed bin split that
! `bins_plan_emit`'s canonicalisation exists for) and check_sym_grad_f (mode
! "f": the uncontracted l = 0..3 set with TRC_BIN_CAP lowered, which prints
! SKIP below TERCO_LMAX = 3). That last combination -- an f basis, a split
! bin and a skeleton plan that cuts the cheaper side -- is the one that can
! present the gradient class filter with an ordering it has no kernel for,
! and the one that proves it cannot (src/CMakeLists.txt's
! _rys_grad_classes).
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_grad_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_bins, only: trc_bin_plan_t
   use trc_eri, only: trc_eri_t
   use trc_pc_kernels, only: CLASS_RADIX
   use trc_rys_grad_kernels, only: rysg_has_class
   use trc_error, only: error_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run, trc_pair_eri_cache_t
   use trc_scf_grad, only: trc_scf_gradient, trc_scf_grad_sym_t
   use trc_sph, only: sph_sph_to_cart
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_scf, only: trc_sym_scf_t, sym_scf_grad_symmetrize
   use trc_test_basis, only: read_xyz, build_aux
   implicit none
   private
   public :: collect_sym_grad, g_mode

   character(len=8) :: g_mode = "main"
      !! "main", "quick" or "f", set from argv before run_testsuite.

   real(dp), parameter :: THRESH = 1.0e-10_dp
   real(dp), parameter :: TOL_REL = 1.0e-10_dp
      !! Allowed |sym(grad_skel) - grad_full| relative to
      !! max(1, |grad_full|_max). The design prototype measured 3.1e-15 (water
      !! C2v), 1.7e-14 (B2H6 D2h), 2.1e-12 (SF6 D2h) and 1.0e-11 (benzene)
      !! relative, the last being the double-precision summation noise of
      !! ~9.4M accumulated segment contributions in a prototype that
      !! accumulated about TWICE the full build's work items (it applied the
      !! orbit weight by repetition); a real weighted kernel accumulates
      !! 3.5x FEWER than the full build. 1e-10 leaves a factor of ten on the
      !! noisiest measured case. A symmetry error is O(1) here, not marginal:
      !! raw_case measures the unprojected skeleton at 0.3 to 2.
   real(dp), parameter :: TOL_INV = 1.0e-10_dp
      !! Translational invariance of the RAW skeleton gradient, the same gate
      !! check_full_grad_fd and check_rys_grad_fd use for the full one.
   real(dp), parameter :: TOL_C1 = 1.0e-9_dp
      !! |production symmetric gradient - C1 gradient|, as the per-group
      !! check_sym_* gradient cases use.
   real(dp), parameter :: TOL_EQV = 1.0e-14_dp
      !! Equivariance and forbidden components of a projected gradient.
   real(dp), parameter :: FD_H = 1.0e-3_dp
   real(dp), parameter :: FD_RELTOL = 1.0e-5_dp
   real(dp), parameter :: FD_ABSTOL = 5.0e-7_dp
      !! check_full_grad_fd's own gates, for the same reason: at H = 1e-3 the
      !! central difference of a reconverged SCF energy has a ~1e-7 absolute
      !! truncation floor.

   integer(int64), parameter :: WORK_WATER(2) = [3081_int64, 1868_int64]
      !! Water/cc-pVDZ, full and skeleton work items of the gradient plan, on
      !! the default split-when-it-pays basis (14 627 and 9 968 are the
      !! fully split basis, TRC_DECON_SHELL_COST=0).
      !! A 1.65x reduction; Cartesian and spherical alike.
   integer(int64), parameter :: WORK_B2H6(2) = [45108_int64, 12690_int64]
      !! B2H6/6-31G*, the same way: a 3.55x reduction.

   type :: ctx_t
      !! One molecule in one AO space, with its symmetric engine.
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_sym_scf_t) :: sc
      type(trc_eri_t) :: eri
      integer, allocatable :: elem(:)
      logical :: sph = .false.
      integer :: nc = 0, ns = 0
      character(len=96) :: label = " "
      real(dp) :: worst = 0.0_dp
      real(dp) :: worst_raw = 0.0_dp
      integer(int64) :: want_full = -1, want_skel = -1
      !! When set (>= 0), the exact work items of the full and the skeleton
      !! plan: `<` alone would let the reduction regress unnoticed.
   end type ctx_t

contains

   subroutine collect_sym_grad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == "f") then
         testsuite = [ &
                     new_unittest("uncontracted l = 0..3, split bins: water C2v, B2H6 D2h", test_lmax3) &
                     ]
      else if (g_mode == "quick") then
         testsuite = [ &
                     new_unittest("water C2v, B2H6 D2h: skeleton vs full gradient", test_quick) &
                     ]
      else
         testsuite = [ &
                     new_unittest("Cs CHFCl2", test_cs), &
                     new_unittest("Cs HOCl: every pair orbit a singleton", test_cs_planar), &
                     new_unittest("C2 H2O2", test_c2), &
                     new_unittest("C2v water", test_c2v), &
                     new_unittest("C2h trans-N2H2", test_c2h), &
                     new_unittest("D2 twisted ethylene", test_d2), &
                     new_unittest("Ci N2H2F2", test_ci), &
                     new_unittest("D2h B2H6 (nonplanar)", test_b2h6), &
                     new_unittest("D2h SF6", test_sf6), &
                     new_unittest("D6h benzene as D2h", test_benzene), &
                     new_unittest("the raw skeleton gradient is NOT a gradient", test_raw), &
                     new_unittest("production path: skeleton + projection vs C1", test_production), &
                     new_unittest("a DECLINED skeleton reuses the SCF's engine and still projects", test_declined), &
                     new_unittest("finite difference of the reconverged SCF energy", test_fd) &
                     ]
      end if
   end subroutine collect_sym_grad

   ! ------------------------------------------------------------------ setup

   subroutine setup_ctx(ctx, xyz, group, basis_json, auxl, sph, error)
      !! Detect the group, build the basis on its standard frame (from a JSON
      !! file, or the uncontracted l = 0..auxl set when auxl >= 0), the
      !! production symmetrizer of the caller's AO space, and the engine with
      !! the point group (which is what the skeleton plan needs).
      type(ctx_t), intent(inout) :: ctx
      character(len=*), intent(in) :: xyz, group, basis_json
      integer, intent(in) :: auxl
      logical, intent(in) :: sph
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), ident(:, :), zat(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer, allocatable :: sh_l(:), sh_np(:)
      integer :: nat, i, nsh, maxnp

      call read_xyz(trim(xyz)//".xyz", nat, ctx%elem, r)
      ctx%geom = sym_geom_detect(r, ctx%elem, err, basis_id=ctx%elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, ctx%geom%group%name, trim(group), "group of "//trim(xyz))
      if (allocated(error)) return
      if (auxl >= 0) then
         call build_aux(nat, ctx%elem, ctx%geom%r_std, auxl, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
         allocate (zat(nat))
         zat = real(ctx%elem, dp)
         call ctx%bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, ctx%geom%r_std, maxnp)
      else
         call trc_basis_from_json(trim(basis_json), nat, ctx%elem, ctx%geom%r_std, ctx%bas, err, force_cart=.not. sph)
         call check(error, .not. err%has_error(), "basis: "//err%get_message())
         if (allocated(error)) return
      end if
      ctx%bas%spherical = sph
      ctx%sph = sph
      call ctx%bas%to_device()
      ctx%nc = ctx%bas%nao
      ctx%ns = merge(ctx%bas%nao_sph, ctx%bas%nao, sph)
      ctx%label = trim(xyz)//" "//trim(group)//merge(" sph ", " cart", sph)

      ctx%ao = sym_ao_build(ctx%bas, ctx%geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      ! The production symmetrizer needs its device data, made by `prepare`
      ! from an overlap; the identity gives every SO block a unit metric.
      ctx%sc%ao = ctx%ao
      allocate (ident(ctx%ns, ctx%ns))
      ident = 0.0_dp
      do i = 1, ctx%ns
         ident(i, i) = 1.0_dp
      end do
      call ctx%sc%prepare(ident, 1, 0.0_dp, err)
      call check(error, .not. err%has_error(), "sym_scf prepare: "//err%get_message())
      if (allocated(error)) return

      call ctx%eri%build(ctx%bas, THRESH, sym=ctx%geom)
   end subroutine setup_ctx

   subroutine release_ctx(ctx)
      type(ctx_t), intent(inout) :: ctx
      call ctx%eri%release()
      call ctx%sc%release()
      call ctx%bas%release()
   end subroutine release_ctx

   ! -------------------------------------------------------------- densities

   subroutine random_density(ctx, seed, scale, ds)
      !! A symmetric matrix of pseudo-random numbers in the caller's AO space,
      !! projected onto the group by the production symmetrizer. A linear
      !! congruential generator, so the test does not depend on the
      !! compiler's random_number.
      type(ctx_t), intent(inout) :: ctx
      integer, intent(in) :: seed
      real(dp), intent(in) :: scale
      real(dp), intent(out) :: ds(ctx%ns, ctx%ns)
      real(dp), allocatable :: w(:, :)
      integer(int64) :: st
      integer :: i, j, n
      n = ctx%ns
      st = int(seed, int64)*2654435761_int64 + 12345_int64
      do j = 1, n
         do i = 1, j
            st = mod(st*6364136223846793005_int64 + 1442695040888963407_int64, 2_int64**62)
            ds(i, j) = scale*(real(mod(st/1024_int64, 20001_int64), dp)/10000.0_dp - 1.0_dp)
            ds(j, i) = ds(i, j)
         end do
      end do
      allocate (w(n, n))
      call ctx%sc%symmetrize(n, ds, w)
   end subroutine random_density

   subroutine scf_density(ctx, ds, error)
      !! The converged density of the symmetric SCF in this basis.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(out) :: ds(ctx%ns, ctx%ns)
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: rs
      type(trc_sym_scf_t) :: scr
      integer :: nocc
      nocc = sum(ctx%elem)/2
      ! Just above the accumulation-order noise floor:
      ! the OpenMP and device builds accumulate F in a run-dependent order, so identical runs differ by 1e-13 to 1e-12 Eh.
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      scr%ao = ctx%ao
      call trc_scf_run(ctx%bas, nocc, nocc, opts, rs, sym=scr)
      call check(error, rs%converged, trim(ctx%label)//": the SCF did not converge: "//trim(rs%message))
      if (allocated(error)) return
      call scr%release()
      ds = rs%dmat(:, :, 1)
   end subroutine scf_density

   subroutine to_cart(ctx, ds, dc)
      !! The invariant density, in the Cartesian basis the integrals are
      !! formed in. The fold is where the spherical axis differs from the
      !! Cartesian one at all.
      type(ctx_t), intent(in) :: ctx
      real(dp), intent(in) :: ds(ctx%ns, ctx%ns)
      real(dp), intent(out) :: dc(ctx%nc, ctx%nc)
      if (ctx%sph) then
         call sph_sph_to_cart(ctx%bas, ds, dc)
      else
         dc = ds
      end if
   end subroutine to_cart

   ! ------------------------------------------------------------ comparisons

   subroutine equivariance(geom, g, asym, forbidden)
      !! g(image) = s g(atom) for every operation, s the sign of the
      !! component; on an atom that is its own image a component with
      !! s = -1 must vanish.
      type(trc_sym_geom_t), intent(in) :: geom
      real(dp), intent(in) :: g(:, :)
      real(dp), intent(out) :: asym, forbidden
      integer :: a, k, ic
      real(dp) :: sg
      asym = 0.0_dp
      forbidden = 0.0_dp
      do a = 1, size(g, 2)
         do k = 2, geom%group%order
            do ic = 1, 3
               sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
               asym = max(asym, abs(g(ic, geom%perm(a, k)) - sg*g(ic, a)))
               if (geom%perm(a, k) == a .and. sg < 0.0_dp) forbidden = max(forbidden, abs(g(ic, a)))
            end do
         end do
      end do
   end subroutine equivariance

   subroutine canon_case(ctx, reordered, error)
      !! PLAN LEVEL, no integrals, so it costs nothing. Both plans the
      !! gradient can build must be ANGULAR-MOMENTUM CANONICAL --
      !! (sLA, sLB) >= (sLC, sLD) on every segment -- and every class they
      !! carry must be one this build generated a gradient kernel for. That is
      !! the whole claim of bins_plan_emit's `canon`, and above d it is the
      !! difference between a launch and grad_classes' `error stop`: the
      !! filter in src/CMakeLists.txt builds no f-containing class with
      !! (lc, ld) > (la, lb) (30 of them, listed there).
      !!
      !! The third count is what makes this test non-vacuous: a segment whose
      !! sub-range sits on the KET side (sNB < sSB) is one where
      !! bins_skel_cut_first chose the far bin -- the cheaper side -- and the
      !! canonicalisation then moved that sub-range to the ket rather than
      !! swapping the class. Without `canon` it would be the bra side and the
      !! count would be zero, which is the assertion that fails without the
      !! change. CLASS_RADIX is trc_grad's own GCLASS_RADIX (both 11).
      type(ctx_t), intent(inout) :: ctx
      logical, intent(in) :: reordered
      !! At least one segment must have had its two sides reordered. False
      !! for a molecule all of whose pair orbits are singletons (a planar Cs
      !! molecule, see `test_cs_planar`): every section then covers its whole
      !! bin, so no segment is a sub-range at all and there is nothing to
      !! order.
      type(error_type), allocatable, intent(out) :: error
      type(trc_bin_plan_t) :: pf, ps
      integer :: i, nbad_f, nbad_s, nmiss_f, nmiss_s, nket, nbra

      ! Both plans as grad_bins builds them: only the gradient asks for the
      ! canonical order, the energy path keeps the walk order.
      call pf%build(ctx%eri%bins, ctx%eri%thresh, ctx%eri%quartets%bin_pbase, .false., canonical=.true.)
      call ps%build_skeleton(ctx%eri%bins, ctx%eri%quartets%orb, ctx%eri%thresh, ctx%eri%quartets%bin_pbase, &
                             canonical=.true.)
      nbad_f = 0
      nbad_s = 0
      nmiss_f = 0
      nmiss_s = 0
      nket = 0
      nbra = 0
      do i = 1, pf%nseg
         if (pf%sLC(i)*16 + pf%sLD(i) > pf%sLA(i)*16 + pf%sLB(i)) nbad_f = nbad_f + 1
         if (.not. rysg_has_class(pf%class_key(i, CLASS_RADIX))) nmiss_f = nmiss_f + 1
      end do
      do i = 1, ps%nseg
         if (ps%sLC(i)*16 + ps%sLD(i) > ps%sLA(i)*16 + ps%sLB(i)) nbad_s = nbad_s + 1
         if (.not. rysg_has_class(ps%class_key(i, CLASS_RADIX))) nmiss_s = nmiss_s + 1
         if (ps%sNB(i) < ps%sSB(i)) nket = nket + 1
         if (ps%sA(i) < ps%sSA(i)) nbra = nbra + 1
      end do
      print "(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)", "  ", trim(ctx%label), &
         " plans: full segs ", pf%nseg, " non-canonical ", nbad_f, " without a kernel ", nmiss_f, &
         ";  skeleton segs ", ps%nseg, " non-canonical ", nbad_s, " without a kernel ", nmiss_s, &
         ";  sub-range on the ket ", nket
      call pf%release()
      call ps%release()

      call check(error, nbad_f == 0, trim(ctx%label)//": "//to_string(nbad_f)// &
                 " segments of the FULL gradient plan have (lc, ld) > (la, lb)")
      if (allocated(error)) return
      call check(error, nbad_s == 0, trim(ctx%label)//": "//to_string(nbad_s)// &
                 " segments of the SKELETON gradient plan have (lc, ld) > (la, lb)")
      if (allocated(error)) return
      call check(error, nmiss_f == 0 .and. nmiss_s == 0, trim(ctx%label)//": "// &
                 to_string(nmiss_f)//" full and "//to_string(nmiss_s)// &
                 " skeleton segments carry a class with no gradient kernel")
      if (allocated(error)) return
      if (.not. reordered) return
      call check(error, nket > 0, trim(ctx%label)// &
                 ": no skeleton segment has its sub-range on the ket side, so the cheaper-side cut "// &
                 "never had to be reordered and this check proves nothing")
      if (allocated(error)) return
      call check(error, nbra > 0, trim(ctx%label)//": no skeleton segment has its sub-range on the bra side")
   end subroutine canon_case

   subroutine compare_grad(ctx, dc, ndens, jscale, kscale, omega, tag, smaller, error, graw, gfull)
      !! The skeleton and the full gradient of the same operator on the same
      !! engine: the first projected over the group, the second as it comes.
      !! `graw`/`gfull` come back for the raw-difference and the
      !! two-pass cases.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      integer, intent(in) :: ndens
      real(dp), intent(in) :: jscale, kscale, omega
      character(len=*), intent(in) :: tag
      logical, intent(in) :: smaller
      !! The skeleton plan must have FEWER work items than the full one. False
      !! for a molecule all of whose pair orbits are singletons -- see
      !! `singleton_case`, where it legitimately has more.
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(out), optional :: graw(3, ctx%bas%natm), gfull(3, ctx%bas%natm)
      real(dp), allocatable :: dw(:, :, :), gf(:, :), gs(:, :), g0(:, :)
      real(dp) :: err_sym, err_raw, gmax, inv
      integer(int64) :: nw_full, nw_skel, nk_full, nk_skel
      integer :: natm, id

      natm = ctx%bas%natm
      allocate (dw(ctx%nc, ctx%nc, ndens), gf(3, natm), gs(3, natm), g0(3, natm))
      ! ndens = 2 splits the invariant total density into two equal spin
      ! densities: eri%gradient(Da, Da) must equal eri%gradient(2 Da), which
      ! is what fixes the UHF exchange coefficient (gen_rys_grad.py's
      ! kfac_eff). Both are symmetric under the group because dc is.
      do id = 1, ndens
         dw(:, :, id) = dc/real(ndens, dp)
      end do

      call ctx%eri%gradient(ctx%bas, dw, gf, k_scale=kscale, j_scale=jscale, omega=omega)
      nw_full = ctx%eri%nwork
      nk_full = ctx%eri%nkept
      call ctx%eri%gradient(ctx%bas, dw, gs, k_scale=kscale, j_scale=jscale, omega=omega, skeleton=.true.)
      nw_skel = ctx%eri%nwork
      nk_skel = ctx%eri%nkept

      gmax = maxval(abs(gf))
      err_raw = maxval(abs(gs - gf))
      ! Translational invariance of the RAW skeleton: every quartet's force
      ! contribution sums to zero over its own four centres, so the
      ! unprojected weighted sum does too, whatever the weights are.
      inv = maxval(abs(sum(gs, dim=2)))
      g0 = gs
      call sym_scf_grad_symmetrize(ctx%geom%group, ctx%geom%perm, gs)
      err_sym = maxval(abs(gs - gf))
      ctx%worst = max(ctx%worst, err_sym/max(1.0_dp, gmax))
      ctx%worst_raw = max(ctx%worst_raw, err_raw/max(1.0_dp, gmax))
      print "(a,a,a,a,a,i0,a,es10.3,a,es10.3,a,es9.2,a,es8.1,a,i0,a,i0,a,i0,a,i0)", &
         "  ", trim(ctx%label), " ", trim(tag), " nd", ndens, &
         ": max|g| ", gmax, "  raw err ", err_raw, "  sym err ", err_sym, "  trans ", inv, &
         "  work ", nw_full, " -> ", nw_skel, "  kept ", nk_full, " -> ", nk_skel

      call check(error, err_sym <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)//" "//trim(tag)// &
                 ": projected skeleton gradient differs from the full one by "//to_string(err_sym, "es10.3")// &
                 " (max|g| "//to_string(gmax, "es10.3")//")")
      if (allocated(error)) return
      if (ctx%want_full >= 0) then
         call check(error, nw_full == ctx%want_full .and. nw_skel == ctx%want_skel, trim(ctx%label)//" "// &
                    trim(tag)//": work items "//to_string(nw_full)//" -> "//to_string(nw_skel)//", pinned "// &
                    to_string(ctx%want_full)//" -> "//to_string(ctx%want_skel))
         if (allocated(error)) return
      end if
      ! A silent fall back to the full plan would leave both counts equal.
      if (smaller) then
         call check(error, nw_skel < nw_full, trim(ctx%label)//" "//trim(tag)// &
                    ": the skeleton plan is not smaller than the full one ("//to_string(nw_skel)//" vs "// &
                    to_string(nw_full)//" work items)")
         if (allocated(error)) return
         call check(error, nk_skel < nk_full .and. nk_skel > 0_int64, trim(ctx%label)//" "//trim(tag)// &
                    ": the skeleton screen kept "//to_string(nk_skel)//" quartets against the full "// &
                    to_string(nk_full))
         if (allocated(error)) return
      else
         call check(error, nk_skel > 0_int64, trim(ctx%label)//" "//trim(tag)//": the screen kept nothing")
         if (allocated(error)) return
      end if
      call check(error, inv <= TOL_INV, trim(ctx%label)//" "//trim(tag)// &
                 ": the raw skeleton gradient is not translation invariant: "//to_string(inv, "es10.2"))
      if (allocated(error)) return
      if (present(graw)) graw = g0
      if (present(gfull)) gfull = gf
   end subroutine compare_grad

   subroutine compare_two_pass(ctx, dc, error)
      !! A range-separated hybrid runs eri%gradient TWICE (the full-range
      !! term, then j_scale = 0 with k_scale = exx_lr and omega > 0) and adds
      !! the two. Both passes take the same plan, so both are skeleton, and
      !! the projector is applied ONCE to their sum -- which is what
      !! two_electron_4c does. Linearity of the projector is the claim.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: dw(:, :, :), ga(:, :), gb(:, :), gf(:, :)
      real(dp) :: e, gmax
      integer :: natm
      natm = ctx%bas%natm
      allocate (dw(ctx%nc, ctx%nc, 1), ga(3, natm), gb(3, natm), gf(3, natm))
      dw(:, :, 1) = dc
      call ctx%eri%gradient(ctx%bas, dw, gf, k_scale=0.25_dp)
      call ctx%eri%gradient(ctx%bas, dw, ga, j_scale=0.0_dp, k_scale=0.6_dp, omega=0.3_dp)
      gf = gf + ga
      call ctx%eri%gradient(ctx%bas, dw, ga, k_scale=0.25_dp, skeleton=.true.)
      call ctx%eri%gradient(ctx%bas, dw, gb, j_scale=0.0_dp, k_scale=0.6_dp, omega=0.3_dp, skeleton=.true.)
      ga = ga + gb
      call sym_scf_grad_symmetrize(ctx%geom%group, ctx%geom%perm, ga)
      gmax = maxval(abs(gf))
      e = maxval(abs(ga - gf))
      ctx%worst = max(ctx%worst, e/max(1.0_dp, gmax))
      print "(a,a,a,es10.3,a,es10.3)", "  ", trim(ctx%label), " range-separated, two passes one projection: max|g| ", &
         gmax, "  err ", e
      call check(error, e <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)// &
                 ": the two skeleton passes, projected once, differ from the two full ones by "// &
                 to_string(e, "es10.3"))
   end subroutine compare_two_pass

   subroutine run_densities(ctx, scf, smaller, error)
      !! The operator and spin axes on one or two densities.
      type(ctx_t), intent(inout) :: ctx
      logical, intent(in) :: scf, smaller
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ds(:, :), dc(:, :), g1(:, :), g2(:, :)
      real(dp) :: e
      integer :: id
      character(len=96) :: base

      allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc), g1(3, ctx%bas%natm), g2(3, ctx%bas%natm))
      base = ctx%label
      do id = 1, merge(2, 1, scf)
         if (id == 1) then
            call random_density(ctx, 7, 0.5_dp, ds)
            ctx%label = trim(base)//" rand"
         else
            call scf_density(ctx, ds, error)
            if (allocated(error)) return
            ctx%label = trim(base)//" scf"
         end if
         call to_cart(ctx, ds, dc)
         call compare_grad(ctx, dc, 1, 1.0_dp, 1.0_dp, 0.0_dp, "J+K ", smaller, error, gfull=g1)
         if (allocated(error)) return
         call compare_grad(ctx, dc, 2, 1.0_dp, 1.0_dp, 0.0_dp, "J+K ", smaller, error, gfull=g2)
         if (allocated(error)) return
         ! D = 2 Da is the RHF solution, so the two must be the same gradient.
         e = maxval(abs(g1 - g2))
         call check(error, e <= TOL_REL*max(1.0_dp, maxval(abs(g1))), trim(ctx%label)// &
                    ": the UHF (ndens = 2) gradient of D/2, D/2 differs from the RHF gradient of D by "// &
                    to_string(e, "es10.3"))
         if (allocated(error)) return
         call compare_grad(ctx, dc, 1, 1.0_dp, 0.0_dp, 0.0_dp, "J   ", smaller, error)
         if (allocated(error)) return
         call compare_grad(ctx, dc, 1, 0.0_dp, 1.0_dp, 0.0_dp, "K   ", smaller, error)
         if (allocated(error)) return
         call compare_grad(ctx, dc, 1, 0.0_dp, 0.6_dp, 0.3_dp, "K_lr", smaller, error)
         if (allocated(error)) return
         call compare_grad(ctx, dc, 2, 0.0_dp, 0.6_dp, 0.3_dp, "K_lr", smaller, error)
         if (allocated(error)) return
         if (id == 1) then
            call compare_two_pass(ctx, dc, error)
            if (allocated(error)) return
         end if
      end do
      ctx%label = base
   end subroutine run_densities

   subroutine run_molecule(xyz, group, basis_json, scf, smaller, error, work)
      !! One molecule, Cartesian and spherical.
      character(len=*), intent(in) :: xyz, group, basis_json
      logical, intent(in) :: scf, smaller
      type(error_type), allocatable, intent(out) :: error
      integer(int64), intent(in), optional :: work(2)
      !! The expected work items, full and skeleton (`ctx_t%want_full`).
      type(ctx_t) :: ctx
      integer :: ik, dev
      dev = trc_bind_device(0)
      do ik = 1, 2
         call setup_ctx(ctx, xyz, group, basis_json, -1, ik == 2, error)
         if (allocated(error)) return
         if (present(work)) then
            ctx%want_full = work(1)
            ctx%want_skel = work(2)
         end if
         ! RELEASE ON EVERY PATH OUT. An early return would leave the basis
         ! and the engine mapped; the next molecule's `to_device` then meets a
         ! host address the device present table still holds, with the wrong
         ! size, and nvfortran dumps the present table and dies. Measured:
         ! the GPU build's first failing case took the whole program with it.
         call canon_case(ctx, smaller, error)
         if (.not. allocated(error)) call run_densities(ctx, scf, smaller, error)
         if (.not. allocated(error)) then
            print "(a,a,a,es9.2,a,es9.2)", "  ", trim(ctx%label), ": worst relative error ", ctx%worst, &
               ", worst RAW relative error ", ctx%worst_raw
         end if
         call release_ctx(ctx)
         if (allocated(error)) return
      end do
   end subroutine run_molecule

   ! ------------------------------------------------------------------ tests

   subroutine test_cs(error)
      !! CHFCl2: the mirror EXCHANGES the two chlorines, so the pair orbits
      !! are not all singletons and the skeleton plan is smaller. See
      !! `test_cs_planar` for the molecule where they are.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("chfcl2", "Cs", "basis_sets/6-31gs.json", .false., .true., error)
   end subroutine test_cs

   subroutine test_cs_planar(error)
      !! A PLANAR Cs molecule, where the only operation besides the identity
      !! is the molecular plane and every atom is therefore its own image. No
      !! shell pair has an image but itself, so every pair orbit is a
      !! singleton, every segment weight is one, and the skeleton restriction
      !! drops NOTHING: the raw skeleton gradient already equals the full one
      !! and the plan is slightly LARGER, because a diagonal bin pair becomes
      !! the rectangle reps(A) x A with weight 1/2 instead of the triangle
      !! (measured on HOCl/6-31G*: 9355 work items against 9142, +2.3%). The
      !! skeleton is then a small pessimisation, exactly as it is for the
      !! skeleton FOCK build on such a molecule -- which is why ONE size
      !! rule now governs both builds, and why this molecule is in
      !! `test_declined` as well: the rule declines it (3 atoms in 3
      !! orbits), so the production path never runs the larger plan. This
      !! case drives `run_molecule`, which asks for the skeleton directly
      !! and so reaches it anyway. Pinned here so it is a documented
      !! property and not a surprise, and because it is the one
      !! configuration in which the "the plan must be smaller" assertion
      !! must NOT hold.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("hocl", "Cs", "basis_sets/6-31gs.json", .false., .false., error)
   end subroutine test_cs_planar

   subroutine test_c2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("h2o2", "C2", "basis_sets/cc-pvdz.json", .true., .true., error)
   end subroutine test_c2

   subroutine test_c2v(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("water", "C2v", "basis_sets/cc-pvdz.json", .true., .true., error, WORK_WATER)
   end subroutine test_c2v

   subroutine test_c2h(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("n2h2_trans", "C2h", "basis_sets/6-31gs.json", .true., .true., error)
   end subroutine test_c2h

   subroutine test_d2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("c2h4_twisted", "D2", "basis_sets/6-31gs.json", .true., .true., error)
   end subroutine test_d2

   subroutine test_ci(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("n2h2f2_ci", "Ci", "basis_sets/6-31gs.json", .false., .true., error)
   end subroutine test_ci

   subroutine test_b2h6(error)
      !! Nonplanar D2h: every operation MOVES atoms, so a sign or a
      !! permutation error cannot cancel the way it can on a planar molecule.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("b2h6_d2h", "D2h", "basis_sets/6-31gs.json", .false., .true., error, WORK_B2H6)
   end subroutine test_b2h6

   subroutine test_sf6(error)
      !! The one system the design measured where the orbit sizes are mixed
      !! (2.42x rather than ~3.8x), so several weights occur in one plan.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("tops_sf6", "D2h", "basis_sets/cc-pvdz.json", .false., .true., error)
   end subroutine test_sf6

   subroutine test_benzene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("benzene", "D2h", "basis_sets/6-31gs.json", .false., .true., error)
   end subroutine test_benzene

   subroutine test_quick(error)
      type(error_type), allocatable, intent(out) :: error
      ! Not pinned here: check_sym_grad_cap runs this with TRC_BIN_CAP=40, whose
      ! split bins admit a different count (test_c2v and test_b2h6 pin them).
      call run_molecule("water", "C2v", "basis_sets/cc-pvdz.json", .false., .true., error)
      if (allocated(error)) return
      call run_molecule("b2h6_d2h", "D2h", "basis_sets/6-31gs.json", .false., .true., error)
   end subroutine test_quick

   subroutine test_lmax3(error)
      !! The uncontracted l = 0..3 set, with TRC_BIN_CAP lowered by the ctest
      !! registration so that split_large_bins fires. An f basis plus a split
      !! bin plus a skeleton plan that cuts the cheaper side is the only
      !! combination that can present grad_classes with a class ordering the
      !! filter drops, and this is the test that would `error stop` -- or,
      !! before reaching a launch at all, fail canon_case -- without
      !! bins_plan_emit's canonicalisation.
      !!
      !! The PLAN checks run on both molecules in both AO spaces, because they
      !! cost nothing. The integral comparison runs on water, Cartesian, J+K
      !! only: an f-class gradient over a bin split at 40 pairs is minutes of
      !! host time per evaluation, and what it adds over the plan checks is
      !! only that the kernels index a sub-range correctly, which every case
      !! in the main mode already covers.
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      integer :: im, ik, dev
      character(len=12), parameter :: MOL(2) = [character(len=12) :: "water", "b2h6_d2h"]
      character(len=4), parameter :: GRP(2) = [character(len=4) :: "C2v", "D2h"]
      dev = trc_bind_device(0)
      do im = 1, 2
         do ik = 1, 2
            call setup_ctx(ctx, trim(MOL(im)), trim(GRP(im)), " ", 3, ik == 2, error)
            if (allocated(error)) return
            call canon_case(ctx, .true., error)
            if (.not. allocated(error) .and. im == 1 .and. ik == 1) then
               allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc))
               call random_density(ctx, 3, 0.1_dp, ds)
               call to_cart(ctx, ds, dc)
               call compare_grad(ctx, dc, 1, 1.0_dp, 1.0_dp, 0.0_dp, "J+K ", .true., error)
               deallocate (ds, dc)
            end if
            call release_ctx(ctx)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_lmax3

   subroutine test_raw(error)
      !! The raw skeleton gradient must be WRONG. Without this a test that
      !! accidentally compared two full builds would pass everything above:
      !! the projection is not a correction here, it is the whole of the
      !! symmetry (measured: 0.58 to 1.84
      !! against components of 0.27 to 2.14).
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :), graw(:, :), gfull(:, :)
      real(dp) :: rel
      integer :: dev
      dev = trc_bind_device(0)
      call setup_ctx(ctx, "b2h6_d2h", "D2h", "basis_sets/6-31gs.json", -1, .false., error)
      if (allocated(error)) return
      allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc), graw(3, ctx%bas%natm), gfull(3, ctx%bas%natm))
      call random_density(ctx, 7, 0.5_dp, ds)
      call to_cart(ctx, ds, dc)
      call compare_grad(ctx, dc, 1, 1.0_dp, 1.0_dp, 0.0_dp, "J+K ", .true., error, graw=graw, gfull=gfull)
      if (.not. allocated(error)) then
         rel = maxval(abs(graw - gfull))/max(1.0_dp, maxval(abs(gfull)))
         print "(a,es10.3,a,es10.3)", "  raw skeleton vs full, relative: ", rel, "  max|g| ", maxval(abs(gfull))
         call check(error, rel > 0.1_dp, "the raw skeleton gradient is within 10% of the full one ("// &
                    to_string(rel, "es10.3")//"): the projection is not being tested")
      end if
      call release_ctx(ctx)
   end subroutine test_raw

   subroutine production_case(error, xyz, basis_json, spherical, group)
      !! trc_scf_gradient with `sym_geom`: the four-centre term is a skeleton
      !! build projected inside two_electron_4c, with the one-electron,
      !! overlap and nuclear-repulsion terms on top and the end-of-routine
      !! projection as well. Against the same molecule run in C1.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, group
      logical, intent(in) :: spherical
      type(trc_basis_t) :: bas, bas1
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_sym_scf_t) :: sc
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: rs, r1
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :), gsym(:, :), g1(:, :)
      real(dp) :: dg, asym, forbidden
      integer :: nat, nocc, dev

      dev = trc_bind_device(0)
      call read_xyz(trim(xyz)//".xyz", nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, trim(group), "group of "//trim(xyz))
      if (allocated(error)) return
      call trc_basis_from_json(trim(basis_json), nat, elem, geom%r_std, bas, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = spherical
      call bas%to_device()
      nocc = sum(elem)/2
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100

      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      sc%ao = ao
      call trc_scf_run(bas, nocc, nocc, opts, rs, sym=sc)
      call check(error, rs%converged, trim(xyz)//": the symmetric SCF did not converge")
      if (allocated(error)) return
      ! A second basis object on the SAME standard frame for the C1 run: the
      ! geometry must be identical, only the symmetry handling differs.
      call trc_basis_from_json(trim(basis_json), nat, elem, geom%r_std, bas1, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas1%spherical = spherical
      call bas1%to_device()
      call trc_scf_run(bas1, nocc, nocc, opts, r1)
      call check(error, r1%converged, trim(xyz)//": the C1 SCF did not converge")
      if (allocated(error)) return

      allocate (gsym(3, nat), g1(3, nat))
      call lp%init(bas%nao)
      call trc_scf_gradient(bas, lp, rs%nspin, nocc, nocc, rs%cmo, rs%eps, 1.0e-12_dp, gsym, sym=trc_scf_grad_sym_t(geometry=geom))
      call trc_scf_gradient(bas1, lp, r1%nspin, nocc, nocc, r1%cmo, r1%eps, 1.0e-12_dp, g1)
      call lp%release()
      call sc%release()
      call bas%release()
      call bas1%release()

      dg = maxval(abs(gsym - g1))
      call equivariance(geom, gsym, asym, forbidden)
      print "(a,a,a,es10.2,a,es10.2,a,es10.2)", "  ", trim(xyz), ": |skeleton+projection - C1| ", dg, &
         "  equivariance ", asym, "  forbidden ", forbidden
      call check(error, dg <= TOL_C1, trim(xyz)//": the skeleton symmetric gradient differs from the C1 one by "// &
                 to_string(dg, "es10.2"))
      if (allocated(error)) return
      call check(error, asym <= TOL_EQV, trim(xyz)//": the projected gradient is not equivariant: "// &
                 to_string(asym, "es10.2"))
      if (allocated(error)) return
      call check(error, forbidden <= TOL_EQV, trim(xyz)//": forbidden gradient component "// &
                 to_string(forbidden, "es10.2"))
   end subroutine production_case

   subroutine declined_case(error, xyz, basis_json, spherical, group, expect_skeleton)
      !! The skeleton size rule governs the four-centre GRADIENT as well as
      !! the Fock build, which is what `trc_scf_gradient`'s `sym_skeleton`
      !! carries and what trc_capi.F90 sets from `trc_scf_result_t%sym_skeleton`.
      !!
      !! Two things are pinned, and both are properties of the ENGINE, not
      !! of the numbers -- the numbers agree on either path, which is exactly
      !! why the routing needs its own test.
      !!
      !! 1. A DECLINED skeleton reuses the SCF's engine. The SCF built that
      !!    engine WITHOUT the group (trc_scf_driver calls
      !!    scf_context%setup with no `sym_geom` when the rule declines, so
      !!    `orbits_built` is false), and a skeleton gradient cannot run on
      !!    it -- it would silently take the full plan -- so
      !!    trc_scf_gradient refuses to reuse it and builds its own with the
      !!    group. With the skeleton declined there is nothing to refuse:
      !!    `reused` must come back TRUE, and the whole second engine build
      !!    does not happen. On an ACCEPTED molecule the SCF's engine has
      !!    the orbits and `reused` is true as well, so the assertion is
      !!    `reused` on both, with `expect_skeleton` saying which reason.
      !!
      !! 2. The gradient is still projected. A declined skeleton takes the
      !!    full C1 four-centre term, so two_electron_4c does not project
      !!    and the end-of-routine sym_scf_grad_symmetrize is the only
      !!    projection left -- it has to still be there, or the symmetry
      !!    stops zeroing the forbidden components on exactly the small
      !!    molecules the rule protects.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, group
      logical, intent(in) :: spherical, expect_skeleton
      type(trc_basis_t) :: bas, bas1
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_sym_scf_t) :: sc
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: rs, r1
      type(trc_pair_eri_cache_t) :: cache
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :), gsym(:, :), g1(:, :)
      real(dp) :: dg, asym, forbidden
      integer :: nat, nocc, dev
      logical :: reused

      dev = trc_bind_device(0)
      call read_xyz(trim(xyz)//".xyz", nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, trim(group), "group of "//trim(xyz))
      if (allocated(error)) return
      call trc_basis_from_json(trim(basis_json), nat, elem, geom%r_std, bas, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = spherical
      call bas%to_device()
      nocc = sum(elem)/2
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100

      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      sc%ao = ao
      ! The SCF with its engine kept, exactly as trc_capi.F90 keeps `cx%peri`:
      ! the rule decides inside here, and `rs%sym_skeleton` is its verdict.
      call trc_scf_run(bas, nocc, nocc, opts, rs, sym=sc, pair_eri_cache=cache)
      call check(error, rs%converged, trim(xyz)//": the symmetric SCF did not converge")
      if (allocated(error)) return
      call check(error, rs%sym_skeleton .eqv. expect_skeleton, trim(xyz)//": the size rule's verdict is not the one "// &
                 "this case was written for (sym_atom_orbits = "//to_string(rs%sym_atom_orbits)//")")
      if (allocated(error)) return
      call check(error, cache%orbits_built .eqv. expect_skeleton, trim(xyz)//": the SCF's engine has the orbit "// &
                 "tables when and only when the skeleton build was taken")
      if (allocated(error)) return

      call trc_basis_from_json(trim(basis_json), nat, elem, geom%r_std, bas1, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas1%spherical = spherical
      call bas1%to_device()
      call trc_scf_run(bas1, nocc, nocc, opts, r1)
      call check(error, r1%converged, trim(xyz)//": the C1 SCF did not converge")
      if (allocated(error)) return

      allocate (gsym(3, nat), g1(3, nat))
      call lp%init(bas%nao)
      call trc_scf_gradient(bas, lp, rs%nspin, nocc, nocc, rs%cmo, rs%eps, 1.0e-12_dp, gsym, &
                            pair_eri_cache=cache, reused=reused, sym=trc_scf_grad_sym_t(geometry=geom, skeleton=rs%sym_skeleton))
      call trc_scf_gradient(bas1, lp, r1%nspin, nocc, nocc, r1%cmo, r1%eps, 1.0e-12_dp, g1)
      call lp%release()
      call cache%release()
      call sc%release()
      call bas%release()
      call bas1%release()

      dg = maxval(abs(gsym - g1))
      call equivariance(geom, gsym, asym, forbidden)
      print "(a,a,a,l1,a,l1,a,es10.2,a,es10.2)", "  ", trim(xyz), ": skeleton ", rs%sym_skeleton, &
         "  engine reused ", reused, "  |sym - C1| ", dg, "  forbidden ", forbidden
      call check(error, reused, trim(xyz)//": the gradient did not reuse the SCF's engine, so it rebuilt one")
      if (allocated(error)) return
      call check(error, dg <= TOL_C1, trim(xyz)//": the symmetric gradient differs from the C1 one by "// &
                 to_string(dg, "es10.2"))
      if (allocated(error)) return
      call check(error, asym <= TOL_EQV, trim(xyz)//": the gradient is not equivariant: "// &
                 to_string(asym, "es10.2"))
      if (allocated(error)) return
      call check(error, forbidden <= TOL_EQV, trim(xyz)//": forbidden gradient component "// &
                 to_string(forbidden, "es10.2"))
   end subroutine declined_case

   subroutine test_declined(error)
      !! Water is DECLINED by the skeleton size rule (3 atoms in 2 orbits,
      !! trc_sym_run's sym_run_skeleton_size_rule), HOCl too (3 in 3, and it is
      !! the planar Cs molecule whose skeleton plan is LARGER -- see
      !! `test_cs_planar`), and B2H6 is accepted. All three must give the
      !! same gradient; the first two must do it over the SCF's own engine.
      type(error_type), allocatable, intent(out) :: error
      call declined_case(error, "water", "basis_sets/6-31gs.json", .false., "C2v", .false.)
      if (allocated(error)) return
      call declined_case(error, "hocl", "basis_sets/6-31gs.json", .false., "Cs", .false.)
      if (allocated(error)) return
      call declined_case(error, "b2h6_d2h", "basis_sets/6-31gs.json", .false., "D2h", .true.)
   end subroutine test_declined

   subroutine test_production(error)
      type(error_type), allocatable, intent(out) :: error
      call production_case(error, "water", "basis_sets/6-31gs.json", .false., "C2v")
      if (allocated(error)) return
      call production_case(error, "water", "basis_sets/cc-pvdz.json", .true., "C2v")
      if (allocated(error)) return
      call production_case(error, "b2h6_d2h", "basis_sets/6-31gs.json", .false., "D2h")
   end subroutine test_production

   function scf_energy_at(basis_json, nat, elem, r, nocc, opts) result(e)
      !! The reconverged SCF energy at a displaced geometry, in C1: a
      !! displacement along one Cartesian component breaks the symmetry, so
      !! the finite difference cannot be taken inside the group.
      character(len=*), intent(in) :: basis_json
      integer, intent(in) :: nat
      integer, intent(in) :: elem(nat)
      real(dp), intent(in) :: r(3, nat)
      integer, intent(in) :: nocc
      type(trc_scf_options_t), intent(in) :: opts
      real(dp) :: e
      type(trc_basis_t) :: bas
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      e = 0.0_dp
      call trc_basis_from_json(basis_json, nat, elem, r, bas, err)
      if (err%has_error()) return
      call bas%to_device()
      call trc_scf_run(bas, nocc, nocc, opts, res)
      if (res%converged) e = res%energy
      call bas%release()
   end function scf_energy_at

   subroutine test_fd(error)
      !! The production symmetric (skeleton) gradient against a central
      !! difference of terco's own fully reconverged SCF energy -- the only
      !! check here that would catch an error consistent between the skeleton
      !! and the full four-centre path. Water/6-31G*, 18 displaced SCFs.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_sym_scf_t) :: sc
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: rs
      type(trc_linalg_t) :: lp
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r0(:, :), rt(:, :), grad(:, :)
      real(dp) :: ep, em, fd, worst_abs, worst_rel, denom
      integer :: nat, nocc, a, d, dev
      character(len=*), parameter :: BASFILE = "basis_sets/6-31gs.json"

      dev = trc_bind_device(0)
      call read_xyz("water.xyz", nat, elem, r0)
      geom = sym_geom_detect(r0, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, "C2v", "group of water")
      if (allocated(error)) return
      call trc_basis_from_json(BASFILE, nat, elem, geom%r_std, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      nocc = sum(elem)/2
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 200
      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      sc%ao = ao
      call trc_scf_run(bas, nocc, nocc, opts, rs, sym=sc)
      call check(error, rs%converged, "water: the symmetric SCF did not converge")
      if (allocated(error)) return
      allocate (grad(3, nat))
      call lp%init(bas%nao)
      call trc_scf_gradient(bas, lp, rs%nspin, nocc, nocc, rs%cmo, rs%eps, 1.0e-12_dp, grad, sym=trc_scf_grad_sym_t(geometry=geom))
      call lp%release()
      call sc%release()
      call bas%release()

      allocate (rt(3, nat))
      worst_abs = 0.0_dp
      worst_rel = 0.0_dp
      do a = 1, nat
         do d = 1, 3
            rt = geom%r_std
            rt(d, a) = geom%r_std(d, a) + FD_H
            ep = scf_energy_at(BASFILE, nat, elem, rt, nocc, opts)
            rt = geom%r_std
            rt(d, a) = geom%r_std(d, a) - FD_H
            em = scf_energy_at(BASFILE, nat, elem, rt, nocc, opts)
            call check(error, ep /= 0.0_dp .and. em /= 0.0_dp, "water: a displaced SCF did not converge")
            if (allocated(error)) return
            fd = (ep - em)/(2.0_dp*FD_H)
            denom = max(abs(grad(d, a)), abs(fd), 1.0e-10_dp)
            worst_abs = max(worst_abs, abs(grad(d, a) - fd))
            worst_rel = max(worst_rel, abs(grad(d, a) - fd)/denom)
         end do
      end do
      print "(a,es10.3,a,es10.3)", "  water C2v skeleton gradient vs FD: worst abs ", worst_abs, &
         "  worst rel ", worst_rel
      call check(error, worst_abs <= FD_ABSTOL .or. worst_rel <= FD_RELTOL, &
                 "water: the skeleton symmetric gradient disagrees with the finite difference: abs "// &
                 to_string(worst_abs, "es10.3")//", rel "//to_string(worst_rel, "es10.3"))
   end subroutine test_fd

end module check_sym_grad_suite

program check_sym_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use check_sym_grad_suite, only: collect_sym_grad, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if (trim(g_mode) == "f" .and. LMAX < 3) then
      print "(a,i0,a)", "  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ", LMAX, "."
      stop 0
   end if
   stat = 0
   ! The basis reader is not thread-safe, so the cases run serially.
   call run_testsuite(collect_sym_grad, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_grad
