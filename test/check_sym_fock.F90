!
! The skeleton Fock build: the S14 oracle.
!
! `eri%fock(..., skeleton=.true.)` launches only the orbit representatives of
! one side of each bin pair, weighted by their orbit size; the matrix it
! returns is G_skel, and (1/|G|) sum_g R G_skel R^T, the symmetrization the
! caller applies in its own AO space after the Cartesian -> spherical fold, is
! the full G. This check builds both on the SAME object (bins closed under the
! group, `eri%build(..., sym=geom)`; the full build walks the same pairs in
! another order) and compares them:
!
!   |sym(fold(G_skel)) - fold(G_full)|_max  <=  TOL_REL * max(1, |G_full|_max)
!
! the tolerance relative to the largest element as check_nosym's is: the sums
! have 1e5..1e7 terms and their order differs between the two builds. TOL_REL
! is justified from the measured errors the test prints (see the header of
! `TOL_REL`).
!
! Matrix, per molecule:
!   - Cartesian and spherical caller AO spaces (the symmetrizer is the
!     production one, `trc_sym_scf_t%symmetrize`, on the AO maps of the
!     caller's space; the skeleton is formed in the Cartesian working basis);
!   - J only, K only (the exchange-only launch, the one the exchange tiles
!     prune) and J+K;
!   - a random symmetrized density and a converged one (the symmetric SCF's);
!   - incremental: sym(G_skel(D_ref) + G_skel(dD)) against G(D_ref) + G(dD);
!   - the engine routing of the environment (TRC_RYS_CLASSES none / all runs
!     are separate ctest entries), TRC_BIN_CAP (the orbit-closed split) and
!     TRC_BIN_SORT (the orbit-aware sort).
! Molecules: one or more for every enabled group (Cs CHFCl2, Ci, C2 H2O2,
! C2v water, C2h N2H2, D2 twisted ethylene, D2h B2H6, naphthalene and
! benzene, which is D6h detected as D2h). Pentacene and hexacene (D2h, banded
! density) so that the exchange tiles cut: the K launch with tiles has fewer
! work items than without, in the full and in the skeleton plan alike.
!
! The quartet counts of the plan itself (`eri%nwork`), on a zero density so
! that no integral is evaluated, against the S13 prototype: cc-pVDZ at 1e-10
! gives 14627 / 9968 (water), 360415 / 107763 (B2H6), 4741001 / 1342393
! (benzene), 20671728 / 5451483 (naphthalene): the full count and the
! skeleton's (1.47, 3.34, 3.53, 3.79 times fewer). A silent fall back to the
! full build would show as a ratio of one.
!
! Registered as check_sym_fock (mode "main"), check_sym_fock_quick/_hgp/_rys/
! _cap/_sort (mode "quick" under the environment that names them) and
! check_sym_fock_f (mode "f": the uncontracted l = 0..3 set, which prints SKIP
! when this build's LMAX is below 3).
!
! The McMurchie-Davidson engines (S14.4): their bra is restricted to the orbit
! representatives too, on the engines' own counters (`eri%md%npq`, the
! (P, Q) primitive pairs of the Hermite stage; `eri%mg%nquart`, the
! shell-pair quartets of the exchange windows), which a silent fall back to
! the unrestricted engine would leave unchanged. check_sym_fock_md (mode
! "md", TRC_MDJ=1: the resident entry's J from MD J on the uncontracted
! l = 0..3 set; SKIP below LMAX 3), check_sym_fock_mdrule (mode "mdrule",
! TRC_MDJ unset: a system for which eri_mdj_size_rule picks MD J by itself) and
! check_sym_fock_mdg (mode "mdg": water cc-pVQZ, whose g shell sends the
! classes above TERCO_RYS_G_LTOT to MD exchange and J to MD J; SKIP below
! LMAX 4, where those classes are not generated).
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_fock_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_bin_tiles, only: bin_tiles_set_enabled
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_mdj, only: mdj_skeleton_bra_available
   use trc_parity, only: PARITY_OFF
   use trc_error, only: error_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_test_basis, only: read_xyz, build_aux
   implicit none
   private
   public :: collect_sym_fock, g_mode

   character(len=8) :: g_mode = "main"
      !! "main", "quick" or "f", set from argv before run_testsuite.

   real(dp), parameter :: THRESH = 1.0e-10_dp
   real(dp), parameter :: TOL_REL = 5.0e-11_dp
      !! Allowed |sym(G_skel) - G_full| relative to max(1, |G_full|_max).
      !! Measured (host, serial): 1e-15 for water, 1e-13 for B2H6 and the
      !! twisted ethylene, up to 7e-12 for benzene and the Ci molecule with
      !! its chlorines, whose random-density G has cancelling sums of 1e5..1e6
      !! terms (the S13 prototype's naphthalene was 5e-14 relative on a
      !! converged density). Under OpenMP the accumulation order is run
      !! dependent, so two builds of one matrix differ at this level by
      !! themselves. The 1e-11 the design suggested left under a factor of
      !! two of margin on the worst case; 5e-11 leaves seven, and a symmetry
      !! error is O(1) (the unsymmetrized skeleton is off by 1..80).
   real(dp), parameter :: TOL_REF = 1.0e-8_dp
      !! Allowed |G_full(sym-built bins) - G_full(bins without symmetry)|, relative to
      !! max(1, |G|_max): the symmetry-aware bins screen on the orbit maximum of
      !! Q and dsh, which keeps a handful of
      !! quartets the plain bins drop, each below the 1e-10 threshold.

   type :: ctx_t
      !! One molecule in one AO space, with its engine.
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_sym_scf_t) :: sc
      type(trc_eri_t) :: eri, eri0
      integer, allocatable :: elem(:)
      logical :: sph = .false.
      logical :: has0 = .false.
      integer :: nc = 0, ns = 0
      character(len=96) :: label = " "
      real(dp) :: worst = 0.0_dp
   end type ctx_t

   type :: mol_t
      character(len=20) :: xyz
      character(len=4) :: group
      character(len=32) :: basis
      logical :: scf
   end type mol_t

   type :: count_t
      character(len=20) :: xyz
      character(len=4) :: group
      integer(int64) :: nfull, nskel
   end type count_t

contains

   subroutine collect_sym_fock(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == "md") then
         testsuite = [ &
                     new_unittest("MD J forced, uncontracted l = 0..3: water C2v, B2H6 D2h", test_md_f) &
                     ]
      else if (g_mode == "mdrule") then
         testsuite = [ &
                     new_unittest("MD J by the size rule: naphthalene D2h, d-aug-cc-pVDZ", test_md_rule) &
                     ]
      else if (g_mode == "mdg") then
         testsuite = [ &
                     new_unittest("MD J and MD exchange: water C2v cc-pVQZ", test_mdg), &
                     new_unittest("MD J and MD exchange: uncontracted l = 0..4, B2H6 D2h", test_mdg_aux) &
                     ]
      else if (g_mode == "f") then
         testsuite = [ &
                     new_unittest("uncontracted l = 0..3, water C2v and B2H6 D2h: J, K, J+K, resident", test_lmax3) &
                     ]
      else if (g_mode == "quick") then
         testsuite = [ &
                     new_unittest("water C2v cc-pVDZ, B2H6 D2h 6-31G*: skeleton vs full", test_quick), &
                     new_unittest("quartet counts against the S13 prototype (water, B2H6)", test_counts_small) &
                     ]
      else
         testsuite = [ &
                     new_unittest("Cs CHFCl2", test_cs), &
                     new_unittest("C2 H2O2", test_c2), &
                     new_unittest("C2v water", test_c2v), &
                     new_unittest("C2h trans-N2H2", test_c2h), &
                     new_unittest("D2 twisted ethylene", test_d2), &
                     new_unittest("Ci meso-CHFCl-CHFCl", test_ci), &
                     new_unittest("D2h B2H6 (nonplanar)", test_b2h6), &
                     new_unittest("D2h naphthalene, 6-31G*", test_naphthalene), &
                     new_unittest("D6h benzene as D2h, 6-31G*", test_benzene), &
                     new_unittest("D2h pentacene and hexacene: the exchange tiles cut", test_acenes), &
                     new_unittest("quartet counts against the S13 prototype", test_counts) &
                     ]
      end if
   end subroutine collect_sym_fock

   ! ------------------------------------------------------------------ setup

   subroutine setup_ctx(ctx, xyz, group, basis_json, auxl, sph, with_ref, error)
      !! Detect the group, build the basis on its standard frame (from a JSON
      !! file, or the uncontracted l = 0..auxl set when auxl >= 0), the
      !! production symmetrizer of the caller's AO space, and the engines: one
      !! with the point group (the skeleton build) and, with `with_ref`, one
      !! without.
      type(ctx_t), intent(inout) :: ctx
      character(len=*), intent(in) :: xyz, group, basis_json
      integer, intent(in) :: auxl
      logical, intent(in) :: sph, with_ref
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
      ctx%has0 = with_ref
      if (with_ref) call ctx%eri0%build(ctx%bas, THRESH)
   end subroutine setup_ctx

   subroutine release_ctx(ctx)
      type(ctx_t), intent(inout) :: ctx
      call ctx%eri%release()
      if (ctx%has0) call ctx%eri0%release()
      call ctx%sc%release()
      call ctx%bas%release()
   end subroutine release_ctx

   ! -------------------------------------------------------------- densities

   subroutine random_density(ctx, seed, scale, ds)
      !! A symmetric matrix of pseudo-random numbers in the caller's AO space,
      !! symmetrized over the group. A linear congruential generator, so that
      !! the test does not depend on the compiler's random_number.
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

   subroutine banded_density(ctx, seed, reach, dc)
      !! A Cartesian density that vanishes between atoms farther apart than
      !! `reach` (bohr) and is invariant under the group (the distance is), so
      !! that the exchange tiles have blocks to drop. Cartesian caller space.
      type(ctx_t), intent(inout) :: ctx
      integer, intent(in) :: seed
      real(dp), intent(in) :: reach
      real(dp), intent(out) :: dc(ctx%nc, ctx%nc)
      real(dp), allocatable :: w(:, :)
      integer, allocatable :: at(:)
      integer :: s, a, i, j, n
      real(dp) :: dist
      n = ctx%nc
      allocate (at(n), w(n, n))
      do s = 1, ctx%bas%nshell
         do a = 0, (ctx%bas%sh_l(s) + 1)*(ctx%bas%sh_l(s) + 2)/2 - 1
            at(ctx%bas%sh_ao(s) + a) = ctx%bas%sh_at(s)
         end do
      end do
      call random_density(ctx, seed, 0.4_dp, dc)
      do j = 1, n
         do i = 1, n
            dist = norm2(ctx%bas%at_r(:, at(i)) - ctx%bas%at_r(:, at(j)))
            if (dist > reach) dc(i, j) = 0.0_dp
         end do
      end do
      call ctx%sc%symmetrize(n, dc, w)
   end subroutine banded_density

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
      type(ctx_t), intent(in) :: ctx
      real(dp), intent(in) :: ds(ctx%ns, ctx%ns)
      real(dp), intent(out) :: dc(ctx%nc, ctx%nc)
      if (ctx%sph) then
         call sph_sph_to_cart(ctx%bas, ds, dc)
      else
         dc = ds
      end if
   end subroutine to_cart

   subroutine to_caller(ctx, gc, gs)
      !! The fold of the Fock build: the Cartesian matrix in the caller's space.
      type(ctx_t), intent(in) :: ctx
      real(dp), intent(in) :: gc(ctx%nc, ctx%nc)
      real(dp), intent(out) :: gs(ctx%ns, ctx%ns)
      if (ctx%sph) then
         call sph_cart_to_sph(ctx%bas, gc, gs)
      else
         gs = gc
      end if
   end subroutine to_caller

   ! ------------------------------------------------------------ comparisons

   subroutine compare_build(ctx, dc, jscale, kscale, tag, tiles, error, nfull, nskel)
      !! The skeleton and the full build of the same operator on the same
      !! engine, the first folded and symmetrized, against the second folded.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      real(dp), intent(in) :: jscale, kscale
      character(len=*), intent(in) :: tag
      logical, intent(in) :: tiles
      type(error_type), allocatable, intent(out) :: error
      integer(int64), intent(out), optional :: nfull, nskel
      real(dp), allocatable :: gfc(:, :), gsc(:, :), gf(:, :), gs(:, :), w(:, :)
      real(dp) :: err_sym, err_raw, gmax
      integer(int64) :: nw_full, nw_skel

      allocate (gfc(ctx%nc, ctx%nc), gsc(ctx%nc, ctx%nc), gf(ctx%ns, ctx%ns), gs(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
      call bin_tiles_set_enabled(tiles)
      call ctx%eri%fock(ctx%bas, dc, gfc, opts=trc_fock_options_t(k_scale=kscale, j_scale=jscale))
      nw_full = ctx%eri%nwork
      call ctx%eri%fock(ctx%bas, dc, gsc, opts=trc_fock_options_t(k_scale=kscale, j_scale=jscale, skeleton=.true.))
      nw_skel = ctx%eri%nwork
      call bin_tiles_set_enabled(.true.)
      call to_caller(ctx, gfc, gf)
      call to_caller(ctx, gsc, gs)
      gmax = maxval(abs(gf))
      err_raw = maxval(abs(gs - gf))
      call ctx%sc%symmetrize(ctx%ns, gs, w)
      err_sym = maxval(abs(gs - gf))
      ctx%worst = max(ctx%worst, err_sym/max(1.0_dp, gmax))
      print "(a,a,a,a,a,es10.3,a,es10.3,a,es9.2,a,i0,a,i0)", "  ", trim(ctx%label), " ", trim(tag), &
         ": max|G| ", gmax, "  raw skeleton err ", err_raw, "  symmetrized err ", err_sym, &
         "  quartets ", nw_full, " -> ", nw_skel
      call check(error, err_sym <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)//" "//trim(tag)// &
                 ": symmetrized skeleton differs from the full G by "//to_string(err_sym, "es10.3")// &
                 " (max|G| "//to_string(gmax, "es10.3")//")")
      if (allocated(error)) return
      if (present(nfull)) nfull = nw_full
      if (present(nskel)) nskel = nw_skel
   end subroutine compare_build

   subroutine compare_nosym(ctx, dc, error)
      !! The full build on the symmetry-aware bins against the one on bins
      !! built without a point group: the same operator, screened on the orbit
      !! maximum of Q, and the walk in another order.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: g1(:, :), g0(:, :)
      real(dp) :: e, gmax
      if (.not. ctx%has0) return
      allocate (g1(ctx%nc, ctx%nc), g0(ctx%nc, ctx%nc))
      call ctx%eri%fock(ctx%bas, dc, g1)
      call ctx%eri0%fock(ctx%bas, dc, g0)
      gmax = maxval(abs(g0))
      e = maxval(abs(g1 - g0))
      print "(a,a,a,es10.3,a,es10.3)", "  ", trim(ctx%label), " full on symmetry-aware bins vs plain bins: max|G| ", gmax, &
         "  err ", e
      call check(error, e <= TOL_REF*max(1.0_dp, gmax), trim(ctx%label)// &
                 ": full build on the symmetry-aware bins differs from the plain bins by "//to_string(e, "es10.3"))
   end subroutine compare_nosym

   subroutine compare_incremental(ctx, dref, ddel, error)
      !! G(D_ref) + G(dD): the skeleton matrices accumulated and symmetrized
      !! once, against the sum of the full ones.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dref(ctx%nc, ctx%nc), ddel(ctx%nc, ctx%nc)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ga(:, :), gb(:, :), gf(:, :), gs(:, :), w(:, :)
      real(dp) :: e, gmax
      allocate (ga(ctx%nc, ctx%nc), gb(ctx%nc, ctx%nc), gf(ctx%ns, ctx%ns), gs(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
      call ctx%eri%fock(ctx%bas, dref, ga)
      call ctx%eri%fock(ctx%bas, ddel, gb)
      ga = ga + gb
      call to_caller(ctx, ga, gf)
      call ctx%eri%fock(ctx%bas, dref, ga, opts=trc_fock_options_t(skeleton=.true.))
      call ctx%eri%fock(ctx%bas, ddel, gb, opts=trc_fock_options_t(skeleton=.true.))
      ga = ga + gb
      call to_caller(ctx, ga, gs)
      call ctx%sc%symmetrize(ctx%ns, gs, w)
      gmax = maxval(abs(gf))
      e = maxval(abs(gs - gf))
      ctx%worst = max(ctx%worst, e/max(1.0_dp, gmax))
      print "(a,a,a,es10.3,a,es10.3)", "  ", trim(ctx%label), " incremental G(D_ref) + G(dD): max|G| ", gmax, "  err ", e
      call check(error, e <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)// &
                 ": incremental skeleton differs from the full G by "//to_string(e, "es10.3"))
   end subroutine compare_incremental

   subroutine compare_resident(ctx, dc, error)
      !! The resident entry (which the SCF uses, and which sends J to the MD J
      !! engine in an f basis, restricted like the rest) against the host full build.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gfc(:, :), gsc(:, :), gf(:, :), gs(:, :), w(:, :)
      real(dp) :: e, gmax
      allocate (gfc(ctx%nc, ctx%nc), gsc(ctx%nc, ctx%nc), gf(ctx%ns, ctx%ns), gs(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
      call ctx%eri%fock_resident(ctx%bas, dc, gfc)
      call to_caller(ctx, gfc, gf)
      call ctx%eri%fock_resident(ctx%bas, dc, gsc, opts=trc_fock_options_t(skeleton=.true.))
      call to_caller(ctx, gsc, gs)
      call ctx%sc%symmetrize(ctx%ns, gs, w)
      gmax = maxval(abs(gf))
      e = maxval(abs(gs - gf))
      ctx%worst = max(ctx%worst, e/max(1.0_dp, gmax))
      print "(a,a,a,es10.3,a,es10.3)", "  ", trim(ctx%label), " resident entry: max|G| ", gmax, "  err ", e
      call check(error, e <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)// &
                 ": resident skeleton differs from the full G by "//to_string(e, "es10.3"))
   end subroutine compare_resident

   subroutine compare_md(ctx, dc, jscale, kscale, omega, tag, resident, want_j, want_g, error)
      !! The skeleton and the full build through the entry that reaches the MD
      !! engines, which must have run (their counters are nonzero in the full
      !! build) and run less in the skeleton build (smaller counters).
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      real(dp), intent(in) :: jscale, kscale, omega
      character(len=*), intent(in) :: tag
      logical, intent(in) :: resident
      !! fock_resident, else fock
      logical, intent(in) :: want_j, want_g
      !! MD J / MD exchange must have run
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gfc(:, :), gsc(:, :), gf(:, :), gs(:, :), w(:, :)
      real(dp) :: e, gmax
      integer(int64) :: pq(2), qq(2)

      allocate (gfc(ctx%nc, ctx%nc), gsc(ctx%nc, ctx%nc), gf(ctx%ns, ctx%ns), gs(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
      ctx%eri%md%npq = 0
      ctx%eri%mg%nquart = 0
      if (resident) then
         call ctx%eri%fock_resident(ctx%bas, dc, gfc, &
                                    opts=trc_fock_options_t(k_scale=kscale, j_scale=jscale, omega=omega))
      else
         call ctx%eri%fock(ctx%bas, dc, gfc, opts=trc_fock_options_t(k_scale=kscale, j_scale=jscale, omega=omega))
      end if
      pq(1) = ctx%eri%md%npq
      qq(1) = ctx%eri%mg%nquart
      ctx%eri%md%npq = 0
      ctx%eri%mg%nquart = 0
      if (resident) then
         call ctx%eri%fock_resident(ctx%bas, dc, gsc, &
              opts=trc_fock_options_t(k_scale=kscale, j_scale=jscale, &
                        omega=omega, skeleton=.true.))
      else
         call ctx%eri%fock(ctx%bas, dc, gsc, &
                           opts=trc_fock_options_t(k_scale=kscale, j_scale=jscale, omega=omega, skeleton=.true.))
      end if
      pq(2) = ctx%eri%md%npq
      qq(2) = ctx%eri%mg%nquart
      call to_caller(ctx, gfc, gf)
      call to_caller(ctx, gsc, gs)
      gmax = maxval(abs(gf))
      call ctx%sc%symmetrize(ctx%ns, gs, w)
      e = maxval(abs(gs - gf))
      ctx%worst = max(ctx%worst, e/max(1.0_dp, gmax))
      print "(a,a,a,a,a,es10.3,a,es10.3,a,i0,a,i0,a,i0,a,i0)", "  ", trim(ctx%label), " ", trim(tag), &
         ": max|G| ", gmax, "  symmetrized err ", e, "  MD J (P,Q) pairs ", pq(1), " -> ", pq(2), &
         "  MD exchange quartets ", qq(1), " -> ", qq(2)
      call check(error, e <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)//" "//trim(tag)// &
                 ": symmetrized skeleton differs from the full G by "//to_string(e, "es10.3"))
      if (allocated(error)) return
      if (want_j) then
         call check(error, pq(1) > 0, trim(ctx%label)//" "//trim(tag)//": MD J did not run")
         if (allocated(error)) return
         if (mdj_skeleton_bra_available()) then
            call check(error, pq(2) > 0 .and. pq(2) < pq(1), trim(ctx%label)//" "//trim(tag)// &
                       ": the skeleton MD J did not restrict its bra")
         else
            ! Only the symmetric Hermite kernels in this build: the skeleton
            ! Fock build takes the full J (trc_mdj's fock), so the counts match.
            call check(error, pq(2) == pq(1), trim(ctx%label)//" "//trim(tag)// &
                       ": without gather kernels the skeleton build must take the full MD J")
         end if
         if (allocated(error)) return
      end if
      if (want_g) then
         call check(error, qq(1) > 0, trim(ctx%label)//" "//trim(tag)//": MD exchange did not run")
         if (allocated(error)) return
         call check(error, qq(2) > 0 .and. qq(2) < qq(1), trim(ctx%label)//" "//trim(tag)// &
                    ": the skeleton MD exchange did not restrict its bra")
         if (allocated(error)) return
      end if
   end subroutine compare_md

   subroutine run_densities(ctx, scf, error)
      !! The matrix of one molecule in one AO space: a random and (with
      !! `scf`) a converged density, J+K, J and K, the K launch with and
      !! without tiles, the increment.
      type(ctx_t), intent(inout) :: ctx
      logical, intent(in) :: scf
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ds(:, :), dc(:, :), dd(:, :)
      integer(int64) :: nf, ns
      integer :: id
      character(len=96) :: base

      allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc), dd(ctx%nc, ctx%nc))
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
         call compare_build(ctx, dc, 1.0_dp, 1.0_dp, "J+K", .true., error, nf, ns)
         if (allocated(error)) return
         ! the skeleton must launch fewer quartets: a fall back to the full plan fails here
         call check(error, ns < nf, trim(ctx%label)//": the skeleton plan is not smaller than the full one")
         if (allocated(error)) return
         call compare_build(ctx, dc, 1.0_dp, 0.0_dp, "J  ", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  tiles", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  no tiles", .false., error)
         if (allocated(error)) return
         if (id == 1) then
            call compare_nosym(ctx, dc, error)
            if (allocated(error)) return
            call random_density(ctx, 11, 0.02_dp, ds)
            call to_cart(ctx, ds, dd)
            call compare_incremental(ctx, dc, dd, error)
            if (allocated(error)) return
         else
            call random_density(ctx, 13, 0.01_dp, ds)
            call to_cart(ctx, ds, dd)
            call compare_incremental(ctx, dc, dd, error)
            if (allocated(error)) return
         end if
      end do
      ctx%label = base
   end subroutine run_densities

   subroutine run_molecule(xyz, group, basis_json, scf, error)
      !! One molecule, Cartesian and spherical.
      character(len=*), intent(in) :: xyz, group, basis_json
      logical, intent(in) :: scf
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      integer :: ik, dev
      dev = trc_bind_device(0)
      do ik = 1, 2
         call setup_ctx(ctx, xyz, group, basis_json, -1, ik == 2, ik == 1, error)
         if (allocated(error)) return
         call run_densities(ctx, scf, error)
         if (allocated(error)) return
         print "(a,a,a,es9.2)", "  ", trim(ctx%label), ": worst relative error ", ctx%worst
         call release_ctx(ctx)
      end do
   end subroutine run_molecule

   ! ------------------------------------------------------------------ tests

   subroutine test_cs(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("chfcl2", "Cs", "basis_sets/6-31gs.json", .true., error)
   end subroutine test_cs

   subroutine test_c2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("h2o2", "C2", "basis_sets/cc-pvdz.json", .true., error)
   end subroutine test_c2

   subroutine test_c2v(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("water", "C2v", "basis_sets/cc-pvdz.json", .true., error)
   end subroutine test_c2v

   subroutine test_c2h(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("n2h2_trans", "C2h", "basis_sets/6-31gs.json", .true., error)
   end subroutine test_c2h

   subroutine test_d2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("c2h4_twisted", "D2", "basis_sets/6-31gs.json", .true., error)
   end subroutine test_d2

   subroutine test_ci(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("chfcl2_ci", "Ci", "basis_sets/6-31gs.json", .false., error)
   end subroutine test_ci

   subroutine test_b2h6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("b2h6_d2h", "D2h", "basis_sets/6-31gs.json", .true., error)
   end subroutine test_b2h6

   subroutine test_naphthalene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("naphthalene_d2h", "D2h", "basis_sets/6-31gs.json", .false., error)
   end subroutine test_naphthalene

   subroutine test_benzene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("tops_benzene", "D2h", "basis_sets/6-31gs.json", .false., error)
   end subroutine test_benzene

   subroutine test_quick(error)
      !! The routing variants (TRC_RYS_CLASSES, TRC_BIN_CAP, TRC_BIN_SORT) run this.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("water", "C2v", "basis_sets/cc-pvdz.json", .false., error)
      if (allocated(error)) return
      call run_molecule("b2h6_d2h", "D2h", "basis_sets/6-31gs.json", .false., error)
   end subroutine test_quick

   subroutine test_acenes(error)
      !! Banded densities on the two long D2h molecules: the K launch with
      !! tiles has fewer work items than without, in the full plan and in the
      !! skeleton plan, and every result agrees with the full G.
      type(error_type), allocatable, intent(out) :: error
      character(len=14), parameter :: MOL(2) = [character(len=14) :: "pentacene_d2h", "hexacene_d2h"]
      type(ctx_t) :: ctx
      real(dp), allocatable :: dc(:, :), gfc(:, :), gsc(:, :), gf(:, :), gs(:, :), w(:, :)
      integer(int64) :: nw(4)
      integer :: im, dev
      real(dp) :: e, gmax

      dev = trc_bind_device(0)
      do im = 1, 2
         call setup_ctx(ctx, trim(MOL(im)), "D2h", "basis_sets/6-31g.json", -1, .false., .false., error)
         if (allocated(error)) return
         allocate (dc(ctx%nc, ctx%nc), gfc(ctx%nc, ctx%nc), gsc(ctx%nc, ctx%nc), gf(ctx%ns, ctx%ns), &
                   gs(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
         ! 7 bohr: ends of the molecule do not see each other
         call banded_density(ctx, 5, 14.0_dp, dc)
         call compare_build(ctx, dc, 1.0_dp, 1.0_dp, "J+K", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 1.0_dp, 0.0_dp, "J  ", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  tiles", .true., error, nw(1), nw(2))
         if (allocated(error)) return
         call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  no tiles", .false., error, nw(3), nw(4))
         if (allocated(error)) return
         print "(a,a,a,i0,a,i0,a,i0,a,i0)", "  ", trim(ctx%label), " K-only work items, full plan: tiles ", nw(1), &
            " vs ", nw(3), "; skeleton: tiles ", nw(2), " vs ", nw(4)
         call check(error, nw(1) < nw(3), trim(ctx%label)//": the exchange tiles dropped nothing from the full plan")
         if (allocated(error)) return
         call check(error, nw(2) < nw(4), trim(ctx%label)//": the exchange tiles dropped nothing from the skeleton plan")
         if (allocated(error)) return
         call check(error, nw(2) < nw(1), trim(ctx%label)//": the tiled skeleton plan is not smaller than the tiled full plan")
         if (allocated(error)) return
         ! the tiled skeleton against the untiled full build
         call bin_tiles_set_enabled(.false.)
         call ctx%eri%fock(ctx%bas, dc, gfc, opts=trc_fock_options_t(j_scale=0.0_dp))
         call bin_tiles_set_enabled(.true.)
         call ctx%eri%fock(ctx%bas, dc, gsc, opts=trc_fock_options_t(j_scale=0.0_dp, skeleton=.true.))
         call to_caller(ctx, gfc, gf)
         call to_caller(ctx, gsc, gs)
         call ctx%sc%symmetrize(ctx%ns, gs, w)
         gmax = maxval(abs(gf))
         e = maxval(abs(gs - gf))
         print "(a,a,a,es10.3,a,es10.3)", "  ", trim(ctx%label), " tiled skeleton vs untiled full K: max|G| ", gmax, "  err ", e
         call check(error, e <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)// &
                    ": tiled skeleton K differs from the untiled full K by "//to_string(e, "es10.3"))
         if (allocated(error)) return
         deallocate (dc, gfc, gsc, gf, gs, w)
         call release_ctx(ctx)
      end do
   end subroutine test_acenes

   subroutine count_case(xyz, group, basis_json, nfull, nskel, mixed, error)
      !! The plan's own counts on a zero density (every quartet fails the
      !! density screen, so nothing is evaluated).
      character(len=*), intent(in) :: xyz, group, basis_json
      integer(int64), intent(in) :: nfull, nskel
      logical, intent(in) :: mixed
      !! Part of the molecule's pairs lie in the dead plane of parity skipping and part do not
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      type(trc_eri_t) :: eoff
      real(dp), allocatable :: dc(:, :), gc(:, :)
      integer(int64) :: n1, n2, n1off, n2off
      integer :: dev, st
      character(len=16) :: cap
      dev = trc_bind_device(0)
      call setup_ctx(ctx, xyz, group, basis_json, -1, .false., .false., error)
      if (allocated(error)) return
      allocate (dc(ctx%nc, ctx%nc), gc(ctx%nc, ctx%nc))
      dc = 0.0_dp
      call ctx%eri%fock(ctx%bas, dc, gc)
      n1 = ctx%eri%nwork
      call ctx%eri%fock(ctx%bas, dc, gc, opts=trc_fock_options_t(skeleton=.true.))
      n2 = ctx%eri%nwork
      ! The S13 pins are the plan without parity skipping: the same engine with
      ! the masks off. Parity skipping bins the pairs of the dead plane apart,
      ! which can change which side of a bin pair the skeleton cuts.
      eoff%parity_mode = PARITY_OFF
      call eoff%build(ctx%bas, THRESH, sym=ctx%geom)
      call eoff%fock(ctx%bas, dc, gc)
      n1off = eoff%nwork
      call eoff%fock(ctx%bas, dc, gc, opts=trc_fock_options_t(skeleton=.true.))
      n2off = eoff%nwork
      call eoff%release()
      print "(a,a,a,i0,a,i0,a,f6.3,a,i0)", "  ", trim(xyz), " cc-pVDZ quartets: full ", n1, "  skeleton ", n2, "  (", &
         real(n1, dp)/real(n2, dp), " times fewer); skeleton without parity skipping ", n2off
      call check(error, n1off == nfull, trim(xyz)//": full plan quartets "//to_string(n1off)//", S13: "//to_string(nfull))
      if (allocated(error)) return
      ! Binning the dead-plane pairs apart never makes the skeleton plan larger:
      ! per bin pair it cuts the side with fewer quartets, rep(A)*|B| or |A|*rep(B),
      ! and for the sub-bins A = A1 + A2 and B = B1 + B2 the sum of the minima over
      ! the four sub-pairs is at most the minimum of the sums (a bin against itself:
      ! the triangle R*N becomes two triangles and the cheaper rectangle). The full plan is the
      ! same set of quartets however the bins are cut.
      call check(error, n1 == n1off, trim(xyz)//": the full plan changed with parity skipping: "// &
                 to_string(n1)//" against "//to_string(n1off))
      if (allocated(error)) return
      call check(error, n2 <= n2off, trim(xyz)//": parity skipping made the skeleton plan larger: "// &
                 to_string(n2)//" against "//to_string(n2off))
      if (allocated(error)) return
      ! A molecule with every pair in the plane (or none) has the same bins.
      if (mixed) then
         call check(error, n2 < n2off, trim(xyz)//": a mixed-mask molecule left the skeleton plan unchanged")
      else
         call check(error, n2 == n2off, trim(xyz)//": a uniform-mask molecule changed its skeleton plan")
      end if
      if (allocated(error)) return
      ! A lowered bin cap (TRC_BIN_CAP, the orbit-closed split) changes which
      ! side of a bin pair is cheaper to cut, so only the S13 count is exact.
      call get_environment_variable("TRC_BIN_CAP", cap, status=st)
      if (st == 0) then
         call check(error, n2 < n1, trim(xyz)//": the skeleton plan under a bin cap is not smaller than the full one")
      else
         call check(error, n2off == nskel, trim(xyz)//": skeleton plan quartets "//to_string(n2off)//", S13: "//to_string(nskel))
      end if
      if (allocated(error)) return
      call release_ctx(ctx)
   end subroutine count_case

   subroutine test_counts_small(error)
      type(error_type), allocatable, intent(out) :: error
      call count_case("water", "C2v", "basis_sets/cc-pvdz.json", 14627_int64, 9968_int64, .false., error)
      if (allocated(error)) return
      call count_case("b2h6_d2h", "D2h", "basis_sets/cc-pvdz.json", 360415_int64, 107763_int64, .true., error)
   end subroutine test_counts_small

   subroutine test_counts(error)
      type(error_type), allocatable, intent(out) :: error
      call test_counts_small(error)
      if (allocated(error)) return
      call count_case("tops_benzene", "D2h", "basis_sets/cc-pvdz.json", 4741001_int64, 1342393_int64, .false., error)
      if (allocated(error)) return
      call count_case("naphthalene_d2h", "D2h", "basis_sets/cc-pvdz.json", 20671728_int64, 5451483_int64, .false., error)
   end subroutine test_counts

   subroutine test_md_f(error)
      !! TRC_MDJ=1 (the ctest entry's environment): the resident entry forms J
      !! with MD J on the uncontracted l = 0..3 set. J+K (the bins' exchange
      !! skeleton next to the MD J skeleton), J alone and K alone.
      type(error_type), allocatable, intent(out) :: error
      character(len=12), parameter :: MOL(2) = [character(len=12) :: "water", "b2h6_d2h"]
      character(len=4), parameter :: GRP(2) = [character(len=4) :: "C2v", "D2h"]
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      character(len=8) :: e
      integer :: im, ik, dev
      dev = trc_bind_device(0)
      call get_environment_variable("TRC_MDJ", e)
      call check(error, trim(e) == "1", "check_sym_fock md mode needs TRC_MDJ=1 in the environment")
      if (allocated(error)) return
      do im = 1, 2
         do ik = 1, 2
            call setup_ctx(ctx, trim(MOL(im)), trim(GRP(im)), " ", 3, ik == 2, .false., error)
            if (allocated(error)) return
            allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc))
            call random_density(ctx, 3, 0.1_dp, ds)
            call to_cart(ctx, ds, dc)
            call compare_md(ctx, dc, 1.0_dp, 1.0_dp, 0.0_dp, "J+K", .true., .true., .false., error)
            if (allocated(error)) return
            call compare_md(ctx, dc, 1.0_dp, 0.0_dp, 0.0_dp, "J  ", .true., .true., .false., error)
            if (allocated(error)) return
            call compare_md(ctx, dc, 0.0_dp, 1.0_dp, 0.0_dp, "K  ", .true., .false., .false., error)
            if (allocated(error)) return
            deallocate (ds, dc)
            call release_ctx(ctx)
         end do
      end do
   end subroutine test_md_f

   subroutine test_md_rule(error)
      !! TRC_MDJ unset: eri_mdj_size_rule picks MD J for this extended system,
      !! with no f function and no environment to force it.
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      character(len=8) :: e
      integer :: dev, st
      dev = trc_bind_device(0)
      call get_environment_variable("TRC_MDJ", e, status=st)
      call check(error, st /= 0 .or. len_trim(e) == 0, "check_sym_fock mdrule mode needs TRC_MDJ unset")
      if (allocated(error)) return
      call setup_ctx(ctx, "naphthalene_d2h", "D2h", "basis_sets/d-aug-cc-pvdz.json", -1, .false., .false., error)
      if (allocated(error)) return
      call check(error, ctx%eri%mdj_extended, trim(ctx%label)//": the size rule does not choose MD J here")
      if (allocated(error)) return
      allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc))
      call random_density(ctx, 3, 0.1_dp, ds)
      call to_cart(ctx, ds, dc)
      call compare_md(ctx, dc, 1.0_dp, 1.0_dp, 0.0_dp, "J+K", .true., .true., .false., error)
      if (allocated(error)) return
      deallocate (ds, dc)
      call release_ctx(ctx)
   end subroutine test_md_rule

   subroutine test_mdg(error)
      !! Water cc-pVQZ: the g shell on O. `fock` sends J to MD J and the
      !! classes above TERCO_RYS_G_LTOT to MD exchange; with omega > 0 MD J
      !! stands aside and MD exchange forms J too.
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      integer :: ik, dev
      dev = trc_bind_device(0)
      do ik = 1, 2
         call setup_ctx(ctx, "water", "C2v", "basis_sets/cc-pvqz-pyscf17.json", -1, ik == 2, .false., error)
         if (allocated(error)) return
         allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc))
         call random_density(ctx, 3, 0.1_dp, ds)
         call to_cart(ctx, ds, dc)
         call compare_md(ctx, dc, 1.0_dp, 1.0_dp, 0.0_dp, "J+K", .false., .true., .true., error)
         if (allocated(error)) return
         call compare_md(ctx, dc, 0.0_dp, 1.0_dp, 0.0_dp, "K  ", .false., .false., .true., error)
         if (allocated(error)) return
         call compare_md(ctx, dc, 1.0_dp, 1.0_dp, 0.3_dp, "J+K omega 0.3", .false., .false., .true., error)
         if (allocated(error)) return
         call compare_md(ctx, dc, 1.0_dp, 1.0_dp, 0.0_dp, "J+K resident", .true., .true., .true., error)
         if (allocated(error)) return
         deallocate (ds, dc)
         call release_ctx(ctx)
      end do
   end subroutine test_mdg

   subroutine test_mdg_aux(error)
      !! The uncontracted l = 0..4 set (Cartesian: 8 minutes with the spherical
      !! run too) on nonplanar D2h B2H6, whose operations
      !! move atoms: the g classes of every pair of atoms, MD exchange and MD J
      !! in all the group's orbit sizes.
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      integer :: dev
      dev = trc_bind_device(0)
      call setup_ctx(ctx, "b2h6_d2h", "D2h", " ", 4, .false., .false., error)
      if (allocated(error)) return
      allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc))
      call random_density(ctx, 3, 0.05_dp, ds)
      call to_cart(ctx, ds, dc)
      call compare_md(ctx, dc, 1.0_dp, 1.0_dp, 0.0_dp, "J+K", .false., .true., .true., error)
      if (allocated(error)) return
      deallocate (ds, dc)
      call release_ctx(ctx)
   end subroutine test_mdg_aux

   subroutine test_lmax3(error)
      !! Uncontracted l = 0..3 shells on every atom: the f classes (Rys only),
      !! and MD J in the resident entry.
      type(error_type), allocatable, intent(out) :: error
      character(len=12), parameter :: MOL(2) = [character(len=12) :: "water", "b2h6_d2h"]
      character(len=4), parameter :: GRP(2) = [character(len=4) :: "C2v", "D2h"]
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      integer :: im, ik, dev
      dev = trc_bind_device(0)
      do im = 1, 2
         do ik = 1, 2
            call setup_ctx(ctx, trim(MOL(im)), trim(GRP(im)), " ", 3, ik == 2, .false., error)
            if (allocated(error)) return
            allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc))
            call random_density(ctx, 3, 0.1_dp, ds)
            call to_cart(ctx, ds, dc)
            call compare_build(ctx, dc, 1.0_dp, 1.0_dp, "J+K", .true., error)
            if (allocated(error)) return
            call compare_build(ctx, dc, 1.0_dp, 0.0_dp, "J  ", .true., error)
            if (allocated(error)) return
            call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  ", .true., error)
            if (allocated(error)) return
            call compare_resident(ctx, dc, error)
            if (allocated(error)) return
            deallocate (ds, dc)
            call release_ctx(ctx)
         end do
      end do
   end subroutine test_lmax3

end module check_sym_fock_suite

program check_sym_fock
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use check_sym_fock_suite, only: collect_sym_fock, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if ((trim(g_mode) == "f" .or. trim(g_mode) == "md") .and. LMAX < 3) then
      print "(a,i0,a)", "  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ", LMAX, "."
      stop 0
   end if
   if (trim(g_mode) == "mdg" .and. LMAX < 4) then
      print "(a,i0,a)", "  SKIP: MD exchange (g classes) needs LMAX = 4 and this build has LMAX = ", LMAX, "."
      stop 0
   end if
   stat = 0
   call run_testsuite(collect_sym_fock, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_fock
