!
! Parity skipping (docs/SYMMETRY_PARITY_SKIP.md): the Fock matrix built with
! the dead-axis masks and the engine frame equals the one built without them.
!
! `ctx%eri` is built with the point group and the dead axis of the molecule
! made the kernels' z (PARITY_FORCE_FRAME: the axis permutation is applied
! even before a masked kernel exists); `ctx%eri0` is the same engine with
! parity skipping off. For each molecule, Cartesian and spherical, with J, K
! and J+K, tiles on and off, an increment, and a random and a converged
! density, the symmetrized skeleton of the first and its full G are compared
! with the full G of the second against the S14.2 gate
!   |.|_max <= 5e-11 max(1, |G|_max).
! Counters (`eri%nwork_masked`, `nwork_variant`, `parity_axis`) show that the
! masks and, once a kernel exists, the variants really ran: a silent
! fall back to the unmasked plan fails on them. Without a point group
! nothing may be masked and the matrix must be bit-identical.
!
! Mode "f" is the uncontracted l = 0..3 set (skipped below LMAX = 3).
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_parity_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_bin_tiles, only: bin_tiles_set_enabled
   use trc_binkernel, only: binkernel_have_variants
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_parity, only: PARITY_OFF, PARITY_AUTO, PARITY_FORCE_FRAME, parity_set_default
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_test_basis, only: read_xyz, build_aux
   implicit none
   private
   public :: collect_sym_parity, g_mode

   character(len=8) :: g_mode = "main"

   real(dp), parameter :: THRESH = 1.0e-10_dp
   real(dp), parameter :: TOL_REL = 5.0e-11_dp
      !! The S14.2 gate, relative to max(1, |G|_max).

   type :: ctx_t
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

contains

   subroutine collect_sym_parity(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == "f") then
         testsuite = [new_unittest("uncontracted l = 0..3: water C2v, B2H6 D2h", test_lmax3)]
      else
         testsuite = [ &
                     new_unittest("planar D2h benzene", test_benzene), &
                     new_unittest("planar D2h naphthalene", test_naphthalene), &
                     new_unittest("planar D2h pentacene, banded density", test_pentacene), &
                     new_unittest("planar Cs HOCl", test_hocl), &
                     new_unittest("planar C2v water and C2h N2H2", test_c2), &
                     new_unittest("mixed masks: D2h B2H6", test_b2h6), &
                     new_unittest("no point group: nothing is masked, bit-identical", test_nosym), &
                     new_unittest("SCF energy with and without parity skipping", test_scf) &
                     ]
      end if
   end subroutine collect_sym_parity

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
      ctx%geom = sym_detect(r, ctx%elem, err, basis_id=ctx%elem)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
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

      ! The engine under test applies the frame even when no masked kernel
      ! exists (so the permutation round trip is tested on its own); eri0 is
      ! the same engine with parity skipping off.
      ctx%eri%parity_mode = PARITY_FORCE_FRAME
      call ctx%eri%build(ctx%bas, THRESH, sym=ctx%geom)
      ctx%eri0%parity_mode = PARITY_OFF
      call ctx%eri0%build(ctx%bas, THRESH, sym=ctx%geom)
      ctx%has0 = .true.
   end subroutine setup_ctx

   subroutine release_ctx(ctx)
      type(ctx_t), intent(inout) :: ctx
      call ctx%eri%release()
      if (ctx%has0) call ctx%eri0%release()
      call ctx%sc%release()
      call ctx%bas%release()
   end subroutine release_ctx

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
      ! Just above the accumulation-order noise floor (briefs, section 8).
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

   subroutine compare_build(ctx, dc, jscale, kscale, tag, tiles, error)
      !! The symmetrized skeleton and the full build of the engine under test
      !! against the full build of the engine without parity skipping.
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dc(ctx%nc, ctx%nc)
      real(dp), intent(in) :: jscale, kscale
      character(len=*), intent(in) :: tag
      logical, intent(in) :: tiles
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gr(:, :), gf(:, :), gsk(:, :), gref(:, :), gfull(:, :), gskel(:, :), w(:, :)
      real(dp) :: e_skel, e_full, gmax

      allocate (gr(ctx%nc, ctx%nc), gf(ctx%nc, ctx%nc), gsk(ctx%nc, ctx%nc))
      allocate (gref(ctx%ns, ctx%ns), gfull(ctx%ns, ctx%ns), gskel(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
      call bin_tiles_set_enabled(tiles)
      call ctx%eri0%fock(ctx%bas, dc, gr, k_scale=kscale, j_scale=jscale)
      call ctx%eri%fock(ctx%bas, dc, gf, k_scale=kscale, j_scale=jscale)
      call ctx%eri%fock(ctx%bas, dc, gsk, k_scale=kscale, j_scale=jscale, skeleton=.true.)
      call bin_tiles_set_enabled(.true.)
      call to_caller(ctx, gr, gref)
      call to_caller(ctx, gf, gfull)
      call to_caller(ctx, gsk, gskel)
      call ctx%sc%symmetrize(ctx%ns, gskel, w)
      gmax = maxval(abs(gref))
      e_skel = maxval(abs(gskel - gref))
      e_full = maxval(abs(gfull - gref))
      ctx%worst = max(ctx%worst, e_skel/max(1.0_dp, gmax), e_full/max(1.0_dp, gmax))
      print "(a,a,a,a,a,es10.3,a,es9.2,a,es9.2)", "  ", trim(ctx%label), " ", trim(tag), &
         ": max|G| ", gmax, "  skeleton err ", e_skel, "  full err ", e_full
      call check(error, e_skel <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)//" "//trim(tag)// &
                 ": symmetrized skeleton with parity skipping differs from the full G by "//to_string(e_skel, "es10.3"))
      if (allocated(error)) return
      call check(error, e_full <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)//" "//trim(tag)// &
                 ": full G with parity skipping differs by "//to_string(e_full, "es10.3"))
   end subroutine compare_build

   subroutine compare_incremental(ctx, dref, ddel, error)
      type(ctx_t), intent(inout) :: ctx
      real(dp), intent(in) :: dref(ctx%nc, ctx%nc), ddel(ctx%nc, ctx%nc)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ga(:, :), gb(:, :), gf(:, :), gs(:, :), w(:, :)
      real(dp) :: e, gmax
      allocate (ga(ctx%nc, ctx%nc), gb(ctx%nc, ctx%nc), gf(ctx%ns, ctx%ns), gs(ctx%ns, ctx%ns), w(ctx%ns, ctx%ns))
      call ctx%eri0%fock(ctx%bas, dref, ga)
      call ctx%eri0%fock(ctx%bas, ddel, gb)
      ga = ga + gb
      call to_caller(ctx, ga, gf)
      call ctx%eri%fock(ctx%bas, dref, ga, skeleton=.true.)
      call ctx%eri%fock(ctx%bas, ddel, gb, skeleton=.true.)
      ga = ga + gb
      call to_caller(ctx, ga, gs)
      call ctx%sc%symmetrize(ctx%ns, gs, w)
      gmax = maxval(abs(gf))
      e = maxval(abs(gs - gf))
      ctx%worst = max(ctx%worst, e/max(1.0_dp, gmax))
      print "(a,a,a,es10.3,a,es10.3)", "  ", trim(ctx%label), " incremental: max|G| ", gmax, "  err ", e
      call check(error, e <= TOL_REL*max(1.0_dp, gmax), trim(ctx%label)// &
                 ": incremental skeleton with parity skipping differs by "//to_string(e, "es10.3"))
   end subroutine compare_incremental

   subroutine check_counters(ctx, axis, min_masked, error)
      !! The masks were computed for the expected axis, the skeleton plan has
      !! at least `min_masked` of its work in masked segments, the engine
      !! without parity skipping has none, and (when this build has masked
      !! kernels) the variants ran for part of the masked work.
      type(ctx_t), intent(inout) :: ctx
      integer, intent(in) :: axis
      real(dp), intent(in) :: min_masked
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: dc(:, :), gc(:, :)
      real(dp) :: frac, ratio
      allocate (dc(ctx%nc, ctx%nc), gc(ctx%nc, ctx%nc))
      dc = 1.0e-3_dp
      call ctx%eri%fock(ctx%bas, dc, gc, skeleton=.true.)
      frac = real(ctx%eri%nwork_masked, dp)/real(max(ctx%eri%nwork, 1_int64), dp)
      print "(a,a,a,i0,a,f7.4,a,i0,a,i0)", "  ", trim(ctx%label), " dead axis ", ctx%eri%parity_axis, &
         "  masked share of the skeleton plan ", frac, "  variant launches ", ctx%eri%nlaunch_variant, &
         "  masked work in variants ", ctx%eri%nwork_variant
      call check(error, ctx%eri%parity_axis, axis, trim(ctx%label)//": dead axis")
      if (allocated(error)) return
      call check(error, frac >= min_masked, trim(ctx%label)//": masked share "//to_string(frac, "f7.4")// &
                 " is below "//to_string(min_masked, "f7.4")//" (the masks were dropped)")
      if (allocated(error)) return
      call ctx%eri0%fock(ctx%bas, dc, gc, skeleton=.true.)
      call check(error, ctx%eri0%nwork_masked == 0 .and. ctx%eri0%nlaunch_variant == 0 .and. ctx%eri0%parity_axis == 0, &
                 trim(ctx%label)//": the engine with parity skipping off masked something")
      if (allocated(error)) return
      if (binkernel_have_variants()) then
         call check(error, ctx%eri%nlaunch_variant > 0 .and. ctx%eri%nwork_variant > 0, &
                    trim(ctx%label)//": no launch went to a masked variant")
         if (allocated(error)) return
         ! The work the variants drop, modelled as the operations of the
         ! kernel text (VRR, accumulation, HRR) of the launches that ran in
         ! a variant: a variant that fell back to the full kernel, or one
         ! that dropped nothing, shows as a ratio near one. The z-masked
         ! kernels measure 0.4 to 0.7 of the unmasked ones (3.6 of the design).
         ratio = ctx%eri%parity_ops_masked/max(ctx%eri%parity_ops_base, 1.0_dp)
         print "(a,a,a,f7.4)", "  ", trim(ctx%label), " modelled operations of the variant launches, masked/unmasked: ", ratio
         call check(error, ratio < 0.8_dp, trim(ctx%label)//": the masked variants drop too little work: "// &
                    to_string(ratio, "f7.4"))
      else
         call check(error, ctx%eri%nlaunch_variant == 0, trim(ctx%label)//": a variant ran but none exists")
      end if
   end subroutine check_counters

   subroutine run_densities(ctx, scf, error)
      type(ctx_t), intent(inout) :: ctx
      logical, intent(in) :: scf
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: ds(:, :), dc(:, :), dd(:, :)
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
         call compare_build(ctx, dc, 1.0_dp, 1.0_dp, "J+K", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 1.0_dp, 0.0_dp, "J  ", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  tiles", .true., error)
         if (allocated(error)) return
         call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  no tiles", .false., error)
         if (allocated(error)) return
         call random_density(ctx, 11 + id, 0.02_dp, ds)
         call to_cart(ctx, ds, dd)
         call compare_incremental(ctx, dc, dd, error)
         if (allocated(error)) return
      end do
      ctx%label = base
   end subroutine run_densities

   subroutine run_molecule(xyz, group, basis_json, scf, axis, min_masked, error)
      character(len=*), intent(in) :: xyz, group, basis_json
      logical, intent(in) :: scf
      integer, intent(in) :: axis
      real(dp), intent(in) :: min_masked
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      integer :: ik, dev
      dev = trc_bind_device(0)
      do ik = 1, 2
         call setup_ctx(ctx, xyz, group, basis_json, -1, ik == 2, .true., error)
         if (allocated(error)) return
         call check_counters(ctx, axis, min_masked, error)
         if (allocated(error)) return
         call run_densities(ctx, scf, error)
         if (allocated(error)) return
         print "(a,a,a,es9.2)", "  ", trim(ctx%label), ": worst relative error ", ctx%worst
         call release_ctx(ctx)
      end do
   end subroutine run_molecule

   ! ------------------------------------------------------------------ tests

   subroutine test_benzene(error)
      !! Benzene is D6h, detected as D2h with the plane normal on x: the
      !! engine frame moves x to z. Every pair of atoms is in the plane.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("benzene", "D2h", "basis_sets/6-31gs.json", .false., 1, 0.95_dp, error)
      if (allocated(error)) return
      ! the same molecule in the other orientation of the test set: the normal on z
      call run_molecule("tops_benzene", "D2h", "basis_sets/6-31gs.json", .false., 3, 0.95_dp, error)
   end subroutine test_benzene

   subroutine test_naphthalene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("naphthalene_d2h", "D2h", "basis_sets/6-31gs.json", .false., 1, 0.95_dp, error)
   end subroutine test_naphthalene

   subroutine test_pentacene(error)
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: dc(:, :)
      integer :: dev
      dev = trc_bind_device(0)
      call setup_ctx(ctx, "pentacene_d2h", "D2h", "basis_sets/6-31g.json", -1, .false., .true., error)
      if (allocated(error)) return
      call check_counters(ctx, 1, 0.95_dp, error)
      if (allocated(error)) return
      allocate (dc(ctx%nc, ctx%nc))
      call banded_density(ctx, 5, 14.0_dp, dc)
      call compare_build(ctx, dc, 1.0_dp, 1.0_dp, "J+K", .true., error)
      if (allocated(error)) return
      call compare_build(ctx, dc, 0.0_dp, 1.0_dp, "K  tiles", .true., error)
      if (allocated(error)) return
      call release_ctx(ctx)
   end subroutine test_pentacene

   subroutine test_hocl(error)
      !! Planar Cs, plane normal already z: no permutation, masks only.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("hocl", "Cs", "basis_sets/cc-pvdz.json", .true., 3, 0.95_dp, error)
   end subroutine test_hocl

   subroutine test_c2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("water", "C2v", "basis_sets/cc-pvdz.json", .false., 1, 0.5_dp, error)
      if (allocated(error)) return
      call run_molecule("n2h2_trans", "C2h", "basis_sets/6-31gs.json", .true., 3, 0.5_dp, error)
   end subroutine test_c2

   subroutine test_b2h6(error)
      !! Six of the eight atoms lie in one mirror plane: part of the plan is masked (doc 3.8: 52 % of the
      !! full plan, 59 % of the skeleton), the rest runs unmasked in the same launch list.
      type(error_type), allocatable, intent(out) :: error
      call run_molecule("b2h6_d2h", "D2h", "basis_sets/6-31gs.json", .true., 1, 0.35_dp, error)
   end subroutine test_b2h6

   subroutine test_nosym(error)
      !! Without a point group there is nothing to mask, and the matrix is
      !! the same bits with parity skipping on and off.
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      type(trc_eri_t) :: ea, eb
      real(dp), allocatable :: ds(:, :), dc(:, :), ga(:, :), gb(:, :)
      integer :: dev
      dev = trc_bind_device(0)
      call setup_ctx(ctx, "tops_benzene", "D2h", "basis_sets/6-31gs.json", -1, .false., .false., error)
      if (allocated(error)) return
      allocate (ds(ctx%ns, ctx%ns), dc(ctx%nc, ctx%nc), ga(ctx%nc, ctx%nc), gb(ctx%nc, ctx%nc))
      call random_density(ctx, 5, 0.5_dp, ds)
      call to_cart(ctx, ds, dc)
      ea%parity_mode = PARITY_AUTO
      call ea%build(ctx%bas, THRESH)
      eb%parity_mode = PARITY_OFF
      call eb%build(ctx%bas, THRESH)
      call ea%fock(ctx%bas, dc, ga)
      call eb%fock(ctx%bas, dc, gb)
      call check(error, ea%parity_axis == 0 .and. ea%nwork_masked == 0 .and. ea%nlaunch_variant == 0, &
                 "a build without a point group masked something")
      if (allocated(error)) return
      call check(error, ea%nwork == eb%nwork .and. ea%nlaunch == eb%nlaunch, "the plan differs with parity skipping on")
      if (allocated(error)) return
      call check(error, all(ga == gb), "G is not bit-identical without a point group")
      call ea%release()
      call eb%release()
      call release_ctx(ctx)
   end subroutine test_nosym

   subroutine test_scf(error)
      !! The SCF (symmetric, skeleton Fock build) with the default parity
      !! skipping against the same run with it off: planar D2h, Cs, C2v, C2h.
      type(error_type), allocatable, intent(out) :: error
      character(len=16), parameter :: MOL(4) = [character(len=16) :: "tops_benzene", "hocl", "water", "n2h2_trans"]
      character(len=4), parameter :: GRP(4) = [character(len=4) :: "D2h", "Cs", "C2v", "C2h"]
      character(len=24), parameter :: BAS(4) = [character(len=24) :: "basis_sets/6-31gs.json", "basis_sets/cc-pvdz.json", &
                                                "basis_sets/cc-pvdz.json", "basis_sets/6-31gs.json"]
      type(ctx_t) :: ctx
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: ron, roff
      type(trc_sym_scf_t) :: sc
      integer :: im, nocc, dev
      dev = trc_bind_device(0)
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      do im = 1, 4
         call setup_ctx(ctx, trim(MOL(im)), trim(GRP(im)), trim(BAS(im)), -1, .false., .false., error)
         if (allocated(error)) return
         nocc = sum(ctx%elem)/2
         sc%ao = ctx%ao
         call parity_set_default(PARITY_AUTO)
         call trc_scf_run(ctx%bas, nocc, nocc, opts, ron, sym=sc)
         call sc%release()
         sc%ao = ctx%ao
         call parity_set_default(PARITY_OFF)
         call trc_scf_run(ctx%bas, nocc, nocc, opts, roff, sym=sc)
         call sc%release()
         call parity_set_default(-1)
         call check(error, ron%converged .and. roff%converged, trim(MOL(im))//": the SCF did not converge")
         if (allocated(error)) return
         print "(a,a,a,f18.10,a,es9.2)", "  ", trim(MOL(im)), " SCF energy ", ron%energy, "  parity on - off ", &
            ron%energy - roff%energy
         call check(error, abs(ron%energy - roff%energy) <= 1.0e-8_dp, trim(MOL(im))// &
                    ": SCF energy with parity skipping differs by "//to_string(ron%energy - roff%energy, "es10.3"))
         if (allocated(error)) return
         call release_ctx(ctx)
      end do
   end subroutine test_scf

   subroutine test_lmax3(error)
      type(error_type), allocatable, intent(out) :: error
      type(ctx_t) :: ctx
      real(dp), allocatable :: ds(:, :), dc(:, :)
      character(len=12), parameter :: MOL(2) = [character(len=12) :: "water", "b2h6_d2h"]
      character(len=4), parameter :: GRP(2) = [character(len=4) :: "C2v", "D2h"]
      integer :: im, ik, dev
      dev = trc_bind_device(0)
      do im = 1, 2
         do ik = 1, 2
            call setup_ctx(ctx, trim(MOL(im)), trim(GRP(im)), " ", 3, ik == 2, .true., error)
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
            deallocate (ds, dc)
            call release_ctx(ctx)
         end do
      end do
   end subroutine test_lmax3

end module check_sym_parity_suite

program check_sym_parity
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use check_sym_parity_suite, only: collect_sym_parity, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if (trim(g_mode) == "f" .and. LMAX < 3) then
      print "(a,i0,a)", "  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ", LMAX, "."
      stop 0
   end if
   stat = 0
   call run_testsuite(collect_sym_parity, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_parity
