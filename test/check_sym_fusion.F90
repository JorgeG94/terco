!
! Weighted-segment fusion in the bin Fock launch (W2 of
! .claude/context/segment-fusion-plan.md): ONE LAUNCH PER (class, mask).
!
! The skeleton plan's segments carry a weight -- the bra pairs' orbit size,
! halved on a diagonal bin pair, so a multiple of one half. The weight used to
! be folded into the launch scalars (`jf = launch%jfac*wt`), which forced
! `binkernel_fock_bins` to break its class run wherever the weight changed:
! one launch per (class, mask, WEIGHT). W2 made it a per-segment `real(dp)`
! array (`plan%sW`) that the kernels read at digestion and dropped the
! `plan%sW(c1 + 1) /= wt` exit, taking a steady-state Fock build from 85 to 33
! launches on benzene, 88 -> 36 naphthalene, 90 -> 36 pentacene and 218 -> 72
! on SF6.
!
! `check_sym_fock` is the oracle for the NUMBERS, and it pins `eri%nwork` --
! which quartets run, something fusion does not change. Nothing pinned the
! fusion. A change that reintroduced a per-weight split would keep every
! energy and every quartet count right and silently hand back the whole
! performance gain. This check is that pin, and it is structural rather than a
! literal 33/36/36/72: it rebuilds the plan the launcher built, counts the
! distinct (class, mask) groups in it, and asserts
!
!     eri%nlaunch == the number of distinct (class, mask) groups,
!
! i.e. that no group was split any further. The literals would be fragile
! against the bin split, the pair sort and the basis; the invariant is not,
! and it is what "fused" means.
!
! The counts come from the plan and not from an instrumented engine because
! `binkernel_fock_bins` releases its plan before it returns. The plan is
! reproducible from the engine's own state -- `build_skeleton`/`build` read
! the bins, the orbits, the threshold and `quartets%bin_pbase` and nothing
! else (the full plan's bin-level density screen is off on this path:
! `launch%use_density` is never set), so this file builds the same plan with
! the same two calls the launcher makes, in the same order, AFTER the Fock
! build, by which time the engine's one-time pair resort (TRC_BIN_SORT) has
! already settled the bins.
!
! Each case also asserts the invariant is not VACUOUS: the skeleton plan must
! hold more distinct (class, mask, weight) triples than (class, mask) groups,
! which is exactly the count the pre-W2 launcher would have reached. Without
! it a test that no group is split passes trivially on a plan with one weight
! per group.
!
! Covered per molecule: the skeleton plan with parity skipping as the
! environment leaves it (masked groups split by `sM`), the skeleton plan with
! parity skipping off (every mask 0, so a group is a bare class), and the
! FULL plan, which has no weights at all and whose launcher fills `plan%sW`
! with ones rather than branch in the kernels -- that path must be fused too.
!
! Molecules: water C2v, B2H6 D2h and benzene (D6h detected as D2h), cc-pVDZ,
! on a zero density so that no integral is evaluated (the plan is built and
! launched; every quartet then fails the kernels' density screen). Measured on
! this host build, skeleton groups against the triples the pre-W2 launcher
! would have launched: water 36 against 71, B2H6 63 against 168 (33 against
! 101 with the masks off), benzene 32 against 82. The whole file runs in under
! a second, because nothing is evaluated.
!
! `nlaunch` AND HIGHER LMAX: `launch_classes` does `nl = nl + 1` before the
! `if (mdg_routes(key))` branch that issues no kernel at all, so the counter
! counts a group whose class is routed to the MD exchange pass (trc_mdg). That
! makes `nlaunch` an overcount of the real kernel launches for any basis with
! g shells, and it would ALSO break the identity asserted here if the
! increment were ever moved after the branch. Both are kept out of the way by
! asserting that no class of these plans routes at all: cc-pVDZ for H, B, C
! and O reaches d, so no class key carries an l of 4 whatever TERCO_LMAX is,
! and the check on `nrouted` says so rather than leaving it to be assumed. A
! g basis (water cc-pVQZ, as check_sym_fock_mdg uses) must not be added here
! without deciding which of the two counts the assertion is about.
!
! TRC_LAUNCH_CAP re-splits a class run every `lcap` quartets. It is 0 by
! default and no test, script or CI job sets it (W0 audited this), so the
! cases below skip their assertions rather than model it if it is ever set.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_fusion_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_bins, only: trc_bin_plan_t
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_mdg, only: mdg_routes
   use trc_parity, only: PARITY_OFF
   use trc_pc_kernels, only: CLASS_RADIX
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sym_fusion

   real(dp), parameter :: THRESH = 1.0e-10_dp

contains

   subroutine collect_sym_fusion(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("one launch per (class, mask): water C2v", test_water), &
                  new_unittest("one launch per (class, mask): B2H6 D2h", test_b2h6), &
                  new_unittest("one launch per (class, mask): benzene D2h", test_benzene) &
                  ]
   end subroutine collect_sym_fusion

   subroutine plan_groups(eri, skel, nseg, ngroup, ntriple, nrouted)
      !! Rebuild the plan `binkernel_fock_bins` built for this engine, sort it
      !! the way the Fock launcher does, and count the distinct (class, mask)
      !! groups, the distinct (class, mask, weight) triples (the pre-W2 launch
      !! count) and the groups whose class is routed to MD exchange.
      type(trc_eri_t), intent(in) :: eri
      logical, intent(in) :: skel
      integer, intent(out) :: nseg, ngroup, ntriple, nrouted
      type(trc_bin_plan_t) :: plan
      integer :: s, key, mask, pkey, pmask
      real(dp) :: w, pw

      if (skel) then
         call plan%build_skeleton(eri%bins, eri%orb, eri%thresh, eri%quartets%bin_pbase)
      else
         call plan%build(eri%bins, eri%thresh, eri%quartets%bin_pbase, .false.)
      end if
      call plan%sort_by_class(CLASS_RADIX, by_mask=.true.)

      nseg = plan%nseg
      ngroup = 0
      ntriple = 0
      nrouted = 0
      pkey = -1
      pmask = -1
      pw = -1.0_dp
      do s = 1, plan%nseg
         key = plan%class_key(s, CLASS_RADIX)
         mask = plan%sM(s)
         w = 1.0_dp
         if (allocated(plan%sW)) w = plan%sW(s)
         ! The sorted key is ((class*2 + mask)*32 + nint(2*sW)), injective in
         ! all three (a mask is 0 or 1, and twice a realised weight is 1..8),
         ! so a change in any of them starts a new run and every group and
         ! every triple is one contiguous run. The weights are exact multiples
         ! of one half, which is why `w /= pw` is an exact comparison here and
         ! an integer `nint(2*sW)` inside the sort.
         if (key /= pkey .or. mask /= pmask) then
            ngroup = ngroup + 1
            ntriple = ntriple + 1
            if (mdg_routes(key)) nrouted = nrouted + 1
         else if (w /= pw) then
            ntriple = ntriple + 1
         end if
         pkey = key
         pmask = mask
         pw = w
      end do
      call plan%release()
   end subroutine plan_groups

   subroutine fusion_case(xyz, group, basis_json, error)
      !! One molecule: the skeleton launch is fused with parity skipping as
      !! configured and with it off, and so is the full launch.
      character(len=*), intent(in) :: xyz, group, basis_json
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_eri_t) :: eri, eoff
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), dc(:, :), gc(:, :)
      integer, allocatable :: elem(:)
      integer :: nat, dev, st, nseg, ngroup, ntriple, nrouted
      character(len=16) :: capenv

      ! TRC_LAUNCH_CAP splits a class run by quartet count, which this check
      ! does not model. Nothing in the tree sets it.
      capenv = " "
      call get_environment_variable("TRC_LAUNCH_CAP", capenv, status=st)
      if (st == 0) then
         if (len_trim(capenv) > 0) then
            print "(a,a,a)", "  ", trim(xyz), ": TRC_LAUNCH_CAP is set, launch fusion not asserted"
            return
         end if
      end if

      dev = trc_bind_device(0)
      call read_xyz(trim(xyz)//".xyz", nat, elem, r)
      geom = sym_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, trim(group), "group of "//trim(xyz))
      if (allocated(error)) return
      call trc_basis_from_json(trim(basis_json), nat, elem, geom%r_std, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = .false.
      call bas%to_device()
      allocate (dc(bas%nao, bas%nao), gc(bas%nao, bas%nao))
      dc = 0.0_dp

      call eri%build(bas, THRESH, sym=geom)
      ! The g classes have no HGP or Rys kernel and would make `nlaunch` count
      ! launches that are never issued; this basis has none, and `nrouted`
      ! below proves it rather than assuming it.
      call check(error, .not. eri%needs_md, trim(xyz)//": the basis needs the MD engines, see this file's header")
      if (allocated(error)) return

      ! Skeleton, parity skipping as configured.
      call eri%fock(bas, dc, gc, skeleton=.true.)
      call plan_groups(eri, .true., nseg, ngroup, ntriple, nrouted)
      print "(a,a,a,i0,a,i0,a,i0,a,i0)", "  ", trim(xyz), " skeleton: segments ", nseg, &
         ", (class, mask) groups ", ngroup, ", (class, mask, weight) triples ", ntriple, &
         ", launches ", eri%nlaunch
      call check(error, nrouted == 0, trim(xyz)//": "//to_string(nrouted)// &
                 " skeleton groups route to MD exchange, so nlaunch is not a group count")
      if (allocated(error)) return
      call check(error, eri%nlaunch == ngroup, trim(xyz)//": the skeleton build launched "// &
                 to_string(eri%nlaunch)//" times for "//to_string(ngroup)// &
                 " (class, mask) groups -- a group was split (the weight split would give "// &
                 to_string(ntriple)//")")
      if (allocated(error)) return
      ! Without this the assertion above would hold on a plan with one weight
      ! per group, where there is nothing to fuse.
      call check(error, ntriple > ngroup, trim(xyz)//": the skeleton plan has one weight per (class, mask) group ("// &
                 to_string(ntriple)//" triples, "//to_string(ngroup)//" groups), so it cannot show fusion")
      if (allocated(error)) return

      ! The full plan over the same bins: no weights at all, and the launcher
      ! hands the kernels an array of ones. One launch per (class, mask) too.
      call eri%fock(bas, dc, gc)
      call plan_groups(eri, .false., nseg, ngroup, ntriple, nrouted)
      print "(a,a,a,i0,a,i0,a,i0)", "  ", trim(xyz), " full:     segments ", nseg, &
         ", (class, mask) groups ", ngroup, ", launches ", eri%nlaunch
      call check(error, eri%nlaunch == ngroup, trim(xyz)//": the full build launched "// &
                 to_string(eri%nlaunch)//" times for "//to_string(ngroup)//" (class, mask) groups")
      if (allocated(error)) return

      ! Parity skipping off: every segment mask is 0, so a group is a bare
      ! class and the weight is the only thing left that could split a run.
      eoff%parity_mode = PARITY_OFF
      call eoff%build(bas, THRESH, sym=geom)
      call eoff%fock(bas, dc, gc, skeleton=.true.)
      call plan_groups(eoff, .true., nseg, ngroup, ntriple, nrouted)
      print "(a,a,a,i0,a,i0,a,i0,a,i0)", "  ", trim(xyz), " skeleton, no parity masks: segments ", nseg, &
         ", classes ", ngroup, ", (class, weight) pairs ", ntriple, ", launches ", eoff%nlaunch
      call check(error, nrouted == 0, trim(xyz)//": "//to_string(nrouted)//" unmasked groups route to MD exchange")
      if (allocated(error)) return
      call check(error, eoff%nlaunch == ngroup, trim(xyz)//": the unmasked skeleton build launched "// &
                 to_string(eoff%nlaunch)//" times for "//to_string(ngroup)//" classes (the weight split would give "// &
                 to_string(ntriple)//")")
      if (allocated(error)) return
      call check(error, ntriple > ngroup, trim(xyz)//": the unmasked skeleton plan has one weight per class, "// &
                 "so it cannot show fusion")
      if (allocated(error)) return

      call eoff%release()
      call eri%release()
      call bas%release()
      deallocate (dc, gc, r, elem)
   end subroutine fusion_case

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      call fusion_case("water", "C2v", "basis_sets/cc-pvdz.json", error)
   end subroutine test_water

   subroutine test_b2h6(error)
      type(error_type), allocatable, intent(out) :: error
      call fusion_case("b2h6_d2h", "D2h", "basis_sets/cc-pvdz.json", error)
   end subroutine test_b2h6

   subroutine test_benzene(error)
      type(error_type), allocatable, intent(out) :: error
      call fusion_case("tops_benzene", "D2h", "basis_sets/cc-pvdz.json", error)
   end subroutine test_benzene

end module check_sym_fusion_suite

program check_sym_fusion
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_fusion_suite, only: collect_sym_fusion
   implicit none
   integer :: stat
   stat = 0
   ! Serial: test-drive's default runs the cases under an OpenMP parallel do
   ! and the basis JSON reader is not thread-safe.
   call run_testsuite(collect_sym_fusion, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_fusion
