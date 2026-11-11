!
! The skeleton Fock build's SIZE RULE (trc_sym_run's sym_run_skeleton_size_rule,
! read through sym_run_skeleton_mode and applied by trc_scf_run): which symmetric
! runs build the Fock matrix from the orbit representatives and which build it
! in full and symmetrize it afterwards.
!
! WHAT IS PINNED. The rule is "at least two atoms per symmetry-inequivalent
! atom", so the counts it reads are part of it: this test pins the atom and
! orbit counts of the four systems its comment was measured on (tops_sf6,
! benzene, naphthalene_d2h, pentacene_d2h) and of the two other D2h molecules
! of the suite, together with the verdict for each. Detection alone, no basis
! and no integrals, so it is cheap.
!
! AND WHAT THE RULE DECIDES FOR EVERY OTHER SYMMETRY ORACLE. The rule also
! silently settles which Fock build each per-group SCF oracle exercises, so
! `test_oracle_molecules` pins its verdict for every molecule those oracles
! run, group by group. See that routine's comment.
!
! AND THAT DECLINING CHANGES NOTHING BUT THE WORK. Three molecules are run
! with the rule in charge, with the skeleton build forced on and with it
! forced off: HOCl (declined, 3 atoms in 3 orbits -- a planar Cs molecule has
! no atom orbit bigger than one), B2H6 (accepted, 8 atoms in 3 orbits) and
! SF6 (declined, 7 atoms in 4 orbits, the rule's own measured loss and the
! only case here with degenerate levels). The energy, the orbital energies,
! the density, the irrep decomposition, the dropped-function count and the
! iteration count must agree across the three.
!
! WHAT "THE SAME IRREP LABELS" MEANS, AND WHY IT IS NOT PER-ORBITAL. The two
! Fock builds accumulate G in a different order, so an orbital energy moves in
! its last bits, and `trc_sym_scf`'s energy merge then breaks a near-tie
! between degenerate partners the other way round. Measured here on SF6
! (cc-pVDZ, D2h, 102 orbitals): 18 of 102 entries move, as 9 TRANSPOSITIONS,
! each one between two orbitals whose energies agree to ~1e-13 Eh. (SF6 in D2h
! has twofold and threefold degenerate levels, and a transposition inside a
! threefold level leaves its third member alone, so the count of transpositions
! is not a count of levels.) That is pre-existing between the two paths and is
! not a defect -- src/api/trc_capi.F90's 1.25 note documents it -- so a
! per-orbital label comparison is not an invariant of the Fock path and must
! not be asserted as one.
!
! What IS invariant is the irrep DECOMPOSITION: the multiset of labels
! inside each degenerate block, which `irrep_blocks_agree` compares. On a
! non-degenerate level that block has one member and the check reduces to the
! per-orbital label, so a genuine relabelling is still caught wherever the
! label is determined; summed over blocks it is the per-irrep orbital count.
!
! THE ENVIRONMENT OVERRIDE is tested by the ctest entries, which run this
! program three times: without TRC_SYM_SKELETON, with 1 and with 0. The
! expectations below follow whatever the environment says, so the override is
! checked end to end, through the SCF and not only through the mode function.
!
! WHAT EACH ctest VARIANT CAN FAIL. Only `_auto` compares the two Fock builds:
! there the forced runs take opposite paths and the agreement above is a real
! two-path comparison. Under `_on` and `_off` the environment forces all three
! runs onto one path, so their agreement is near-vacuous and is NOT what those
! entries are for. What they pin, and can fail, is the override plumbing:
! TRC_SYM_SKELETON must reach every run whatever that run requested, must
! override a request to the contrary, and must leave the size rule unconsulted
! (`sym_atom_orbits` = 0). The energy agreement across the three runs is also
! not asserted bit-for-bit anywhere: see `run_three`.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_skeleton_rule_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_sorting, only: sort
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect, sym_geom_atom_orbits
   use trc_sym_run, only: sym_run_skeleton_mode, sym_run_skeleton_size_rule, &
                          SYM_SKEL_AUTO, SYM_SKEL_ON, SYM_SKEL_OFF, SYM_SKEL_MIN_ORBIT
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sym_skeleton_rule

   real(dp), parameter :: TOL_E = 1.0e-10_dp   !! the project's energy gate, skeleton against full
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp  !! orbital energies
   real(dp), parameter :: TOL_D = 1.0e-7_dp    !! density

contains

   subroutine collect_sym_skeleton_rule(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("the measured systems' atom orbits and verdicts", test_counts), &
                  new_unittest("every per-group SCF oracle's molecules, and what the rule does with them", &
                               test_oracle_molecules), &
                  new_unittest("the threshold itself, at and below the boundary", test_boundary), &
                  new_unittest("TRC_SYM_SKELETON overrides the request", test_env), &
                  new_unittest("HOCl cc-pVDZ (declined): the same SCF three ways", test_declined), &
                  new_unittest("B2H6 6-31G* (accepted): the same SCF three ways", test_accepted), &
                  new_unittest("SF6 cc-pVDZ (declined, degenerate levels): the same SCF three ways", test_sf6) &
                  ]
   end subroutine collect_sym_skeleton_rule

   subroutine one_count(error, xyz, group, natm, norbit, verdict)
      !! Detect `xyz`'s group, count its atom orbits and check the rule's
      !! verdict against what sym_run_skeleton_size_rule's table says.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, group
      integer, intent(in) :: natm, norbit
      logical, intent(in) :: verdict
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer :: nat, norb
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)

      call read_xyz(xyz, nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), xyz//": sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, group, xyz//": group")
      if (allocated(error)) return
      call check(error, nat, natm, xyz//": atoms")
      if (allocated(error)) return
      norb = sym_geom_atom_orbits(geom%perm)
      call check(error, norb, norbit, xyz//": atoms inequivalent under "//group)
      if (allocated(error)) return
      call check(error, sym_run_skeleton_size_rule(nat, norb) .eqv. verdict, &
                 xyz//": the size rule wants the "//trim(merge("skeleton", "full    ", verdict))// &
                 " build at "//to_string(nat)//" atoms in "//to_string(norb)//" orbits")
   end subroutine one_count

   subroutine test_counts(error)
      !! The rows of sym_run_skeleton_size_rule's measured table, plus the two
      !! other D2h molecules of the suite. SF6 is the only measured loss: its
      !! sulfur is fixed by all eight operations and each fluorine by four, so
      !! its orbits are small although it carries 102 cc-pVDZ functions.
      type(error_type), allocatable, intent(out) :: error
      call one_count(error, "tops_sf6.xyz", "D2h", 7, 4, .false.)
      if (allocated(error)) return
      call one_count(error, "benzene.xyz", "D2h", 12, 4, .true.)
      if (allocated(error)) return
      call one_count(error, "naphthalene_d2h.xyz", "D2h", 18, 5, .true.)
      if (allocated(error)) return
      call one_count(error, "pentacene_d2h.xyz", "D2h", 36, 10, .true.)
      if (allocated(error)) return
      call one_count(error, "b2h6_d2h.xyz", "D2h", 8, 3, .true.)
      if (allocated(error)) return
      call one_count(error, "c2h4_d2h.xyz", "D2h", 6, 2, .true.)
      if (allocated(error)) return
      ! A planar Cs molecule: the mirror is the molecular plane, so every atom
      ! is its own image and no shell pair has a nontrivial orbit. The
      ! skeleton build would launch every quartet at weight one and save
      ! nothing, which is what the rule says.
      call one_count(error, "hocl.xyz", "Cs", 3, 3, .false.)
   end subroutine test_counts

   subroutine test_oracle_molecules(error)
      !! One row per molecule that a per-group SCF oracle runs, with the
      !! rule's verdict. This exists because the rule QUIETLY DECIDES WHAT
      !! THOSE ORACLES COVER: a molecule it accepts is run on the skeleton
      !! build by default, one it declines on the full build, and nothing in
      !! the oracle's own output says which. Two whole groups come out
      !! declined -- every C2v molecule of check_sym_c2v and every Cs molecule
      !! of check_sym_scf_cs -- so those two files now force the skeleton
      !! build on a case each and assert the path (see their headers). The
      !! other five groups keep it by this verdict alone, which is why the
      !! verdict is pinned here: a change to the rule or to a geometry that
      !! took the skeleton build away from Ci, C2, C2h, D2 or D2h would turn
      !! this test red instead of silently narrowing what they test.
      type(error_type), allocatable, intent(out) :: error
      ! Ci, C2, C2h, D2: accepted, so their oracles exercise the skeleton
      ! build on their default path.
      call one_count(error, "chfcl2_ci.xyz", "Ci", 8, 4, .true.)
      if (allocated(error)) return
      call one_count(error, "n2h2f2_ci.xyz", "Ci", 6, 3, .true.)
      if (allocated(error)) return
      call one_count(error, "h2o2.xyz", "C2", 4, 2, .true.)
      if (allocated(error)) return
      call one_count(error, "n2h4.xyz", "C2", 6, 3, .true.)
      if (allocated(error)) return
      call one_count(error, "n2h2_trans.xyz", "C2h", 4, 2, .true.)
      if (allocated(error)) return
      call one_count(error, "glyoxal_trans.xyz", "C2h", 6, 3, .true.)
      if (allocated(error)) return
      call one_count(error, "clh2c_ch2cl_anti.xyz", "C2h", 8, 3, .true.)
      if (allocated(error)) return
      call one_count(error, "c2h4_twisted.xyz", "D2", 6, 2, .true.)
      if (allocated(error)) return
      call one_count(error, "c2f4_twisted.xyz", "D2", 6, 2, .true.)
      if (allocated(error)) return
      ! C2v and Cs: every molecule declined. The C2 axis and the mirrors fix
      ! most of their atoms, and a planar Cs molecule fixes all of them.
      call one_count(error, "water.xyz", "C2v", 3, 2, .false.)
      if (allocated(error)) return
      call one_count(error, "nh2.xyz", "C2v", 3, 2, .false.)
      if (allocated(error)) return
      call one_count(error, "formaldehyde.xyz", "C2v", 4, 3, .false.)
      if (allocated(error)) return
      call one_count(error, "ch2cl2.xyz", "C2v", 5, 3, .false.)
      if (allocated(error)) return
      call one_count(error, "hco.xyz", "Cs", 3, 3, .false.)
      if (allocated(error)) return
      call one_count(error, "chfcl2.xyz", "Cs", 5, 4, .false.)
   end subroutine test_oracle_molecules

   subroutine test_boundary(error)
      !! The threshold is a mean orbit size of SYM_SKEL_MIN_ORBIT, exactly.
      type(error_type), allocatable, intent(out) :: error
      integer :: n
      n = SYM_SKEL_MIN_ORBIT
      call check(error, sym_run_skeleton_size_rule(5*n, 5), "5 orbits of exactly the threshold size are accepted")
      if (allocated(error)) return
      call check(error, .not. sym_run_skeleton_size_rule(5*n - 1, 5), "one atom below the threshold is declined")
      if (allocated(error)) return
      call check(error, sym_run_skeleton_size_rule(5*n + 1, 5), "one atom above the threshold is accepted")
      if (allocated(error)) return
      ! A count the detection never produces must not be read as a verdict.
      call check(error, .not. sym_run_skeleton_size_rule(0, 0), "an empty molecule is declined, not accepted")
   end subroutine test_boundary

   subroutine test_env(error)
      !! sym_run_skeleton_mode returns the request unless TRC_SYM_SKELETON names a
      !! mode, and then that mode whatever was requested.
      type(error_type), allocatable, intent(out) :: error
      integer :: forced
      character(len=16) :: e
      e = " "
      call get_environment_variable("TRC_SYM_SKELETON", e)
      forced = sym_run_skeleton_mode(SYM_SKEL_AUTO)
      select case (trim(e))
      case ("1")
         call check(error, forced, SYM_SKEL_ON, "TRC_SYM_SKELETON=1 must force the skeleton build on")
         if (allocated(error)) return
         call check(error, sym_run_skeleton_mode(SYM_SKEL_OFF), SYM_SKEL_ON, "and must override a request to decline")
      case ("0")
         call check(error, forced, SYM_SKEL_OFF, "TRC_SYM_SKELETON=0 must force the skeleton build off")
         if (allocated(error)) return
         call check(error, sym_run_skeleton_mode(SYM_SKEL_ON), SYM_SKEL_OFF, "and must override a request to use it")
      case default
         call check(error, forced, SYM_SKEL_AUTO, "without TRC_SYM_SKELETON the request stands")
         if (allocated(error)) return
         call check(error, sym_run_skeleton_mode(SYM_SKEL_ON), SYM_SKEL_ON, "and a forced request is passed through")
      end select
   end subroutine test_env

   logical function irrep_blocks_agree(eps_ref, irr_ref, irr_other, tol) result(ok)
      !! Whether two runs assign the same irreps to the same orbitals, read as
      !! an invariant that survives degeneracy: the orbitals are cut into
      !! blocks of levels no further apart than `tol`, and the two runs must
      !! carry the same MULTISET of labels in every block. A block of one is
      !! the per-orbital label, so nothing is weakened where the label is
      !! determined; a block of several is where it is not, because the Fock
      !! path decides which partner the energy sort puts first (see this
      !! file's header). `eps_ref` must be ascending, as `trc_scf_result_t%eps`
      !! is.
      real(dp), intent(in) :: eps_ref(:)
      integer, intent(in) :: irr_ref(:), irr_other(:)
      real(dp), intent(in) :: tol
      integer :: i, j, n
      integer, allocatable :: a(:), b(:)

      ok = .true.
      n = size(eps_ref)
      i = 1
      do while (i <= n)
         j = i
         do while (j < n)
            if (eps_ref(j + 1) - eps_ref(j) > tol) exit
            j = j + 1
         end do
         a = irr_ref(i:j)
         b = irr_other(i:j)
         call sort(a)
         call sort(b)
         if (any(a /= b)) then
            ok = .false.
            return
         end if
         i = j + 1
      end do
   end function irrep_blocks_agree

   subroutine run_three(error, xyz, basis_json, label, rule_wants)
      !! `xyz` run with the rule in charge, with the skeleton build forced on
      !! and with it forced off. Checks which path each run took, then that
      !! the three agree on everything the user is told.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, label
      logical, intent(in) :: rule_wants
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_sym_scf_t) :: sca, scon, scoff
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: ra, ron, roff
      type(error_t) :: err
      integer :: nat, nelec, nalpha, nbeta, dev, norb, nmo, sp
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      logical :: want_auto
      real(dp) :: de, deps, dd

      dev = trc_bind_device(0)
      call read_xyz(xyz, nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), xyz//": sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call trc_basis_from_json(basis_json, nat, elem, geom%r_std, bas, err)
      call check(error, .not. err%has_error(), xyz//": basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), xyz//": sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      norb = sym_geom_atom_orbits(geom%perm)

      nelec = sum(elem)
      nalpha = nelec/2
      nbeta = nelec - nalpha
      opts%conv_energy = 1.0e-10_dp
      opts%max_iteration = 100

      ! One object per run: `prepare` builds per-irrep workspace into it, and
      ! the three runs must not share any of it.
      sca%ao = ao
      sca%geom = geom
      sca%skeleton = SYM_SKEL_AUTO
      call trc_scf_run(bas, nalpha, nbeta, opts, ra, sym=sca)
      call sca%release()
      scon%ao = ao
      scon%geom = geom
      scon%skeleton = SYM_SKEL_ON
      call trc_scf_run(bas, nalpha, nbeta, opts, ron, sym=scon)
      call scon%release()
      scoff%ao = ao
      scoff%geom = geom
      scoff%skeleton = SYM_SKEL_OFF
      call trc_scf_run(bas, nalpha, nbeta, opts, roff, sym=scoff)
      call scoff%release()
      call bas%release()

      call check(error, ra%converged .and. ron%converged .and. roff%converged, &
                 label//": an SCF did not converge: "//trim(ra%message)//trim(ron%message)//trim(roff%message))
      if (allocated(error)) return

      ! Which path each run took. TRC_SYM_SKELETON, when set, decides all
      ! three; otherwise the forced runs are forced and the rule decides the
      ! first one.
      select case (sym_run_skeleton_mode(SYM_SKEL_AUTO))
      case (SYM_SKEL_ON)
         want_auto = .true.
         call check(error, ron%sym_skeleton .and. roff%sym_skeleton, label//": TRC_SYM_SKELETON=1 must reach every run")
         if (allocated(error)) return
         ! The environment settled the mode, so the rule was never asked, and
         ! `sym_atom_orbits` must say so. This is the assertion the forced
         ! variants carry that they can actually fail.
         call check(error, ra%sym_atom_orbits, 0, &
                    label//": TRC_SYM_SKELETON=1 forced the mode, so the size rule must not have been consulted")
      case (SYM_SKEL_OFF)
         want_auto = .false.
         call check(error, .not. (ron%sym_skeleton .or. roff%sym_skeleton), &
                    label//": TRC_SYM_SKELETON=0 must reach every run")
         if (allocated(error)) return
         call check(error, ra%sym_atom_orbits, 0, &
                    label//": TRC_SYM_SKELETON=0 forced the mode, so the size rule must not have been consulted")
      case default
         want_auto = rule_wants
         call check(error, ron%sym_skeleton, label//": the skeleton build was asked for and not used")
         if (allocated(error)) return
         call check(error, .not. roff%sym_skeleton, label//": the full build was asked for and the skeleton was used")
         if (allocated(error)) return
         call check(error, ra%sym_atom_orbits, norb, label//": the rule read the wrong orbit count")
         if (allocated(error)) return
         call check(error, ron%sym_atom_orbits == 0 .and. roff%sym_atom_orbits == 0, &
                    label//": a forced request must not consult the size rule")
      end select
      if (allocated(error)) return
      call check(error, ra%sym_skeleton .eqv. want_auto, &
                 label//": the run in the rule's charge took the "//trim(merge("skeleton", "full    ", ra%sym_skeleton))// &
                 " build, expected the "//trim(merge("skeleton", "full    ", want_auto))//" one")
      if (allocated(error)) return

      ! And the three runs report the same molecule. The irrep labels, the
      ! group and the SO dimensions come from trc_sym_ao, which neither path
      ! touches; the energies and orbitals are the gate.
      call check(error, ra%n_lindep == roff%n_lindep .and. ron%n_lindep == roff%n_lindep, &
                 label//": the dropped-function count moved with the Fock build")
      if (allocated(error)) return
      call check(error, ra%iterations == roff%iterations .and. ron%iterations == roff%iterations, &
                 label//": the iteration count moved with the Fock build")
      if (allocated(error)) return
      nmo = size(roff%eps, 1) - roff%n_lindep
      ! The irrep decomposition, block by degenerate block rather than orbital
      ! by orbital: on SF6 nine pairs of orbitals swap labels inside their own
      ! degenerate level between the two Fock paths, so the per-orbital label
      ! is not an invariant of the path
      ! (this file's header has the measurement). TOL_EPS is also what the
      ! orbital energies themselves are compared to below, so a block here is
      ! exactly a set of levels this test cannot tell apart.
      do sp = 1, size(roff%mo_irrep, 2)
         call check(error, irrep_blocks_agree(roff%eps(1:nmo, sp), roff%mo_irrep(1:nmo, sp), &
                                              ra%mo_irrep(1:nmo, sp), TOL_EPS), &
                    label//": the orbital symmetries moved with the Fock build (rule's run against the full build)")
         if (allocated(error)) return
         call check(error, irrep_blocks_agree(roff%eps(1:nmo, sp), roff%mo_irrep(1:nmo, sp), &
                                              ron%mo_irrep(1:nmo, sp), TOL_EPS), &
                    label//": the orbital symmetries moved with the Fock build (skeleton against the full build)")
         if (allocated(error)) return
      end do
      de = max(abs(ra%energy - roff%energy), abs(ron%energy - roff%energy))
      deps = max(maxval(abs(ra%eps(1:nmo, :) - roff%eps(1:nmo, :))), &
                 maxval(abs(ron%eps(1:nmo, :) - roff%eps(1:nmo, :))))
      dd = max(maxval(abs(ra%dmat - roff%dmat)), maxval(abs(ron%dmat - roff%dmat)))
      print '(a,a,a,es10.2,a,es10.2,a,es10.2,a,i0)', '  ', label, ': |dE| ', de, '  max|d eps| ', deps, &
         '  max|dD| ', dd, '  iterations ', roff%iterations
      call check(error, de <= TOL_E, label//": the Fock build moved the energy by "//to_string(de, "es10.2"))
      if (allocated(error)) return
      call check(error, deps <= TOL_EPS, label//": orbital energies differ by "//to_string(deps, "es10.2"))
      if (allocated(error)) return
      call check(error, dd <= TOL_D, label//": densities differ by "//to_string(dd, "es10.2"))
      ! NOT asserted: that a declined run equals the forced-full run BIT FOR
      ! BIT because it is the same code path. From a fresh process it does
      ! (`terco tops_sf6.xyz --symmetry auto` and `--no-sym-skeleton` print the
      ! same 15 digits), but inside one process the three SCFs here do not
      ! reproduce each other exactly: two repeats of the identical forced-full
      ! run on SF6 came out 9.1e-12 Eh apart in this host build, the noise
      ! floor for a ~1000 Eh system (run-dependent accumulation order, 1e-13 to 1e-12 Eh).
      ! Comparing two runs bitwise is not meaningful at that level; the structural
      ! claim is carried by the `sym_skeleton` and `sym_atom_orbits` checks
      ! above instead, and the numbers by the gates above.
   end subroutine run_three

   subroutine test_declined(error)
      type(error_type), allocatable, intent(out) :: error
      call run_three(error, "hocl.xyz", "basis_sets/cc-pvdz.json", "HOCl Cs declined", .false.)
   end subroutine test_declined

   subroutine test_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      call run_three(error, "b2h6_d2h.xyz", "basis_sets/6-31gs.json", "B2H6 D2h accepted", .true.)
   end subroutine test_accepted

   subroutine test_sf6(error)
      !! The rule's own measured loss, and the only case here whose spectrum
      !! is degenerate: its sulfur is fixed by all eight operations of D2h and
      !! each fluorine by four, so 7 atoms sit in 4 orbits and the skeleton
      !! build is declined. It is in the cases because the degenerate levels
      !! are what makes the irrep comparison above non-trivial -- nine pairs of
      !! orbitals swap labels inside their own level between the two Fock
      !! paths, 18 of its 102 entries.
      type(error_type), allocatable, intent(out) :: error
      call run_three(error, "tops_sf6.xyz", "basis_sets/cc-pvdz.json", "SF6 D2h declined", .false.)
   end subroutine test_sf6

end module check_sym_skeleton_rule_suite

program check_sym_skeleton_rule
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_skeleton_rule_suite, only: collect_sym_skeleton_rule
   implicit none
   integer :: stat
   stat = 0
   ! Serial: the SCF cases build device-resident engines.
   call run_testsuite(collect_sym_skeleton_rule, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_skeleton_rule
