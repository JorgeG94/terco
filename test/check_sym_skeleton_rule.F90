!
! The skeleton Fock build's SIZE RULE (trc_sym_run's sym_skeleton_size_rule,
! read through sym_skeleton_mode and applied by trc_scf_run): which symmetric
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
! AND THAT DECLINING CHANGES NOTHING BUT THE WORK. HOCl is declined (3 atoms,
! 3 orbits: a planar Cs molecule has no atom orbit bigger than one) and
! B2H6 is accepted (8 atoms, 3 orbits). Each is run with the rule in charge,
! with the skeleton build forced on and with it forced off, and the energy, the
! orbital energies, the density, the irrep labels, the group, the SO block
! dimensions and the iteration count must be the same every time. The run the
! rule declines must also be bit-for-bit the run that was asked to decline:
! it is the same code path.
!
! THE ENVIRONMENT OVERRIDE is tested by the ctest entries, which run this
! program three times: without TRC_SYM_SKELETON, with 1 and with 0. The
! expectations below follow whatever the environment says, so the override is
! checked end to end, through the SCF and not only through the mode function.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_skeleton_rule_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect, sym_geom_atom_orbits
   use trc_sym_run, only: sym_skeleton_mode, sym_skeleton_size_rule, &
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
                  new_unittest("the threshold itself, at and below the boundary", test_boundary), &
                  new_unittest("TRC_SYM_SKELETON overrides the request", test_env), &
                  new_unittest("HOCl cc-pVDZ (declined): the same SCF three ways", test_declined), &
                  new_unittest("B2H6 6-31G* (accepted): the same SCF three ways", test_accepted) &
                  ]
   end subroutine collect_sym_skeleton_rule

   subroutine one_count(error, xyz, group, natm, norbit, verdict)
      !! Detect `xyz`'s group, count its atom orbits and check the rule's
      !! verdict against what sym_skeleton_size_rule's table says.
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
      geom = sym_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), xyz//": sym_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, group, xyz//": group")
      if (allocated(error)) return
      call check(error, nat, natm, xyz//": atoms")
      if (allocated(error)) return
      norb = sym_geom_atom_orbits(geom%perm)
      call check(error, norb, norbit, xyz//": atoms inequivalent under "//group)
      if (allocated(error)) return
      call check(error, sym_skeleton_size_rule(nat, norb) .eqv. verdict, &
                 xyz//": the size rule wants the "//trim(merge("skeleton", "full    ", verdict))// &
                 " build at "//to_string(nat)//" atoms in "//to_string(norb)//" orbits")
   end subroutine one_count

   subroutine test_counts(error)
      !! The rows of sym_skeleton_size_rule's measured table, plus the two
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

   subroutine test_boundary(error)
      !! The threshold is a mean orbit size of SYM_SKEL_MIN_ORBIT, exactly.
      type(error_type), allocatable, intent(out) :: error
      integer :: n
      n = SYM_SKEL_MIN_ORBIT
      call check(error, sym_skeleton_size_rule(5*n, 5), "5 orbits of exactly the threshold size are accepted")
      if (allocated(error)) return
      call check(error, .not. sym_skeleton_size_rule(5*n - 1, 5), "one atom below the threshold is declined")
      if (allocated(error)) return
      call check(error, sym_skeleton_size_rule(5*n + 1, 5), "one atom above the threshold is accepted")
      if (allocated(error)) return
      ! A count the detection never produces must not be read as a verdict.
      call check(error, .not. sym_skeleton_size_rule(0, 0), "an empty molecule is declined, not accepted")
   end subroutine test_boundary

   subroutine test_env(error)
      !! sym_skeleton_mode returns the request unless TRC_SYM_SKELETON names a
      !! mode, and then that mode whatever was requested.
      type(error_type), allocatable, intent(out) :: error
      integer :: forced
      character(len=16) :: e
      e = " "
      call get_environment_variable("TRC_SYM_SKELETON", e)
      forced = sym_skeleton_mode(SYM_SKEL_AUTO)
      select case (trim(e))
      case ("1")
         call check(error, forced, SYM_SKEL_ON, "TRC_SYM_SKELETON=1 must force the skeleton build on")
         if (allocated(error)) return
         call check(error, sym_skeleton_mode(SYM_SKEL_OFF), SYM_SKEL_ON, "and must override a request to decline")
      case ("0")
         call check(error, forced, SYM_SKEL_OFF, "TRC_SYM_SKELETON=0 must force the skeleton build off")
         if (allocated(error)) return
         call check(error, sym_skeleton_mode(SYM_SKEL_ON), SYM_SKEL_OFF, "and must override a request to use it")
      case default
         call check(error, forced, SYM_SKEL_AUTO, "without TRC_SYM_SKELETON the request stands")
         if (allocated(error)) return
         call check(error, sym_skeleton_mode(SYM_SKEL_ON), SYM_SKEL_ON, "and a forced request is passed through")
      end select
   end subroutine test_env

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
      integer :: nat, nelec, nalpha, nbeta, dev, norb, nmo
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      logical :: want_auto
      real(dp) :: de, deps, dd

      dev = trc_bind_device(0)
      call read_xyz(xyz, nat, elem, r)
      geom = sym_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), xyz//": sym_detect: "//err%get_message())
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
      select case (sym_skeleton_mode(SYM_SKEL_AUTO))
      case (SYM_SKEL_ON)
         want_auto = .true.
         call check(error, ron%sym_skeleton .and. roff%sym_skeleton, label//": TRC_SYM_SKELETON=1 must reach every run")
      case (SYM_SKEL_OFF)
         want_auto = .false.
         call check(error, .not. (ron%sym_skeleton .or. roff%sym_skeleton), &
                    label//": TRC_SYM_SKELETON=0 must reach every run")
      case default
         want_auto = rule_wants
         call check(error, ron%sym_skeleton, label//": the skeleton build was asked for and not used")
         if (allocated(error)) return
         call check(error, .not. roff%sym_skeleton, label//": the full build was asked for and the skeleton was used")
         if (allocated(error)) return
         call check(error, ra%sym_atom_orbits, norb, label//": the rule read the wrong orbit count")
      end select
      if (allocated(error)) return
      call check(error, ra%sym_skeleton .eqv. want_auto, &
                 label//": the run in the rule's charge took the "//trim(merge("skeleton", "full    ", ra%sym_skeleton))// &
                 " build, expected the "//trim(merge("skeleton", "full    ", want_auto))//" one")
      if (allocated(error)) return

      ! And the three runs report the same molecule. The irrep labels, the
      ! group and the SO dimensions come from trc_sym_ao, which neither path
      ! touches; the energies and orbitals are the gate.
      call check(error, all(ra%mo_irrep == roff%mo_irrep) .and. all(ron%mo_irrep == roff%mo_irrep), &
                 label//": the orbital symmetries moved with the Fock build")
      if (allocated(error)) return
      call check(error, ra%n_lindep == roff%n_lindep .and. ron%n_lindep == roff%n_lindep, &
                 label//": the dropped-function count moved with the Fock build")
      if (allocated(error)) return
      call check(error, ra%iterations == roff%iterations .and. ron%iterations == roff%iterations, &
                 label//": the iteration count moved with the Fock build")
      if (allocated(error)) return
      nmo = size(roff%eps, 1) - roff%n_lindep
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
      if (allocated(error)) return
      ! A declined run is not merely close to the run that was told to
      ! decline: it is the same code path on the same engine.
      if (.not. ra%sym_skeleton) then
         call check(error, ra%energy == roff%energy, &
                    label//": a declined skeleton build is the full build and must give its energy exactly")
      end if
   end subroutine run_three

   subroutine test_declined(error)
      type(error_type), allocatable, intent(out) :: error
      call run_three(error, "hocl.xyz", "basis_sets/cc-pvdz.json", "HOCl Cs declined", .false.)
   end subroutine test_declined

   subroutine test_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      call run_three(error, "b2h6_d2h.xyz", "basis_sets/6-31gs.json", "B2H6 D2h accepted", .true.)
   end subroutine test_accepted

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
