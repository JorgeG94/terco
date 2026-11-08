module trc_sym_run
   !! The user-facing side of symmetry: what a run asks for (`off`, `auto`, or
   !! a named group), which group that means for a detected geometry, and the
   !! per-irrep occupations a user can fix. The detection, the SO basis and
   !! the SCF are `trc_sym_geom`, `trc_sym_ao` and `trc_sym_scf`; this module
   !! only decides and parses, so the CLI, the Fortran API and the C ABI
   !! share one set of rules and one set of messages.
   !!
   !! RULES. `off` never reaches the symmetry code. `auto` uses
   !! the largest enabled subgroup of the detected group (`trc_sym_geom`
   !! applies the cap). `default` is `auto` asked for by nobody: the request
   !! of a fresh C-ABI context and of a CLI run without `--symmetry` (S17).
   !! It differs from an explicit `auto` in two ways: the results are
   !! reported in the user's input frame, and a path the symmetric SCF does
   !! not support runs in C1 instead of being refused (`by_default`). The
   !! Fortran library entry `trc_scf_run` has no default: `sym` is absent
   !! unless the caller passes it, and `trc_sym_request_t` itself starts at
   !! `SYM_OFF`. A named group must be the detected group or one of
   !! its subgroups, else the run is refused. The caller passes the named
   !! group to `sym_detect(..., requested=)`, which reports it in pyscf's
   !! frame for that subgroup; a group the build does not
   !! support yet (`SYM_GROUP_ENABLED`) is refused with its own message, after
   !! the subgroup test so that a wrong group is reported as wrong first.
   !! `C1` is always allowed and means "no symmetry in the SCF".
   !!
   !! IRREP OCCUPATIONS are given as electrons per irrep, in pyscf's
   !! `irrep_nelec` spirit: "A'=6,A\"=2" for a restricted run (even counts,
   !! each orbital holds two), "A'=4/3,A\"=1/1" for an unrestricted one
   !! (alpha/beta). Irreps that are not listed hold no electrons. The text is
   !! parsed against the irrep names of the group actually used, so a typo is
   !! reported with the valid names.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_sym_group, only: trc_sym_group_t, sym_group_from_name, SYM_NGROUPS, sym_group_names
   use trc_sym_geom, only: trc_sym_geom_t, SYM_GROUP_ENABLED, SYM_DEFAULT_TOL
   implicit none
   private

   public :: trc_sym_request_t, SYM_OFF, SYM_AUTO, SYM_EXPLICIT
   public :: SYM_SKEL_AUTO, SYM_SKEL_ON, SYM_SKEL_OFF
   public :: sym_parse_request, sym_choose_group, sym_is_subgroup, sym_request_text
   public :: sym_parse_irrep_nelec, sym_aufbau_matches, sym_fixed_occupation, sym_occupied_counts
   public :: sym_skeleton_mode, sym_skeleton_size_rule, SYM_SKEL_MIN_ORBIT

   integer, parameter :: SYM_OFF = 0       !! no symmetry
   integer, parameter :: SYM_AUTO = 1      !! the detected group, capped to what is enabled
   integer, parameter :: SYM_EXPLICIT = 2  !! a named group

   integer, parameter :: SYM_SKEL_AUTO = 0  !! the skeleton Fock build where `sym_skeleton_size_rule` says it pays
   integer, parameter :: SYM_SKEL_ON = 1    !! the skeleton Fock build whatever the rule says
   integer, parameter :: SYM_SKEL_OFF = 2   !! never the skeleton Fock build: the full build, then symmetrize

   integer, parameter :: SYM_SKEL_MIN_ORBIT = 2
      !! The threshold of `sym_skeleton_size_rule`: the smallest mean atom
      !! orbit size (atoms per symmetry-inequivalent atom) at which the
      !! skeleton Fock build is used. See that routine for the measurements.

   type :: trc_sym_request_t
      !! What the user asked for.
      integer :: mode = SYM_OFF
      character(len=3) :: group = "C1"  !! canonical name, with `SYM_EXPLICIT`
      real(dp) :: tol = SYM_DEFAULT_TOL  !! geometry tolerance, Bohr
      character(len=:), allocatable :: irrep_nelec  !! electrons per irrep, unparsed; absent = aufbau
      integer :: skeleton = SYM_SKEL_AUTO
      !! Whether the Fock build is the skeleton one (trc_sym_scf_t%skeleton):
      !! SYM_SKEL_AUTO leaves it to `sym_skeleton_size_rule`, SYM_SKEL_ON and
      !! SYM_SKEL_OFF force it either way.
      logical :: by_default = .false.  !! `auto` nobody asked for (`default`): input-frame report, C1 fallback
   end type trc_sym_request_t

contains

   pure function lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, c
      t = s
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) t(i:i) = achar(c + 32)
      end do
   end function lower

   pure function sym_request_text(req) result(text)
      !! The request as the user wrote it: "off", "auto", "default" or the group name.
      type(trc_sym_request_t), intent(in) :: req
      character(len=:), allocatable :: text
      select case (req%mode)
      case (SYM_AUTO)
         text = merge("default", "auto   ", req%by_default)
         text = trim(text)
      case (SYM_EXPLICIT)
         text = trim(req%group)
      case default
         text = "off"
      end select
   end function sym_request_text

   subroutine sym_parse_request(text, req, error)
      !! `off`, `auto`, `default` or a group name (any case) into
      !! `req%mode`/`req%group`/`req%by_default`. Leaves `req%tol` and
      !! `req%irrep_nelec` alone.
      character(len=*), intent(in) :: text
      type(trc_sym_request_t), intent(inout) :: req
      type(error_t), intent(out) :: error
      type(trc_sym_group_t) :: g
      character(len=:), allocatable :: t

      t = lower(trim(adjustl(text)))
      req%by_default = .false.
      select case (t)
      case ("off", "none", "")
         req%mode = SYM_OFF
         req%group = "C1"
      case ("auto")
         req%mode = SYM_AUTO
         req%group = "C1"
      case ("default")
         req%mode = SYM_AUTO
         req%group = "C1"
         req%by_default = .true.
      case default
         g = sym_group_from_name(t, error)
         if (error%has_error()) then
            call error%set(ERROR_VALIDATION, "symmetry must be off, auto or a point group (D2h, C2h, C2v, D2, Cs, Ci, "// &
                           "C2, C1), got '"//trim(adjustl(text))//"'")
            return
         end if
         req%mode = SYM_EXPLICIT
         req%group = g%name
      end select
   end subroutine sym_parse_request

   pure logical function sym_is_subgroup(sub, group) result(ok)
      !! Whether the group called `sub` is `group` or a subgroup of it, among
      !! the Abelian subgroups of D2h.
      character(len=*), intent(in) :: sub, group
      character(len=3) :: s
      s = trim(sub)
      select case (trim(group))
      case ("D2h")
         ok = any(s == sym_group_names)
      case ("C2h")
         ok = any(s == ["C2h", "Cs ", "Ci ", "C2 ", "C1 "])
      case ("C2v")
         ok = any(s == ["C2v", "Cs ", "C2 ", "C1 "])
      case ("D2")
         ok = any(s == ["D2 ", "C2 ", "C1 "])
      case ("Cs", "Ci", "C2")
         ok = s == trim(group) .or. s == "C1"
      case default
         ok = s == "C1"
      end select
   end function sym_is_subgroup

   pure logical function group_enabled(name) result(on)
      character(len=*), intent(in) :: name
      integer :: ig
      on = .false.
      do ig = 1, SYM_NGROUPS
         if (trim(sym_group_names(ig)) == trim(name)) on = SYM_GROUP_ENABLED(ig)
      end do
   end function group_enabled

   function sym_choose_group(req, geom, error) result(name)
      !! The group to run in, for the geometry `geom` that `sym_detect` found:
      !! the reported group for `auto`, the requested one otherwise (after the
      !! subgroup and support checks). "C1" means no symmetry in the SCF.
      type(trc_sym_request_t), intent(in) :: req
      type(trc_sym_geom_t), intent(in) :: geom
      type(error_t), intent(out) :: error
      character(len=3) :: name
      character(len=:), allocatable :: why

      name = "C1"
      select case (req%mode)
      case (SYM_AUTO)
         name = geom%group%name
      case (SYM_EXPLICIT)
         if (.not. sym_is_subgroup(req%group, geom%detected)) then
            why = ""
            if (geom%degenerate) why = " (a symmetric or spherical top, linear molecule or atom: requested "// &
                                       "subgroups of these are not supported, use symmetry = auto)"
            call error%set(ERROR_VALIDATION, "symmetry "//trim(req%group)//": the molecule's point group is "// &
                           trim(geom%detected)//trim(why)//", and "//trim(req%group)// &
                           " is not that group or one of its subgroups; use symmetry = auto or a subgroup")
            return
         end if
         if (.not. group_enabled(req%group)) then
            call error%set(ERROR_VALIDATION, "symmetry "//trim(req%group)//" is not yet supported (the molecule is "// &
                           trim(geom%detected)//"; supported so far: "//enabled_list()//"; symmetry = auto uses the "// &
                           "largest supported subgroup)")
            return
         end if
         name = req%group
         ! `sym_detect(..., requested=req%group)` reports a requested subgroup that has a
         ! frame rule in `trc_sym_geom`. A requested group that it did not report
         ! (no rule for this pair yet) can only be C1, which runs without symmetry.
         if (name /= geom%group%name .and. name /= "C1") then
            call error%set(ERROR_VALIDATION, "symmetry "//trim(req%group)//" cannot be used for this molecule "// &
                           "(the molecule is "//trim(geom%detected)//", the group found is "// &
                           trim(geom%group%name)//", and no frame is defined for "//trim(req%group)// &
                           " as a subgroup of "//trim(geom%detected)//"; pass the request to sym_detect)")
         end if
      case default
         name = "C1"
      end select
   end function sym_choose_group

   subroutine sym_parse_irrep_nelec(text, group, nspin, nocc, error)
      !! `text` ("A'=6,A\"=2", or "A'=4/3,A\"=1/1" with two spin channels)
      !! into the orbitals per irrep and spin, `nocc(nirrep, nspin)`. Irreps
      !! that are not listed get 0.
      character(len=*), intent(in) :: text
      type(trc_sym_group_t), intent(in) :: group
      integer, intent(in) :: nspin
      integer, allocatable, intent(out) :: nocc(:, :)
      type(error_t), intent(out) :: error
      integer :: start, comma, eq, slash, ir, s, ios, na, nb, n
      logical, allocatable :: seen(:)
      character(len=:), allocatable :: item, label, val

      allocate (nocc(group%order, nspin), seen(group%order))
      nocc = 0
      seen = .false.
      start = 1
      do while (start <= len(text))
         comma = index(text(start:), ",")
         if (comma == 0) then
            item = trim(adjustl(text(start:)))
            start = len(text) + 1
         else
            item = trim(adjustl(text(start:start + comma - 2)))
            start = start + comma
         end if
         if (len(item) == 0) cycle
         eq = index(item, "=")
         if (eq < 2) then
            call error%set(ERROR_VALIDATION, "irrep_nelec: '"//item//"' is not LABEL=ELECTRONS")
            return
         end if
         label = trim(adjustl(item(1:eq - 1)))
         val = trim(adjustl(item(eq + 1:)))
         ir = 0
         do s = 1, group%order
            if (lower(trim(group%irrep_names(s))) == lower(label)) ir = s
         end do
         if (ir == 0) then
            call error%set(ERROR_VALIDATION, "irrep_nelec: '"//label//"' is not an irrep of "//trim(group%name)// &
                           " (the irreps are "//names_list(group)//")")
            return
         end if
         if (seen(ir)) then
            call error%set(ERROR_VALIDATION, "irrep_nelec: irrep "//label//" is given twice")
            return
         end if
         seen(ir) = .true.
         slash = index(val, "/")
         if (nspin == 1) then
            if (slash > 0) then
               call error%set(ERROR_VALIDATION, "irrep_nelec: a restricted run takes one electron count per irrep "// &
                              "(A'=6), not alpha/beta ('"//item//"')")
               return
            end if
            read (val, *, iostat=ios) n
            if (ios /= 0 .or. n < 0 .or. mod(n, 2) /= 0) then
               call error%set(ERROR_VALIDATION, "irrep_nelec: '"//item//"' must be a non-negative even number of "// &
                              "electrons in a restricted run")
               return
            end if
            nocc(ir, 1) = n/2
         else
            if (slash < 2) then
               call error%set(ERROR_VALIDATION, "irrep_nelec: an unrestricted run takes alpha/beta electrons per "// &
                              "irrep (A'=4/3), got '"//item//"'")
               return
            end if
            read (val(1:slash - 1), *, iostat=ios) na
            if (ios == 0) read (val(slash + 1:), *, iostat=ios) nb
            if (ios /= 0 .or. na < 0 .or. nb < 0) then
               call error%set(ERROR_VALIDATION, "irrep_nelec: '"//item//"' must be non-negative alpha/beta counts")
               return
            end if
            nocc(ir, 1) = na
            nocc(ir, 2) = nb
         end if
      end do
   end subroutine sym_parse_irrep_nelec

   function enabled_list() result(text)
      !! The enabled groups (`SYM_GROUP_ENABLED`), comma separated, for messages.
      character(len=:), allocatable :: text
      integer :: ig
      text = "C1"
      do ig = SYM_NGROUPS - 1, 1, -1
         if (SYM_GROUP_ENABLED(ig)) text = text//", "//trim(sym_group_names(ig))
      end do
   end function enabled_list

   function names_list(group) result(text)
      type(trc_sym_group_t), intent(in) :: group
      character(len=:), allocatable :: text
      integer :: i
      text = ""
      do i = 1, group%order
         if (i > 1) text = text//", "
         text = text//trim(group%irrep_names(i))
      end do
   end function names_list

   pure logical function sym_aufbau_matches(mo_irrep, nocc, irrep_nocc) result(ok)
      !! Whether fixed occupations `irrep_nocc(nirrep, nspin)` are the aufbau
      !! ones of the orbitals labelled `mo_irrep(nmo, nspin)`: the lowest
      !! `nocc(s)` orbitals then hold exactly `irrep_nocc(:, s)` per irrep.
      !! The gradient rebuilds its density from the first `nocc` orbitals, so
      !! it is right only when this is true.
      integer, intent(in) :: mo_irrep(:, :), nocc(:), irrep_nocc(:, :)
      integer :: s, ir
      ok = .true.
      do s = 1, size(nocc)
         do ir = 1, size(irrep_nocc, 1)
            if (count(mo_irrep(1:nocc(s), s) == ir) /= irrep_nocc(ir, s)) ok = .false.
         end do
      end do
   end function sym_aufbau_matches

   pure function sym_fixed_occupation(mo_irrep, irrep_nocc) result(occupied)
      !! Which of the orbitals of one spin (labels `mo_irrep(nmo)`, in energy
      !! order) are occupied when each irrep holds its lowest `irrep_nocc(ir)`.
      integer, intent(in) :: mo_irrep(:), irrep_nocc(:)
      logical :: occupied(size(mo_irrep))
      integer :: m, ir
      integer :: seen(size(irrep_nocc))
      seen = 0
      do m = 1, size(mo_irrep)
         ir = mo_irrep(m)
         occupied(m) = .false.
         if (ir < 1) cycle
         seen(ir) = seen(ir) + 1
         occupied(m) = seen(ir) <= irrep_nocc(ir)
      end do
   end function sym_fixed_occupation

   pure function sym_occupied_counts(mo_irrep, occupied, nirrep) result(counts)
      !! Occupied orbitals per irrep, from the labels and an occupation mask.
      integer, intent(in) :: mo_irrep(:)
      integer, intent(in) :: nirrep
      logical, intent(in) :: occupied(:)
      integer :: counts(nirrep)
      integer :: ir
      do ir = 1, nirrep
         counts(ir) = count(occupied .and. mo_irrep == ir)
      end do
   end function sym_occupied_counts

   integer function sym_skeleton_mode(request) result(mode)
      !! The skeleton Fock build's mode in force: `request`
      !! (`trc_sym_request_t%skeleton`, i.e. what the user asked for), unless
      !! TRC_SYM_SKELETON names one. TRC_SYM_SKELETON=1 forces the skeleton
      !! build on and 0 forces it off, so both sides of
      !! `sym_skeleton_size_rule` stay measurable from one binary without
      !! rebuilding; `auto` puts the rule back in charge. This mirrors TRC_MDJ
      !! over the J engine's own size rule (trc_eri's mdj_requested).
      integer, intent(in) :: request
      character(len=8) :: e
      mode = request
      e = " "
      call get_environment_variable("TRC_SYM_SKELETON", e)
      select case (trim(e))
      case ("1", "on")
         mode = SYM_SKEL_ON
      case ("0", "off")
         mode = SYM_SKEL_OFF
      case ("auto")
         mode = SYM_SKEL_AUTO
      case default
         ! Unset, blank, or a word nobody documented: the request stands.
         ! A typo must not silently turn a measured default into a forced run.
      end select
   end function sym_skeleton_mode

   pure logical function sym_skeleton_size_rule(natm, norbit) result(res)
      !! Whether the skeleton Fock build (docs/SYMMETRY_FOCK.md) is worth
      !! using on a molecule of `natm` atoms sitting in `norbit` orbits of the
      !! group (`sym_geom_atom_orbits` of trc_sym_geom). It is, when the mean
      !! atom orbit size is at least SYM_SKEL_MIN_ORBIT, i.e. when
      !!   natm >= 2 * norbit.
      !! Nothing else about the molecule enters: not the basis, not the atom
      !! count on its own.
      !!
      !! WHY ORBIT SIZE AND NOT SIZE. The skeleton build computes one shell
      !! quartet per orbit and weights it, so what it SAVES is the orbit
      !! reduction, and what it COSTS is the extra launches the weights and
      !! `bins_skel_cut_first`'s bra/ket swap create. The mean atom orbit size
      !! is a cheap proxy for the saving, available before the engine is built
      !! (the decision has to be made then: the engine is built with the group
      !! or without it). It tracks the measured scheme-2a quartet reduction
      !! across every system the project has counted -- docs/SYMMETRY_FOCK.md
      !! section 2.2, cc-pVDZ, tau = 1e-10, density-screened:
      !!
      !!   system       group  natm/norbit  s2a quartet reduction
      !!   water         C2v      1.50          1.47x
      !!   tops_sf6      D2h      1.75          2.44x
      !!   b2h6_d2h      D2h      2.67          3.34x
      !!   benzene       D2h      3.00          3.53x
      !!   naphthalene   D2h      3.60          3.79x
      !!   pentacene     D2h      3.60          3.85x
      !!
      !! A molecule whose atoms sit on many symmetry elements at once has
      !! small orbits however many atoms or functions it has: SF6's sulfur is
      !! fixed by all eight operations of D2h and each fluorine by four, and a
      !! planar molecule's own plane fixes every atom (which is why planar D2h
      !! reaches an effective order of 4 and not 8). That is also why the
      !! basis does not enter: S15 measured benzene at 0.46x (cc-pVDZ, 114 bf)
      !! and 0.48x (cc-pVTZ, 264 bf) and naphthalene at 0.98x and 0.99x across
      !! the same 2.3x increase, while naphthalene at 412 bf was a wash and
      !! pentacene at 378 bf a 2.14x win.
      !!
      !! Measured ERI time (cc-pVDZ, --econv 1e-10, Tesla V100-DGXS-32GB
      !! (cc70), nvfortran 26.5, -stdpar=gpu -acc -gpu=cc70,mem:separate,
      !! TERCO_LMAX=3, per steady-state Fock build, baseline the same build
      !! with --symmetry off, W2 of .claude/context/segment-fusion-plan.md
      !! landed). `*` marks what this rule chooses:
      !!
      !!   system            natm  norbit  natm/norbit  ERI on/off
      !!   tops_sf6             7       4     1.75        0.44x   (full*)
      !!   benzene             12       4     3.00        1.41x   (skeleton*)
      !!   naphthalene_d2h     18       5     3.60        2.07x   (skeleton*)
      !!   pentacene_d2h       36      10     3.60        2.78x   (skeleton*)
      !!
      !! Energies are unaffected either way: |dE| between the symmetric and
      !! the C1 run on the same symmetrized geometry was 2.9e-13 to 1.4e-12 Eh
      !! against a 1e-10 gate, with identical iteration counts (S15).
      !!
      !! WHERE THE THRESHOLD COMES FROM. SF6 loses with a 2.44x reduction and
      !! benzene wins with 3.53x, so break-even is a reduction somewhere in
      !! (2.44x, 3.53x), i.e. a mean orbit size somewhere in (1.75, 3.00).
      !! The threshold is set at the LOW end of that interval, so the rule
      !! declines only what was measured to lose and keeps the skeleton build
      !! wherever it might still pay -- an over-eager rule would now cost time,
      !! since post-W2 benzene at 12 atoms is a win. B2H6 (2.67, 3.34x) is the
      !! one system inside the interval, and the run that would settle the
      !! threshold; docs/SYMMETRY_FOCK.md section 4.2a lists it and the rest.
      !!
      !! REJECTED PREDICTORS, on this same data:
      !!  - the atom count: the four ERI rows allow a threshold anywhere in
      !!    (7, 12], and every such threshold also declines 8-atom b2h6_d2h,
      !!    whose 3.34x reduction is closer to benzene's than to SF6's. The
      !!    pre-W2 crossover near twenty atoms is stale: benzene is now a win.
      !!  - the number of inequivalent atoms, `norbit`: it does not separate
      !!    the data at all (SF6 and benzene both have 4).
      !!  - predicted quartets per launch against the device's fill threshold
      !!    (~164k work items, 80 SMs at 128-thread blocks): W0 measured
      !!    benzene and SF6 at the SAME ~46% of quartets in launches too small
      !!    to fill the card after fusion, yet benzene gains 1.41x and SF6
      !!    loses 2.3x, so fill is not what separates them. The reduction
      !!    (3.55x against 2.44x) and the launch inflation over the C1 plan
      !!    (1.57x against 3.43x) are, and both follow from the orbit sizes.
      !!
      !! The groups of order 2 and 4 (Cs, Ci, C2, C2v, C2h) are declined
      !! wherever an atom sits on a symmetry element, which is most of them --
      !! water at 1.47x, and a planar Cs molecule has no nontrivial shell-pair
      !! orbit at all. That answers docs/SYMMETRY_FOCK.md section 6's open
      !! question 2 from measurement: a ~1.5x reduction cannot pay for the
      !! launches when 2.44x does not.
      !!
      !! TRC_SYM_SKELETON=0/1 forces either path.
      integer, intent(in) :: natm, norbit
      res = .false.
      if (natm < 1 .or. norbit < 1) return
      res = natm >= SYM_SKEL_MIN_ORBIT*norbit
   end function sym_skeleton_size_rule

end module trc_sym_run
