submodule (trc_run) trc_run_opt
   !! The `--opt` / NEB support of the `terco` command: option parsing and
   !! checks, the optimiser summary and the xyz writers (bodies of `trc_run`).
   use pic_types, only: default_int
   use pic_error, only: ERROR_VALIDATION, ERROR_IO, ERROR_PARSE
   use pic_strings, only: to_string, strip, find
   use pic_string_type, only: string_type, char
   use pic_tokenizer, only: split, tokenize, parse_int
   use trc_app, only: app_element_symbol, app_atomic_number, BOHR_PER_ANGSTROM, app_write_opt_xyz_frame, &
                      app_write_traj_xyz_frames
   use trc_vib, only: trc_vib_rotor
   use trc_dlfind_bridge, only: TRC_OPT_STEP_RULE
   use trc_constants, only: trc_atomic_mass, trc_atomic_mass_average, trc_covalent_radius
   implicit none

contains

   !
   ! Per-atom masses for --hessian's vibrational analysis: "isotope"
   ! (trc_atomic_mass, the most-abundant isotope -- the CLI's own
   ! default) or "average" (trc_atomic_mass_average, the isotope-averaged
   ! standard atomic weight). `kind` is already lower-cased by the
   ! caller; anything else is a programming error (the CLI validates
   ! --masses before this is ever called).
   !
   module subroutine run_build_masses(z, kind, mass, err)
      integer, intent(in) :: z(:)
      character(len=*), intent(in) :: kind
      real(dp), allocatable, intent(out) :: mass(:)
      type(mass_error_t), intent(inout) :: err
      integer :: i
      allocate (mass(size(z)))
      do i = 1, size(z)
         if (kind == "average") then
            mass(i) = trc_atomic_mass_average(z(i), err)
         else
            mass(i) = trc_atomic_mass(z(i), err)
         end if
         if (err%has_error()) return
      end do
   end subroutine run_build_masses

   !
   ! Whether this run's method has a gradient at all, checked BEFORE any
   ! SCF runs: every one of these conditions is already known from the
   ! command line alone, so there is no need to spend an SCF on a probe
   ! gradient just to find out. The same set trc_run_gradient's own
   ! gradient_capability_check refuses after a converged SCF
   ! (src/api/trc_capi.F90). Returns "" when a gradient is available, the
   ! refusal message otherwise.
   !
   module function run_opt_capability_check(r) result(refusal)
      use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_FAMILY_MGGA
      type(run_t), intent(in) :: r
      character(len=:), allocatable :: refusal
      type(trc_xc_functional_t) :: func
      refusal = ""
      if (len(r%func) > 0) then
         if (xc_functional_by_name(r%func, func)) then
            if (func%family == XC_FAMILY_MGGA) then
               refusal = "--opt needs a gradient at every step; meta-GGA ('"//trim(r%func)// &
                        "') exchange-correlation gradient not implemented"
               return
            end if
         end if
      end if
      if (r%rimp2 .or. r%dh) then
         refusal = "--opt needs a gradient at every step; RI-MP2/double-hybrid correlation has none"
         return
      end if
      if (len(r%pcm) > 0) then
         refusal = "--opt needs a gradient at every step; PCM gradient not implemented"
         return
      end if
   end function run_opt_capability_check

   !
   ! Whether the start geometry is linear, for --opt-zero-modes' own
   ! auto-default (5 for a linear molecule's Cartesian P-RFO/NR/dimer
   ! search, 6 otherwise -- DL-FIND's own zero-mode count for the rigid
   ! translations/rotations a Cartesian search carries). Reuses
   ! trc_vib's own moment-of-inertia rotor test (trc_vib_rotor,
   ! MOMENT_LINEAR_TOL_AMUA2) rather than a separate geometric criterion,
   ! so "linear" means the same thing everywhere in terco: trc_thermo's
   ! own is_atom/is_linear picks its rotor formulas from the SAME
   ! threshold, and this and that can never disagree on the same
   ! geometry. trc_vib_rotor needs only masses and Cartesian coordinates,
   ! no Hessian -- what this is available with at --opt's own start.
   !
   ! This runs ONCE, on the start geometry: DL-FIND's own interface takes
   ! `nzero` as a fixed parameter at cb_get_params (called once, before
   ! the search begins, see trc_dlfind_bridge), so a search that bends a
   ! started-linear molecule (or straightens a bent one) cannot have its
   ! zero-mode count re-evaluated mid-search through that interface.
   !
   module function run_opt_geometry_is_linear(natm, z, xyz, err) result(linear)
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: xyz(3, natm)
      type(mass_error_t), intent(inout) :: err
      logical :: linear
      real(dp), allocatable :: mass(:)
      integer :: n_tr, info
      linear = .false.
      call run_build_masses(z, "isotope", mass, err)
      if (err%has_error()) return
      call trc_vib_rotor(mass, xyz, linear, n_tr, info)
   end function run_opt_geometry_is_linear

   !
   ! Connected components of the covalent-bond graph (Cordero radii,
   ! trc_constants%trc_covalent_radius, bond threshold 1.3x the sum of two
   ! atoms' radii), 1-based, for --opt-coords hdlc's `residues`: one
   ! residue per covalently bonded molecule, so HDLC builds internal
   ! coordinates within each and Cartesians between them. --opt-coords
   ! cartesian/dlc never calls this -- trc_dlfind_bridge's driver takes
   ! all-ones directly, which is correct for cartesian and for a single
   ! bonded molecule in dlc (a multi-fragment cluster in pure dlc fails
   ! inside DL-FIND itself with "cyclic failure at residue 1" -- fixing
   ! that would need a real multi-residue partition for dlc too, out of
   ! scope here).
   !
   module subroutine run_build_opt_residues(natm, z, xyz, residues, err)
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: xyz(3, natm)
      integer, intent(out) :: residues(natm)
      type(error_t), intent(inout) :: err
      real(dp), parameter :: BOND_SCALE = 1.3_dp
      real(dp) :: radius(natm)
      type(mass_error_t) :: rerr
      integer :: i, j, comp, sp, a
      integer :: stack(natm)
      logical :: visited(natm)
      do i = 1, natm
         radius(i) = trc_covalent_radius(z(i), rerr)
         if (rerr%has_error()) then
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-coords hdlc: "//rerr%get_message())
            return
         end if
      end do
      residues = 0
      visited = .false.
      comp = 0
      do i = 1, natm
         if (visited(i)) cycle
         comp = comp + 1
         sp = 1
         stack(1) = i
         visited(i) = .true.
         do while (sp > 0)
            a = stack(sp)
            sp = sp - 1
            residues(a) = comp
            do j = 1, natm
               if (visited(j) .or. j == a) cycle
               if (norm2(xyz(:, a) - xyz(:, j)) <= BOND_SCALE*(radius(a) + radius(j))) then
                  visited(j) = .true.
                  sp = sp + 1
                  stack(sp) = j
               end if
            end do
         end do
      end do
   end subroutine run_build_opt_residues

   !
   ! --opt-seed and --opt-global-*: refused unless --opt-algorithm is
   ! stochastic or genetic (r%opt_is_global, set above), so a command
   ! line that asks for a global-search parameter without the search
   ! that reads it never means something different from what it says --
   ! the same convention --opt-algorithm itself follows against
   ! --opt-target saddle. Validated here (positive, and distribution in
   ! DL-FIND's own {1,2,3}) rather than left for trc_dlfind_bridge to
   ! pass straight through to libdlfind. --opt-seed alone is ALSO
   ! accepted for --opt-saddle-method dimer: trc_dlfind_bridge seeds the
   ! Fortran intrinsic generator dlf_dimer.f90's own random_number(
   ! glob%xcoords2) reads, the only way to make that initial direction
   ! reproducible (that module's own header).
   !
   module subroutine run_parse_opt_global_options(r)
      type(run_t), intent(inout) :: r
      integer, parameter :: ALG_STOCHASTIC = 51, ALG_GENETIC = 52
         !! trc_dlfind_bridge's own TRC_OPT_ALG_STOCHASTIC/GENETIC values,
         !! not imported directly (this module builds without
         !! trc_dlfind_bridge under TERCO_ENABLE_DLFIND=OFF) -- the same
         !! "plain integer, both modules define it identically" convention
         !! run_t%opt_algorithm_kind's own comment gives.
      integer, parameter :: SADDLE_DIMER = 1
         !! trc_dlfind_bridge's own TRC_OPT_SADDLE_DIMER, same convention.
      character(len=24), parameter :: GLOBAL_OPTS(*) = [character(len=24) :: &
         "opt-global-pop-size", "opt-global-maxcycle", "opt-global-init-pop-size", &
         "opt-global-nsave", "opt-global-distribution", "opt-global-reset", "opt-global-radius", &
         "opt-global-contraction", "opt-global-tolerance-r", "opt-global-tolerance-g", &
         "opt-global-mutation-rate", "opt-global-death-rate"]
      logical :: seed_ok
      integer :: i

      seed_ok = r%opt_is_global .or. r%opt_saddle_method_kind == SADDLE_DIMER
      if (r%cli%is_set("opt-seed") .and. .not. seed_ok) then
         call run_die(r, "--opt-seed needs --opt-algorithm stochastic/genetic or --opt-saddle-method dimer")
      end if

      if (.not. r%opt_is_global) then
         do i = 1, size(GLOBAL_OPTS)
            if (r%cli%is_set(trim(GLOBAL_OPTS(i)))) call run_die(r, "--"//trim(GLOBAL_OPTS(i))// &
                     " needs --opt-algorithm stochastic or genetic")
         end do
         if (r%cli%is_set("opt-global-scalefac")) call run_die(r, "--opt-global-scalefac needs "// &
                  "--opt-algorithm stochastic or genetic")
         if (r%cli%is_set("opt-seed")) then
            call r%cli%get("opt-seed", r%opt_seed, r%err)
            call run_check(r, "--opt-seed")
         end if
         return
      end if

      if (r%cli%is_set("opt-seed")) then
         call r%cli%get("opt-seed", r%opt_seed, r%err)
         call run_check(r, "--opt-seed")
      end if
      if (r%cli%is_set("opt-global-pop-size")) then
         call r%cli%get("opt-global-pop-size", r%opt_global_pop_size, r%err)
         call run_check(r, "--opt-global-pop-size")
         if (r%opt_global_pop_size < 1) call run_die(r, "--opt-global-pop-size must be at least 1")
      end if
      if (r%cli%is_set("opt-global-maxcycle")) then
         call r%cli%get("opt-global-maxcycle", r%opt_global_maxcycle, r%err)
         call run_check(r, "--opt-global-maxcycle")
         if (r%opt_global_maxcycle < 1) call run_die(r, "--opt-global-maxcycle must be at least 1")
      end if
      if (r%cli%is_set("opt-global-init-pop-size")) then
         call r%cli%get("opt-global-init-pop-size", r%opt_global_init_pop_size, r%err)
         call run_check(r, "--opt-global-init-pop-size")
         if (r%opt_global_init_pop_size < 1) call run_die(r, "--opt-global-init-pop-size must be at least 1")
         if (r%opt_algorithm_kind /= ALG_GENETIC) call run_die(r, "--opt-global-init-pop-size "// &
                  "needs --opt-algorithm genetic")
      end if
      if (r%cli%is_set("opt-global-nsave")) then
         call r%cli%get("opt-global-nsave", r%opt_global_nsave, r%err)
         call run_check(r, "--opt-global-nsave")
         if (r%opt_global_nsave < 1) call run_die(r, "--opt-global-nsave must be at least 1")
      end if
      if (r%cli%is_set("opt-global-distribution")) then
         call r%cli%get("opt-global-distribution", r%opt_global_distribution, r%err)
         call run_check(r, "--opt-global-distribution")
         if (r%opt_global_distribution < 1 .or. r%opt_global_distribution > 3) then
            call run_die(r, "--opt-global-distribution must be 1 (box), 2 (sphere) or 3 (Gaussian)")
         end if
         if (r%opt_algorithm_kind /= ALG_STOCHASTIC) call run_die(r, "--opt-global-distribution "// &
                  "needs --opt-algorithm stochastic")
      end if
      if (r%cli%is_set("opt-global-reset")) then
         call r%cli%get("opt-global-reset", r%opt_global_reset, r%err)
         call run_check(r, "--opt-global-reset")
         if (r%opt_global_reset < 1) call run_die(r, "--opt-global-reset must be at least 1")
         if (r%opt_algorithm_kind /= ALG_GENETIC) call run_die(r, "--opt-global-reset "// &
                  "needs --opt-algorithm genetic")
      end if
      if (r%cli%is_set("opt-global-radius")) then
         call r%cli%get("opt-global-radius", r%opt_global_radius, r%err)
         call run_check(r, "--opt-global-radius")
         if (r%opt_global_radius <= 0.0_dp) call run_die(r, "--opt-global-radius must be positive")
      end if
      if (r%cli%is_set("opt-global-contraction")) then
         call r%cli%get("opt-global-contraction", r%opt_global_contraction, r%err)
         call run_check(r, "--opt-global-contraction")
         if (r%opt_global_contraction <= 0.0_dp .or. r%opt_global_contraction >= 1.0_dp) then
            call run_die(r, "--opt-global-contraction must be in (0, 1)")
         end if
         if (r%opt_algorithm_kind /= ALG_STOCHASTIC) call run_die(r, "--opt-global-contraction "// &
                  "needs --opt-algorithm stochastic")
      end if
      if (r%cli%is_set("opt-global-tolerance-r")) then
         call r%cli%get("opt-global-tolerance-r", r%opt_global_tolerance_r, r%err)
         call run_check(r, "--opt-global-tolerance-r")
         if (r%opt_global_tolerance_r <= 0.0_dp) call run_die(r, "--opt-global-tolerance-r must be positive")
         if (r%opt_algorithm_kind /= ALG_STOCHASTIC) call run_die(r, "--opt-global-tolerance-r "// &
                  "needs --opt-algorithm stochastic")
      end if
      if (r%cli%is_set("opt-global-tolerance-g")) then
         call r%cli%get("opt-global-tolerance-g", r%opt_global_tolerance_g, r%err)
         call run_check(r, "--opt-global-tolerance-g")
         if (r%opt_global_tolerance_g <= 0.0_dp) call run_die(r, "--opt-global-tolerance-g must be positive")
      end if
      if (r%cli%is_set("opt-global-mutation-rate")) then
         call r%cli%get("opt-global-mutation-rate", r%opt_global_mutation_rate, r%err)
         call run_check(r, "--opt-global-mutation-rate")
         if (r%opt_global_mutation_rate <= 0.0_dp .or. r%opt_global_mutation_rate >= 1.0_dp) then
            call run_die(r, "--opt-global-mutation-rate must be in (0, 1)")
         end if
         if (r%opt_algorithm_kind /= ALG_GENETIC) call run_die(r, "--opt-global-mutation-rate "// &
                  "needs --opt-algorithm genetic")
      end if
      if (r%cli%is_set("opt-global-death-rate")) then
         call r%cli%get("opt-global-death-rate", r%opt_global_death_rate, r%err)
         call run_check(r, "--opt-global-death-rate")
         if (r%opt_global_death_rate <= 0.0_dp .or. r%opt_global_death_rate >= 1.0_dp) then
            call run_die(r, "--opt-global-death-rate must be in (0, 1)")
         end if
         if (r%opt_algorithm_kind /= ALG_GENETIC) call run_die(r, "--opt-global-death-rate "// &
                  "needs --opt-algorithm genetic")
      end if
      if (r%cli%is_set("opt-global-scalefac")) then
         call r%cli%get("opt-global-scalefac", r%opt_global_scalefac, r%err)
         call run_check(r, "--opt-global-scalefac")
         if (r%opt_global_scalefac <= 0.0_dp) call run_die(r, "--opt-global-scalefac must be positive")
         if (r%opt_algorithm_kind /= ALG_STOCHASTIC) call run_die(r, "--opt-global-scalefac "// &
                  "needs --opt-algorithm stochastic")
      end if
   end subroutine run_parse_opt_global_options

   !
   ! --opt-frozen "1-3,7": a comma-separated list of 1-based atom indices
   ! and ranges, into a sorted array with no duplicates. Validated here
   ! (range, at least one atom) so an out-of-range index is refused before
   ! trc_dlfind_optimize ever runs, the same convention every other --opt-*
   ! validation in this file follows.
   !
   module subroutine run_parse_frozen_list(spec, natm, frozen, err)
      character(len=*), intent(in) :: spec
      integer, intent(in) :: natm
      integer, allocatable, intent(out) :: frozen(:)
      type(error_t), intent(inout) :: err
      type(string_type), allocatable :: toks(:), bounds(:)
      type(error_t) :: perr
      integer :: i, j, lo, hi, cnt
      character(len=:), allocatable :: tok
      logical :: seen(natm)
      logical :: bad, have_lo, have_hi

      seen = .false.
      cnt = 0
      toks = split(strip(spec), ",")
      do i = 1, size(toks)
         tok = strip(char(toks(i)))
         if (len(tok) == 0) cycle
         bad = .false.
         have_lo = .false.
         have_hi = .false.
         if (find(tok, "-") > 0) then
            bounds = split(tok, "-")
            if (size(bounds) == 2) then
               call parse_int(bounds(1), lo, perr)
               have_lo = .not. perr%has_error()
               call parse_int(bounds(2), hi, perr)
               have_hi = .not. perr%has_error()
            end if
            bad = .not. (have_lo .and. have_hi)
         else
            call parse_int(tok, lo, perr)
            bad = perr%has_error()
            hi = lo
         end if
         if (bad .or. lo < 1 .or. hi < lo .or. hi > natm) then
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-frozen: '"//tok// &
                         "' is not a valid 1-based atom index or range for "//to_string(natm)//" atoms")
            return
         end if
         do j = lo, hi
            if (.not. seen(j)) then
               seen(j) = .true.
               cnt = cnt + 1
            end if
         end do
      end do
      if (cnt == 0) then
         call err%set(int(ERROR_VALIDATION, default_int), "--opt-frozen: no atom indices found in '"//trim(spec)//"'")
         return
      end if
      allocate (frozen(cnt))
      cnt = 0
      do i = 1, natm
         if (seen(i)) then
            cnt = cnt + 1
            frozen(cnt) = i
         end if
      end do
   end subroutine run_parse_frozen_list

   !
   ! --opt-constrain "angle 1 2 3; bond 2 4": a ';'-separated list (pic_cli's
   ! own options are single-valued -- its header comment, "the last
   ! occurrence wins" -- so a repeatable CLI flag is not how this bridge
   ! takes more than one constraint). Each entry is "<type> <atoms...>":
   ! bond i j; angle i j k (j is the vertex); torsion i j k l (j-k is the
   ! central bond); cartesian i (expanded to all three components by
   ! trc_dlfind_bridge%build_icons); bond-difference i j k (bond(i,j) -
   ! bond(j,k) -- VERIFIED from libdlfind's dlf_hdlc_constraint.f90%
   ! assign_cons directly: this type uses three atoms, not four). The
   ! parallel arrays returned are trc_dlfind_bridge's own TRC_OPT_CONS_*
   ! values (1..5) by construction -- see run_t%opt_cons_kind's own comment
   ! for why this module cannot just import that enum.
   !
   module subroutine run_parse_constraints(spec, natm, cons_kind, cons_a, cons_b, cons_c, cons_d, err)
      character(len=*), intent(in) :: spec
      integer, intent(in) :: natm
      integer, allocatable, intent(out) :: cons_kind(:), cons_a(:), cons_b(:), cons_c(:), cons_d(:)
      type(error_t), intent(inout) :: err
      integer, parameter :: MAX_CONSTRAINTS = 64   !! generous fixed cap, checked below
      integer, parameter :: N_CONS_KINDS = 5   !! TRC_OPT_CONS_{BOND,ANGLE,TORSION,CARTESIAN,BOND_DIFF}
      integer, parameter :: NSLOTS(N_CONS_KINDS) = [2, 3, 4, 1, 3]
         !! atom slots TRC_OPT_CONS_{BOND,ANGLE,TORSION,CARTESIAN,BOND_DIFF}
         !! (1..5, found_kind's own values) actually use -- the rest of
         !! atom_index(:) stays the 0 this subroutine's own reset put there.
      character(len=:), allocatable :: entry, kind_tok
      type(string_type), allocatable :: entries(:), words(:)
      type(error_t) :: perr
      integer :: i, ncons, k, k2, nslots_here
      integer :: found_kind(MAX_CONSTRAINTS), found_a(MAX_CONSTRAINTS), found_b(MAX_CONSTRAINTS)
      integer :: found_c(MAX_CONSTRAINTS), found_d(MAX_CONSTRAINTS)
      integer :: atom_index(4)
      logical :: bad

      ncons = 0
      entries = split(strip(spec), ";")
      do i = 1, size(entries)
         entry = strip(char(entries(i)))
         if (len(entry) == 0) cycle
         ncons = ncons + 1
         if (ncons > size(found_kind)) then
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: more than "// &
                         to_string(size(found_kind))//" constraints in one run is not supported")
            return
         end if
         words = tokenize(entry)
         if (size(words) < 1) then
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: '"//entry// &
                         "' has no atom indices")
            return
         end if
         call run_to_lower_into(char(words(1)), kind_tok)
         select case (kind_tok)
         case ("bond")
            found_kind(ncons) = 1
         case ("angle")
            found_kind(ncons) = 2
         case ("torsion")
            found_kind(ncons) = 3
         case ("cartesian")
            found_kind(ncons) = 4
         case ("bond-difference")
            found_kind(ncons) = 5
         case default
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: unknown type '"//kind_tok// &
                         "': use bond, angle, torsion, cartesian or bond-difference")
            return
         end select
         nslots_here = NSLOTS(found_kind(ncons))
         atom_index = 0
         bad = (size(words) - 1 < nslots_here)
         if (.not. bad) then
            do k = 1, nslots_here
               call parse_int(words(k + 1), atom_index(k), perr)
               bad = bad .or. perr%has_error()
            end do
         end if
         if (bad) then
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: '"//entry// &
                         "' does not have the atom indices '"//kind_tok//"' needs")
            return
         end if
         found_a(ncons) = atom_index(1)
         found_b(ncons) = atom_index(2)
         found_c(ncons) = atom_index(3)
         found_d(ncons) = atom_index(4)
         ! Only the first NSLOTS(kind) entries of atom_index are atoms THIS type
         ! uses; the rest are the 0 this subroutine's own `atom_index = 0` reset
         ! put there and must stay 0 (never a slot a caller filled in).
         ! Checking atom_index(k) /= 0 here (the bug this replaced) let a TYPED 0
         ! through as if it meant "unused": DL-FIND's own assign_cons then
         ! reports "crosses residue boundary" and `error stop`s, since
         ! atom 0 belongs to no residue at all.
         do k = 1, NSLOTS(found_kind(ncons))
            if (atom_index(k) < 1 .or. atom_index(k) > natm) then
               call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: '"//entry// &
                            "': atom index out of range for "//to_string(natm)//" atoms")
               return
            end if
         end do
         do k = 1, NSLOTS(found_kind(ncons))
            do k2 = k + 1, NSLOTS(found_kind(ncons))
               if (atom_index(k) == atom_index(k2)) then
                  call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: '"//entry// &
                               "': the same atom appears twice in one constraint")
                  return
               end if
            end do
         end do
      end do
      if (ncons == 0) then
         call err%set(int(ERROR_VALIDATION, default_int), "--opt-constrain: no constraints found in '"// &
                      trim(spec)//"'")
         return
      end if
      allocate (cons_kind(ncons), source=found_kind(1:ncons))
      allocate (cons_a(ncons), source=found_a(1:ncons))
      allocate (cons_b(ncons), source=found_b(1:ncons))
      allocate (cons_c(ncons), source=found_c(1:ncons))
      allocate (cons_d(ncons), source=found_d(1:ncons))
   end subroutine run_parse_constraints

   !
   ! Checks every --opt-constrain entry against the HDLC residue
   ! partition and the --opt-frozen list, both only known once the
   ! geometry-dependent residues are built (trc_constants%trc_covalent_
   ! radius-based, build_opt_residues) -- run_parse_constraints above only
   ! has `natm`, so a syntactically valid entry can still name two atoms
   ! DL-FIND's own assign_cons would refuse once it ran. `die`s on the
   ! first bad entry, the same convention every other --opt-* validation
   ! in this file follows.
   !
   module subroutine run_check_opt_constrain_residues(r, opt_residues)
      type(run_t), intent(in) :: r
      integer, intent(in) :: opt_residues(:)
      integer, parameter :: N_CONS_KINDS = 5   !! TRC_OPT_CONS_{BOND,ANGLE,TORSION,CARTESIAN,BOND_DIFF}
      integer, parameter :: NSLOTS(N_CONS_KINDS) = [2, 3, 4, 1, 3]
         !! run_parse_constraints' own NSLOTS -- which of atoms(1:4) a given
         !! TRC_OPT_CONS_* kind actually uses.
      integer :: ic, k, na, res1
      integer :: atoms(4)

      do ic = 1, size(r%opt_cons_kind)
         na = NSLOTS(r%opt_cons_kind(ic))
         atoms = [r%opt_cons_a(ic), r%opt_cons_b(ic), r%opt_cons_c(ic), r%opt_cons_d(ic)]
         do k = 1, na
            if (allocated(r%opt_frozen_atoms)) then
               if (any(r%opt_frozen_atoms == atoms(k))) then
                  call run_die(r, "--opt-constrain: atom "//to_string(atoms(k))//" is frozen (--opt-frozen); "// &
                           "a constraint cannot touch a frozen atom")
               end if
            end if
         end do
         res1 = opt_residues(atoms(1))
         do k = 2, na
            if (opt_residues(atoms(k)) /= res1) then
               call run_die(r, "--opt-constrain: atoms "//to_string(atoms(1))//" and "//to_string(atoms(k))// &
                        " are in different HDLC residues; DL-FIND cannot constrain across a residue boundary")
            end if
         end do
      end do
   end subroutine run_check_opt_constrain_residues

   !
   ! --opt-neb <endpoint.xyz>: the product geometry must have the SAME
   ! atoms in the SAME order as the start geometry -- DL-FIND's own
   ! coords2 is a flat Cartesian array with no atom labels of its own, so a
   ! silently reordered or substituted atom would just move the wrong
   ! nucleus. Checked element-by-element rather than just atom count.
   !
   module subroutine run_check_neb_endpoint(natm, z, z2, natm2, err, what)
      integer, intent(in) :: natm, natm2
      integer, intent(in) :: z(natm), z2(natm2)
      type(error_t), intent(inout) :: err
      character(len=*), intent(in), optional :: what
         !! What to call the file that failed in the message -- "the
         !! endpoint geometry" (--opt-neb's own default) or "--opt-neb-via
         !! frame N" (run_read_xyz_frames' own caller, below).
      character(len=:), allocatable :: label
      integer :: i
      label = "the endpoint geometry"
      if (present(what)) label = what
      if (natm2 /= natm) then
         call err%set(int(ERROR_VALIDATION, default_int), "--opt-neb: "//label//" has "// &
                      to_string(natm2)//" atoms, the start geometry has "//to_string(natm))
         return
      end if
      do i = 1, natm
         if (z2(i) /= z(i)) then
            call err%set(int(ERROR_VALIDATION, default_int), "--opt-neb: atom "//to_string(i)// &
                         " differs between the start geometry and "//label//" (atom order must match)")
            return
         end if
      end do
   end subroutine run_check_neb_endpoint

   !
   ! --opt-neb-via <file.xyz>: one or more xyz FRAMES, standard multi-frame
   ! xyz (each frame is its own "natm / title / natm atom lines" block,
   ! concatenated -- the same convention write_traj_xyz already writes),
   ! read until end of file. Each frame's own atom count/order is checked
   ! against the start geometry by the caller (run_check_neb_endpoint, above),
   ! not here -- this only parses.
   !
   module subroutine run_read_xyz_frames(path, natm, z, nframes, xyz, err)
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
         !! the start geometry's own atomic numbers -- every frame is
         !! checked against this directly (run_check_neb_endpoint's own
         !! element-by-element comparison), same atoms, same order.
      integer, intent(out) :: nframes
      real(dp), allocatable, intent(out) :: xyz(:, :, :)   !! (3, natm, nframes), Bohr
      type(error_t), intent(inout) :: err
      integer :: u, ios, ia, natm_f
      character(len=1024) :: line
      character(len=16) :: tok
      integer :: z_f(natm)
      real(dp) :: r(3)
      real(dp), allocatable :: frames_read(:, :, :)

      nframes = 0
      open (newunit=u, file=path, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         call err%set(int(ERROR_IO, default_int), "--opt-neb-via: cannot open '"//trim(path)//"'")
         return
      end if
      allocate (frames_read(3, natm, 16))   ! grown below if more frames turn up
      do
         read (u, *, iostat=ios) natm_f
         if (ios /= 0) exit   ! clean end of file: done
         if (natm_f /= natm) then
            call err%set(int(ERROR_PARSE, default_int), "--opt-neb-via: '"//trim(path)//"': frame "// &
                         to_string(nframes + 1)//" has "//to_string(natm_f)//" atoms, the start "// &
                         "geometry has "//to_string(natm))
            close (u)
            return
         end if
         read (u, "(a)", iostat=ios) line   ! title line, discarded
         if (ios /= 0) then
            call err%set(int(ERROR_PARSE, default_int), "--opt-neb-via: '"//trim(path)//"': frame "// &
                         to_string(nframes + 1)//" ends after its atom count, before its title line")
            close (u)
            return
         end if
         nframes = nframes + 1
         if (nframes > size(frames_read, 3)) then
            block
               real(dp), allocatable :: grown(:, :, :)
               allocate (grown(3, natm, 2*size(frames_read, 3)))
               grown(:, :, 1:size(frames_read, 3)) = frames_read
               call move_alloc(grown, frames_read)
            end block
         end if
         do ia = 1, natm
            read (u, "(a)", iostat=ios) line
            if (ios /= 0) then
               call err%set(int(ERROR_PARSE, default_int), "--opt-neb-via: '"//trim(path)//"': frame "// &
                            to_string(nframes)//" ends after "//to_string(ia - 1)//" of "// &
                            to_string(natm)//" atoms")
               close (u)
               return
            end if
            read (line, *, iostat=ios) tok, r
            if (ios /= 0) then
               call err%set(int(ERROR_PARSE, default_int), "--opt-neb-via: '"//trim(path)//"': frame "// &
                            to_string(nframes)//", atom "//to_string(ia)//" is not 'symbol x y z': '"// &
                            trim(line)//"'")
               close (u)
               return
            end if
            z_f(ia) = app_atomic_number(tok)
            frames_read(:, ia, nframes) = r*BOHR_PER_ANGSTROM
         end do
         call run_check_neb_endpoint(natm, z, z_f, natm, err, what="--opt-neb-via frame "//to_string(nframes))
         if (err%has_error()) then
            close (u)
            return
         end if
      end do
      close (u)
      if (nframes == 0) then
         call err%set(int(ERROR_PARSE, default_int), "--opt-neb-via: '"//trim(path)//"' has no frames")
         return
      end if
      allocate (xyz(3, natm, nframes), source=frames_read(:, :, 1:nframes))
   end subroutine run_read_xyz_frames

   !
   ! `<geometry base>_neb.xyz`: the converged NEB path, one frame per image
   ! (both endpoints included), Angstrom, comment line carries the image
   ! index and energy, with the climbing/highest-energy image marked, so
   ! the whole path can be inspected or animated without re-running DL-FIND.
   !
   module subroutine run_write_neb_path_xyz(r, path)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      integer :: u, ia, k
      character(len=:), allocatable :: mark
      open (newunit=u, file=path, status="replace", action="write")
      do k = 1, r%opt_neb_n_images
         write (u, "(i0)") r%natm
         mark = ""
         if (k == r%opt_neb_climb_image) mark = " <- climbing/highest-energy image"
         write (u, "(a,i0,a,i0,a,f20.12,a,a)") "image ", k, " of ", r%opt_neb_n_images, "   E = ", &
            r%opt_neb_image_energy(k), " Eh", mark
         do ia = 1, r%natm
            write (u, "(a2,3f20.12)") app_element_symbol(r%z(ia)), r%opt_neb_image_xyz(:, ia, k)/BOHR_PER_ANGSTROM
         end do
      end do
      close (u)
   end subroutine run_write_neb_path_xyz

   !
   ! The optimisation's own final summary: converged or not, steps,
   ! energy evaluations and the final energy.
   ! Per-step lines (energy, max/rms |g|, SCF iterations, wall time) are
   ! logged from trc_dlfind_bridge%cb_put_coords itself as each step is
   ! accepted, not replayed here.
   !
   module subroutine run_opt_summary(r)
      type(run_t), intent(in) :: r
      integer :: n_imag
      ! The step table's own closing rule (app/trc_dlfind_bridge.F90's
      ! consider_step prints the header/rows; this is the one place that
      ! knows the whole run is over) -- skipped for --opt-neb, which
      ! never calls consider_step at all (cb_put_coords's own NEB branch:
      ! no header was ever printed, so no rule belongs here either).
      if (r%opt_steps > 0 .and. .not. r%opt_neb_have) call global_logger%info(TRC_OPT_STEP_RULE)
      call global_logger%info("")
      call global_logger%info("   geometry optimisation (--opt)")
      call global_logger%info("   ----------------------------------------------------------")
      block
         character(len=:), allocatable :: logline
         logline = "   target: "//opt_target_name(r%opt_target_kind)
         if (r%opt_target_kind == 1) then
            logline = logline//"   saddle method: "//opt_saddle_method_name(r%opt_saddle_method_kind)
         else
            logline = logline//"   algorithm: "//opt_algorithm_name(r%opt_algorithm_kind)
         end if
         if (r%opt_connect) logline = logline//"   (connect)"
         call global_logger%info(logline)
      end block
      block
         character(len=1024) :: logline
         write (logline, "(a,a,a,i0,a,i0,a,a,a,es10.3,a)") "   ", &
            trim(merge("converged    ", "NOT converged", r%opt_converged)), &
            "   steps: ", r%opt_steps, "   energy evaluations: ", r%opt_evals, &
            "   coordinates: ", trim(opt_coords_name(r%opt_coords_kind)), &
            "   final max|g|: ", r%opt_gmax, " Eh/Bohr"
         call global_logger%info(trim(logline))
      end block
      if (r%opt_hess_evals > 0) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a)") "   terco analytic-Hessian callback (cb_get_hessian): called ", &
               r%opt_hess_evals, " time(s)"
            call global_logger%info(trim(logline))
         end block
      end if
      block
         character(len=1024) :: logline
         write (logline, "(a,f22.12,a)") "   final energy     : ", r%etot, " Eh"
         call global_logger%info(trim(logline))
      end block
      ! The TS check: for a saddle target, state how many imaginary modes
      ! the final geometry has -- exactly one is a genuine first-order
      ! saddle point. Only available together with --hessian (which runs
      ! AFTER --opt, at the geometry --opt reports -- trc_main's own
      ! ordering), since that is what computes r%freq_cm1 at all.
      if (r%opt_target_kind == 1 .and. r%hess_ok) then
         n_imag = count(r%freq_cm1 < 0.0_dp)
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a)") "   saddle run_check: ", n_imag, " imaginary mode(s) (expect exactly 1 "// &
               "for a transition state)"
            call global_logger%info(trim(logline))
            if (n_imag >= 1) then
               write (logline, "(a,f12.2,a)") "   imaginary frequency (lowest): ", minval(r%freq_cm1), " i cm^-1"
               call global_logger%info(trim(logline))
            end if
         end block
      end if
      if (r%opt_connect) then
         if (r%opt_have_ts) then
            block
               character(len=1024) :: logline
               write (logline, "(a,f22.12,a,a,a,es10.3,a)") "   transition state : ", r%opt_ts_energy, " Eh   ", &
                  trim(merge("converged    ", "NOT converged", r%opt_ts_converged)), "   max|g|: ", r%opt_ts_gmax, &
                  " Eh/Bohr"
               call global_logger%info(trim(logline))
            end block
         end if
         if (r%opt_have_min_minus) then
            block
               character(len=1024) :: logline
               write (logline, "(a,f22.12,a,a,a,es10.3,a)") "   minimum (-distort): ", r%opt_min_minus_energy, " Eh   ", &
                  trim(merge("converged    ", "NOT converged", r%opt_min_minus_converged)), "   max|g|: ", &
                  r%opt_min_minus_gmax, " Eh/Bohr"
               call global_logger%info(trim(logline))
            end block
         end if
         if (r%opt_have_min_plus) then
            block
               character(len=1024) :: logline
               write (logline, "(a,f22.12,a,a,a,es10.3,a)") "   minimum (+distort): ", r%opt_min_plus_energy, " Eh   ", &
                  trim(merge("converged    ", "NOT converged", r%opt_min_plus_converged)), "   max|g|: ", &
                  r%opt_min_plus_gmax, " Eh/Bohr"
               call global_logger%info(trim(logline))
            end block
         end if
      end if
      if (r%opt_neb_have) then
         block
            character(len=1024) :: logline
            character(len=:), allocatable :: climb_mark, image_mark
            integer :: k
            climb_mark = ""
            if (r%opt_neb_climb) climb_mark = "   (climbing image)"
            write (logline, "(a,i0,a,a,a)") "   NEB: ", r%opt_neb_n_images, " images, endpoints "// &
               trim(r%opt_neb_endpoints_arg), climb_mark
            call global_logger%info(trim(logline))
            do k = 1, r%opt_neb_n_images
               image_mark = ""
               if (k == r%opt_neb_climb_image) image_mark = " <- highest energy"
               write (logline, "(a,i0,a,f22.12,a,a)") "     image ", k, "  E = ", r%opt_neb_image_energy(k), &
                  " Eh", image_mark
               call global_logger%info(trim(logline))
            end do
         end block
      end if
   end subroutine run_opt_summary

   pure module function opt_coords_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name
      select case (kind)
      case (1)
         name = "hdlc"
      case (2)
         name = "hdlc-tc"
      case (3)
         name = "dlc"
      case (4)
         name = "dlc-tc"
      case default
         name = "cartesian"
      end select
   end function opt_coords_name

   pure module function opt_target_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name
      if (kind == 1) then
         name = "saddle"
      else
         name = "minimum"
      end if
   end function opt_target_name

   pure module function opt_saddle_method_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name
      if (kind == 1) then
         name = "dimer"
      else
         name = "prfo"
      end if
   end function opt_saddle_method_name

   pure module function opt_algorithm_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name
      select case (kind)
      case (0)
         name = "sd"
      case (1)
         name = "cg-auto"
      case (2)
         name = "cg"
      case (10)
         name = "prfo"
      case (20)
         name = "nr"
      case (30)
         name = "damped"
      case (51)
         name = "stochastic"
      case (52)
         name = "genetic"
      case default
         name = "lbfgs"
      end select
   end function opt_algorithm_name

   !
   ! `<geometry base>_opt.xyz`: the optimised geometry in Angstrom, same
   ! atom order as the input. Written whether or not the run converged
   ! (the comment line says which), so a failed optimisation's last
   ! geometry is never silently lost.
   !
   module subroutine run_write_opt_xyz(r, path, xyz_bohr)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: xyz_bohr(3, r%natm)
      call app_write_opt_xyz_frame(path, r%natm, r%z, real(r%etot, dp), xyz_bohr, &
                               trim(merge("converged    ", "NOT converged", r%opt_converged)))
   end subroutine run_write_opt_xyz

   !
   ! `<geometry base>_ts.xyz`/`_minus.xyz`/`_plus.xyz`: --opt-connect's own
   ! three structures (the transition state and its two downhill minima),
   ! one labelled energy each -- write_opt_xyz's own format, but with an
   ! energy that is NOT r%etot (the context is re-converged at the TS,
   ! final_coords, after --opt returns; the two minima's own energies came
   ! straight from DL-FIND's cb_put_coords and were never re-run).
   !
   module subroutine run_write_labelled_xyz(r, path, xyz_bohr, e, label)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: xyz_bohr(3, r%natm)
      real(dp), intent(in) :: e
      character(len=*), intent(in) :: label
      integer :: u, ia
      open (newunit=u, file=path, status="replace", action="write")
      write (u, "(i0)") r%natm
      write (u, "(a,a,f0.12)") label, "   E = ", e
      do ia = 1, r%natm
         write (u, "(a2,3f20.12)") app_element_symbol(r%z(ia)), xyz_bohr(:, ia)/BOHR_PER_ANGSTROM
      end do
      close (u)
   end subroutine run_write_labelled_xyz

   !
   ! `<geometry base>_traj.xyz`: one frame per ACCEPTED step (never a
   ! line-search trial -- trc_dlfind_bridge%cb_put_coords's own
   ! distinction), Angstrom, comment line carries the step index and
   ! energy.
   !
   module subroutine run_write_traj_xyz(r, path, traj_xyz_bohr, traj_energy)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: traj_xyz_bohr(3, r%natm, r%opt_n_traj)
      real(dp), intent(in) :: traj_energy(r%opt_n_traj)
      call app_write_traj_xyz_frames(path, r%natm, r%z, r%opt_n_traj, traj_xyz_bohr, traj_energy)
   end subroutine run_write_traj_xyz

end submodule trc_run_opt
