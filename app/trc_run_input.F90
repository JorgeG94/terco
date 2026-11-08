module trc_run_input
   !! The TOML input file (docs/INPUT.md): `terco run.toml [flags]`.
   !!
   !! A file is turned into the command-line arguments it stands for and
   !! parsed by the same `cli_t` as a command line, so every check `main`
   !! makes on a flag applies to its key unchanged; there is no second
   !! validation path to drift. The arguments from the file come first and the real
   !! command line after them, and `cli_t` keeps the last value of an option,
   !! which is what makes flags override the file.
   !!
   !! `RUN_INPUT_KEYS` is the whole format: one row per key, naming its flag
   !! and how its value becomes the argument of that flag. docs/INPUT.md is its
   !! documentation; keep the two in step.
   !!
   !! The `"input"` block of `--json` is written from the PARSED command line
   !! (`run_input_echo`), not from the file: it shows what the run used,
   !! overrides and declared defaults included. A key whose flag has no
   !! declared default and was not given (the charge, the guess, ...) is left
   !! out, because its value is decided later from the system.
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_ascii, only: to_lower
   use pic_string_type, only: string_type, assignment(=)
   use pic_cli, only: cli_t
   use json_module, only: json_core, json_value
   use tomlf, only: toml_table, toml_array, toml_key, toml_error, toml_context, toml_load, get_value, len
   implicit none
   private

   public :: run_input_is_toml, run_input_args, run_input_echo

   integer, parameter :: KIND_STRING = 1
      !! A string, passed through.
   integer, parameter :: KIND_INT = 2
   integer, parameter :: KIND_FLOAT = 3
      !! A float; an integer is accepted too (`epsilon = 78`).
   integer, parameter :: KIND_BOOL = 4
      !! A switch: true gives the flag, false gives nothing.
   integer, parameter :: KIND_PATH = 5
      !! A file path, relative to the directory of the TOML file.
   integer, parameter :: KIND_FILE_OR_NAME = 6
      !! A basis: a path relative to the TOML file if that file exists, else
      !! a name looked up as on the command line.
   integer, parameter :: KIND_INT_LIST = 7
      !! An integer array, joined with commas (`--opt-frozen`).
   integer, parameter :: KIND_STRING_LIST = 8
      !! A string array, joined with semicolons (`--opt-constrain`).
   integer, parameter :: KIND_METHOD = 9
      !! `[correlation] method`: "rimp2" gives `--rimp2`, "none" nothing.
   integer, parameter :: KIND_BOOL_ON = 10
      !! A switch that is on by default: true gives the flag, false gives
      !! its `--no-` form (`incremental = false` is `--no-incremental`).
   integer, parameter :: KIND_BOOL_TRI = 11
      !! A switch with three states: true gives the flag, false gives its
      !! `--no-` form, and absent means neither, because the default is
      !! something the code decides (`symmetry_skeleton` and its size rule).
      !! The `--json` echo leaves it out when neither flag was given, so
      !! writing the echo back as a TOML file keeps that third state.

   type :: run_input_key_t
      !! One key of the format.
      character(len=16) :: section
         !! The table it lives in, dotted for a sub-table (`opt.neb`).
      character(len=20) :: key
      character(len=26) :: flag
         !! Its command-line option, without the dashes; blank for the
         !! geometry, which is the positional argument.
      integer :: kind
   end type run_input_key_t

   integer, parameter :: NKEYS = 91
   type(run_input_key_t), parameter :: RUN_INPUT_KEYS(NKEYS) = [ &
      run_input_key_t("molecule", "geometry", "", KIND_PATH), &
      run_input_key_t("molecule", "charge", "chrg", KIND_INT), &
      run_input_key_t("molecule", "unpaired", "uhf", KIND_INT), &
      run_input_key_t("basis", "main", "basis", KIND_FILE_OR_NAME), &
      run_input_key_t("basis", "auxiliary", "aux", KIND_FILE_OR_NAME), &
      run_input_key_t("basis", "jk_auxiliary", "jkfit", KIND_FILE_OR_NAME), &
      run_input_key_t("basis", "cartesian", "cart", KIND_BOOL), &
      run_input_key_t("basis", "linear_dependence", "lindep", KIND_FLOAT), &
      run_input_key_t("scf", "functional", "func", KIND_STRING), &
      run_input_key_t("scf", "guess", "guess", KIND_STRING), &
      run_input_key_t("scf", "max_iterations", "iterations", KIND_INT), &
      run_input_key_t("scf", "energy_tol", "econv", KIND_FLOAT), &
      run_input_key_t("scf", "commutator_tol", "gconv", KIND_FLOAT), &
      run_input_key_t("scf", "screening", "thresh", KIND_FLOAT), &
      run_input_key_t("scf", "incremental", "incremental", KIND_BOOL_ON), &
      run_input_key_t("scf", "density_fitted", "df", KIND_BOOL), &
      run_input_key_t("scf", "restart", "restart", KIND_PATH), &
      run_input_key_t("scf", "symmetry", "symmetry", KIND_STRING), &
      run_input_key_t("scf", "symmetry_tol", "symmetry-tol", KIND_FLOAT), &
      run_input_key_t("scf", "irrep_nelec", "irrep-nelec", KIND_STRING), &
      run_input_key_t("scf", "symmetry_skeleton", "sym-skeleton", KIND_BOOL_TRI), &
      run_input_key_t("scf.convergence", "diis_size", "diis", KIND_INT), &
      run_input_key_t("scf.convergence", "adiis", "adiis", KIND_BOOL), &
      run_input_key_t("scf.convergence", "level_shift", "level-shift", KIND_FLOAT), &
      run_input_key_t("scf.convergence", "damping", "damp", KIND_FLOAT), &
      run_input_key_t("scf.convergence", "damping_steps", "damp-steps", KIND_INT), &
      run_input_key_t("scf.convergence", "diis_start", "diis-start", KIND_FLOAT), &
      run_input_key_t("scf.convergence", "soscf", "soscf", KIND_BOOL), &
      run_input_key_t("scf.convergence", "soscf_start", "soscf-start", KIND_FLOAT), &
      run_input_key_t("dft", "grid", "grid", KIND_INT), &
      run_input_key_t("dft", "nlc_grid", "nlc-grid", KIND_INT), &
      run_input_key_t("solvation", "model", "pcm", KIND_STRING), &
      run_input_key_t("solvation", "epsilon", "eps", KIND_FLOAT), &
      run_input_key_t("correlation", "method", "rimp2", KIND_METHOD), &
      run_input_key_t("correlation", "frozen_core", "frozen", KIND_INT), &
      run_input_key_t("correlation", "aux_block", "aux-block", KIND_INT), &
      run_input_key_t("driver", "gradient", "grad", KIND_BOOL), &
      run_input_key_t("driver", "hessian", "hessian", KIND_BOOL), &
      run_input_key_t("driver", "optimize", "opt", KIND_BOOL), &
      run_input_key_t("analysis", "charges", "charges", KIND_STRING), &
      run_input_key_t("analysis", "fukui", "fukui", KIND_STRING), &
      run_input_key_t("hessian", "step", "step", KIND_FLOAT), &
      run_input_key_t("thermo", "temperature", "temperature", KIND_FLOAT), &
      run_input_key_t("thermo", "pressure", "pressure", KIND_FLOAT), &
      run_input_key_t("thermo", "symmetry_number", "symmetry-number", KIND_INT), &
      run_input_key_t("thermo", "masses", "masses", KIND_STRING), &
      run_input_key_t("opt", "coordinates", "opt-coords", KIND_STRING), &
      run_input_key_t("opt", "max_steps", "opt-maxiter", KIND_INT), &
      run_input_key_t("opt", "gradient_tol", "opt-gtol", KIND_FLOAT), &
      run_input_key_t("opt", "energy_tol", "opt-etol", KIND_FLOAT), &
      run_input_key_t("opt", "scf_guess", "opt-guess", KIND_STRING), &
      run_input_key_t("opt", "algorithm", "opt-algorithm", KIND_STRING), &
      run_input_key_t("opt", "target", "opt-target", KIND_STRING), &
      run_input_key_t("opt", "saddle_method", "opt-saddle-method", KIND_STRING), &
      run_input_key_t("opt", "connect", "opt-connect", KIND_BOOL), &
      run_input_key_t("opt", "initial_hessian", "opt-init-hessian", KIND_STRING), &
      run_input_key_t("opt", "hessian_update", "opt-hessian-update", KIND_STRING), &
      run_input_key_t("opt", "max_hessian_updates", "opt-max-updates", KIND_INT), &
      run_input_key_t("opt", "zero_modes", "opt-zero-modes", KIND_INT), &
      run_input_key_t("opt", "soft_modes", "opt-soft-modes", KIND_FLOAT), &
      run_input_key_t("opt", "max_step", "opt-max-step", KIND_FLOAT), &
      run_input_key_t("opt", "line_search", "opt-linesearch", KIND_STRING), &
      run_input_key_t("opt", "step_scale", "opt-scalestep", KIND_FLOAT), &
      run_input_key_t("opt", "lbfgs_memory", "opt-lbfgs-memory", KIND_INT), &
      run_input_key_t("opt", "print_level", "opt-print-level", KIND_INT), &
      run_input_key_t("opt", "frozen_atoms", "opt-frozen", KIND_INT_LIST), &
      run_input_key_t("opt", "constraints", "opt-constrain", KIND_STRING_LIST), &
      run_input_key_t("opt", "seed", "opt-seed", KIND_INT), &
      run_input_key_t("opt.neb", "product", "opt-neb", KIND_PATH), &
      run_input_key_t("opt.neb", "images", "opt-neb-images", KIND_INT), &
      run_input_key_t("opt.neb", "spring", "opt-neb-spring", KIND_FLOAT), &
      run_input_key_t("opt.neb", "climb", "opt-neb-climb", KIND_BOOL), &
      run_input_key_t("opt.neb", "endpoints", "opt-neb-endpoints", KIND_STRING), &
      run_input_key_t("opt.neb", "via", "opt-neb-via", KIND_PATH), &
      run_input_key_t("opt.global", "population", "opt-global-pop-size", KIND_INT), &
      run_input_key_t("opt.global", "initial_population", "opt-global-init-pop-size", KIND_INT), &
      run_input_key_t("opt.global", "max_cycles", "opt-global-maxcycle", KIND_INT), &
      run_input_key_t("opt.global", "keep_best", "opt-global-nsave", KIND_INT), &
      run_input_key_t("opt.global", "distribution", "opt-global-distribution", KIND_INT), &
      run_input_key_t("opt.global", "reset_after", "opt-global-reset", KIND_INT), &
      run_input_key_t("opt.global", "radius", "opt-global-radius", KIND_FLOAT), &
      run_input_key_t("opt.global", "contraction", "opt-global-contraction", KIND_FLOAT), &
      run_input_key_t("opt.global", "radius_tol", "opt-global-tolerance-r", KIND_FLOAT), &
      run_input_key_t("opt.global", "gradient_tol", "opt-global-tolerance-g", KIND_FLOAT), &
      run_input_key_t("opt.global", "mutation_rate", "opt-global-mutation-rate", KIND_FLOAT), &
      run_input_key_t("opt.global", "death_rate", "opt-global-death-rate", KIND_FLOAT), &
      run_input_key_t("opt.global", "scale_factor", "opt-global-scalefac", KIND_FLOAT), &
      run_input_key_t("output", "json", "json", KIND_PATH), &
      run_input_key_t("output", "save_scf", "save-scf", KIND_PATH), &
      run_input_key_t("output", "silent", "silent", KIND_BOOL), &
      run_input_key_t("output", "log_level", "logger", KIND_STRING)]

   type :: input_state_t
      !! What walking a file builds up: the flags so far, the geometry (the
      !! positional, kept apart so it can go first) and the first error.
      type(string_type), allocatable :: flags(:)
      integer :: nflags = 0
      character(len=:), allocatable :: geometry
      character(len=:), allocatable :: message
   end type input_state_t

   integer, parameter :: NSECTIONS = 15
   character(len=16), parameter :: SECTIONS(NSECTIONS) = [character(len=16) :: &
      "molecule", "basis", "scf", "scf.convergence", "dft", "solvation", "correlation", "driver", &
      "analysis", "hessian", "thermo", "opt", "opt.neb", "opt.global", "output"]

   character(len=*), parameter :: REAL_FORMAT = "es24.16e3"
      !! Every float goes through this, not the g0 default of pic, whose digits
      !! differ by compiler: es24.16 reads back to the same double.

contains

   pure function run_input_is_toml(arg) result(is_toml)
      !! Whether the first command-line argument names an input file.
      character(len=*), intent(in) :: arg
      logical :: is_toml
      integer :: n
      n = len_trim(arg)
      is_toml = .false.
      if (n > 5) is_toml = to_lower(arg(n - 4:n)) == ".toml"
   end function run_input_is_toml

   subroutine run_input_args(path, args, message)
      !! The command-line arguments `path` stands for: the geometry first
      !! (the positional), then one `--flag=value` per key. `message` is empty
      !! on success, else what is wrong and where in the file.
      character(len=*), intent(in) :: path
      type(string_type), allocatable, intent(out) :: args(:)
      character(len=:), allocatable, intent(out) :: message
      type(toml_table), allocatable :: table
      type(toml_error), allocatable :: error
      type(toml_context) :: context
      type(input_state_t) :: state
      character(len=:), allocatable :: dir
      integer :: slash

      message = ""
      allocate (state%flags(16))
      state%geometry = ""
      state%message = ""
      call toml_load(table, path, context=context, error=error)
      if (allocated(error)) then
         message = error%message
         return
      end if
      slash = index(path, "/", back=.true.)
      dir = path(1:slash)
      call walk(table, "", dir, context, state)
      if (len(state%message) > 0) then
         message = state%message
         return
      end if
      if (len(state%geometry) == 0) then
         message = path//": [molecule] geometry is required"
         return
      end if
      allocate (args(state%nflags + 1))
      args(1) = state%geometry
      args(2:) = state%flags(1:state%nflags)
   end subroutine run_input_args

   recursive subroutine walk(table, prefix, dir, context, state)
      !! Every key of `table` (section `prefix`), turned into flags; recurses
      !! into the sub-tables the format has.
      type(toml_table), intent(inout) :: table
      character(len=*), intent(in) :: prefix, dir
      type(toml_context), intent(in) :: context
      type(input_state_t), intent(inout) :: state
      type(toml_key), allocatable :: keys(:)
      type(toml_table), pointer :: child
      character(len=:), allocatable :: name
      integer :: i, row, stat

      call table%get_keys(keys)
      do i = 1, size(keys)
         if (len(prefix) == 0) then
            name = keys(i)%key
         else
            name = prefix//"."//keys(i)%key
         end if
         call get_value(table, keys(i)%key, child, requested=.false., stat=stat)
         if (associated(child)) then
            if (.not. any(SECTIONS == name)) then
               state%message = context%report("unknown section ["//name//"] (docs/INPUT.md lists them)", &
                                        keys(i)%origin)
               return
            end if
            call walk(child, name, dir, context, state)
            if (len(state%message) > 0) return
            cycle
         end if
         if (len(prefix) == 0) then
            state%message = context%report("'"//name//"' is outside any section; keys go under a [section]", &
                                     keys(i)%origin)
            return
         end if
         row = find_key(prefix, keys(i)%key)
         if (row == 0) then
            state%message = context%report("unknown key '"//keys(i)%key//"' in ["//prefix//"] (docs/INPUT.md)", &
                                     keys(i)%origin)
            return
         end if
         call convert(table, keys(i), RUN_INPUT_KEYS(row), dir, context, state)
         if (len(state%message) > 0) return
      end do
   end subroutine walk

   subroutine convert(table, key, spec, dir, context, state)
      !! The value of one key, checked against its kind and appended as its flag.
      type(toml_table), intent(inout) :: table
      type(toml_key), intent(in) :: key
      type(run_input_key_t), intent(in) :: spec
      character(len=*), intent(in) :: dir
      type(toml_context), intent(in) :: context
      type(input_state_t), intent(inout) :: state
      character(len=:), allocatable :: sval, flag, joined
      type(toml_array), pointer :: array
      integer(int64) :: ival
      real(dp) :: rval
      logical :: lval, exists
      integer :: stat, j

      flag = "--"//trim(spec%flag)
      select case (spec%kind)
      case (KIND_STRING, KIND_PATH, KIND_FILE_OR_NAME, KIND_METHOD)
         call get_value(table, key%key, sval, stat=stat)
         if (stat /= 0 .or. .not. allocated(sval)) then
            state%message = context%report("'"//key%key//"' must be a string", key%origin)
            return
         end if
         if (spec%kind == KIND_PATH) then
            sval = resolve(dir, sval)
         else if (spec%kind == KIND_FILE_OR_NAME) then
            inquire (file=resolve(dir, sval), exist=exists)
            if (exists) sval = resolve(dir, sval)
         else if (spec%kind == KIND_METHOD) then
            select case (to_lower(sval))
            case ("rimp2")
               call push(state, flag)
            case ("none")
            case default
               state%message = context%report("[correlation] method must be 'rimp2' or 'none', got '"//sval//"'", &
                                        key%origin)
            end select
            return
         end if
         if (len_trim(spec%flag) == 0) then
            state%geometry = sval
         else
            call push(state, flag//"="//sval)
         end if
      case (KIND_INT)
         call get_value(table, key%key, ival, stat=stat)
         if (stat /= 0) then
            state%message = context%report("'"//key%key//"' must be an integer", key%origin)
            return
         end if
         call push(state, flag//"="//to_string(ival))
      case (KIND_FLOAT)
         call get_value(table, key%key, rval, stat=stat)
         if (stat /= 0) then
            call get_value(table, key%key, ival, stat=stat)
            rval = real(ival, dp)
         end if
         if (stat /= 0) then
            state%message = context%report("'"//key%key//"' must be a number", key%origin)
            return
         end if
         call push(state, flag//"="//trim(adjustl(to_string(rval, REAL_FORMAT))))
      case (KIND_BOOL)
         call get_value(table, key%key, lval, stat=stat)
         if (stat /= 0) then
            state%message = context%report("'"//key%key//"' must be true or false", key%origin)
            return
         end if
         if (lval) call push(state, flag)
      case (KIND_BOOL_ON, KIND_BOOL_TRI)
         call get_value(table, key%key, lval, stat=stat)
         if (stat /= 0) then
            state%message = context%report("'"//key%key//"' must be true or false", key%origin)
            return
         end if
         if (lval) then
            call push(state, flag)
         else
            call push(state, "--no-"//trim(spec%flag))
         end if
      case (KIND_INT_LIST, KIND_STRING_LIST)
         call get_value(table, key%key, array, requested=.false., stat=stat)
         if (.not. associated(array)) then
            state%message = context%report("'"//key%key//"' must be an array", key%origin)
            return
         end if
         joined = ""
         do j = 1, len(array)
            if (spec%kind == KIND_INT_LIST) then
               call get_value(array, j, ival, stat=stat)
               if (stat /= 0) then
                  state%message = context%report("'"//key%key//"' must be an array of integers", key%origin)
                  return
               end if
               sval = to_string(ival)
            else
               call get_value(array, j, sval, stat=stat)
               if (stat /= 0 .or. .not. allocated(sval)) then
                  state%message = context%report("'"//key%key//"' must be an array of strings", key%origin)
                  return
               end if
            end if
            if (j > 1) joined = joined//merge(",", ";", spec%kind == KIND_INT_LIST)
            joined = joined//sval
         end do
         if (len(joined) > 0) call push(state, flag//"="//joined)
      case default
         state%message = "trc_run_input: key kind "//to_string(spec%kind)//" has no conversion"
      end select
   end subroutine convert

   subroutine run_input_echo(cli, json, root)
      !! The `"input"` object of `--json`: every key as the parsed command
      !! line has it (file, overrides and declared defaults together), in
      !! the sections of the file format, so it can be written back as a
      !! TOML file and rerun.
      type(cli_t), intent(in) :: cli
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(inout) :: root
      type(json_value), pointer :: input, section
      character(len=:), allocatable :: sval
      type(run_input_key_t) :: spec
      integer :: i
      real(dp) :: rval
      integer(int64) :: ival
      integer :: stat

      call json%create_object(input, "input")
      call json%add(root, input)
      do i = 1, NKEYS
         spec = RUN_INPUT_KEYS(i)
         if (spec%kind == KIND_BOOL) then
            call section_object(json, input, trim(spec%section), section)
            call json%add(section, trim(spec%key), cli%is_set(trim(spec%flag)))
            cycle
         end if
         if (spec%kind == KIND_BOOL_ON) then
            call section_object(json, input, trim(spec%section), section)
            call json%add(section, trim(spec%key), .not. cli%is_set("no-"//trim(spec%flag)))
            cycle
         end if
         if (spec%kind == KIND_BOOL_TRI) then
            ! Neither flag set is the third state, and echoing either value
            ! would turn it into a forced one on a rerun of the echo.
            if (.not. (cli%is_set(trim(spec%flag)) .or. cli%is_set("no-"//trim(spec%flag)))) cycle
            call section_object(json, input, trim(spec%section), section)
            call json%add(section, trim(spec%key), cli%is_set(trim(spec%flag)))
            cycle
         end if
         if (spec%kind == KIND_METHOD) then
            call section_object(json, input, trim(spec%section), section)
            call json%add(section, trim(spec%key), trim(merge("rimp2", "none ", cli%is_set(trim(spec%flag)))))
            cycle
         end if
         if (len_trim(spec%flag) == 0) then
            call cli%get("geometry", sval)
         else
            call cli%get(trim(spec%flag), sval)
         end if
         if (.not. allocated(sval)) cycle
         if (len(sval) == 0) cycle
         call section_object(json, input, trim(spec%section), section)
         select case (spec%kind)
         case (KIND_INT)
            read (sval, *, iostat=stat) ival
            if (stat == 0) then
               call json%add(section, trim(spec%key), int(ival))
            else
               call json%add(section, trim(spec%key), sval)
            end if
         case (KIND_FLOAT)
            read (sval, *, iostat=stat) rval
            if (stat == 0) then
               call json%add(section, trim(spec%key), rval)
            else
               call json%add(section, trim(spec%key), sval)
            end if
         case (KIND_INT_LIST)
            call add_int_list(json, section, trim(spec%key), sval)
         case (KIND_STRING_LIST)
            call add_string_list(json, section, trim(spec%key), sval)
         case default
            call json%add(section, trim(spec%key), sval)
         end select
      end do
   end subroutine run_input_echo

   recursive subroutine section_object(json, input, name, section)
      !! The object for section `name` (dotted for a sub-table) under
      !! `input`, made on first use.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(inout) :: input
      character(len=*), intent(in) :: name
      type(json_value), pointer, intent(out) :: section
      type(json_value), pointer :: parent
      logical :: found
      integer :: dot

      call json%get(input, name, section, found)
      if (found) return
      dot = index(name, ".", back=.true.)
      if (dot == 0) then
         parent => input
      else
         call section_object(json, input, name(1:dot - 1), parent)
      end if
      call json%create_object(section, name(dot + 1:))
      call json%add(parent, section)
   end subroutine section_object

   subroutine add_int_list(json, section, key, text)
      !! The `--opt-frozen` string "1-3,7" as [1, 2, 3, 7].
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(inout) :: section
      character(len=*), intent(in) :: key, text
      integer, allocatable :: atoms(:)
      integer :: start, comma, dash, lo, hi, stat1, stat2, k
      character(len=:), allocatable :: piece

      allocate (atoms(0))
      start = 1
      do while (start <= len(text))
         comma = index(text(start:), ",")
         if (comma == 0) then
            piece = trim(adjustl(text(start:)))
            start = len(text) + 1
         else
            piece = trim(adjustl(text(start:start + comma - 2)))
            start = start + comma
         end if
         dash = index(piece, "-")
         if (dash > 1) then
            read (piece(1:dash - 1), *, iostat=stat1) lo
            read (piece(dash + 1:), *, iostat=stat2) hi
         else
            read (piece, *, iostat=stat1) lo
            hi = lo
            stat2 = 0
         end if
         if (stat1 /= 0 .or. stat2 /= 0) then
            call json%add(section, key, text)
            return
         end if
         atoms = [atoms, (k, k=lo, hi)]
      end do
      call json%add(section, key, atoms)
   end subroutine add_int_list

   subroutine add_string_list(json, section, key, text)
      !! The `--opt-constrain` string "a;b" as ["a", "b"].
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(inout) :: section
      character(len=*), intent(in) :: key, text
      type(json_value), pointer :: array
      integer :: start, semi

      call json%create_array(array, key)
      start = 1
      do while (start <= len(text))
         semi = index(text(start:), ";")
         if (semi == 0) then
            call json%add(array, "", trim(adjustl(text(start:))))
            exit
         end if
         call json%add(array, "", trim(adjustl(text(start:start + semi - 2))))
         start = start + semi
      end do
      call json%add(section, array)
   end subroutine add_string_list

   pure function find_key(section, key) result(row)
      !! Row of RUN_INPUT_KEYS for `key` in `section`, 0 if there is none.
      character(len=*), intent(in) :: section, key
      integer :: row
      integer :: i
      row = 0
      do i = 1, NKEYS
         if (trim(RUN_INPUT_KEYS(i)%section) == section .and. trim(RUN_INPUT_KEYS(i)%key) == key) then
            row = i
            return
         end if
      end do
   end function find_key

   pure function resolve(dir, path) result(resolved)
      !! `path` relative to the directory of the TOML file, unless absolute.
      character(len=*), intent(in) :: dir, path
      character(len=:), allocatable :: resolved
      if (len(path) > 0) then
         if (path(1:1) == "/") then
            resolved = path
            return
         end if
      end if
      resolved = dir//path
   end function resolve

   subroutine push(state, flag)
      !! Append one argument, growing the list by doubling.
      type(input_state_t), intent(inout) :: state
      character(len=*), intent(in) :: flag
      type(string_type), allocatable :: grown(:)
      if (state%nflags == size(state%flags)) then
         allocate (grown(2*size(state%flags)))
         grown(1:state%nflags) = state%flags(1:state%nflags)
         call move_alloc(grown, state%flags)
      end if
      state%nflags = state%nflags + 1
      state%flags(state%nflags) = flag
   end subroutine push

end module trc_run_input
