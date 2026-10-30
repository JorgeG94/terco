module trc_app
   !! What the `terco` command needs that the library does not: reading a
   !! geometry file, finding a basis set by name, and talking to the C API's
   !! NUL-terminated strings.
   !!
   !! Support for the `terco` command line program (app/main.F90).
   !!
   !! The program drives the library through its context API (trc_capi), the
   !! same surface a C or Python host uses, so every guess and option the
   !! context offers is reachable from the command line without a second
   !! implementation of any of it. What lives here is only what a context
   !! does not do: turn an xyz file into nuclear charges and Bohr
   !! coordinates, and turn `--basis def2-svp` into a file path.
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp, default_int
   use pic_error, only: error_t, ERROR_IO, ERROR_PARSE, ERROR_VALIDATION
   use pic_ascii, only: to_lower
   use pic_logger, only: global_logger, debug_level, verbose_level, large_info_level, info_level, &
                         performance_level, warning_level, error_level, knowledge_level
   use pic_strings, only: to_string, ends_with
   use trc_constants, only: ANGSTROM_TO_BOHR
   use json_module, only: json_core, json_value
   implicit none
   private

   public :: app_read_xyz, app_resolve_basis, app_element_symbol, app_atomic_number, app_c_string, app_f_string, app_set_env
   public :: BOHR_PER_ANGSTROM
   public :: app_help_for_topic, app_set_log_level
   public :: app_write_opt_xyz_frame, app_write_traj_xyz_frames, app_write_opt_failure_json

   real(dp), parameter :: BOHR_PER_ANGSTROM = ANGSTROM_TO_BOHR
      !! CODATA 2018 (trc_constants' ANGSTROM_TO_BOHR), the constant every
      !! reference terco is checked against uses. Kept under this name
      !! rather than renamed at the use site.

   integer, parameter :: NELEM = 118
   character(len=2), parameter :: SYMBOLS(NELEM) = [character(len=2) :: &
      "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
      "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", &
      "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
      "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
      "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
      "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
      "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
      "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", &
      "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", &
      "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

   interface
      function c_setenv(name, value, overwrite) result(rc) bind(c, name="setenv")
         import :: c_char, c_int
         implicit none
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: rc
      end function c_setenv
   end interface

#ifndef TERCO_BASIS_DIR
#define TERCO_BASIS_DIR ""
#endif
   character(len=*), parameter :: INSTALLED_BASIS_DIR = TERCO_BASIS_DIR
      !! Where `cmake --install` puts the basis sets; empty in an fpm build.
#ifndef TERCO_SOURCE_BASIS_DIR
#define TERCO_SOURCE_BASIS_DIR ""
#endif
   character(len=*), parameter :: SOURCE_BASIS_DIR = TERCO_SOURCE_BASIS_DIR
      !! The source tree's test/basis_sets, searched last: a binary run from
      !! its build tree, never installed, still finds the data files terco
      !! ships there (minao.json, hirshfeld_proatoms.json, ...) when
      !! $TERCO_BASIS_PATH names only someone else's basis directory.

contains

   pure function app_element_symbol(z) result(sym)
      !! The symbol for nuclear charge `z`, or "?" outside 1..118.
      integer, intent(in) :: z
      character(len=2) :: sym
      sym = "? "
      if (z >= 1 .and. z <= NELEM) sym = SYMBOLS(z)
   end function app_element_symbol

   pure integer function app_atomic_number(tok) result(z)
      !! Nuclear charge from an element symbol in any case, or from an
      !! integer (some xyz writers use Z). Zero if it is neither.
      character(len=*), intent(in) :: tok
      integer :: i, ios
      character(len=2) :: s
      z = 0
      read (tok, *, iostat=ios) i
      if (ios == 0) then
         if (i >= 1 .and. i <= NELEM) z = i
         return
      end if
      if (len_trim(tok) > 2) return
      s = to_lower(adjustl(tok))
      do i = 1, NELEM
         if (to_lower(SYMBOLS(i)) == s) then
            z = i
            return
         end if
      end do
   end function app_atomic_number

   subroutine app_read_xyz(path, natm, z, xyz, title, err)
      !! A standard xyz file: atom count, a title line, then one
      !! `symbol x y z` line per atom in Angstrom. Coordinates come back in
      !! Bohr, atom-major (3, natm), which is what trc_set_molecule takes.
      character(len=*), intent(in) :: path
      integer, intent(out) :: natm
      integer, allocatable, intent(out) :: z(:)
      real(dp), allocatable, intent(out) :: xyz(:, :)
      character(len=:), allocatable, intent(out) :: title
      type(error_t), intent(inout) :: err
      integer :: u, ios, i
      character(len=1024) :: line
      character(len=16) :: tok
      real(dp) :: r(3)

      natm = 0
      title = ""
      open (newunit=u, file=path, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         call err%set(int(ERROR_IO, default_int), "cannot open geometry file '"//trim(path)//"'")
         return
      end if
      read (u, *, iostat=ios) natm
      if (ios /= 0 .or. natm < 1) then
         call err%set(int(ERROR_PARSE, default_int), trim(path)// &
                      ": the first line must be the number of atoms (an xyz file)")
         close (u)
         return
      end if
      read (u, "(a)", iostat=ios) line
      title = trim(line)
      allocate (z(natm), xyz(3, natm))
      do i = 1, natm
         read (u, "(a)", iostat=ios) line
         if (ios /= 0) then
            call err%set(int(ERROR_PARSE, default_int), trim(path)//": expected "// &
                         to_string(natm)//" atoms, the file ends after "//to_string(i - 1))
            close (u)
            return
         end if
         read (line, *, iostat=ios) tok, r
         if (ios /= 0) then
            call err%set(int(ERROR_PARSE, default_int), trim(path)//": atom "//to_string(i)// &
                         " is not 'symbol x y z': '"//trim(line)//"'")
            close (u)
            return
         end if
         z(i) = app_atomic_number(tok)
         if (z(i) == 0) then
            call err%set(int(ERROR_PARSE, default_int), trim(path)//": atom "//to_string(i)// &
                         ": unknown element '"//trim(tok)//"'")
            close (u)
            return
         end if
         xyz(:, i) = r*BOHR_PER_ANGSTROM
      end do
      close (u)
   end subroutine app_read_xyz

   subroutine app_resolve_basis(name, path, err)
      !! A basis set given on the command line, as a file path.
      !!
      !! Anything that names an existing file is used as is. Otherwise it is a
      !! name: lower-cased, with ".json" added unless it is already there, and
      !! looked up in each directory of $TERCO_BASIS_PATH (colon-separated) in
      !! order, then in the installed basis directory, then in the source
      !! tree's test/basis_sets. The files are Basis Set
      !! Exchange JSON, as trc_basis_from_json reads them.
      character(len=*), intent(in) :: name
      character(len=:), allocatable, intent(out) :: path
      type(error_t), intent(inout) :: err
      character(len=:), allocatable :: file, dirs, tried
      character(len=4096) :: env
      integer :: n, i0, i1, stat

      path = ""
      if (exists(name)) then
         path = trim(name)
         return
      end if
      file = to_lower(trim(name))
      if (.not. ends_with(file, ".json")) file = file//".json"

      call get_environment_variable("TERCO_BASIS_PATH", env, n, stat)
      dirs = ""
      if (stat == 0 .and. n > 0) dirs = trim(env)
      if (len_trim(INSTALLED_BASIS_DIR) > 0) then
         if (len(dirs) > 0) then
            dirs = dirs//":"//INSTALLED_BASIS_DIR
         else
            dirs = INSTALLED_BASIS_DIR
         end if
      end if
      if (len_trim(SOURCE_BASIS_DIR) > 0) then
         if (len(dirs) > 0) then
            dirs = dirs//":"//SOURCE_BASIS_DIR
         else
            dirs = SOURCE_BASIS_DIR
         end if
      end if

      tried = ""
      i0 = 1
      do while (i0 <= len(dirs))
         i1 = index(dirs(i0:), ":")
         if (i1 == 0) then
            i1 = len(dirs) + 1
         else
            i1 = i0 + i1 - 1
         end if
         if (i1 > i0) then
            if (exists(dirs(i0:i1 - 1)//"/"//file)) then
               path = dirs(i0:i1 - 1)//"/"//file
               return
            end if
            tried = tried//new_line("a")//"    "//dirs(i0:i1 - 1)//"/"//file
         end if
         i0 = i1 + 1
      end do
      if (len(tried) == 0) then
         call err%set(int(ERROR_VALIDATION, default_int), "basis '"//trim(name)// &
                      "' is not a file, and there is nowhere to look it up: set TERCO_BASIS_PATH "// &
                      "to a directory of Basis Set Exchange JSON files, or pass a path")
      else
         call err%set(int(ERROR_VALIDATION, default_int), "basis '"//trim(name)// &
                      "' is not a file and was not found as:"//tried)
      end if
   end subroutine app_resolve_basis

   logical function exists(path) result(res)
      character(len=*), intent(in) :: path
      inquire (file=trim(path), exist=res)
   end function exists

   pure function app_c_string(s) result(cs)
      !! `s` NUL-terminated, for a `character(kind=c_char) :: x(*)` dummy.
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function app_c_string

   pure function app_f_string(buffer) result(s)
      !! The text of a NUL-terminated C buffer.
      character(kind=c_char), intent(in) :: buffer(:)
      character(len=:), allocatable :: s
      integer :: i, n
      n = size(buffer)
      do i = 1, size(buffer)
         if (buffer(i) == c_null_char) then
            n = i - 1
            exit
         end if
      end do
      allocate (character(len=n) :: s)
      do i = 1, n
         s(i:i) = buffer(i)
      end do
   end function app_f_string

   subroutine app_set_env(name, value)
      !! POSIX setenv. The library reads a few settings (the SAP fit file
      !! among them) only from the environment.
      character(len=*), intent(in) :: name, value
      integer(c_int) :: rc
      rc = c_setenv(app_c_string(name), app_c_string(value), 1_c_int)
   end subroutine app_set_env

   subroutine app_write_opt_xyz_frame(path, natm, z, energy, xyz_bohr, label)
      !! `<geometry base>_opt.xyz`: one geometry in Angstrom, atom order
      !! unchanged, comment line carrying `label` and the energy. Shared by
      !! an ordinary --opt finish (app/main.F90's own write_opt_xyz) and
      !! trc_dlfind_bridge's unrecoverable-failure path, which has no other
      !! chance to write this file: DL-FIND's own api.f90 calls `error stop`
      !! right after a gradient evaluation reports status=1, so nothing runs
      !! after api_dl_find returns on that path.
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: xyz_bohr(3, natm)
      character(len=*), intent(in) :: label
      integer :: u, ia
      open (newunit=u, file=path, status="replace", action="write")
      write (u, "(i0)") natm
      write (u, "(a,a,f0.12)") trim(label), "   E = ", energy
      do ia = 1, natm
         write (u, "(a2,3f20.12)") app_element_symbol(z(ia)), xyz_bohr(:, ia)/BOHR_PER_ANGSTROM
      end do
      close (u)
   end subroutine app_write_opt_xyz_frame

   subroutine app_write_traj_xyz_frames(path, natm, z, n_traj, traj_xyz_bohr, traj_energy)
      !! `<geometry base>_traj.xyz`: one frame per geometry handed in, in
      !! Angstrom, comment line carrying the step index and energy. The
      !! caller decides what belongs in `traj_xyz_bohr` -- app/main.F90
      !! passes the ACCEPTED trajectory a finished run collected;
      !! trc_dlfind_bridge passes the same trajectory so far when DL-FIND
      !! is about to `error stop` on a failed evaluation.
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm, n_traj
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: traj_xyz_bohr(3, natm, n_traj)
      real(dp), intent(in) :: traj_energy(n_traj)
      integer :: u, ia, k
      open (newunit=u, file=path, status="replace", action="write")
      do k = 1, n_traj
         write (u, "(i0)") natm
         write (u, "(a,i0,a,f0.12)") "step ", k, "   E = ", traj_energy(k)
         do ia = 1, natm
            write (u, "(a2,3f20.12)") app_element_symbol(z(ia)), traj_xyz_bohr(:, ia, k)/BOHR_PER_ANGSTROM
         end do
      end do
      close (u)
   end subroutine app_write_traj_xyz_frames

   subroutine app_write_opt_failure_json(path, geometry, reason, steps, energy_evaluations, &
                                     final_energy, final_gradient_max, final_gradient_rms, &
                                     coords_name, guess_name, traj_energy)
      !! A minimal "optimization" JSON block for the one outcome
      !! app/main.F90's own add_optimization never gets a chance to write:
      !! DL-FIND's `error stop` right after a failed evaluation
      !! (trc_dlfind_bridge%cb_get_gradient never returns to anything that
      !! could call add_optimization). Smaller than add_optimization's own
      !! block on purpose -- no gradient_tolerance, no final geometry array
      !! -- but carries the two fields a caller checking a failed run needs:
      !! "failed" and "failure_reason".
      character(len=*), intent(in) :: path, geometry, reason, coords_name, guess_name
      integer, intent(in) :: steps, energy_evaluations
      real(dp), intent(in) :: final_energy, final_gradient_max, final_gradient_rms
      real(dp), intent(in) :: traj_energy(:)
      type(json_core) :: json
      type(json_value), pointer :: root, oobj
      call json%initialize()
      call json%create_object(root, "")
      call json%add(root, "program", "terco")
      call json%add(root, "geometry", geometry)
      call json%create_object(oobj, "optimization")
      call json%add(oobj, "converged", .false.)
      call json%add(oobj, "failed", .true.)
      call json%add(oobj, "failure_reason", reason)
      call json%add(oobj, "steps", steps)
      call json%add(oobj, "energy_evaluations", energy_evaluations)
      call json%add(oobj, "final_energy", final_energy)
      call json%add(oobj, "final_gradient_max", final_gradient_max)
      call json%add(oobj, "final_gradient_rms", final_gradient_rms)
      call json%add(oobj, "coordinates", coords_name)
      call json%add(oobj, "guess", guess_name)
      if (size(traj_energy) > 0) call json%add(oobj, "trajectory_energies", traj_energy)
      call json%add(root, oobj)
      call json%print(root, path)
      call json%destroy(root)
   end subroutine app_write_opt_failure_json

   function app_help_for_topic(text, topic) result(out)
      !! `pic_cli`'s help text cut down to one topic. Without a topic the
      !! `--opt-*` tuning options (about fifty lines) are dropped, except
      !! `--opt-freq`, which is a headline switch; with `--help opt` only the
      !! options named `opt` or `opt-*` are kept. `pic_cli` has no notion of
      !! option groups, so this filters its rendered lines on the long name.
      character(len=*), intent(in) :: text, topic
      character(len=:), allocatable :: out
      character(len=:), allocatable :: line, name, t
      integer :: pos, nl, eq
      logical :: is_option, in_group

      t = to_lower_alloc_app(trim(adjustl(topic)))
      out = ""
      pos = 1
      do while (pos <= len(text))
         nl = index(text(pos:), new_line("a"))
         if (nl == 0) then
            line = text(pos:)
            pos = len(text) + 1
         else
            line = text(pos:pos + nl - 2)
            pos = pos + nl
         end if
         eq = index(line, "--")
         is_option = eq > 0 .and. eq <= 8 .and. index(line(:eq - 1), "[") == 0 .and. line(1:2) == "  "
         if (is_option) then
            name = line(eq + 2:)
            name = name(:index(name//" ", " ") - 1)
            if (len(t) > 0) then
               in_group = name == t .or. (len(name) > len(t) .and. name(:len(t) + 1) == t//"-")
               if (.not. in_group) cycle
            else
               if (len(name) > 4 .and. name(:4) == "opt-" .and. name /= "opt-freq") cycle
            end if
         else if (len(t) > 0 .and. (index(line, "Usage:") == 1 .or. index(line, "Options:") == 1 &
                                    .or. index(line, "Arguments:") == 1 .or. len_trim(line) == 0)) then
            if (index(line, "Usage:") /= 1 .and. index(line, "Options:") /= 1) cycle
         end if
         out = out//line//new_line("a")
      end do
      if (len(t) == 0) then
         out = out//new_line("a")//"Run 'terco --help opt' for the geometry-optimisation (--opt-*) options."// &
               new_line("a")
      end if
   end function app_help_for_topic

   subroutine app_set_log_level(name, ok)
      !! Apply `--logger <name>` to the global logger. Names are pic_logger's
      !! own levels, quietest first: knowledge, error, warning, performance,
      !! info, large_info, verbose, debug.
      character(len=*), intent(in) :: name
      logical, intent(out) :: ok

      ok = .true.
      select case (to_lower_alloc_app(trim(adjustl(name))))
      case ("knowledge")
      call global_logger%configure(knowledge_level)
      case ("error")
      call global_logger%configure(error_level)
      case ("warning")
      call global_logger%configure(warning_level)
      case ("performance")
      call global_logger%configure(performance_level)
      case ("info")
      call global_logger%configure(info_level)
      case ("large_info")
      call global_logger%configure(large_info_level)
      case ("verbose")
      call global_logger%configure(verbose_level)
      case ("debug")
      call global_logger%configure(debug_level)
      case default
      ok = .false.
      end select
   end subroutine app_set_log_level

   function to_lower_alloc_app(s) result(lower_s)
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: lower_s
      integer :: i
      lower_s = s
      do i = 1, len(s)
         lower_s(i:i) = to_lower(s(i:i))
      end do
   end function to_lower_alloc_app

end module trc_app
