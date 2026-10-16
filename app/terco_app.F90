module terco_app
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
   implicit none
   private

   public :: read_xyz, resolve_basis, element_symbol, c_string, f_string, set_env, lower

   real(dp), parameter, public :: BOHR_PER_ANGSTROM = 1.0_dp/0.529177210903_dp
      !! CODATA 2018, the constant every reference terco is checked against uses.

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

contains

   pure function lower(s) result(r)
      !! ASCII lower case.
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, c
      r = s
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) r(i:i) = achar(c + 32)
      end do
   end function lower

   pure function element_symbol(z) result(sym)
      !! The symbol for nuclear charge `z`, or "?" outside 1..118.
      integer, intent(in) :: z
      character(len=2) :: sym
      sym = "? "
      if (z >= 1 .and. z <= NELEM) sym = SYMBOLS(z)
   end function element_symbol

   pure integer function atomic_number(tok) result(z)
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
      s = lower(adjustl(tok))
      do i = 1, NELEM
         if (lower(SYMBOLS(i)) == s) then
            z = i
            return
         end if
      end do
   end function atomic_number

   subroutine read_xyz(path, natm, z, xyz, title, err)
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
                         itoa(natm)//" atoms, the file ends after "//itoa(i - 1))
            close (u)
            return
         end if
         read (line, *, iostat=ios) tok, r
         if (ios /= 0) then
            call err%set(int(ERROR_PARSE, default_int), trim(path)//": atom "//itoa(i)// &
                         " is not 'symbol x y z': '"//trim(line)//"'")
            close (u)
            return
         end if
         z(i) = atomic_number(tok)
         if (z(i) == 0) then
            call err%set(int(ERROR_PARSE, default_int), trim(path)//": atom "//itoa(i)// &
                         ": unknown element '"//trim(tok)//"'")
            close (u)
            return
         end if
         xyz(:, i) = r*BOHR_PER_ANGSTROM
      end do
      close (u)
   end subroutine read_xyz

   subroutine resolve_basis(name, path, err)
      !! A basis set given on the command line, as a file path.
      !!
      !! Anything that names an existing file is used as is. Otherwise it is a
      !! name: lower-cased, with ".json" added unless it is already there, and
      !! looked up in each directory of $TERCO_BASIS_PATH (colon-separated) in
      !! order, then in the installed basis directory. The files are Basis Set
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
      file = lower(trim(name))
      if (len(file) < 5) then
         file = file//".json"
      else if (file(len(file) - 4:) /= ".json") then
         file = file//".json"
      end if

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
   end subroutine resolve_basis

   logical function exists(path)
      character(len=*), intent(in) :: path
      inquire (file=trim(path), exist=exists)
   end function exists

   pure function c_string(s) result(cs)
      !! `s` NUL-terminated, for a `character(kind=c_char) :: x(*)` dummy.
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function c_string

   pure function f_string(buf) result(s)
      !! The text of a NUL-terminated C buffer.
      character(kind=c_char), intent(in) :: buf(:)
      character(len=:), allocatable :: s
      integer :: i, n
      n = size(buf)
      do i = 1, size(buf)
         if (buf(i) == c_null_char) then
            n = i - 1
            exit
         end if
      end do
      allocate (character(len=n) :: s)
      do i = 1, n
         s(i:i) = buf(i)
      end do
   end function f_string

   subroutine set_env(name, value)
      !! POSIX setenv. The library reads a few settings (the SAP fit file
      !! among them) only from the environment.
      character(len=*), intent(in) :: name, value
      integer(c_int) :: rc
      rc = c_setenv(c_string(name), c_string(value), 1_c_int)
   end subroutine set_env

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: b
      write (b, "(i0)") i
      s = trim(b)
   end function itoa

end module terco_app
