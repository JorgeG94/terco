program terco_main
   !! The `terco` command: an SCF on an xyz file, from the shell.
   !!
   !! `terco molecule.xyz [options]`, in the spirit of xtb's command line.
   !!
   !! Everything below the argument handling is the library's context API
   !! (trc_capi): the guesses, the SCF, RI-MP2 and the energies are exactly
   !! what a C or Python host gets. `terco --help` lists the options.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_ptr
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, int64
   use pic_types, only: dp, default_int
   use pic_error, only: error_t
   use pic_cli, only: cli_t
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_aux_json, trc_set_method, trc_set_convergence, trc_set_screening, &
                       trc_set_guess, trc_set_verbose, trc_set_rimp2, trc_set_force_cartesian, &
                       trc_run_scf, trc_run_rimp2, trc_nao, trc_is_spherical, trc_nspin, trc_energy, &
                       trc_energy_parts, trc_converged, trc_iterations, trc_rimp2_energy, trc_message, &
                       TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_CORE, TRC_GUESS_GWH, TRC_GUESS_SAD, &
                       TRC_GUESS_SADQ, TRC_GUESS_SAP
   use trc_api, only: trc_bind_device
   use terco_app, only: read_xyz, resolve_basis, element_symbol, c_string, f_string, set_env, lower
   implicit none

#ifndef TERCO_VERSION
#define TERCO_VERSION "unknown"
#endif
   character(len=*), parameter :: VERSION = TERCO_VERSION

   type(cli_t) :: cli
   type(error_t) :: err
   type(c_ptr) :: h
   character(len=:), allocatable :: geom, basis_arg, basis_path, aux_path, func, guess, title, json
   character(len=:), allocatable :: sval
   integer, allocatable :: z(:)
   real(dp), allocatable :: xyz(:, :)
   integer :: natm, chrg, uhf, nelec, mult, grid, maxiter, nfrozen, guess_kind, ia, dev
   integer(c_int) :: st, nao, nspin, conv, niter, sph
   real(dp) :: econv, gconv, thresh
   real(c_double) :: etot, enuc, eone, etwo, exc, eos, ess
   logical :: rimp2, cart
   integer(int64) :: t0, t1, rate

   call system_clock(t0, rate)

   ! ---- the interface ---------------------------------------------------
   call cli%set_program("terco", "Hartree-Fock and Kohn-Sham DFT, in standard Fortran")
   call cli%add_positional("geometry", "molecule in xyz format, Angstrom", required=.false.)
   call cli%add_option("basis", "basis set: a Basis Set Exchange JSON file, or a name looked up "// &
                       "in $TERCO_BASIS_PATH and then the installed basis sets (def2-svp)", short="b", default="def2-svp")
   call cli%add_option("chrg", "molecular charge (default: a .CHRG file here, else 0)", short="c")
   call cli%add_option("uhf", "number of unpaired electrons (default: a .UHF file here, else 0 "// &
                       "or 1 by electron count)", short="u")
   call cli%add_option("func", "hf (the default), or a functional: lda (svwn), pbe, blyp, pbe0, b3lyp, b3lyp5", &
                       short="f", default="hf")
   call cli%add_option("guess", "sad, sadq, sap, gwh or core (default: sad; sadq for an ion)", short="g")
   call cli%add_flag("cart", "Cartesian AOs even when the basis is defined spherical")
   call cli%add_option("grid", "DFT integration grid level (3)", default="3")
   call cli%add_option("iterations", "maximum SCF iterations (100)", default="100")
   call cli%add_option("econv", "energy convergence, Eh (1e-10)", default="1e-10")
   call cli%add_option("gconv", "convergence on the DIIS error norm |FDS - SDF| (1e-5)", default="1e-5")
   call cli%add_option("thresh", "integral screening threshold (1e-12)", default="1e-12")
   call cli%add_flag("rimp2", "RI-MP2 correlation energy after a closed-shell SCF")
   call cli%add_option("aux", "auxiliary basis for --rimp2 (cc-pvdz-rifit)", default="cc-pvdz-rifit")
   call cli%add_option("frozen", "frozen core orbitals for --rimp2; -1, the default, counts them by element", default="-1")
   call cli%add_option("json", "also write a JSON summary to this file")
   call cli%add_flag("silent", "only the summary: no SCF iterations, timings or guess details", short="s")
   call cli%add_flag("version", "print the version and stop")

   call cli%parse(err)
   if (err%has_error()) call die(err%get_message()//new_line("a")//"try 'terco --help'")
   if (cli%help_requested()) then
      write (output_unit, "(a)") cli%help_text()
      stop
   end if
   if (cli%is_set("version")) then
      write (output_unit, "(a)") "terco version "//VERSION
      stop
   end if
   call cli%get("geometry", geom, err)
   if (err%has_error() .or. len_trim(geom) == 0) &
      call die("no geometry given: terco molecule.xyz [options]"//new_line("a")//"try 'terco --help'")

   ! ---- inputs ----------------------------------------------------------
   call read_xyz(geom, natm, z, xyz, title, err)
   if (err%has_error()) call die(err%get_message())

   call cli%get("basis", basis_arg, err)
   call resolve_basis(basis_arg, basis_path, err)
   if (err%has_error()) call die(err%get_message())

   chrg = int_option("chrg", ".CHRG", 0)
   nelec = sum(z) - chrg
   if (nelec < 1) call die("charge "//itoa(chrg)//" leaves "//itoa(nelec)//" electrons")
   uhf = int_option("uhf", ".UHF", mod(nelec, 2))
   if (uhf < 0 .or. uhf > nelec .or. mod(nelec - uhf, 2) /= 0) &
      call die("--uhf "//itoa(uhf)//" is impossible with "//itoa(nelec)//" electrons: the unpaired "// &
               "count must be between 0 and the electron count, with the same parity")
   mult = uhf + 1

   call cli%get("func", sval, err)
   func = lower(trim(sval))
   if (func == "hf" .or. func == "rhf" .or. func == "uhf") func = ""

   guess = "sad"
   if (chrg /= 0) guess = "sadq"
   if (cli%is_set("guess")) then
      call cli%get("guess", sval, err)
      guess = lower(trim(sval))
   end if
   select case (guess)
   case ("sad"); guess_kind = TRC_GUESS_SAD
   case ("sadq"); guess_kind = TRC_GUESS_SADQ
   case ("sap"); guess_kind = TRC_GUESS_SAP
   case ("gwh"); guess_kind = TRC_GUESS_GWH
   case ("core", "hcore"); guess_kind = TRC_GUESS_CORE
   case default
      call die("unknown guess '"//guess//"': use sad, sadq, sap, gwh or core")
   end select

   call cli%get("grid", grid, err); call check("--grid")
   call cli%get("iterations", maxiter, err); call check("--iterations")
   call cli%get("econv", econv, err); call check("--econv")
   call cli%get("gconv", gconv, err); call check("--gconv")
   call cli%get("thresh", thresh, err); call check("--thresh")
   call cli%get("frozen", nfrozen, err); call check("--frozen")
   rimp2 = cli%is_set("rimp2")
   cart = cli%is_set("cart")
   if (rimp2) then
      ! RI-MP2 is closed-shell and, today, Cartesian only: say so here rather
      ! than after an SCF it cannot use.
      if (mult /= 1) call die("--rimp2 needs a closed-shell reference, and this molecule has "// &
                              itoa(uhf)//" unpaired electrons")
      cart = .true.
      call cli%get("aux", sval, err)
      call resolve_basis(sval, aux_path, err)
      if (err%has_error()) call die("--aux: "//err%get_message())
   end if
   json = ""
   if (cli%is_set("json")) call cli%get("json", json, err)

   ! SAP's fitted atomic potentials ship as a basis-format file; find it the
   ! way a basis is found, unless the caller has already said where it is.
   if (guess_kind == TRC_GUESS_SAP) then
      block
         character(len=1024) :: e
         integer :: n, s
         call get_environment_variable("TRC_SAP_BASIS", e, n, s)
         if (s /= 0 .or. n == 0) then
            call resolve_basis("sap_helfem_large", sval, err)
            if (err%has_error()) call die("--guess sap needs its fitted potentials: "//err%get_message())
            call set_env("TRC_SAP_BASIS", sval)
         end if
      end block
   end if

   call banner()

   ! ---- the calculation -------------------------------------------------
   ! One process, one device: the first. Before anything reaches the device,
   ! as trc_bind_device asks; a host build has nothing to bind.
   dev = trc_bind_device(0)
   st = trc_create(h)
   call api(trc_set_verbose(h, merge(0_c_int, 1_c_int, cli%is_set("silent"))))
   call api(trc_set_molecule(h, int(natm, c_int), real(z, c_double), xyz, int(chrg, c_int), int(mult, c_int)))
   call api(trc_set_force_cartesian(h, merge(1_c_int, 0_c_int, cart)))
   call api(trc_set_basis_json(h, c_string(basis_path)))
   if (rimp2) call api(trc_set_aux_json(h, c_string(aux_path)))
   call api(trc_set_method(h, c_string(func), int(grid, c_int)))
   call api(trc_set_convergence(h, real(econv, c_double), real(gconv, c_double), int(maxiter, c_int), 10_c_int))
   call api(trc_set_screening(h, real(thresh, c_double)))
   call api(trc_set_guess(h, int(guess_kind, c_int), c_null_ptr, 1_c_int))
   if (rimp2) call api(trc_set_rimp2(h, int(nfrozen, c_int), 0_c_int))

   st = trc_run_scf(h)
   if (st /= TRC_OK .and. st /= TRC_ERR_NOCONV) call api(st)
   call api(trc_nao(h, nao))
   call api(trc_is_spherical(h, sph))
   call api(trc_nspin(h, nspin))
   call api(trc_energy(h, etot))
   call api(trc_energy_parts(h, enuc, eone, etwo, exc))
   call api(trc_converged(h, conv))
   call api(trc_iterations(h, niter))
   eos = 0; ess = 0
   if (rimp2 .and. conv /= 0) then
      call api(trc_run_rimp2(h))
      call api(trc_rimp2_energy(h, eos, ess))
   end if

   call system_clock(t1)
   call summary()
   if (len(json) > 0) call write_json(json)
   st = trc_destroy(h)

   if (conv == 0) then
      write (error_unit, "(a)") "terco: the SCF did not converge in "//itoa(int(niter))//" iterations"
      stop 2, quiet=.true.
   end if
   write (output_unit, "(/,a)") " normal termination of terco"

contains

   subroutine die(msg)
      character(len=*), intent(in) :: msg
      write (error_unit, "(a)") "terco: "//msg
      stop 1, quiet=.true.
   end subroutine die

   subroutine check(what)
      character(len=*), intent(in) :: what
      if (err%has_error()) call die(what//": "//err%get_message())
   end subroutine check

   subroutine api(status)
      !! Stop on a failed library call, with the library's own explanation.
      integer(c_int), intent(in) :: status
      character(kind=c_char) :: buf(1024)
      integer(c_int) :: s
      if (status == TRC_OK) return
      s = trc_message(h, buf, int(size(buf), c_int))
      call die(f_string(buf))
   end subroutine api

   integer function int_option(name, file, default) result(v)
      !! An integer option, else the first integer in `file` in the working
      !! directory (xtb's .CHRG/.UHF convention), else `default`.
      character(len=*), intent(in) :: name, file
      integer, intent(in) :: default
      integer :: u, ios
      logical :: there
      v = default
      if (cli%is_set(name)) then
         call cli%get(name, v, err)
         call check("--"//name)
         return
      end if
      inquire (file=file, exist=there)
      if (.not. there) return
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      read (u, *, iostat=ios) v
      close (u)
      if (ios /= 0) call die(file//" does not start with an integer")
   end function int_option

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: b
      write (b, "(i0)") i
      s = trim(b)
   end function itoa

   function method() result(m)
      character(len=:), allocatable :: m
      if (len(func) == 0) then
         m = merge("RHF", "UHF", mult == 1)
      else
         m = merge("RKS", "UKS", mult == 1)//"/"//func
      end if
   end function method

   subroutine banner()
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") "      -----------------------------------------------------------"
      write (output_unit, "(a)") "     |                  t e r c o   "//VERSION
      write (output_unit, "(a)") "     |   Hartree-Fock and DFT through do concurrent, in Fortran"
      write (output_unit, "(a)") "      -----------------------------------------------------------"
      write (output_unit, "(a)") ""
      write (output_unit, "(a,a)") "   geometry         : ", geom
      if (len_trim(title) > 0) write (output_unit, "(a,a)") "   title            : ", trim(title)
      write (output_unit, "(a,i0,a)", advance="no") "   atoms            : ", natm, "  ("
      call formula()
      write (output_unit, "(a)") ")"
      write (output_unit, "(a,i0,a,i0,a,i0)") "   charge           : ", chrg, "   unpaired: ", uhf, &
         "   electrons: ", nelec
      write (output_unit, "(a,a)") "   method           : ", method()
      write (output_unit, "(a,a)") "   basis            : ", basis_path
      if (rimp2) write (output_unit, "(a,a)") "   auxiliary basis  : ", aux_path
      if (rimp2 .and. .not. cli%is_set("cart")) &
         write (output_unit, "(a)") "   AO space         : Cartesian (RI-MP2 does not support spherical AOs yet)"
      write (output_unit, "(a,a)") "   guess            : ", guess
      write (output_unit, "(a,es8.1,a,es8.1,a,es8.1)") "   convergence      : dE ", econv, &
         "   |FDS-SDF| ", gconv, "   screening ", thresh
      if (len(func) > 0) write (output_unit, "(a,i0)") "   grid level       : ", grid
   end subroutine banner

   subroutine formula()
      !! Hill order: with carbon, C then H then the rest alphabetically;
      !! without, everything alphabetically.
      integer :: k, n, order(118)
      n = 0
      if (any(z == 6)) then
         n = n + 1; order(n) = 6
         if (any(z == 1)) then
            n = n + 1; order(n) = 1
         end if
      end if
      k = n
      do ia = 1, 118
         if (any(z == ia) .and. .not. any(order(:k) == ia)) then
            n = n + 1; order(n) = ia
         end if
      end do
      call sort_by_symbol(order(k + 1:n))
      do k = 1, n
         write (output_unit, "(a)", advance="no") trim(element_symbol(order(k)))
         if (count(z == order(k)) > 1) write (output_unit, "(i0)", advance="no") count(z == order(k))
      end do
   end subroutine formula

   subroutine sort_by_symbol(a)
      integer, intent(inout) :: a(:)
      integer :: i, j, t
      do i = 2, size(a)
         t = a(i); j = i - 1
         do while (j >= 1)
            if (element_symbol(a(j)) <= element_symbol(t)) exit
            a(j + 1) = a(j); j = j - 1
         end do
         a(j + 1) = t
      end do
   end subroutine sort_by_symbol

   subroutine summary()
      character(len=*), parameter :: bar = "          ---------------------------------------------------"
      write (output_unit, "(a)") ""
      write (output_unit, "(a,i0,a,a)") "   AO basis         : ", nao, " functions, ", &
         merge("spherical", "Cartesian", sph /= 0)
      write (output_unit, "(a,i0,a)") "   SCF              : ", niter, merge(" iterations, converged    ", &
                                                                            " iterations, NOT converged", conv /= 0)
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") bar
      call line("TOTAL ENERGY", etot)
      if (rimp2 .and. conv /= 0) then
         call line("MP2 CORRELATION", eos + ess)
         call line("TOTAL MP2 ENERGY", etot + eos + ess)
      end if
      write (output_unit, "(a)") bar
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") "   ::::::::::::::::::::::::::::::::::::::::::::::::::::::"
      write (output_unit, "(a)") "   ::"//repeat(" ", 21)//"SUMMARY"//repeat(" ", 22)//"::"
      write (output_unit, "(a)") "   ::::::::::::::::::::::::::::::::::::::::::::::::::::::"
      call row("total energy", etot)
      call row("nuclear repulsion", enuc)
      call row("one-electron energy", eone)
      call row("two-electron energy", etwo)
      if (len(func) > 0) call row("exchange-correlation", exc)
      if (rimp2 .and. conv /= 0) then
         call row("MP2 opposite spin", eos)
         call row("MP2 same spin", ess)
      end if
      write (output_unit, "(a)") "   ::::::::::::::::::::::::::::::::::::::::::::::::::::::"
      write (output_unit, "(a,f12.3,a)") "   wall time        : ", real(t1 - t0, dp)/real(rate, dp), " s"
   end subroutine summary

   subroutine line(label, e)
      character(len=*), intent(in) :: label
      real(c_double), intent(in) :: e
      write (output_unit, "(a,a20,f22.12,a)") "          | ", label, e, " Eh    |"
   end subroutine line

   subroutine row(label, e)
      character(len=*), intent(in) :: label
      real(c_double), intent(in) :: e
      write (output_unit, "(a,a22,f22.12,a)") "   :: ", label, e, " Eh  ::"
   end subroutine row

   pure function json_escape(t) result(r)
      !! `t` as the inside of a JSON string: backslash and quote escaped.
      character(len=*), intent(in) :: t
      character(len=:), allocatable :: r
      integer :: i
      r = ""
      do i = 1, len(t)
         ! achar(92), not a literal: some compilers read a backslash in a
         ! string as an escape unless told otherwise.
         if (t(i:i) == '"' .or. t(i:i) == achar(92)) then
            r = r//achar(92)//t(i:i)
         else
            r = r//t(i:i)
         end if
      end do
   end function json_escape

   subroutine write_json(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      open (newunit=u, file=path, status="replace", action="write", iostat=ios)
      if (ios /= 0) call die("cannot write "//path)
      write (u, "(a)") "{"
      write (u, "(a)") '  "program": "terco",'
      write (u, "(a)") '  "version": "'//VERSION//'",'
      write (u, "(a)") '  "geometry": "'//json_escape(geom)//'",'
      write (u, "(a)") '  "method": "'//json_escape(method())//'",'
      write (u, "(a)") '  "basis": "'//json_escape(basis_path)//'",'
      write (u, "(a)") '  "guess": "'//json_escape(guess)//'",'
      write (u, "(a,i0,a)") '  "charge": ', chrg, ','
      write (u, "(a,i0,a)") '  "unpaired electrons": ', uhf, ','
      write (u, "(a,i0,a)") '  "number of basis functions": ', nao, ','
      write (u, "(a,a,a)") '  "spherical": ', trim(merge("true ", "false", sph /= 0)), ','
      write (u, "(a,i0,a)") '  "iterations": ', niter, ','
      write (u, "(a,a,a)") '  "converged": ', trim(merge("true ", "false", conv /= 0)), ','
      write (u, "(a,es23.15e3,a)") '  "nuclear repulsion energy": ', enuc, ','
      write (u, "(a,es23.15e3,a)") '  "one-electron energy": ', eone, ','
      write (u, "(a,es23.15e3,a)") '  "two-electron energy": ', etwo, ','
      write (u, "(a,es23.15e3,a)") '  "exchange-correlation energy": ', exc, ','
      if (rimp2 .and. conv /= 0) then
         write (u, "(a,es23.15e3,a)") '  "mp2 opposite-spin energy": ', eos, ','
         write (u, "(a,es23.15e3,a)") '  "mp2 same-spin energy": ', ess, ','
      end if
      write (u, "(a,es23.15e3)") '  "total energy": ', etot
      write (u, "(a)") "}"
      close (u)
   end subroutine write_json

end program terco_main
