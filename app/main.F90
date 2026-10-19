module terco_run
   !! Everything `terco_main` (below) used to keep as internal (`contains`)
   !! procedures, moved out to genuine module procedures that take the
   !! run's state as an explicit argument (`run_t`) instead of reaching it
   !! by host association.
   !!
   !! This is not a style preference: an internal procedure ANYWHERE in this
   !! executable -- even one that is never called, and even one with no
   !! state in common with the SCF call at all -- corrupts the two-electron
   !! part of the Fock build the moment `trc_run_scf` runs collectively
   !! (more than one rank, this build's NVHPC 26.5 + this OpenMPI/HPCX):
   !! `conv` and `niter` still look fine, but the two-electron energy comes
   !! back as exactly 0 and the total is wrong, deterministically, on every
   !! rank alike. Measured by elimination: the program (or procedure) that
   !! makes the collective call must have NO sibling internal procedure at
   !! all -- not "none that touch the same state", none, full stop -- and a
   !! genuine module procedure (this file's old `contains` turned into
   !! `terco_run`'s) does not trip it. Reproduced down to an EMPTY, never
   !! invoked `contains subroutine noop(); end subroutine` in the same
   !! program; unaffected by optimisation level, by `-acc`/`-stdpar=gpu` on
   !! the caller, by MPI transport (ob1 vs UCX), or by UCX's memtype cache.
   !! Looks like an NVHPC/OpenMPI codegen interaction specific to this
   !! toolchain; worth a closer look (and a vendor report) on its own.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use pic_types, only: dp, int64
   use pic_error, only: error_t
   use pic_cli, only: cli_t
   use pic_mpi_lib, only: comm_t, abort_comm, pic_mpi_finalize
   use trc_capi, only: trc_message, TRC_OK
   use terco_app, only: element_symbol, f_string
   implicit none
   private
   public :: run_t, VERSION, die, quit, check, api, int_option, itoa, banner, summary, write_json, &
             gradient_summary

#ifndef TERCO_VERSION
#define TERCO_VERSION "unknown"
#endif
   character(len=*), parameter :: VERSION = TERCO_VERSION

   type :: run_t
      !! Every setting and result `terco_main`'s procedures need, passed
      !! explicitly -- see this module's own header comment for why.
      type(cli_t) :: cli
      type(error_t) :: err
      type(c_ptr) :: h
      character(len=:), allocatable :: geom, basis_arg, basis_path, aux_path, jkfit_path, func, guess, title, json
      character(len=:), allocatable :: pcm
      character(len=:), allocatable :: save_scf, restart   !! "" when not asked for
      integer, allocatable :: z(:)
      integer :: natm, chrg, uhf, nelec, mult, grid
      integer(c_int) :: nao, conv, niter, sph, npcm
      real(dp) :: econv, gconv, thresh, eps
      real(c_double) :: etot, enuc, eone, etwo, exc, eos, ess, epcm, qpcm
      logical :: rimp2, df
      logical :: dh = .false.   !! --func names a double hybrid: rimp2 runs even without --rimp2
      real(c_double) :: pt2_cos = 0.0_c_double, pt2_css = 0.0_c_double, etot_dh = 0.0_c_double
      integer(int64) :: t0, t1, rate
      type(comm_t) :: wcomm
      integer :: nranks
      logical :: rank0
      logical :: grad
      real(dp), allocatable :: gradient(:, :)   !! (3, natm), Eh/Bohr, after --grad
   end type run_t

contains

   subroutine die(r, msg)
      !! Every caller of `die` in this file reaches it the same way on
      !! every rank -- the same argv, the same files, the same collective
      !! library status -- so printing on rank 0 alone says it once, not
      !! `nranks` times, and every other rank is about to call `quit` with
      !! the same `code` regardless.
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: msg
      if (r%rank0) write (error_unit, "(a)") "terco: "//msg
      call quit(r, 1)
   end subroutine die

   subroutine quit(r, code)
      !! End the process with exit status `code`, printing nothing more.
      !! `stop code, quiet=.true.` is the standard way to say this, but
      !! nvfortran rejects it, and a plain `stop code` prints the code.
      !!
      !! More than one rank and a non-zero code: `abort_comm` rather than a
      !! clean `pic_mpi_finalize`, because a caller of `die` cannot promise
      !! every OTHER rank reached the same conclusion the same way -- and
      !! `MPI_Abort` tears down the whole job either way, so there is
      !! nothing to wait for. One rank keeps today's path exactly: no pic-mpi
      !! call before it that was not already a no-op there.
      type(run_t), intent(in) :: r
      integer, intent(in) :: code
      interface
         subroutine c_exit(status) bind(c, name="exit")
            import :: c_int
            implicit none
            integer(c_int), value :: status
         end subroutine c_exit
      end interface
      flush (output_unit)
      flush (error_unit)
      if (code /= 0 .and. r%nranks > 1) call abort_comm(r%wcomm, code)
      call pic_mpi_finalize()
      call c_exit(int(code, c_int))
   end subroutine quit

   subroutine check(r, what)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: what
      if (r%err%has_error()) call die(r, what//": "//r%err%get_message())
   end subroutine check

   subroutine api(r, status)
      !! Stop on a failed library call, with the library's own explanation.
      type(run_t), intent(in) :: r
      integer(c_int), intent(in) :: status
      integer, parameter :: MSG_LEN = 1024
      character(kind=c_char) :: buf(MSG_LEN)
      integer(c_int) :: s
      if (status == TRC_OK) return
      s = trc_message(r%h, buf, int(size(buf), c_int))
      call die(r, f_string(buf))
   end subroutine api

   integer function int_option(r, name, file, default) result(v)
      !! An integer option, else the first integer in `file` in the working
      !! directory (xtb's .CHRG/.UHF convention), else `default`.
      type(run_t), intent(inout) :: r
      character(len=*), intent(in) :: name, file
      integer, intent(in) :: default
      integer :: u, ios
      logical :: there
      v = default
      if (r%cli%is_set(name)) then
         call r%cli%get(name, v, r%err)
         call check(r, "--"//name)
         return
      end if
      inquire (file=file, exist=there)
      if (.not. there) return
      open (newunit=u, file=file, status="old", action="read", iostat=ios)
      if (ios /= 0) return
      read (u, *, iostat=ios) v
      close (u)
      if (ios /= 0) call die(r, file//" does not start with an integer")
   end function int_option

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: b
      write (b, "(i0)") i
      s = trim(b)
   end function itoa

   function method(r) result(m)
      type(run_t), intent(in) :: r
      character(len=:), allocatable :: m
      if (len(r%func) == 0) then
         m = merge("RHF", "UHF", r%mult == 1)
      else
         m = merge("RKS", "UKS", r%mult == 1)//"/"//r%func
      end if
   end function method

   subroutine banner(r)
      type(run_t), intent(in) :: r
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") "      -----------------------------------------------------------"
      write (output_unit, "(a)") "     |                  t e r c o   "//VERSION
      write (output_unit, "(a)") "     |   Hartree-Fock and DFT through do concurrent, in Fortran"
      write (output_unit, "(a)") "      -----------------------------------------------------------"
      write (output_unit, "(a)") ""
      write (output_unit, "(a,a)") "   geometry         : ", r%geom
      if (len_trim(r%title) > 0) write (output_unit, "(a,a)") "   title            : ", trim(r%title)
      write (output_unit, "(a,i0,a)", advance="no") "   atoms            : ", r%natm, "  ("
      call formula(r%z, r%natm)
      write (output_unit, "(a)") ")"
      write (output_unit, "(a,i0,a,i0,a,i0)") "   charge           : ", r%chrg, "   unpaired: ", r%uhf, &
         "   electrons: ", r%nelec
      if (r%nranks > 1) write (output_unit, "(a,i0,a)") "   ranks            : ", r%nranks, " (one GPU each)"
      write (output_unit, "(a,a)") "   method           : ", method(r)
      write (output_unit, "(a,a)") "   basis            : ", r%basis_path
      if (r%rimp2) write (output_unit, "(a,a)") "   auxiliary basis  : ", r%aux_path
      ! The AO space itself (spherical or Cartesian) is not known until the
      ! basis is loaded, below -- `summary`'s own "AO basis" line reports
      ! it, for every run, once it is.
      if (len(r%restart) > 0) then
         write (output_unit, "(a,a,a)") "   SCF              : read from ", r%restart, ", not run"
      else
         write (output_unit, "(a,a)") "   guess            : ", r%guess
      end if
      write (output_unit, "(a,es8.1,a,es8.1,a,es8.1)") "   convergence      : dE ", r%econv, &
         "   |FDS-SDF| ", r%gconv, "   screening ", r%thresh
      if (len(r%func) > 0) write (output_unit, "(a,i0)") "   grid level       : ", r%grid
      if (len(r%pcm) > 0) write (output_unit, "(a,a,a,f0.4)") "   solvation        : ", r%pcm, ", eps ", r%eps
   end subroutine banner

   subroutine formula(z, natm)
      !! Hill order: with carbon, C then H then the rest alphabetically;
      !! without, everything alphabetically.
      !!
      !! Counted once into a per-element histogram with plain loops. The
      !! previous `any(z == ia)` test came back true for every element on
      !! nvfortran/aarch64 (a GH200, 1155-atom rna9), which printed all 118
      !! symbols while `count` still gave the right numbers.
      integer, intent(in) :: natm
      integer, intent(in) :: z(:)
      integer, parameter :: NELEM = 118
      integer :: ia, k, n
      integer :: cnt(NELEM), order(NELEM)
      cnt = 0
      do ia = 1, natm
         if (z(ia) >= 1 .and. z(ia) <= NELEM) cnt(z(ia)) = cnt(z(ia)) + 1
      end do
      n = 0
      if (cnt(6) > 0) then
         n = n + 1
         order(n) = 6
         if (cnt(1) > 0) then
            n = n + 1
            order(n) = 1
         end if
      end if
      k = n
      do ia = 1, NELEM
         if (cnt(ia) == 0) cycle
         if (k > 0 .and. (ia == 6 .or. (ia == 1 .and. k == 2))) cycle
         n = n + 1
         order(n) = ia
      end do
      call sort_by_symbol(order(k + 1:n))
      do ia = 1, n
         write (output_unit, "(a)", advance="no") trim(element_symbol(order(ia)))
         if (cnt(order(ia)) > 1) write (output_unit, "(i0)", advance="no") cnt(order(ia))
      end do
   end subroutine formula

   subroutine sort_by_symbol(a)
      integer, intent(inout) :: a(:)
      integer :: i, j, t
      do i = 2, size(a)
         t = a(i)
         j = i - 1
         do while (j >= 1)
            if (element_symbol(a(j)) <= element_symbol(t)) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = t
      end do
   end subroutine sort_by_symbol

   subroutine summary(r)
      type(run_t), intent(in) :: r
      character(len=*), parameter :: bar = "          ---------------------------------------------------"
      write (output_unit, "(a)") ""
      write (output_unit, "(a,i0,a,a)") "   AO basis         : ", r%nao, " functions, ", &
         merge("spherical", "Cartesian", r%sph /= 0)
      write (output_unit, "(a,i0,a)", advance="no") "   SCF              : ", r%niter, &
         trim(merge(" iterations, converged    ", " iterations, NOT converged", r%conv /= 0))
      if (len(r%restart) > 0) then
         write (output_unit, "(a,a,a)") " (read from ", r%restart, ")"
      else
         write (output_unit, "(a)") ""
      end if
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") bar
      ! "TOTAL ENERGY"/"total energy" means the SAME thing in text and
      ! JSON: the double hybrid's weighted total when one ran, the SCF
      ! total otherwise -- JSON already worked this way; this merge
      ! makes the text banner/summary agree with it, rather than
      ! renaming either key (every external reader of terco's output,
      ! metalquicha included, looks for "total energy"/"TOTAL ENERGY"
      ! by that name). "SCF ENERGY"/"scf energy" is new: always
      ! present, always the SCF-level value alone, so a caller wanting
      ! that specifically no longer has to infer it from context.
      call line("TOTAL ENERGY", merge(r%etot_dh, r%etot, r%dh .and. r%conv /= 0))
      call line("SCF ENERGY", r%etot)
      if (r%rimp2 .and. r%conv /= 0) then
         if (r%dh) then
            call line("PT2 CORRELATION", r%eos + r%ess)
            call line("DOUBLE HYBRID TOTAL", r%etot_dh)
         else
            call line("MP2 CORRELATION", r%eos + r%ess)
            call line("TOTAL MP2 ENERGY", r%etot + r%eos + r%ess)
         end if
      end if
      write (output_unit, "(a)") bar
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") "   ::::::::::::::::::::::::::::::::::::::::::::::::::::::"
      write (output_unit, "(a)") "   ::"//repeat(" ", 21)//"SUMMARY"//repeat(" ", 22)//"::"
      write (output_unit, "(a)") "   ::::::::::::::::::::::::::::::::::::::::::::::::::::::"
      call row("total energy", merge(r%etot_dh, r%etot, r%dh .and. r%conv /= 0))
      call row("scf energy", r%etot)
      call row("nuclear repulsion", r%enuc)
      call row("one-electron energy", r%eone)
      call row("two-electron energy", r%etwo)
      if (len(r%func) > 0) call row("exchange-correlation", r%exc)
      if (len(r%pcm) > 0) call row("solvation (E_diel)", r%epcm)
      if (r%rimp2 .and. r%conv /= 0) then
         call row(merge("PT2 opposite spin", "MP2 opposite spin", r%dh), r%eos)
         call row(merge("PT2 same spin    ", "MP2 same spin    ", r%dh), r%ess)
         if (r%dh) call row("double hybrid total", r%etot_dh)
      end if
      write (output_unit, "(a)") "   ::::::::::::::::::::::::::::::::::::::::::::::::::::::"
      write (output_unit, "(a,f12.3,a)") "   wall time        : ", real(r%t1 - r%t0, dp)/real(r%rate, dp), " s"
   end subroutine summary

   subroutine gradient_summary(r)
      !! Per-atom analytic nuclear gradient (Eh/Bohr), plus the RMS and the
      !! largest per-atom norm -- the two scalars an optimiser's own
      !! convergence check usually gates on. `r%gradient` is already
      !! filled by the time this is called (terco_main's own job, right
      !! after `trc_gradient`).
      type(run_t), intent(in) :: r
      real(dp) :: rms, gmax
      integer :: ia
      write (output_unit, "(a)") ""
      write (output_unit, "(a)") "   gradient (Eh/Bohr)"
      write (output_unit, "(a)") "   ----------------------------------------------------------"
      do ia = 1, r%natm
         write (output_unit, "(a,i4,2x,a2,3f16.8)") "   ", ia, element_symbol(r%z(ia)), r%gradient(:, ia)
      end do
      rms = sqrt(sum(r%gradient**2)/real(3*r%natm, dp))
      gmax = sqrt(maxval(sum(r%gradient**2, dim=1)))
      write (output_unit, "(a)") "   ----------------------------------------------------------"
      write (output_unit, "(a,es12.4)") "   RMS gradient     : ", rms
      write (output_unit, "(a,es12.4)") "   max |gradient|   : ", gmax
   end subroutine gradient_summary

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

   pure function json_escape(t) result(esc)
      !! `t` as the inside of a JSON string: backslash and quote escaped.
      character(len=*), intent(in) :: t
      character(len=:), allocatable :: esc
      integer :: i
      esc = ""
      do i = 1, len(t)
         ! achar(92), not a literal: some compilers read a backslash in a
         ! string as an escape unless told otherwise.
         if (t(i:i) == '"' .or. t(i:i) == achar(92)) then
            esc = esc//achar(92)//t(i:i)
         else
            esc = esc//t(i:i)
         end if
      end do
   end function json_escape

   subroutine write_json(r, path)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      integer :: u, ios
      open (newunit=u, file=path, status="replace", action="write", iostat=ios)
      if (ios /= 0) call die(r, "cannot write "//path)
      write (u, "(a)") "{"
      write (u, "(a)") '  "program": "terco",'
      write (u, "(a)") '  "version": "'//VERSION//'",'
      write (u, "(a)") '  "geometry": "'//json_escape(r%geom)//'",'
      write (u, "(a)") '  "method": "'//json_escape(method(r))//'",'
      write (u, "(a)") '  "basis": "'//json_escape(r%basis_path)//'",'
      write (u, "(a)") '  "guess": "'//json_escape(r%guess)//'",'
      write (u, "(a,i0,a)") '  "charge": ', r%chrg, ","
      write (u, "(a,i0,a)") '  "unpaired electrons": ', r%uhf, ","
      write (u, "(a,i0,a)") '  "number of basis functions": ', r%nao, ","
      write (u, "(a,a,a)") '  "spherical": ', trim(merge("true ", "false", r%sph /= 0)), ","
      write (u, "(a,i0,a)") '  "iterations": ', r%niter, ","
      write (u, "(a,a,a)") '  "converged": ', trim(merge("true ", "false", r%conv /= 0)), ","
      write (u, "(a,es23.15e3,a)") '  "nuclear repulsion energy": ', r%enuc, ","
      write (u, "(a,es23.15e3,a)") '  "one-electron energy": ', r%eone, ","
      write (u, "(a,es23.15e3,a)") '  "two-electron energy": ', r%etwo, ","
      write (u, "(a,es23.15e3,a)") '  "exchange-correlation energy": ', r%exc, ","
      if (len(r%pcm) > 0) then
         write (u, "(a,a,a)") '  "solvation model": "', r%pcm, '",'
         write (u, "(a,es23.15e3,a)") '  "dielectric constant": ', r%eps, ","
         write (u, "(a,es23.15e3,a)") '  "solvation energy": ', r%epcm, ","
         write (u, "(a,i0,a)") '  "solvation surface points": ', r%npcm, ","
      end if
      if (r%rimp2 .and. r%conv /= 0) then
         write (u, "(a,es23.15e3,a)") '  "mp2 opposite-spin energy": ', r%eos, ","
         write (u, "(a,es23.15e3,a)") '  "mp2 same-spin energy": ', r%ess, ","
         if (r%dh) then
            write (u, "(a,es23.15e3,a)") '  "pt2 opposite-spin fraction": ', r%pt2_cos, ","
            write (u, "(a,es23.15e3,a)") '  "pt2 same-spin fraction": ', r%pt2_css, ","
            write (u, "(a,es23.15e3,a)") '  "scf energy": ', r%etot, ","
            write (u, "(a,es23.15e3)") '  "total energy": ', r%etot_dh
            write (u, "(a)") "}"
            close (u)
            return
         end if
      end if
      write (u, "(a,es23.15e3,a)") '  "scf energy": ', r%etot, ","
      if (r%grad .and. allocated(r%gradient)) then
         write (u, "(a,es23.15e3,a)") '  "total energy": ', r%etot, ","
         call write_json_gradient(u, r)
      else
         write (u, "(a,es23.15e3)") '  "total energy": ', r%etot
      end if
      write (u, "(a)") "}"
      close (u)
   end subroutine write_json

   subroutine write_json_gradient(u, r)
      !! The analytic gradient as a JSON array of [x, y, z] triples, one per
      !! atom in `r%gradient`'s own (and `trc_set_molecule`'s) atom order,
      !! plus the RMS and max per-atom norm `gradient_summary` also prints.
      integer, intent(in) :: u
      type(run_t), intent(in) :: r
      integer :: ia
      real(dp) :: rms, gmax
      character(len=2) :: closing
      write (u, "(a)") '  "gradient": ['
      do ia = 1, r%natm
         closing = merge("] ", "],", ia == r%natm)
         write (u, "(a,es23.15e3,a,es23.15e3,a,es23.15e3,a)") "    [", r%gradient(1, ia), ", ", &
            r%gradient(2, ia), ", ", r%gradient(3, ia), trim(closing)
      end do
      write (u, "(a)") "  ],"
      rms = sqrt(sum(r%gradient**2)/real(3*r%natm, dp))
      gmax = sqrt(maxval(sum(r%gradient**2, dim=1)))
      write (u, "(a,es23.15e3,a)") '  "gradient rms": ', rms, ","
      write (u, "(a,es23.15e3)") '  "gradient max": ', gmax
   end subroutine write_json_gradient

end module terco_run

program terco_main
   !! The `terco` command: an SCF on an xyz file, from the shell.
   !!
   !! `terco molecule.xyz [options]`, in the spirit of xtb's command line.
   !!
   !! Everything below the argument handling is the library's context API
   !! (trc_capi): the guesses, the SCF, RI-MP2 and the energies are exactly
   !! what a C or Python host gets. `terco --help` lists the options.
   !!
   !! This program itself has NO internal (`contains`) procedures -- see
   !! terco_run's header comment for why that is load-bearing, not a style
   !! choice, the moment more than one rank is involved.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init
   use pic_knowledge, only: get_knowledge
   use pic_profiler, only: profiler_init
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_aux_json, trc_set_method, trc_set_nlc_grid, trc_set_convergence, trc_set_screening, &
                       trc_set_incremental, trc_set_df, trc_set_jkfit_json, &
                       trc_set_guess, trc_set_verbose, trc_set_rimp2, trc_set_force_cartesian, trc_set_comm, &
                       trc_run_scf, trc_run_rimp2, trc_nao, trc_is_spherical, trc_nspin, trc_energy, &
                       trc_energy_parts, trc_converged, trc_iterations, trc_rimp2_energy, &
                       trc_set_pcm, trc_pcm_energy, trc_save_scf, trc_load_scf, &
                       trc_run_gradient, trc_gradient, &
                       TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_CORE, TRC_GUESS_GWH, TRC_GUESS_SAD, &
                       TRC_GUESS_SADQ, TRC_GUESS_SAP
   use trc_api, only: trc_bind_device
   use trc_xc_functional, only: XC_AVAILABLE_FUNCTIONALS, trc_xc_functional_t, xc_functional_by_name, &
                                xc_is_double_hybrid
   use terco_app, only: read_xyz, resolve_basis, c_string, set_env, lower
   use terco_run, only: run_t, VERSION, die, quit, check, api, int_option, itoa, banner, summary, write_json, &
                       gradient_summary
   implicit none

   type(run_t) :: r
   character(len=:), allocatable :: sval
   real(dp), allocatable :: xyz(:, :)
   real(c_double), allocatable :: gradbuf(:)
   integer :: nlc_grid, maxiter, nfrozen, guess_kind, dev, aux_block
   integer(c_int) :: st, nspin
   logical :: cart

   ! MPI, in every build: pic-mpi's single-rank backend (TERCO_ENABLE_MPI=OFF)
   ! makes this a no-op, so the rest of the program need not branch on it.
   ! `r%wcomm` and `r%nranks` answer everything the rest of this file needs
   ! to know about ranks; `r%rank0` is what gates every print below so N
   ! ranks produce ONE copy of the output, not N.
   call pic_mpi_init()
   ! pic's region profiler, which is what puts the library's named NVTX
   ! ranges on an Nsight Systems timeline (a TERCO_NVTX build); without NVTX
   ! it only keeps host wall times nobody reads, so it is always on.
   call profiler_init()
   r%wcomm = comm_world()
   r%nranks = r%wcomm%size()
   r%rank0 = r%wcomm%leader()

   call system_clock(r%t0, r%rate)

   ! ---- the interface ---------------------------------------------------
   call r%cli%set_program("terco", "Hartree-Fock and Kohn-Sham DFT, in standard Fortran")
   call r%cli%add_positional("geometry", "molecule in xyz format, Angstrom", required=.false.)
   call r%cli%add_option("basis", "basis set: a Basis Set Exchange JSON file, or a name looked up "// &
                         "in $TERCO_BASIS_PATH and then the installed basis sets (def2-svp)", short="b", &
                         default="def2-svp")
   call r%cli%add_option("chrg", "molecular charge (default: a .CHRG file here, else 0)", short="c")
   call r%cli%add_option("uhf", "number of unpaired electrons (default: a .UHF file here, else 0 "// &
                         "or 1 by electron count)", short="u")
   call r%cli%add_option("func", "hf (the default), or a functional: "//XC_AVAILABLE_FUNCTIONALS, &
                         short="f", default="hf")
   call r%cli%add_option("guess", "sad, sadq, sap, gwh or core (default: sad; sadq for an ion)", short="g")
   call r%cli%add_flag("cart", "Cartesian AOs even when the basis is defined spherical")
   call r%cli%add_option("grid", "DFT integration grid level (3)", default="3")
   call r%cli%add_option("nlc-grid", "grid level for VV10 nonlocal correlation, e.g. wb97x-v (1)", default="1")
   call r%cli%add_option("iterations", "maximum SCF iterations (100)", default="100")
   call r%cli%add_option("econv", "energy convergence, Eh (1e-10)", default="1e-10")
   call r%cli%add_option("gconv", "convergence on the DIIS error norm |FDS - SDF| (1e-5)", default="1e-5")
   call r%cli%add_option("thresh", "integral screening threshold (1e-12)", default="1e-12")
   call r%cli%add_flag("incremental", "build later Fock matrices from the density difference (closed shell)")
   call r%cli%add_flag("df", "density-fitted J and K instead of the four-centre engine (restricted, "// &
                       "Cartesian only; fits both, never DF-J with 4-centre K)")
   call r%cli%add_option("jkfit", "fitting basis for --df (def2-universal-jkfit)", default="def2-universal-jkfit")
   call r%cli%add_option("pcm", "polarisable continuum solvation: cpcm or iefpcm (pyscf's SWIG "// &
                         "cavity, 302 points per sphere, radii 1.2 x Bondi); needs --eps")
   call r%cli%add_option("eps", "solvent dielectric constant for --pcm, e.g. 78.3553 for water")
   call r%cli%add_flag("rimp2", "RI-MP2 correlation energy after the SCF (restricted or "// &
                       "unrestricted; unrestricted refuses more than one rank)")
   call r%cli%add_option("aux", "auxiliary basis for --rimp2 (cc-pvdz-rifit)", default="cc-pvdz-rifit")
   call r%cli%add_option("frozen", "frozen core orbitals for --rimp2; -1, the default, counts them by element", &
                         default="-1")
   call r%cli%add_option("aux-block", "auxiliary functions per three-index block for --rimp2 (default: as many "// &
                         "as keep the block under 2 GB)")
   call r%cli%add_option("save-scf", "write the SCF result (orbitals, density, energies) to this file, "// &
                         "for a later --restart")
   call r%cli%add_option("restart", "read the SCF result from a --save-scf file instead of running the "// &
                         "SCF, e.g. to go straight to --rimp2; the molecule, basis and method must match")
   call r%cli%add_option("json", "also write a JSON summary to this file")
   call r%cli%add_flag("grad", "print the analytic nuclear gradient (Eh/Bohr) after the SCF -- "// &
                       "RHF/UHF, or KS-DFT for an LDA/GGA --func (global and range-separated "// &
                       "hybrids included); refused for a VV10 --func, --rimp2, --df, --pcm, "// &
                       "or more than one rank")
   call r%cli%add_flag("silent", "only the summary: no SCF iterations, timings or guess details", short="s")
   call r%cli%add_flag("version", "print the version and stop")

   call r%cli%parse(r%err)
   if (r%err%has_error()) call die(r, r%err%get_message()//new_line("a")//"try 'terco --help'")
   if (r%cli%help_requested()) then
      if (r%rank0) write (output_unit, "(a)") r%cli%help_text()
      call quit(r, 0)
   end if
   if (r%cli%is_set("version")) then
      if (r%rank0) write (output_unit, "(a)") "terco version "//VERSION
      call quit(r, 0)
   end if
   call r%cli%get("geometry", r%geom, r%err)
   if (r%err%has_error() .or. len_trim(r%geom) == 0) then
     call die(r, "no geometry given: terco molecule.xyz [options]"//new_line("a")//"try 'terco --help'")
   end if

   ! ---- inputs ----------------------------------------------------------
   call read_xyz(r%geom, r%natm, r%z, xyz, r%title, r%err)
   if (r%err%has_error()) call die(r, r%err%get_message())

   call r%cli%get("basis", r%basis_arg, r%err)
   call resolve_basis(r%basis_arg, r%basis_path, r%err)
   if (r%err%has_error()) call die(r, r%err%get_message())

   r%chrg = int_option(r, "chrg", ".CHRG", 0)
   r%nelec = sum(r%z) - r%chrg
   if (r%nelec < 1) call die(r, "charge "//itoa(r%chrg)//" leaves "//itoa(r%nelec)//" electrons")
   r%uhf = int_option(r, "uhf", ".UHF", mod(r%nelec, 2))
   if (r%uhf < 0 .or. r%uhf > r%nelec .or. mod(r%nelec - r%uhf, 2) /= 0) then
     call die(r, "--uhf "//itoa(r%uhf)//" is impossible with "//itoa(r%nelec)//" electrons: the unpaired "// &
               "count must be between 0 and the electron count, with the same parity")
   end if
   r%mult = r%uhf + 1

   call r%cli%get("func", sval, r%err)
   r%func = lower(trim(sval))
   if (r%func == "hf" .or. r%func == "rhf" .or. r%func == "uhf") r%func = ""

   ! A double hybrid's PT2 part is always RI-MP2 on the converged KS
   ! orbitals (no canonical, non-fitted MP2 path exists in terco at all):
   ! --func b2plyp runs it without needing --rimp2 too, with --aux's usual
   ! default (cc-pvdz-rifit) unless overridden below.
   if (len(r%func) > 0) then
      block
         type(trc_xc_functional_t) :: func_spec
         if (xc_functional_by_name(r%func, func_spec)) r%dh = xc_is_double_hybrid(func_spec)
         if (r%dh) then
            r%pt2_cos = real(func_spec%pt2_cos, c_double)
            r%pt2_css = real(func_spec%pt2_css, c_double)
         end if
      end block
   end if

   r%guess = "sad"
   if (r%chrg /= 0) r%guess = "sadq"
   if (r%cli%is_set("guess")) then
      call r%cli%get("guess", sval, r%err)
      r%guess = lower(trim(sval))
   end if
   select case (r%guess)
   case ("sad")
   guess_kind = TRC_GUESS_SAD
   case ("sadq")
   guess_kind = TRC_GUESS_SADQ
   case ("sap")
   guess_kind = TRC_GUESS_SAP
   case ("gwh")
   guess_kind = TRC_GUESS_GWH
   case ("core", "hcore")
   guess_kind = TRC_GUESS_CORE
   case default
      call die(r, "unknown guess '"//r%guess//"': use sad, sadq, sap, gwh or core")
   end select

   call r%cli%get("grid", r%grid, r%err)
   call check(r, "--grid")
   call r%cli%get("nlc-grid", nlc_grid, r%err)
   call check(r, "--nlc-grid")
   call r%cli%get("iterations", maxiter, r%err)
   call check(r, "--iterations")
   call r%cli%get("econv", r%econv, r%err)
   call check(r, "--econv")
   call r%cli%get("gconv", r%gconv, r%err)
   call check(r, "--gconv")
   call r%cli%get("thresh", r%thresh, r%err)
   call check(r, "--thresh")
   call r%cli%get("frozen", nfrozen, r%err)
   call check(r, "--frozen")
   r%rimp2 = r%cli%is_set("rimp2") .or. r%dh
   r%df = r%cli%is_set("df")
   r%grad = r%cli%is_set("grad")
   cart = r%cli%is_set("cart")
   r%pcm = ""
   r%eps = 0.0_dp
   if (r%cli%is_set("pcm")) then
      call r%cli%get("pcm", sval, r%err)
      r%pcm = lower(trim(sval))
      if (r%pcm /= "cpcm" .and. r%pcm /= "iefpcm") call die(r, "unknown --pcm '"//r%pcm//"': use cpcm or iefpcm")
      if (.not. r%cli%is_set("eps")) call die(r, "--pcm needs the solvent's dielectric constant: "// &
                                              "--eps 78.3553 for water; nothing is assumed")
      call r%cli%get("eps", r%eps, r%err)
      call check(r, "--eps")
   else if (r%cli%is_set("eps")) then
      call die(r, "--eps is the dielectric constant for --pcm, which was not given")
   end if
   if (r%rimp2) then
      ! RI-MP2 (restricted or unrestricted, UHF/UKS) now follows the
      ! basis's own AO space, spherical or Cartesian, same as a plain
      ! SCF -- trc_run_rimp2 folds the MO coefficients into the Cartesian
      ! space its three-centre tensor runs in itself, see trc_capi.F90's
      ! own header. `--cart` is still the explicit override. An open-shell
      ! --rimp2 refuses more than one rank, deeper in (trc_run_rimp2) --
      ! that needs the communicator this parses later, not knowable yet.
      call r%cli%get("aux", sval, r%err)
      call resolve_basis(sval, r%aux_path, r%err)
      if (r%err%has_error()) call die(r, "--aux: "//r%err%get_message())
   end if
   if (r%df) then
      if (r%cli%is_set("incremental")) call die(r, "--df with --incremental is refused -- a "// &
                                               "density-fitted build is a full build every "// &
                                               "iteration, and the increment buys nothing")
      ! `--df` follows the basis's own AO space too, same as --rimp2 above
      ! (trc_scf_df.F90's df_jk_setup only ever touches the Cartesian
      ! d_c/gmat_c trc_scf_fock.F90's to_cart/to_sph fold at the boundary;
      ! it never reads b%spherical itself). `--cart` is still the explicit
      ! override.
      call r%cli%get("jkfit", sval, r%err)
      call resolve_basis(sval, r%jkfit_path, r%err)
      if (r%err%has_error()) call die(r, "--jkfit: "//r%err%get_message())
   end if
   r%json = ""
   if (r%cli%is_set("json")) call r%cli%get("json", r%json, r%err)
   r%save_scf = ""
   if (r%cli%is_set("save-scf")) call r%cli%get("save-scf", r%save_scf, r%err)
   r%restart = ""
   if (r%cli%is_set("restart")) call r%cli%get("restart", r%restart, r%err)

   ! SAP's fitted atomic potentials ship as a basis-format file; find it the
   ! way a basis is found, unless the caller has already said where it is.
   if (guess_kind == TRC_GUESS_SAP) then
      block
         character(len=1024) :: e
         integer :: n, s
         call get_environment_variable("TRC_SAP_BASIS", e, n, s)
         if (s /= 0 .or. n == 0) then
            call resolve_basis("sap_helfem_large", sval, r%err)
            if (r%err%has_error()) call die(r, "--guess sap needs its fitted potentials: "//r%err%get_message())
            call set_env("TRC_SAP_BASIS", sval)
         end if
      end block
   end if

   if (r%rank0) call banner(r)

   ! ---- the calculation -------------------------------------------------
   ! One process, one device: the first, or (more than one rank) this
   ! rank's own, by NODE-LOCAL rank -- the same split trc_build_comm does
   ! inside the library, done here too because device binding must happen
   ! before anything reaches the device, and the guess the library builds
   ! starts doing that before it gets to call trc_build_comm itself. A host
   ! build has nothing to bind either way.
   if (r%nranks > 1) then
      block
         type(comm_t) :: node_comm
         node_comm = r%wcomm%split()
         dev = trc_bind_device(node_comm%rank())
         call node_comm%finalize()
      end block
   else
      dev = trc_bind_device(0)
   end if
   st = trc_create(r%h)
   ! Every rank must make this call, with the same handle, before trc_run_scf
   ! -- trc_set_comm's own contract. One rank keeps today's path exactly:
   ! the context's fcomm defaults to -1 and is never touched.
   if (r%nranks > 1) call api(r, trc_set_comm(r%h, int(r%wcomm%fortran_handle(), c_int)))
   call api(r, trc_set_verbose(r%h, merge(0_c_int, 1_c_int, r%cli%is_set("silent"))))
   call api(r, trc_set_molecule(r%h, int(r%natm, c_int), real(r%z, c_double), xyz, int(r%chrg, c_int), &
                                int(r%mult, c_int)))
   call api(r, trc_set_force_cartesian(r%h, merge(1_c_int, 0_c_int, cart)))
   call api(r, trc_set_basis_json(r%h, c_string(r%basis_path)))
   if (r%rimp2) call api(r, trc_set_aux_json(r%h, c_string(r%aux_path)))
   if (r%df) then
      call api(r, trc_set_jkfit_json(r%h, c_string(r%jkfit_path)))
      call api(r, trc_set_df(r%h, 1_c_int))
   end if
   call api(r, trc_set_method(r%h, c_string(r%func), int(r%grid, c_int)))
   call api(r, trc_set_nlc_grid(r%h, int(nlc_grid, c_int)))
   if (len(r%pcm) > 0) call api(r, trc_set_pcm(r%h, c_string(r%pcm), real(r%eps, c_double), 0_c_int, 0.0_c_double))
   call api(r, trc_set_convergence(r%h, real(r%econv, c_double), real(r%gconv, c_double), int(maxiter, c_int), &
                                   10_c_int))
   call api(r, trc_set_screening(r%h, real(r%thresh, c_double)))
   ! On unless asked otherwise: the library's default is off, so a host
   ! code keeps the build it had, but a command-line run gets the speed.
   call api(r, trc_set_incremental(r%h, merge(1_c_int, 0_c_int, r%cli%is_set("incremental"))))
   call api(r, trc_set_guess(r%h, int(guess_kind, c_int), c_null_ptr, 1_c_int))
   if (r%rimp2) then
      ! The three-index tensor (mu nu|P) is nao^2 per auxiliary function and
      ! is held one block of them at a time; the library's own default is
      ! one block of all of them, which at a thousand functions and a
      ! cc-pVDZ-RIFIT basis is 30 GB. 2 GB a block leaves the device for the
      ! transformed tensor, which is the one RI-MP2 cannot split.
      if (r%cli%is_set("aux-block")) then
         call r%cli%get("aux-block", aux_block, r%err)
         call check(r, "--aux-block")
         if (aux_block < 1) call die(r, "--aux-block must be at least 1")
      else
         call api(r, trc_nao(r%h, r%nao))
         aux_block = int(max(16.0_dp, min(real(huge(1), dp), 2.0e9_dp/(8.0_dp*real(r%nao, dp)**2))))
      end if
      call api(r, trc_set_rimp2(r%h, int(nfrozen, c_int), int(aux_block, c_int)))
   end if

   if (len(r%restart) > 0) then
      st = trc_load_scf(r%h, c_string(r%restart))
   else
      st = trc_run_scf(r%h)
   end if
   if (st /= TRC_OK .and. st /= TRC_ERR_NOCONV) call api(r, st)
   ! Rank 0 alone: every rank holds the same result.
   if (len(r%save_scf) > 0 .and. r%rank0) call api(r, trc_save_scf(r%h, c_string(r%save_scf)))
   call api(r, trc_nao(r%h, r%nao))
   call api(r, trc_is_spherical(r%h, r%sph))
   call api(r, trc_nspin(r%h, nspin))
   call api(r, trc_energy(r%h, r%etot))
   call api(r, trc_energy_parts(r%h, r%enuc, r%eone, r%etwo, r%exc))
   call api(r, trc_converged(r%h, r%conv))
   call api(r, trc_iterations(r%h, r%niter))
   call api(r, trc_pcm_energy(r%h, r%epcm, r%qpcm, r%npcm))
   r%eos = 0
   r%ess = 0
   if (r%rimp2 .and. r%conv /= 0) then
      call api(r, trc_run_rimp2(r%h))
      call api(r, trc_rimp2_energy(r%h, r%eos, r%ess))
      ! The double hybrid's own weighted total comes from trc_energy
      ! itself (now that RI-MP2 has run, cx%rimp2_ok is set and its
      ! pt2_cos/pt2_css branch applies), not recomputed here from
      ! r%pt2_cos/r%pt2_css -- those stay around only to report the
      ! fractions themselves (write_json), never to duplicate the
      ! weighting trc_energy already does from the SAME table.
      if (r%dh) call api(r, trc_energy(r%h, r%etot_dh))
   end if
   if (r%grad .and. r%conv /= 0) then
      ! trc_run_gradient's own refusals (a VV10 --func, RI-MP2, --df,
      ! PCM, multi-rank) come back through `api`, which prints the
      ! library's own message and stops -- the same clear-refusal path
      ! every other unsupported combination in this CLI already takes.
      ! A plain LDA/GGA --func (global and range-separated hybrids
      ! included) is NOT refused: trc_run_gradient builds the KS-DFT
      ! gradient at the SCF's own functional and grid.
      call api(r, trc_run_gradient(r%h))
      allocate (gradbuf(3*r%natm), r%gradient(3, r%natm))
      call api(r, trc_gradient(r%h, gradbuf))
      r%gradient = reshape(real(gradbuf, dp), [3, r%natm])
      deallocate (gradbuf)
   end if

   call system_clock(r%t1)
   if (r%rank0) call summary(r)
   if (r%rank0 .and. r%grad .and. r%conv /= 0) call gradient_summary(r)
   if (r%rank0 .and. len(r%json) > 0) call write_json(r, r%json)
   st = trc_destroy(r%h)

   if (r%conv == 0) then
      if (r%rank0) write (error_unit, "(a)") "terco: the SCF did not converge in "//itoa(int(r%niter))//" iterations"
      call quit(r, 2)
   end if
   if (r%rank0) write (output_unit, "(/,a)") " normal termination of terco"
   ! A line of pic's collected lore, as metalquicha prints at the end of a
   ! successful run: rank 0 only, and not under --silent, whose promise is
   ! the summary and nothing else.
   if (r%rank0 .and. .not. r%cli%is_set("silent")) call get_knowledge()
   call quit(r, 0)

end program terco_main
