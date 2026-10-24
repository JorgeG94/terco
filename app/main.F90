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
   use pic_types, only: dp, int64, default_int
   use pic_error, only: error_t, ERROR_VALIDATION
   use pic_cli, only: cli_t
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_mpi_lib, only: comm_t, abort_comm, pic_mpi_finalize
   use json_module, only: json_core, json_value
   use trc_capi, only: trc_message, TRC_OK
   use terco_app, only: element_symbol, f_string, BOHR_PER_ANGSTROM, write_opt_xyz_frame, write_traj_xyz_frames
   use trc_vib, only: trc_vib_log_table, trc_vib_rotor
   use trc_thermo, only: trc_thermo_result_t, trc_thermo_log_summary
   use trc_constants, only: trc_atomic_mass, trc_atomic_mass_average, trc_covalent_radius
   use trc_error, only: mass_error_t => error_t
   implicit none
   private
   public :: run_t, VERSION, die, quit, check, api, int_option, banner, summary, write_json, &
             gradient_summary, hessian_summary, build_masses
   public :: opt_capability_check, build_opt_residues, opt_summary, write_opt_xyz, write_traj_xyz
   public :: opt_geometry_is_linear, write_labelled_xyz, to_lower_alloc

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
      integer :: natm, chrg, uhf, nelec, mult, grid, nlc_grid
      integer(c_int) :: nao, conv, niter, sph, npcm
      real(dp) :: econv, gconv, thresh, eps
      real(c_double) :: etot, enuc, eone, etwo, exc, eos, ess, epcm, qpcm
      logical :: rimp2, df
      logical :: dh = .false.   !! --func names a double hybrid: rimp2 runs even without --rimp2
      logical :: vv10 = .false.  !! --func has VV10 nonlocal correlation, on its own --nlc-grid
      real(c_double) :: pt2_cos = 0.0_c_double, pt2_css = 0.0_c_double, etot_dh = 0.0_c_double
      integer(int64) :: t0, t1, rate
      type(comm_t) :: wcomm
      integer :: nranks
      logical :: rank0
      logical :: grad
      real(dp), allocatable :: gradient(:, :)   !! (3, natm), Eh/Bohr, after --grad
      ! --hessian: the FD Hessian (trc_run_hessian/trc_hessian) and the
      ! vibrational analysis/thermochemistry built on it (trc_vib/
      ! trc_thermo).
      logical :: hessian = .false.
      real(dp) :: hess_step = 0.005_dp, temperature_k = 298.15_dp, pressure_atm = 1.0_dp
      integer :: sym_number = 1
      character(len=:), allocatable :: masses_kind   !! "isotope" or "average"
      real(dp), allocatable :: hess(:, :)            !! (3*natm, 3*natm), Eh/Bohr^2
      real(dp), allocatable :: freq_cm1(:), redmass_amu(:), force_const(:)
      logical :: hess_ok = .false.
      logical :: is_linear = .false.
      integer :: n_tr = 0
      type(trc_thermo_result_t) :: thermo
      ! --opt: geometry optimisation through app/trc_dlfind_bridge.F90's
      ! trc_dlfind_optimize (TERCO_ENABLE_DLFIND=ON) or its stub (OFF, the
      ! default -- always refused with a message). opt_coords_kind/
      ! opt_guess_kind are trc_dlfind_bridge's own TRC_OPT_COORDS_*/
      ! TRC_OPT_GUESS_* enums, resolved from --opt-coords/--opt-guess here
      ! so main.F90 need not import that module's constants under every
      ! build configuration.
      logical :: opt = .false.
      character(len=:), allocatable :: opt_coords_arg, opt_guess_arg
      integer :: opt_coords_kind = 0, opt_guess_kind = 0
      integer :: opt_maxiter = 100
      real(dp) :: opt_gtol = 4.5e-4_dp
      real(dp) :: opt_max_step = -1.0_dp   !! negative = DL-FIND's own default (trc_opt_options_t%max_step)
      logical :: opt_converged = .false.
      integer :: opt_steps = 0, opt_evals = 0, opt_hess_evals = 0
      real(dp) :: opt_gmax = 0.0_dp, opt_grms = 0.0_dp
      integer :: opt_n_traj = 0
      real(dp), allocatable :: opt_traj_energy(:)
      real(dp), allocatable :: opt_final_xyz(:, :)   !! (3, natm), Bohr -- the optimised geometry
      ! Algorithm, target (minimum/saddle), the two saddle methods,
      ! `connect`, and the external-Hessian control -- every one of these
      ! is also one of
      ! trc_dlfind_bridge's own TRC_OPT_* enums (app/trc_dlfind_bridge.F90,
      ! app/trc_dlfind_stub.F90), stored here as the plain integer both
      ! modules define identically, same reason opt_coords_kind/
      ! opt_guess_kind above need no import of that module.
      character(len=:), allocatable :: opt_algorithm_arg, opt_target_arg, opt_saddle_method_arg, &
                                       opt_init_hessian_arg, opt_hessian_update_arg, opt_linesearch_arg
      integer :: opt_algorithm_kind = 3   !! TRC_OPT_ALG_LBFGS
      integer :: opt_target_kind = 0      !! TRC_OPT_TARGET_MINIMUM
      integer :: opt_saddle_method_kind = 0   !! TRC_OPT_SADDLE_PRFO
      logical :: opt_connect = .false.
      logical :: opt_is_linear = .false.
      integer :: opt_hessian_source_kind = 0  !! TRC_OPT_HESS_TERCO
      integer :: opt_hessian_update_kind = -1, opt_max_updates = -1, opt_zero_modes = -1
      real(dp) :: opt_soft_modes = -1.0_dp, opt_scalestep = -1.0_dp
      integer :: opt_linesearch_kind = -1   !! TRC_OPT_LINE_AUTO
      integer :: opt_lbfgs_memory = -1
      real(dp) :: opt_etol = -1.0_dp
      integer :: opt_print_level = -1
      ! --opt-connect's own three structures (the TS and both downhill
      ! minima); populated only when opt_connect is true and the run did
      ! not fail.
      logical :: opt_have_ts = .false., opt_have_min_plus = .false., opt_have_min_minus = .false.
      real(dp) :: opt_ts_energy = 0.0_dp, opt_min_plus_energy = 0.0_dp, opt_min_minus_energy = 0.0_dp
      real(dp) :: opt_ts_gmax = 0.0_dp, opt_min_plus_gmax = 0.0_dp, opt_min_minus_gmax = 0.0_dp
      logical :: opt_ts_converged = .false., opt_min_plus_converged = .false., opt_min_minus_converged = .false.
      real(dp), allocatable :: opt_ts_xyz(:, :), opt_min_plus_xyz(:, :), opt_min_minus_xyz(:, :)
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
      if (r%rank0) call global_logger%error("terco: "//msg)
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

   function to_lower_alloc(s) result(lower_s)
      !! A non-elemental, allocatable-result ASCII lowercaser for every CLI
      !! normalisation site in `terco_main` below (`--func`, `--guess`,
      !! `--masses`, `--pcm`, the `--opt-*` string options, and the geometry
      !! file-extension check). It replaces pic's `to_lower`
      !! (`pic_ascii.f90`): nvfortran 26.5 and 26.9, host and GPU builds
      !! alike, with `terco_main`'s real `-fast -O3` flags, miscompile
      !! pic's `to_lower` -- `pure elemental`, automatic-length (not
      !! allocatable) result -- at those call sites. The result comes back
      !! with the right length but blank (all-space) content, so e.g.
      !! `--func pbe` silently falls back to a plain SCF with no exchange-
      !! correlation term. A non-elemental function with an allocatable
      !! result, called the same way, is unaffected; see
      !! `.claude/context/gpu-pitfalls.md` for the trap. `terco_main` has no
      !! internal procedures (see this module's header comment), so this
      !! lives here instead.
      character(len=*), intent(in) :: s
      character(len=:), allocatable :: lower_s
      integer :: i, c
      allocate (character(len=len(s)) :: lower_s)
      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            lower_s(i:i) = achar(c + 32)
         else
            lower_s(i:i) = s(i:i)
         end if
      end do
   end function to_lower_alloc

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
      call global_logger%info("")
      call global_logger%info("      -----------------------------------------------------------")
      call global_logger%info("     |                  t e r c o   "//VERSION)
      call global_logger%info("     |         GPU accelerated Quantum Chemistry ")
      call global_logger%info("     |     (c) 2026, Jorge Luis Galvez Vallejo, NCI  ")
      call global_logger%info("     |               Powered by Fortran        ")
      call global_logger%info("      -----------------------------------------------------------")
      call global_logger%info("")
      call global_logger%info("   geometry         : "//r%geom)
      if (len_trim(r%title) > 0) call global_logger%info("   title            : "//trim(r%title))
      call global_logger%info("   atoms            : "//to_string(r%natm)//"  ("//formula(r%z, r%natm)//")")
      block
         character(len=1024) :: logline
         write (logline, "(a,i0,a,i0,a,i0)") "   charge           : ", r%chrg, "   unpaired: ", r%uhf, &
            "   electrons: ", r%nelec
         call global_logger%info(trim(logline))
      end block
      if (r%nranks > 1) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a)") "   ranks            : ", r%nranks, " (one GPU each)"
            call global_logger%info(trim(logline))
         end block
      end if
      call global_logger%info("   method           : "//method(r))
      call global_logger%info("   basis            : "//r%basis_path)
      if (r%rimp2) call global_logger%info("   auxiliary basis  : "//r%aux_path)
      ! The AO space itself (spherical or Cartesian) is not known until the
      ! basis is loaded, below -- `summary`'s own "AO basis" line reports
      ! it, for every run, once it is.
      if (len(r%restart) > 0) then
         call global_logger%info("   SCF              : read from "//r%restart//", not run")
      else
         call global_logger%info("   guess            : "//r%guess)
      end if
      block
         character(len=1024) :: logline
         write (logline, "(a,es8.1,a,es8.1,a,es8.1)") "   convergence      : dE ", r%econv, &
            "   |FDS-SDF| ", r%gconv, "   screening ", r%thresh
         call global_logger%info(trim(logline))
      end block
      if (len(r%func) > 0) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0)") "   grid level       : ", r%grid
            if (r%vv10) write (logline, "(a,i0,a,i0)") "   grid level       : ", r%grid, &
               "   VV10: ", r%nlc_grid
            call global_logger%info(trim(logline))
         end block
      end if
      if (len(r%pcm) > 0) then
         block
            character(len=1024) :: logline
            write (logline, "(a,a,a,f0.4)") "   solvation        : ", r%pcm, ", eps ", r%eps
            call global_logger%info(trim(logline))
         end block
      end if
   end subroutine banner

   function formula(z, natm) result(f)
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
      character(len=:), allocatable :: f
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
      f = ""
      do ia = 1, n
         f = f//trim(element_symbol(order(ia)))
         if (cnt(order(ia)) > 1) f = f//to_string(cnt(order(ia)))
      end do
   end function formula

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
      call global_logger%info("")
      block
         character(len=1024) :: logline
         write (logline, "(a,i0,a,a)") "   AO basis         : ", r%nao, " functions, ", &
            merge("spherical", "Cartesian", r%sph /= 0)
         call global_logger%info(trim(logline))
      end block
      if (len(r%restart) > 0) then
         call global_logger%info("   SCF              : "//to_string(r%niter) &
                                 //trim(merge(" iterations, converged    ", " iterations, NOT converged", r%conv /= 0)) &
                                 //" (read from "//r%restart//")")
      else
         call global_logger%info("   SCF              : "//to_string(r%niter) &
                                 //trim(merge(" iterations, converged    ", " iterations, NOT converged", r%conv /= 0)))
      end if
      call global_logger%info("")
      call global_logger%info(bar)
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
      call global_logger%info(bar)
      call global_logger%info("")
      call global_logger%info("   ::::::::::::::::::::::::::::::::::::::::::::::::::::::")
      call global_logger%info("   ::"//repeat(" ", 21)//"SUMMARY"//repeat(" ", 22)//"::")
      call global_logger%info("   ::::::::::::::::::::::::::::::::::::::::::::::::::::::")
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
      call global_logger%info("   ::::::::::::::::::::::::::::::::::::::::::::::::::::::")
      block
         character(len=1024) :: logline
         write (logline, "(a,f12.3,a)") "   wall time        : ", real(r%t1 - r%t0, dp)/real(r%rate, dp), " s"
         call global_logger%info(trim(logline))
      end block
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
      call global_logger%info("")
      call global_logger%info("   gradient (Eh/Bohr)")
      call global_logger%info("   ----------------------------------------------------------")
      do ia = 1, r%natm
         block
            character(len=1024) :: logline
            write (logline, "(a,i4,2x,a2,3f16.8)") "   ", ia, element_symbol(r%z(ia)), r%gradient(:, ia)
            call global_logger%info(trim(logline))
         end block
      end do
      rms = sqrt(sum(r%gradient**2)/real(3*r%natm, dp))
      gmax = sqrt(maxval(sum(r%gradient**2, dim=1)))
      call global_logger%info("   ----------------------------------------------------------")
      block
         character(len=1024) :: logline
         write (logline, "(a,es12.4)") "   RMS gradient     : ", rms
         call global_logger%info(trim(logline))
      end block
      block
         character(len=1024) :: logline
         write (logline, "(a,es12.4)") "   max |gradient|   : ", gmax
         call global_logger%info(trim(logline))
      end block
   end subroutine gradient_summary

   !
   ! Per-atom masses for --hessian's vibrational analysis: "isotope"
   ! (trc_atomic_mass, the most-abundant isotope -- the CLI's own
   ! default) or "average" (trc_atomic_mass_average, the isotope-averaged
   ! standard atomic weight). `kind` is already lower-cased by the
   ! caller; anything else is a programming error (the CLI validates
   ! --masses before this is ever called).
   !
   subroutine build_masses(z, kind, mass, err)
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
   end subroutine build_masses

   !
   ! Whether this run's method has a gradient at all, checked BEFORE any
   ! SCF runs: every one of these conditions is already known from the
   ! command line alone, so there is no need to spend an SCF on a probe
   ! gradient just to find out. The same set trc_run_gradient's own
   ! gradient_capability_check refuses after a converged SCF
   ! (src/api/trc_capi.F90). Returns "" when a gradient is available, the
   ! refusal message otherwise.
   !
   function opt_capability_check(r) result(refusal)
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
   end function opt_capability_check

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
   function opt_geometry_is_linear(natm, z, xyz, err) result(linear)
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: xyz(3, natm)
      type(mass_error_t), intent(inout) :: err
      logical :: linear
      real(dp), allocatable :: mass(:)
      integer :: n_tr, info
      linear = .false.
      call build_masses(z, "isotope", mass, err)
      if (err%has_error()) return
      call trc_vib_rotor(mass, xyz, linear, n_tr, info)
   end function opt_geometry_is_linear

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
   subroutine build_opt_residues(natm, z, xyz, residues, err)
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
   end subroutine build_opt_residues

   !
   ! The optimisation's own final summary: converged or not, steps,
   ! energy evaluations and the final energy.
   ! Per-step lines (energy, max/rms |g|, SCF iterations, wall time) are
   ! logged from trc_dlfind_bridge%cb_put_coords itself as each step is
   ! accepted, not replayed here.
   !
   subroutine opt_summary(r)
      type(run_t), intent(in) :: r
      integer :: n_imag
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
      ! AFTER --opt, at the geometry --opt reports -- terco_main's own
      ! ordering), since that is what computes r%freq_cm1 at all.
      if (r%opt_target_kind == 1 .and. r%hess_ok) then
         n_imag = count(r%freq_cm1 < 0.0_dp)
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a)") "   saddle check: ", n_imag, " imaginary mode(s) (expect exactly 1 "// &
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
   end subroutine opt_summary

   pure function opt_coords_name(kind) result(name)
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

   pure function opt_target_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name
      if (kind == 1) then
         name = "saddle"
      else
         name = "minimum"
      end if
   end function opt_target_name

   pure function opt_saddle_method_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=:), allocatable :: name
      if (kind == 1) then
         name = "dimer"
      else
         name = "prfo"
      end if
   end function opt_saddle_method_name

   pure function opt_algorithm_name(kind) result(name)
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
   subroutine write_opt_xyz(r, path, xyz_bohr)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: xyz_bohr(3, r%natm)
      call write_opt_xyz_frame(path, r%natm, r%z, real(r%etot, dp), xyz_bohr, &
                               trim(merge("converged    ", "NOT converged", r%opt_converged)))
   end subroutine write_opt_xyz

   !
   ! `<geometry base>_ts.xyz`/`_minus.xyz`/`_plus.xyz`: --opt-connect's own
   ! three structures (the transition state and its two downhill minima),
   ! one labelled energy each -- write_opt_xyz's own format, but with an
   ! energy that is NOT r%etot (the context is re-converged at the TS,
   ! final_coords, after --opt returns; the two minima's own energies came
   ! straight from DL-FIND's cb_put_coords and were never re-run).
   !
   subroutine write_labelled_xyz(r, path, xyz_bohr, e, label)
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
         write (u, "(a2,3f20.12)") element_symbol(r%z(ia)), xyz_bohr(:, ia)/BOHR_PER_ANGSTROM
      end do
      close (u)
   end subroutine write_labelled_xyz

   !
   ! `<geometry base>_traj.xyz`: one frame per ACCEPTED step (never a
   ! line-search trial -- trc_dlfind_bridge%cb_put_coords's own
   ! distinction), Angstrom, comment line carries the step index and
   ! energy.
   !
   subroutine write_traj_xyz(r, path, traj_xyz_bohr, traj_energy)
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: traj_xyz_bohr(3, r%natm, r%opt_n_traj)
      real(dp), intent(in) :: traj_energy(r%opt_n_traj)
      call write_traj_xyz_frames(path, r%natm, r%z, r%opt_n_traj, traj_xyz_bohr, traj_energy)
   end subroutine write_traj_xyz

   !
   ! The frequency table and thermochemistry summary for --hessian,
   ! through trc_vib's/trc_thermo's own logger helpers -- `r%freq_cm1`
   ! etc. and `r%thermo` are already filled by terco_main itself, right
   ! after trc_vibrational_analysis/trc_thermochemistry succeed.
   !
   subroutine hessian_summary(r)
      type(run_t), intent(in) :: r
      call global_logger%info("")
      call global_logger%info("   vibrational analysis and thermochemistry (--hessian)")
      call global_logger%info("   ----------------------------------------------------------")
      block
         character(len=1024) :: logline
         write (logline, "(a,f0.4,a,a,a,i0,a,i0)") "   step ", r%hess_step, " Bohr   masses: ", r%masses_kind, &
            "   symmetry number: ", r%sym_number, " (not determined automatically)   mult: ", r%mult
         call global_logger%info(trim(logline))
      end block
      call trc_vib_log_table(r%freq_cm1, r%redmass_amu, r%force_const)
      call trc_thermo_log_summary(r%thermo)
   end subroutine hessian_summary

   subroutine line(label, e)
      character(len=*), intent(in) :: label
      real(c_double), intent(in) :: e
      block
         character(len=1024) :: logline
         write (logline, "(a,a20,f22.12,a)") "          | ", label, e, " Eh    |"
         call global_logger%info(trim(logline))
      end block
   end subroutine line

   subroutine row(label, e)
      character(len=*), intent(in) :: label
      real(c_double), intent(in) :: e
      block
         character(len=1024) :: logline
         write (logline, "(a,a22,f22.12,a)") "   :: ", label, e, " Eh  ::"
         call global_logger%info(trim(logline))
      end block
   end subroutine row

   subroutine write_json(r, path)
      !! The run's JSON summary, built with json-fortran, which escapes the
      !! strings and writes the reals at full precision. The key names and
      !! their order are an interface: scripts and the check_cli_grad tests
      !! read them.
      type(run_t), intent(in) :: r
      character(len=*), intent(in) :: path
      type(json_core) :: json
      type(json_value), pointer :: root
      call json%initialize()
      call json%create_object(root, "")
      call json%add(root, "program", "terco")
      call json%add(root, "version", VERSION)
      call json%add(root, "geometry", r%geom)
      call json%add(root, "method", method(r))
      call json%add(root, "basis", r%basis_path)
      call json%add(root, "guess", r%guess)
      call json%add(root, "charge", r%chrg)
      call json%add(root, "unpaired electrons", r%uhf)
      call json%add(root, "number of basis functions", int(r%nao))
      call json%add(root, "spherical", r%sph /= 0)
      call json%add(root, "iterations", int(r%niter))
      call json%add(root, "converged", r%conv /= 0)
      call json%add(root, "nuclear repulsion energy", real(r%enuc, dp))
      call json%add(root, "one-electron energy", real(r%eone, dp))
      call json%add(root, "two-electron energy", real(r%etwo, dp))
      call json%add(root, "exchange-correlation energy", real(r%exc, dp))
      if (len(r%pcm) > 0) then
         call json%add(root, "solvation model", r%pcm)
         call json%add(root, "dielectric constant", r%eps)
         call json%add(root, "solvation energy", real(r%epcm, dp))
         call json%add(root, "solvation surface points", int(r%npcm))
      end if
      if (r%rimp2 .and. r%conv /= 0) then
         call json%add(root, "mp2 opposite-spin energy", real(r%eos, dp))
         call json%add(root, "mp2 same-spin energy", real(r%ess, dp))
      end if
      if (r%rimp2 .and. r%conv /= 0 .and. r%dh) then
         call json%add(root, "pt2 opposite-spin fraction", real(r%pt2_cos, dp))
         call json%add(root, "pt2 same-spin fraction", real(r%pt2_css, dp))
         call json%add(root, "scf energy", real(r%etot, dp))
         call json%add(root, "total energy", real(r%etot_dh, dp))
      else
         call json%add(root, "scf energy", real(r%etot, dp))
         call json%add(root, "total energy", real(r%etot, dp))
         if (r%grad .and. allocated(r%gradient)) call add_gradient(json, root, r)
      end if
      if (r%hessian .and. r%hess_ok) call add_hessian(json, root, r)
      if (r%opt) call add_optimization(json, root, r)
      call json%print(root, path)
      if (json%failed()) then
         call json%destroy(root)
         call die(r, "cannot write "//path)
      end if
      call json%destroy(root)
   end subroutine write_json

   subroutine add_gradient(json, root, r)
      !! The analytic gradient as an array of [x, y, z] triples, one per atom
      !! in `r%gradient`'s own (and `trc_set_molecule`'s) atom order, plus the
      !! RMS and max per-atom norm `gradient_summary` also prints.
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: rows
      integer :: ia
      call json%add(root, "gradient status", "validated: single rank (rank 0 alone on a multi-rank run)")
      call json%create_array(rows, "gradient")
      do ia = 1, r%natm
         call json%add(rows, "", r%gradient(:, ia))
      end do
      call json%add(root, rows)
      call json%add(root, "gradient rms", sqrt(sum(r%gradient**2)/real(3*r%natm, dp)))
      call json%add(root, "gradient max", sqrt(maxval(sum(r%gradient**2, dim=1))))
   end subroutine add_gradient

   !
   ! --hessian's JSON: three blocks, as asked for -- "hessian" (the
   ! matrix itself and the step it was built at), "frequencies"
   ! (trc_vibrational_analysis's own per-mode output) and
   ! "thermochemistry" (trc_thermo_result_t, trc_thermo.F90's own
   ! fields, nothing recomputed here). Row/column k = 3*(atom-1)+
   ! {x,y,z}, the same order `trc_hessian` and `trc_set_molecule` use.
   !
   subroutine add_hessian(json, root, r)
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: hobj, hrows, fobj, tobj, tsub
      integer :: ia, ndof

      ndof = 3*r%natm

      call json%create_object(hobj, "hessian")
      call json%add(hobj, "step_bohr", r%hess_step)
      call json%create_array(hrows, "matrix")
      do ia = 1, ndof
         call json%add(hrows, "", r%hess(ia, :))
      end do
      call json%add(hobj, hrows)
      call json%add(root, hobj)

      call json%create_object(fobj, "frequencies")
      call json%add(fobj, "n_modes", size(r%freq_cm1))
      call json%add(fobj, "n_imaginary", r%thermo%n_imaginary)
      call json%add(fobj, "is_linear", r%is_linear)
      call json%add(fobj, "frequencies_cm1", r%freq_cm1)
      call json%add(fobj, "reduced_masses_amu", r%redmass_amu)
      call json%add(fobj, "force_constants_hartree_bohr2", r%force_const)
      call json%add(root, fobj)

      call json%create_object(tobj, "thermochemistry")
      call json%add(tobj, "temperature_K", r%thermo%temperature_k)
      call json%add(tobj, "pressure_atm", r%pressure_atm)
      call json%add(tobj, "symmetry_number", r%sym_number)
      call json%add(tobj, "spin_multiplicity", r%mult)
      call json%add(tobj, "masses", r%masses_kind)
      call json%add(tobj, "zero_point_energy_hartree", r%thermo%zpe_hartree)

      call json%create_object(tsub, "thermal_corrections_hartree")
      call json%add(tsub, "to_energy", r%thermo%delta_e_hartree)
      call json%add(tsub, "to_enthalpy", r%thermo%delta_h_hartree)
      call json%add(tsub, "to_gibbs", r%thermo%delta_g_hartree)
      call json%add(tobj, tsub)

      call json%create_object(tsub, "total_energies_hartree")
      call json%add(tsub, "electronic", r%thermo%e_elec_hartree)
      call json%add(tsub, "electronic_plus_zpe", r%thermo%e_zpe_hartree)
      call json%add(tsub, "electronic_plus_thermal_E", r%thermo%e_tot_hartree)
      call json%add(tsub, "electronic_plus_thermal_H", r%thermo%h_tot_hartree)
      call json%add(tsub, "electronic_plus_thermal_G", r%thermo%g_tot_hartree)
      call json%add(tobj, tsub)

      call json%create_object(tsub, "entropy_cal_mol_K")
      call json%add(tsub, "translational", r%thermo%s_trans_calmolk)
      call json%add(tsub, "rotational", r%thermo%s_rot_calmolk)
      call json%add(tsub, "vibrational", r%thermo%s_vib_calmolk)
      call json%add(tsub, "electronic", r%thermo%s_elec_calmolk)
      call json%add(tsub, "total", r%thermo%s_tot_calmolk)
      call json%add(tobj, tsub)

      call json%create_object(tsub, "heat_capacity_cal_mol_K")
      call json%add(tsub, "translational", r%thermo%cv_trans_calmolk)
      call json%add(tsub, "rotational", r%thermo%cv_rot_calmolk)
      call json%add(tsub, "vibrational", r%thermo%cv_vib_calmolk)
      call json%add(tsub, "electronic", r%thermo%cv_elec_calmolk)
      call json%add(tsub, "total", r%thermo%cv_tot_calmolk)
      call json%add(tobj, tsub)

      call json%add(root, tobj)
   end subroutine add_hessian

   !
   ! --opt's JSON block: converged, step/evaluation counts, the final
   ! energy and gradient, the coordinate system, the trajectory's
   ! per-step energies and the final geometry. The final geometry is in
   ! Angstrom, named accordingly so a reader never has to guess the unit
   ! the way it would from a bare "final_geometry".
   !
   subroutine add_optimization(json, root, r)
      type(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: root
      type(run_t), intent(in) :: r
      type(json_value), pointer :: oobj, rows, csub
      integer :: ia

      call json%create_object(oobj, "optimization")
      call json%add(oobj, "converged", r%opt_converged)
      call json%add(oobj, "steps", r%opt_steps)
      call json%add(oobj, "energy_evaluations", r%opt_evals)
      call json%add(oobj, "final_energy", real(r%etot, dp))
      call json%add(oobj, "final_gradient_max", r%opt_gmax)
      call json%add(oobj, "final_gradient_rms", r%opt_grms)
      call json%add(oobj, "gradient_tolerance", r%opt_gtol)
      call json%add(oobj, "coordinates", opt_coords_name(r%opt_coords_kind))
      call json%add(oobj, "guess", trim(r%opt_guess_arg))
      call json%add(oobj, "target", opt_target_name(r%opt_target_kind))
      if (r%opt_target_kind == 1) then
         call json%add(oobj, "saddle_method", opt_saddle_method_name(r%opt_saddle_method_kind))
      else
         call json%add(oobj, "algorithm", opt_algorithm_name(r%opt_algorithm_kind))
      end if
      call json%add(oobj, "connect", r%opt_connect)
      call json%add(oobj, "hessian_evaluations", r%opt_hess_evals)
      if (allocated(r%opt_traj_energy)) then
         call json%add(oobj, "trajectory_energies", r%opt_traj_energy)
      end if
      if (allocated(r%opt_final_xyz)) then
         call json%create_array(rows, "final_geometry_angstrom")
         do ia = 1, r%natm
            call json%add(rows, "", r%opt_final_xyz(:, ia)/BOHR_PER_ANGSTROM)
         end do
         call json%add(oobj, rows)
      end if
      if (r%opt_connect) then
         ! Each sub-search's OWN convergence verdict and gradient, never
         ! another sub-search's -- the top-level converged/final_gradient_max
         ! above are the TS search's own (trc_dlfind_bridge's own override),
         ! but a reader who wants the two downhill minima's separately
         ! finds them here instead of having to infer them.
         call json%add(oobj, "ts_found", r%opt_have_ts)
         if (r%opt_have_ts) then
            call json%create_object(csub, "ts")
            call json%add(csub, "energy", r%opt_ts_energy)
            call json%add(csub, "converged", r%opt_ts_converged)
            call json%add(csub, "final_gradient_max", r%opt_ts_gmax)
            call json%add(oobj, csub)
         end if
         call json%add(oobj, "minimum_minus_found", r%opt_have_min_minus)
         if (r%opt_have_min_minus) then
            call json%create_object(csub, "minimum_minus")
            call json%add(csub, "energy", r%opt_min_minus_energy)
            call json%add(csub, "converged", r%opt_min_minus_converged)
            call json%add(csub, "final_gradient_max", r%opt_min_minus_gmax)
            call json%add(oobj, csub)
         end if
         call json%add(oobj, "minimum_plus_found", r%opt_have_min_plus)
         if (r%opt_have_min_plus) then
            call json%create_object(csub, "minimum_plus")
            call json%add(csub, "energy", r%opt_min_plus_energy)
            call json%add(csub, "converged", r%opt_min_plus_converged)
            call json%add(csub, "final_gradient_max", r%opt_min_plus_gmax)
            call json%add(oobj, csub)
         end if
      end if
      call json%add(root, oobj)
   end subroutine add_optimization

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
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init
   use pic_knowledge, only: get_knowledge
   use pic_profiler, only: profiler_init
   use pic_logger, only: global_logger
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                       trc_set_aux_json, trc_set_method, trc_set_nlc_grid, trc_set_convergence, trc_set_screening, &
                       trc_set_incremental, trc_set_df, trc_set_jkfit_json, &
                       trc_set_guess, trc_set_verbose, trc_set_rimp2, trc_set_force_cartesian, trc_set_comm, &
                       trc_run_scf, trc_run_rimp2, trc_nao, trc_is_spherical, trc_nspin, trc_energy, &
                       trc_energy_parts, trc_converged, trc_iterations, trc_rimp2_energy, &
                       trc_set_pcm, trc_pcm_energy, trc_save_scf, trc_load_scf, &
                       trc_run_gradient, trc_gradient, trc_set_geometry, trc_set_hessian_step, trc_run_hessian, &
                       trc_hessian, &
                       TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_CORE, TRC_GUESS_GWH, TRC_GUESS_SAD, &
                       TRC_GUESS_SADQ, TRC_GUESS_SAP
   use trc_api, only: trc_bind_device, trc_bind_device_env
   use trc_xc_functional, only: XC_AVAILABLE_FUNCTIONALS, trc_xc_functional_t, xc_functional_by_name, &
                                xc_is_double_hybrid
   use trc_vib, only: trc_vibrational_analysis
   use trc_thermo, only: trc_thermochemistry
   use trc_constants, only: ATM_TO_PA
   use trc_error, only: hess_err_t => error_t
   use pic_logger, only: global_logger
   use terco_app, only: read_xyz, resolve_basis, c_string, set_env
   use pic_strings, only: to_string
   use terco_run, only: run_t, VERSION, die, quit, check, api, int_option, banner, summary, write_json, &
                       gradient_summary, hessian_summary, build_masses, &
                       opt_capability_check, build_opt_residues, opt_summary, write_opt_xyz, write_traj_xyz, &
                       opt_geometry_is_linear, write_labelled_xyz, to_lower_alloc
   use trc_dlfind_bridge, only: trc_dlfind_available, trc_opt_options_t, trc_opt_result_t, trc_dlfind_optimize, &
                                TRC_OPT_COORDS_CARTESIAN, TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_HDLC_TC, &
                                TRC_OPT_COORDS_DLC, TRC_OPT_COORDS_DLC_TC, &
                                TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, &
                                TRC_OPT_ALG_SD, TRC_OPT_ALG_CG_AUTO, TRC_OPT_ALG_CG, TRC_OPT_ALG_LBFGS, &
                                TRC_OPT_ALG_PRFO, TRC_OPT_ALG_NR, TRC_OPT_ALG_DAMPED, &
                                TRC_OPT_TARGET_MINIMUM, TRC_OPT_TARGET_SADDLE, &
                                TRC_OPT_SADDLE_PRFO, TRC_OPT_SADDLE_DIMER, &
                                TRC_OPT_HESS_TERCO, TRC_OPT_HESS_DLFIND, &
                                TRC_OPT_HUPD_NONE, TRC_OPT_HUPD_POWELL, TRC_OPT_HUPD_BOFILL, TRC_OPT_HUPD_BFGS, &
                                TRC_OPT_LINE_AUTO, TRC_OPT_LINE_SCALING, TRC_OPT_LINE_ENERGY_TRUST, &
                                TRC_OPT_LINE_GRADIENT_TRUST
   implicit none

   type(run_t) :: r
   character(len=:), allocatable :: sval
   real(dp), allocatable :: xyz(:, :)
   real(c_double), allocatable :: gradbuf(:), hessbuf(:)
   real(dp), allocatable :: mass(:), displacement(:, :, :)
   type(hess_err_t) :: herr
   integer :: maxiter, nfrozen, guess_kind, dev, aux_block, ndof
   integer(c_int) :: st, nspin
   logical :: cart
   ! --opt
   character(len=:), allocatable :: opt_base, opt_refusal
   integer :: opt_base_slash
   logical :: opt_hessian_relevant, opt_zero_soft_relevant
   integer, allocatable :: opt_residues(:)
   type(trc_opt_options_t) :: opt_opts
   type(trc_opt_result_t) :: opt_res

   ! MPI, in every build: pic-mpi's single-rank backend (TERCO_ENABLE_MPI=OFF)
   ! makes this a no-op, so the rest of the program need not branch on it.
   ! `r%wcomm` and `r%nranks` answer everything the rest of this file needs
   ! to know about ranks; `r%rank0` is what gates every print below so N
   ! ranks produce ONE copy of the output, not N.
   !
   ! The device first, from the launcher's node-local rank: a GPU-aware MPI
   ! records the current device in MPI_Init, and switching afterwards costs
   ! Cray MPICH its CUDA IPC (trc_bind_device_env). The split-communicator
   ! bind below still runs, and is a no-op when it lands on this device.
   dev = trc_bind_device_env()
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
   call r%cli%add_flag("grad", "Print the analytic nuclear gradient (Eh/Bohr) after the SCF -- "// &
                       "RHF/UHF, or KS-DFT for an LDA/GGA --func (global and range-separated "// &
                       "hybrids and VV10 included); with --df, the fitted RHF/RKS gradient; "// &
                       "refused for a meta-GGA, --rimp2, a double hybrid or --pcm; with more than "// &
                       "one rank it is computed on rank 0 alone, with a warning")
   call r%cli%add_flag("hessian", "finite-difference Cartesian Hessian (central differences of the "// &
                       "analytic gradient) after the SCF, then harmonic vibrational analysis and RRHO "// &
                       "thermochemistry; refused wherever --grad is (meta-GGA, --rimp2, a double "// &
                       "hybrid, --pcm), checked before any displacement runs; single rank only for now")
   call r%cli%add_option("step", "--hessian's displacement, Bohr (0.005)", default="0.005")
   call r%cli%add_option("temperature", "temperature for --hessian's thermochemistry, K (298.15)", &
                         default="298.15")
   call r%cli%add_option("pressure", "pressure for --hessian's thermochemistry, atm (1.0)", default="1.0")
   call r%cli%add_option("symmetry-number", "rotational symmetry number for --hessian's thermochemistry "// &
                         "(1; not determined automatically -- a warning names the quantities it affects)", &
                         default="1")
   call r%cli%add_option("masses", "isotope (trc_atomic_mass, the most-abundant isotope) or average "// &
                         "(trc_atomic_mass_average, the isotope-averaged standard atomic weight), for "// &
                         "--hessian's vibrational analysis (isotope)", default="isotope")
   call r%cli%add_flag("opt", "optimise the geometry with DL-FIND (needs a TERCO_ENABLE_DLFIND=ON build; "// &
                       "refused otherwise). Warm-starts the SCF between steps; writes <geometry>_opt.xyz "// &
                       "and <geometry>_traj.xyz into the CURRENT WORKING DIRECTORY (named from the input "// &
                       "geometry's own base name), and an 'optimization' block to --json. "// &
                       "Combined with --hessian, the Hessian runs at the OPTIMISED geometry. If a step's "// &
                       "SCF fails to converge even after one retry, DL-FIND has no recovery: terco writes "// &
                       "the files/JSON above (marked failed, with the last ACCEPTED geometry) and the "// &
                       "process then stops abnormally from inside DL-FIND, not through terco's own exit path")
   call r%cli%add_option("opt-coords", "internal coordinates for --opt: cartesian, hdlc, hdlc-tc, dlc "// &
                         "or dlc-tc (cartesian)", default="cartesian")
   call r%cli%add_option("opt-maxiter", "maximum number of --opt steps (100)", default="100")
   call r%cli%add_option("opt-gtol", "--opt convergence: this run's own verdict is the max Cartesian "// &
                         "gradient component of the last accepted step, Eh/Bohr (4.5e-4) -- NOT "// &
                         "DL-FIND's own five-criterion internal test, which in hdlc/dlc is over internal "// &
                         "coordinates and is not this number", default="4.5e-4")
   call r%cli%add_option("opt-guess", "SCF warm start between --opt steps: density (the previous step's "// &
                         "converged density) or projected (its orbitals re-orthonormalised in the new "// &
                         "overlap) (density)", default="density")
   call r%cli%add_option("opt-algorithm", "--opt's minimiser: sd, cg, cg-auto, lbfgs, nr, damped or prfo "// &
                         "(lbfgs); refused together with --opt-target saddle, which picks its own "// &
                         "algorithm from --opt-saddle-method", default="lbfgs")
   call r%cli%add_option("opt-target", "minimum or saddle (minimum)", default="minimum")
   call r%cli%add_option("opt-saddle-method", "prfo (uses a Hessian, --opt-init-hessian) or dimer (no "// &
                         "Hessian); needs --opt-target saddle (prfo)", default="prfo")
   call r%cli%add_flag("opt-connect", "after a saddle search (--opt-target saddle), displace along the "// &
                       "imaginary mode both ways and minimise each to a downhill minimum (DL-FIND task "// &
                       "1011); needs --opt-target saddle. Writes <geometry>_ts.xyz, <geometry>_minus.xyz "// &
                       "and <geometry>_plus.xyz in addition to the usual --opt output")
   call r%cli%add_option("opt-init-hessian", "for --opt-target saddle --opt-saddle-method prfo (or "// &
                         "--opt-algorithm nr): terco (terco's own FD Hessian, trc_run_hessian, via "// &
                         "cb_get_hessian) or dlfind (DL-FIND's own two-point FD of the gradient) (terco)", &
                         default="terco")
   call r%cli%add_option("opt-hessian-update", "Hessian update for prfo/nr: none, powell, bofill or bfgs "// &
                         "(DL-FIND's own default, bofill)")
   call r%cli%add_option("opt-max-updates", "maximum Hessian updates before prfo/nr recomputes one "// &
                         "(DL-FIND's own default, 50)")
   call r%cli%add_option("opt-zero-modes", "modes assumed to have zero vibrational frequency, for "// &
                         "prfo/nr's Hessian (nzero); default 6 with --opt-coords cartesian and "// &
                         "--opt-algorithm prfo/nr (minimum or saddle) or --opt-saddle-method dimer (5 "// &
                         "if the start geometry is linear -- otherwise the search climbs a rotation), "// &
                         "DL-FIND's own default (0) otherwise")
   call r%cli%add_option("opt-soft-modes", "eigenmodes of prfo/nr's Hessian below this value (Eh/Bohr^2) "// &
                         "are treated as soft (DL-FIND's own default: none)")
   call r%cli%add_option("opt-max-step", "DL-FIND's own ceiling on a trust-radius step, Bohr -- also the "// &
                         "ceiling the trust radius may grow to across accepted steps (DL-FIND's own "// &
                         "default, 0.5, unless given); raising it lets the optimiser attempt a larger, "// &
                         "possibly-rejected trial sooner")
   call r%cli%add_option("opt-linesearch", "trust radius/line search: auto (DL-FIND's own default, which "// &
                         "depends on --opt-algorithm), scaling, energy-trust or gradient-trust (auto)", &
                         default="auto")
   call r%cli%add_option("opt-scalestep", "step scale factor for --opt-linesearch scaling (DL-FIND's "// &
                         "own default, 1.0)")
   call r%cli%add_option("opt-lbfgs-memory", "L-BFGS history length for --opt-algorithm lbfgs "// &
                         "(DL-FIND's own default, min(3*natm, 50))")
   call r%cli%add_option("opt-etol", "--opt's energy-change convergence criterion, Eh (DL-FIND's own "// &
                         "default, --opt-gtol/450)")
   call r%cli%add_option("opt-print-level", "DL-FIND's own stdout verbosity (printl): 0, 2, 4 or 6 "// &
                         "(follows this process's own logger level: 0 normally, 2 under --verbose, 4 "// &
                         "under --debug)")
   call r%cli%add_flag("silent", "only the summary: no SCF iterations, timings or guess details", short="s")
   call r%cli%add_flag("version", "print the version and stop")

   call r%cli%parse(r%err)
   if (r%err%has_error()) call die(r, r%err%get_message()//new_line("a")//"try 'terco --help'")
   if (r%cli%help_requested()) then
      if (r%rank0) call global_logger%info(r%cli%help_text())
      call quit(r, 0)
   end if
   if (r%cli%is_set("version")) then
      if (r%rank0) call global_logger%info("terco version "//VERSION)
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
   if (r%nelec < 1) call die(r, "charge "//to_string(r%chrg)//" leaves "//to_string(r%nelec)//" electrons")
   r%uhf = int_option(r, "uhf", ".UHF", mod(r%nelec, 2))
   if (r%uhf < 0 .or. r%uhf > r%nelec .or. mod(r%nelec - r%uhf, 2) /= 0) then
     call die(r, "--uhf "//to_string(r%uhf)//" is impossible with "//to_string(r%nelec)//" electrons: the unpaired "// &
               "count must be between 0 and the electron count, with the same parity")
   end if
   r%mult = r%uhf + 1

   call r%cli%get("func", sval, r%err)
   r%func = to_lower_alloc(trim(sval))
   if (r%func == "hf" .or. r%func == "rhf" .or. r%func == "uhf") r%func = ""

   ! A double hybrid's PT2 part is always RI-MP2 on the converged KS
   ! orbitals (no canonical, non-fitted MP2 path exists in terco at all):
   ! --func b2plyp runs it without needing --rimp2 too, with --aux's usual
   ! default (cc-pvdz-rifit) unless overridden below.
   if (len(r%func) > 0) then
      block
         type(trc_xc_functional_t) :: func_spec
         if (xc_functional_by_name(r%func, func_spec)) then
            r%dh = xc_is_double_hybrid(func_spec)
            r%vv10 = func_spec%nlc_b /= 0.0_dp
         end if
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
      r%guess = to_lower_alloc(trim(sval))
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
   call r%cli%get("nlc-grid", r%nlc_grid, r%err)
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
   r%hessian = r%cli%is_set("hessian")
   cart = r%cli%is_set("cart")
   if (r%hessian) then
      call r%cli%get("step", r%hess_step, r%err)
      call check(r, "--step")
      if (r%hess_step <= 0.0_dp) call die(r, "--step must be positive, got "//to_string(r%hess_step))
      call r%cli%get("temperature", r%temperature_k, r%err)
      call check(r, "--temperature")
      if (r%temperature_k <= 0.0_dp) call die(r, "--temperature must be positive")
      call r%cli%get("pressure", r%pressure_atm, r%err)
      call check(r, "--pressure")
      if (r%pressure_atm <= 0.0_dp) call die(r, "--pressure must be positive")
      call r%cli%get("symmetry-number", r%sym_number, r%err)
      call check(r, "--symmetry-number")
      if (r%sym_number < 1) call die(r, "--symmetry-number must be at least 1")
      if (r%rank0) call global_logger%warning("--hessian: the rotational symmetry number ("// &
         to_string(r%sym_number)//") is not determined automatically -- it affects the rotational "// &
         "entropy and so the Gibbs free energy; pass --symmetry-number if the molecule's point group "// &
         "has one above 1 (water: 2, ammonia: 3, ...)")
      call r%cli%get("masses", sval, r%err)
      call check(r, "--masses")
      r%masses_kind = to_lower_alloc(trim(sval))
      if (r%masses_kind /= "isotope" .and. r%masses_kind /= "average") then
         call die(r, "--masses must be isotope or average, got '"//r%masses_kind//"'")
      end if
   end if
   r%opt = r%cli%is_set("opt")
   r%pcm = ""
   r%eps = 0.0_dp
   if (r%cli%is_set("pcm")) then
      call r%cli%get("pcm", sval, r%err)
      r%pcm = to_lower_alloc(trim(sval))
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
   ! Every --opt-* option is parsed inside `if (r%opt) then` below, so one
   ! given without --opt itself would otherwise do nothing silently --
   ! refused here instead, before that block, so a command line that sets
   ! one always means what it says.
   if (.not. r%opt) then
      block
         character(len=*), parameter :: opt_options(*) = [character(len=18) :: &
            "opt-coords", "opt-guess", "opt-maxiter", "opt-gtol", "opt-target", "opt-saddle-method", &
            "opt-connect", "opt-algorithm", "opt-init-hessian", "opt-hessian-update", "opt-max-updates", &
            "opt-zero-modes", "opt-soft-modes", "opt-max-step", "opt-linesearch", "opt-scalestep", &
            "opt-lbfgs-memory", "opt-etol", "opt-print-level"]
         integer :: k
         do k = 1, size(opt_options)
            if (r%cli%is_set(trim(opt_options(k)))) call die(r, "--"//trim(opt_options(k))//" needs --opt")
         end do
      end block
   end if
   ! Validated LAST among the feature flags above: opt_capability_check
   ! reads r%func/r%rimp2/r%dh/r%pcm, every one of which must already be
   ! resolved (not just r%cli%is_set) for the check to see what the run
   ! will actually do -- a PCM run that reached DL-FIND with this check
   ! run too early crashed with `error stop` instead of this message
   ! (this is exactly why opt_capability_check exists: DL-FIND itself
   ! has no graceful way to decline a run it cannot do).
   if (r%opt) then
      if (.not. trc_dlfind_available()) call die(r, "--opt needs libdlfind: this terco was built with "// &
                                                  "-DTERCO_ENABLE_DLFIND=OFF (the default). Reconfigure "// &
                                                  "with -DTERCO_ENABLE_DLFIND=ON.")
      if (r%nranks > 1) call die(r, "--opt is single-rank only for now (DL-FIND itself is serial; it is "// &
                                   "driven from the CLI, never run collectively)")
      if (len(r%restart) > 0) call die(r, "--opt cannot start from --restart: restart loads one SCF "// &
                                       "result, and --opt needs to move the geometry and build a new "// &
                                       "one at every step")
      opt_refusal = opt_capability_check(r)
      if (len(opt_refusal) > 0) call die(r, opt_refusal)
      call r%cli%get("opt-coords", sval, r%err)
      r%opt_coords_arg = to_lower_alloc(trim(sval))
      select case (r%opt_coords_arg)
      case ("cartesian")
         r%opt_coords_kind = TRC_OPT_COORDS_CARTESIAN
      case ("hdlc")
         r%opt_coords_kind = TRC_OPT_COORDS_HDLC
      case ("hdlc-tc")
         r%opt_coords_kind = TRC_OPT_COORDS_HDLC_TC
      case ("dlc")
         r%opt_coords_kind = TRC_OPT_COORDS_DLC
      case ("dlc-tc")
         r%opt_coords_kind = TRC_OPT_COORDS_DLC_TC
      case default
         call die(r, "--opt-coords must be cartesian, hdlc, hdlc-tc, dlc or dlc-tc, got '"//r%opt_coords_arg//"'")
      end select
      call r%cli%get("opt-guess", sval, r%err)
      r%opt_guess_arg = to_lower_alloc(trim(sval))
      select case (r%opt_guess_arg)
      case ("density")
         r%opt_guess_kind = TRC_OPT_GUESS_DENSITY
      case ("projected")
         r%opt_guess_kind = TRC_OPT_GUESS_PROJECTED
      case default
         call die(r, "--opt-guess must be density or projected, got '"//r%opt_guess_arg//"'")
      end select
      call r%cli%get("opt-maxiter", r%opt_maxiter, r%err)
      call check(r, "--opt-maxiter")
      if (r%opt_maxiter < 1) call die(r, "--opt-maxiter must be at least 1")
      call r%cli%get("opt-gtol", r%opt_gtol, r%err)
      call check(r, "--opt-gtol")
      if (r%opt_gtol <= 0.0_dp) call die(r, "--opt-gtol must be positive")

      ! ---- target, saddle method, connect -------------------------------
      call r%cli%get("opt-target", sval, r%err)
      r%opt_target_arg = to_lower_alloc(trim(sval))
      select case (r%opt_target_arg)
      case ("minimum")
         r%opt_target_kind = TRC_OPT_TARGET_MINIMUM
      case ("saddle")
         r%opt_target_kind = TRC_OPT_TARGET_SADDLE
      case default
         call die(r, "--opt-target must be minimum or saddle, got '"//r%opt_target_arg//"'")
      end select
      r%opt_connect = r%cli%is_set("opt-connect")
      if (r%opt_target_kind == TRC_OPT_TARGET_MINIMUM) then
         if (r%cli%is_set("opt-saddle-method")) call die(r, "--opt-saddle-method needs --opt-target saddle")
         if (r%opt_connect) call die(r, "--opt-connect needs --opt-target saddle")
         if (r%cli%is_set("opt-algorithm")) then
            call r%cli%get("opt-algorithm", sval, r%err)
            r%opt_algorithm_arg = to_lower_alloc(trim(sval))
            select case (r%opt_algorithm_arg)
            case ("sd")
               r%opt_algorithm_kind = TRC_OPT_ALG_SD
            case ("cg-auto")
               r%opt_algorithm_kind = TRC_OPT_ALG_CG_AUTO
            case ("cg")
               r%opt_algorithm_kind = TRC_OPT_ALG_CG
            case ("lbfgs")
               r%opt_algorithm_kind = TRC_OPT_ALG_LBFGS
            case ("prfo")
               r%opt_algorithm_kind = TRC_OPT_ALG_PRFO   ! P-RFO as a plain minimiser: valid, just unusual
            case ("nr")
               r%opt_algorithm_kind = TRC_OPT_ALG_NR
            case ("damped")
               r%opt_algorithm_kind = TRC_OPT_ALG_DAMPED
            case default
               call die(r, "--opt-algorithm must be sd, cg, cg-auto, lbfgs, nr, damped or prfo, got '"// &
                        r%opt_algorithm_arg//"'")
            end select
         else
            r%opt_algorithm_kind = TRC_OPT_ALG_LBFGS
         end if
      else
         ! target = saddle: the saddle method fixes DL-FIND's own
         ! algorithm/coordinate family (cb_get_params) -- --opt-algorithm
         ! is refused rather than silently ignored, so a command line
         ! that sets both never means something different from what it
         ! says.
         if (r%cli%is_set("opt-algorithm")) call die(r, "--opt-algorithm is not used with --opt-target "// &
                                                     "saddle: --opt-saddle-method (prfo or dimer) fixes "// &
                                                     "DL-FIND's own algorithm/coordinate family")
         call r%cli%get("opt-saddle-method", sval, r%err)
         r%opt_saddle_method_arg = to_lower_alloc(trim(sval))
         select case (r%opt_saddle_method_arg)
         case ("prfo")
            r%opt_saddle_method_kind = TRC_OPT_SADDLE_PRFO
         case ("dimer")
            r%opt_saddle_method_kind = TRC_OPT_SADDLE_DIMER
         case default
            call die(r, "--opt-saddle-method must be prfo or dimer, got '"//r%opt_saddle_method_arg//"'")
         end select
         ! --opt-connect with P-RFO crashes INSIDE libdlfind itself, not in
         ! this bridge: confirmed with a debug (-fcheck=all) build of
         ! libdlfind, the fault is a bad array copy at
         ! dlf_formstep.f90:2362 (`array(:)=tsmode_r(:)` in
         ! dlf_formstep_get_ra, called from dlf_task.f90's task=1011 case
         ! right after the TS sub-search converges) -- the size check one
         ! line above it passes, so `tsmode_r` itself is a corrupted
         ! allocation by the time the downhill phase asks for it, not a
         ! size mismatch this bridge could catch or avoid. The dimer
         ! method's own TS-mode bookkeeping does not go through
         ! dlf_formstep_set_tsmode's mode 10/11 path at all and does not
         ! reproduce this: a connect run from the SAME geometry with
         ! --opt-saddle-method dimer finds the TS and both downhill minima
         ! correctly (test/check_opt_connect.F90). Refused here rather
         ! than left to crash the process.
         if (r%opt_connect .and. r%opt_saddle_method_kind == TRC_OPT_SADDLE_PRFO) then
            call die(r, "--opt-connect with --opt-saddle-method prfo is refused: it crashes inside "// &
                     "libdlfind itself (a corrupted TS-mode array in dlf_formstep.f90, confirmed with "// &
                     "a debug build -- not something this bridge can work around). Use "// &
                     "--opt-saddle-method dimer with --opt-connect instead.")
         end if
      end if

      ! ---- Hessian handling (meaningful for prfo/nr; harmless, DL-FIND
      ! simply ignores them, for every other algorithm) -------------------
      call r%cli%get("opt-init-hessian", sval, r%err)
      r%opt_init_hessian_arg = to_lower_alloc(trim(sval))
      select case (r%opt_init_hessian_arg)
      case ("terco")
         r%opt_hessian_source_kind = TRC_OPT_HESS_TERCO
      case ("dlfind")
         r%opt_hessian_source_kind = TRC_OPT_HESS_DLFIND
      case default
         call die(r, "--opt-init-hessian must be terco or dlfind, got '"//r%opt_init_hessian_arg//"'")
      end select
      if (r%cli%is_set("opt-hessian-update")) then
         call r%cli%get("opt-hessian-update", sval, r%err)
         r%opt_hessian_update_arg = to_lower_alloc(trim(sval))
         select case (r%opt_hessian_update_arg)
         case ("none")
            r%opt_hessian_update_kind = TRC_OPT_HUPD_NONE
         case ("powell")
            r%opt_hessian_update_kind = TRC_OPT_HUPD_POWELL
         case ("bofill")
            r%opt_hessian_update_kind = TRC_OPT_HUPD_BOFILL
         case ("bfgs")
            r%opt_hessian_update_kind = TRC_OPT_HUPD_BFGS
         case default
            call die(r, "--opt-hessian-update must be none, powell, bofill or bfgs, got '"// &
                     r%opt_hessian_update_arg//"'")
         end select
      end if
      if (r%cli%is_set("opt-max-updates")) then
         call r%cli%get("opt-max-updates", r%opt_max_updates, r%err)
         call check(r, "--opt-max-updates")
         ! 0 is accepted by DL-FIND itself (dlf_formstep.f90's own
         ! iupd>=maxupd test, iupd starting at 0): it recomputes the
         ! Hessian fresh every step rather than ever updating it, not a
         ! degenerate "no Hessian" case.
         if (r%opt_max_updates < 0) call die(r, "--opt-max-updates must be at least 0")
      end if
      if (r%cli%is_set("opt-zero-modes")) then
         call r%cli%get("opt-zero-modes", r%opt_zero_modes, r%err)
         call check(r, "--opt-zero-modes")
         if (r%opt_zero_modes < 0 .or. r%opt_zero_modes > 3*r%natm) then
            call die(r, "--opt-zero-modes must be between 0 and 3*natm ("//to_string(3*r%natm)// &
                     "), got "//to_string(r%opt_zero_modes))
         end if
      end if
      if (r%cli%is_set("opt-soft-modes")) then
         call r%cli%get("opt-soft-modes", r%opt_soft_modes, r%err)
         call check(r, "--opt-soft-modes")
         if (r%opt_soft_modes <= 0.0_dp) call die(r, "--opt-soft-modes must be positive")
      end if

      ! --opt-init-hessian/--opt-hessian-update/--opt-max-updates:
      ! meaningful only for P-RFO or NR (cb_get_params' own
      ! is_prfo_or_nr branch; DL-FIND never reads inithessian/update/
      ! maxupd for any other algorithm) -- refused rather than silently
      ! ignored, so a command line that sets one means what it says.
      opt_hessian_relevant = (r%opt_target_kind == TRC_OPT_TARGET_SADDLE .and. &
                              r%opt_saddle_method_kind == TRC_OPT_SADDLE_PRFO) .or. &
                             (r%opt_target_kind == TRC_OPT_TARGET_MINIMUM .and. &
                              (r%opt_algorithm_kind == TRC_OPT_ALG_PRFO .or. r%opt_algorithm_kind == TRC_OPT_ALG_NR))
      if (.not. opt_hessian_relevant) then
         if (r%cli%is_set("opt-init-hessian")) call die(r, "--opt-init-hessian needs --opt-algorithm prfo/nr "// &
                  "or --opt-target saddle --opt-saddle-method prfo")
         if (r%cli%is_set("opt-hessian-update")) call die(r, "--opt-hessian-update needs --opt-algorithm "// &
                  "prfo/nr or --opt-target saddle --opt-saddle-method prfo")
         if (r%cli%is_set("opt-max-updates")) call die(r, "--opt-max-updates needs --opt-algorithm prfo/nr "// &
                  "or --opt-target saddle --opt-saddle-method prfo")
      end if

      ! --opt-zero-modes/--opt-soft-modes: meaningful for P-RFO, NR or
      ! the dimer method (cb_get_params' own nzero/soft handling) --
      ! refused the same way for every algorithm that never reads them.
      opt_zero_soft_relevant = (r%opt_target_kind == TRC_OPT_TARGET_SADDLE) .or. &
                               (r%opt_target_kind == TRC_OPT_TARGET_MINIMUM .and. &
                                (r%opt_algorithm_kind == TRC_OPT_ALG_PRFO .or. r%opt_algorithm_kind == TRC_OPT_ALG_NR))
      if (.not. opt_zero_soft_relevant) then
         if (r%cli%is_set("opt-zero-modes")) call die(r, "--opt-zero-modes needs --opt-algorithm prfo/nr "// &
                  "or --opt-target saddle")
         if (r%cli%is_set("opt-soft-modes")) call die(r, "--opt-soft-modes needs --opt-algorithm prfo/nr "// &
                  "or --opt-target saddle")
      end if

      ! A single atom has no internal degrees of freedom: P-RFO/NR/the
      ! dimer method in Cartesian coordinates default nzero to 6 (5 if
      ! linear, never reached with one atom), which already exceeds
      ! nvar=3 before any eigenvector-following step could even be
      ! attempted -- refused here rather than left to DL-FIND, which has
      ! no graceful way to decline a run it cannot do.
      if (r%natm == 1 .and. r%opt_coords_kind == TRC_OPT_COORDS_CARTESIAN .and. opt_zero_soft_relevant) then
         call die(r, "--opt-algorithm prfo/nr and --opt-target saddle need more than one atom in Cartesian "// &
                  "coordinates: a single atom has no internal degrees of freedom to search over")
      end if
      if (r%cli%is_set("opt-zero-modes") .and. r%opt_zero_modes > 3*r%natm - 1) then
         call die(r, "--opt-zero-modes "//to_string(r%opt_zero_modes)//" leaves fewer than one true "// &
                  "degree of freedom for a "//to_string(r%natm)//"-atom molecule (3*natm = "// &
                  to_string(3*r%natm)//")")
      end if

      ! ---- step/trust controls -------------------------------------------
      if (r%cli%is_set("opt-max-step")) then
         call r%cli%get("opt-max-step", r%opt_max_step, r%err)
         call check(r, "--opt-max-step")
         if (r%opt_max_step <= 0.0_dp) call die(r, "--opt-max-step must be positive")
      end if
      call r%cli%get("opt-linesearch", sval, r%err)
      r%opt_linesearch_arg = to_lower_alloc(trim(sval))
      select case (r%opt_linesearch_arg)
      case ("auto")
         r%opt_linesearch_kind = TRC_OPT_LINE_AUTO
      case ("scaling")
         r%opt_linesearch_kind = TRC_OPT_LINE_SCALING
      case ("energy-trust")
         r%opt_linesearch_kind = TRC_OPT_LINE_ENERGY_TRUST
      case ("gradient-trust")
         r%opt_linesearch_kind = TRC_OPT_LINE_GRADIENT_TRUST
      case default
         call die(r, "--opt-linesearch must be auto, scaling, energy-trust or gradient-trust, got '"// &
                  r%opt_linesearch_arg//"'")
      end select
      if (r%cli%is_set("opt-scalestep")) then
         call r%cli%get("opt-scalestep", r%opt_scalestep, r%err)
         call check(r, "--opt-scalestep")
         if (r%opt_scalestep <= 0.0_dp) call die(r, "--opt-scalestep must be positive")
      end if
      if (r%cli%is_set("opt-lbfgs-memory")) then
         call r%cli%get("opt-lbfgs-memory", r%opt_lbfgs_memory, r%err)
         call check(r, "--opt-lbfgs-memory")
         if (r%opt_lbfgs_memory < 1) call die(r, "--opt-lbfgs-memory must be at least 1")
      end if
      if (r%cli%is_set("opt-etol")) then
         call r%cli%get("opt-etol", r%opt_etol, r%err)
         call check(r, "--opt-etol")
         if (r%opt_etol <= 0.0_dp) call die(r, "--opt-etol must be positive")
      end if
      if (r%cli%is_set("opt-print-level")) then
         call r%cli%get("opt-print-level", r%opt_print_level, r%err)
         call check(r, "--opt-print-level")
         if (all(r%opt_print_level /= [0, 2, 4, 6])) call die(r, "--opt-print-level must be 0, 2, 4 or 6")
      end if
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
   ! build has nothing to bind either way. When trc_bind_device_env above
   ! already bound this device before MPI_Init, this picks the same one and
   ! changes nothing.
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
   call api(r, trc_set_nlc_grid(r%h, int(r%nlc_grid, c_int)))
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

   if (r%opt) then
      allocate (opt_residues(r%natm))
      if (r%opt_coords_kind == TRC_OPT_COORDS_HDLC .or. r%opt_coords_kind == TRC_OPT_COORDS_HDLC_TC) then
         call build_opt_residues(r%natm, r%z, xyz, opt_residues, r%err)
         if (r%err%has_error()) call die(r, r%err%get_message())
      else
         opt_residues = 1   ! cartesian/dlc/dlc-tc: trc_dlfind_bridge's own driver needs no real partition
      end if
      r%opt_is_linear = opt_geometry_is_linear(r%natm, r%z, xyz, herr)
      if (herr%has_error()) call die(r, "--opt: "//herr%get_message())

      opt_opts%coords_kind = r%opt_coords_kind
      opt_opts%guess_kind = r%opt_guess_kind
      opt_opts%fallback_guess_kind = int(guess_kind, c_int)
      opt_opts%maxiter = r%opt_maxiter
      opt_opts%gtol = r%opt_gtol
      opt_opts%algorithm = r%opt_algorithm_kind
      opt_opts%target_kind = r%opt_target_kind
      opt_opts%saddle_method = r%opt_saddle_method_kind
      opt_opts%connect = r%opt_connect
      opt_opts%is_linear = r%opt_is_linear
      opt_opts%hessian_source = r%opt_hessian_source_kind
      opt_opts%hessian_update = r%opt_hessian_update_kind
      opt_opts%max_updates = r%opt_max_updates
      opt_opts%zero_modes = r%opt_zero_modes
      opt_opts%soft_modes = r%opt_soft_modes
      opt_opts%max_step = r%opt_max_step
      opt_opts%linesearch = r%opt_linesearch_kind
      opt_opts%scalestep = r%opt_scalestep
      opt_opts%lbfgs_memory = r%opt_lbfgs_memory
      opt_opts%etol = r%opt_etol
      opt_opts%print_level = r%opt_print_level

      ! <geometry base>_opt.xyz/_traj.xyz's own path, and --json's (if
      ! given): trc_dlfind_optimize needs these BEFORE it runs, not after,
      ! because an unrecoverable mid-run SCF failure writes them from
      ! inside a DL-FIND callback (trc_dlfind_bridge%fail_eval) -- DL-FIND
      ! itself stops the process right after, so nothing after this call
      ! returns would ever get the chance.
      opt_base = r%geom
      if (len(opt_base) > 4) then
         if (to_lower_alloc(opt_base(len(opt_base) - 3:)) == ".xyz") opt_base = opt_base(1:len(opt_base) - 4)
      end if
      ! Written into the CURRENT WORKING DIRECTORY, not next to the input
      ! geometry: a directory part in r%geom (an absolute path, or one
      ! outside the CWD) must not steer these writes somewhere a test's
      ! own WORKING_DIRECTORY (ctest's convention, test/CMakeLists.txt)
      ! does not expect -- DL-FIND's own path*.xyz already go to the CWD
      ! (dl-find.f90), so this matches.
      opt_base_slash = scan(opt_base, "/", back=.true.)
      if (opt_base_slash > 0) opt_base = opt_base(opt_base_slash + 1:)
      opt_opts%opt_xyz_path = opt_base//"_opt.xyz"
      opt_opts%traj_xyz_path = opt_base//"_traj.xyz"
      opt_opts%json_path = r%json
      opt_opts%geometry_path = r%geom
      opt_opts%coords_name = r%opt_coords_arg
      opt_opts%guess_name = r%opt_guess_arg

      call trc_dlfind_optimize(r%h, r%natm, r%z, opt_residues, xyz, opt_opts, opt_res)
      if (opt_res%failed) call die(r, "--opt: "//opt_res%message)

      r%opt_converged = opt_res%converged
      r%opt_steps = opt_res%steps
      r%opt_evals = opt_res%energy_evaluations
      r%opt_hess_evals = opt_res%hessian_evaluations
      r%opt_gmax = opt_res%final_gmax
      r%opt_grms = opt_res%final_grms
      r%opt_n_traj = opt_res%n_traj
      if (allocated(opt_res%traj_energy)) r%opt_traj_energy = opt_res%traj_energy
      r%opt_have_ts = opt_res%have_ts
      if (opt_res%have_ts) then
         r%opt_ts_energy = opt_res%ts_energy
         r%opt_ts_xyz = opt_res%ts_coords
         r%opt_ts_gmax = opt_res%ts_gmax
         r%opt_ts_converged = opt_res%ts_converged
      end if
      r%opt_have_min_plus = opt_res%have_minimum_plus
      if (opt_res%have_minimum_plus) then
         r%opt_min_plus_energy = opt_res%minimum_plus_energy
         r%opt_min_plus_xyz = opt_res%minimum_plus_coords
         r%opt_min_plus_gmax = opt_res%minimum_plus_gmax
         r%opt_min_plus_converged = opt_res%minimum_plus_converged
      end if
      r%opt_have_min_minus = opt_res%have_minimum_minus
      if (opt_res%have_minimum_minus) then
         r%opt_min_minus_energy = opt_res%minimum_minus_energy
         r%opt_min_minus_xyz = opt_res%minimum_minus_coords
         r%opt_min_minus_gmax = opt_res%minimum_minus_gmax
         r%opt_min_minus_converged = opt_res%minimum_minus_converged
      end if
      xyz = opt_res%final_coords
      r%opt_final_xyz = xyz
      r%grad = .true.   ! so the existing --grad path below reports the final gradient too

      ! Re-converge the context at the reported final (ACCEPTED) geometry:
      ! cb_get_gradient calls trc_set_geometry/trc_run_scf for EVERY
      ! evaluation, including a rejected trial that happens after the last
      ! accepted step (e.g. --opt-maxiter stops the run right after one),
      ! so by the time api_dl_find returns, the context may be sitting at
      ! that later trial's geometry/SCF, not opt_res%final_coords's own.
      ! Every getter below must report exactly the geometry --opt writes
      ! out.
      call api(r, trc_set_geometry(r%h, reshape(real(xyz, c_double), [3*r%natm])))
      st = trc_run_scf(r%h)
   else if (len(r%restart) > 0) then
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
      ! trc_run_gradient's own refusals (RI-MP2, PCM) come
      ! back through `api`, which prints the library's own message and
      ! stops -- the same clear-refusal path every other unsupported
      ! combination in this CLI already takes. An LDA/GGA --func (global
      ! and range-separated hybrids and VV10 included) is NOT refused:
      ! trc_run_gradient builds the KS-DFT gradient at the SCF's own
      ! functional and grids, and with --df the fitted one.
      call api(r, trc_run_gradient(r%h))
      allocate (gradbuf(3*r%natm), r%gradient(3, r%natm))
      call api(r, trc_gradient(r%h, gradbuf))
      r%gradient = reshape(real(gradbuf, dp), [3, r%natm])
      deallocate (gradbuf)
   end if
   if (r%hessian .and. r%conv /= 0) then
      ! trc_run_hessian's own refusals (meta-GGA, --rimp2, a double
      ! hybrid, --pcm, --df outside its own scope, more than one rank)
      ! come back through `api`, exactly as --grad's do -- checked
      ! before any displacement, so a refused run costs nothing beyond
      ! the SCF this CLI already ran for --grad/the energy anyway.
      call api(r, trc_set_hessian_step(r%h, real(r%hess_step, c_double)))
      call api(r, trc_run_hessian(r%h))
      ndof = 3*r%natm
      allocate (hessbuf(ndof*ndof), r%hess(ndof, ndof))
      call api(r, trc_hessian(r%h, hessbuf))
      r%hess = reshape(real(hessbuf, dp), [ndof, ndof])
      deallocate (hessbuf)

      call build_masses(r%z, r%masses_kind, mass, herr)
      if (herr%has_error()) call die(r, "--masses: "//herr%get_message())
      call trc_vibrational_analysis(mass, xyz, r%hess, r%freq_cm1, r%redmass_amu, r%force_const, &
                                    displacement, r%is_linear, r%n_tr, herr)
      if (herr%has_error()) call die(r, "--hessian: vibrational analysis: "//herr%get_message())
      call trc_thermochemistry(mass, xyz, r%freq_cm1, real(r%etot, dp), r%sym_number, r%mult, r%thermo, herr, &
                               temperature_k=r%temperature_k, pressure_pa=r%pressure_atm*ATM_TO_PA)
      if (herr%has_error()) call die(r, "--hessian: thermochemistry: "//herr%get_message())
      r%hess_ok = .true.
   end if

   call system_clock(r%t1)
   if (r%rank0) call summary(r)
   if (r%rank0 .and. r%grad .and. r%conv /= 0) call gradient_summary(r)
   if (r%rank0 .and. r%opt) call opt_summary(r)
   if (r%rank0 .and. r%hessian .and. r%hess_ok) call hessian_summary(r)
   if (r%rank0 .and. len(r%json) > 0) call write_json(r, r%json)
   if (r%rank0 .and. r%opt) then
      ! opt_base (and opt_opts%opt_xyz_path/traj_xyz_path built from it)
      ! was already computed before trc_dlfind_optimize ran, above.
      call write_opt_xyz(r, opt_base//"_opt.xyz", xyz)
      if (r%opt_n_traj > 0) call write_traj_xyz(r, opt_base//"_traj.xyz", opt_res%traj_coords, r%opt_traj_energy)
      call global_logger%info("   wrote "//opt_base//"_opt.xyz and "//opt_base//"_traj.xyz")
      if (r%opt_connect) then
         if (r%opt_have_ts) then
            call write_labelled_xyz(r, opt_base//"_ts.xyz", r%opt_ts_xyz, r%opt_ts_energy, "transition state")
         end if
         if (r%opt_have_min_minus) then
            call write_labelled_xyz(r, opt_base//"_minus.xyz", r%opt_min_minus_xyz, r%opt_min_minus_energy, &
                                    "downhill minimum (-distort)")
         end if
         if (r%opt_have_min_plus) then
            call write_labelled_xyz(r, opt_base//"_plus.xyz", r%opt_min_plus_xyz, r%opt_min_plus_energy, &
                                    "downhill minimum (+distort)")
         end if
         call global_logger%info("   opt-connect: wrote "//opt_base//"_ts.xyz, "//opt_base//"_minus.xyz and "// &
                                 opt_base//"_plus.xyz")
      end if
   end if
   st = trc_destroy(r%h)

   if (r%conv == 0) then
      if (r%rank0) call global_logger%error("terco: the SCF did not converge in "//to_string(int(r%niter))//" iterations")
      call quit(r, 2)
   end if
   if (r%opt .and. .not. r%opt_converged) then
      ! Matches the SCF's own non-convergence above: a CLI run that did
      ! not reach the thing it was asked to do must not exit 0 and claim
      ! "normal termination" -- opt_summary already printed "NOT
      ! converged" above, this just makes the exit code agree with it.
      if (r%rank0) call global_logger%error("terco: --opt did not converge in "// &
                                             to_string(r%opt_steps)//" accepted step(s)")
      call quit(r, 2)
   end if
   if (r%rank0) then
      call global_logger%info("")
      call global_logger%info(" normal termination of terco")
   end if
   ! A line of pic's collected lore, as metalquicha prints at the end of a
   ! successful run: rank 0 only, and not under --silent, whose promise is
   ! the summary and nothing else.
   if (r%rank0 .and. .not. r%cli%is_set("silent")) call get_knowledge()
   call quit(r, 0)

end program terco_main
