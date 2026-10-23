!
! terco's side of DL-FIND: fills libdlfind's `spec`/parameter block, turns
! its seven callbacks into calls on a terco context handle (trc_set_geometry
! / trc_run_scf / trc_run_gradient), and turns the result back into a
! trc_opt_result_t. The seven-callback contract, the `spec` layout and its
! traps below come from reading libdlfind's own api.f90 and
! dlf_parameter_module.f90 directly; each one is explained inline where it
! applies rather than pointed at elsewhere.
!
! THIS MODULE IS NEVER COMPILED INTO libterco. It lives in app/ and is
! linked only into the terco_cli executable (app/CMakeLists.txt), which is
! what keeps DL-FIND's LGPL-3.0-or-later licence and its `error stop` on a
! failure (dlf_error, below) out of the library -- see
! cmake/modules/Findlibdlfind.cmake and the root CMakeLists.txt's
! TERCO_ENABLE_DLFIND option.
!
! ONE OPTIMIZATION AT A TIME. DL-FIND's callbacks are plain function
! pointers with no user-data argument (api_dl_find takes them as
! type(c_funptr), nothing else), so the system being optimized has to live
! in module `save` state where the callbacks can reach it: the context
! handle, counters and latest result below, all reset at the start of each
! run so nothing from a previous optimisation in the same process leaks
! into the next one. `active` is the guard against two running at once.
!
! CONVERGENCE AUTHORITY. This bridge reports its own
! verdict -- the max Cartesian gradient component of the LAST accepted
! step against `opts%gtol` -- not DL-FIND's internal one. DL-FIND tests
! five criteria (max/RMS gradient, max/RMS step, energy change) and, in
! HDLC/DLC, tests the INTERNAL-coordinate gradient, which is not the
! number this CLI prints or writes to --json either way; reading it back
! out reliably across coordinate systems would need reaching into
! dlf_global state this bridge does not otherwise touch. The Cartesian
! max-gradient check is also what `--opt-gtol` documents it means.
!
module trc_dlfind_bridge
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_funptr, c_funloc, c_loc, &
                                          c_null_ptr, c_associated
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   use trc_capi, only: trc_set_geometry, trc_run_scf, trc_run_gradient, trc_gradient, &
                       trc_energy, trc_converged, trc_iterations, trc_nspin, trc_nao, &
                       trc_set_guess, trc_density, trc_message, &
                       TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_GIVEN, TRC_GUESS_PROJECTED
   use terco_app, only: f_string
   implicit none
   private

   public :: trc_dlfind_available
   public :: TRC_OPT_COORDS_CARTESIAN, TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_DLC
   public :: TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, TRC_OPT_GUESS_NONE
   public :: trc_opt_options_t, trc_opt_result_t
   public :: trc_dlfind_optimize

   integer, parameter :: TRC_OPT_COORDS_CARTESIAN = 0, TRC_OPT_COORDS_HDLC = 1, TRC_OPT_COORDS_DLC = 2
      !! This module's own coordinate-system enum, in app/main.F90's
      !! --opt-coords; translated to DL-FIND's icoord numbering (which is
      !! not contiguous -- DLC is 3, not 2) only inside cb_get_params.

   integer, parameter :: TRC_OPT_GUESS_DENSITY = 0, TRC_OPT_GUESS_PROJECTED = 1, TRC_OPT_GUESS_NONE = 2
      !! --opt-guess: the previous step's converged density (TRC_GUESS_GIVEN,
      !! fetched with trc_density before the geometry moves) or its
      !! orbitals re-orthonormalised in the new overlap (TRC_GUESS_PROJECTED,
      !! a338056's warm start -- needs no density fetch at all, see
      !! cb_get_gradient's own comment). TRC_OPT_GUESS_NONE re-runs
      !! `fallback_guess_kind` (ordinarily SAD/SADQ) at every step, the
      !! mqc baseline this PR's warm-start measurement compares the other
      !! two against; not reachable from --opt-guess (only "density" and
      !! "projected" are), used only by the measurement itself.

   ! DL-FIND's own numbering (api_dl_find's icoord/iopt/switch) -- this
   ! module is the only place terco needs to know it; see libdlfind's
   ! src/api.f90 and src/dlf_parameter_module.f90 for the full table this
   ! is a small subset of (no NEB, no dimer, no P-RFO: out of scope here).
   integer(c_int), parameter :: DLF_COORDS_CARTESIAN = 0, DLF_COORDS_HDLC = 1, DLF_COORDS_DLC = 3
   integer(c_int), parameter :: DLF_OPT_LBFGS = 3
      !! DL-FIND's own engine default (dlf_default_set's inline default
      !! for iopt); this PR offers no --opt-algorithm, so it is the
      !! only one wired up. The code is shaped so a later PR can add one
      !! (trc_opt_options_t gains a field, cb_get_params's `iopt` branches
      !! on it) without touching anything else here.
   integer(c_int), parameter :: DLF_INITHESS_TWO_POINT = 2
      !! No analytic Hessian is offered to DL-FIND (cb_get_hessian always
      !! declines) -- only P-RFO/NR read this, and LBFGS is the only
      !! algorithm this bridge drives, so it is set but never consulted.
   integer(c_int), parameter :: DLF_PUT_STEP = 1, DLF_PUT_TS = 3
   integer, parameter :: MSG_LEN = 1024
      !! trc_message's own buffer length convention (app/main.F90's `api`)

   type :: trc_opt_options_t
      integer :: coords_kind = TRC_OPT_COORDS_CARTESIAN
      integer :: guess_kind = TRC_OPT_GUESS_DENSITY
      integer(c_int) :: fallback_guess_kind = 2_c_int
         !! trc_capi's TRC_GUESS_* (default 2 = TRC_GUESS_SAD) -- what a
         !! warm-started SCF that failed to converge retries from once,
         !! before this step is declared failed. The CLI sets this to
         !! whatever --guess already resolved (sad/sadq/gwh/core/sap), so
         !! the retry lands on the same guess a non-optimisation run of
         !! this molecule would have used.
      integer :: maxiter = 100
      real(dp) :: gtol = 4.5e-4_dp
   end type trc_opt_options_t

   type :: trc_opt_result_t
      logical :: converged = .false.
      logical :: failed = .false.
      character(len=:), allocatable :: message
      integer :: steps = 0               !! accepted geometries (cb_put_coords, switch 1)
      integer :: energy_evaluations = 0  !! every cb_get_gradient call, line-search trials included
      real(dp) :: final_energy = 0.0_dp
      real(dp) :: final_gmax = 0.0_dp
      real(dp) :: final_grms = 0.0_dp
      real(dp), allocatable :: final_coords(:, :)     !! (3, natm), Bohr
      integer :: n_traj = 0
      real(dp), allocatable :: traj_energy(:)
      real(dp), allocatable :: traj_coords(:, :, :)   !! (3, natm, n_traj), Bohr
      integer :: total_scf_iterations = 0
         !! Sum of trc_iterations over every ACCEPTED step -- the single
         !! number the warm-start measurement (--opt-guess density vs
         !! projected vs the mqc baseline, trc_opt_options_t's own
         !! TRC_OPT_GUESS_NONE) compares.
   end type trc_opt_result_t

   interface
      subroutine api_dl_find(nvarin, nvarin2, nspec, master, &
                             c_dlf_error, c_dlf_get_gradient, c_dlf_get_hessian, &
                             c_dlf_get_multistate_gradients, c_dlf_get_params, &
                             c_dlf_put_coords, c_dlf_update) bind(C, name="api_dl_find")
         !! libdlfind's own single entry point, declared here by hand
         !! rather than taken from its .mod -- see
         !! cmake/modules/Findlibdlfind.cmake's own header for why (the
         !! two projects may be built by different Fortran compilers).
         import :: c_int, c_funptr
         implicit none
         integer(c_int), value :: nvarin   !! 3*natm
         integer(c_int), value :: nvarin2  !! second coordinate array's length; 0, nothing offered
         integer(c_int), value :: nspec    !! length of the integer `spec` array
         integer(c_int), value :: master   !! 1: this task leads (and is the only task in) the run
         type(c_funptr), value :: c_dlf_error, c_dlf_get_gradient, c_dlf_get_hessian, &
                                  c_dlf_get_multistate_gradients, c_dlf_get_params, &
                                  c_dlf_put_coords, c_dlf_update
      end subroutine api_dl_find
   end interface

   ! --- module `save` state the seven callbacks share -----------------------
   logical, save :: active = .false.
   type(c_ptr), save :: ctx = c_null_ptr
   integer, save :: natm = 0
   integer, allocatable, save :: residues(:), znuc(:)
   real(dp), allocatable, save :: start_coords(:, :)
   type(trc_opt_options_t), save :: cur_opts
   logical, save :: have_prev_scf = .false.
   integer(c_int), save :: cached_nspin = 1_c_int
   logical, save :: eval_failed = .false.
   character(len=512), save :: eval_fail_msg = ""
   integer, save :: n_evals = 0, n_steps = 0, total_scf_iters = 0
   logical, save :: have_result = .false.
   real(dp), save :: latest_energy = 0.0_dp
   real(dp), allocatable, save :: latest_coords(:, :)
   real(dp), save :: last_gmax = 0.0_dp, last_grms = 0.0_dp
   integer, save :: last_niter = 0
   type(timer_type), save :: step_clock
   logical, save :: step_clock_running = .false.
   integer, save :: traj_cap = 0, traj_n = 0
   real(dp), allocatable, save :: traj_e(:)
   real(dp), allocatable, save :: traj_xyz(:, :, :)

contains

   pure function trc_dlfind_available() result(ok)
      !! Whether THIS build can optimise a geometry -- always true here;
      !! app/trc_dlfind_stub.F90's own answer (the TERCO_ENABLE_DLFIND=OFF
      !! build) is always false, which is all app/main.F90 ever asks.
      logical :: ok
      ok = .true.
   end function trc_dlfind_available

   subroutine trc_dlfind_optimize(handle, n, z, residues_in, coords, opts, res)
      !! Minimise the molecule already set up on `handle` (molecule, basis
      !! and method must already match `coords`/`z` -- this routine only
      !! ever moves atoms on it, through trc_set_geometry, never builds
      !! anything). The optimised geometry, energy and trajectory come
      !! back in `res`; `coords` itself is read-only (the starting point
      !! cb_get_params hands DL-FIND).
      type(c_ptr), intent(in) :: handle
      integer, intent(in) :: n
      integer, intent(in) :: z(n)
      integer, intent(in) :: residues_in(n)
         !! 1-based residue number per atom (connected components for
         !! --opt-coords hdlc, all 1 for cartesian/dlc -- app/main.F90's
         !! own job to build, trc_constants%trc_covalent_radius-based).
      real(dp), intent(in) :: coords(3, n)
      type(trc_opt_options_t), intent(in) :: opts
      type(trc_opt_result_t), intent(out) :: res

      integer(c_int) :: nvar, nvar2, nspec

      if (active) then
         res%failed = .true.
         res%message = "trc_dlfind_optimize: an optimisation is already running on this process"
         return
      end if
      active = .true.

      ctx = handle
      natm = n
      if (allocated(znuc)) deallocate (znuc)
      if (allocated(residues)) deallocate (residues)
      if (allocated(start_coords)) deallocate (start_coords)
      allocate (znuc(n), source=z)
      allocate (residues(n), source=residues_in)
      allocate (start_coords(3, n), source=coords)
      cur_opts = opts
      have_prev_scf = .false.
      cached_nspin = 1
      eval_failed = .false.
      eval_fail_msg = ""
      n_evals = 0
      n_steps = 0
      total_scf_iters = 0
      have_result = .false.
      latest_energy = 0.0_dp
      if (allocated(latest_coords)) deallocate (latest_coords)
      allocate (latest_coords(3, n), source=coords)
      last_gmax = 0.0_dp
      last_grms = 0.0_dp
      last_niter = 0
      step_clock_running = .false.

      ! Upper bound on accepted steps: DL-FIND itself enforces maxcycle,
      ! so this never has to grow.
      traj_cap = max(opts%maxiter + 4, 8)
      if (allocated(traj_e)) deallocate (traj_e)
      if (allocated(traj_xyz)) deallocate (traj_xyz)
      allocate (traj_e(traj_cap), traj_xyz(3, n, traj_cap))
      traj_n = 0

      nvar = int(3*n, c_int)
      nvar2 = 0_c_int   !! no second structure (no NEB endpoint, no dimer)
      nspec = int(3*n, c_int)   !! [residue(n), znuc(n), micspec(n)]: no constraints, so no icons block

      call api_dl_find(nvar, nvar2, nspec, 1_c_int, &
                       c_funloc(cb_error), c_funloc(cb_get_gradient), c_funloc(cb_get_hessian), &
                       c_funloc(cb_get_multistate_gradients), c_funloc(cb_get_params), &
                       c_funloc(cb_put_coords), c_funloc(cb_update))

      res%steps = n_steps
      res%energy_evaluations = n_evals
      res%total_scf_iterations = total_scf_iters
      if (eval_failed) then
         res%failed = .true.
         res%message = "geometry optimisation stopped: "//trim(eval_fail_msg)
      else if (.not. have_result) then
         ! DL-FIND reports a geometry through dlf_put_coords only when
         ! printf >= 1 -- leaving it at DL-FIND's own default of 0 lets a
         ! run converge with no geometry ever reported back, which is why
         ! cb_get_params always sets it to 1 below. Reaching here with
         ! have_result still false means api_dl_find returned without
         ! ever calling it.
         res%failed = .true.
         res%message = "DL-FIND returned without reporting a geometry"
      else
         res%final_energy = latest_energy
         res%final_coords = latest_coords
         res%final_gmax = last_gmax
         res%final_grms = last_grms
         res%converged = last_gmax <= opts%gtol
      end if
      res%n_traj = traj_n
      if (traj_n > 0) then
         res%traj_energy = traj_e(1:traj_n)
         res%traj_coords = traj_xyz(:, :, 1:traj_n)
      end if

      active = .false.
      ctx = c_null_ptr
      if (allocated(znuc)) deallocate (znuc)
      if (allocated(residues)) deallocate (residues)
      if (allocated(start_coords)) deallocate (start_coords)
      if (allocated(latest_coords)) deallocate (latest_coords)
      if (allocated(traj_e)) deallocate (traj_e)
      if (allocated(traj_xyz)) deallocate (traj_xyz)
      natm = 0
   end subroutine trc_dlfind_optimize

   subroutine cb_get_params( &
      nvar, nvar2, nspec, coords, coords2, spec, ierr, tolerance, printl, maxcycle, maxene, tatoms, icoord, iopt, iline, &
      maxstep, scalestep, lbfgs_mem, nimage, nebk, dump, restart, nz, ncons, nconn, update, maxupd, delta, soft, &
      inithessian, carthessian, tsrel, maxrot, tolrot, nframe, nmass, nweight, timestep, fric0, fricfac, fricp, &
      imultistate, state_i, state_j, pf_c1, pf_c2, gp_c3, gp_c4, ln_t1, ln_t2, printf, tolerance_e, distort, massweight, &
      minstep, maxdump, task, temperature, po_pop_size, po_radius, po_contraction, po_tolerance_r, po_tolerance_g, &
      po_distribution, po_maxcycle, po_init_pop_size, po_reset, po_mutation_rate, po_death_rate, po_scalefac, po_nsave, &
      ntasks, tdlf_farm, n_po_scaling, neb_climb_test, neb_freeze_test, nzero, coupled_states, qtsflag, imicroiter, &
      maxmicrocycle, micro_esp_fit) bind(c)
      !! DL-FIND asking for the start geometry and every setting this run
      !! cares about. Every argument arrives already holding DL-FIND's own
      !! default (dlf_default_init ran first), so only what this bridge
      !! actually sets below is touched.
      integer(c_int), intent(in), value :: nvar, nvar2, nspec
      real(c_double), intent(inout) :: coords(nvar), coords2(nvar2)
      integer(c_int), intent(inout) :: spec(nspec)
      integer(c_int), intent(out) :: ierr
      real(c_double), intent(inout) :: tolerance
      integer(c_int), intent(inout) :: printl, maxcycle, maxene, tatoms, icoord, iopt, iline
      real(c_double), intent(inout) :: maxstep, scalestep
      integer(c_int), intent(inout) :: lbfgs_mem, nimage
      real(c_double), intent(inout) :: nebk
      integer(c_int), intent(inout) :: dump, restart, nz, ncons, nconn, update, maxupd
      real(c_double), intent(inout) :: delta, soft
      integer(c_int), intent(inout) :: inithessian, carthessian, tsrel, maxrot
      real(c_double), intent(inout) :: tolrot
      integer(c_int), intent(inout) :: nframe, nmass, nweight
      real(c_double), intent(inout) :: timestep, fric0, fricfac, fricp
      integer(c_int), intent(inout) :: imultistate, state_i, state_j
      real(c_double), intent(inout) :: pf_c1, pf_c2, gp_c3, gp_c4, ln_t1, ln_t2
      integer(c_int), intent(inout) :: printf
      real(c_double), intent(inout) :: tolerance_e, distort
      integer(c_int), intent(inout) :: massweight
      real(c_double), intent(inout) :: minstep
      integer(c_int), intent(inout) :: maxdump, task
      real(c_double), intent(inout) :: temperature
      integer(c_int), intent(inout) :: po_pop_size
      real(c_double), intent(inout) :: po_radius, po_contraction, po_tolerance_r, po_tolerance_g
      integer(c_int), intent(inout) :: po_distribution, po_maxcycle, po_init_pop_size, po_reset
      real(c_double), intent(inout) :: po_mutation_rate, po_death_rate, po_scalefac
      integer(c_int), intent(inout) :: po_nsave, ntasks, tdlf_farm, n_po_scaling
      real(c_double), intent(inout) :: neb_climb_test, neb_freeze_test
      integer(c_int), intent(inout) :: nzero, coupled_states, qtsflag
      integer(c_int), intent(inout) :: imicroiter, maxmicrocycle, micro_esp_fit

      integer :: i

      ierr = 0_c_int

      coords(1:nvar) = reshape(real(start_coords, c_double), [3*natm])

      ! `spec` = [residue(natm), znuc(natm), micspec(natm)] -- no constraints
      ! block, so no 5*ncons entries between them. The FIRST block is
      ! residue numbers, never a frozen/free flag: a positive entry is the
      ! residue the atom belongs to (all 1 means "ignore residues", which
      ! Cartesian/DLC do; HDLC needs the real partition -- an all-zero
      ! residue block breaks HDLC outright, avoided here because
      ! residues(:) is never zero by construction, app/main.F90's own
      ! build_opt_residues).
      do i = 1, natm
         spec(i) = int(residues(i), c_int)
      end do
      do i = 1, natm
         spec(natm + i) = int(znuc(i), c_int)
      end do
      spec(2*natm + 1:3*natm) = 1_c_int   !! micspec: DL-FIND requires 0 or 1 here; 1 everywhere

      nz = int(natm, c_int)       !! lets HDLC/DLC perceive connectivity on their own
      ncons = 0_c_int
      nconn = 0_c_int             !! DL-FIND reads nconn and micspec at inconsistent offsets when both are present; 0 avoids that
      tatoms = 1_c_int            !! the variables are atoms, not arbitrary DOF -- HDLC needs to know

      ! Arrive uninitialised, or at -1 and then checked positive, so each
      ! needs an explicit value even though it looks optional.
      nframe = 0_c_int
      nmass = 0_c_int
      nweight = 0_c_int
      ntasks = 1_c_int

      select case (cur_opts%coords_kind)
      case (TRC_OPT_COORDS_HDLC)
         icoord = DLF_COORDS_HDLC
      case (TRC_OPT_COORDS_DLC)
         icoord = DLF_COORDS_DLC
      case default
         icoord = DLF_COORDS_CARTESIAN
      end select
      iopt = DLF_OPT_LBFGS
      inithessian = DLF_INITHESS_TWO_POINT   !! never consulted: LBFGS holds no Hessian

      tolerance = real(cur_opts%gtol, c_double)
      maxcycle = int(cur_opts%maxiter, c_int)

      ! printf >= 1 is what makes DL-FIND call dlf_put_coords at all --
      ! left at DL-FIND's own default of 0, a run converges silently with
      ! no geometry ever reported back. printl follows this process's own
      ! logger level, so --silent also quiets DL-FIND's stdout convergence
      ! table.
      printf = 1_c_int
      printl = dlfind_print_level()

      dump = 0_c_int      !! no checkpoint files -- two optimisations sharing a CWD must not read each other's
      restart = 0_c_int

   end subroutine cb_get_params

   subroutine cb_get_gradient(nvar, coords, energy, gradient, iimage, kiter, status) bind(c)
      !! DL-FIND asking for the energy and gradient at a geometry -- an
      !! accepted step OR a line-search trial point; cb_put_coords, not
      !! this, is where an ACCEPTED step is recorded (line-search trials
      !! hit this callback too, so counting evaluations here would
      !! overcount accepted steps).
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: energy
      real(c_double), intent(out) :: gradient(nvar)
      integer(c_int), intent(in), value :: iimage   !! image index; always 1, no NEB here
      integer(c_int), intent(in), value :: kiter    !! microiteration flag; unused
      integer(c_int), intent(out) :: status

      real(dp) :: xyz(3, natm)
      real(c_double), allocatable, target :: dens_buf(:)
      integer(c_int) :: st, conv, niter, nao_c, nspin_c
      real(c_double) :: e
      real(c_double), allocatable :: gbuf(:)
      logical :: converged_now
      character(kind=c_char) :: msg_buf(MSG_LEN)

      n_evals = n_evals + 1
      xyz = reshape(real(coords, dp), [3, natm])

      ! The warm-start density has to be read BEFORE the geometry moves:
      ! trc_set_geometry invalidates the context's scf_ok, and trc_density
      ! refuses once that happens. TRC_GUESS_PROJECTED needs no such fetch
      ! -- it reads the context's own last-converged orbitals internally
      ! (a338056's trc_scf_context%project_density), which invalidate_runs
      ! never clears.
      if (have_prev_scf .and. cur_opts%guess_kind == TRC_OPT_GUESS_DENSITY) then
         st = trc_nao(ctx, nao_c)
         allocate (dens_buf(int(nao_c)*int(nao_c)*cached_nspin))
         st = trc_density(ctx, dens_buf)
         if (st /= TRC_OK) then
            call fail_eval(nvar, energy, gradient, status, "could not read the previous step's density for the warm start")
            return
         end if
      end if

      st = trc_set_geometry(ctx, reshape(real(xyz, c_double), [3*natm]))
      if (st /= TRC_OK) then
         call fail_eval(nvar, energy, gradient, status, "trc_set_geometry failed mid-optimisation")
         return
      end if

      if (have_prev_scf .and. cur_opts%guess_kind /= TRC_OPT_GUESS_NONE) then
         if (cur_opts%guess_kind == TRC_OPT_GUESS_DENSITY) then
            st = trc_set_guess(ctx, TRC_GUESS_GIVEN, c_loc(dens_buf), cached_nspin)
         else
            st = trc_set_guess(ctx, TRC_GUESS_PROJECTED, c_null_ptr, cached_nspin)
         end if
         if (st /= TRC_OK) then
            call fail_eval(nvar, energy, gradient, status, "trc_set_guess (warm start) failed mid-optimisation")
            return
         end if
      else if (have_prev_scf) then
         ! TRC_OPT_GUESS_NONE re-runs the ordinary guess every step with
         ! no density, orbitals or DIIS history carried over -- the
         ! baseline the warm-start measurement (--opt-guess density vs
         ! projected) compares against.
         st = trc_set_guess(ctx, cur_opts%fallback_guess_kind, c_null_ptr, 1_c_int)
         if (st /= TRC_OK) then
            call fail_eval(nvar, energy, gradient, status, "trc_set_guess (baseline) failed mid-optimisation")
            return
         end if
      end if

      st = trc_run_scf(ctx)
      converged_now = .false.
      if (st == TRC_OK .or. st == TRC_ERR_NOCONV) then
         st = trc_converged(ctx, conv)
         converged_now = (conv /= 0)
      end if

      if (.not. converged_now .and. have_prev_scf) then
         ! One retry from the molecule's ordinary (non-warm-started) guess
         ! before giving up -- the warm start is a speed optimisation, not
         ! something this step's success should depend on.
         call global_logger%warning("--opt: the warm-started SCF did not converge at step "// &
                                    to_string(n_steps + 1)//"; retrying from the standard guess")
         st = trc_set_guess(ctx, cur_opts%fallback_guess_kind, c_null_ptr, 1_c_int)
         if (st == TRC_OK) then
            st = trc_run_scf(ctx)
            if (st == TRC_OK .or. st == TRC_ERR_NOCONV) then
               st = trc_converged(ctx, conv)
               converged_now = (conv /= 0)
            end if
         end if
      end if

      if (.not. converged_now) then
         call fail_eval(nvar, energy, gradient, status, "the SCF did not converge at this geometry")
         return
      end if
      ! From here on there IS a converged SCF on `ctx` to warm-start the
      ! NEXT step's guess from (density or projected orbitals, per
      ! cur_opts%guess_kind) -- set once true and left there for the rest
      ! of this optimisation.
      have_prev_scf = .true.

      st = trc_run_gradient(ctx)
      if (st /= TRC_OK) then
         st = trc_message(ctx, msg_buf, int(size(msg_buf), c_int))
         call fail_eval(nvar, energy, gradient, status, "trc_run_gradient: "//f_string(msg_buf))
         return
      end if

      allocate (gbuf(nvar))
      st = trc_gradient(ctx, gbuf)
      if (st /= TRC_OK) then
         call fail_eval(nvar, energy, gradient, status, "trc_gradient failed after a successful trc_run_gradient")
         deallocate (gbuf)
         return
      end if
      st = trc_energy(ctx, e)
      st = trc_nspin(ctx, nspin_c)
      cached_nspin = nspin_c
      st = trc_iterations(ctx, niter)
      last_niter = int(niter)
      last_gmax = sqrt(maxval(sum(reshape(real(gbuf, dp), [3, natm])**2, dim=1)))
      last_grms = sqrt(sum(real(gbuf, dp)**2)/real(nvar, dp))

      energy = e
      gradient = gbuf
      status = 0_c_int
      deallocate (gbuf)
      if (allocated(dens_buf)) deallocate (dens_buf)

   end subroutine cb_get_gradient

   subroutine fail_eval(nvar, energy, gradient, status, why)
      !! One place to mark an evaluation failed and hand DL-FIND a zero
      !! energy/gradient and status = 1 -- DL-FIND's own contract for a
      !! failed evaluation. There is no recovery inside DL-FIND for this;
      !! trc_dlfind_optimize's caller reports the run failed once
      !! api_dl_find returns.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(out) :: energy
      real(c_double), intent(out) :: gradient(nvar)
      integer(c_int), intent(out) :: status
      character(len=*), intent(in) :: why
      eval_failed = .true.
      eval_fail_msg = why
      call global_logger%error("--opt: "//trim(why))
      energy = 0.0_c_double
      gradient = 0.0_c_double
      status = 1_c_int
   end subroutine fail_eval

   subroutine cb_put_coords(nvar, switch, energy, coords, iam) bind(c)
      !! DL-FIND reporting an ACCEPTED geometry -- switch 1 (an ordinary
      !! step) here; this bridge offers no P-RFO/TS mode, so switch 3
      !! (the transition structure) is accepted as an ordinary step too
      !! rather than ignored, and nothing else (switch 2/4/5: a
      !! transition mode or a `connect` run's downhill minima) is ever
      !! sent, since neither feature is wired up on the DL-FIND side
      !! (iopt is always DLF_OPT_LBFGS).
      integer(c_int), intent(in), value :: nvar
      integer(c_int), intent(in), value :: switch
      real(c_double), intent(in), value :: energy
      real(c_double), intent(in) :: coords(nvar)
      integer(c_int), intent(in), value :: iam   !! task id; only task 0 ever reports (ntasks = 1 here anyway)

      real(dp) :: wall

      if (iam /= 0_c_int) return
      if (switch /= DLF_PUT_STEP .and. switch /= DLF_PUT_TS) return
      if (nvar /= 3*natm) return

      latest_coords = reshape(real(coords, dp), [3, natm])
      latest_energy = real(energy, dp)
      have_result = .true.
      n_steps = n_steps + 1
      total_scf_iters = total_scf_iters + last_niter

      wall = 0.0_dp
      if (step_clock_running) then
         call step_clock%stop()
         wall = step_clock%get_elapsed_time()
      end if
      block
         character(len=1024) :: line
         write (line, "(a,i5,a,f20.12,a,es11.4,a,es11.4,a,i5,a,f10.2,a)") &
            "   opt step ", n_steps, "   E = ", latest_energy, "   max|g| = ", last_gmax, &
            "   rms|g| = ", last_grms, "   SCF iters ", last_niter, "   t = ", wall, " s"
         call global_logger%info(trim(line))
      end block
      call step_clock%start()
      step_clock_running = .true.

      if (traj_n < traj_cap) then
         traj_n = traj_n + 1
         traj_e(traj_n) = latest_energy
         traj_xyz(:, :, traj_n) = latest_coords
      end if

   end subroutine cb_put_coords

   subroutine cb_get_hessian(nvar, coords, hessian, status) bind(c)
      !! Second derivatives for P-RFO/NR -- always declined. A non-zero
      !! status is a SUPPORTED answer: DL-FIND falls back to a two-point
      !! finite difference of the gradients it is already requesting.
      !! Never reached anyway, since iopt is always
      !! DLF_OPT_LBFGS, but every one of the seven callbacks must exist.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: hessian(nvar, nvar)
      integer(c_int), intent(out) :: status
      hessian = 0.0_c_double
      status = 1_c_int
   end subroutine cb_get_hessian

   subroutine cb_get_multistate_gradients(nvar, coords, energy, gradient, coupling, &
                                          needcoupling, iimage, status) bind(c)
      !! Conical-intersection search -- not offered.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: energy(2)
      real(c_double), intent(out) :: gradient(nvar, 2)
      real(c_double), intent(out) :: coupling(nvar)
      integer(c_int), intent(in), value :: needcoupling, iimage
      integer(c_int), intent(out) :: status
      energy = 0.0_c_double
      gradient = 0.0_c_double
      coupling = 0.0_c_double
      status = 1_c_int
   end subroutine cb_get_multistate_gradients

   subroutine cb_error() bind(c)
      !! DL-FIND's last word before libdlfind's own api.f90 executes
      !! `error stop` -- there is no recovering from this, which is why
      !! this bridge is linked only into terco_cli, never into libterco.
      !! Say what was being optimised while there is still somewhere to
      !! say it.
      call global_logger%error("--opt: DL-FIND failed internally and is about to stop the process "// &
                               "(atoms: "//to_string(natm)//", accepted steps so far: "//to_string(n_steps)// &
                               ", energy evaluations: "//to_string(n_evals)//")")
   end subroutine cb_error

   subroutine cb_update() bind(c)
      !! Called between steps; nothing here needs to happen then.
   end subroutine cb_update

   function dlfind_print_level() result(printl)
      !! DL-FIND's own `printl`, from this process's logger level --
      !! DL-FIND writes straight to stdout and cannot be routed into
      !! pic_logger, so --silent/--verbose at least still mutes it.
      use pic_logger, only: debug_level, verbose_level
      integer(c_int) :: printl
      if (global_logger%log_level >= debug_level) then
         printl = 4_c_int
      else if (global_logger%log_level >= verbose_level) then
         printl = 2_c_int
      else
         printl = 0_c_int
      end if
   end function dlfind_print_level

end module trc_dlfind_bridge
