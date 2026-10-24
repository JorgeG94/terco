!
! terco's side of DL-FIND: fills libdlfind's `spec`/parameter block, turns
! its seven callbacks into calls on a terco context handle (trc_set_geometry
! / trc_run_scf / trc_run_gradient / trc_run_hessian / trc_hessian), and
! turns the result back into a trc_opt_result_t. The seven-callback
! contract, the `spec` layout and its traps below come from reading
! libdlfind's own api.f90, dlf_parameter_module.f90, dl-find.f90 and
! dlf_task.f90 directly; each one is explained inline where it applies
! rather than pointed at elsewhere.
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
! WHAT THIS BRIDGE EXPOSES (see app/main.F90's --opt-* options for the CLI
! names) and what it does NOT, with the reason:
!   - exposed: every single-state minimiser (sd, cg-auto, cg, lbfgs,
!     damped, nr) and P-RFO, both as a minimiser and for a saddle search;
!     all five icoord families (cartesian, hdlc, hdlc-tc, dlc, dlc-tc);
!     the dimer method for a saddle with no Hessian; the `connect` task
!     (1011: TS then both downhill minima); the external-Hessian callback
!     (cb_get_hessian, via trc_run_hessian/trc_hessian) and DL-FIND's own
!     Hessian update/zero-mode/soft-mode/trust/line-search/L-BFGS-memory
!     controls.
!   - NOT exposed, and why:
!       * imultistate (conical intersections, Lagrange-Newton, iopt 40) --
!         needs excited-state gradients/couplings terco has no reference
!         for; cb_get_multistate_gradients always declines (status 1).
!       * imicroiter (microiterative QM/MM) -- needs a second, cheaper MM
!         engine terco does not have.
!       * NEB/climbing image (icoord 100-199), frozen atoms, internal-
!         coordinate constraints, the parallel/genetic/stochastic global
!         searches (iopt 51/52), qTS/instanton (iopt 11/12/13) -- not yet
!         exposed, or out of scope for a single-state QC engine (global
!         search belongs to a different kind of surface).
!       * ntasks > 1 (parallel multi-image tasks), dump/restart
!         checkpointing, the `distort` override for `connect` (DL-FIND's
!         own default of 0.1 is used), a second coordinate set (nvar2;
!         only meaningful for NEB/dimer endpoints, and the dimer method
!         here uses DL-FIND's own random/internal initial direction, not
!         a supplied one) -- no CLI need for them yet; dimer's own
!         maxrot/tolrot (rotation count/tolerance) are left at DL-FIND's
!         defaults for the same reason (not asked for, and the dimer
!         rotation cost is already dominated by the gradient callback,
!         not the rotation count).
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
                       trc_set_guess, trc_density, trc_message, trc_run_hessian, trc_hessian, &
                       TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_GIVEN, TRC_GUESS_PROJECTED
   use terco_app, only: f_string, write_opt_xyz_frame, write_traj_xyz_frames, write_opt_failure_json
   implicit none
   private

   public :: trc_dlfind_available
   public :: TRC_OPT_COORDS_CARTESIAN, TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_HDLC_TC, &
             TRC_OPT_COORDS_DLC, TRC_OPT_COORDS_DLC_TC
   public :: TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, TRC_OPT_GUESS_NONE
   public :: TRC_OPT_ALG_SD, TRC_OPT_ALG_CG_AUTO, TRC_OPT_ALG_CG, TRC_OPT_ALG_LBFGS, &
             TRC_OPT_ALG_PRFO, TRC_OPT_ALG_NR, TRC_OPT_ALG_DAMPED
   public :: TRC_OPT_TARGET_MINIMUM, TRC_OPT_TARGET_SADDLE
   public :: TRC_OPT_SADDLE_PRFO, TRC_OPT_SADDLE_DIMER
   public :: TRC_OPT_HESS_TERCO, TRC_OPT_HESS_DLFIND
   public :: TRC_OPT_HUPD_NONE, TRC_OPT_HUPD_POWELL, TRC_OPT_HUPD_BOFILL, TRC_OPT_HUPD_BFGS
   public :: TRC_OPT_LINE_AUTO, TRC_OPT_LINE_SCALING, TRC_OPT_LINE_ENERGY_TRUST, TRC_OPT_LINE_GRADIENT_TRUST
   public :: trc_opt_options_t, trc_opt_result_t
   public :: trc_dlfind_optimize

   ! --opt-coords: terco's own enum (main.F90), equal by construction to
   ! DL-FIND's own icoord base (0-4) so cb_get_params needs no translation
   ! table -- only the dimer/NEB offsets (+100, +220, ...) are added on
   ! top, below.
   integer, parameter :: TRC_OPT_COORDS_CARTESIAN = 0, TRC_OPT_COORDS_HDLC = 1, TRC_OPT_COORDS_HDLC_TC = 2, &
                         TRC_OPT_COORDS_DLC = 3, TRC_OPT_COORDS_DLC_TC = 4

   integer, parameter :: TRC_OPT_GUESS_DENSITY = 0, TRC_OPT_GUESS_PROJECTED = 1, TRC_OPT_GUESS_NONE = 2
      !! --opt-guess: the previous step's converged density (TRC_GUESS_GIVEN,
      !! fetched with trc_density before the geometry moves) or its
      !! orbitals re-orthonormalised in the new overlap (TRC_GUESS_PROJECTED,
      !! needs no density fetch at all, see scf_at_geometry's own comment).
      !! TRC_OPT_GUESS_NONE re-runs `fallback_guess_kind` (ordinarily
      !! SAD/SADQ) at every step, the mqc baseline the warm-start
      !! measurement compares the other two against; not reachable from
      !! --opt-guess (only "density" and "projected" are), used only by
      !! that measurement itself.

   ! --opt-algorithm: equal to DL-FIND's own iopt (dl-find.f90's own
   ! inline table; iopt==9/51/52/40 are never offered -- see this
   ! module's header comment's "NOT exposed" list).
   integer, parameter :: TRC_OPT_ALG_SD = 0, TRC_OPT_ALG_CG_AUTO = 1, TRC_OPT_ALG_CG = 2, TRC_OPT_ALG_LBFGS = 3, &
                         TRC_OPT_ALG_PRFO = 10, TRC_OPT_ALG_NR = 20, TRC_OPT_ALG_DAMPED = 30

   integer, parameter :: TRC_OPT_TARGET_MINIMUM = 0, TRC_OPT_TARGET_SADDLE = 1
   integer, parameter :: TRC_OPT_SADDLE_PRFO = 0, TRC_OPT_SADDLE_DIMER = 1
   integer, parameter :: TRC_OPT_HESS_TERCO = 0, TRC_OPT_HESS_DLFIND = 1
      !! --opt-init-hessian: terco's own analytic-from-FD Cartesian
      !! Hessian through cb_get_hessian (inithessian 0), or DL-FIND's own
      !! two-point finite difference of the gradients it is already
      !! requesting (inithessian 2; DL-FIND's own fallback whenever
      !! cb_get_hessian declines, so "dlfind" here and a decline are the
      !! same inithessian value -- the CLI value is honoured by never
      !! calling trc_run_hessian at all when it is set, not just by
      !! declining after the fact).
   integer, parameter :: TRC_OPT_HUPD_NONE = 0, TRC_OPT_HUPD_POWELL = 1, TRC_OPT_HUPD_BOFILL = 2, &
                         TRC_OPT_HUPD_BFGS = 3
      !! --opt-hessian-update: equal to DL-FIND's own `update`.
   integer, parameter :: TRC_OPT_LINE_AUTO = -1, TRC_OPT_LINE_SCALING = 0, TRC_OPT_LINE_ENERGY_TRUST = 1, &
                         TRC_OPT_LINE_GRADIENT_TRUST = 2
      !! --opt-linesearch: equal to DL-FIND's own `iline` (AUTO leaves
      !! DL-FIND's own default, which depends on iopt).

   ! DL-FIND's own numbering (api_dl_find's icoord/iopt/task/switch) that
   ! this module itself adds on top of the CLI-facing enums above -- see
   ! libdlfind's src/dl-find.f90 and src/dlf_task.f90.
   integer(c_int), parameter :: DLF_COORDS_DIMER_BASE = 220_c_int
      !! Dimer method, rotation by line search with extrapolation (the
      !! family DL-FIND's own dlf_task.f90 picks for an unspecified TS
      !! search): icoord = 220 + the inner coordinate system's own base
      !! (0 cartesian ... 4 dlc-tc).
   integer(c_int), parameter :: DLF_INITHESS_EXTERNAL = 0, DLF_INITHESS_TWO_POINT = 2
   integer(c_int), parameter :: DLF_PUT_STEP = 1, DLF_PUT_TS = 3, DLF_PUT_MIN_PLUS = 4, DLF_PUT_MIN_MINUS = 5
      !! cb_put_coords's own `switch` (dlf_task.f90's task-1011 case):
      !! 1 an ordinary accepted step; 3 the transition structure (also
      !! the final answer of a plain, non-connect P-RFO/dimer search); 4
      !! and 5 the two downhill minima a `connect` run adds afterwards.
   integer(c_int), parameter :: DLF_TASK_CONNECT = 1011_c_int
   integer, parameter :: MSG_LEN = 1024
      !! trc_message's own buffer length convention (app/main.F90's `api`)
   integer, parameter :: ACCEPT_ALWAYS = 0, ACCEPT_ENERGY = 1
      !! This bridge's own replica of DL-FIND's default `iline` rule --
      !! see accept_rule's own comment, below.

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

      ! --- algorithm, saddle-search target and Hessian handling ---------
      integer :: algorithm = TRC_OPT_ALG_LBFGS
      integer :: target_kind = TRC_OPT_TARGET_MINIMUM
      integer :: saddle_method = TRC_OPT_SADDLE_PRFO
      logical :: connect = .false.
      logical :: is_linear = .false.
         !! Geometric test on the START geometry (main.F90's own
         !! collinearity check -- no mass or Hessian needed yet); used
         !! only to pick zero_modes' own default below.
      integer :: hessian_source = TRC_OPT_HESS_TERCO
      integer :: hessian_update = -1    !! sentinel: leave DL-FIND's own default (Bofill)
      integer :: max_updates = -1       !! sentinel: leave DL-FIND's own default (50)
      integer :: zero_modes = -1        !! sentinel: see cb_get_params' own auto-default
      real(dp) :: soft_modes = -1.0_dp  !! sentinel: leave DL-FIND's own default (disabled)
      real(dp) :: max_step = -1.0_dp    !! sentinel: leave DL-FIND's own default (0.5 Bohr)
      integer :: linesearch = TRC_OPT_LINE_AUTO
      real(dp) :: scalestep = -1.0_dp   !! sentinel: leave DL-FIND's own default (1.0)
      integer :: lbfgs_memory = -1      !! sentinel: leave DL-FIND's own default
      real(dp) :: etol = -1.0_dp        !! sentinel: leave DL-FIND's own default (tolerance/450)
      integer :: print_level = -1       !! sentinel: follow this process's own logger level
      character(len=:), allocatable :: opt_xyz_path   !! unallocated/"" = do not write
      character(len=:), allocatable :: traj_xyz_path  !! unallocated/"" = do not write
      character(len=:), allocatable :: json_path      !! unallocated/"" = no --json requested
      character(len=:), allocatable :: geometry_path  !! for the failure JSON's "geometry" field
      character(len=:), allocatable :: coords_name    !! "cartesian"/"hdlc"/.../"dlc-tc", for the failure JSON
      character(len=:), allocatable :: guess_name     !! "density"/"projected", for the failure JSON
   end type trc_opt_options_t

   type :: trc_opt_result_t
      logical :: converged = .false.
      logical :: failed = .false.
      character(len=:), allocatable :: message
      integer :: steps = 0               !! ACCEPTED geometries (cb_put_coords, this module's own
                                          !! replica of DL-FIND's test_acceptance/the always-accept
                                          !! rule for the algorithm in use -- see cb_put_coords)
      integer :: energy_evaluations = 0  !! every cb_get_gradient call, rejected trials included
      integer :: hessian_evaluations = 0  !! every cb_get_hessian call that actually ran trc_run_hessian
      real(dp) :: final_energy = 0.0_dp
      real(dp) :: final_gmax = 0.0_dp
      real(dp) :: final_grms = 0.0_dp
      real(dp), allocatable :: final_coords(:, :)     !! (3, natm), Bohr -- the last ACCEPTED geometry;
                                                       !! a connect run's TS when one was found (below)
      integer :: n_traj = 0
      real(dp), allocatable :: traj_energy(:)
      real(dp), allocatable :: traj_coords(:, :, :)   !! (3, natm, n_traj), Bohr -- ACCEPTED steps only
      integer :: total_scf_iterations = 0
         !! Sum of trc_iterations over EVERY evaluation, rejected trials
         !! included (cb_get_gradient's own `st = trc_iterations` runs for
         !! each one) -- the single number the warm-start measurement
         !! (--opt-guess density vs projected vs the mqc baseline,
         !! trc_opt_options_t's own TRC_OPT_GUESS_NONE) compares, since a
         !! rejected trial's SCF iterations are real work too.
      ! --opt-connect (task 1011): the transition state cb_put_coords's
      ! switch 3 reported, and the two downhill minima switches 4/5
      ! reported -- ALL THREE in addition to final_energy/final_coords
      ! above, which for a connect run stay the transition structure
      ! (the primary object of the search; the minima are confirmation).
      ! Each carries its OWN convergence verdict and gradient, snapshotted
      ! when that sub-search's own final point arrived -- never the other
      ! two sub-searches' (see cb_put_coords).
      logical :: have_ts = .false.
      real(dp) :: ts_energy = 0.0_dp
      real(dp) :: ts_gmax = 0.0_dp, ts_grms = 0.0_dp
      logical :: ts_converged = .false.
      real(dp), allocatable :: ts_coords(:, :)
      logical :: have_minimum_plus = .false., have_minimum_minus = .false.
      real(dp) :: minimum_plus_energy = 0.0_dp, minimum_minus_energy = 0.0_dp
      real(dp) :: minimum_plus_gmax = 0.0_dp, minimum_plus_grms = 0.0_dp
      real(dp) :: minimum_minus_gmax = 0.0_dp, minimum_minus_grms = 0.0_dp
      logical :: minimum_plus_converged = .false., minimum_minus_converged = .false.
      real(dp), allocatable :: minimum_plus_coords(:, :), minimum_minus_coords(:, :)
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
   integer, save :: n_evals = 0, n_steps = 0, n_hess_evals = 0, total_scf_iters = 0
   logical, save :: have_result = .false.
   real(dp), save :: latest_energy = 0.0_dp
   real(dp), allocatable, save :: latest_coords(:, :)
   real(dp), save :: last_gmax = 0.0_dp, last_grms = 0.0_dp
      !! THIS evaluation's own max/rms |g| (cb_get_gradient), whether or
      !! not it turns out accepted -- cb_put_coords reads these to decide
      !! what to snapshot into accepted_gmax/accepted_grms below. Never
      !! read directly by trc_dlfind_optimize's own result assembly: a run
      !! that stops right after a REJECTED trial would otherwise report
      !! that trial's gradient against an ACCEPTED geometry.
   integer, save :: last_niter = 0
   logical, save :: have_accepted_energy = .false.
   real(dp), save :: last_accepted_energy = 0.0_dp
      !! DL-FIND's own `glob%oldenergy`: an evaluation in an energy-based
      !! phase (accept_rule == ACCEPT_ENERGY, below) is accepted iff it is
      !! the first evaluation of its phase, or its energy is lower than
      !! this. Reset to "no accepted energy yet" at the start of each
      !! phase a `connect` run adds (the TS search, then each downhill
      !! minimisation DL-FIND resets its own trust radius/oldenergy for --
      !! see cb_put_coords' own comment on switches 3 and 4).
   real(dp), save :: accepted_gmax = 0.0_dp, accepted_grms = 0.0_dp
      !! The max/rms |g| of the last ACCEPTED evaluation IN THE CURRENT
      !! PHASE -- what trc_dlfind_optimize reports as final_gmax/final_grms
      !! (or, for a connect run, ts_gmax/minimum_plus_gmax/minimum_minus_gmax)
      !! and what the CLI's convergence verdict is tested against, never
      !! last_gmax/last_grms above.
   integer, save :: accept_rule = ACCEPT_ALWAYS
      !! This run's own replica of DL-FIND's default `iline` rule
      !! (dlf_default_set's inline table in dl-find.f90, reproduced in
      !! cb_get_params below): ACCEPT_ENERGY for L-BFGS on a non-dimer
      !! coordinate system (DL-FIND's iline=1, test_acceptance's energy
      !! comparison); ACCEPT_ALWAYS for P-RFO, NR, damped or the dimer
      !! method (iline=0 -- these never reject a trial at this level; a
      !! trust radius internal to dlf_formstep governs the step instead,
      !! which test_acceptance/test_acceptance_g never run for them) AND,
      !! approximated the same way for lack of access to DL-FIND's own
      !! internal-coordinate step/gradient history, for sd/cg/cg-auto
      !! (iline=2, test_acceptance_g's gradient-projection rule -- this
      !! bridge cannot replicate it from a Cartesian gradient alone, so it
      !! never rejects a trial for these either; not a regression, since
      !! no earlier version of this bridge modelled rejection for them).
      !! A connect run's downhill phases force ACCEPT_ENERGY regardless
      !! (dlf_task.f90's task==1011 case sets iopt=3/iline=1 for both
      !! unconditionally, no matter what the TS search itself used).
   type(timer_type), save :: step_clock
   logical, save :: step_clock_running = .false.
   integer, save :: traj_cap = 0, traj_n = 0
   real(dp), allocatable, save :: traj_e(:)
   real(dp), allocatable, save :: traj_xyz(:, :, :)
   real(dp), allocatable, save :: last_scf_xyz(:, :)
      !! The geometry scf_at_geometry last ran AND CONVERGED at -- lets
      !! cb_get_hessian skip re-running the SCF when DL-FIND asks for a
      !! Hessian at the point cb_get_gradient just converged (the common
      !! case for P-RFO/NR: the Hessian callback runs right after the step
      !! that needs it), rather than repeat the identical SCF.
   ! --opt-connect's own three structures (cb_put_coords's own switches
   ! 3/4/5, dlf_task.f90's task==1011 case): snapshotted as they arrive,
   ! independent of the ordinary step/trajectory bookkeeping above, which
   ! keeps running underneath each of the connect run's three sub-optimisations.
   logical, save :: have_ts = .false., have_min_plus = .false., have_min_minus = .false.
   real(dp), save :: ts_energy = 0.0_dp, min_plus_energy = 0.0_dp, min_minus_energy = 0.0_dp
   real(dp), save :: ts_gmax = 0.0_dp, ts_grms = 0.0_dp
   real(dp), save :: min_plus_gmax = 0.0_dp, min_plus_grms = 0.0_dp
   real(dp), save :: min_minus_gmax = 0.0_dp, min_minus_grms = 0.0_dp
   real(dp), allocatable, save :: ts_coords(:, :), min_plus_coords(:, :), min_minus_coords(:, :)

contains

   pure function trc_dlfind_available() result(ok)
      !! Whether THIS build can optimise a geometry -- always true here;
      !! app/trc_dlfind_stub.F90's own answer (the TERCO_ENABLE_DLFIND=OFF
      !! build) is always false, which is all app/main.F90 ever asks.
      logical :: ok
      ok = .true.
   end function trc_dlfind_available

   subroutine trc_dlfind_optimize(handle, n, z, residues_in, coords, opts, res)
      !! Minimise (or find a saddle on) the molecule already set up on
      !! `handle` (molecule, basis and method must already match
      !! `coords`/`z` -- this routine only ever moves atoms on it,
      !! through trc_set_geometry, never builds anything). The optimised
      !! geometry, energy and trajectory come back in `res`; `coords`
      !! itself is read-only (the starting point cb_get_params hands
      !! DL-FIND).
      type(c_ptr), intent(in) :: handle
      integer, intent(in) :: n
      integer, intent(in) :: z(n)
      integer, intent(in) :: residues_in(n)
         !! 1-based residue number per atom (connected components for
         !! --opt-coords hdlc/hdlc-tc, all 1 for cartesian/dlc/dlc-tc --
         !! app/main.F90's own job to build, trc_constants%trc_covalent_radius-based).
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
      n_hess_evals = 0
      total_scf_iters = 0
      have_result = .false.
      latest_energy = 0.0_dp
      if (allocated(latest_coords)) deallocate (latest_coords)
      allocate (latest_coords(3, n), source=coords)
      last_gmax = 0.0_dp
      last_grms = 0.0_dp
      last_niter = 0
      have_accepted_energy = .false.
      last_accepted_energy = 0.0_dp
      accepted_gmax = 0.0_dp
      accepted_grms = 0.0_dp
      accept_rule = ACCEPT_ALWAYS
      step_clock_running = .false.
      if (allocated(last_scf_xyz)) deallocate (last_scf_xyz)
      have_ts = .false.
      have_min_plus = .false.
      have_min_minus = .false.
      ts_energy = 0.0_dp
      min_plus_energy = 0.0_dp
      min_minus_energy = 0.0_dp
      ts_gmax = 0.0_dp
      ts_grms = 0.0_dp
      min_plus_gmax = 0.0_dp
      min_plus_grms = 0.0_dp
      min_minus_gmax = 0.0_dp
      min_minus_grms = 0.0_dp
      if (allocated(ts_coords)) deallocate (ts_coords)
      if (allocated(min_plus_coords)) deallocate (min_plus_coords)
      if (allocated(min_minus_coords)) deallocate (min_minus_coords)

      ! A STARTING capacity, not a hard upper bound: HDLC/DLC's own
      ! `trerun_energy` evaluations (dl-find.f90's own COORDS TRANSFORM
      ! section) don't advance `stat%ccycle`, so more than maxiter
      ! ACCEPTED steps per sub-optimisation is possible even though
      ! DL-FIND enforces maxcycle on cycles, not accepted steps as such --
      ! grow_traj (record_step) doubles these with move_alloc rather than
      ! drop a frame past this.
      traj_cap = max(opts%maxiter*merge(3, 1, opts%connect) + 8, 8)
      if (allocated(traj_e)) deallocate (traj_e)
      if (allocated(traj_xyz)) deallocate (traj_xyz)
      allocate (traj_e(traj_cap), traj_xyz(3, n, traj_cap))
      traj_n = 0

      nvar = int(3*n, c_int)
      nvar2 = 0_c_int   ! no second structure (no NEB endpoint, no supplied dimer direction)
      nspec = int(3*n, c_int)   ! [residue(n), znuc(n), micspec(n)]: no constraints, so no icons block

      call api_dl_find(nvar, nvar2, nspec, 1_c_int, &
                       c_funloc(cb_error), c_funloc(cb_get_gradient), c_funloc(cb_get_hessian), &
                       c_funloc(cb_get_multistate_gradients), c_funloc(cb_get_params), &
                       c_funloc(cb_put_coords), c_funloc(cb_update))

      res%steps = n_steps
      res%energy_evaluations = n_evals
      res%hessian_evaluations = n_hess_evals
      res%total_scf_iterations = total_scf_iters
      if (eval_failed) then
         res%failed = .true.
         res%message = "geometry optimisation stopped: "//trim(eval_fail_msg)
      else if (.not. have_result) then
         ! DL-FIND reports a geometry through dlf_put_coords only when
         ! printf >= 1 -- leaving it at DL-FIND's own default of 0 lets a
         ! run converge with no geometry ever reported back, which is why
         ! cb_get_params always sets it to at least 1 below. Reaching here
         ! with have_result still false means api_dl_find returned
         ! without ever calling it.
         res%failed = .true.
         res%message = "DL-FIND returned without reporting a geometry"
      else
         ! latest_coords/latest_energy and accepted_gmax/accepted_grms are
         ! the last ACCEPTED evaluation's own (cb_put_coords), never
         ! whatever cb_get_gradient last happened to compute -- a run that
         ! stops right after a rejected trial must still report the
         ! accepted point.
         res%final_energy = latest_energy
         res%final_coords = latest_coords
         res%final_gmax = accepted_gmax
         res%final_grms = accepted_grms
         res%converged = accepted_gmax <= opts%gtol
      end if
      res%n_traj = traj_n
      if (traj_n > 0) then
         res%traj_energy = traj_e(1:traj_n)
         res%traj_coords = traj_xyz(:, :, 1:traj_n)
      end if

      res%have_ts = have_ts
      if (have_ts) then
         res%ts_energy = ts_energy
         res%ts_coords = ts_coords
         res%ts_gmax = ts_gmax
         res%ts_grms = ts_grms
         res%ts_converged = ts_gmax <= opts%gtol
         ! A connect run's primary answer is the TS, not whichever
         ! sub-optimisation happened to run last (the second downhill
         ! minimisation) -- override final_* with the TS SEARCH's own
         ! point AND its own convergence verdict/gradient, not the last
         ! downhill minimisation's (what final_gmax/final_grms/converged
         ! were set to above): a TS sub-search that never reached
         ! --opt-gtol must not report "converged" just because the two
         ! downhill minimisations after it did.
         if (opts%connect .and. .not. res%failed) then
            res%final_energy = ts_energy
            res%final_coords = ts_coords
            res%final_gmax = ts_gmax
            res%final_grms = ts_grms
            res%converged = res%ts_converged
         end if
      end if
      res%have_minimum_plus = have_min_plus
      if (have_min_plus) then
         res%minimum_plus_energy = min_plus_energy
         res%minimum_plus_coords = min_plus_coords
         res%minimum_plus_gmax = min_plus_gmax
         res%minimum_plus_grms = min_plus_grms
         res%minimum_plus_converged = min_plus_gmax <= opts%gtol
      end if
      res%have_minimum_minus = have_min_minus
      if (have_min_minus) then
         res%minimum_minus_energy = min_minus_energy
         res%minimum_minus_coords = min_minus_coords
         res%minimum_minus_gmax = min_minus_gmax
         res%minimum_minus_grms = min_minus_grms
         res%minimum_minus_converged = min_minus_gmax <= opts%gtol
      end if

      active = .false.
      ctx = c_null_ptr
      if (allocated(znuc)) deallocate (znuc)
      if (allocated(residues)) deallocate (residues)
      if (allocated(start_coords)) deallocate (start_coords)
      if (allocated(latest_coords)) deallocate (latest_coords)
      if (allocated(traj_e)) deallocate (traj_e)
      if (allocated(traj_xyz)) deallocate (traj_xyz)
      if (allocated(ts_coords)) deallocate (ts_coords)
      if (allocated(min_plus_coords)) deallocate (min_plus_coords)
      if (allocated(min_minus_coords)) deallocate (min_minus_coords)
      if (allocated(last_scf_xyz)) deallocate (last_scf_xyz)
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
      !! default (dlf_default_init/dlf_default_set ran first), so only
      !! what this bridge actually sets below is touched -- see this
      !! module's own header comment for what is left untouched on
      !! purpose.
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
      integer(c_int) :: coord_base
      logical :: is_saddle, is_dimer, is_prfo_or_nr

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
      spec(2*natm + 1:3*natm) = 1_c_int   ! micspec: DL-FIND requires 0 or 1 here; 1 everywhere

      nz = int(natm, c_int)       ! lets HDLC/DLC perceive connectivity on their own
      ncons = 0_c_int
      nconn = 0_c_int             ! DL-FIND reads nconn and micspec at inconsistent offsets when both are present; 0 avoids that
      tatoms = 1_c_int            ! the variables are atoms, not arbitrary DOF -- HDLC needs to know

      ! Arrive uninitialised, or at -1 and then checked positive, so each
      ! needs an explicit value even though it looks optional.
      nframe = 0_c_int
      nmass = 0_c_int
      nweight = 0_c_int
      ntasks = 1_c_int

      select case (cur_opts%coords_kind)
      case (TRC_OPT_COORDS_HDLC)
         coord_base = 1_c_int
      case (TRC_OPT_COORDS_HDLC_TC)
         coord_base = 2_c_int
      case (TRC_OPT_COORDS_DLC)
         coord_base = 3_c_int
      case (TRC_OPT_COORDS_DLC_TC)
         coord_base = 4_c_int
      case default
         coord_base = 0_c_int
      end select

      is_saddle = (cur_opts%target_kind == TRC_OPT_TARGET_SADDLE)
      is_dimer = (is_saddle .and. cur_opts%saddle_method == TRC_OPT_SADDLE_DIMER)

      if (is_dimer) then
         ! Dimer method, rotation by line search with extrapolation --
         ! DL-FIND's own choice when a TS search does not pick P-RFO
         ! (dlf_task.f90's task==1011 case; picked directly here too, so
         ! a standalone --opt-saddle-method dimer run -- no --opt-connect
         ! -- behaves the same way). No Hessian is used; iopt is fixed
         ! (the dimer's translation step is DL-FIND's own L-BFGS, not a
         ! choice this bridge exposes -- --opt-algorithm is refused
         ! together with --opt-target saddle in app/main.F90 for exactly
         ! this reason).
         icoord = DLF_COORDS_DIMER_BASE + coord_base
         iopt = TRC_OPT_ALG_LBFGS
      else if (is_saddle) then
         icoord = coord_base
         iopt = TRC_OPT_ALG_PRFO
      else
         icoord = coord_base
         iopt = int(cur_opts%algorithm, c_int)
      end if

      is_prfo_or_nr = (iopt == TRC_OPT_ALG_PRFO .or. iopt == TRC_OPT_ALG_NR)

      if (is_prfo_or_nr) then
         ! cb_get_hessian below actually supplies one only when this is
         ! TRC_OPT_HESS_TERCO; "dlfind" here means never calling it in
         ! the first place, not just declining once asked (same effect
         ! on DL-FIND's side either way: inithessian = 2, its own
         ! two-point FD of the gradients it already requests).
         if (cur_opts%hessian_source == TRC_OPT_HESS_DLFIND) then
            inithessian = DLF_INITHESS_TWO_POINT
         else
            inithessian = DLF_INITHESS_EXTERNAL
         end if
         if (cur_opts%hessian_update >= 0) update = int(cur_opts%hessian_update, c_int)
         if (cur_opts%max_updates >= 0) maxupd = int(cur_opts%max_updates, c_int)
      else
         inithessian = DLF_INITHESS_TWO_POINT   ! never consulted by a Hessian-free algorithm
      end if

      ! --opt-zero-modes' own auto default: ANY Cartesian search that
      ! follows a Hessian eigenvector -- P-RFO (minimum or saddle), NR, or
      ! the dimer method -- chases the LOWEST-curvature direction and,
      ! left with translation/rotation still mixed in among the true
      ! vibrational modes, chases one of those instead. Measured for
      ! Cartesian P-RFO on HCN<->HNC: 60 steps with no convergence vs 13
      ! with nzero=6; separately measured for NR-as-a-minimiser on water,
      ! which does not reach --opt-gtol in 100 steps at
      ! nzero=0 but does with the same default; the dimer method's own
      ! rotation step is exposed to the same failure. Internal
      ! coordinates already remove translation/rotation, so they keep
      ! DL-FIND's own already-resolved default (0) unless the CLI gave
      ! one.
      if (cur_opts%zero_modes >= 0) then
         nzero = int(cur_opts%zero_modes, c_int)
      else if ((iopt == TRC_OPT_ALG_PRFO .or. iopt == TRC_OPT_ALG_NR .or. is_dimer) .and. &
              cur_opts%coords_kind == TRC_OPT_COORDS_CARTESIAN) then
         nzero = merge(5_c_int, 6_c_int, cur_opts%is_linear)
      end if

      if (cur_opts%soft_modes > 0.0_dp) soft = real(cur_opts%soft_modes, c_double)
      if (cur_opts%max_step > 0.0_dp) maxstep = real(cur_opts%max_step, c_double)
      if (cur_opts%linesearch /= TRC_OPT_LINE_AUTO) iline = int(cur_opts%linesearch, c_int)

      ! This bridge's own replica of dlf_default_set's defaulting rule for
      ! `iline` (dl-find.f90): cb_put_coords has no other way to learn
      ! which acceptance test DL-FIND is actually using for this run's
      ! algorithm, since `iline` itself is left at DL-FIND's own sentinel
      ! here unless --opt-linesearch overrode it just above.
      if (iline >= 0_c_int) then
         accept_rule = merge(ACCEPT_ENERGY, ACCEPT_ALWAYS, iline == 1_c_int)
      else if (iopt < 3_c_int) then
         ! sd/cg/cg-auto: DL-FIND's own default is iline=2
         ! (test_acceptance_g's gradient-projection rule) -- approximated
         ! as ACCEPT_ALWAYS, see accept_rule's own comment.
         accept_rule = ACCEPT_ALWAYS
      else if (iopt == TRC_OPT_ALG_LBFGS .and. icoord < 10_c_int) then
         ! L-BFGS on cartesian/hdlc/hdlc-tc/dlc/dlc-tc: DL-FIND's own
         ! default is iline=1 (test_acceptance's energy comparison).
         accept_rule = ACCEPT_ENERGY
      else
         ! P-RFO, NR, damped, or L-BFGS+dimer: DL-FIND's own default is
         ! iline=0 -- no energy-based rejection at this level; a trust
         ! radius internal to dlf_formstep governs the step instead.
         accept_rule = ACCEPT_ALWAYS
      end if

      if (cur_opts%scalestep > 0.0_dp) scalestep = real(cur_opts%scalestep, c_double)
      if (cur_opts%lbfgs_memory > 0) lbfgs_mem = int(cur_opts%lbfgs_memory, c_int)
      if (cur_opts%etol > 0.0_dp) tolerance_e = real(cur_opts%etol, c_double)

      tolerance = real(cur_opts%gtol, c_double)
      maxcycle = int(cur_opts%maxiter, c_int)

      if (cur_opts%connect) then
         task = DLF_TASK_CONNECT
         ! dlf_task.f90's task==1011 case writes TS.xyz/minimum_+-.xyz and
         ! calls dlf_put_coords for switches 3/4/5 ONLY under printf>=4.
         printf = 4_c_int
      else
         ! printf >= 1 is what makes DL-FIND call dlf_put_coords at all
         ! for an ordinary run -- left at DL-FIND's own default of 0, a
         ! run converges silently with no geometry ever reported back.
         printf = 1_c_int
      end if

      if (cur_opts%print_level >= 0) then
         printl = int(cur_opts%print_level, c_int)
      else
         ! printl follows this process's own logger level, so --silent
         ! also quiets DL-FIND's stdout convergence table.
         printl = dlfind_print_level()
      end if

      dump = 0_c_int      ! no checkpoint files -- two optimisations sharing a CWD must not read each other's
      restart = 0_c_int

   end subroutine cb_get_params

   !
   ! Converge the SCF at `xyz`, warm-started from the previous accepted
   ! point if there is one, same contract for every caller that needs a
   ! converged SCF on `ctx` at a given geometry (cb_get_gradient AND
   ! cb_get_hessian -- factored out so the warm-start logic exists once).
   ! Does NOT touch module `save` step/trajectory bookkeeping (last_niter,
   ! n_steps, ...): that stays the gradient callback's own job, since only
   ! an ACCEPTED gradient evaluation's iteration count is what
   ! cb_put_coords later adds into total_scf_iterations.
   !
   subroutine scf_at_geometry(xyz, niter, ok, why)
      real(dp), intent(in) :: xyz(3, natm)
      integer, intent(out) :: niter
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why

      real(c_double), allocatable, target :: dens_buf(:)
      integer(c_int) :: st, conv, niter_c, nao_c, nspin_c
      logical :: converged_now

      ok = .false.
      why = ""
      niter = 0

      ! The warm-start density has to be read BEFORE the geometry moves:
      ! trc_set_geometry invalidates the context's scf_ok, and trc_density
      ! refuses once that happens. TRC_GUESS_PROJECTED needs no such fetch
      ! -- it reads the context's own last-converged orbitals internally,
      ! which invalidate_runs never clears.
      if (have_prev_scf .and. cur_opts%guess_kind == TRC_OPT_GUESS_DENSITY) then
         st = trc_nao(ctx, nao_c)
         allocate (dens_buf(int(nao_c)*int(nao_c)*cached_nspin))
         st = trc_density(ctx, dens_buf)
         if (st /= TRC_OK) then
            why = "could not read the previous step's density for the warm start"
            return
         end if
      end if

      st = trc_set_geometry(ctx, reshape(real(xyz, c_double), [3*natm]))
      if (st /= TRC_OK) then
         why = "trc_set_geometry failed mid-optimisation"
         return
      end if

      if (have_prev_scf .and. cur_opts%guess_kind /= TRC_OPT_GUESS_NONE) then
         if (cur_opts%guess_kind == TRC_OPT_GUESS_DENSITY) then
            st = trc_set_guess(ctx, TRC_GUESS_GIVEN, c_loc(dens_buf), cached_nspin)
         else
            st = trc_set_guess(ctx, TRC_GUESS_PROJECTED, c_null_ptr, cached_nspin)
         end if
         if (st /= TRC_OK) then
            why = "trc_set_guess (warm start) failed mid-optimisation"
            return
         end if
      else if (have_prev_scf) then
         ! TRC_OPT_GUESS_NONE re-runs the ordinary guess every evaluation
         ! with no density, orbitals or DIIS history carried over -- the
         ! baseline the warm-start measurement (--opt-guess density vs
         ! projected) compares against.
         st = trc_set_guess(ctx, cur_opts%fallback_guess_kind, c_null_ptr, 1_c_int)
         if (st /= TRC_OK) then
            why = "trc_set_guess (baseline) failed mid-optimisation"
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
         ! something this evaluation's success should depend on.
         call global_logger%warning("--opt: a warm-started SCF did not converge; retrying from the standard guess")
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
         why = "the SCF did not converge at this geometry"
         return
      end if
      ! From here on there IS a converged SCF on `ctx` to warm-start the
      ! NEXT evaluation's guess from (density or projected orbitals, per
      ! cur_opts%guess_kind) -- set once true and left there for the rest
      ! of this optimisation.
      have_prev_scf = .true.
      st = trc_nspin(ctx, nspin_c)
      cached_nspin = nspin_c
      st = trc_iterations(ctx, niter_c)
      niter = int(niter_c)
      ok = .true.
      last_scf_xyz = xyz   ! cb_get_hessian's own redundant-SCF check, below
      if (allocated(dens_buf)) deallocate (dens_buf)
   end subroutine scf_at_geometry

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
      integer(c_int) :: st, nspin_c
      real(c_double) :: e
      real(c_double), allocatable :: gbuf(:)
      logical :: ok
      character(len=512) :: why
      character(kind=c_char) :: msg_buf(MSG_LEN)

      n_evals = n_evals + 1
      xyz = reshape(real(coords, dp), [3, natm])

      ! Test-only hook: TRC_OPT_TEST_FAIL_AFTER_EVALS, read only here, lets
      ! a test force a failing evaluation at a chosen point deterministically
      ! (never from a real SCF non-convergence, which would be a flaky
      ! trigger) to exercise fail_eval's write-before-error-stop path.
      ! Unset or non-positive: no effect, every build.
      if (test_fail_after_evals() > 0 .and. n_evals >= test_fail_after_evals()) then
         call fail_eval(nvar, energy, gradient, status, &
                        "TRC_OPT_TEST_FAIL_AFTER_EVALS forced a failure at evaluation "//to_string(n_evals))
         return
      end if

      call scf_at_geometry(xyz, last_niter, ok, why)
      if (.not. ok) then
         call fail_eval(nvar, energy, gradient, status, trim(why))
         return
      end if

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
      ! The max CARTESIAN COMPONENT, matching DL-FIND's own test
      ! (dlf_convergence.f90's convergence_set_info: valg=maxval(abs(
      ! gradient))) and what --opt-gtol's own help/the JSON document --
      ! NOT the largest per-atom gradient NORM, which is a different,
      ! larger number for the same gradient whenever more than one
      ! component of an atom's gradient is non-negligible.
      last_gmax = maxval(abs(real(gbuf, dp)))
      last_grms = sqrt(sum(real(gbuf, dp)**2)/real(nvar, dp))

      energy = e
      gradient = gbuf
      status = 0_c_int
      deallocate (gbuf)

   end subroutine cb_get_gradient

   subroutine fail_eval(nvar, energy, gradient, status, why)
      !! One place to mark an evaluation failed and hand DL-FIND a zero
      !! energy/gradient and status = 1 -- DL-FIND's own contract for a
      !! failed evaluation. There is no recovery inside DL-FIND for this:
      !! dl-find.f90's own EANDG block turns a non-zero status straight
      !! into `dlf_fail("Energy evaluation failed")`, which api.f90's
      !! `dlf_error` ends with `error stop` -- api_dl_find NEVER returns on
      !! this path, so trc_dlfind_optimize's own post-call reporting
      !! (res%failed, its message) never runs either. Whatever a finished
      !! run writes to disk (app/main.F90's write_opt_xyz/write_traj_xyz/
      !! add_optimization) has to be written from HERE instead, before
      !! `status = 1` is handed back, while there still is somewhere to
      !! write it.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(out) :: energy
      real(c_double), intent(out) :: gradient(nvar)
      integer(c_int), intent(out) :: status
      character(len=*), intent(in) :: why
      eval_failed = .true.
      eval_fail_msg = why
      call global_logger%error("--opt: "//trim(why))
      call global_logger%error("--opt: DL-FIND has no recovery from a failed evaluation and is about "// &
                               "to stop this process -- writing the last accepted geometry, the "// &
                               "trajectory so far and the JSON optimization block (if requested) now")
      call write_failure_outputs(why)
      energy = 0.0_c_double
      gradient = 0.0_c_double
      status = 1_c_int
   end subroutine fail_eval

   subroutine write_failure_outputs(reason)
      !! See fail_eval's own header for why this runs from inside a
      !! callback rather than after trc_dlfind_optimize returns. Always the
      !! last ACCEPTED geometry/energy (latest_coords/latest_energy) and
      !! the ACCEPTED trajectory collected so far -- never this failed
      !! trial's own coordinates, which cb_get_gradient's caller (DL-FIND)
      !! is about to discard anyway.
      character(len=*), intent(in) :: reason
      if (opt_str_set(cur_opts%opt_xyz_path)) then
         call write_opt_xyz_frame(cur_opts%opt_xyz_path, natm, znuc, latest_energy, latest_coords, &
                                  "NOT CONVERGED / FAILED: "//trim(reason))
      end if
      if (opt_str_set(cur_opts%traj_xyz_path) .and. traj_n > 0) then
         call write_traj_xyz_frames(cur_opts%traj_xyz_path, natm, znuc, traj_n, traj_xyz(:, :, 1:traj_n), &
                                    traj_e(1:traj_n))
      end if
      if (opt_str_set(cur_opts%json_path)) then
         call write_opt_failure_json(cur_opts%json_path, opt_str_or_blank(cur_opts%geometry_path), reason, &
                                     n_steps, n_evals, latest_energy, accepted_gmax, accepted_grms, &
                                     opt_str_or_blank(cur_opts%coords_name), opt_str_or_blank(cur_opts%guess_name), &
                                     traj_e(1:traj_n))
      end if
   end subroutine write_failure_outputs

   pure function opt_str_set(s) result(set)
      !! Whether an optional trc_opt_options_t path field was actually
      !! given a value -- unallocated (a caller that never set it, e.g.
      !! test/bench_opt_warmstart.F90) and "" (set, but asking for no file)
      !! both mean "do not write", and `len()` on an unallocated deferred-
      !! length allocatable is undefined, so this checks `allocated` first.
      character(len=:), allocatable, intent(in) :: s
      logical :: set
      set = allocated(s)
      if (set) set = len(s) > 0
   end function opt_str_set

   pure function opt_str_or_blank(s) result(text)
      character(len=:), allocatable, intent(in) :: s
      character(len=:), allocatable :: text
      if (allocated(s)) then
         text = s
      else
         text = ""
      end if
   end function opt_str_or_blank

   subroutine cb_put_coords(nvar, switch, energy, coords, iam) bind(c)
      !! DL-FIND reporting the geometry it just evaluated.
      !!   switch 1: an ordinary evaluation -- called for EVERY one,
      !!     accepted or not (dl-find.f90 calls it right after
      !!     cb_get_gradient, before test_acceptance/test_acceptance_g
      !!     decides), so whether THIS one is accepted is decided below,
      !!     replicating DL-FIND's own rule for the algorithm in use
      !!     (accept_rule's own comment) -- everything from n_steps down
      !!     is touched only when it is.
      !!   switch 3: the transition structure -- a SEPARATE, redundant
      !!     report of the point the last switch-1 call already recorded
      !!     (dlf_formstep.f90's set_tsmode, called once a saddle search
      !!     converges), so it only snapshots ts_*/ts_gmax/ts_grms from
      !!     that already-accepted state; it never counts as a new step.
      !!     For a `connect` run, the two downhill minimisations that
      !!     follow reset the per-phase acceptance state here, because
      !!     DL-FIND itself resets its own trust radius/oldenergy
      !!     (`toldenergy`) at the start of each one (dlf_task.f90's
      !!     task==1011 case) and unconditionally forces iopt=3/iline=1
      !!     (energy-based acceptance) for both, regardless of what the TS
      !!     search itself used.
      !!   switch 4/5: `connect`'s two downhill minima -- likewise a
      !!     redundant report of each sub-optimisation's own last
      !!     switch-1 point, snapshotted the same way switch 3 is; switch
      !!     5 needs no further phase reset since the run ends after it.
      integer(c_int), intent(in), value :: nvar
      integer(c_int), intent(in), value :: switch
      real(c_double), intent(in), value :: energy
      real(c_double), intent(in) :: coords(nvar)
      integer(c_int), intent(in), value :: iam   !! task id; only task 0 ever reports (ntasks = 1 here anyway)

      if (iam /= 0_c_int) return
      if (nvar /= 3*natm) return

      select case (switch)
      case (DLF_PUT_STEP)
         call consider_step(coords, energy)
      case (DLF_PUT_TS)
         have_ts = .true.
         ts_energy = latest_energy
         ts_coords = latest_coords
         ts_gmax = accepted_gmax
         ts_grms = accepted_grms
         if (cur_opts%connect) then
            accept_rule = ACCEPT_ENERGY
            have_accepted_energy = .false.
         end if
      case (DLF_PUT_MIN_PLUS)
         have_min_plus = .true.
         min_plus_energy = real(energy, dp)
         min_plus_coords = reshape(real(coords, dp), [3, natm])
         min_plus_gmax = accepted_gmax
         min_plus_grms = accepted_grms
         call global_logger%info("   opt-connect: minimum (+distort)   E = "//to_string(min_plus_energy)//" Eh")
         have_accepted_energy = .false.   ! DL-FIND resets toldenergy before the second downhill run too
      case (DLF_PUT_MIN_MINUS)
         have_min_minus = .true.
         min_minus_energy = real(energy, dp)
         min_minus_coords = reshape(real(coords, dp), [3, natm])
         min_minus_gmax = accepted_gmax
         min_minus_grms = accepted_grms
         call global_logger%info("   opt-connect: minimum (-distort)   E = "//to_string(min_minus_energy)//" Eh")
      case default
         ! switch 2 (a TS-mode report, task 1001 -- not offered here) or
         ! anything else DL-FIND might send: nothing to record.
      end select

   end subroutine cb_put_coords

   subroutine consider_step(coords, energy)
      !! An ordinary evaluation (switch 1): total_scf_iters counts it
      !! regardless (trc_opt_result_t%total_scf_iterations's own header --
      !! a rejected trial's SCF iterations are real work too), then
      !! accept_rule decides whether it is ALSO an accepted step --
      !! everything else here (the latest point, the step counter, the
      !! trajectory) only updates when it is.
      real(c_double), intent(in) :: coords(3*natm)
      real(c_double), intent(in), value :: energy
      real(dp) :: wall, e
      logical :: accepted

      total_scf_iters = total_scf_iters + last_niter

      e = real(energy, dp)
      if (accept_rule == ACCEPT_ENERGY) then
         accepted = (.not. have_accepted_energy) .or. (e < last_accepted_energy)
      else
         accepted = .true.
      end if
      if (.not. accepted) then
         call global_logger%verbose("   opt trial rejected: E = "//to_string(e)// &
                                    " Eh (last accepted "//to_string(last_accepted_energy)//" Eh)")
         return
      end if
      have_accepted_energy = .true.
      last_accepted_energy = e

      latest_coords = reshape(real(coords, dp), [3, natm])
      latest_energy = e
      accepted_gmax = last_gmax
      accepted_grms = last_grms
      have_result = .true.
      n_steps = n_steps + 1

      wall = 0.0_dp
      if (step_clock_running) then
         call step_clock%stop()
         wall = step_clock%get_elapsed_time()
      end if
      block
         character(len=1024) :: line
         write (line, "(a,i5,a,f20.12,a,es11.4,a,es11.4,a,i5,a,f10.2,a)") &
            "   opt step ", n_steps, "   E = ", latest_energy, "   max|g| = ", accepted_gmax, &
            "   rms|g| = ", accepted_grms, "   SCF iters ", last_niter, "   t = ", wall, " s"
         call global_logger%info(trim(line))
      end block
      call step_clock%start()
      step_clock_running = .true.

      if (traj_n >= traj_cap) call grow_traj()
      traj_n = traj_n + 1
      traj_e(traj_n) = latest_energy
      traj_xyz(:, :, traj_n) = latest_coords
   end subroutine consider_step

   subroutine grow_traj()
      !! Doubles traj_e/traj_xyz's capacity with move_alloc rather than
      !! drop a frame past traj_cap -- HDLC/DLC's `trerun_energy`
      !! evaluations (dl-find.f90's own COORDS TRANSFORM section) do not
      !! advance `stat%ccycle`, so more accepted steps than
      !! trc_dlfind_optimize's own starting traj_cap is possible in those
      !! coordinate systems.
      real(dp), allocatable :: new_e(:)
      real(dp), allocatable :: new_xyz(:, :, :)
      integer :: new_cap
      new_cap = traj_cap*2
      allocate (new_e(new_cap))
      allocate (new_xyz(3, natm, new_cap))
      new_e(1:traj_cap) = traj_e
      new_xyz(:, :, 1:traj_cap) = traj_xyz
      call move_alloc(new_e, traj_e)
      call move_alloc(new_xyz, traj_xyz)
      traj_cap = new_cap
   end subroutine grow_traj

   subroutine cb_get_hessian(nvar, coords, hessian, status) bind(c)
      !! Second derivatives for P-RFO/NR. When --opt-init-hessian is
      !! "terco" (TRC_OPT_HESS_TERCO), converges the SCF at `coords` (the
      !! same warm-start logic cb_get_gradient uses) and calls terco's own
      !! trc_run_hessian/trc_hessian -- a Cartesian finite-difference
      !! Hessian of terco's OWN analytic gradient, 6*natm+1 SCF+gradient
      !! evaluations, warm-started internally from the SCF this call just
      !! converged (trc_run_hessian's own docstring: "restores the
      !! context" to that reference state once done, so cb_get_gradient's
      !! warm start for the NEXT step is unaffected).
      !!
      !! A non-zero status is a SUPPORTED decline either way:
      !! --opt-init-hessian dlfind never calls trc_run_hessian at all
      !! (inithessian is set to 2 in cb_get_params instead); a convergence
      !! failure at this geometry declines too, so DL-FIND falls back to
      !! its own two-point FD of the gradients it is already requesting
      !! rather than stopping the whole optimisation over a Hessian it
      !! can live without.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: hessian(nvar, nvar)
      integer(c_int), intent(out) :: status

      real(dp) :: xyz(3, natm)
      real(c_double), allocatable :: hbuf(:)
      integer(c_int) :: st
      integer :: niter, ndof
      logical :: ok, already_converged_here
      character(len=512) :: why

      hessian = 0.0_c_double

      if (cur_opts%hessian_source /= TRC_OPT_HESS_TERCO) then
         status = 1_c_int
         return
      end if

      xyz = reshape(real(coords, dp), [3, natm])
      ! DL-FIND asks for the Hessian right after accepting the step that
      ! needed it, at the SAME geometry cb_get_gradient (via
      ! scf_at_geometry) just converged -- skip re-running an identical
      ! SCF when that is exactly what happened, rather than pay for it
      ! twice. Bit-for-bit, not a tolerance: anything else means some
      ! other evaluation ran in between and the warm start is now stale.
      already_converged_here = .false.
      if (allocated(last_scf_xyz)) then
         if (size(last_scf_xyz) == size(xyz)) already_converged_here = all(xyz == last_scf_xyz)
      end if
      if (already_converged_here) then
         niter = 0
         ok = .true.
      else
         call scf_at_geometry(xyz, niter, ok, why)
      end if
      if (.not. ok) then
         call global_logger%warning("--opt: terco's analytic-Hessian callback could not converge the SCF "// &
                                    "at this geometry ("//trim(why)//"); falling back to DL-FIND's own "// &
                                    "finite-difference Hessian")
         status = 1_c_int
         return
      end if

      st = trc_run_hessian(ctx)
      if (st /= TRC_OK) then
         call global_logger%warning("--opt: trc_run_hessian failed at this geometry; falling back to "// &
                                    "DL-FIND's own finite-difference Hessian")
         status = 1_c_int
         return
      end if

      ndof = 3*natm
      allocate (hbuf(ndof*ndof))
      st = trc_hessian(ctx, hbuf)
      if (st /= TRC_OK) then
         call global_logger%warning("--opt: trc_hessian failed after a successful trc_run_hessian; "// &
                                    "falling back to DL-FIND's own finite-difference Hessian")
         status = 1_c_int
         deallocate (hbuf)
         return
      end if

      hessian = reshape(real(hbuf, c_double), [ndof, ndof])
      n_hess_evals = n_hess_evals + 1
      status = 0_c_int
      deallocate (hbuf)
   end subroutine cb_get_hessian

   subroutine cb_get_multistate_gradients(nvar, coords, energy, gradient, coupling, &
                                          needcoupling, iimage, status) bind(c)
      !! Conical-intersection search -- not offered (this module's own
      !! header comment's "NOT exposed" list).
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
                               ", energy evaluations: "//to_string(n_evals)// &
                               ", Hessian evaluations: "//to_string(n_hess_evals)//")")
   end subroutine cb_error

   subroutine cb_update() bind(c)
      !! Called between steps; nothing here needs to happen then.
   end subroutine cb_update

   function test_fail_after_evals() result(n)
      !! TRC_OPT_TEST_FAIL_AFTER_EVALS, cached after the first read (one
      !! `get_environment_variable` call per process, not per evaluation).
      !! -1 (the cache's own "not read yet" AND "unset/unparseable" value)
      !! means disabled -- see cb_get_gradient's own comment for what this
      !! is for.
      integer :: n
      integer, save :: cached = -2   !! -2: not read yet; -1: read, disabled
      character(len=32) :: buf
      integer :: ios, slen, stat
      if (cached == -2) then
         call get_environment_variable("TRC_OPT_TEST_FAIL_AFTER_EVALS", buf, slen, stat)
         cached = -1
         if (stat == 0 .and. slen > 0) then
            read (buf, *, iostat=ios) cached
            if (ios /= 0) cached = -1
         end if
      end if
      n = cached
   end function test_fail_after_evals

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
