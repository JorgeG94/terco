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
! THIS MODULE IS NEVER COMPILED INTO libterco. It lives in src/opt/ and is
! built into its OWN shared library, libterco_opt (src/opt/CMakeLists.txt,
! TERCO_ENABLE_DLFIND=ON only), which is what keeps DL-FIND's
! LGPL-3.0-or-later licence and its `error stop` on a failure (dlf_error,
! below) out of libterco itself -- see cmake/modules/Findlibdlfind.cmake and
! the root CMakeLists.txt's TERCO_ENABLE_DLFIND option. libterco_opt links
! libterco (for trc_capi, below) and libdlfind; terco's own CLI (app/
! main.F90) and any other caller of trc_opt_capi.F90's C ABI (include/
! terco_opt.h) share this exact module -- ONE optimiser driver, not a copy
! each. With TERCO_ENABLE_DLFIND=OFF this file is not built at all;
! app/trc_dlfind_stub.F90 (still in app/, compiled directly into terco_cli)
! is the ONLY implementation of this module's public interface in that
! build, and it is what lets terco_cli refuse --opt with a message instead
! of failing to link.
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
! ACCEPTED-STEP BOOKKEEPING, per icoord/iopt family (consider_step and
! cb_put_coords below implement these three rules, not a single generic
! one):
!   - Single-state (sd/cg/cg-auto/lbfgs/damped/nr/P-RFO/dimer): one
!     dlf_put_coords(switch=1) per DL-FIND cycle, accepted or not
!     (accept_rule's own comment) -- "a step" here means a CYCLE.
!   - NEB: no single accepted step at all -- dlf_put_coords(switch=1)
!     fires once per IMAGE per cycle (cb_get_gradient's own `iimage`),
!     so counting it the single-state way would count every image
!     evaluation as its own step. trc_dlfind_optimize counts an
!     n_steps/n_traj frame once per PATH CYCLE instead (the NEB branch
!     in cb_put_coords, keyed off `last_image`, every image's own latest
!     geometry overwritten in place -- the last one standing once
!     api_dl_find returns IS the converged path).
!   - Global (stochastic/genetic, iopt 51/52): population-based, so
!     every dlf_put_coords(switch=1) IS accepted (accept_rule stays
!     ACCEPT_ALWAYS); but the ENERGY/COORDS it reports belong to a
!     candidate dlf_parallel_opt.f90 evaluated in an EARLIER cycle, not
!     the gradient this bridge's own `last_gmax` holds (the single most
!     recently evaluated candidate, almost never the same one) --
!     measured directly: check_opt_global.F90 reported "NOT converged"
!     from the stale gradient while DL-FIND's own dlf_stoch_endcycle had
!     already judged the ACTUAL reported point converged. consider_step
!     recovers the right one from a small ring buffer keyed by energy
!     (push_global_eval/lookup_global_eval, below) -- exact, since
!     dlf_parallel_opt.f90 carries an energy through its own arrays
!     unchanged, never recomputing it.
!
! WHAT THIS BRIDGE EXPOSES (see app/main.F90's --opt-* options for the CLI
! names) and what it does NOT, with the reason. The single-state
! minimisers, P-RFO/dimer saddle search and `connect` came first; NEB,
! frozen atoms and internal-coordinate constraints below are later
! additions on top of that.
!   - exposed: every single-state minimiser (sd, cg-auto, cg, lbfgs,
!     damped, nr) and P-RFO, both as a minimiser and for a saddle search;
!     all five icoord families (cartesian, hdlc, hdlc-tc, dlc, dlc-tc);
!     the dimer method for a saddle with no Hessian; the `connect` task
!     (1011: TS then both downhill minima); the external-Hessian callback
!     (cb_get_hessian, via trc_run_hessian/trc_hessian) and DL-FIND's own
!     Hessian update/zero-mode/soft-mode/trust/line-search/L-BFGS-memory
!     controls.
!     --- NEB, frozen atoms, constraints and global search ---
!     * NEB and climbing-image NEB (icoord 100-199: 10X free endpoints, 11X
!       endpoints moving perpendicular to the path tangent -- DL-FIND's own
!       default, 12X frozen endpoints; X = the icoord base, so every
!       --opt-coords family NEB offers works for NEB too), with a second
!       geometry (--opt-neb) as DL-FIND's own coords2/nframe=1 product
!       endpoint, --opt-neb-images (nimage, resolved by THIS bridge before
!       api_dl_find runs, not left to DL-FIND's own post-hoc default, so
!       the per-image arrays below can be sized first), --opt-neb-spring
!       (nebk) and --opt-neb-climb (there is no separate on/off switch in
!       DL-FIND's own interface for the climbing image; this bridge forces
!       neb_climb_test to 0.0 -- practically unreachable -- when not asked
!       for, and leaves DL-FIND's own default, 3.0, otherwise). A converged
!       NEB path (trc_opt_result_t%neb_image_energy/neb_image_coords, both
!       endpoints included) and its highest-energy image
!       (neb_climb_image -- with --opt-neb-climb this IS the refined
!       climbing image, so it is also this run's own TS estimate) are
!       always reported, climbing or not.
!     * A warm-start density PER IMAGE (cb_get_gradient's own `iimage`,
!       1-based, one call per image per cycle plus DL-FIND's own scratch
!       slot n_images_dlf once climbing starts): the shared trc context is
!       reused across every image in turn, so by the time an image comes
!       round again `ctx` holds a DIFFERENT image's converged state --
!       reading "the context's last density" (the single-state warm start)
!       would warm-start from the wrong image. test/bench_opt_neb_warmstart.F90
!       measures this on the HCN<->HNC NEB check below: per-image warm
!       start against a single shared density across every image
!       (trc_opt_options_t%neb_shared_density, bench-only).
!       --opt-guess projected is accepted but has no well-defined per-image
!       meaning (it reads the context's own last-converged ORBITALS, one
!       slot, not one per image) and is silently treated as density for
!       --opt-neb, with a one-time warning.
!     * Frozen atoms (--opt-frozen, 1-based, ranges like 1-3,7): spec(i) =
!       -1, already the exact mechanism the part-1 bridge's residue block
!       carries (no bridge code changed at all -- only app/main.F90's own
!       validation, which refuses it together with pure DLC/DLC-TC,
!       dlf_coords.f90's own `minval(spec)<0` check, confirmed by reading
!       it directly).
!     * Internal-coordinate constraints (--opt-constrain, repeatable via a
!       ';'-separated list -- pic_cli's own options are single-valued):
!       bond/angle/torsion/cartesian/bond-difference, DL-FIND's own icons
!       block (type in slot 1, atom indices in slots 2-5 -- verified by
!       reading dlf_hdlc_constraint.f90%assign_cons directly; a CARTESIAN
!       entry is expanded into three icons rows, one per component).
!       HDLC/DLC only (app/main.F90 refuses Cartesian coordinates, and
!       CARTESIAN-type constraints under DLC/DLC-TC, matching
!       dlf_hdlc_interface.f90's own hdlc%internal=(icoord>=3) check).
!       There is NO target-value slot: dlf_hdlc_interface.f90's own caller
!       of assign_cons hardcodes vconstr=0.0 and mconstr=1, so the
!       constrained value is always whatever the START geometry has --
!       not a choice this bridge makes, the only thing DL-FIND's public
!       interface can do.
!     * Global optimisation, iopt 51/52 (dlf_parallel_opt.f90):
!       --opt-algorithm stochastic/genetic, with --opt-seed and the
!       --opt-global-* population/cycle/radius/mutation parameters
!       dlf_parallel_opt.f90 reads (trc_opt_options_t%po_*, dispatched in
!       cb_get_params below). Confirmed single-task-safe by reading
!       dlf_parallel_opt.f90 directly ("call this subroutine even if
!       glob%ntasks == 1", iimage fixed at 1 -- no MPI task farm needed).
!       UNTIL cmake/modules/Findlibdlfind.cmake's pin moved past upstream
!       4167998, every one of dlf_parallel_opt.f90's own three calls to
!       dlf_get_gradient was missing the `kiter` argument that api.f90's
!       own dlf_get_gradient subroutine requires (every OTHER call site
!       in libdlfind, including dl-find.f90's own main loop, passes all
!       seven) -- Positional Fortran argument binding one argument short
!       shifted every later dummy by one slot, so `status` was never
!       written and the real status was read from whatever memory
!       followed the actual argument list, segfaulting before the first
!       energy evaluation on every attempt, independent of molecule,
!       basis or population settings. The pin's own terco-patches branch
!       fixes it (86ca316). The same commit also lets DLF_PO_SEED, an
!       integer environment variable this module sets (terco_app%set_env,
!       trc_dlfind_optimize below) right before api_dl_find runs, replace
!       dlf_parallel_opt.f90's own wall-clock reseed on every call -- with
!       no such variable in the environment DL-FIND's behaviour is
!       unchanged (clock-seeded, as before the patch).
!     * --opt-connect with --opt-saddle-method prfo: used to be refused
!       outright (a bad array copy at dlf_formstep.f90:2362, confirmed
!       with a debug (-fcheck=all) build -- glob%xcoords2 is allocated
!       only when the caller supplies a second geometry, which a connect
!       run never does, so task 1011's own `xcoords2(:,:,1)=tsmode_r`
!       wrote through an unallocated array; the dimer method was
!       unaffected because dlf_dimer allocates xcoords2 for its own
!       direction first). The pin's terco-patches branch fixes this too
!       (3c27977: allocates xcoords2 and sets tcoords2 before task 1011
!       needs them, the same way dlf_dimer already did for the dimer
!       path). test/check_opt_connect.F90 runs connect with both methods.
!   - NOT exposed, and why:
!       * imultistate (conical intersections, Lagrange-Newton, iopt 40) --
!         needs excited-state gradients/couplings terco has no reference
!         for; cb_get_multistate_gradients always declines (status 1).
!       * imicroiter (microiterative QM/MM) -- needs a second, cheaper MM
!         engine terco does not have.
!       * qTS/instanton (iopt 11/12/13, icoord 190) -- investigated
!         directly (dlf_qts.f90, dl-find.f90's icoord==190 special
!         cases): needs a FULL Hessian at every image up front (even the
!         cheapest hessian_mode, `qts%vhessian(varperimage,varperimage,
!         nimage)`), typically 10-40 images, so 10-40 of terco's own FD
!         Hessians (6*natm+1 SCF+gradient evaluations each) before the
!         first qTS step -- not a clean fit for an FD-Hessian-only engine,
!         left out.
!       * ntasks > 1 (parallel multi-image tasks -- NEB runs single-task
!         here, confirmed safe for that), dump/restart checkpointing, the
!         `distort` override for `connect` (DL-FIND's own default of 0.1
!         is used), a second coordinate set (nvar2; only meaningful for
!         NEB/dimer endpoints, and the dimer method here uses DL-FIND's own
!         random/internal initial direction, not a supplied one) -- no CLI
!         need for them yet; dimer's own
!         maxrot/tolrot (rotation count/tolerance) are left at DL-FIND's
!         defaults for the same reason (not asked for, and the dimer
!         rotation cost is already dominated by the gradient callback,
!         not the rotation count); a caller-supplied constraint TARGET
!         value -- not reachable through DL-FIND's own public interface at
!         all (see above).
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
   ! trc_opt_io, not terco_app (app/terco_app.F90): this module is built into
   ! libterco_opt, which must not depend on app/ -- CODE_STRUCTURE.md's
   ! "entry points are leaves" rule, and the only way a caller outside
   ! terco's own build (cyclone) could ever link this at all. trc_opt_io
   ! carries its OWN small copies of these (see that file's own header for
   ! why duplicating three file-writers was preferred over moving them and
   ! leaving app/main.F90's non-opt CLI code without f_string/set_env in a
   ! TERCO_ENABLE_DLFIND=OFF build).
   use trc_opt_io, only: f_string, set_env, write_opt_xyz_frame, write_traj_xyz_frames, write_opt_failure_json
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
   public :: TRC_OPT_ALG_STOCHASTIC, TRC_OPT_ALG_GENETIC
   public :: TRC_OPT_NEB_FREE, TRC_OPT_NEB_PERP, TRC_OPT_NEB_FROZEN
   public :: TRC_OPT_CONS_BOND, TRC_OPT_CONS_ANGLE, TRC_OPT_CONS_TORSION, TRC_OPT_CONS_CARTESIAN, TRC_OPT_CONS_BOND_DIFF
   public :: trc_opt_constraint_t
   public :: trc_opt_options_t, trc_opt_result_t
   public :: trc_dlfind_optimize
   public :: TRC_OPT_STEP_HEADER, TRC_OPT_STEP_RULE

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

   ! --opt-algorithm stochastic/genetic: equal to DL-FIND's own iopt, same
   ! convention as TRC_OPT_ALG_SD etc above (dlf_parallel_opt.f90 -- iopt/10
   ! == 5 is its own dispatch, confirmed single-task-safe: "call this
   ! subroutine even if glob%ntasks == 1", iimage fixed at 1).
   integer, parameter :: TRC_OPT_ALG_STOCHASTIC = 51, TRC_OPT_ALG_GENETIC = 52

   ! --opt-neb-endpoints: the icoord "tens" digit for a NEB run (icoord =
   ! 100 + 10*this + coord_base, verified against dl-find.f90's own
   ! icoord/10==1X dispatch and dlf_default's icoord=-1 fallback, which
   ! picks 110 -- "perpendicular" -- whenever nimage>1, so that is this
   ! enum's own default too).
   integer, parameter :: TRC_OPT_NEB_FREE = 0, TRC_OPT_NEB_PERP = 1, TRC_OPT_NEB_FROZEN = 2

   ! --opt-constrain's own type codes: DL-FIND's `icons(1,:)` (verified by
   ! reading libdlfind's dlf_hdlc_constraint.f90%assign_cons directly --
   ! type in slot 1, atom indices in slots 2-5, NOT the other way around).
   ! CARTESIAN's slot 4 carries the NEGATED component (x=1,y=2,z=3 here;
   ! assign_cons reads mcomp = -icons(4,:)); BOND_DIFF uses only 3 atoms
   ! (i,j,k: bond(i,j) - bond(j,k) held fixed), slot 5 unused.
   integer, parameter :: TRC_OPT_CONS_BOND = 1, TRC_OPT_CONS_ANGLE = 2, TRC_OPT_CONS_TORSION = 3, &
                         TRC_OPT_CONS_CARTESIAN = 4, TRC_OPT_CONS_BOND_DIFF = 5

   type :: trc_opt_constraint_t
      !! One --opt-constrain entry. Atom indices are 1-based, into the SAME
      !! atom order as the geometry; unused slots stay 0. There is no
      !! target-value slot: libdlfind's own assign_cons always takes the
      !! constrained value from the START geometry (vconstr is hardcoded to
      !! 0.0 inside dlf_hdlc_interface.f90's own caller of assign_cons, with
      !! no path from api_dl_find's `spec` to override it) -- "the
      !! constrained value is the start geometry's" is not a choice this
      !! bridge makes, it is the only thing DL-FIND's own public interface
      !! can do.
      integer :: kind = 0   !! TRC_OPT_CONS_*
      integer :: a = 0, b = 0, c = 0, d = 0
         !! bond: a,b. angle: a,b,c (b is the vertex). torsion: a,b,c,d
         !! (b-c is the central bond). cartesian: a = atom, d = component
         !! (1=x,2=y,3=z; negated when written into icons -- see above).
         !! bond-difference: a,b,c (bond(a,b) - bond(b,c)), d unused.
   end type trc_opt_constraint_t

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
   integer(c_int), parameter :: DLF_PUT_GLOBAL_BEST = -1_c_int
      !! --opt-algorithm genetic's own final answer (dlf_parallel_opt.f90's
      !! po_nsave-slot report after its cycles loop, switch = -i, lowest
      !! energy first) -- stochastic reports through DLF_PUT_STEP instead
      !! (dlf_stoch_endcycle calls dlf_put_coords with switch 1 directly,
      !! mid-loop, once it converges or printf>=2).
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

      ! --- NEB / climbing-image NEB (icoord 100-199) ------------------------
      logical :: neb = .false.
      integer :: neb_images = -1        !! sentinel: DL-FIND's own default (10); its own `nimage` directly
      real(dp) :: neb_spring = -1.0_dp  !! sentinel: DL-FIND's own default (0.01 Eh/Bohr^2); its own `nebk` directly
      logical :: neb_climb = .false.
      integer :: neb_endpoints = TRC_OPT_NEB_PERP
      real(dp), allocatable :: neb_product(:, :)   !! (3, natm) Bohr: the second (product) endpoint
      real(dp), allocatable :: neb_via(:, :, :)
         !! (3, natm, n_via) Bohr: intermediate frames handed to DL-FIND
         !! as EXTRA known points along the path (coords2/nframe > 1,
         !! dlf_neb.f90's own "do iframe=1,min(nframe,nimage-1)" loop,
         !! verified directly), in order, BEFORE neb_product -- a bent
         !! via point for two collinear endpoints (--opt-neb-via), so
         !! DL-FIND's own linear-transit guess does not pass one atom
         !! straight through another (see this module's own header).
         !! Unallocated or zero-size: no via points, exactly today's
         !! single-product-endpoint behaviour.
      logical :: neb_shared_density = .false.
         !! Benchmark-only (bench_opt_neb_warmstart.F90; not reachable from
         !! the CLI): warm-start every image from whatever density the
         !! context happens to hold at the time, the way a single-state run
         !! would, instead of each image keeping its OWN stored density --
         !! kept only so the two can be measured against each other.

      ! --- internal-coordinate constraints -----------------------------------
      type(trc_opt_constraint_t), allocatable :: constraints(:)

      ! --- global optimisation (iopt 51/52) ----------------------------------
      ! Sentinels throughout: -1/-1.0 leaves DL-FIND's own default (dl-
      ! find.f90's own inline defaults, not a separate dlf_default_set
      ! routine -- see this module's own header for the numbers).
      integer :: po_pop_size = -1, po_maxcycle = -1, po_init_pop_size = -1, po_nsave = -1, po_distribution = -1
      integer :: po_reset = -1
      real(dp) :: po_radius = -1.0_dp, po_contraction = -1.0_dp, po_tolerance_r = -1.0_dp, po_tolerance_g = -1.0_dp, &
                  po_mutation_rate = -1.0_dp, po_death_rate = -1.0_dp, po_scalefac = -1.0_dp
      integer :: seed = huge(0)
         !! --opt-seed: DLF_PO_SEED, set in THIS process's own environment
         !! (terco_app%set_env) right before api_dl_find runs below --
         !! dlf_parallel_opt.f90's own only path to a reproducible seed
         !! (see this module's own header). Sentinel huge(0): no call to
         !! set_env at all, so a run that did not ask for a seed leaves
         !! DLF_PO_SEED alone -- an unrelated value already in the
         !! environment (another terco_cli invocation's leftover, in the
         !! same shell) IS then read by dlf_parallel_opt.f90, which reads it
         !! unconditionally once it is set; only this bridge's choice to
         !! call set_env or not controls whether that happens.
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
      ! --opt-neb's own path: the final (converged) chain, BOTH endpoints
      ! included, the scratch climbing-image slot (DL-FIND's own nimage-th
      ! position) excluded -- see trc_dlfind_bridge's header for why.
      logical :: have_neb = .false.
      integer :: neb_n_images = 0                        !! DL-FIND's nimage - 1
      real(dp), allocatable :: neb_image_energy(:)       !! (neb_n_images)
      real(dp), allocatable :: neb_image_coords(:, :, :)  !! (3, natm, neb_n_images)
      integer :: neb_climb_image = 0
         !! argmax(neb_image_energy) -- the highest-energy point on the
         !! converged path; with --opt-neb-climb this IS the climbing
         !! image DL-FIND refined, so it is also this run's own TS estimate.
      integer, allocatable :: neb_scf_iterations_per_image(:)   !! (neb_n_images): the warm-start measurement
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
   integer, parameter :: GLOBAL_EVAL_HIST = 4096
      !! Ring buffer capacity for the energy -> gradient lookup below --
      !! generous for --opt-global-pop-size * --opt-global-maxcycle
      !! (just memory, 4096*24 bytes).
   real(dp), save :: global_eval_energy(GLOBAL_EVAL_HIST) = 0.0_dp
   real(dp), save :: global_eval_gmax(GLOBAL_EVAL_HIST) = 0.0_dp
   real(dp), save :: global_eval_grms(GLOBAL_EVAL_HIST) = 0.0_dp
   integer, save :: global_eval_count = 0, global_eval_next = 1
      !! --opt-algorithm stochastic/genetic only (cb_get_gradient's own
      !! guard, below). dlf_parallel_opt.f90 reports its best/saved
      !! candidate through dlf_put_coords with an ENERGY and COORDS it
      !! already has from an EARLIER cb_get_gradient call, not a fresh
      !! one -- by the time that report fires, this bridge's own
      !! last_gmax/last_grms (the single most recent evaluation) almost
      !! always belong to a LATER, different population member, not the
      !! one being reported. Measured: check_opt_global.F90's own
      !! gtol=1e-4 run reported "NOT converged" at max|g|=9.5e-4 from
      !! the stale last-evaluated candidate, while DL-FIND's own
      !! dlf_stoch_endcycle (read directly) had ALREADY judged the
      !! ACTUAL best candidate converged by its own po_tolerance_g test
      !! before ever calling dlf_put_coords. A ring buffer keyed by
      !! ENERGY (carried through dlf_parallel_opt.f90's own arrays bit
      !! for bit, never recomputed) recovers the RIGHT gradient --
      !! push_global_eval/lookup_global_eval below.
   character(len=*), parameter :: TRC_OPT_STEP_HEADER = &
      "      step              E(total)             dE        max|g|        rms|g|  SCF it    t (s)"
      !! consider_step's own step-table header (printed once, mqc's own
      !! style: a label row, then bare aligned columns, never a label
      !! repeated on every data row) -- column widths match consider_
      !! step's own write format exactly. Public so app/main.F90's
      !! opt_summary can close the table with a rule of the same width
      !! once trc_dlfind_optimize returns (consider_step itself has no
      !! way to know which call is the last one).
   character(len=len(TRC_OPT_STEP_HEADER)), parameter :: TRC_OPT_STEP_RULE = &
      "   -----------------------------------------------------------------------------------------"
   type(timer_type), save :: step_clock
      !! consider_step's own "t (s)" column -- (re)started at the top of
      !! every trc_dlfind_optimize call (below), so step 1's own line
      !! times the FIRST evaluation too, not just the gap between later
      !! steps (the previous %start() sat at the BOTTOM of consider_step,
      !! so step 1 always read 0.0 s -- never started before it ran).
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

   ! --opt-neb's own module state: a per-DL-FIND-image warm-start density
   ! (indices 1..n_images_dlf -- the scratch climbing-image slot included,
   ! since cb_get_gradient is called with iimage=n_images_dlf too once
   ! climbing starts, verified against dl-find.f90/dlf_neb.f90), and the
   ! latest accepted geometry/energy per image. cb_put_coords carries no
   ! iimage of its own (verified: dl-find.f90 calls it once per image,
   ! right after that image's own cb_get_gradient, so `last_image` -- set
   ! by the gradient callback -- correctly attributes the put_coords call
   ! that follows it to the same image).
   logical, save :: opt_is_neb = .false.
   integer, save :: n_images_dlf = 0      !! DL-FIND's own nimage
   integer, save :: neb_n_via = 0         !! --opt-neb-via's own frame count (0: none)
   integer, save :: neb_n_frames = 0      !! coords2's own nframe: neb_n_via + 1 (the product), 0 outside --opt-neb
   integer, save :: last_image = 1        !! set by cb_get_gradient, read by cb_put_coords
   integer, save :: neb_prev_image = 0
      !! The previous cb_put_coords(switch=1) call's own `last_image`, for
      !! the NEB branch there to detect a new PATH CYCLE (the per-image
      !! index wrapping back down) rather than count every image's own
      !! report as its own accepted step -- see this module's own header,
      !! "ACCEPTED-STEP BOOKKEEPING".
   logical, allocatable, save :: have_prev_scf_img(:)
   integer(c_int), allocatable, save :: cached_nspin_img(:)
   real(c_double), allocatable, save :: image_densities(:, :)
      !! (nao*nao*nspin, n_images_dlf) -- grown from an initial size of 1
      !! the first time it is actually written (scf_at_geometry's own
      !! comment for why it cannot be sized any earlier).
   real(dp), allocatable, save :: neb_last_coords(:, :, :)    !! (3, natm, n_images_dlf)
   real(dp), allocatable, save :: neb_last_energy(:)          !! (n_images_dlf)
   real(dp), allocatable, save :: neb_last_gmax(:), neb_last_grms(:)   !! (n_images_dlf): that image's own last gradient
   real(dp), allocatable, save :: neb_last_grad(:, :, :)        !! (3, natm, n_images_dlf): that image's own raw gradient
   logical, allocatable, save :: neb_have_image(:)             !! (n_images_dlf)
   integer, allocatable, save :: neb_scf_iters(:)
      !! (n_images_dlf): running total per image, the warm-start measurement
   integer, allocatable, save :: icons_flat(:)                 !! --opt-constrain's own icons block, build_icons' output

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
      integer :: ncons_total, nimg

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
      global_eval_count = 0
      global_eval_next = 1
      neb_prev_image = 0
      call step_clock%start()
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
      ! grow_traj (consider_step) doubles these with move_alloc rather than
      ! drop a frame past this.
      traj_cap = max(opts%maxiter*merge(3, 1, opts%connect) + 8, 8)
      if (allocated(traj_e)) deallocate (traj_e)
      if (allocated(traj_xyz)) deallocate (traj_xyz)
      allocate (traj_e(traj_cap), traj_xyz(3, n, traj_cap))
      traj_n = 0

      ! --- --opt-constrain: expand into DL-FIND's own icons block ----------
      ! (type, atom1..atom4), a CARTESIAN entry expanded into three (one per
      ! component) -- dlf_hdlc_constraint.f90's own assign_cons, read
      ! directly rather than taken on trust: slot 1 is the type, slots 2-5
      ! the atom indices (NOT the other way around).
      call build_icons(opts%constraints, ncons_total)

      ! --- --opt-neb: a per-image warm-start density needs n_images_dlf
      ! fixed BEFORE api_dl_find runs (cb_get_gradient/cb_put_coords use it
      ! on every call), so this bridge resolves DL-FIND's own `nimage`
      ! default (10, dl-find.f90's own inline default block) itself here
      ! rather than leaving it a sentinel for DL-FIND to resolve AFTER
      ! cb_get_params returns, too late for this module's own arrays.
      opt_is_neb = opts%neb
      if (opt_is_neb) then
         nimg = opts%neb_images
         if (nimg < 3) nimg = 10
         n_images_dlf = nimg
         if (cur_opts%guess_kind == TRC_OPT_GUESS_PROJECTED .and. .not. opts%neb_shared_density) then
            ! TRC_GUESS_PROJECTED reads the CONTEXT's own last-converged
            ! orbitals (scf_at_geometry's own comment) -- a single shared
            ! slot, not one per image. For a genuine per-image warm start
            ! that is the wrong previous state once a second image has run,
            ! so NEB always warm-starts by density instead; --opt-guess
            ! projected is accepted but has no effect here.
            call global_logger%warning("--opt-neb: --opt-guess projected has no well-defined per-image "// &
                                       "meaning; warm-starting by density instead")
            cur_opts%guess_kind = TRC_OPT_GUESS_DENSITY
         end if
      else
         n_images_dlf = 1
      end if
      ! image_densities starts at size 1 and is grown the first time a
      ! density is actually persisted (scf_at_geometry, below): trc_nspin
      ! refuses (TRC_ERR_STATE) until the FIRST SCF on this context has
      ! run (cx%scf_ok), which is always still false here -- no SCF has
      ! run yet when --opt starts -- so querying nao*nao*nspin up front
      ! and trusting it silently kept every per-image density at its
      ! size-1 zero placeholder (never matching, never persisted), which
      ! fed trc_set_guess a zero density every time and forced a SAD
      ! retry on nearly every evaluation. Measured by elimination,
      ! comparing bench_opt_neb_warmstart.F90's own two runs: fixed here
      ! rather than ever querying the context before it has a result.
      call reset_image_state(n, n_images_dlf, 1)

      neb_n_via = 0
      if (opt_is_neb .and. allocated(opts%neb_via)) neb_n_via = size(opts%neb_via, 3)
      neb_n_frames = merge(1 + neb_n_via, 0, opt_is_neb)
         !! coords2's own frame count: the via points, in order, THEN the
         !! product endpoint -- dlf_neb.f90's own "do iframe=1,min(nframe,
         !! nimage-1)" loop treats every frame as a KNOWN point along the
         !! path, verified directly (not "nframe-1 via points plus a
         !! special last one": DL-FIND does not distinguish them, this
         !! bridge's own ordering is what makes the LAST frame the product).

      nvar = int(3*n, c_int)
      nvar2 = int(neb_n_frames*3*n, c_int)
         !! dl-find.f90's own check: nvarin2 == nframe*nat*3 (+ nweight +
         !! nmass + n_po_scaling, all 0 here); 0 outside --opt-neb -- no
         !! second structure for any other mode this bridge offers (no
         !! supplied dimer direction).
      nspec = int(3*n + 5*ncons_total, c_int)
         !! [residue(n), znuc(n), icons(5*ncons_total), micspec(n)] --
         !! dl-find.f90's own check: nspec == nat + nz + 5*ncons + 2*nconn + nat (nconn always 0 here).

      res%failed = .false.
      if (opt_is_neb) then
         ! A cheap pre-check, before DL-FIND ever runs: linear Cartesian
         ! interpolation between the start, any via points and the
         ! product can put one atom on top of another for two COLLINEAR
         ! endpoints (HCN<->HNC: H slides straight through C and N) --
         ! the overlap at that image is then singular, which used to
         ! reach LAPACK's own eigensolver and `error stop` there (src/
         ! core/trc_linalg.F90, src/scf/trc_scf_context.F90%sym_orthog --
         ! now a clean error either way, but this check means the user
         ! learns about the geometry, not a diagonalisation info code).
         call check_neb_initial_path(n, znuc, start_coords, opts%neb_via, opts%neb_product, n_images_dlf, &
                                     res%failed, res%message)
      end if

      if (.not. res%failed) then
         ! --opt-seed: the pinned libdlfind's own dlf_parallel_opt.f90
         ! reads DLF_PO_SEED from THIS process's environment (the
         ! terco-patches branch's 86ca316) in place of its wall-clock
         ! reseed, on every call -- the only path a caller outside
         ! libdlfind has to a reproducible global search. Set right
         ! before api_dl_find, not in cb_get_params: dlf_default_init/
         ! dlf_default_set run before the first callback, but the
         ! environment is read inside dlf_parallel_opt.f90 itself, at
         ! the FIRST seeding call, which only --opt-algorithm stochastic/
         ! genetic ever reaches, so setting it unconditionally here is
         ! harmless for every other algorithm too.
         !
         ! The SAME seed also reseeds the Fortran intrinsic generator:
         ! dlf_dimer.f90's own initial-direction guess is
         ! random_number(glob%xcoords2) with no random_seed call of its
         ! own (confirmed by reading it), so it draws from whatever state
         ! the intrinsic generator is already in. libdlfind and this
         ! bridge share one runtime (built by the same compiler, one
         ! process), so seeding it HERE, right before api_dl_find and
         ! after every earlier random_number call this run could have
         ! made, is the only way to put the dimer method's own direction
         ! under --opt-seed too -- there is no "supplied initial
         ! direction" in DL-FIND's public interface (this module's own
         ! header). Harmless for every method that never calls
         ! random_number.
         if (opts%seed /= huge(0)) then
            call set_env("DLF_PO_SEED", to_string(opts%seed))
            block
               integer :: nseed, i
               integer, allocatable :: seed_values(:)
               call random_seed(size=nseed)
               allocate (seed_values(nseed))
               seed_values = [(opts%seed + i, i=1, nseed)]
               call random_seed(put=seed_values)
            end block
         end if
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
            ! printf >= 1 -- leaving it at DL-FIND's own default of 0 lets
            ! a run converge with no geometry ever reported back, which
            ! is why cb_get_params always sets it to at least 1 below.
            ! Reaching here with have_result still false means
            ! api_dl_find returned without ever calling it.
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
      end if

      if (opt_is_neb .and. .not. res%failed) then
         ! The converged path, BOTH endpoints included, the scratch
         ! climbing-image slot (DL-FIND's own image n_images_dlf) excluded
         ! -- dlf_neb.f90's own real chain is indices 1..nimage-1.
         block
            integer :: ni, k
            logical :: have_cimg
            ni = n_images_dlf - 1
            res%have_neb = .true.
            res%neb_n_images = ni
            allocate (res%neb_image_energy(ni), res%neb_image_coords(3, n, ni), &
                     res%neb_scf_iterations_per_image(ni))
            do k = 1, ni
               res%neb_image_energy(k) = neb_last_energy(k)
               res%neb_image_coords(:, :, k) = neb_last_coords(:, :, k)
               res%neb_scf_iterations_per_image(k) = neb_scf_iters(k)
            end do
            res%neb_climb_image = maxloc(res%neb_image_energy, dim=1)
            ! Once climbing actually spawns, DL-FIND freezes the regular
            ! image nearest it (its own cb_get_gradient's iimage argument
            ! stops being called for that index) and instead keeps
            ! refining a SEPARATE point at the scratch slot n_images_dlf
            ! (cb_get_gradient's own docstring: "plus once more at image
            ! n_images_dlf ... once climbing starts") -- excluded from
            ! the loop above on purpose, since it is not part of the
            ! nimage-1 chain. Report THAT refined point as the climbing
            ! image once it has been evaluated at all: the frozen
            ! regular snapshot is measured (--opt-neb-via smoke test,
            ! HCN->HNC) to sit ~5e-7 Eh / several mBohr short of the
            ! saddle that the scratch slot actually converges to.
            have_cimg = cur_opts%neb_climb .and. neb_scf_iters(n_images_dlf) > 0
            if (have_cimg) then
               res%neb_image_energy(res%neb_climb_image) = neb_last_energy(n_images_dlf)
               res%neb_image_coords(:, :, res%neb_climb_image) = neb_last_coords(:, :, n_images_dlf)
               ! Carry the scratch slot's own (small, converged) gradient
               ! into the regular slot too, so the final_gmax/converged
               ! block below judges the climbing image on the point it is
               ! actually reporting, not the frozen snapshot's stale,
               ! pre-climb raw QC gradient.
               neb_last_gmax(res%neb_climb_image) = neb_last_gmax(n_images_dlf)
               neb_last_grms(res%neb_climb_image) = neb_last_grms(n_images_dlf)
               neb_last_grad(:, :, res%neb_climb_image) = neb_last_grad(:, :, n_images_dlf)
            end if
            ! A connect-style override does not apply to NEB (mutually
            ! exclusive at the CLI, app/main.F90): report the climbing
            ! image as the primary final_energy/final_coords, the same way
            ! --opt-connect reports its own TS instead of the last
            ! sub-optimisation's geometry.
            res%final_energy = res%neb_image_energy(res%neb_climb_image)
            res%final_coords = res%neb_image_coords(:, :, res%neb_climb_image)
            block
               ! A chain image's RAW gradient never goes to zero even once
               ! converged: every interior image keeps a component ALONG
               ! the path (the spring force balances it, not the QC
               ! gradient), so judging convergence on neb_last_gmax
               ! directly reports "NOT converged" by construction, every
               ! time, independent of how well the path has actually
               ! settled (measured: check_opt_neb.F90's own climbing run,
               ! raw max|g| ~4.8e-2 at a chain DL-FIND itself had already
               ! frozen). This mirrors dlf_neb_improved_tangent_neb's own
               ! fper (dlf_neb.f90): project out the along-path component
               ! before taking the max, using a central-difference tangent
               ! over the FINAL converged chain (res%neb_image_coords,
               ! climbing override already applied above) -- not a plain
               ! central difference (measured too crude: residual
               ! perpendicular components 2-5e-3 on images a plain
               ! difference called "frozen", against DL-FIND's own
               ! fper ~1e-4 there), but dlf_neb_improved_tangent_neb's
               ! OWN energy-ordered Tau+/Tau-/blend rule (dlf_neb.f90,
               ! read directly): path(k2)=coords(k2+1)-coords(k2) (ONE-
               ! SIDED, forward difference); an image strictly between
               ! its neighbours in energy takes the one-sided tangent
               ! pointing the way energy increases (Tau+ going uphill,
               ! Tau- going downhill); a local extremum (image energy
               ! higher OR lower than both neighbours) blends path(k2)
               ! and path(k2-1), weighted by the SMALLER/LARGER of the
               ! two neighbouring energy gaps (emin/emax), matching the
               ! higher-energy neighbour's own path more. Not DL-FIND's
               ! own climbing-image tangent specifically (a distance-
               ! weighted blend of its two nearest REGULAR images' own
               ! tau, which this reuses via the regular k2==climb_image
               ! case instead -- close enough once the path has settled,
               ! not chased further here).
               !
               ! Frozen endpoints never move at all, so they are excluded
               ! outright (as before). Free endpoints (TRC_OPT_NEB_FREE)
               ! move by the FULL, unprojected gradient in DL-FIND's own
               ! dlf_neb.f90 (icoord's own "0: free" mode), so their raw
               ! gmax/grms, not a perpendicular one, is the right
               ! criterion for them.
               integer :: lo, hi, k2
               real(dp) :: tnorm, emin, emax
               real(dp) :: tvec(3, n), gperp(3, n)
               real(dp) :: path_fwd(3, n, ni)
               real(dp) :: perp_gmax(ni), perp_grms(ni)
               logical :: endpoints_free
               lo = 1
               hi = ni
               if (cur_opts%neb_endpoints == TRC_OPT_NEB_FROZEN .and. ni >= 3) then
                  lo = 2
                  hi = ni - 1
               end if
               endpoints_free = (cur_opts%neb_endpoints == TRC_OPT_NEB_FREE)
               do k2 = 1, ni - 1
                  path_fwd(:, :, k2) = res%neb_image_coords(:, :, k2 + 1) - res%neb_image_coords(:, :, k2)
               end do
               do k2 = lo, hi
                  if (endpoints_free .and. (k2 == 1 .or. k2 == ni)) then
                     perp_gmax(k2) = neb_last_gmax(k2)
                     perp_grms(k2) = neb_last_grms(k2)
                     cycle
                  end if
                  if (k2 == 1) then
                     tvec = path_fwd(:, :, 1)
                  else if (k2 == ni) then
                     tvec = path_fwd(:, :, ni - 1)
                  else if (res%neb_image_energy(k2 + 1) > res%neb_image_energy(k2) .and. &
                          res%neb_image_energy(k2) > res%neb_image_energy(k2 - 1)) then
                     tvec = path_fwd(:, :, k2)
                  else if (res%neb_image_energy(k2 + 1) < res%neb_image_energy(k2) .and. &
                          res%neb_image_energy(k2) < res%neb_image_energy(k2 - 1)) then
                     tvec = path_fwd(:, :, k2 - 1)
                  else
                     emin = min(abs(res%neb_image_energy(k2 + 1) - res%neb_image_energy(k2)), &
                                abs(res%neb_image_energy(k2 - 1) - res%neb_image_energy(k2)))
                     emax = max(abs(res%neb_image_energy(k2 + 1) - res%neb_image_energy(k2)), &
                                abs(res%neb_image_energy(k2 - 1) - res%neb_image_energy(k2)))
                     if (res%neb_image_energy(k2 + 1) < res%neb_image_energy(k2 - 1)) then
                        ! interchange them (dlf_neb.f90's own order: swap
                        ! BEFORE the emin=emax=0 catch below)
                        tnorm = emin
                        emin = emax
                        emax = tnorm
                     end if
                     if (emin == emax .and. emin == 0.0_dp) then
                        emin = 1.0_dp
                        emax = 1.0_dp
                     end if
                     tvec = emax*path_fwd(:, :, k2) + emin*path_fwd(:, :, k2 - 1)
                  end if
                  tnorm = sqrt(sum(tvec**2))
                  if (tnorm > 0.0_dp) then
                     tvec = tvec/tnorm
                     gperp = neb_last_grad(:, :, k2) - sum(neb_last_grad(:, :, k2)*tvec)*tvec
                  else
                     gperp = neb_last_grad(:, :, k2)
                  end if
                  perp_gmax(k2) = maxval(abs(gperp))
                  perp_grms(k2) = sqrt(sum(gperp**2)/real(3*n, dp))
               end do
               res%final_gmax = maxval(perp_gmax(lo:hi))
               res%final_grms = sqrt(sum(perp_grms(lo:hi)**2)/real(hi - lo + 1, dp))
               res%converged = res%final_gmax <= opts%gtol
            end block
         end block
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
      call reset_image_state(0, 0, 1)
      opt_is_neb = .false.
      natm = 0
   end subroutine trc_dlfind_optimize

   subroutine reset_image_state(n, nimg, density_size)
      !! (Re)allocate the per-image warm-start/path bookkeeping for `nimg`
      !! DL-FIND images (1 for every non-NEB run, so the rest of this
      !! module can always index these arrays by image without a branch),
      !! each with room for a density of `density_size` reals (every image
      !! shares the same basis, so this is the same number for all of
      !! them). Called with nimg=0 at the end of trc_dlfind_optimize to
      !! release everything between runs.
      integer, intent(in) :: n, nimg
      integer, intent(in) :: density_size
      if (allocated(have_prev_scf_img)) deallocate (have_prev_scf_img)
      if (allocated(cached_nspin_img)) deallocate (cached_nspin_img)
      if (allocated(image_densities)) deallocate (image_densities)
      if (allocated(neb_last_coords)) deallocate (neb_last_coords)
      if (allocated(neb_last_energy)) deallocate (neb_last_energy)
      if (allocated(neb_last_gmax)) deallocate (neb_last_gmax)
      if (allocated(neb_last_grms)) deallocate (neb_last_grms)
      if (allocated(neb_last_grad)) deallocate (neb_last_grad)
      if (allocated(neb_have_image)) deallocate (neb_have_image)
      if (allocated(neb_scf_iters)) deallocate (neb_scf_iters)
      if (nimg <= 0) return
      allocate (have_prev_scf_img(nimg), source=.false.)
      allocate (cached_nspin_img(nimg), source=1_c_int)
      allocate (image_densities(max(1, density_size), nimg), source=0.0_c_double)
      allocate (neb_last_coords(3, max(n, 1), nimg), source=0.0_dp)
      allocate (neb_last_energy(nimg), source=0.0_dp)
      allocate (neb_last_gmax(nimg), source=0.0_dp)
      allocate (neb_last_grms(nimg), source=0.0_dp)
      allocate (neb_last_grad(3, max(n, 1), nimg), source=0.0_dp)
         !! The raw Cartesian gradient at the last cb_get_gradient call for
         !! each image, kept (not just its gmax/grms scalars) so the
         !! perpendicular-to-path projection below (trc_dlfind_optimize's
         !! own NEB convergence block) can be computed from the FINAL
         !! converged path, after every image has its own ultimate
         !! position -- the path tangent is only meaningful relative to
         !! where images actually settled, not where they started.
      allocate (neb_have_image(nimg), source=.false.)
      allocate (neb_scf_iters(nimg), source=0)
   end subroutine reset_image_state

   subroutine build_icons(constraints, ncons_total)
      !! Flatten --opt-constrain's own list into DL-FIND's icons block
      !! (module `save` icons_flat, consumed by cb_get_params), expanding
      !! a CARTESIAN entry into three (x,y,z) icons rows -- DL-FIND has no
      !! single "freeze this atom's position" constraint type of its own
      !! (dlf_hdlc_constraint.f90's assign_cons takes one Cartesian
      !! COMPONENT at a time).
      type(trc_opt_constraint_t), allocatable, intent(in) :: constraints(:)
      integer, intent(out) :: ncons_total
      integer :: k, comp, p, nc

      nc = 0
      if (allocated(constraints)) nc = size(constraints)
      ncons_total = 0
      do k = 1, nc
         if (constraints(k)%kind == TRC_OPT_CONS_CARTESIAN) then
            ncons_total = ncons_total + 3
         else
            ncons_total = ncons_total + 1
         end if
      end do

      if (allocated(icons_flat)) deallocate (icons_flat)
      allocate (icons_flat(max(1, 5*ncons_total)))
      icons_flat = 0

      p = 0
      do k = 1, nc
         if (constraints(k)%kind == TRC_OPT_CONS_CARTESIAN) then
            do comp = 1, 3
               icons_flat(p + 1) = TRC_OPT_CONS_CARTESIAN
               icons_flat(p + 2) = constraints(k)%a
               icons_flat(p + 3) = 0
               icons_flat(p + 4) = -comp
               icons_flat(p + 5) = 0
               p = p + 5
            end do
         else
            icons_flat(p + 1) = constraints(k)%kind
            icons_flat(p + 2) = constraints(k)%a
            icons_flat(p + 3) = constraints(k)%b
            icons_flat(p + 4) = constraints(k)%c
            icons_flat(p + 5) = constraints(k)%d
            p = p + 5
         end if
      end do
   end subroutine build_icons

   subroutine check_neb_initial_path(n, znuc_in, start, via, product_, nimage, failed, message)
      !! A cheap, APPROXIMATE replay of dlf_neb.f90's own linear-transit
      !! initial-path construction (dlf_neb_definepath, read directly):
      !! every known point (the start, any --opt-neb-via frames, the
      !! product) is placed at an image slot proportional to its
      !! cumulative Euclidean distance along the chain, and the gaps are
      !! filled by linear interpolation between consecutive known points.
      !! Not a bit-exact replica (DL-FIND can additionally re-map/drop a
      !! known point when several land on the same slot, which this never
      !! needs for the handful of points --opt-neb-via realistically
      !! supplies), but it reproduces the ONE failure mode that matters:
      !! two atoms landing on top of each other in an INTERPOLATED image.
      !! Two COLLINEAR endpoints with no via point (HCN<->HNC: the bridging
      !! H slides straight through C and N) is exactly this -- the
      !! resulting image's overlap matrix is singular, which used to run
      !! into LAPACK's own eigensolver and `error stop` there (src/core/
      !! trc_linalg.F90, src/scf/trc_scf_context.F90%sym_orthog -- a clean
      !! error either way now, but naming the image and the atoms here, up
      !! front, is far more useful than a diagonalisation info code).
      integer, intent(in) :: n
      integer, intent(in) :: znuc_in(n)
      real(dp), intent(in) :: start(3, n)
      real(dp), allocatable, intent(in) :: via(:, :, :)
      real(dp), intent(in) :: product_(3, n)
      integer, intent(in) :: nimage   !! DL-FIND's own nimage; the chain has nimage-1 real images
      logical, intent(out) :: failed
      character(len=:), allocatable, intent(out) :: message
      real(dp), parameter :: MIN_SEP = 0.5_dp   !! Bohr

      integer :: n_known, k, img, ia, ja, nv, ni
      real(dp), allocatable :: known(:, :, :), dist_cum(:), path(:, :, :)
      integer, allocatable :: slot(:)
      real(dp) :: total_dist, svar, sep

      failed = .false.
      message = ""
      ni = nimage - 1
      if (ni < 2) return   ! --opt-neb-images' own floor is validated elsewhere; nothing to interpolate here

      nv = 0
      if (allocated(via)) nv = size(via, 3)
      n_known = nv + 2
      allocate (known(3, n, n_known))
      known(:, :, 1) = start
      do k = 1, nv
         known(:, :, 1 + k) = via(:, :, k)
      end do
      known(:, :, n_known) = product_

      allocate (dist_cum(n_known))
      dist_cum(1) = 0.0_dp
      do k = 2, n_known
         dist_cum(k) = dist_cum(k - 1) + norm2(reshape(known(:, :, k) - known(:, :, k - 1), [3*n]))
      end do
      total_dist = dist_cum(n_known)

      allocate (slot(n_known))
      slot(1) = 1
      slot(n_known) = ni
      do k = 2, n_known - 1
         if (total_dist > 0.0_dp) then
            slot(k) = nint(1.0_dp + dist_cum(k)/total_dist*real(ni - 1, dp))
         else
            slot(k) = slot(k - 1) + 1
         end if
         slot(k) = min(ni - 1, max(slot(k - 1) + 1, slot(k)))
      end do

      allocate (path(3, n, ni))
      do k = 1, n_known
         path(:, :, slot(k)) = known(:, :, k)
      end do
      do k = 1, n_known - 1
         do img = slot(k) + 1, slot(k + 1) - 1
            svar = real(img - slot(k), dp)/real(slot(k + 1) - slot(k), dp)
            path(:, :, img) = (1.0_dp - svar)*known(:, :, k) + svar*known(:, :, k + 1)
         end do
      end do

      do img = 1, ni
         do ia = 1, n - 1
            do ja = ia + 1, n
               sep = norm2(path(:, ia, img) - path(:, ja, img))
               if (sep < MIN_SEP) then
                  failed = .true.
                  message = "--opt-neb: the predicted initial path puts atoms "//to_string(ia)//" (Z="// &
                            to_string(znuc_in(ia))//") and "//to_string(ja)//" (Z="//to_string(znuc_in(ja))// &
                            ") only "//to_string(sep)//" Bohr apart at image "//to_string(img)//" of "// &
                            to_string(ni)//" (below "//to_string(MIN_SEP)//" Bohr) -- DL-FIND's own linear "// &
                            "interpolation between the endpoints has put one atom on top of another there "// &
                            "(a collinear pair of endpoints with no intermediate frame is the usual cause: "// &
                            "a bridging atom slides straight through another along the straight-line path). "// &
                            "Add an intermediate geometry with --opt-neb-via that keeps every atom clear of "// &
                            "every other one along the path."
                  return
               end if
            end do
         end do
      end do
   end subroutine check_neb_initial_path

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

      ! `spec` = [residue(natm), znuc(natm), icons(5*ncons), micspec(natm)]
      ! (dl-find.f90's own layout, verified directly). The FIRST block is
      ! residue numbers, never a frozen/free flag: a positive entry is the
      ! residue the atom belongs to (all 1 means "ignore residues", which
      ! Cartesian/DLC do; HDLC needs the real partition -- an all-zero
      ! residue block breaks HDLC outright, avoided here because
      ! residues(:) is never zero by construction, app/main.F90's own
      ! build_opt_residues). A FROZEN atom (app/main.F90's --opt-frozen)
      ! is residues(i) = -1, already carried straight through here with
      ! no translation needed.
      do i = 1, natm
         spec(i) = int(residues(i), c_int)
      end do
      do i = 1, natm
         spec(natm + i) = int(znuc(i), c_int)
      end do
      block
         integer :: ncons_here
         logical :: have_cons
         ncons = 0_c_int
         have_cons = .false.
         if (allocated(icons_flat)) then
            ncons_here = size(icons_flat)/5
            have_cons = (ncons_here > 0) .and. any(icons_flat /= 0)
            if (have_cons) then
               spec(2*natm + 1:2*natm + 5*ncons_here) = int(icons_flat, c_int)
               ncons = int(ncons_here, c_int)
            end if
         end if
         spec(2*natm + 5*int(ncons) + 1:3*natm + 5*int(ncons)) = 1_c_int   ! micspec: DL-FIND requires 0 or 1 here; 1 everywhere
      end block

      nz = int(natm, c_int)       ! lets HDLC/DLC perceive connectivity on their own
      nconn = 0_c_int             ! DL-FIND reads nconn and micspec at inconsistent offsets when both are present; 0 avoids that
      tatoms = 1_c_int            ! the variables are atoms, not arbitrary DOF -- HDLC needs to know

      ! Arrive uninitialised, or at -1 and then checked positive, so each
      ! needs an explicit value even though it looks optional.
      nframe = int(neb_n_frames, c_int)
      nmass = 0_c_int
      nweight = 0_c_int
      ntasks = 1_c_int

      if (opt_is_neb) then
         ! --opt-neb's own via points (if any), THEN the product endpoint
         ! -- dlf_neb.f90's own "do iframe=1,min(nframe,nimage-1)" loop
         ! treats every frame as a known point along the path in order,
         ! so this bridge's own ordering is what makes the LAST frame the
         ! product (this module's own header on trc_opt_options_t%neb_via).
         block
            integer :: f, off
            off = 0
            do f = 1, neb_n_via
               coords2(off + 1:off + 3*natm) = reshape(real(cur_opts%neb_via(:, :, f), c_double), [3*natm])
               off = off + 3*natm
            end do
            coords2(off + 1:off + 3*natm) = reshape(real(cur_opts%neb_product, c_double), [3*natm])
         end block
      end if

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

      if (opt_is_neb) then
         ! icoord = 100 + 10*endpoint_mode + coord_base -- dl-find.f90's own
         ! icoord/10 dispatch (10=free, 11=perpendicular, 12=frozen),
         ! verified directly; app/main.F90 refuses --opt-neb together with
         ! --opt-target saddle/--opt-connect, so this branch and the two
         ! below are mutually exclusive by construction.
         icoord = 100_c_int + 10_c_int*int(cur_opts%neb_endpoints, c_int) + coord_base
         iopt = int(cur_opts%algorithm, c_int)
         nimage = int(n_images_dlf, c_int)
         if (cur_opts%neb_spring >= 0.0_dp) nebk = real(cur_opts%neb_spring, c_double)
         ! --opt-neb-climb: dlf_neb.f90 spawns the climbing image on its
         ! own once every (non-end) image's perpendicular gradient drops
         ! below neb_climb_test*tolerance -- there is no separate on/off
         ! switch in DL-FIND's own interface. Leaving neb_climb_test at its
         ! sentinel (-1, DL-FIND's own default 3.0) lets that happen;
         ! forcing it to 0.0 makes the test's own threshold unreachable in
         ! practice (a gradient component is never exactly zero before
         ! this run's OWN convergence test ends it anyway), which is this
         ! bridge's way of turning climbing off when not asked for.
         if (.not. cur_opts%neb_climb) neb_climb_test = 0.0_c_double
      else if (is_dimer) then
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
      else if (iopt == TRC_OPT_ALG_STOCHASTIC .or. iopt == TRC_OPT_ALG_GENETIC) then
         ! dlf_stoch_endcycle's own report of the best point
         ! (DLF_PUT_GLOBAL_BEST/DLF_PUT_STEP, this module's own cb_put_
         ! coords) is itself gated on printf>=2 (dlf_parallel_opt.f90,
         ! independent of the printf>=1 the rest of DL-FIND uses) --
         ! left at 1 like every other algorithm, stochastic would run to
         ! completion and report "DL-FIND returned without reporting a
         ! geometry" regardless of how it converged. Genetic's own final
         ! report (the po_nsave-slot loop) is unconditional, so 2 is
         ! harmless for it too.
         printf = 2_c_int
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

      ! --opt-algorithm stochastic/genetic (iopt 51/52): dlf_parallel_opt.f90
      ! reads these po_* variables directly -- sentinels here leave
      ! dl-find.f90's own inline defaults (po_pop_size=25, po_radius=1.0,
      ! po_contraction=0.9, po_tolerance_r=1e-8, po_tolerance_g=1e-3,
      ! po_distribution=3, po_maxcycle=10000, po_init_pop_size=2*
      ! po_pop_size, po_reset=500 (genetic only), po_mutation_rate=0.15,
      ! po_death_rate=0.5, po_scalefac=10.0, po_nsave=10). The reproducible
      ! seed (--opt-seed, DLF_PO_SEED) is set in trc_dlfind_optimize,
      ! before api_dl_find runs, not here -- this callback only forwards
      ! the population/cycle/radius/mutation settings.
      if (cur_opts%po_pop_size > 0) po_pop_size = int(cur_opts%po_pop_size, c_int)
      if (cur_opts%po_maxcycle > 0) po_maxcycle = int(cur_opts%po_maxcycle, c_int)
      if (cur_opts%po_init_pop_size > 0) po_init_pop_size = int(cur_opts%po_init_pop_size, c_int)
      if (cur_opts%po_nsave > 0) po_nsave = int(cur_opts%po_nsave, c_int)
      if (cur_opts%po_distribution > 0) po_distribution = int(cur_opts%po_distribution, c_int)
      if (cur_opts%po_reset > 0) po_reset = int(cur_opts%po_reset, c_int)   ! genetic-only; dl-find.f90's own default 500
      if (cur_opts%po_radius > 0.0_dp) po_radius = real(cur_opts%po_radius, c_double)
      if (cur_opts%po_contraction > 0.0_dp) po_contraction = real(cur_opts%po_contraction, c_double)
      if (cur_opts%po_tolerance_r > 0.0_dp) po_tolerance_r = real(cur_opts%po_tolerance_r, c_double)
      if (cur_opts%po_tolerance_g > 0.0_dp) po_tolerance_g = real(cur_opts%po_tolerance_g, c_double)
      if (cur_opts%po_mutation_rate > 0.0_dp) po_mutation_rate = real(cur_opts%po_mutation_rate, c_double)
      if (cur_opts%po_death_rate > 0.0_dp) po_death_rate = real(cur_opts%po_death_rate, c_double)
      if (cur_opts%po_scalefac > 0.0_dp) po_scalefac = real(cur_opts%po_scalefac, c_double)

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
   subroutine scf_at_geometry(xyz, image, niter, ok, why)
      !! `image` is the DL-FIND image this geometry belongs to (always 1
      !! outside --opt-neb). For --opt-neb (and not opts%neb_shared_density,
      !! the benchmark-only override), the warm-start density/orbital-count
      !! bookkeeping below is kept ONE SLOT PER IMAGE
      !! (have_prev_scf_img/cached_nspin_img/image_densities), because the
      !! shared trc context is reused across every image in turn: by the
      !! time this image comes round again, `ctx` holds whatever image ran
      !! most recently, not this one's own last converged state. The
      !! scalar have_prev_scf/cached_nspin pair below is what every other
      !! caller (the single-state path, and --opt-neb with
      !! neb_shared_density for the warm-start measurement) still uses --
      !! exactly today's pre-NEB behaviour.
      real(dp), intent(in) :: xyz(3, natm)
      integer, intent(in) :: image
      integer, intent(out) :: niter
      logical, intent(out) :: ok
      character(len=*), intent(out) :: why

      real(c_double), allocatable, target :: dens_buf(:)
      integer(c_int) :: st, conv, niter_c, nao_c, nspin_c, guess_nspin
      logical :: converged_now, per_image, has_prev

      ok = .false.
      why = ""
      niter = 0

      per_image = opt_is_neb .and. .not. cur_opts%neb_shared_density
      if (per_image) then
         has_prev = have_prev_scf_img(image)
         guess_nspin = cached_nspin_img(image)
      else
         has_prev = have_prev_scf
         guess_nspin = cached_nspin
      end if

      ! The warm-start density has to be read BEFORE the geometry moves:
      ! trc_set_geometry invalidates the context's scf_ok, and trc_density
      ! refuses once that happens. TRC_GUESS_PROJECTED needs no such fetch
      ! -- it reads the context's own last-converged orbitals internally,
      ! which invalidate_runs never clears (never reached when per_image:
      ! trc_dlfind_optimize forces guess_kind to density for --opt-neb).
      if (has_prev .and. cur_opts%guess_kind == TRC_OPT_GUESS_DENSITY) then
         if (per_image) then
            ! THIS image's OWN previously converged density (stored at the
            ! bottom of this routine, the last time it converged) -- not
            ! whatever `ctx` currently holds, which belongs to whichever
            ! image was evaluated last.
            allocate (dens_buf(size(image_densities, 1)))
            dens_buf = image_densities(:, image)
         else
            st = trc_nao(ctx, nao_c)
            allocate (dens_buf(int(nao_c)*int(nao_c)*cached_nspin))
            st = trc_density(ctx, dens_buf)
            if (st /= TRC_OK) then
               why = "could not read the previous step's density for the warm start"
               return
            end if
         end if
      end if

      st = trc_set_geometry(ctx, reshape(real(xyz, c_double), [3*natm]))
      if (st /= TRC_OK) then
         why = "trc_set_geometry failed mid-optimisation"
         return
      end if

      if (has_prev .and. cur_opts%guess_kind /= TRC_OPT_GUESS_NONE) then
         if (cur_opts%guess_kind == TRC_OPT_GUESS_DENSITY) then
            st = trc_set_guess(ctx, TRC_GUESS_GIVEN, c_loc(dens_buf), guess_nspin)
         else
            st = trc_set_guess(ctx, TRC_GUESS_PROJECTED, c_null_ptr, guess_nspin)
         end if
         if (st /= TRC_OK) then
            why = "trc_set_guess (warm start) failed mid-optimisation"
            return
         end if
      else if (has_prev) then
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

      if (.not. converged_now .and. has_prev) then
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
      ! From here on there IS a converged SCF on `ctx` to warm-start a
      ! LATER evaluation's guess from (density or projected orbitals, per
      ! cur_opts%guess_kind).
      st = trc_nspin(ctx, nspin_c)
      st = trc_iterations(ctx, niter_c)
      niter = int(niter_c)
      ok = .true.
      ! --opt-neb has no single "accepted step" the way consider_step counts
      ! for the single-state path (cb_put_coords never calls it for NEB --
      ! see that routine's own comment), so total_scf_iters is kept valid
      ! for NEB by counting every SCF HERE instead, per_image or shared
      ! alike (neb_scf_iters, below, is the additional PER-IMAGE
      ! breakdown, per_image only).
      if (opt_is_neb) total_scf_iters = total_scf_iters + niter
      if (per_image) then
         have_prev_scf_img(image) = .true.
         cached_nspin_img(image) = nspin_c
         neb_scf_iters(image) = neb_scf_iters(image) + niter
         ! Persist THIS image's own converged density for the NEXT time it
         ! is evaluated -- fetched fresh now, since `dens_buf` above (if
         ! allocated at all) still held the OLD one. nao/nspin are only
         ! KNOWN now (a successful SCF just ran), so image_densities is grown
         ! to the real size the first time this fires -- see
         ! trc_dlfind_optimize's own comment on why it cannot be sized
         ! any earlier.
         st = trc_nao(ctx, nao_c)
         block
            integer :: real_size
            real_size = int(nao_c)*int(nao_c)*int(nspin_c)
            if (real_size > 0 .and. real_size /= size(image_densities, 1)) then
               block
                  real(c_double), allocatable :: grown(:, :)
                  integer :: copy_n
                  allocate (grown(real_size, n_images_dlf), source=0.0_c_double)
                  copy_n = min(real_size, size(image_densities, 1))
                  grown(1:copy_n, :) = image_densities(1:copy_n, :)
                  call move_alloc(grown, image_densities)
               end block
            end if
            if (real_size == size(image_densities, 1)) st = trc_density(ctx, image_densities(:, image))
         end block
      else
         have_prev_scf = .true.
         cached_nspin = nspin_c
      end if
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
      integer(c_int), intent(in), value :: iimage
         !! image index; always 1 outside --opt-neb. For --opt-neb, DL-FIND
         !! calls this once per image, 1-based, every cycle (verified
         !! against dl-find.f90's own main loop: one dlf_get_gradient per
         !! image, immediately followed by dlf_put_coords for THAT image --
         !! "last_image" below is what lets cb_put_coords attribute its own
         !! call to the right one), plus once more at image n_images_dlf
         !! (DL-FIND's own scratch slot) once climbing starts.
      integer(c_int), intent(in), value :: kiter    !! microiteration flag; unused
      integer(c_int), intent(out) :: status

      real(dp) :: xyz(3, natm)
      integer(c_int) :: st, nspin_c
      real(c_double) :: e
      real(c_double), allocatable :: gbuf(:)
      logical :: ok
      character(len=512) :: why
      character(kind=c_char) :: msg_buf(MSG_LEN)
      integer :: img, niter_here
      real(dp) :: gmax_here, grms_here

      n_evals = n_evals + 1
      xyz = reshape(real(coords, dp), [3, natm])

      img = 1
      if (opt_is_neb) img = max(1, min(int(iimage), n_images_dlf))
      last_image = img

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

      call scf_at_geometry(xyz, img, niter_here, ok, why)
      last_niter = niter_here
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
      gmax_here = maxval(abs(real(gbuf, dp)))
      grms_here = sqrt(sum(real(gbuf, dp)**2)/real(nvar, dp))
      last_gmax = gmax_here
      last_grms = grms_here
      if (opt_is_neb) then
         neb_last_gmax(img) = gmax_here
         neb_last_grms(img) = grms_here
         neb_last_grad(:, :, img) = reshape(real(gbuf, dp), [3, natm])
      end if
      if (cur_opts%algorithm == TRC_OPT_ALG_STOCHASTIC .or. cur_opts%algorithm == TRC_OPT_ALG_GENETIC) then
         call push_global_eval(real(e, dp), gmax_here, grms_here)
      end if

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

      if (opt_is_neb .and. switch == DLF_PUT_STEP) then
         ! --opt-neb: this callback carries no iimage of its own, so
         ! `last_image` (set by the cb_get_gradient call that always
         ! immediately precedes it, for the SAME image -- verified against
         ! dl-find.f90's own main loop) is what this geometry belongs to.
         ! Recorded as "the latest accepted geometry for this image" rather
         ! than through consider_step's ordinary single-trajectory
         ! bookkeeping, which has no per-image notion at all; each
         ! overwrite is harmless because the last one standing, once
         ! api_dl_find returns, IS the converged path.
         neb_last_coords(:, :, last_image) = reshape(real(coords, dp), [3, natm])
         neb_last_energy(last_image) = real(energy, dp)
         neb_have_image(last_image) = .true.
         have_result = .true.
         ! One accepted step per PATH CYCLE, not per image (this module's
         ! own header, "ACCEPTED-STEP BOOKKEEPING"): a new cycle is
         ! detected here as the per-image index wrapping back down from
         ! the previous call (cb_get_gradient's own docstring: one
         ! dlf_get_gradient/dlf_put_coords pair per image, every cycle,
         ! plus one more at the scratch climbing slot once climbing
         ! starts -- that scratch index is always LARGER than any regular
         ! image, so it never itself looks like a wrap).
         if (last_image <= neb_prev_image) n_steps = n_steps + 1
         neb_prev_image = last_image
         return
      end if

      select case (switch)
      case (DLF_PUT_STEP, DLF_PUT_GLOBAL_BEST)
         ! DLF_PUT_GLOBAL_BEST (-1): dlf_parallel_opt.f90's own genetic
         ! algorithm never calls dlf_put_coords with switch 1 the way
         ! stochastic's dlf_stoch_endcycle does -- its only report is the
         ! po_nsave-slot loop after the cycles loop ends, switch = -i for
         ! the i-th saved minimum, LOWEST energy first (dlf_sort runs
         ! before each slot is saved). -1 is that lowest one; -2 and
         ! below (a second, third, ... local minimum the population also
         ! found) are not this run's own answer and are ignored below.
         ! accept_rule stays ACCEPT_ALWAYS for both global algorithms, so
         ! consider_step records it unconditionally.
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
      real(dp) :: wall, e, de
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
      ! consider_step's own dE column: the PREVIOUS accepted step's energy
      ! -- exactly what have_accepted_energy/last_accepted_energy already
      ! track, read here before this step overwrites them.
      de = 0.0_dp
      if (have_accepted_energy) de = e - last_accepted_energy
      have_accepted_energy = .true.
      last_accepted_energy = e

      latest_coords = reshape(real(coords, dp), [3, natm])
      latest_energy = e
      accepted_gmax = last_gmax
      accepted_grms = last_grms
      if (cur_opts%algorithm == TRC_OPT_ALG_STOCHASTIC .or. cur_opts%algorithm == TRC_OPT_ALG_GENETIC) then
         ! The reported point's OWN gradient, not whatever cb_get_gradient
         ! last happened to compute (this subroutine's own header, the
         ! global-search case) -- see global_eval_count's own comment.
         call lookup_global_eval(e, accepted_gmax, accepted_grms)
      end if
      have_result = .true.
      n_steps = n_steps + 1

      call step_clock%stop()
      wall = step_clock%get_elapsed_time()
      ! Header (and its underline) once, right before the first accepted
      ! step's own row -- not up front in trc_dlfind_optimize, since
      ! --opt-neb never calls this routine at all (cb_put_coords's own
      ! NEB branch) and would otherwise get a table header for a table
      ! that never appears. mqc's own style (read read-only from
      ! metalquicha): a label row, bare aligned columns under it, never
      ! a label repeated on every data row -- app/main.F90's own
      ! opt_summary prints the matching closing rule once the whole run
      ! (trc_dlfind_optimize) returns, since this routine has no way to
      ! know which call is the last one.
      if (n_steps == 1) then
         call global_logger%info(TRC_OPT_STEP_HEADER)
         call global_logger%info(TRC_OPT_STEP_RULE)
      end if
      block
         character(len=1024) :: line
         write (line, "(i7,f22.12,es15.4,es14.4,es14.4,i8,f9.2)") &
            n_steps, latest_energy, de, accepted_gmax, accepted_grms, last_niter, wall
         call global_logger%info("   "//trim(line))
      end block
      call step_clock%start()

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

   subroutine push_global_eval(energy, gmax, grms)
      !! Records one cb_get_gradient evaluation's own (energy, gmax, grms)
      !! into the ring buffer lookup_global_eval reads back from -- see
      !! global_eval_count's own comment for why. Overwrites the oldest
      !! entry once full; a population search re-evaluates the same
      !! candidates every cycle, so the entry this run actually needs is
      !! always recent.
      real(dp), intent(in) :: energy, gmax, grms
      global_eval_energy(global_eval_next) = energy
      global_eval_gmax(global_eval_next) = gmax
      global_eval_grms(global_eval_next) = grms
      global_eval_next = global_eval_next + 1
      if (global_eval_next > GLOBAL_EVAL_HIST) global_eval_next = 1
      global_eval_count = min(global_eval_count + 1, GLOBAL_EVAL_HIST)
   end subroutine push_global_eval

   subroutine lookup_global_eval(energy, gmax, grms)
      !! Finds the MOST RECENT push_global_eval entry whose energy
      !! matches `energy` bit for bit (dlf_parallel_opt.f90 carries an
      !! evaluated candidate's energy through its own arrays unchanged,
      !! never recomputing it, so an exact match is the right test, not
      !! a tolerance) and overwrites gmax/grms with its own. Leaves
      !! gmax/grms untouched (the caller's own last_gmax/last_grms
      !! fallback) when no entry matches -- defensive only; every energy
      !! dlf_put_coords reports here should have come from a push first.
      real(dp), intent(in) :: energy
      real(dp), intent(inout) :: gmax, grms
      integer :: k, j
      do k = 1, global_eval_count
         j = global_eval_next - k
         if (j < 1) j = j + GLOBAL_EVAL_HIST
         if (global_eval_energy(j) == energy) then
            gmax = global_eval_gmax(j)
            grms = global_eval_grms(j)
            return
         end if
      end do
   end subroutine lookup_global_eval

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
         ! image 1: never reached for --opt-neb, which refuses
         ! --opt-init-hessian terco together with --opt-neb at the CLI.
         call scf_at_geometry(xyz, 1, niter, ok, why)
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
