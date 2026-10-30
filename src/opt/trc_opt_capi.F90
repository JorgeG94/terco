!
! C ABI for the geometry optimiser (include/terco_opt.h, generated from this
! file by tools/gen_capi_headers.py exactly as include/terco.h is generated
! from src/api/trc_capi.F90 -- see that script's own header).
!
! This is the ONE place terco's own CLI (app/main.F90) and any other caller
! (cyclone, linking only libterco_opt.so and this header) reach
! trc_dlfind_bridge%trc_dlfind_optimize: main.F90 calls the Fortran module
! procedure directly (it is built by the same compiler, into the same
! library); a C or C++ caller goes through trc_opt_run below instead. Both
! land on the SAME trc_dlfind_optimize call -- one optimiser driver, not a
! copy per caller.
!
! HANDLES. Three kinds, all opaque `trc_handle_t` (void*) at the C boundary,
! same convention as terco.h:
!   - the context handle trc_create (terco.h) already returned -- this file
!     adds no context operations of its own, it only hands the same handle
!     on to trc_dlfind_optimize. The caller sets up molecule, basis and
!     method through terco.h first, exactly as it would for an energy or a
!     gradient; trc_opt_run only ever moves the atoms on it.
!   - an options handle (trc_opt_options_create/destroy, trc_opt_set_*):
!     wraps trc_dlfind_bridge%trc_opt_options_t, which is not itself
!     interoperable (allocatable components), plus the one extra field
!     (frozen atom list) that CLI-side code builds OUTSIDE that type
!     (app/main.F90's own --opt-frozen handling, folded into the residue
!     array trc_dlfind_optimize takes, not into trc_opt_options_t). A null
!     options handle means "every default", same as a freshly constructed
!     trc_opt_options_t.
!   - a result handle (trc_opt_result_create/destroy, trc_opt_result_*):
!     wraps trc_opt_result_t. Fetched field by field rather than as a
!     struct, for the same interoperability reason, and because the
!     trajectory is variable length (trc_opt_result_trajectory_count first,
!     then the caller's own correctly sized buffers).
!
! WHAT THIS v1 ABI DOES NOT YET EXPOSE. trc_opt_options_t also carries NEB,
! `--opt-connect`, internal-coordinate constraints and the stochastic/
! genetic global search (trc_dlfind_bridge.F90's own header lists what each
! one is for) -- none of those have a setter here yet, only what was asked
! for: task, algorithm, coordinates, step/tolerance limits, the initial
! Hessian source and frozen atoms. trc_opt_run always takes the
! single-state path (never --opt-neb/--opt-connect); terco's own CLI still
! reaches the rest directly through trc_opt_options_t/trc_dlfind_optimize
! (app/main.F90), which is why those fields stay on that type rather than
! being removed. Extending this ABI is adding setters/getters here, not
! touching trc_dlfind_bridge.F90 at all.
!
! DL-FIND'S `error stop`. A failed SCF evaluation mid-optimisation writes
! what it can (the last accepted geometry/trajectory/JSON, if the caller
! asked for any of those paths) and then DL-FIND's OWN api.f90 calls
! `error stop` -- trc_opt_run, like app/main.F90's own call, can simply
! never return on that path. This is libdlfind's contract, not a bug this
! ABI works around (see trc_dlfind_bridge.F90's own fail_eval); a caller
! that cannot tolerate the HOST PROCESS stopping must not call trc_opt_run
! in that process at all.
!
module trc_opt_capi
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_ptr, c_null_char, c_associated, &
                                          c_loc, c_f_pointer
   use pic_types, only: dp
   use trc_capi, only: TRC_OK, TRC_ERR_NULL, TRC_ERR_BADARG, TRC_ERR_NOCONV, TRC_ERR_STATE
   use trc_dlfind_bridge, only: trc_opt_options_t, trc_opt_result_t, trc_dlfind_optimize, trc_dlfind_available, &
      coords_cartesian_ => TRC_OPT_COORDS_CARTESIAN, coords_hdlc_ => TRC_OPT_COORDS_HDLC, &
      coords_hdlc_tc_ => TRC_OPT_COORDS_HDLC_TC, coords_dlc_ => TRC_OPT_COORDS_DLC, &
      coords_dlc_tc_ => TRC_OPT_COORDS_DLC_TC, &
      guess_density_ => TRC_OPT_GUESS_DENSITY, guess_projected_ => TRC_OPT_GUESS_PROJECTED, &
      guess_none_ => TRC_OPT_GUESS_NONE, &
      alg_sd_ => TRC_OPT_ALG_SD, alg_cg_auto_ => TRC_OPT_ALG_CG_AUTO, alg_cg_ => TRC_OPT_ALG_CG, &
      alg_lbfgs_ => TRC_OPT_ALG_LBFGS, alg_prfo_ => TRC_OPT_ALG_PRFO, alg_nr_ => TRC_OPT_ALG_NR, &
      alg_damped_ => TRC_OPT_ALG_DAMPED, &
      target_minimum_ => TRC_OPT_TARGET_MINIMUM, target_saddle_ => TRC_OPT_TARGET_SADDLE, &
      saddle_prfo_ => TRC_OPT_SADDLE_PRFO, saddle_dimer_ => TRC_OPT_SADDLE_DIMER, &
      hess_terco_ => TRC_OPT_HESS_TERCO, hess_dlfind_ => TRC_OPT_HESS_DLFIND
   use trc_constants, only: trc_covalent_radius
   use trc_error, only: error_t
   implicit none
   private

   public :: trc_opt_available
   public :: trc_opt_options_create, trc_opt_options_destroy
   public :: trc_opt_set_target, trc_opt_set_algorithm, trc_opt_set_coords, trc_opt_set_saddle_method
   public :: trc_opt_set_max_steps, trc_opt_set_gtol, trc_opt_set_max_step
   public :: trc_opt_set_hessian_source, trc_opt_set_guess, trc_opt_set_frozen_atoms
   public :: trc_opt_result_create, trc_opt_result_destroy
   public :: trc_opt_result_converged, trc_opt_result_failed, trc_opt_result_steps
   public :: trc_opt_result_energy_evaluations, trc_opt_result_energy
   public :: trc_opt_result_gradient_max, trc_opt_result_gradient_rms
   public :: trc_opt_result_geometry, trc_opt_result_trajectory_count, trc_opt_result_trajectory
   public :: trc_opt_result_message
   public :: trc_opt_run

   ! --opt-coords (TRC_OPT_COORDS_*): equal to trc_dlfind_bridge's own
   ! icoord base, re-exported here, by value, the same way trc_capi.F90
   ! re-exports trc_ao_norm's TRC_AO_NORM_SHELL/UNIT -- see that file's own
   ! header, and tools/gen_capi_headers.py's `find_module_constant`, which
   ! is what resolves these for include/terco_opt.h/trc_opt_c_interfaces.f90.
   integer(c_int), parameter, public :: TRC_OPT_COORDS_CARTESIAN = coords_cartesian_
   integer(c_int), parameter, public :: TRC_OPT_COORDS_HDLC = coords_hdlc_
   integer(c_int), parameter, public :: TRC_OPT_COORDS_HDLC_TC = coords_hdlc_tc_
   integer(c_int), parameter, public :: TRC_OPT_COORDS_DLC = coords_dlc_
   integer(c_int), parameter, public :: TRC_OPT_COORDS_DLC_TC = coords_dlc_tc_

   ! --opt-guess (TRC_OPT_GUESS_*).
   integer(c_int), parameter, public :: TRC_OPT_GUESS_DENSITY = guess_density_
   integer(c_int), parameter, public :: TRC_OPT_GUESS_PROJECTED = guess_projected_
   integer(c_int), parameter, public :: TRC_OPT_GUESS_NONE = guess_none_

   ! --opt-algorithm (TRC_OPT_ALG_*), the single-state minimisers only --
   ! stochastic/genetic (iopt 51/52) are not offered through this ABI yet
   ! (this module's own header).
   integer(c_int), parameter, public :: TRC_OPT_ALG_SD = alg_sd_
   integer(c_int), parameter, public :: TRC_OPT_ALG_CG_AUTO = alg_cg_auto_
   integer(c_int), parameter, public :: TRC_OPT_ALG_CG = alg_cg_
   integer(c_int), parameter, public :: TRC_OPT_ALG_LBFGS = alg_lbfgs_
   integer(c_int), parameter, public :: TRC_OPT_ALG_PRFO = alg_prfo_
   integer(c_int), parameter, public :: TRC_OPT_ALG_NR = alg_nr_
   integer(c_int), parameter, public :: TRC_OPT_ALG_DAMPED = alg_damped_

   ! --opt-target / --opt-saddle-method.
   integer(c_int), parameter, public :: TRC_OPT_TARGET_MINIMUM = target_minimum_
   integer(c_int), parameter, public :: TRC_OPT_TARGET_SADDLE = target_saddle_
   integer(c_int), parameter, public :: TRC_OPT_SADDLE_PRFO = saddle_prfo_
   integer(c_int), parameter, public :: TRC_OPT_SADDLE_DIMER = saddle_dimer_

   ! --opt-init-hessian.
   integer(c_int), parameter, public :: TRC_OPT_HESS_TERCO = hess_terco_
   integer(c_int), parameter, public :: TRC_OPT_HESS_DLFIND = hess_dlfind_

   type :: trc_opt_options_box
      type(trc_opt_options_t) :: opts
      integer, allocatable :: frozen_atoms(:)
         !! 1-based, into trc_opt_run's own atom order -- app/main.F90's
         !! --opt-frozen (build_opt_residues' own comment: folded into the
         !! residue array AT RUN TIME, residues(i) = -1, never into
         !! trc_opt_options_t itself).
   end type trc_opt_options_box

   type :: trc_opt_result_box
      type(trc_opt_result_t) :: res
      integer :: natm = 0
         !! Captured by trc_opt_run, checked against every getter's own
         !! `natm`/`n_traj` argument below so a caller's mismatched buffer
         !! is TRC_ERR_BADARG, never an out-of-bounds Fortran copy.
   end type trc_opt_result_box

contains

   function trc_opt_available() result(available) bind(c, name="trc_opt_available")
      !! Whether THIS build of libterco_opt can optimise a geometry at
      !! all -- always true: with TERCO_ENABLE_DLFIND=OFF the library
      !! itself is not built, so there is no trc_opt_available to call in
      !! the first place. Kept for the same reason app/main.F90 still
      !! calls trc_dlfind_bridge%trc_dlfind_available() itself rather than
      !! assuming: a caller that detects this library at runtime (dlopen)
      !! rather than at link time can ask instead of guessing.
      integer(c_int) :: available
      available = merge(1_c_int, 0_c_int, trc_dlfind_available())
   end function trc_opt_available

   function trc_opt_options_create(handle) result(status) bind(c, name="trc_opt_options_create")
      !! A fresh options handle, every field at trc_opt_options_t's own
      !! default (app/main.F90's unoptimised run is --opt-algorithm lbfgs,
      !! --opt-coords cartesian, --opt-gtol 4.5e-4, ... -- trc_dlfind_bridge.F90's
      !! own type definition).
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      allocate (ob)
      handle = c_loc(ob)
      status = TRC_OK
   end function trc_opt_options_create

   function trc_opt_options_destroy(handle) result(status) bind(c, name="trc_opt_options_destroy")
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      deallocate (ob)
      status = TRC_OK
   end function trc_opt_options_destroy

   function trc_opt_set_target(handle, target_kind) result(status) bind(c, name="trc_opt_set_target")
      !! TRC_OPT_TARGET_MINIMUM or TRC_OPT_TARGET_SADDLE (--opt-target).
      type(c_ptr), value :: handle
      integer(c_int), value :: target_kind
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      status = TRC_ERR_BADARG
      if (target_kind /= TRC_OPT_TARGET_MINIMUM .and. target_kind /= TRC_OPT_TARGET_SADDLE) return
      ob%opts%target_kind = int(target_kind)
      status = TRC_OK
   end function trc_opt_set_target

   function trc_opt_set_algorithm(handle, algorithm) result(status) bind(c, name="trc_opt_set_algorithm")
      !! TRC_OPT_ALG_* (--opt-algorithm). Only meaningful for a minimum
      !! search (--opt-target minimum); a saddle search picks its own
      !! algorithm from --opt-saddle-method instead, same as the CLI.
      type(c_ptr), value :: handle
      integer(c_int), value :: algorithm
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      ob%opts%algorithm = int(algorithm)
      status = TRC_OK
   end function trc_opt_set_algorithm

   function trc_opt_set_coords(handle, coords_kind) result(status) bind(c, name="trc_opt_set_coords")
      !! TRC_OPT_COORDS_* (--opt-coords).
      type(c_ptr), value :: handle
      integer(c_int), value :: coords_kind
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      ob%opts%coords_kind = int(coords_kind)
      status = TRC_OK
   end function trc_opt_set_coords

   function trc_opt_set_saddle_method(handle, method) result(status) bind(c, name="trc_opt_set_saddle_method")
      !! TRC_OPT_SADDLE_PRFO or TRC_OPT_SADDLE_DIMER (--opt-saddle-method).
      type(c_ptr), value :: handle
      integer(c_int), value :: method
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      ob%opts%saddle_method = int(method)
      status = TRC_OK
   end function trc_opt_set_saddle_method

   function trc_opt_set_max_steps(handle, maxiter) result(status) bind(c, name="trc_opt_set_max_steps")
      !! --opt-max-steps: the cycle cap (DL-FIND's own `maxcycle`).
      type(c_ptr), value :: handle
      integer(c_int), value :: maxiter
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      status = TRC_ERR_BADARG
      if (maxiter < 1) return
      ob%opts%maxiter = int(maxiter)
      status = TRC_OK
   end function trc_opt_set_max_steps

   function trc_opt_set_gtol(handle, gtol) result(status) bind(c, name="trc_opt_set_gtol")
      !! --opt-gtol, Eh/Bohr: the max CARTESIAN COMPONENT of the gradient
      !! at the last accepted step, same test trc_dlfind_bridge's own
      !! convergence verdict uses (never DL-FIND's internal one -- that
      !! file's own header).
      type(c_ptr), value :: handle
      real(c_double), value :: gtol
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      status = TRC_ERR_BADARG
      if (gtol <= 0.0_c_double) return
      ob%opts%gtol = real(gtol, dp)
      status = TRC_OK
   end function trc_opt_set_gtol

   function trc_opt_set_max_step(handle, max_step) result(status) bind(c, name="trc_opt_set_max_step")
      !! --opt-max-step, Bohr: the trust radius / max Cartesian step
      !! DL-FIND's own `maxstep` enforces per cycle. A non-positive value
      !! leaves DL-FIND's own default (0.5 Bohr, trc_opt_options_t's own
      !! sentinel for "not set").
      type(c_ptr), value :: handle
      real(c_double), value :: max_step
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      ob%opts%max_step = real(max_step, dp)
      status = TRC_OK
   end function trc_opt_set_max_step

   function trc_opt_set_hessian_source(handle, source) result(status) bind(c, name="trc_opt_set_hessian_source")
      !! TRC_OPT_HESS_TERCO or TRC_OPT_HESS_DLFIND (--opt-init-hessian),
      !! for P-RFO/NR only -- ignored by every other algorithm, same as
      !! the CLI.
      type(c_ptr), value :: handle
      integer(c_int), value :: source
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      status = TRC_ERR_BADARG
      if (source /= TRC_OPT_HESS_TERCO .and. source /= TRC_OPT_HESS_DLFIND) return
      ob%opts%hessian_source = int(source)
      status = TRC_OK
   end function trc_opt_set_hessian_source

   function trc_opt_set_guess(handle, guess_kind) result(status) bind(c, name="trc_opt_set_guess")
      !! TRC_OPT_GUESS_* (--opt-guess): the warm start between steps.
      type(c_ptr), value :: handle
      integer(c_int), value :: guess_kind
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      ob%opts%guess_kind = int(guess_kind)
      status = TRC_OK
   end function trc_opt_set_guess

   function trc_opt_set_frozen_atoms(handle, n, atoms) result(status) bind(c, name="trc_opt_set_frozen_atoms")
      !! --opt-frozen: `atoms` holds `n` 1-based indices into trc_opt_run's
      !! own atom order. Replaces whatever a previous call set (not
      !! additive): a caller that wants to grow the list reads it back
      !! itself first -- there is no getter for it, cheap enough to track
      !! on the caller's own side instead.
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      integer(c_int), intent(in) :: atoms(*)
      integer(c_int) :: status
      type(trc_opt_options_box), pointer :: ob
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, ob)
      status = TRC_ERR_BADARG
      if (n < 0) return
      if (allocated(ob%frozen_atoms)) deallocate (ob%frozen_atoms)
      allocate (ob%frozen_atoms(n), source=int(atoms(1:n)))
      status = TRC_OK
   end function trc_opt_set_frozen_atoms

   function trc_opt_result_create(handle) result(status) bind(c, name="trc_opt_result_create")
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      allocate (rb)
      handle = c_loc(rb)
      status = TRC_OK
   end function trc_opt_result_create

   function trc_opt_result_destroy(handle) result(status) bind(c, name="trc_opt_result_destroy")
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      deallocate (rb)
      status = TRC_OK
   end function trc_opt_result_destroy

   function trc_opt_result_converged(handle, converged) result(status) bind(c, name="trc_opt_result_converged")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: converged
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      converged = merge(1_c_int, 0_c_int, rb%res%converged)
      status = TRC_OK
   end function trc_opt_result_converged

   function trc_opt_result_failed(handle, failed) result(status) bind(c, name="trc_opt_result_failed")
      !! Whether this run stopped on a recoverable failure (an SCF that
      !! would not converge mid-optimisation) BEFORE DL-FIND's own
      !! `error stop` could fire -- see this module's own header for why
      !! that is not the only way an evaluation failure ends a run.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: failed
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      failed = merge(1_c_int, 0_c_int, rb%res%failed)
      status = TRC_OK
   end function trc_opt_result_failed

   function trc_opt_result_steps(handle, steps) result(status) bind(c, name="trc_opt_result_steps")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: steps
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      steps = int(rb%res%steps, c_int)
      status = TRC_OK
   end function trc_opt_result_steps

   function trc_opt_result_energy_evaluations(handle, n) &
      result(status) bind(c, name="trc_opt_result_energy_evaluations")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: n
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      n = int(rb%res%energy_evaluations, c_int)
      status = TRC_OK
   end function trc_opt_result_energy_evaluations

   function trc_opt_result_energy(handle, energy) result(status) bind(c, name="trc_opt_result_energy")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: energy
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      energy = real(rb%res%final_energy, c_double)
      status = TRC_OK
   end function trc_opt_result_energy

   function trc_opt_result_gradient_max(handle, gmax) result(status) bind(c, name="trc_opt_result_gradient_max")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: gmax
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      gmax = real(rb%res%final_gmax, c_double)
      status = TRC_OK
   end function trc_opt_result_gradient_max

   function trc_opt_result_gradient_rms(handle, grms) result(status) bind(c, name="trc_opt_result_gradient_rms")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: grms
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      grms = real(rb%res%final_grms, c_double)
      status = TRC_OK
   end function trc_opt_result_gradient_rms

   function trc_opt_result_geometry(handle, natm, xyz) result(status) bind(c, name="trc_opt_result_geometry")
      !! The last ACCEPTED geometry, atom-major Bohr (x1 y1 z1 x2 ...),
      !! same convention as terco.h's trc_set_geometry/trc_gradient --
      !! `xyz` is the caller's own buffer, length 3*natm.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm
      real(c_double), intent(out) :: xyz(*)
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      status = TRC_ERR_BADARG
      if (int(natm) /= rb%natm .or. .not. allocated(rb%res%final_coords)) return
      xyz(1:3*rb%natm) = reshape(real(rb%res%final_coords, c_double), [3*rb%natm])
      status = TRC_OK
   end function trc_opt_result_geometry

   function trc_opt_result_trajectory_count(handle, n_traj) &
      result(status) bind(c, name="trc_opt_result_trajectory_count")
      !! The number of ACCEPTED-step frames available from
      !! trc_opt_result_trajectory below (0: none, or the caller asked for
      !! no trajectory to be kept -- there is no such switch in this v1
      !! ABI, the trajectory trc_dlfind_optimize always collects is always
      !! cheap enough to keep, same as app/main.F90's own --opt).
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: n_traj
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      n_traj = int(rb%res%n_traj, c_int)
      status = TRC_OK
   end function trc_opt_result_trajectory_count

   function trc_opt_result_trajectory(handle, natm, n_traj, xyz, energy) &
      result(status) bind(c, name="trc_opt_result_trajectory")
      !! `xyz` is the caller's own buffer, length 3*natm*n_traj, frame
      !! `k` (1-based) at offset 3*natm*(k-1); `energy` is length n_traj.
      !! Both `natm` and `n_traj` must match trc_opt_run's own atom count
      !! and trc_opt_result_trajectory_count's own answer exactly.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm, n_traj
      real(c_double), intent(out) :: xyz(*)
      real(c_double), intent(out) :: energy(*)
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      status = TRC_ERR_BADARG
      if (int(natm) /= rb%natm .or. int(n_traj) /= rb%res%n_traj) return
      if (rb%res%n_traj > 0) then
         if (.not. allocated(rb%res%traj_coords) .or. .not. allocated(rb%res%traj_energy)) return
         xyz(1:3*rb%natm*rb%res%n_traj) = reshape(real(rb%res%traj_coords, c_double), &
                                                   [3*rb%natm*rb%res%n_traj])
         energy(1:rb%res%n_traj) = real(rb%res%traj_energy, c_double)
      end if
      status = TRC_OK
   end function trc_opt_result_trajectory

   function trc_opt_result_message(handle, buffer, buflen) result(status) bind(c, name="trc_opt_result_message")
      !! The run's own message (why it failed, or "" for an ordinary
      !! finish), NUL-terminated into `buffer` of `buflen` chars -- same
      !! contract as terco.h's trc_message.
      type(c_ptr), value :: handle
      character(kind=c_char), intent(out) :: buffer(*)
      integer(c_int), value :: buflen
      integer(c_int) :: status
      type(trc_opt_result_box), pointer :: rb
      integer :: i, n
      character(len=:), allocatable :: msg
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, rb)
      msg = ""
      if (allocated(rb%res%message)) msg = rb%res%message
      n = min(len_trim(msg), int(buflen) - 1)
      do i = 1, n
         buffer(i) = msg(i:i)
      end do
      if (buflen > 0) buffer(n + 1) = c_null_char
      status = TRC_OK
   end function trc_opt_result_message

   function trc_opt_run(handle, natm, z, xyz, options, result_handle) result(status) bind(c, name="trc_opt_run")
      !! Optimise (or find a saddle on) the molecule already set up on
      !! `handle` through terco.h (molecule, basis and method, matching
      !! `z`/`xyz` -- same requirement as trc_dlfind_bridge%trc_dlfind_optimize,
      !! which this calls directly). `z` (nuclear charges) and `xyz`
      !! (atom-major Bohr, length 3*natm) are the caller's OWN copy of what
      !! it already gave trc_set_molecule -- there is no getter for them on
      !! the context handle (terco.h has none), so the caller hands them in
      !! again rather than this ABI reading them back out of an opaque
      !! object.
      !!
      !! `options` may be a null handle (every default); `result_handle`
      !! must already exist (trc_opt_result_create) and is overwritten in
      !! place. Returns TRC_OK when converged, TRC_ERR_NOCONV when the run
      !! finished without converging (the usual "ran out of steps"),
      !! TRC_ERR_STATE when an evaluation failed outright (the message is
      !! on result_handle either way, trc_opt_result_message) -- but see
      !! this module's own header: DL-FIND's `error stop` can still end the
      !! calling PROCESS before this function ever returns at all.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm
      integer(c_int), intent(in) :: z(*)
      real(c_double), intent(in) :: xyz(*)
      type(c_ptr), value :: options
      type(c_ptr), value :: result_handle
      integer(c_int) :: status

      type(trc_opt_options_box), pointer :: ob
      type(trc_opt_result_box), pointer :: rb
      type(trc_opt_options_t) :: opts
      integer, allocatable :: residues(:), z_in(:)
      real(dp), allocatable :: xyz_in(:, :)
      type(error_t) :: err

      status = TRC_ERR_NULL
      if (.not. c_associated(handle) .or. .not. c_associated(result_handle)) return
      call c_f_pointer(result_handle, rb)

      status = TRC_ERR_BADARG
      if (int(natm) < 1) return
      rb%natm = int(natm)

      if (c_associated(options)) then
         call c_f_pointer(options, ob)
         opts = ob%opts
      end if

      allocate (z_in(natm), source=int(z(1:natm)))
      allocate (xyz_in(3, natm), source=reshape(real(xyz(1:3*natm), dp), [3, int(natm)]))

      allocate (residues(natm))
      if (opts%coords_kind == TRC_OPT_COORDS_HDLC .or. opts%coords_kind == TRC_OPT_COORDS_HDLC_TC) then
         call build_residues(int(natm), z_in, xyz_in, residues, err)
         if (err%has_error()) return
      else
         residues = 1   ! cartesian/dlc/dlc-tc: trc_dlfind_bridge needs no real partition
      end if
      if (c_associated(options)) then
         if (allocated(ob%frozen_atoms)) then
            status = TRC_ERR_BADARG
            if (any(ob%frozen_atoms < 1) .or. any(ob%frozen_atoms > int(natm))) return
            residues(ob%frozen_atoms) = -1
         end if
      end if

      call trc_dlfind_optimize(handle, int(natm), z_in, residues, xyz_in, opts, rb%res)

      if (rb%res%failed) then
         status = TRC_ERR_STATE
      else if (rb%res%converged) then
         status = TRC_OK
      else
         status = TRC_ERR_NOCONV
      end if
   end function trc_opt_run

   subroutine build_residues(natm, z, xyz, residues, err)
      !! --opt-coords hdlc/hdlc-tc's own connectivity partition -- the same
      !! bonded-component graph app/main.F90's own build_opt_residues
      !! builds (BOND_SCALE 1.3x the sum of covalent radii), duplicated
      !! here for the reason this module's own header gives for
      !! trc_opt_io.F90: this library cannot reach app/main.F90's module.
      integer, intent(in) :: natm
      integer, intent(in) :: z(natm)
      real(dp), intent(in) :: xyz(3, natm)
      integer, intent(out) :: residues(natm)
      type(error_t), intent(inout) :: err
      real(dp), parameter :: BOND_SCALE = 1.3_dp
      real(dp) :: radius(natm)
      integer :: i, j, comp, sp, a
      integer :: stack(natm)
      logical :: visited(natm)
      do i = 1, natm
         radius(i) = trc_covalent_radius(z(i), err)
         if (err%has_error()) return
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
   end subroutine build_residues

end module trc_opt_capi
