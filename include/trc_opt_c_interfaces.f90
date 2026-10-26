!
! The authoritative Fortran declaration of libterco_opt's C ABI --
! GENERATED, do not edit by hand.
!
! tools/gen_capi_headers.py parses src/opt/trc_opt_capi.F90's own
! `bind(c, name=...)` procedures and public TRC_OPT_* parameters and
! writes this file and include/terco_opt.h from the same pass -- see
! that pair's own header for the full rationale. Regenerate with that
! script; do not edit the output.
!
! Only built/installed with -DTERCO_ENABLE_DLFIND=ON.
!
module trc_opt_c_interfaces
   use, intrinsic :: iso_c_binding
   implicit none
   public

   integer(c_int), parameter :: TRC_OPT_COORDS_CARTESIAN = 0
   integer(c_int), parameter :: TRC_OPT_COORDS_HDLC = 1
   integer(c_int), parameter :: TRC_OPT_COORDS_HDLC_TC = 2
   integer(c_int), parameter :: TRC_OPT_COORDS_DLC = 3
   integer(c_int), parameter :: TRC_OPT_COORDS_DLC_TC = 4
   integer(c_int), parameter :: TRC_OPT_GUESS_DENSITY = 0
   integer(c_int), parameter :: TRC_OPT_GUESS_PROJECTED = 1
   integer(c_int), parameter :: TRC_OPT_GUESS_NONE = 2
   integer(c_int), parameter :: TRC_OPT_ALG_SD = 0
   integer(c_int), parameter :: TRC_OPT_ALG_CG_AUTO = 1
   integer(c_int), parameter :: TRC_OPT_ALG_CG = 2
   integer(c_int), parameter :: TRC_OPT_ALG_LBFGS = 3
   integer(c_int), parameter :: TRC_OPT_ALG_PRFO = 10
   integer(c_int), parameter :: TRC_OPT_ALG_NR = 20
   integer(c_int), parameter :: TRC_OPT_ALG_DAMPED = 30
   integer(c_int), parameter :: TRC_OPT_TARGET_MINIMUM = 0
   integer(c_int), parameter :: TRC_OPT_TARGET_SADDLE = 1
   integer(c_int), parameter :: TRC_OPT_SADDLE_PRFO = 0
   integer(c_int), parameter :: TRC_OPT_SADDLE_DIMER = 1
   integer(c_int), parameter :: TRC_OPT_HESS_TERCO = 0
   integer(c_int), parameter :: TRC_OPT_HESS_DLFIND = 1

   interface

      integer(c_int) function trc_opt_available() bind(c) result(res)
         !! Whether THIS build of libterco_opt can optimise a geometry at
         !! all -- always true: with TERCO_ENABLE_DLFIND=OFF the library
         !! itself is not built, so there is no trc_opt_available to call in
         !! the first place. Kept for the same reason app/main.F90 still
         !! calls trc_dlfind_bridge%trc_dlfind_available() itself rather than
         !! assuming: a caller that detects this library at runtime (dlopen)
         !! rather than at link time can ask instead of guessing.
         import :: c_int
         implicit none
      end function trc_opt_available

      integer(c_int) function trc_opt_options_create(handle) bind(c) result(res)
         !! A fresh options handle, every field at trc_opt_options_t's own
         !! default (app/main.F90's unoptimised run is --opt-algorithm lbfgs,
         !! --opt-coords cartesian, --opt-gtol 4.5e-4, ... -- trc_dlfind_bridge.F90's
         !! own type definition).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(out) :: handle
      end function trc_opt_options_create

      integer(c_int) function trc_opt_options_destroy(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_opt_options_destroy

      integer(c_int) function trc_opt_set_target(handle, target_kind) bind(c) result(res)
         !! TRC_OPT_TARGET_MINIMUM or TRC_OPT_TARGET_SADDLE (--opt-target).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: target_kind
      end function trc_opt_set_target

      integer(c_int) function trc_opt_set_algorithm(handle, algorithm) bind(c) result(res)
         !! TRC_OPT_ALG_* (--opt-algorithm). Only meaningful for a minimum
         !! search (--opt-target minimum); a saddle search picks its own
         !! algorithm from --opt-saddle-method instead, same as the CLI.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: algorithm
      end function trc_opt_set_algorithm

      integer(c_int) function trc_opt_set_coords(handle, coords_kind) bind(c) result(res)
         !! TRC_OPT_COORDS_* (--opt-coords).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: coords_kind
      end function trc_opt_set_coords

      integer(c_int) function trc_opt_set_saddle_method(handle, method) bind(c) result(res)
         !! TRC_OPT_SADDLE_PRFO or TRC_OPT_SADDLE_DIMER (--opt-saddle-method).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: method
      end function trc_opt_set_saddle_method

      integer(c_int) function trc_opt_set_max_steps(handle, maxiter) bind(c) result(res)
         !! --opt-max-steps: the cycle cap (DL-FIND's own `maxcycle`).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: maxiter
      end function trc_opt_set_max_steps

      integer(c_int) function trc_opt_set_gtol(handle, gtol) bind(c) result(res)
         !! --opt-gtol, Eh/Bohr: the max CARTESIAN COMPONENT of the gradient
         !! at the last accepted step, same test trc_dlfind_bridge's own
         !! convergence verdict uses (never DL-FIND's internal one -- that
         !! file's own header).
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), value :: gtol
      end function trc_opt_set_gtol

      integer(c_int) function trc_opt_set_max_step(handle, max_step) bind(c) result(res)
         !! --opt-max-step, Bohr: the trust radius / max Cartesian step
         !! DL-FIND's own `maxstep` enforces per cycle. A non-positive value
         !! leaves DL-FIND's own default (0.5 Bohr, trc_opt_options_t's own
         !! sentinel for "not set").
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), value :: max_step
      end function trc_opt_set_max_step

      integer(c_int) function trc_opt_set_hessian_source(handle, source) bind(c) result(res)
         !! TRC_OPT_HESS_TERCO or TRC_OPT_HESS_DLFIND (--opt-init-hessian),
         !! for P-RFO/NR only -- ignored by every other algorithm, same as
         !! the CLI.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: source
      end function trc_opt_set_hessian_source

      integer(c_int) function trc_opt_set_guess(handle, guess_kind) bind(c) result(res)
         !! TRC_OPT_GUESS_* (--opt-guess): the warm start between steps.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: guess_kind
      end function trc_opt_set_guess

      integer(c_int) function trc_opt_set_frozen_atoms(handle, n, atoms) bind(c) result(res)
         !! --opt-frozen: `atoms` holds `n` 1-based indices into trc_opt_run's
         !! own atom order. Replaces whatever a previous call set (not
         !! additive): a caller that wants to grow the list reads it back
         !! itself first -- there is no getter for it, cheap enough to track
         !! on the caller's own side instead.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         integer(c_int), intent(in) :: atoms(*)
      end function trc_opt_set_frozen_atoms

      integer(c_int) function trc_opt_result_create(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(out) :: handle
      end function trc_opt_result_create

      integer(c_int) function trc_opt_result_destroy(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_opt_result_destroy

      integer(c_int) function trc_opt_result_converged(handle, converged) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: converged
      end function trc_opt_result_converged

      integer(c_int) function trc_opt_result_failed(handle, failed) bind(c) result(res)
         !! Whether this run stopped on a recoverable failure (an SCF that
         !! would not converge mid-optimisation) BEFORE DL-FIND's own
         !! `error stop` could fire -- see this module's own header for why
         !! that is not the only way an evaluation failure ends a run.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: failed
      end function trc_opt_result_failed

      integer(c_int) function trc_opt_result_steps(handle, steps) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: steps
      end function trc_opt_result_steps

      integer(c_int) function trc_opt_result_energy_evaluations(handle, n) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: n
      end function trc_opt_result_energy_evaluations

      integer(c_int) function trc_opt_result_energy(handle, energy) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: energy
      end function trc_opt_result_energy

      integer(c_int) function trc_opt_result_gradient_max(handle, gmax) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: gmax
      end function trc_opt_result_gradient_max

      integer(c_int) function trc_opt_result_gradient_rms(handle, grms) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: grms
      end function trc_opt_result_gradient_rms

      integer(c_int) function trc_opt_result_geometry(handle, natm, xyz) bind(c) result(res)
         !! The last ACCEPTED geometry, atom-major Bohr (x1 y1 z1 x2 ...),
         !! same convention as terco.h's trc_set_geometry/trc_gradient --
         !! `xyz` is the caller's own buffer, length 3*natm.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natm
         real(c_double), intent(out) :: xyz(*)
      end function trc_opt_result_geometry

      integer(c_int) function trc_opt_result_trajectory_count(handle, n_traj) bind(c) result(res)
         !! The number of ACCEPTED-step frames available from
         !! trc_opt_result_trajectory below (0: none, or the caller asked for
         !! no trajectory to be kept -- there is no such switch in this v1
         !! ABI, the trajectory trc_dlfind_optimize always collects is always
         !! cheap enough to keep, same as app/main.F90's own --opt).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: n_traj
      end function trc_opt_result_trajectory_count

      integer(c_int) function trc_opt_result_trajectory(handle, natm, n_traj, xyz, energy) bind(c) result(res)
         !! `xyz` is the caller's own buffer, length 3*natm*n_traj, frame
         !! `k` (1-based) at offset 3*natm*(k-1); `energy` is length n_traj.
         !! Both `natm` and `n_traj` must match trc_opt_run's own atom count
         !! and trc_opt_result_trajectory_count's own answer exactly.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natm
         integer(c_int), value :: n_traj
         real(c_double), intent(out) :: xyz(*)
         real(c_double), intent(out) :: energy(*)
      end function trc_opt_result_trajectory

      integer(c_int) function trc_opt_result_message(handle, buf, buflen) bind(c) result(res)
         !! The run's own message (why it failed, or "" for an ordinary
         !! finish), NUL-terminated into `buf` of `buflen` chars -- same
         !! contract as terco.h's trc_message.
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(out) :: buf(*)
         integer(c_int), value :: buflen
      end function trc_opt_result_message

      integer(c_int) function trc_opt_run(handle, natm, z, xyz, options, result_handle) bind(c) result(res)
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
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natm
         integer(c_int), intent(in) :: z(*)
         real(c_double), intent(in) :: xyz(*)
         type(c_ptr), value :: options
         type(c_ptr), value :: result_handle
      end function trc_opt_run

   end interface

end module trc_opt_c_interfaces
