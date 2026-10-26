/*
 * libterco_opt's C ABI -- GENERATED, do not edit by hand.
 *
 * tools/gen_capi_headers.py parses src/opt/trc_opt_capi.F90's own
 * `bind(c, name=...)` procedures and public TRC_OPT_* parameters and
 * writes this file and include/trc_opt_c_interfaces.f90 from the same
 * pass, exactly as it does terco.h/trc_c_interfaces.f90 from
 * src/api/trc_capi.F90 -- see that pair for the full rationale.
 * Regenerate with that script; do not edit the output.
 *
 * Needs terco.h for `trc_handle_t` and the TRC_OK/TRC_ERR_* status
 * codes: this ABI defines neither again, it reuses terco's own (the
 * handle trc_opt_run below takes is the SAME context handle
 * trc_create returned). Only built/installed with
 * -DTERCO_ENABLE_DLFIND=ON (libterco_opt links libdlfind,
 * LGPL-3.0-or-later -- see the root CMakeLists.txt's option).
 */
#ifndef TERCO_OPT_H
#define TERCO_OPT_H

#include "terco.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Enum values (see each function's own doc for which one). */
#define TRC_OPT_COORDS_CARTESIAN 0
#define TRC_OPT_COORDS_HDLC 1
#define TRC_OPT_COORDS_HDLC_TC 2
#define TRC_OPT_COORDS_DLC 3
#define TRC_OPT_COORDS_DLC_TC 4
#define TRC_OPT_GUESS_DENSITY 0
#define TRC_OPT_GUESS_PROJECTED 1
#define TRC_OPT_GUESS_NONE 2
#define TRC_OPT_ALG_SD 0
#define TRC_OPT_ALG_CG_AUTO 1
#define TRC_OPT_ALG_CG 2
#define TRC_OPT_ALG_LBFGS 3
#define TRC_OPT_ALG_PRFO 10
#define TRC_OPT_ALG_NR 20
#define TRC_OPT_ALG_DAMPED 30
#define TRC_OPT_TARGET_MINIMUM 0
#define TRC_OPT_TARGET_SADDLE 1
#define TRC_OPT_SADDLE_PRFO 0
#define TRC_OPT_SADDLE_DIMER 1
#define TRC_OPT_HESS_TERCO 0
#define TRC_OPT_HESS_DLFIND 1

/**
 * Whether THIS build of libterco_opt can optimise a geometry at
 * all -- always true: with TERCO_ENABLE_DLFIND=OFF the library
 * itself is not built, so there is no trc_opt_available to call in
 * the first place. Kept for the same reason app/main.F90 still
 * calls trc_dlfind_bridge%trc_dlfind_available() itself rather than
 * assuming: a caller that detects this library at runtime (dlopen)
 * rather than at link time can ask instead of guessing.
 */
int trc_opt_available(void);

/**
 * A fresh options handle, every field at trc_opt_options_t's own
 * default (app/main.F90's unoptimised run is --opt-algorithm lbfgs,
 * --opt-coords cartesian, --opt-gtol 4.5e-4, ... -- trc_dlfind_bridge.F90's
 * own type definition).
 */
int trc_opt_options_create(trc_handle_t* handle);

int trc_opt_options_destroy(trc_handle_t handle);

/**
 * TRC_OPT_TARGET_MINIMUM or TRC_OPT_TARGET_SADDLE (--opt-target).
 */
int trc_opt_set_target(trc_handle_t handle, int target_kind);

/**
 * TRC_OPT_ALG_* (--opt-algorithm). Only meaningful for a minimum
 * search (--opt-target minimum); a saddle search picks its own
 * algorithm from --opt-saddle-method instead, same as the CLI.
 */
int trc_opt_set_algorithm(trc_handle_t handle, int algorithm);

/**
 * TRC_OPT_COORDS_* (--opt-coords).
 */
int trc_opt_set_coords(trc_handle_t handle, int coords_kind);

/**
 * TRC_OPT_SADDLE_PRFO or TRC_OPT_SADDLE_DIMER (--opt-saddle-method).
 */
int trc_opt_set_saddle_method(trc_handle_t handle, int method);

/**
 * --opt-max-steps: the cycle cap (DL-FIND's own `maxcycle`).
 */
int trc_opt_set_max_steps(trc_handle_t handle, int maxiter);

/**
 * --opt-gtol, Eh/Bohr: the max CARTESIAN COMPONENT of the gradient
 * at the last accepted step, same test trc_dlfind_bridge's own
 * convergence verdict uses (never DL-FIND's internal one -- that
 * file's own header).
 */
int trc_opt_set_gtol(trc_handle_t handle, double gtol);

/**
 * --opt-max-step, Bohr: the trust radius / max Cartesian step
 * DL-FIND's own `maxstep` enforces per cycle. A non-positive value
 * leaves DL-FIND's own default (0.5 Bohr, trc_opt_options_t's own
 * sentinel for "not set").
 */
int trc_opt_set_max_step(trc_handle_t handle, double max_step);

/**
 * TRC_OPT_HESS_TERCO or TRC_OPT_HESS_DLFIND (--opt-init-hessian),
 * for P-RFO/NR only -- ignored by every other algorithm, same as
 * the CLI.
 */
int trc_opt_set_hessian_source(trc_handle_t handle, int source);

/**
 * TRC_OPT_GUESS_* (--opt-guess): the warm start between steps.
 */
int trc_opt_set_guess(trc_handle_t handle, int guess_kind);

/**
 * --opt-frozen: `atoms` holds `n` 1-based indices into trc_opt_run's
 * own atom order. Replaces whatever a previous call set (not
 * additive): a caller that wants to grow the list reads it back
 * itself first -- there is no getter for it, cheap enough to track
 * on the caller's own side instead.
 */
int trc_opt_set_frozen_atoms(trc_handle_t handle, int n, const int* atoms);

int trc_opt_result_create(trc_handle_t* handle);

int trc_opt_result_destroy(trc_handle_t handle);

int trc_opt_result_converged(trc_handle_t handle, int* converged);

/**
 * Whether this run stopped on a recoverable failure (an SCF that
 * would not converge mid-optimisation) BEFORE DL-FIND's own
 * `error stop` could fire -- see this module's own header for why
 * that is not the only way an evaluation failure ends a run.
 */
int trc_opt_result_failed(trc_handle_t handle, int* failed);

int trc_opt_result_steps(trc_handle_t handle, int* steps);

int trc_opt_result_energy_evaluations(trc_handle_t handle, int* n);

int trc_opt_result_energy(trc_handle_t handle, double* energy);

int trc_opt_result_gradient_max(trc_handle_t handle, double* gmax);

int trc_opt_result_gradient_rms(trc_handle_t handle, double* grms);

/**
 * The last ACCEPTED geometry, atom-major Bohr (x1 y1 z1 x2 ...),
 * same convention as terco.h's trc_set_geometry/trc_gradient --
 * `xyz` is the caller's own buffer, length 3*natm.
 */
int trc_opt_result_geometry(trc_handle_t handle, int natm, double* xyz);

/**
 * The number of ACCEPTED-step frames available from
 * trc_opt_result_trajectory below (0: none, or the caller asked for
 * no trajectory to be kept -- there is no such switch in this v1
 * ABI, the trajectory trc_dlfind_optimize always collects is always
 * cheap enough to keep, same as app/main.F90's own --opt).
 */
int trc_opt_result_trajectory_count(trc_handle_t handle, int* n_traj);

/**
 * `xyz` is the caller's own buffer, length 3*natm*n_traj, frame
 * `k` (1-based) at offset 3*natm*(k-1); `energy` is length n_traj.
 * Both `natm` and `n_traj` must match trc_opt_run's own atom count
 * and trc_opt_result_trajectory_count's own answer exactly.
 */
int trc_opt_result_trajectory(trc_handle_t handle, int natm, int n_traj, double* xyz, double* energy);

/**
 * The run's own message (why it failed, or "" for an ordinary
 * finish), NUL-terminated into `buf` of `buflen` chars -- same
 * contract as terco.h's trc_message.
 */
int trc_opt_result_message(trc_handle_t handle, char* buf, int buflen);

/**
 * Optimise (or find a saddle on) the molecule already set up on
 * `handle` through terco.h (molecule, basis and method, matching
 * `z`/`xyz` -- same requirement as trc_dlfind_bridge%trc_dlfind_optimize,
 * which this calls directly). `z` (nuclear charges) and `xyz`
 * (atom-major Bohr, length 3*natm) are the caller's OWN copy of what
 * it already gave trc_set_molecule -- there is no getter for them on
 * the context handle (terco.h has none), so the caller hands them in
 * again rather than this ABI reading them back out of an opaque
 * object.
 *
 * `options` may be a null handle (every default); `result_handle`
 * must already exist (trc_opt_result_create) and is overwritten in
 * place. Returns TRC_OK when converged, TRC_ERR_NOCONV when the run
 * finished without converging (the usual "ran out of steps"),
 * TRC_ERR_STATE when an evaluation failed outright (the message is
 * on result_handle either way, trc_opt_result_message) -- but see
 * this module's own header: DL-FIND's `error stop` can still end the
 * calling PROCESS before this function ever returns at all.
 */
int trc_opt_run(trc_handle_t handle, int natm, const int* z, const double* xyz, trc_handle_t options, trc_handle_t result_handle);

#ifdef __cplusplus
}
#endif

#endif /* TERCO_OPT_H */
