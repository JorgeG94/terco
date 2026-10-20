/*
 * terco's C ABI -- GENERATED, do not edit by hand.
 *
 * tools/gen_capi_headers.py parses src/api/trc_capi.F90's own
 * `bind(c, name=...)` procedures and public TRC_* parameters and
 * writes this file and include/trc_c_interfaces.f90 from the same
 * pass, so the two cannot drift from each other or from the
 * library. Regenerate with that script; do not edit the output.
 */
#ifndef TERCO_H
#define TERCO_H

#ifdef __cplusplus
extern "C" {
#endif

/* An opaque handle to any object this API hands back -- a basis, a
 * shell-pair list, an ERI context, or a whole calculation context.
 * The Fortran side owns the memory; a caller only ever holds this
 * and passes it back. */
typedef void* trc_handle_t;

/* Status codes. Zero is success; everything else is a reason. */
#define TRC_OK 0
#define TRC_ERR_NULL 1
#define TRC_ERR_BADARG 2
#define TRC_ERR_UNSUPPORTED 3
#define TRC_ERR_NOCONV 4
#define TRC_ERR_STATE 5
#define TRC_GUESS_CORE 0
#define TRC_GUESS_GWH 1
#define TRC_GUESS_SAD 2
#define TRC_GUESS_GIVEN 3
#define TRC_GUESS_SADQ 4
#define TRC_GUESS_SAP 5
#define TRC_AO_NORM_SHELL 0
#define TRC_AO_NORM_UNIT 1
#define TRC_ABI_MAJOR 1
#define TRC_ABI_MINOR 0

int trc_basis_create(int nshell, int maxnp, const int* l, const int* nprim, const double* exps, const double* coefs, const double* centres, int natm, const double* zatm, const double* ratm, trc_handle_t* handle);

int trc_basis_destroy(trc_handle_t handle);

/**
 * Number of atomic orbitals, so a caller can size its own output buffers.
 */
int trc_basis_nao(trc_handle_t handle, int* nao);

int trc_pairs_create(trc_handle_t basis, double thresh, trc_handle_t* handle);

int trc_pairs_destroy(trc_handle_t handle);

int trc_pairs_count(trc_handle_t handle, int* npair);

/**
 * COVARIANT under `basis`'s own ao_norm_kind (set through
 * trc_set_ao_norm on a context that lent this handle via
 * trc_context_basis; a raw trc_basis_create* handle is always
 * TRC_AO_NORM_SHELL and this is then a no-op): S_out = N S N, same
 * for T and V.
 */
int trc_compute_1e(trc_handle_t basis, trc_handle_t pairs, double* smat, double* tmat, double* vmat);

int trc_compute_df2c(trc_handle_t aux, double* jmat);

int trc_compute_df3c(trc_handle_t basis, trc_handle_t pairs, trc_handle_t aux, double* tens);

int trc_basis_create_libcint(const int* atm, int natm, const int* bas, int nshell, const double* env, int nenv, int cartesian, trc_handle_t* handle);

/**
 * Highest angular momentum in the basis, so a caller can decide whether to
 * use this library at all before it builds anything expensive. The
 * four-centre kernels stop at d and the caller has to fall back rather than
 * get a wrong answer.
 */
int trc_basis_maxl(trc_handle_t handle, int* maxl);

/**
 * The angular-momentum ceilings THIS build was compiled for, so a host
 * can check a basis against them instead of hard-coding them: the
 * four-centre ERIs / SCF (TERCO_LMAX: 2 by default, 3 for f), and the
 * density-fitting orbital and auxiliary bases (TERCO_DF_LMAX; the
 * auxiliary one is one higher). Compare with trc_basis_maxl.
 */
int trc_max_angular_momentum(int* eri_lmax, int* df_lmax, int* df_lmax_aux);

/**
 * Number of Cartesian multipole components, so a caller can size its buffer
 * without hardcoding 39 and drifting when the operator order changes.
 * Named for what it counts rather than `trc_nmult`, which collides with the
 * parameter of that name it returns.
 */
int trc_multipole_count(void);

/**
 * COVARIANT under `basis`'s own ao_norm_kind, one component at a
 * time (see trc_compute_1e).
 */
int trc_compute_multipoles(trc_handle_t basis, trc_handle_t pairs, const double* origin, double* mmat);

int trc_eri_create(trc_handle_t basis, double thresh, trc_handle_t* handle);

int trc_eri_destroy(trc_handle_t handle);

/**
 * CONTRAVARIANT in (dmat), COVARIANT out (gmat), under `basis`'s own
 * ao_norm_kind: this is GAMESS's own main caller (its Fock/J/K
 * build, run_Fockbuild's role in hpcchem_driver.src), so it is wired
 * rather than left at TRC_AO_NORM_SHELL like the entries the "not
 * wired" list below still names. See fock_scale_in/fock_scale_out
 * for why SHELL mode still costs nothing extra.
 */
int trc_fock(trc_handle_t eri, trc_handle_t basis, const double* dmat, double* gmat, double jfac, double kfac, int dscreen, double omega);

/**
 * CONTRAVARIANT in, COVARIANT out, under `basis`'s own ao_norm_kind
 * -- see capi_fock; the density-index-fastest layout means the
 * UNIT-mode scaling here is its own small routines
 * (fock_many_scale_in/out) rather than a reuse of fock_scale_in/out.
 */
int trc_fock_many(trc_handle_t eri, trc_handle_t basis, int ndens, const double* dmats, double* gmats, double jfac, double kfac, int dscreen, double omega);

/**
 * CONTRAVARIANT in, COVARIANT out, exactly as capi_fock -- the
 * diagonal congruence ao_scale_cov/fock_scale_in apply does not
 * care whether `dmat` is symmetric, so nothing here differs from
 * the symmetric entry beyond which shared routines it reuses.
 */
int trc_fock_nosym(trc_handle_t eri, trc_handle_t basis, const double* dmat, double* gmat, int dscreen);

int trc_scf(trc_handle_t basis, int nalpha, int nbeta, const char* functional, int grid_level, double conv_energy, double conv_density, int max_iter, trc_handle_t dguess, int verbose, double* energy, double* e_xc, double* dmat, double* eps, int* niter);

int trc_scf_mpi(int fcomm, trc_handle_t basis, int nalpha, int nbeta, const char* functional, int grid_level, double conv_energy, double conv_density, int max_iter, trc_handle_t dguess, int verbose, double* energy, double* e_xc, double* dmat, double* eps, int* niter);

int trc_rimp2(int fcomm, trc_handle_t basis, trc_handle_t aux, int nfrozen, int aux_block, double* e_os, double* e_ss);

int trc_create(trc_handle_t* handle);

int trc_destroy(trc_handle_t handle);

/**
 * Atoms in Bohr, nuclear charges as doubles, atom-major coordinates.
 */
int trc_set_molecule(trc_handle_t handle, int natm, const double* z, const double* xyz, int charge, int multiplicity);

/**
 * The orbital basis from libcint's packed `atm`/`bas`/`env`.
 * `cartesian` nonzero builds terco's own Cartesian AO space; zero
 * tags it spherical (libcint/pyscf order, m = -l..l except p as x,
 * y, z) -- trc_nao and every SCF-facing entry then size and index
 * at nao_sph, exactly as for a spherical trc_set_basis_json. A
 * mismatch against an already-set auxiliary or jkfit basis's own
 * convention is refused, not coerced.
 */
int trc_set_basis_libcint(trc_handle_t handle, const int* atm, int natm, const int* bas, int nshell, const double* env, int nenv, int cartesian);

/**
 * The auxiliary (RI-MP2 rifit) basis, same layout and `cartesian`
 * convention as trc_set_basis_libcint; refused rather than coerced
 * if it disagrees with the orbital basis already set.
 */
int trc_set_aux_libcint(trc_handle_t handle, const int* atm, int natm, const int* bas, int nshell, const double* env, int nenv, int cartesian);

int trc_set_basis_arrays(trc_handle_t handle, int nshell, int maxnp, const int* l, const int* nprim, const double* exps, const double* coefs, const double* centres);

int trc_set_aux_arrays(trc_handle_t handle, int nshell, int maxnp, const int* l, const int* nprim, const double* exps, const double* coefs, const double* centres);

int trc_set_basis_json(trc_handle_t handle, const char* path);

int trc_set_aux_json(trc_handle_t handle, const char* path);

/**
 * An empty functional name is Hartree-Fock.
 */
int trc_set_method(trc_handle_t handle, const char* functional, int grid_level);

/**
 * Grid level for VV10's nonlocal correlation (default 1, coarser than
 * the semilocal grid: VV10 is a double sum over its points). Only
 * used by a functional with VV10, such as wB97X-V.
 */
int trc_set_nlc_grid(trc_handle_t handle, int grid_level);

int trc_set_pcm(trc_handle_t handle, const char* model, double eps, int lebedev, double radii_scale);

/**
 * The continuum's part of the last SCF: E_diel = q.v/2 (already in
 * trc_energy's total), the total apparent surface charge, and the
 * surface points the cavity kept. Zeros when no continuum was set.
 */
int trc_pcm_energy(trc_handle_t handle, double* e_diel, double* q_total, int* npoints);

int trc_set_convergence(trc_handle_t handle, double conv_energy, double conv_diis, int max_iter, int ndiis);

int trc_set_screening(trc_handle_t handle, double thresh);

int trc_set_incremental(trc_handle_t handle, int on);

int trc_set_df(trc_handle_t handle, int on);

int trc_set_jkfit_json(trc_handle_t handle, const char* path);

/**
 * The --df fitting basis from libcint's packed arrays, the same
 * shape and `cartesian` convention as trc_set_aux_libcint; refused
 * rather than coerced if it disagrees with the orbital basis
 * already set (pyscf's auxmol.cart = mol.cart rule).
 */
int trc_set_jkfit_libcint(trc_handle_t handle, const int* atm, int natm, const int* bas, int nshell, const double* env, int nenv, int cartesian);

int trc_set_batching(trc_handle_t handle, int res);

/**
 * Force the Cartesian AO space on every basis this context loads
 * from here on, regardless of what a JSON file's function_type says
 * (trc_scf_options_t%force_cartesian). `flag` non-zero sets it.
 *
 * NOT RETROACTIVE: a basis already loaded (trc_set_basis_json) kept
 * whichever space it was built in at that call -- the same way
 * trc_set_ao_norm's own comment describes for a TRC_GUESS_GIVEN
 * density already copied. Call this BEFORE trc_set_basis_json/
 * trc_set_aux_json for it to take effect on that basis; a caller
 * wanting it to also apply to a basis already loaded has to load
 * the file again.
 */
int trc_set_force_cartesian(trc_handle_t handle, int flag);

/**
 * `kind` one of TRC_GUESS_*. For TRC_GUESS_GIVEN, `dguess` is
 * (nao, nao, nspin) column-major and is COPIED here; the caller's array
 * is free afterwards. The basis must be set first for that copy.
 *
 * CONTRAVARIANT under trc_set_ao_norm: a density given in
 * TRC_AO_NORM_UNIT is converted to terco's own shell-normalised
 * convention here, once, at the copy -- D_terco = N D_given N (see
 * trc_ao_norm.F90's header; this is the "guess coming IN" direction,
 * the covariant congruence run forward). Everything downstream of
 * this copy, including trc_run_scf's own iteration, only ever sees
 * terco's native convention.
 */
int trc_set_guess(trc_handle_t handle, int kind, trc_handle_t dguess, int nspin);

/**
 * How this context's AO-indexed C API entries normalise a Cartesian
 * shell's components: TRC_AO_NORM_SHELL (0, the default -- terco's
 * own convention, and libcint's) or TRC_AO_NORM_UNIT (1 -- every
 * component individually normalised to 1, GAMESS's convention; see
 * trc_ao_norm.F90's header for the closed form and which entries
 * convert). Defaulting to SHELL means no caller that never calls
 * this sees any difference from before this entry existed.
 *
 * Order-independent of trc_set_basis* and trc_set_aux*: settable before
 * or after either, and takes effect on whichever basis is present now
 * and on any set afterwards. NOT retroactive on an already-copied
 * TRC_GUESS_GIVEN density, though -- that conversion happened once,
 * at trc_set_guess's own call, under whatever mode was active then;
 * call this before trc_set_guess if the two must agree.
 */
int trc_set_ao_norm(trc_handle_t handle, int kind);

/**
 * The AO normalisation this context is set to -- TRC_AO_NORM_SHELL
 * unless trc_set_ao_norm was called with TRC_AO_NORM_UNIT.
 */
int trc_get_ao_norm(trc_handle_t handle, int* kind);

/**
 * `fcomm` is the Fortran handle of ANY intra-communicator the caller
 * already has open and collective over -- MPI_Comm_c2f from C, or
 * the integer a Fortran caller holds (mpi_f08's MPI_VAL, or the
 * legacy `mpi` module's handle directly); -1 means one rank, no MPI.
 * It need not be MPI_COMM_WORLD: GAMESS, for example, hands this
 * DDI's compute communicator, whose data-server ranks are not part
 * of it and never call trc_run_scf/trc_run_rimp2 at all. terco
 * duplicates whatever `fcomm` names (MPI_Comm_dup) the first time a
 * run needs it, so every rank of THAT communicator must reach the
 * run together; the caller keeps ownership of its own communicator,
 * since only terco's duplicate is ever freed.
 */
int trc_set_comm(trc_handle_t handle, int fcomm);

int trc_set_verbose(trc_handle_t handle, int level);

/**
 * `nfrozen` -1 counts the core from the elements; `aux_block` <= 0 is
 * the whole auxiliary basis at once.
 */
int trc_set_rimp2(trc_handle_t handle, int nfrozen, int aux_block);

/**
 * Needs the molecule and the basis. Builds the guess it was asked for,
 * runs, keeps everything in the context.
 */
int trc_run_scf(trc_handle_t handle);

/**
 * ======================================================================
 * !!! EXPERIMENTAL !!! Validated against pyscf and finite differences,
 * but NOT performance-tuned: the --df gradient can cost far more than
 * the energy, and only one rank is supported. Do not rely on its cost
 * yet.
 * ======================================================================
 *
 * The full analytic RHF/UHF nuclear gradient at the converged SCF
 * this context already holds (trc_scf_grad%trc_scf_gradient: the
 * one-electron terms, eri%gradient's two-electron part, and nuclear
 * repulsion -- see that module's own header for the assembly order
 * and the factor of 1/2). When the orbital basis is spherical,
 * `trc_scf_gradient` itself folds D and W to the Cartesian space the
 * derivative kernels need (T D T^T) before touching an integral;
 * nothing here has to know which case it is in.
 *
 * Reuses `cx%la` (initialised once per context by `trc_run_scf`,
 * released only at `trc_destroy`) rather than a `trc_linalg_t` of
 * its own: a second handle here would need its own release before
 * going out of scope (trc_linalg_t has no finalizer), for no reason
 * -- the gradient's GEMMs are the same size `cx%la` was already
 * sized for.
 *
 * NOT refused any more: a functional with VV10 nonlocal correlation
 * (wB97X-V). `trc_scf_gradient` adds trc_vv10_grad%vv10_gradient_add's
 * term on top of the semilocal one, on VV10's own grid
 * (`cx%opts%nlc_grid_level`, the SAME level the SCF itself converged
 * it on). wB97M-V is meta-GGA and cannot reach this point at all (see
 * below), so it is not a case this routine has to decide about.
 *
 * REFUSED, clearly, rather than attempted:
 *   - meta-GGA or a double hybrid: `trc_xc_functional_t` has no
 *     meta-GGA family and no MP2-fraction field on this branch, so
 *     neither can reach this point at all -- `trc_set_method` already
 *     refused it (func%missing) before a functional name naming one
 *     could be stored in `cx%opts%functional`. Nothing to check here.
 *   - RI-MP2 (or, once it exists, a double hybrid): no correlation
 *     gradient exists; `trc_scf_gradient` is the SCF reference term
 *     only.
 *   - --df outside the DF-SCF's own scope (unrestricted, spherical,
 *     range-separated -- the SCF refuses those first), or with the fit
 *     gone. A --df gradient otherwise runs the FITTED J/K gradient
 *     (trc_scf_grad_df.F90) over `cx%dfjk`, the fit the SCF used.
 *   - PCM: no cavity/surface-charge derivative exists.
 *   - a multi-rank context (`cx%fcomm /= -1`). `trc_scf_gradient`
 *     takes no communicator: it rebuilds the pair list and ERI
 *     engine from `cx%bb%b` alone, so on a collective run it would
 *     rebuild the FULL (unsplit) pair list, Schwarz screen and ERI
 *     structures on EVERY rank's own device -- not merely redundant
 *     work, but re-materialising, on one GPU, exactly the memory a
 *     multi-rank SCF split across ranks to avoid. `eri%gradient`
 *     itself DOES support a rank split (`this%rank`/`this%nranks`,
 *     `sum_ranks`'s allreduce) but only when its own `eri%build` is
 *     given a `comm` -- a plumbing `trc_scf_gradient` does not do
 *     today, found by reading trc_eri.F90's `eri_gradient` and
 *     comparing it against `eri_build`'s own comm handling.
 *
 * A plain LDA/GGA functional, a GLOBAL hybrid of one (b3lyp and
 * friends), or a RANGE-SEPARATED hybrid (wb97x and friends: nonzero
 * `func%omega`/`exx_lr`, `nlc_b`/`nlc_c` still zero) is NOT refused:
 * `func` is looked up again here (same `xc_functional_by_name`
 * trc_set_method already validated the name with) and handed to
 * `trc_scf_gradient`, along with the SAME grid settings
 * (`cx%opts%grid_level`, `grid_max_pts`, `rho_tol`) the SCF itself
 * converged with -- the gradient's quadrature has to be the SCF's,
 * or the two are not derivatives of the same energy.
 * `trc_scf_gradient` itself builds the erf(omega r)/r long-range
 * exchange term (a second eri%gradient pass, its own header) and
 * trc_xc_grad%xc_gradient_add computes the functional's semilocal
 * GGA part exactly as for any other GGA -- see that module's own
 * header on why `omega`/`exx_lr` do not affect it.
 *
 * GRID REUSE. `cx%kgrid`, the SCF's own grid(s) kept alive since
 * `trc_run_scf` (`context_box%kgrid`'s own header), is handed
 * straight to `trc_scf_gradient` when it is built AND matches
 * every setting the gradient's quadrature has to match the SCF's
 * on (`reuse_grid`/`reuse_nlc`, below) -- the same list this
 * routine's own header already promises (grid level, max points
 * per batch, rho_tol, nlc grid level); geometry and basis are
 * covered for free, because ANY change to either would have
 * cleared `cx%scf_ok` (`invalidate_runs`) before this point could
 * ever be reached. Otherwise `trc_scf_gradient` rebuilds exactly as
 * it always has -- this is a pure performance path, never a
 * correctness one: `grid`/`xg` are read-only there, so a mismatch
 * caught here is the only thing standing between "reuse" and
 * "quietly wrong".
 */
int trc_run_gradient(trc_handle_t handle);

/**
 * (3, natm), Eh/Bohr, column-major (x, y, z fastest, then atom, the
 * same atom order `trc_set_molecule` was given) -- after a
 * successful `trc_run_gradient`. NOT AO-indexed data (see
 * trc_capi.F90's own "AO NORMALISATION" header): a nuclear gradient
 * is a derivative of the SCF energy, which does not depend on which
 * AO normalisation convention computed it, so `trc_set_ao_norm`
 * plays no role here.
 */
int trc_gradient(trc_handle_t handle, double* grad);

/**
 * Needs a converged SCF in the context (restricted, or unrestricted --
 * UHF or UKS, cx%res%nspin == 2) and an auxiliary basis. An
 * unrestricted reference refuses more than one rank below (see
 * trc_urimp2_run's own header for why the restricted ring is not a
 * clean generalisation here). The orbital basis may be spherical:
 * `cx%res%cmo`/`eps` (at whatever space the SCF ran in) are folded
 * into the Cartesian, zero-padded space `trc_rimp2_run`/
 * `trc_urimp2_run` themselves run in by `rimp2_mo_to_cart`
 * (trc_rimp2_driver.F90's own header says why the pad is exact).
 */
int trc_run_rimp2(trc_handle_t handle);

/**
 * Write the context's SCF result -- orbitals, orbital energies,
 * density, energies -- to `path` (NUL-terminated), with the molecule,
 * basis and method it belongs to; trc_load_scf reads it back. One
 * rank is enough: every rank holds the same result.
 */
int trc_save_scf(trc_handle_t handle, const char* path);

/**
 * In place of trc_run_scf: take the SCF result from a file
 * trc_save_scf wrote, so a correlated step (trc_run_rimp2) can follow
 * without repeating the SCF. The molecule, basis, method and solvent
 * must be set as for trc_run_scf, and must be the ones the file was
 * made with; anything else is refused, naming the difference. Every
 * rank reads the file itself. TRC_ERR_NOCONV if the saved SCF had not
 * converged, as trc_run_scf would have said.
 */
int trc_load_scf(trc_handle_t handle, const char* path);

/**
 * 0 (lexicographic) or 1 (gamess): the within-shell Cartesian order every
 * AO-indexed array this API exchanges uses, fixed at compile time
 * (TRC_CART_ORDER) and the same for every handle -- no handle argument,
 * and it always succeeds.
 */
int trc_cart_order(void);

/**
 * This build's ABI version (TRC_ABI_MAJOR, TRC_ABI_MINOR above), so
 * a host can compare it against the constants the header it
 * compiled against carries before trusting any other entry point.
 * Bump TRC_ABI_MINOR for an additive change (a new entry point, a
 * new enum value appended at the end); bump TRC_ABI_MAJOR, and
 * reset TRC_ABI_MINOR to 0, for anything that changes an existing
 * signature or an existing enum value's meaning. No handle; always
 * succeeds.
 */
int trc_abi_version(int* major, int* minor);

/**
 * Look up a functional by name (case-insensitive, trc_xc_functional's
 * own table) without a context: `flags` is a bitmask --
 *
 *   1  known      -- runnable today (trc_set_method would accept it)
 *   2  hybrid     -- exx or exx_lr nonzero (some exact exchange)
 *   4  range-separated -- omega nonzero
 *   8  VV10       -- nlc_b nonzero (nonlocal correlation)
 *   16 meta-GGA   -- needs a kinetic-energy-density kernel; always
 *      paired with bit 64 today, since terco has none translated yet
 *   32 double hybrid -- xc_is_double_hybrid (pt2_cos or pt2_css nonzero)
 *   64 listed but not runnable -- trc_xc_functional_t%missing is set;
 *      bit 1 is then NOT set, and `exx`/`c_os`/`c_ss`/the other bits
 *      reflect only whatever the table's own entry set before
 *      giving up (most "missing" entries set nothing at all, so
 *      every other bit and output reads zero -- that is the table's
 *      gap, not this entry's)
 *
 * `exx` is the fraction of ordinary (full-range) exact exchange;
 * `c_os`/`c_ss` are the opposite-/same-spin RI-MP2 weights a double
 * hybrid adds on top (both 0 for anything that is not one). An
 * unknown name returns a non-OK status with flags = 0 and every
 * output zeroed, rather than TRC_OK with nothing found.
 */
int trc_functional_info(const char* name, int* flags, double* exx, double* c_os, double* c_ss);

/**
 * XC_AVAILABLE_FUNCTIONALS -- every name trc_set_method/
 * trc_functional_info can look up that terco can actually run today
 * -- space-separated and NUL-terminated into `buf` (capacity
 * `buflen` bytes, NUL included). `reqlen` is ALWAYS set to the
 * number of bytes that call needed, NUL included, whether or not it
 * fit: a caller queries it with `buflen = 0` (`buf` is never
 * touched then) and allocates exactly `reqlen` bytes for a second,
 * filling call. Returns TRC_ERR_BADARG, leaving `buf` untouched, if
 * `buflen` is positive but smaller than `reqlen`.
 */
int trc_functional_list(char* buf, int buflen, int* reqlen);

int trc_nao(trc_handle_t handle, int* nao);

/**
 * Non-zero if the orbital basis this context holds is running in
 * the spherical AO space (nao_sph of trc_nao, rather than the
 * Cartesian nao) -- the basis's own JSON function_type unless
 * force_cartesian overrode it. For a caller checking which space
 * `trc_nao`'s number, and every AO-indexed array it sizes, is in.
 */
int trc_is_spherical(trc_handle_t handle, int* flag);

int trc_nspin(trc_handle_t handle, int* nspin);

/**
 * The SCF total, plus the RI-MP2 correlation if that ran: the FULL
 * E_os + E_ss for a plain post-SCF RI-MP2 (any functional, or
 * Hartree-Fock), or c_os E_os + c_ss E_ss for a double hybrid --
 * trc_rimp2_energy below always reports the raw, unscaled E_os/E_ss
 * either way, so a caller after the PT2 part alone never sees this
 * weighting. A caller asking for a double hybrid's total needs
 * nothing beyond the ordinary trc_run_scf/trc_run_rimp2/trc_energy
 * sequence any RI-MP2 run already uses -- the weighting is applied
 * here, once, from the functional's own pt2_cos/pt2_css.
 */
int trc_energy(trc_handle_t handle, double* energy);

int trc_energy_parts(trc_handle_t handle, double* e_nuc, double* e_one, double* e_two, double* e_xc);

int trc_converged(trc_handle_t handle, int* flag);

int trc_iterations(trc_handle_t handle, int* niter);

/**
 * (nao, nao, nspin), column-major; restricted holds the total.
 *
 * CONTRAVARIANT under trc_set_ao_norm: D_out = N^-1 D_terco N^-1 per
 * spin slice when UNIT. Scaled on a COPY -- `cx%res%dmat` itself
 * stays in terco's own convention, because trc_run_rimp2 and a
 * second trc_run_scf both read it as such internally.
 */
int trc_density(trc_handle_t handle, double* dmat);

/**
 * CONTRAVARIANT under trc_set_ao_norm, on the AO (row) index only:
 * C_out = N^-1 C_terco per spin slice when UNIT -- see trc_density
 * for why this scales a copy, not `cx%res%cmo` itself.
 */
int trc_mo_coefficients(trc_handle_t handle, double* cmo);

int trc_mo_energies(trc_handle_t handle, double* eps);

int trc_rimp2_energy(trc_handle_t handle, double* e_os, double* e_ss);

/**
 * The last message, NUL-terminated into `buf` of `buflen` chars.
 */
int trc_message(trc_handle_t handle, char* buf, int buflen);

/**
 * The basis handle inside the context: BORROWED, never destroyed by
 * the caller, valid until the context is destroyed or the basis reset.
 */
int trc_context_basis(trc_handle_t handle, trc_handle_t* basis);

/**
 * The ERI context at the context's screening threshold, built on first
 * use. BORROWED, as trc_context_basis.
 */
int trc_context_eri(trc_handle_t handle, trc_handle_t* eri);

#ifdef __cplusplus
}
#endif

#endif /* TERCO_H */
