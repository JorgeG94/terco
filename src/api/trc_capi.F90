!
! C ABI.
!
! The Fortran API is the real one; this is a thin `bind(c)` skin over it so a
! caller that is not Fortran can use the library. It owns no logic -- every
! entry point unpacks C types, calls the Fortran routine, and returns a status
! code.
!
! HANDLES
! -------
! Objects cross the boundary as opaque `void*`. Each handle is a `c_ptr` to a
! heap-allocated Fortran derived type; the caller never sees the layout and
! never allocates it. That is the same contract cuEST offers and it is what
! lets the Fortran side change shape without breaking callers.
!
! WHAT IS DELIBERATELY ABSENT
! ---------------------------
! No workspace-query protocol. cuEST has one (`*CreateWorkspaceQuery` fills a
! descriptor, the caller allocates persistent and temporary buffers and hands
! them back) because a C library cannot own device memory on the caller's
! behalf without an allocator contract. We are Fortran underneath with
! allocatable components and `release`, so the objects own their own memory.
! If a caller ever needs to place that memory itself, the protocol belongs
! here, at this boundary, and not in the core.
!
! ERRORS
! ------
! Status codes, not exceptions or aborts. A library that stops the process on
! bad input is unusable inside someone else's SCF.
!
! COMMUNICATORS
! -------------
! Every entry that takes an `fcomm` (trc_scf_mpi, trc_rimp2, trc_set_comm)
! accepts the Fortran handle of ANY intra-communicator the caller already
! has open and is willing to be collective over -- not just MPI_COMM_WORLD.
! -1 means one rank and no MPI, in any build. terco duplicates whatever the
! handle names (MPI_Comm_dup, so every rank of it must call together) into
! its own comm_t and frees only that duplicate afterwards; the caller's own
! communicator is untouched and stays open for its own use throughout. This
! is what lets GAMESS hand in DDI's compute communicator -- whose data
! servers sit in MPI_COMM_WORLD but never call into terco, so terco must
! never wait on MPI_COMM_WORLD itself -- and what will let a fragment
! method hand in a per-fragment sub-communicator. See trc_build_comm.
!
! CARTESIAN COMPONENT ORDER
! -------------------------
! Every AO-indexed array this API hands back or takes in (S/T/V, the Fock
! and density matrices, multipoles, DF tensors, `trc_set_guess`'s density)
! uses the WITHIN-SHELL Cartesian order the library was BUILT with
! (TRC_CART_ORDER at compile time: lexicographic by default, or gamess).
! `trc_cart_order` below answers which one a caller is linked against, so a
! caller need not know it out of band; it is a build-wide constant, not
! something that varies per handle, so it takes none.
!
! AO NORMALISATION
! ----------------
! Separate from the order above, and a RUNTIME choice, not a build one:
! trc_set_ao_norm(handle, TRC_AO_NORM_SHELL | TRC_AO_NORM_UNIT) on a
! context, default TRC_AO_NORM_SHELL (terco's and libcint's own convention,
! so nothing changes for a caller that never calls it). TRC_AO_NORM_UNIT is
! GAMESS's: every Cartesian component individually normalised to 1. See
! trc_ao_norm.F90's header for the closed form (n = sqrt((2l-1)!! /
! product over x,y,z of (2*power-1)!!), FURTHER divided by terco's own
! libcint-inherited `common_fac_sp` scale, which is 1 for d and up but not
! for s or p -- also in that header, found the hard way) and the
! covariant/contravariant transform it implies -- M_GAMESS = N M_terco N
! for a two-AO-index matrix, D_GAMESS = N^-1 D_terco N^-1 for a density,
! C_GAMESS = N^-1 C_terco (one index only: the other is an orbital label,
! not an AO).
!
! The setting lives on `basis_box` (propagated from a context's own choice
! by `take_basis`, so a handle borrowed through trc_context_basis carries
! it), not only on the context, because that is the object every AO-array
! entry below actually holds.
!
! SPHERICAL BASES
! ----------------
! A context's orbital basis may be SPHERICAL (trc_api.F90's
! `trc_basis_t%spherical`, set by trc_basis_json.F90 from a JSON basis's
! own function_type, "gto_spherical" meaning yes; `trc_set_force_cartesian`
! overrides it to Cartesian regardless). When it is, every SCF-FACING
! AO-indexed entry below -- trc_nao, trc_density, trc_mo_coefficients,
! trc_mo_energies, trc_set_guess's `dguess`, and `trc_run_scf` itself --
! is sized at `nao_sph`, not the Cartesian `nao`; `trc_nao` is the single
! place that size is computed (`scf_space_nao`) and every other entry
! reads it from there, so a caller sizing its own buffers from `trc_nao`
! never has to know which case it is in. `trc_is_spherical` answers
! directly for a caller that wants to know anyway.
!
! The raw-basis-handle entries (trc_compute_1e, trc_compute_multipoles,
! trc_compute_df2c/df3c, capi_fock and its _many and _nosym siblings,
! capi_scf/capi_scf_mpi through scf_entry) are NOT affected and stay at
! the Cartesian `nao` of `trc_basis_nao` unconditionally: a basis built
! through `trc_basis_create*` has no JSON file to ever set `spherical`
! from, by construction (trc_api.F90's `basis_build` always leaves it
! `.false.`), so there is nothing for these entries to switch on.
! RI-MP2 (trc_run_rimp2) and `--df` (trc_scf_df.F90) both support a
! spherical orbital basis now: their three-centre tensor stays Cartesian
! (every kernel does), but `trc_run_rimp2` folds `cx%res%cmo`/`eps` into
! that space first (`rimp2_mo_to_cart`, zero-padded past the SCF's own
! spherical MO count -- trc_rimp2_driver.F90's own header), and `--df`'s
! `d_c`/`gmat_c` arrive already Cartesian through `trc_scf_fock.F90`'s
! `to_cart`/`to_sph`, which is every AO index `trc_scf_df.F90` ever reads.
! The AUXILIARY basis (rifit for RI-MP2, jkfit for `--df`) follows the
! ORBITAL basis's own convention (pyscf's `auxmol.cart = mol.cart`), not
! necessarily its own JSON tag -- see `basis_json`'s and
! `trc_set_jkfit_json`'s own comments -- and is folded on its own
! auxiliary (P) index inside `trc_scf_df.F90`/`trc_rimp2_driver.F90` when
! it ends up spherical.
!
! TRC_AO_NORM_UNIT AND A SPHERICAL BASIS DO NOT COMBINE. UNIT rescales
! each CARTESIAN component of a shell individually (GAMESS's convention);
! a spherical AO is already a fixed linear combination of several
! Cartesian components, each with its own ratio, so "the UNIT-normalised
! version of a spherical d_z2" is not a well-defined request under
! anything this module implements. `trc_set_ao_norm(..., TRC_AO_NORM_UNIT)`
! refuses if the context's current basis is spherical, and loading a
! spherical basis while UNIT mode is already set is refused the other
! way, at `trc_set_basis_json`/`trc_set_aux_json`. Neither direction
! silently falls back to one interpretation of the two.
!
! EVERY entry that returns or accepts an AO-indexed array, classified:
!
!   WIRED (scales when the basis/context involved is TRC_AO_NORM_UNIT):
!     trc_compute_1e        covariant, S and T and V each
!     trc_compute_multipoles covariant, one multipole component at a time
!     trc_density            contravariant (context)
!     trc_mo_coefficients    contravariant, AO index only (context)
!     trc_set_guess           covariant-IN for TRC_GUESS_GIVEN, once, at the
!                             copy (context) -- see that entry's own comment
!     capi_fock (trc_fock), trc_fock_many, trc_fock_nosym
!         dmat contravariant in, gmat covariant out. This is GAMESS's own
!         MAIN caller -- its Fock/J/K build, run_Fockbuild's role in
!         hpcchem_driver.src -- so it is wired even though it costs the
!         one exception to this module's usual zero-copy Fock path: in
!         TRC_AO_NORM_UNIT, `dmat` is copied to a host scratch array and
!         scaled to terco's own convention THERE (never touching the
!         caller's buffer -- see fock_scale_in), and `gmat` is scaled in
!         place AFTER the build, since it is terco's own to write anyway
!         (fock_scale_out). One O(n^2) copy against an O(n^4) Fock build.
!         TRC_AO_NORM_SHELL, the default, still takes the original
!         sequence-associated NO-COPY path unchanged (see fock_shaped's
!         own comment), so no existing caller pays anything for this.
!         trc_fock_many's `nd` densities (density index fastest) use
!         their own small scale-in/out pair for that layout.
!     capi_scf (trc_scf), capi_scf_mpi (trc_scf_mpi), through scf_entry
!         dguess contravariant in, dmat contravariant out, eps invariant;
!         wired for the same reason as the Fock family -- a UNIT-mode
!         caller can reach this with a context-borrowed `basis` exactly
!         as it can reach trc_fock. `res%dmat`/`bb%cmo`/`bb%eps` stay in
!         terco's own convention internally (trc_run_rimp2 reads them
!         that way); only the copies handed back across the C boundary
!         are scaled.
!
!   INVARIANT, untouched by design (see trc_ao_norm.F90's header for why):
!     trc_mo_energies, trc_energy, trc_energy_parts, trc_converged,
!     trc_iterations, trc_rimp2_energy, trc_nao, trc_nmo, trc_nspin,
!     trc_mo_occupations (MO-indexed, like trc_mo_energies, not AO-indexed),
!     trc_message, trc_set_pcm, trc_pcm_energy, trc_gradient, trc_polarizability
!     (a nuclear-coordinate array and a Cartesian response tensor, neither an
!     AO one -- an SCF energy derivative or response does not depend on
!     which AO normalisation convention computed it)
!
!   NOT AO-INDEXED DATA (build or borrow a HANDLE; nothing to scale here,
!   the handle's own ao_norm_kind is what a later entry reads):
!     trc_basis_create*, trc_pairs_create, trc_eri_create,
!     trc_context_basis, trc_context_eri, trc_set_ao_norm, trc_get_ao_norm
!
!   COVARIANT/CONTRAVARIANT BY THE SAME RULE, BUT NOT WIRED IN THIS RELEASE
!   -- always behave as TRC_AO_NORM_SHELL regardless of what a context they
!   were reached through is set to, which is safe (identical to every
!   caller before this setting existed) but not yet convention-aware:
!     trc_compute_df2c, trc_compute_df3c
!         covariant, one N per AO index (df3c's auxiliary index scales by
!         the AUXILIARY basis's own N, not the orbital basis's -- it is a
!         basis function too, just a different basis). Unreached by
!         trc_set_ao_norm in practice: context never lends out `ab` (no
!         trc_context_aux exists), so every `aux` handle these two ever see
!         was built by trc_basis_create* directly and is always
!         TRC_AO_NORM_SHELL.
!
module trc_capi
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_f_pointer, c_int, c_loc, &
                                          c_null_char, c_null_ptr, c_ptr
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, trc_bind_device_env, &
                        trc_1e, trc_df_2c, trc_df_3c, trc_multipoles, &
                        TRC_NMULT
   use trc_eri, only: trc_eri_t, eri_lmax_refusal
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run, trc_grid_cache_t, trc_pair_eri_cache_t
   use trc_scf_grad, only: trc_scf_gradient
   use trc_scf_settings, only: scf_spin_rule, eri_thresh_env
   use trc_scf_df, only: df_jk_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_sad, only: trc_sad_build, trc_sadq_build
   use trc_sap, only: trc_sap_build
   use trc_minao, only: trc_minao_build
   use trc_error, only: error_t, ERROR_IO
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, xc_is_double_hybrid, &
                                XC_AVAILABLE_FUNCTIONALS, XC_FAMILY_MGGA
   use trc_tables, only: TRC_ERI_LMAX => LMAX
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use pic_mpi_lib, only: comm_t, comm_from_handle, bcast
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   use trc_linalg, only: trc_linalg_t, trc_release_local_mem_reservation
   use trc_rimp2_driver, only: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, rimp2_mo_to_cart
   use trc_scf_file, only: trc_scf_save, trc_scf_load
   use trc_cart_perm, only: trc_cart_order_id
   use trc_ao_norm, only: ao_norm_shell_ => TRC_AO_NORM_SHELL, ao_norm_unit_ => TRC_AO_NORM_UNIT, &
                          ao_norm_vector, ao_scale_cov, ao_scale_contra, ao_scale_contra_row
   use trc_sph, only: sph_cart_to_sph, sph_sph_to_cart
   use trc_hirshfeld, only: trc_hirshfeld_populations
   use trc_fukui, only: trc_fukui_indices
   use trc_pcm, only: pcm_model_from_name, PCM_MODEL_NONE, PCM_DEFAULT_LEBEDEV, PCM_XI_ORDERS
   use trc_pcm_radii, only: PCM_DEFAULT_RADII_SCALE
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_scf_polarizability
   use trc_cpks, only: trc_scf_polarizability_rks
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_constants, only: BOHR_TO_ANGSTROM, HARTREE_TO_EV, HARTREE_TO_KCALMOL, AMU_TO_ME, KB_SI, NA_SI
   implicit none
   private

   public :: trc_basis_create, trc_basis_destroy, trc_basis_nao
   public :: trc_basis_create_libcint, trc_basis_maxl
   public :: trc_basis_shell_count, trc_basis_shell_max_primitives
   public :: trc_basis_shell_info, trc_basis_shell_primitives
   public :: trc_basis_is_spherical
   public :: trc_physical_constant
   public :: trc_pairs_create, trc_pairs_destroy, trc_pairs_count
   public :: trc_compute_1e, trc_compute_df2c, trc_compute_df3c
   public :: trc_compute_multipoles, trc_multipole_count
   public :: trc_eri_create, trc_eri_destroy
   ! A binding label is a global identifier, and so is a module name, and
   ! the standard forbids the two to coincide. The modules behind the
   ! `trc_fock` and `trc_scf` entries are therefore named trc_eri and
   ! trc_scf_driver: with a module called trc_fock in scope, gfortran folded
   ! every call into that module onto the `trc_fock` label and the entry
   ! called itself until the stack ran out, which nvfortran happened to
   ! tolerate. The Fortran function names are capi_* for the same reason;
   ! the C symbols are the ABI and are unchanged.
   public :: capi_fock, trc_fock_many, trc_fock_nosym
   public :: trc_create, trc_destroy, trc_set_molecule, trc_set_charge, trc_set_basis_libcint, trc_set_basis_arrays, &
             trc_set_basis_json, trc_set_aux_libcint, trc_set_aux_arrays, trc_set_aux_json, trc_set_method, trc_set_nlc_grid, &
             trc_set_unrestricted, trc_set_level_shift, &
             trc_set_convergence, trc_set_lindep, trc_set_screening, trc_set_incremental, trc_set_df, trc_set_jkfit_json, &
             trc_set_jkfit_libcint, &
             trc_set_batching, trc_set_guess, trc_set_comm, trc_set_verbose, &
             trc_set_rimp2, trc_run_scf, trc_run_rimp2, trc_nao, trc_nspin, trc_energy, trc_energy_parts, &
             trc_converged, trc_iterations, trc_density, trc_mo_coefficients, trc_mo_energies, &
             trc_rimp2_energy, trc_message, trc_context_basis, trc_context_eri
   public :: trc_run_gradient, trc_gradient
   public :: trc_set_geometry
   public :: trc_set_hessian_step, trc_hessian_step, trc_run_hessian, trc_hessian
   public :: trc_save_scf, trc_load_scf
   public :: trc_set_pcm, trc_set_pcm_radii, trc_pcm_energy
   public :: trc_set_ao_norm, trc_get_ao_norm
   public :: trc_set_force_cartesian, trc_is_spherical
   public :: trc_max_angular_momentum
   public :: capi_scf, capi_scf_mpi, capi_rimp2
   public :: trc_cart_order
   public :: trc_functional_info, trc_functional_list
   public :: trc_abi_version
   public :: trc_run_polarizability, trc_polarizability
   public :: trc_nmo, trc_mo_occupations
   public :: trc_hirshfeld_charges
   public :: trc_fukui_hirshfeld

   ! Status codes. Zero is success; everything else is a reason.
   integer(c_int), parameter, public :: TRC_OK            = 0
   integer(c_int), parameter, public :: TRC_ERR_NULL      = 1
   integer(c_int), parameter, public :: TRC_ERR_BADARG    = 2
   integer(c_int), parameter, public :: TRC_ERR_UNSUPPORTED = 3
   integer(c_int), parameter, public :: TRC_ERR_NOCONV = 4   !! the SCF ran out of iterations
   integer(c_int), parameter, public :: TRC_ERR_STATE = 5    !! called before what it needs was set
   ! Guess kinds for trc_set_guess.
   integer(c_int), parameter, public :: TRC_GUESS_CORE = 0, TRC_GUESS_GWH = 1, TRC_GUESS_SAD = 2, &
                                        TRC_GUESS_GIVEN = 3, TRC_GUESS_SADQ = 4, &
                                        TRC_GUESS_SAP = 5, TRC_GUESS_PROJECTED = 6, &
                                        TRC_GUESS_MINAO = 7
   !! One declaration, so one FORD block; each kind that needs more than
   !! its name has its own paragraph below.
   !!
   !! **TRC_GUESS_MINAO**: the minimal-ANO density (trc_minao's own
   !! trc_minao_build) projected into this run's basis -- no free-atom
   !! SCF, so nothing cached; refused for an element above Kr or a
   !! Cartesian basis with occupied d (trc_minao's own checks). Always
   !! NEUTRAL atoms, as PySCF's init_guess_by_minao: for an ion the guess
   !! carries the neutral electron count, unlike TRC_GUESS_SADQ, and the
   !! first Roothaan step fixes the count.
   !!
   !! **TRC_GUESS_PROJECTED**: the density from THIS context's own last
   !! orbitals (cx%bb%cmo, stored by any run with at least one iteration,
   !! not only a converged one), re-orthonormalised in the CURRENT
   !! geometry's overlap (trc_scf_context's `project_density`) -- a warm
   !! start for a SEQUENCE of nearby geometries at the SAME molecule and
   !! basis (a geometry optimiser's steps, or a finite-difference Hessian's
   !! displacements) that is exact when the geometry has not moved at all
   !! and degrades gracefully as it does, unlike reusing the raw density
   !! unprojected. Refused (TRC_ERR_STATE) if no previous SCF at this nao
   !! has run yet, if it was run with fractional occupation (`project_density`
   !! only builds integer occupations), or at a size/spin mismatch. See
   !! trc_run_scf's own dispatch; test/bench_warmstart.F90 measured it
   !! against the plain density warm start at water/def2-SVP, steps
   !! 0.005-1.0 Bohr: the two warm starts needed the same iteration count
   !! at every step (HF 10/13/15, B3LYP 8/10/12).
   !! Appended, not inserted: 0-3 are what they have always been, so a
   !! caller built against the old header keeps working.
   ! AO normalisation kinds for trc_set_ao_norm -- re-exported from
   ! trc_ao_norm (renamed on import, above, only so this c_int redeclaration
   ! does not clash with the plain-integer names it reads) so a C caller
   ! needs only this header, not that module's own Fortran-only one. See
   ! trc_ao_norm's header for the convention itself and this module's own
   ! "AO NORMALISATION" section below for which entries it applies to.
   integer(c_int), parameter, public :: TRC_AO_NORM_SHELL = ao_norm_shell_
   integer(c_int), parameter, public :: TRC_AO_NORM_UNIT = ao_norm_unit_

   ! This ABI's own version, read back by trc_abi_version so a host can
   ! compare it against the constants its own compiled-against header
   ! carries before trusting anything else in this module. MINOR for an
   ! additive change, MAJOR (with MINOR reset to 0) for anything that
   ! changes an existing signature or enum meaning -- see trc_abi_version's
   ! own comment.
   integer(c_int), parameter, public :: TRC_ABI_MAJOR = 1
   integer(c_int), parameter, public :: TRC_ABI_MINOR = 11

   ! trc_physical_constant's own selector: which of trc_constants' values
   ! a caller is asking for. Appended at the end, like TRC_GUESS_*, so a
   ! caller built against an older header keeps working.
   integer(c_int), parameter, public :: TRC_CONST_BOHR_TO_ANGSTROM = 0
   integer(c_int), parameter, public :: TRC_CONST_HARTREE_TO_EV = 1
   integer(c_int), parameter, public :: TRC_CONST_HARTREE_TO_KCALMOL = 2
   integer(c_int), parameter, public :: TRC_CONST_AMU_TO_ME = 3
   integer(c_int), parameter, public :: TRC_CONST_KB_SI = 4
   integer(c_int), parameter, public :: TRC_CONST_NA_SI = 5

   ! Wrappers so a bare derived type can be pointed at from C.
   type :: basis_box
      type(trc_basis_t) :: b
      ! The orbitals of the last trc_scf on this basis, kept for a
      ! correlated step that follows: trc_scf hands the caller a density and
      ! eigenvalues, not orbitals, and its signature is not changing under
      ! callers already written to it.
      real(dp), allocatable :: cmo(:, :, :), eps(:, :)
      integer :: nalpha = -1, nbeta = -1
      ! TRC_AO_NORM_SHELL by default: an independently-created basis (this
      ! box built through trc_basis_create*, never touched by
      ! trc_set_ao_norm) always reports and accepts terco's own convention,
      ! which is what makes the new setting opt-in. A context propagates its
      ! own choice here (see `take_basis`) so that a handle borrowed through
      ! trc_context_basis carries it to whatever raw entry receives it.
      integer :: ao_norm_kind = TRC_AO_NORM_SHELL
   end type basis_box

   type :: pairs_box
      type(trc_pairlist_t) :: p
   end type pairs_box

   type :: eri_box
      type(trc_eri_t) :: e
   end type eri_box

   !
   ! THE CONTEXT. One handle for a whole calculation: the molecule, the
   ! basis and the auxiliary basis, every setting, the guess, the
   ! communicator, and after a run its results. Set up in steps, in any
   ! order, then run. The basis and ERI objects inside are the same boxes the
   ! three-handle entries use, so the Fock family can borrow them through
   ! trc_context_basis / trc_context_eri without a second copy.
   !
   type :: context_box
      ! molecule
      integer :: natm = 0, charge = 0, multiplicity = 1
      real(dp), allocatable :: at_z(:), at_r(:, :)
      integer, allocatable :: elements(:)
      !! Element per atom (1..118), trc_set_elements's own array -- defaults
      !! to nint(at_z) (trc_set_molecule resets it there). Drives JSON basis
      !! lookup, DFT/NLC grid sizes and guess shape; `at_z` (0 for a ghost)
      !! stays what drives the physics. See trc_set_elements.
      logical :: have_mol = .false.
      ! bases; `bb` also keeps the last SCF's orbitals, as the old entries do
      type(basis_box) :: bb, ab
      logical :: have_basis = .false., have_aux = .false.
      ! the density-fitting (jkfit) basis for `--df`, separate from `ab`
      ! (RI-MP2's rifit) because one run can do DF-SCF and RI-MP2 with two
      ! different auxiliary sets (see trc_scf_df.F90's own header).
      type(basis_box) :: jb
      logical :: have_jkfit = .false.
      ! the density-fitted J/K object itself: fitted once per (molecule,
      ! basis, jkfit) and carried across every `trc_run_scf` on this
      ! handle, so a chained RI-MP2 step could reuse `bpk` instead of
      ! redoing the metric/Cholesky/3c/trsm (not wired up yet -- see
      ! trc_scf_df.F90's `df_jk_dims`). `df_fit_ok` is cleared by anything
      ! that invalidates the fit (a new basis, molecule or jkfit set).
      type(df_jk_t) :: dfjk
      logical :: df_fit_ok = .false.
      ! The file paths trc_set_aux_json/trc_set_jkfit_json were last given:
      ! trc_run_rimp2's only way to tell whether RI-MP2's own auxiliary
      ! basis is the SAME one `dfjk` was fitted on (same JSON path/
      ! contents), the condition the chained path (skip metric/Cholesky/
      ! 3c/trsm, read `dfjk%packed_b()` straight) is restricted to.
      character(len=1024) :: aux_path = "", jkfit_path = ""
      ! settings
      type(trc_scf_options_t) :: opts
      integer :: guess_kind = TRC_GUESS_SAD
      real(dp), allocatable :: dguess(:, :, :)
      integer :: fcomm = -1
      integer :: nfrozen = -1, aux_block = 256
      ! TRC_AO_NORM_SHELL by default: terco's own per-l convention, so no
      ! existing caller sees a difference. trc_set_ao_norm sets this AND,
      ! if a basis is already present, `bb`/`ab`'s own copies (see
      ! `take_basis`, which propagates it into every basis built afterwards
      ! too); the entries that actually convert read `bb`/`ab`'s copy, not
      ! this one, because they are the ones a caller can also reach without
      ! a context (trc_context_basis's borrowed handle).
      integer :: ao_norm_kind = TRC_AO_NORM_SHELL
      ! the ERI context, only if something borrows it
      type(eri_box), allocatable :: eb
      ! cuBLAS/cuSOLVER handles, created at this context's first run (after
      ! the device is bound) and carried through every run after it: the
      ! SAD atoms, the SCF, XC and RI-MP2 all use these, never their own.
      ! Destroyed with the context.
      type(trc_linalg_t) :: la
      ! results
      type(trc_scf_result_t) :: res
      integer :: nalpha = 0, nbeta = 0
      logical :: scf_ok = .false., rimp2_ok = .false.
      real(dp) :: e_os = 0.0_dp, e_ss = 0.0_dp
      ! trc_run_gradient's own result: (3, natm), Eh/Bohr, in terco's atom
      ! order. `grad_ok` is cleared by anything that invalidates `scf_ok`
      ! too (invalidate_runs) -- a gradient at a geometry/density this
      ! context no longer holds is not a gradient of anything.
      real(dp), allocatable :: grad(:, :)
      logical :: grad_ok = .false.
      ! trc_run_polarizability's own result: the static dipole
      ! polarizability (a.u.) at the converged SCF, CPHF (RHF/UHF) or
      ! CPKS (RKS). Cleared by invalidate_runs, same as grad_ok.
      real(dp) :: alpha(3, 3) = 0.0_dp
      logical :: polar_ok = .false.
      ! `trc_run_scf`'s own DFT grid(s), kept alive (device-resident) for
      ! `trc_run_gradient` to reuse instead of rebuilding -- see
      ! trc_scf_context.F90's `trc_grid_cache_t` and trc_scf_driver.F90's
      ! `trc_scf_run`'s own `grid_cache` argument. `kgrid_level`/
      ! `kgrid_max_pts`/`kgrid_rho_tol`/`kgrid_nlc_level` are the settings
      ! it was built at, snapshotted the moment `trc_run_scf` succeeds, so
      ! `trc_run_gradient` can check them against `cx%opts`'s CURRENT
      ! values rather than trust `scf_ok` alone. Released (and these reset)
      ! by `invalidate_runs`, by `trc_load_scf` (a loaded result has no
      ! grid of its own -- see that function's own header) and by
      ! `trc_destroy`.
      type(trc_grid_cache_t) :: kgrid
      integer :: kgrid_level = -1, kgrid_max_pts = -1, kgrid_nlc_level = -1
      real(dp) :: kgrid_rho_tol = -1.0_dp
      ! `trc_run_scf`'s own pair list and (outside `--df`) four-centre
      ! engine, kept alive the same way and for the same reason as
      ! `kgrid` right above -- see trc_scf_context.F90's
      ! `trc_pair_eri_cache_t`. No settings snapshot of its own: its
      ! `thresh_built` IS the setting to check (trc_scf_gradient's own
      ! `pair_eri_cache` header), compared bit for bit against whatever
      ! threshold the gradient call resolves at the time. Released (and
      ! reset) wherever `kgrid` is.
      type(trc_pair_eri_cache_t) :: peri
      ! trc_run_hessian's own result: the Cartesian Hessian (3*natm,
      ! 3*natm), row/column k = 3*(atom-1)+{x,y,z}, Eh/Bohr^2, at the
      ! geometry/density this context held converged when it was called
      ! (the context is restored to that same, reconverged state before
      ! trc_run_hessian returns -- see that function's own header).
      ! `hess_step` is the displacement h, Bohr, trc_set_hessian_step's
      ! own default 0.005. `hess_ok` is cleared by invalidate_runs, same
      ! as `grad_ok`: a Hessian at a geometry/density this context no
      ! longer holds is not a Hessian of anything.
      real(dp) :: hess_step = 0.005_dp
      real(dp), allocatable :: hess(:, :)
      logical :: hess_ok = .false.
      character(len=256) :: message = ""
   end type context_box

contains

   !
   ! The SIZE of the AO space the SCF-facing entries of this file read and
   ! write -- trc_nao, trc_density, trc_mo_coefficients, trc_mo_energies,
   ! trc_set_guess, trc_run_scf. `nao_sph` when the basis this box holds is
   ! spherical, `nao` (Cartesian) otherwise -- the same choice
   ! trc_scf_driver itself makes (trc_scf_driver.F90's own `sph`/`nao` split), so
   ! this is the one place that choice is restated at the C boundary rather
   ! than copied into every entry that needs it. The raw-basis-handle
   ! family (trc_basis_nao, trc_compute_1e, capi_fock and friends) does NOT
   ! use this: those operate on whatever `trc_basis_create*` built
   ! directly, which has no file to read a function_type from and so is
   ! never spherical (see trc_api.F90's `basis_build`, which always leaves
   ! `spherical` at its default .false.) -- `bb%b%nao` is already right
   ! for them, unconditionally.
   !
   pure integer function scf_space_nao(bb) result(n)
      type(basis_box), intent(in) :: bb
      n = merge(bb%b%nao_sph, bb%b%nao, bb%b%spherical)
   end function scf_space_nao

   !
   ! Build a basis from flat per-shell arrays.
   !
   ! `centres` is (3, nshell) and `exps`/`coefs` are (maxnp, nshell), both in
   ! Fortran order -- which for C means column-major, i.e. the fastest index is
   ! the first. A C caller writing exps[shell][prim] has it transposed; the
   ! header comment says so rather than the library guessing.
   !
   function trc_basis_create(nshell, maxnp, l, nprim, exps, coefs, centres, &
                               natm, zatm, ratm, handle) &
      result(status) bind(c, name="trc_basis_create")
      integer(c_int), value :: nshell, maxnp, natm
      integer(c_int), intent(in) :: l(nshell), nprim(nshell)
      real(c_double), intent(in) :: exps(maxnp, nshell), coefs(maxnp, nshell)
      real(c_double), intent(in) :: centres(3, nshell)
      real(c_double), intent(in) :: zatm(natm), ratm(3, natm)
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status

      type(basis_box), pointer :: box
      integer :: i

      status = TRC_ERR_BADARG
      handle = c_null_ptr
      if (nshell <= 0 .or. natm <= 0 .or. maxnp <= 0) return
      do i = 1, nshell
         if (nprim(i) < 1 .or. nprim(i) > maxnp) return
         if (l(i) < 0) return
      end do

      allocate (box)
      call box%b%build(int(nshell), int(l), int(nprim), &
                       real(exps, dp), real(coefs, dp), real(centres, dp), &
                       int(natm), real(zatm, dp), real(ratm, dp), int(maxnp))
      call box%b%to_device()
      handle = c_loc(box)
      status = TRC_OK
   end function trc_basis_create

   function trc_basis_destroy(handle) result(status) &
      bind(c, name="trc_basis_destroy")
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(basis_box), pointer :: box
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      call box%b%release()
      deallocate (box)
      status = TRC_OK
   end function trc_basis_destroy

   function trc_basis_nao(handle, nao) result(status) &
      bind(c, name="trc_basis_nao")
      !! Number of atomic orbitals, so a caller can size its own output buffers.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: nao
      integer(c_int) :: status
      type(basis_box), pointer :: box
      nao = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      nao = int(box%b%nao, c_int)
      status = TRC_OK
   end function trc_basis_nao

   function trc_pairs_create(basis, thresh, handle) result(status) &
      bind(c, name="trc_pairs_create")
      type(c_ptr), value :: basis
      real(c_double), value :: thresh
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(pairs_box), pointer :: pb
      handle = c_null_ptr
      status = TRC_ERR_NULL
      if (.not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      allocate (pb)
      call pb%p%build(bb%b, real(thresh, dp))
      call pb%p%to_device()
      handle = c_loc(pb)
      status = TRC_OK
   end function trc_pairs_create

   function trc_pairs_destroy(handle) result(status) &
      bind(c, name="trc_pairs_destroy")
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(pairs_box), pointer :: pb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, pb)
      call pb%p%release()
      deallocate (pb)
      status = TRC_OK
   end function trc_pairs_destroy

   function trc_pairs_count(handle, npair) result(status) &
      bind(c, name="trc_pairs_count")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: npair
      integer(c_int) :: status
      type(pairs_box), pointer :: pb
      npair = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, pb)
      npair = int(pb%p%npair, c_int)
      status = TRC_OK
   end function trc_pairs_count

   !
   ! S, T and V into caller-provided host buffers of nao*nao each.
   !
   ! The buffers are host memory and the library copies the result back. A
   ! device-pointer variant belongs here later for a caller already holding
   ! device memory; it is a second entry point, not a flag, because getting
   ! the two confused is a segfault rather than a wrong answer.
   !
   function trc_compute_1e(basis, pairs, smat, tmat, vmat) result(status) &
      bind(c, name="trc_compute_1e")
      !! COVARIANT under `basis`'s own ao_norm_kind (set through
      !! trc_set_ao_norm on a context that lent this handle via
      !! trc_context_basis; a raw trc_basis_create* handle is always
      !! TRC_AO_NORM_SHELL and this is then a no-op): S_out = N S N, same
      !! for T and V.
      type(c_ptr), value :: basis, pairs
      real(c_double), intent(out) :: smat(*), tmat(*), vmat(*)
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(pairs_box), pointer :: pb
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :), nrm(:)
      integer :: n

      status = TRC_ERR_NULL
      if (.not. c_associated(basis) .or. .not. c_associated(pairs)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(pairs, pb)
      n = bb%b%nao
      allocate (s(n, n), t(n, n), v(n, n))
      !$acc enter data create(s, t, v)
      call trc_1e(bb%b, pb%p, s, t, v)
      !$acc update self(s, t, v)
      !$acc exit data delete(s, t, v)
      if (bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         allocate (nrm(n))
         call ao_norm_vector(bb%b%nshell, bb%b%sh_l, bb%b%sh_ao, n, bb%ao_norm_kind, nrm)
         call ao_scale_cov(nrm, s)
         call ao_scale_cov(nrm, t)
         call ao_scale_cov(nrm, v)
         deallocate (nrm)
      end if
      smat(1:n*n) = reshape(s, [n*n])
      tmat(1:n*n) = reshape(t, [n*n])
      vmat(1:n*n) = reshape(v, [n*n])
      deallocate (s, t, v)
      status = TRC_OK
   end function trc_compute_1e

   function trc_compute_df2c(aux, jmat) result(status) &
      bind(c, name="trc_compute_df2c")
      type(c_ptr), value :: aux
      real(c_double), intent(out) :: jmat(*)
      integer(c_int) :: status
      type(basis_box), pointer :: ab
      real(dp), allocatable :: j(:, :)
      integer :: n

      status = TRC_ERR_NULL
      if (.not. c_associated(aux)) return
      call c_f_pointer(aux, ab)
      n = ab%b%nao
      allocate (j(n, n))
      !$acc enter data create(j)
      call trc_df_2c(ab%b, j)
      !$acc update self(j)
      !$acc exit data delete(j)
      jmat(1:n*n) = reshape(j, [n*n])
      deallocate (j)
      status = TRC_OK
   end function trc_compute_df2c

   function trc_compute_df3c(basis, pairs, aux, tens) result(status) &
      bind(c, name="trc_compute_df3c")
      type(c_ptr), value :: basis, pairs, aux
      real(c_double), intent(out) :: tens(*)
      integer(c_int) :: status
      type(basis_box), pointer :: bb, ab
      type(pairs_box), pointer :: pb
      real(dp), allocatable :: t(:, :, :)
      integer :: n, m

      status = TRC_ERR_NULL
      if (.not. c_associated(basis) .or. .not. c_associated(pairs) &
          .or. .not. c_associated(aux)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(pairs, pb)
      call c_f_pointer(aux, ab)
      n = bb%b%nao
      m = ab%b%nao
      allocate (t(n, n, m))
      !$acc enter data create(t)
      call trc_df_3c(bb%b, pb%p, ab%b, t)
      !$acc update self(t)
      !$acc exit data delete(t)
      tens(1:int(n, int64)*n*m) = reshape(t, [n*n*m])
      deallocate (t)
      status = TRC_OK
   end function trc_compute_df3c

   !
   ! Build a basis straight from libcint's packed `atm`/`bas`/`env`.
   !
   ! THIS IS THE ONE A HOST CODE ACTUALLY WANTS. Anything with a libcint
   ! backend already holds these three arrays, so the alternative -- unpacking
   ! them into per-shell arrays on the caller's side and handing those to
   ! `trc_basis_create` -- is a conversion written twice, once here and once in
   ! every caller, with two chances to disagree about normalisation.
   !
   ! CARTESIAN ONLY. libcint's `bas` carries the angular form in KAPPA_OF and
   ! terco has no spherical transform; a spherical basis silently reinterpreted
   ! as Cartesian gives wrong numbers rather than an error, so the caller is
   ! required to state it and is refused if it is not Cartesian.
   !
   function trc_basis_create_libcint(atm, natm, bas, nshell, env, nenv, &
                                     cartesian, handle) result(status) &
      bind(c, name="trc_basis_create_libcint")
      integer(c_int), value :: natm, nshell, nenv, cartesian
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status
      type(basis_box), pointer :: box

      handle = c_null_ptr
      status = TRC_ERR_BADARG
      if (natm <= 0 .or. nshell <= 0 .or. nenv <= 0) return
      if (cartesian == 0) then
         status = TRC_ERR_UNSUPPORTED
         return
      end if
      allocate (box)
      call box%b%from_libcint(int(atm(1:6*natm)), int(natm), &
                              int(bas(1:8*nshell)), int(nshell), &
                              real(env(1:nenv), dp))
      call box%b%to_device()
      handle = c_loc(box)
      status = TRC_OK
   end function trc_basis_create_libcint

   function trc_basis_maxl(handle, maxl) result(status) &
      bind(c, name="trc_basis_maxl")
      !! Highest angular momentum in the basis, so a caller can decide whether to
      !! use this library at all before it builds anything expensive. The
      !! four-centre kernels stop at d and the caller has to fall back rather than
      !! get a wrong answer.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: maxl
      integer(c_int) :: status
      type(basis_box), pointer :: box
      maxl = -1
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      maxl = int(maxval(box%b%sh_l), c_int)
      status = TRC_OK
   end function trc_basis_maxl

   !
   ! Shell-data accessors.
   !
   ! A caller that already has its own basis-construction path (cyclone's
   ! GAMESS-to-terco converter is the motivating one) can otherwise check
   ! only the AO count and the final energy against its own input; these
   ! four let it read back the shells THIS handle actually holds, so a
   ! conversion bug shows up as a wrong exponent or coefficient instead of
   ! only a wrong energy several steps later. They work on ANY basis
   ! handle -- one built directly through trc_basis_create*/
   ! trc_basis_create_libcint, or borrowed from a context through
   ! trc_context_basis -- since both are the same `basis_box`.
   !
   ! Caller-provided buffers with a size query first, as trc_hessian/
   ! trc_multipole_count already do: trc_basis_shell_count for NSHELL,
   ! then trc_basis_shell_info's three NSHELL-long buffers, and
   ! trc_basis_shell_max_primitives for MAXNP to size
   ! trc_basis_shell_primitives' two MAXNP*NSHELL buffers.
   !
   function trc_basis_shell_count(handle, nshell) result(status) &
      bind(c, name="trc_basis_shell_count")
      !! Number of shells, so a caller can size trc_basis_shell_info's
      !! buffers. NOT split by decontraction (src/basis/trc_decontract.F90):
      !! that split lives inside trc_eri_t, built from a copy of this
      !! basis, and never touches the trc_basis_t this handle itself
      !! holds -- this is always the basis as the caller (or
      !! trc_basis_from_json) built it.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: nshell
      integer(c_int) :: status
      type(basis_box), pointer :: box
      nshell = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      nshell = int(box%b%nshell, c_int)
      status = TRC_OK
   end function trc_basis_shell_count

   function trc_basis_shell_max_primitives(handle, maxnp) result(status) &
      bind(c, name="trc_basis_shell_max_primitives")
      !! The second dimension trc_basis_shell_primitives' exps/coefs
      !! buffers need (MAXNP*NSHELL each): the widest shell's own
      !! primitive count, box%b%maxnp -- the same bound
      !! trc_basis_create's own `maxnp` argument sets on the way in.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: maxnp
      integer(c_int) :: status
      type(basis_box), pointer :: box
      maxnp = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      maxnp = int(box%b%maxnp, c_int)
      status = TRC_OK
   end function trc_basis_shell_max_primitives

   function trc_basis_shell_info(handle, atom, l, nprim) result(status) &
      bind(c, name="trc_basis_shell_info")
      !! Per shell, into three caller buffers of trc_basis_shell_count's
      !! own NSHELL: ATOM is the 1-based index into trc_set_molecule's own
      !! atom order (0 for a shell on no atom -- a ghost/probe function;
      !! trc_api.F90's `sh_at`, matched by an exact coordinate hit against
      !! the atom list, same convention as every other 1-based index this
      !! library's own Fortran side uses, e.g. trc_pairlist_t's `p_i`/
      !! `p_j`). L is the angular momentum. NPRIM is the number of
      !! primitives actually retained: `basis_build` drops a primitive
      !! whose coefficient is exactly zero (a dead column of a general-
      !! contraction table) once, when the basis is built, so this is the
      !! count trc_basis_shell_primitives below actually fills, which can
      !! be smaller than whatever the caller's own input table listed for
      !! that shell.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: atom(*), l(*), nprim(*)
      integer(c_int) :: status
      type(basis_box), pointer :: box
      integer :: i, n
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      n = box%b%nshell
      do i = 1, n
         atom(i) = int(box%b%sh_at(i), c_int)
         l(i) = int(box%b%sh_l(i), c_int)
         nprim(i) = int(box%b%sh_np(i), c_int)
      end do
      status = TRC_OK
   end function trc_basis_shell_info

   function trc_basis_shell_primitives(handle, exps, coefs) result(status) &
      bind(c, name="trc_basis_shell_primitives")
      !! Per-shell exponents and contraction coefficients, each
      !! (MAXNP, NSHELL) in Fortran column-major order -- the SAME layout
      !! trc_basis_create's own `exps`/`coefs` arguments take -- sized
      !! from trc_basis_shell_max_primitives and trc_basis_shell_count. A
      !! primitive past a shell's own NPRIM (trc_basis_shell_info) is
      !! zero, exactly as a caller building a basis would pad a shorter
      !! shell in the same buffer.
      !!
      !! COEFFICIENTS: `basis_build` folds libcint's common_fac_sp into the
      !! STORED coefficient exactly once, on the way in
      !! (0.282094791773878143 for an s shell, 0.488602511902919921 for a
      !! p shell, 1 for l >= 2 -- AGENTS.md's "libcint's common_fac_sp
      !! normalisation is applied exactly once, in basis construction");
      !! this entry divides that one factor back out before handing the
      !! coefficient back, the same exact inverse
      !! test/trc_test_basis.F90's `unpack_basis` already relies on.
      !!
      !! WHAT THAT RECOVERS DEPENDS ON HOW THIS BASIS WAS BUILT.
      !! trc_basis_create/trc_basis_create_libcint call `build` directly
      !! on the caller's own coefficients with NOTHING else done to them,
      !! so for a handle built either way this is EXACTLY the caller's own
      !! input, unchanged -- the case cyclone's GAMESS-to-terco conversion
      !! check (the motivating one) is in: it builds its own arrays and
      !! can compare them here bit for bit.
      !!
      !! trc_set_basis_json/trc_basis_from_json is NOT that case:
      !! trc_basis_json.F90's own `read_element` renormalises every
      !! shell's coefficients to unit primitive AND unit contraction
      !! self-overlap (its own `normalise`, BSE's/pyscf's own convention)
      !! BEFORE handing them to `build`, which folds common_fac_sp in on
      !! top of THAT. Undoing only the common_fac_sp fold then does not
      !! reach the JSON file's own literal number -- a caller on that path
      !! sees terco's fully renormalised coefficient instead (exponents
      !! are untouched either way, since `normalise` only rescales `c`).
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: exps(*), coefs(*)
      integer(c_int) :: status
      type(basis_box), pointer :: box
      integer :: i, k, n, mnp, off
      real(dp) :: unfold
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      n = box%b%nshell
      mnp = box%b%maxnp
      do i = 1, n
         select case (box%b%sh_l(i))
         case (0)
            unfold = 1.0_dp/0.282094791773878143_dp
         case (1)
            unfold = 1.0_dp/0.488602511902919921_dp
         case default
            unfold = 1.0_dp
         end select
         off = (i - 1)*mnp
         do k = 1, mnp
            exps(off + k) = real(box%b%sh_e(k, i), c_double)
            coefs(off + k) = real(box%b%sh_c(k, i)*unfold, c_double)
         end do
      end do
      status = TRC_OK
   end function trc_basis_shell_primitives

   function trc_basis_is_spherical(handle, spherical) result(status) &
      bind(c, name="trc_basis_is_spherical")
      !! Non-zero if this basis handle's shells are spherical rather than
      !! Cartesian -- box%b%spherical, set by trc_basis_from_json from a
      !! JSON basis's own function_type, the same field trc_is_spherical
      !! reads off a context's current basis (see that entry's own
      !! comment for what it changes downstream). This is the basis-
      !! handle-level answer, for a caller that holds the handle directly
      !! rather than the context it came from -- trc_context_basis hands
      !! out exactly that handle.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: spherical
      integer(c_int) :: status
      type(basis_box), pointer :: box
      spherical = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, box)
      spherical = merge(1_c_int, 0_c_int, box%b%spherical)
      status = TRC_OK
   end function trc_basis_is_spherical

   function trc_max_angular_momentum(eri_lmax, df_lmax, df_lmax_aux) result(status) &
      bind(c, name="trc_max_angular_momentum")
      !! The angular-momentum ceilings THIS build was compiled for, so a host
      !! can check a basis against them instead of hard-coding them: the
      !! four-centre ERIs / SCF (TERCO_LMAX: 2 by default, 3 for f), and the
      !! density-fitting orbital and auxiliary bases (TERCO_DF_LMAX; the
      !! auxiliary one is one higher). Compare with trc_basis_maxl.
      integer(c_int), intent(out) :: eri_lmax, df_lmax, df_lmax_aux
      integer(c_int) :: status
      eri_lmax = int(TRC_ERI_LMAX, c_int)
      df_lmax = int(DF_LMAX, c_int)
      df_lmax_aux = int(DF_LMAX_AUX, c_int)
      status = TRC_OK
   end function trc_max_angular_momentum

   function trc_multipole_count() result(n) &
      bind(c, name="trc_multipole_count")
      !! Number of Cartesian multipole components, so a caller can size its buffer
      !! without hardcoding 39 and drifting when the operator order changes.
      !! Named for what it counts rather than `trc_nmult`, which collides with the
      !! parameter of that name it returns.
      integer(c_int) :: n
      n = int(TRC_NMULT, c_int)
   end function trc_multipole_count

   !
   ! Cartesian multipoles through the octupole about `origin`, into a caller
   ! buffer of nao*nao*trc_nmult(). Component index is slowest, matching the
   ! Fortran (nao, nao, ncomp).
   !
   function trc_compute_multipoles(basis, pairs, origin, mmat) result(status) &
      bind(c, name="trc_compute_multipoles")
      !! COVARIANT under `basis`'s own ao_norm_kind, one component at a
      !! time (see trc_compute_1e).
      type(c_ptr), value :: basis, pairs
      real(c_double), intent(in)  :: origin(3)
      real(c_double), intent(out) :: mmat(*)
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(pairs_box), pointer :: pb
      real(dp), allocatable :: m(:, :, :), nrm(:)
      integer :: n, k
      integer(int64) :: ntot

      status = TRC_ERR_NULL
      if (.not. c_associated(basis) .or. .not. c_associated(pairs)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(pairs, pb)
      n = bb%b%nao
      ntot = int(n, int64)*int(n, int64)*int(TRC_NMULT, int64)
      allocate (m(n, n, TRC_NMULT))
      !$acc enter data create(m)
      call trc_multipoles(bb%b, pb%p, real(origin, dp), m)
      !$acc update self(m)
      !$acc exit data delete(m)
      if (bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         allocate (nrm(n))
         call ao_norm_vector(bb%b%nshell, bb%b%sh_l, bb%b%sh_ao, n, bb%ao_norm_kind, nrm)
         do k = 1, TRC_NMULT
            call ao_scale_cov(nrm, m(:, :, k))
         end do
         deallocate (nrm)
      end if
      mmat(1:ntot) = reshape(m, [ntot])
      deallocate (m)
      status = TRC_OK
   end function trc_compute_multipoles

   !
   ! The four-centre ERI object: Schwarz bounds, the binned pair list and the
   ! primitive-pair data, built once per geometry and reused by every Fock
   ! build. This is the expensive handle -- creating one per SCF iteration
   ! would throw away the whole point of the container.
   !
   function trc_eri_create(basis, thresh, handle) result(status) &
      bind(c, name="trc_eri_create")
      !! The four-centre engine for `basis` at screening threshold `thresh`.
      !! TRC_ERR_UNSUPPORTED, and no handle, when the basis has a shell
      !! above this build's four-centre ceiling (trc_max_angular_momentum).
      type(c_ptr), value :: basis
      real(c_double), value :: thresh
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      handle = c_null_ptr
      status = TRC_ERR_NULL
      if (.not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      ! No context to carry a message: the status alone says the basis is
      ! above this build's four-centre kernels (trc_max_angular_momentum).
      status = TRC_ERR_UNSUPPORTED
      if (len(eri_lmax_refusal(bb%b)) > 0) return
      allocate (eb)
      call eb%e%build(bb%b, real(thresh, dp))
      handle = c_loc(eb)
      status = TRC_OK
   end function trc_eri_create

   function trc_eri_destroy(handle) result(status) &
      bind(c, name="trc_eri_destroy")
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(eri_box), pointer :: eb
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, eb)
      call eb%e%release()
      deallocate (eb)
      status = TRC_OK
   end function trc_eri_destroy

   !
   ! G = j*J - k*K/2 for one density. Host buffers of nao*nao each.
   !
   ! `jfac` and `kfac` are passed by value and are NOT optional here: an
   ! optional argument has no C representation, and a caller that wants plain
   ! Hartree-Fock passes 1 and 1. Making them mandatory at the boundary also
   ! means a Kohn-Sham caller cannot forget them and silently get HF.
   !
   function capi_fock(eri, basis, dmat, gmat, jfac, kfac, dscreen, omega) &
      result(status) bind(c, name="trc_fock")
      !! CONTRAVARIANT in (dmat), COVARIANT out (gmat), under `basis`'s own
      !! ao_norm_kind: this is GAMESS's own main caller (its Fock/J/K
      !! build, run_Fockbuild's role in hpcchem_driver.src), so it is wired
      !! rather than left at TRC_AO_NORM_SHELL like the entries the "not
      !! wired" list below still names. See fock_scale_in/fock_scale_out
      !! for why SHELL mode still costs nothing extra.
      type(c_ptr), value :: eri, basis
      real(c_double), intent(in)  :: dmat(*)
      real(c_double), intent(out) :: gmat(*)
      real(c_double), value :: jfac, kfac
      integer(c_int), value :: dscreen
      !! Nonzero to weight the Schwarz bound by the density (an SCF), zero for
      !! the plain Schwarz screen (a coupled-perturbed solve, whose trial
      !! vector the solver drives to zero -- a density-keyed screen tightens
      !! as it does and the operator stops being a fixed linear map).
      !! An int rather than a logical: Fortran's logical has no guaranteed C
      !! representation and c_bool is one byte against C's four-byte int.
      real(c_double), value :: omega
      !! Range-separation parameter, bohr^-1; zero is the ordinary Coulomb
      !! operator. See eri%fock.
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      integer :: n
      real(dp), allocatable :: nrm(:), dscratch(:, :)

      status = TRC_ERR_NULL
      if (.not. c_associated(eri) .or. .not. c_associated(basis)) return
      if (.not. (omega >= 0.0_c_double)) then
         ! A negative omega reaches trc_eri.F90's `error stop`, which kills
         ! the host process; a NaN omega returns TRC_OK with an unattenuated
         ! K (the `ow < 0` guard there never fires for NaN). Catch both at
         ! the boundary instead. No context handle here (only `eri`/`basis`
         ! boxes), so -- as with `trc_fock_many`'s own `ndens < 1` check,
         ! below -- there is nowhere to stash a message; the status code is
         ! the whole report.
         status = TRC_ERR_BADARG
         return
      end if
      call c_f_pointer(basis, bb)
      call c_f_pointer(eri, eb)
      n = bb%b%nao
      if (bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         call fock_scale_in(bb, n, dmat, nrm, dscratch)
         call fock_shaped(bb%b, eb%e, n, dscratch, gmat, real(jfac, dp), &
                          real(kfac, dp), dscreen /= 0, real(omega, dp))
         call fock_scale_out(nrm, n, gmat)
      else
         ! The zero-copy path, untouched: `dmat` and `gmat` are sequence
         ! associated straight into the caller's own memory (see the
         ! comment on `fock_shaped` below), which is exactly why SHELL
         ! mode -- the default, so every existing caller -- pays nothing
         ! for this entry existing.
         call fock_shaped(bb%b, eb%e, n, dmat, gmat, real(jfac, dp), &
                          real(kfac, dp), dscreen /= 0, real(omega, dp))
      end if
      status = TRC_OK
   end function capi_fock

   subroutine fock_scale_in(bb, n, dmat, nrm, dscratch)
      !! UNIT mode only, shared by capi_fock and trc_fock_nosym (identical
      !! (n,n) shape): copies `dmat` -- never mutated, it is the caller's --
      !! into a host scratch array and scales THAT to terco's own convention,
      !! D_terco = N D_GAMESS N. One O(n^2) copy against an O(n^4) Fock build.
      type(basis_box), intent(in) :: bb
      integer, intent(in) :: n
      real(c_double), intent(in) :: dmat(*)
      real(dp), allocatable, intent(out) :: nrm(:), dscratch(:, :)
      allocate (nrm(n), dscratch(n, n))
      call ao_norm_vector(bb%b%nshell, bb%b%sh_l, bb%b%sh_ao, n, bb%ao_norm_kind, nrm)
      dscratch = reshape(real(dmat(1:n*n), dp), [n, n])
      call ao_scale_cov(nrm, dscratch)
   end subroutine fock_scale_in

   subroutine fock_scale_out(nrm, n, gmat)
      !! UNIT mode only: `gmat` is terco's own to write (intent(out) at the C
      !! boundary, so whatever was in the caller's buffer before this call is
      !! already gone), so scaling it after the build, in place, costs nothing
      !! the build itself did not already cost -- G_GAMESS = N G_terco N.
      real(dp), intent(in) :: nrm(:)
      integer, intent(in) :: n
      real(c_double), intent(inout) :: gmat(*)
      real(dp) :: g(n, n)
      g = reshape(real(gmat(1:n*n), dp), [n, n])
      call ao_scale_cov(nrm, g)
      gmat(1:n*n) = real(reshape(g, [n*n]), c_double)
   end subroutine fock_scale_out

   !
   ! NO COPY in TRC_AO_NORM_SHELL, the default. `dmat` and `gmat` arrive as
   ! assumed-size, and handing an assumed-size actual to an explicit-shape
   ! `(n, n)` dummy is sequence association -- the callee addresses the
   ! caller's memory directly.
   !
   ! The first version reshaped both into allocatable locals, which at 2773
   ! functions is four passes over 61 MB per Fock build for nothing. Reshaping
   ! is the obvious way to write this and the obvious way is wrong: the bytes
   ! are already contiguous and already in the right order, so all the copy
   ! bought was a second name for them. TRC_AO_NORM_UNIT is the one case that
   ! DOES need its own copy (fock_scale_in/fock_scale_out, above, called from
   ! capi_fock and trc_fock_nosym before/after this): an O(n^2) price only a
   ! caller who opted into UNIT mode pays, against an O(n^4) build.
   !
   subroutine fock_shaped(b, e, n, d, g, jfac, kfac, dscreen, omega)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: e
      integer, intent(in) :: n
      real(dp), intent(in)  :: d(n, n)
      real(dp), intent(out) :: g(n, n)
      real(dp), intent(in) :: jfac, kfac
      logical, intent(in) :: dscreen
      real(dp), intent(in) :: omega
      call e%fock(b, d, g, k_scale=kfac, j_scale=jfac, density_screen=dscreen, &
                 omega=omega)
   end subroutine fock_shaped

   !
   ! N densities in one pass over the integrals.
   !
   ! `dmats` and `gmats` are (ndens, nao, nao) with the DENSITY INDEX FASTEST,
   ! which is the layout the kernel wants and not the one a C caller would
   ! guess. Phase 2 measured the alternative: with the density index slowest
   ! the N atomic updates for one (mu,nu) sit nao^2 apart and the batch stops
   ! paying past N = 4. The layout is part of the contract, so it is stated
   ! here rather than transposed silently.
   !
   function trc_fock_many(eri, basis, ndens, dmats, gmats, jfac, kfac, dscreen, omega) &
      result(status) bind(c, name="trc_fock_many")
      !! CONTRAVARIANT in, COVARIANT out, under `basis`'s own ao_norm_kind
      !! -- see capi_fock; the density-index-fastest layout means the
      !! UNIT-mode scaling here is its own small routines
      !! (fock_many_scale_in/out) rather than a reuse of fock_scale_in/out.
      type(c_ptr), value :: eri, basis
      integer(c_int), value :: ndens
      real(c_double), intent(in)  :: dmats(*)
      real(c_double), intent(out) :: gmats(*)
      real(c_double), value :: jfac, kfac
      integer(c_int), value :: dscreen
      real(c_double), value :: omega
      !! See capi_fock.
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      integer :: n
      real(dp), allocatable :: nrm(:), dscratch(:, :, :)

      status = TRC_ERR_NULL
      if (.not. c_associated(eri) .or. .not. c_associated(basis)) return
      if (ndens < 1) then
         status = TRC_ERR_BADARG
         return
      end if
      if (.not. (omega >= 0.0_c_double)) then
         ! See capi_fock's own omega check.
         status = TRC_ERR_BADARG
         return
      end if
      call c_f_pointer(basis, bb)
      call c_f_pointer(eri, eb)
      n = bb%b%nao
      if (bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         call fock_many_scale_in(bb, int(ndens), n, dmats, nrm, dscratch)
         call fock_many_shaped(bb%b, eb%e, int(ndens), n, dscratch, gmats, &
                               real(jfac, dp), real(kfac, dp), dscreen /= 0, real(omega, dp))
         call fock_many_scale_out(nrm, int(ndens), n, gmats)
      else
         call fock_many_shaped(bb%b, eb%e, int(ndens), n, dmats, gmats, &
                               real(jfac, dp), real(kfac, dp), dscreen /= 0, real(omega, dp))
      end if
      status = TRC_OK
   end function trc_fock_many

   subroutine fock_many_shaped(b, e, nd, n, d, g, jfac, kfac, dscreen, omega)
      !! As `fock_shaped`, and for the same reason -- the batch is the case where
      !! the copy would have hurt most, being N times the size (in SHELL mode;
      !! UNIT mode already pays for a copy every entry that wires it, so this
      !! one is no different in kind, only in the extra `nd` loop).
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: e
      integer, intent(in) :: nd, n
      real(dp), intent(in)  :: d(nd, n, n)
      real(dp), intent(out) :: g(nd, n, n)
      real(dp), intent(in) :: jfac, kfac
      logical, intent(in) :: dscreen
      real(dp), intent(in) :: omega
      call e%fock_many(b, nd, d, g, k_scale=kfac, j_scale=jfac, &
                       density_screen=dscreen, omega=omega)
   end subroutine fock_many_shaped

   subroutine fock_many_scale_in(bb, nd, n, dmats, nrm, dscratch)
      !! UNIT mode only: as `fock_scale_in`, but `dmats` is (ndens, nao, nao)
      !! with the DENSITY INDEX FASTEST (trc_fock_many's own layout comment
      !! above), so every one of the `nd` densities gets the SAME `nrm(i)*
      !! nrm(j)` -- one basis, one scale, `nd` copies of it in the array.
      type(basis_box), intent(in) :: bb
      integer, intent(in) :: nd, n
      real(c_double), intent(in) :: dmats(*)
      real(dp), allocatable, intent(out) :: nrm(:), dscratch(:, :, :)
      integer :: i, j
      allocate (nrm(n), dscratch(nd, n, n))
      call ao_norm_vector(bb%b%nshell, bb%b%sh_l, bb%b%sh_ao, n, bb%ao_norm_kind, nrm)
      dscratch = reshape(real(dmats(1:nd*n*n), dp), [nd, n, n])
      do j = 1, n
         do i = 1, n
            dscratch(:, i, j) = dscratch(:, i, j)*nrm(i)*nrm(j)
         end do
      end do
   end subroutine fock_many_scale_in

   subroutine fock_many_scale_out(nrm, nd, n, gmats)
      !! UNIT mode only: as `fock_scale_out`, `nd` densities at a time.
      real(dp), intent(in) :: nrm(:)
      integer, intent(in) :: nd, n
      real(c_double), intent(inout) :: gmats(*)
      real(dp), allocatable :: g(:, :, :)
      integer :: i, j
      allocate (g(nd, n, n))
      g = reshape(real(gmats(1:nd*n*n), dp), [nd, n, n])
      do j = 1, n
         do i = 1, n
            g(:, i, j) = g(:, i, j)*nrm(i)*nrm(j)
         end do
      end do
      gmats(1:nd*n*n) = real(reshape(g, [nd*n*n]), c_double)
   end subroutine fock_many_scale_out

   !
   ! The same build for a density that is NOT symmetric.
   !
   ! A separate entry point rather than a flag on `trc_fock`, for the reason
   ! the Fortran side gives: the folded eight-fold digestion assumes
   ! D(mu,nu) = D(nu,mu), and a caller that passes an antisymmetric response
   ! density to the folded path gets a plausible wrong answer with no
   ! diagnostic. Two names cannot be confused by accident.
   !
   function trc_fock_nosym(eri, basis, dmat, gmat, dscreen) result(status) &
      bind(c, name="trc_fock_nosym")
      !! CONTRAVARIANT in, COVARIANT out, exactly as capi_fock -- the
      !! diagonal congruence ao_scale_cov/fock_scale_in apply does not
      !! care whether `dmat` is symmetric, so nothing here differs from
      !! the symmetric entry beyond which shared routines it reuses.
      type(c_ptr), value :: eri, basis
      real(c_double), intent(in)  :: dmat(*)
      real(c_double), intent(out) :: gmat(*)
      integer(c_int), value :: dscreen
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      integer :: n
      real(dp), allocatable :: nrm(:), dscratch(:, :)

      status = TRC_ERR_NULL
      if (.not. c_associated(eri) .or. .not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(eri, eb)
      n = bb%b%nao
      if (bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         call fock_scale_in(bb, n, dmat, nrm, dscratch)
         call fock_nosym_shaped(bb%b, eb%e, n, dscratch, gmat, dscreen /= 0)
         call fock_scale_out(nrm, n, gmat)
      else
         call fock_nosym_shaped(bb%b, eb%e, n, dmat, gmat, dscreen /= 0)
      end if
      status = TRC_OK
   end function trc_fock_nosym

   subroutine fock_nosym_shaped(b, e, n, d, g, dscreen)
      !! As `fock_shaped`.
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: e
      integer, intent(in) :: n
      real(dp), intent(in)  :: d(n, n)
      real(dp), intent(out) :: g(n, n)
      logical, intent(in) :: dscreen
      call e%fock_nosym(b, d, g, density_screen=dscreen)
   end subroutine fock_nosym_shaped


   !
   ! The whole SCF: Hartree-Fock when `functional` is empty, Kohn-Sham
   ! otherwise; restricted when nalpha == nbeta, unrestricted otherwise, so
   ! `nspin` is 1 or 2 accordingly and the caller sizes `dmat` and `eps` as
   ! (nao, nao, nspin) and (nao, nspin). `dguess` may be NULL for the core
   ! guess, or point at a (nao, nao, nspin) density -- for a restricted
   ! case the total density, for an unrestricted one alpha then beta.
   !
   ! `functional` is a NUL-terminated C string; names are the ones
   ! trc_xc_functional knows. `grid_level` is 0 to 9 in metalquicha's sense.
   ! `verbose` non-zero prints one line per iteration -- energy, its change
   ! and the RMS density change -- to standard output, since the iteration
   ! is otherwise invisible to a caller that only sees the final energy.
   ! A non-converged SCF returns TRC_ERR_NOCONV with the last energy and
   ! density filled in, so the caller can decide what that is worth.
   !
   function capi_scf(basis, nalpha, nbeta, functional, grid_level, conv_energy, &
                     conv_density, max_iter, dguess, verbose, energy, e_xc, dmat, eps, &
                     niter) result(status) bind(c, name="trc_scf")
      type(c_ptr), value :: basis
      integer(c_int), value :: nalpha, nbeta, grid_level, max_iter, verbose
      character(kind=c_char), intent(in) :: functional(*)
      real(c_double), value :: conv_energy, conv_density
      type(c_ptr), value :: dguess
      real(c_double), intent(out) :: energy, e_xc
      real(c_double), intent(out) :: dmat(*), eps(*)
      integer(c_int), intent(out) :: niter
      integer(c_int) :: status
      status = scf_entry(basis, nalpha, nbeta, functional, grid_level, conv_energy, conv_density, &
                         max_iter, dguess, verbose, energy, e_xc, dmat, eps, niter)
   end function capi_scf

   !
   ! The same SCF run collectively by every rank of a communicator: the
   ! Fock build, the XC batches and the grid partition are split across
   ! the ranks and every rank returns the identical result. `fcomm` is the
   ! Fortran handle of the communicator (MPI_Comm_c2f from C, or the
   ! integer a Fortran caller holds), and every rank of it must call this
   ! together, because building terco's own communicator from it duplicates
   ! it. Each rank binds its own device before anything else -- the caller
   ! does not.
   !
   ! ANY intra-communicator is accepted, not just MPI_COMM_WORLD: this is
   ! what lets GAMESS hand in DDI's compute communicator, whose data-server
   ! ranks never call here at all, and what lets a fragment method hand in
   ! a per-fragment sub-communicator later. -1 is trc_scf on one rank in
   ! any build; in a build without MPI every other handle is also the
   ! single rank there is, because that is what pic-mpi's serial backend's
   ! comm_from_handle gives back.
   !
   function capi_scf_mpi(fcomm, basis, nalpha, nbeta, functional, grid_level, conv_energy, &
                         conv_density, max_iter, dguess, verbose, energy, e_xc, dmat, eps, &
                         niter) result(status) bind(c, name="trc_scf_mpi")
      integer(c_int), value :: fcomm
      type(c_ptr), value :: basis
      integer(c_int), value :: nalpha, nbeta, grid_level, max_iter, verbose
      character(kind=c_char), intent(in) :: functional(*)
      real(c_double), value :: conv_energy, conv_density
      type(c_ptr), value :: dguess
      real(c_double), intent(out) :: energy, e_xc
      real(c_double), intent(out) :: dmat(*), eps(*)
      integer(c_int), intent(out) :: niter
      integer(c_int) :: status
      type(comm_t) :: comm
      logical :: collective
      call trc_build_comm(fcomm, comm, collective, status)
      if (status /= TRC_OK) then
         energy = 0.0_c_double
         e_xc = 0.0_c_double
         niter = 0
         return
      end if
      if (.not. collective) then
         ! fcomm was -1: one rank, no communicator: trc_scf.
         status = scf_entry(basis, nalpha, nbeta, functional, grid_level, conv_energy, conv_density, &
                            max_iter, dguess, verbose, energy, e_xc, dmat, eps, niter)
         return
      end if
      status = scf_entry(basis, nalpha, nbeta, functional, grid_level, conv_energy, conv_density, &
                         max_iter, dguess, verbose, energy, e_xc, dmat, eps, niter, comm)
      call comm%finalize()
   end function capi_scf_mpi

   !
   ! RI-MP2 on the orbitals of the last trc_scf (or trc_scf_mpi) on
   ! `basis`, which must have been restricted. `aux` is the auxiliary basis,
   ! made like any other through trc_basis_create*; `nfrozen` doubly
   ! occupied orbitals are left out of the correlation; `aux_block` bounds
   ! the depth of the three-index tensor held at once (0 for all of it).
   ! E_os and E_ss come back separately: MP2 is their sum, SCS and SOS are
   ! weights on them. `fcomm` is ANY intra-communicator's Fortran handle
   ! (see capi_scf_mpi); the occupied orbitals are split over its ranks and
   ! every rank gets the same numbers. -1 runs on one rank in any build.
   !
   function capi_rimp2(fcomm, basis, aux, nfrozen, aux_block, e_os, e_ss) &
      result(status) bind(c, name="trc_rimp2")
      integer(c_int), value :: fcomm
      type(c_ptr), value :: basis, aux
      integer(c_int), value :: nfrozen, aux_block
      real(c_double), intent(out) :: e_os, e_ss
      integer(c_int) :: status
      type(basis_box), pointer :: bb, ab
      type(trc_pairlist_t) :: pl
      type(trc_rimp2_result_t) :: res
      type(comm_t) :: comm
      integer :: nblk
      logical :: collective

      e_os = 0.0_c_double
      e_ss = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(basis) .or. .not. c_associated(aux)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(aux, ab)
      status = TRC_ERR_BADARG
      if (.not. allocated(bb%cmo)) return           ! no SCF on this basis yet
      if (bb%nalpha /= bb%nbeta) return             ! restricted reference only
      if (nfrozen < 0 .or. nfrozen >= bb%nalpha) return
      call trc_build_comm(fcomm, comm, collective, status)
      if (status /= TRC_OK) return
      nblk = int(aux_block)
      if (nblk <= 0) nblk = ab%b%nao
      if (.not. ab%b%on_device) call ab%b%to_device()
      call pl%build(bb%b, 1.0e-12_dp)
      call pl%to_device()
      if (collective) then
         call trc_rimp2_run(bb%b, ab%b, pl, bb%nalpha, bb%cmo(:, :, 1), bb%eps(:, 1), res, &
                            nfrozen=int(nfrozen), aux_block=nblk, comm=comm)
         call comm%finalize()
      else
         call trc_rimp2_run(bb%b, ab%b, pl, bb%nalpha, bb%cmo(:, :, 1), bb%eps(:, 1), res, &
                            nfrozen=int(nfrozen), aux_block=nblk)
      end if
      call pl%release()
      if (len_trim(res%message) > 0) then
         status = TRC_ERR_UNSUPPORTED
         return
      end if
      e_os = real(res%e_os, c_double)
      e_ss = real(res%e_ss, c_double)
      status = TRC_OK
   end function capi_rimp2

   function scf_entry(basis, nalpha, nbeta, functional, grid_level, conv_energy, &
                      conv_density, max_iter, dguess, verbose, energy, e_xc, dmat, eps, &
                      niter, comm) result(status)
      !! `dguess` CONTRAVARIANT in, `dmat` CONTRAVARIANT out, `eps`
      !! invariant, under `basis`'s own ao_norm_kind -- wired for the same
      !! reason capi_fock is: this is the entry capi_scf/capi_scf_mpi both
      !! funnel through, and a UNIT-mode caller can reach it with a
      !! context-borrowed `basis` exactly as it can reach the Fock family.
      type(c_ptr), value :: basis
      integer(c_int), value :: nalpha, nbeta, grid_level, max_iter, verbose
      character(kind=c_char), intent(in) :: functional(*)
      real(c_double), value :: conv_energy, conv_density
      type(c_ptr), value :: dguess
      real(c_double), intent(out) :: energy, e_xc
      real(c_double), intent(out) :: dmat(*), eps(*)
      integer(c_int), intent(out) :: niter
      type(comm_t), intent(in), optional :: comm
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      real(c_double), pointer :: dg(:, :, :)
      real(dp), allocatable :: nrm(:), dgscratch(:, :, :)
      integer :: n, nspin, i, s, k
      logical :: unit_mode

      status = TRC_ERR_NULL
      energy = 0.0_c_double
      e_xc = 0.0_c_double
      niter = 0
      if (.not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      n = bb%b%nao
      status = TRC_ERR_BADARG
      if (nalpha < 0 .or. nbeta < 0 .or. nalpha + nbeta > 2*n) return
      if (grid_level < 0 .or. max_iter < 1) return
      nspin = merge(2, 1, nalpha /= nbeta)
      unit_mode = bb%ao_norm_kind == TRC_AO_NORM_UNIT
      if (unit_mode) then
         allocate (nrm(n))
         call ao_norm_vector(bb%b%nshell, bb%b%sh_l, bb%b%sh_ao, n, bb%ao_norm_kind, nrm)
      end if

      opts%functional = ""
      do i = 1, len(opts%functional)
         if (functional(i) == c_null_char) exit
         opts%functional(i:i) = functional(i)
      end do
      opts%grid_level = int(grid_level)
      opts%conv_energy = real(conv_energy, dp)
      opts%conv_density = real(conv_density, dp)
      opts%max_iter = int(max_iter)
      opts%verbose = verbose /= 0

      if (c_associated(dguess)) then
         call c_f_pointer(dguess, dg, [n, n, nspin])
         if (unit_mode) then
            ! `dg` points at the CALLER's guess; copy before scaling so
            ! that guess is never touched, D_terco = N D_GAMESS N.
            allocate (dgscratch(n, n, nspin))
            dgscratch = real(dg, dp)
            do s = 1, nspin
               call ao_scale_cov(nrm, dgscratch(:, :, s))
            end do
            call trc_scf_run(bb%b, int(nalpha), int(nbeta), opts, res, dguess=dgscratch, comm=comm)
            deallocate (dgscratch)
         else
            call trc_scf_run(bb%b, int(nalpha), int(nbeta), opts, res, dguess=real(dg, dp), comm=comm)
         end if
      else
         call trc_scf_run(bb%b, int(nalpha), int(nbeta), opts, res, comm=comm)
      end if
      if (len_trim(res%message) > 0 .and. res%iterations == 0) then
         status = TRC_ERR_UNSUPPORTED   ! an unknown functional, or a grid that would not build
         return
      end if

      energy = real(res%energy, c_double)
      e_xc = real(res%e_xc, c_double)
      niter = int(res%iterations, c_int)
      ! `res%dmat` is terco's own SCF result, in terco's own convention --
      ! scale a COPY (D_GAMESS = N^-1 D_terco N^-1) for the caller-facing
      ! `dmat` and keep `res%dmat` itself untouched for the orbitals kept
      ! below (bb%cmo/bb%eps stay terco-native; trc_run_rimp2 needs them
      ! that way).
      if (unit_mode) then
         block
            real(dp), allocatable :: dscratch(:, :, :)
            allocate (dscratch(n, n, nspin))
            dscratch = res%dmat
            do s = 1, nspin
               call ao_scale_contra(nrm, dscratch(:, :, s))
            end do
            k = 0
            do s = 1, nspin
               do i = 1, n*n
                  k = k + 1
                  dmat(k) = real(dscratch(mod(i - 1, n) + 1, (i - 1)/n + 1, s), c_double)
               end do
            end do
            deallocate (dscratch)
         end block
      else
         k = 0
         do s = 1, nspin
            do i = 1, n*n
               k = k + 1
               dmat(k) = real(res%dmat(mod(i - 1, n) + 1, (i - 1)/n + 1, s), c_double)
            end do
         end do
      end if
      k = 0
      do s = 1, nspin
         do i = 1, n
            k = k + 1
            eps(k) = real(res%eps(i, s), c_double)
         end do
      end do
      ! Keep the orbitals for a correlated step; the copies above are done.
      if (allocated(bb%cmo)) deallocate (bb%cmo, bb%eps)
      call move_alloc(res%cmo, bb%cmo)
      call move_alloc(res%eps, bb%eps)
      bb%nalpha = int(nalpha)
      bb%nbeta = int(nbeta)
      status = merge(TRC_OK, TRC_ERR_NOCONV, res%converged)
   end function scf_entry


   ! ======================================================================
   ! The context surface
   ! ======================================================================

   function trc_create(handle) result(status) bind(c, name="trc_create")
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      allocate (cx)
      handle = c_loc(cx)
      status = TRC_OK
   end function trc_create

   function trc_destroy(handle) result(status) bind(c, name="trc_destroy")
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      call drop_eri(cx)
      if (cx%df_fit_ok) call cx%dfjk%release()
      call cx%la%release()
      call cx%kgrid%release()
      call cx%peri%release()
      if (cx%have_basis) call cx%bb%b%release()
      if (cx%have_aux) call cx%ab%b%release()
      if (cx%have_jkfit) call cx%jb%b%release()
      deallocate (cx)
      status = TRC_OK
   end function trc_destroy

   subroutine drop_eri(cx)
      type(context_box), intent(inout) :: cx
      if (allocated(cx%eb)) then
         call cx%eb%e%release()
         deallocate (cx%eb)
      end if
   end subroutine drop_eri

   subroutine refuse(cx, why)
      !! What a setter would not take, and what it takes instead: the reason
      !! lives in the context so trc_message answers "what do I change?"
      !! whichever driver asked.
      type(context_box), intent(inout) :: cx
      character(len=*), intent(in) :: why
      cx%message = why
   end subroutine refuse

   pure function itoa(i) result(t)
      integer, intent(in) :: i
      character(len=16) :: t
      write (t, "(i0)") i
   end function itoa

   pure function elements_refusal(cx) result(why)
      !! What trc_run_scf refuses about `cx%elements`/`cx%at_z` together,
      !! or an empty string if there is nothing to refuse. Element 0 at
      !! any atom builds an element-0 grid/guess that silently does the
      !! wrong thing rather than erroring (the bug ghost atoms would have
      !! reopened had this not been here); every atom a ghost leaves
      !! nothing to compute. Ghosts with PCM are refused until a cavity
      !! model is chosen: pyscf (and mqc's CPU backend) give a ghost a full
      !! sphere, the ghost-free cavity below gives it none, and the two
      !! differ. trc_pcm keeps the no-sphere mask, so lifting this is one
      !! test.
      type(context_box), intent(in) :: cx
      character(len=:), allocatable :: why
      integer :: i
      logical :: any_real, any_ghost
      why = ""
      any_real = .false.
      any_ghost = .false.
      do i = 1, cx%natm
         if (cx%elements(i) < 1) then
            why = "atom "//trim(itoa(i))//" has element "//trim(itoa(cx%elements(i)))// &
                  ", which is not a periodic table entry (trc_set_elements takes 1..118)"
            return
         end if
         if (cx%at_z(i) /= 0.0_dp) any_real = .true.
         if (cx%at_z(i) == 0.0_dp) any_ghost = .true.
      end do
      if (.not. any_real) then
         why = "every atom is a ghost (charge 0) -- there is no nucleus and no electron "// &
               "left to run an SCF on"
      else if (any_ghost .and. len_trim(cx%opts%pcm_model) > 0) then
         why = "PCM with ghost atoms is not supported yet: whether a ghost gets a cavity "// &
               "sphere is undecided (pyscf gives it one, terco's ghost-free cavity would not)"
      end if
   end function elements_refusal

   subroutine invalidate_runs(cx)
      !! A setting the results depended on changed.
      type(context_box), intent(inout) :: cx
      cx%scf_ok = .false.
      cx%rimp2_ok = .false.
      cx%grad_ok = .false.
      cx%polar_ok = .false.
      cx%hess_ok = .false.
      call cx%kgrid%release()
      call cx%peri%release()
      cx%kgrid_level = -1
      cx%kgrid_max_pts = -1
      cx%kgrid_nlc_level = -1
      cx%kgrid_rho_tol = -1.0_dp
      cx%message = ""
   end subroutine invalidate_runs

   function trc_set_molecule(handle, natm, z, xyz, charge, multiplicity) result(status) &
      bind(c, name="trc_set_molecule")
      !! Atoms in Bohr, nuclear charges as doubles, atom-major coordinates.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm, charge, multiplicity
      real(c_double), intent(in) :: z(*), xyz(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: i
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (natm <= 0) then
         call refuse(cx, "trc_set_molecule: natm must be positive, got "//trim(itoa(int(natm))))
         return
      end if
      if (multiplicity < 1) then
         call refuse(cx, "trc_set_molecule: multiplicity is 2S+1 and must be at least 1, got "// &
                     trim(itoa(int(multiplicity))))
         return
      end if
      if (allocated(cx%at_z)) deallocate (cx%at_z, cx%at_r)
      if (allocated(cx%elements)) deallocate (cx%elements)
      allocate (cx%at_z(natm), cx%at_r(3, natm), cx%elements(natm))
      do i = 1, natm
         cx%at_z(i) = real(z(i), dp)
         cx%at_r(:, i) = real(xyz(3*i - 2:3*i), dp)
      end do
      ! Back-compatible default: every atom is its own element, i.e. no
      ! ghosts. A caller wanting one calls trc_set_elements next.
      cx%elements = nint(cx%at_z)
      cx%natm = int(natm)
      cx%charge = int(charge)
      cx%multiplicity = int(multiplicity)
      cx%have_mol = .true.
      call invalidate_runs(cx)
      if (cx%df_fit_ok) call cx%dfjk%release()
      cx%df_fit_ok = .false.
      ! A new (possibly different-sized, differently-ordered) atom set
      ! invalidates any trc_set_pcm_radii override -- unlike trc_set_charge,
      ! which keeps the same atoms and so keeps it.
      if (allocated(cx%opts%pcm_radii_override)) deallocate (cx%opts%pcm_radii_override)
      status = TRC_OK
   end function trc_set_molecule

   function trc_set_elements(handle, natm, elements) result(status) bind(c, name="trc_set_elements")
      !! Place a GHOST atom: one with basis functions but no nucleus and no
      !! electrons. `elements(i)` (1..118) names the element whose basis
      !! trc_set_basis_json looks up and whose shape every atom-driven piece
      !! (the DFT/NLC grid, SAD/SADQ/SAP) follows for atom `i`; `trc_set_molecule`'s
      !! own `z(i)` stays what drives the physics (nuclear repulsion and
      !! attraction, the electron count, PCM's nuclear potential, the
      !! Hellmann-Feynman gradient terms). A ghost is `z(i) == 0` with
      !! `elements(i) >= 1`: a monomer run in its dimer's basis, say,
      !! z(i) = 0 but elements(i) names the atom whose shells still sit
      !! there. A normal atom needs no call here -- trc_set_molecule
      !! already defaults `elements` to nint(z).
      !!
      !! Call this AFTER trc_set_molecule and BEFORE trc_set_basis_json/
      !! trc_set_aux_json/trc_set_jkfit_json: those look the basis up by
      !! element, so changing it after they ran would leave the loaded
      !! shells on the wrong atoms. Refused once any of the three has run;
      !! call trc_set_molecule again (which also resets `elements`) to
      !! start over.
      !!
      !! `natm` must equal trc_set_molecule's own count, in the same atom
      !! order; terco copies `elements`, the caller's array need not
      !! outlive this call. Invalidates the last SCF/gradient/
      !! polarizability result exactly as trc_set_molecule does, resets a
      !! TRC_GUESS_GIVEN guess to SAD (the old density's shape no longer
      !! has to match) and releases any --df fit.
      !!
      !! trc_set_charge and trc_set_geometry keep whatever `elements` this
      !! call (or trc_set_molecule's own default) last set -- neither
      !! changes the atom count or which element sits where.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm
      integer(c_int), intent(in) :: elements(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: i
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_elements: set the molecule first")
         return
      end if
      if (cx%have_basis .or. cx%have_aux .or. cx%have_jkfit) then
         call refuse(cx, "trc_set_elements: a basis is already loaded -- call this before "// &
                     "trc_set_basis_*/trc_set_aux_*/trc_set_jkfit_* (JSON, libcint or arrays), not after")
         return
      end if
      status = TRC_ERR_BADARG
      if (int(natm) /= cx%natm) then
         call refuse(cx, "trc_set_elements: natm "//trim(itoa(int(natm)))// &
                     " does not match the molecule's "//trim(itoa(cx%natm))//" atoms")
         return
      end if
      do i = 1, cx%natm
         if (elements(i) < 1 .or. elements(i) > 118) then
            call refuse(cx, "trc_set_elements: element "//trim(itoa(int(elements(i))))// &
                        " at atom "//trim(itoa(i))//" is not 1..118")
            return
         end if
      end do
      cx%elements = int(elements(1:cx%natm))
      if (cx%guess_kind == TRC_GUESS_GIVEN) then
         cx%guess_kind = TRC_GUESS_SAD
         if (allocated(cx%dguess)) deallocate (cx%dguess)
      end if
      call invalidate_runs(cx)
      if (cx%df_fit_ok) call cx%dfjk%release()
      cx%df_fit_ok = .false.
      status = TRC_OK
   end function trc_set_elements

   function trc_set_charge(handle, charge, multiplicity) result(status) bind(c, name="trc_set_charge")
      !! Change the electron count on the SAME molecule, basis, auxiliary
      !! and jkfit bases -- unlike trc_set_molecule, which always
      !! releases any --df fit, this keeps `cx%dfjk`/`cx%df_fit_ok`: a
      !! fit depends on the geometry and the two bases (the metric, the
      !! three-centre tensor and its Cholesky factor -- see
      !! trc_scf_df.F90's `df_jk_setup`'s own header, "meant to survive
      !! across more than one trc_scf_run on the same geometry"), never
      !! on the electron count, so an ionisation series on one handle
      !! (N, N+1, N-1, ...) fits once and reuses it for every charge
      !! state instead of repeating the metric/Cholesky/3c/trsm for each.
      !!
      !! Same multiplicity check as trc_set_molecule (2S+1 >= 1); whether
      !! charge and multiplicity actually fit the electron count is
      !! trc_run_scf's own check, exactly as it is reached through
      !! trc_set_molecule.
      !!
      !! Invalidates the last SCF/gradient/polarizability result
      !! (invalidate_runs) and any TRC_GUESS_GIVEN density: a density
      !! from the OLD electron count's own occupation is not a guess for
      !! the new one, even where its shape happens to still match, so
      !! the guess kind resets to TRC_GUESS_SAD. Every other guess kind
      !! (CORE/GWH/SAD/SADQ/SAP) carries no state keyed to the old
      !! charge and is left as the caller set it.
      type(c_ptr), value :: handle
      integer(c_int), value :: charge, multiplicity
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_charge: set the molecule first")
         return
      end if
      status = TRC_ERR_BADARG
      if (multiplicity < 1) then
         call refuse(cx, "trc_set_charge: multiplicity is 2S+1 and must be at least 1, got "// &
                     trim(itoa(int(multiplicity))))
         return
      end if
      cx%charge = int(charge)
      cx%multiplicity = int(multiplicity)
      if (cx%guess_kind == TRC_GUESS_GIVEN) then
         cx%guess_kind = TRC_GUESS_SAD
         if (allocated(cx%dguess)) deallocate (cx%dguess)
      end if
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_charge

   subroutine take_basis(cx, aux, b)
      !! Move a freshly built basis into the context, then to the device.
      !! The copy goes up, not `b`: the device present-table is keyed on
      !! host addresses, and a basis moved up before being assigned leaves
      !! the context's arrays at addresses the table has never seen.
      type(context_box), intent(inout) :: cx
      logical, intent(in) :: aux
      type(trc_basis_t), intent(inout) :: b
      if (aux) then
         ! A new auxiliary basis touches the correlated step only; the SCF
         ! stands, which is how a driver sets it after the SCF has run.
         if (cx%have_aux) call cx%ab%b%release()
         cx%ab%b = b
         call cx%ab%b%to_device()
         cx%ab%ao_norm_kind = cx%ao_norm_kind
         cx%have_aux = .true.
         cx%rimp2_ok = .false.
         call b%release()
         return
      end if

      if (cx%have_basis) call cx%bb%b%release()
      cx%bb%b = b
      call cx%bb%b%to_device()
      cx%bb%ao_norm_kind = cx%ao_norm_kind
      if (allocated(cx%bb%cmo)) deallocate (cx%bb%cmo, cx%bb%eps)
      cx%have_basis = .true.
      call drop_eri(cx)
      call b%release()
      call invalidate_runs(cx)
      ! A new orbital basis makes any existing --df fit (built on the old
      ! one) wrong, not merely stale.
      if (cx%df_fit_ok) call cx%dfjk%release()
      cx%df_fit_ok = .false.
   end subroutine take_basis

   subroutine context_from_libcint(cx, b, atm, natm, bas, nshell, env, nenv)
      !! A context's orbital/auxiliary/jkfit basis from libcint's arrays, with
      !! the elements trc_set_elements (or trc_set_molecule's nint(z) default)
      !! recorded. libcint's atm carries only charges, 0 on a ghost, and a ghost
      !! must keep its element for the DFT grid and Becke partition, SAD/SAP
      !! and PCM radii -- without it every ghost was element 0, an error of
      !! 1e-5 to 1e-4 Eh in a ghosted B3LYP energy (water dimer, cc-pVDZ). The
      !! JSON and arrays paths already pass them. Only when the atom lists
      !! agree in length; another atom list keeps nint(charge).
      type(context_box), intent(in) :: cx
      type(trc_basis_t), intent(inout) :: b
      integer(c_int), intent(in) :: natm, nshell, nenv
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      if (allocated(cx%elements) .and. int(natm) == cx%natm) then
         call b%from_libcint(int(atm(1:6*natm)), int(natm), int(bas(1:8*nshell)), int(nshell), &
                             real(env(1:nenv), dp), at_elem=cx%elements)
      else
         call b%from_libcint(int(atm(1:6*natm)), int(natm), int(bas(1:8*nshell)), int(nshell), &
                             real(env(1:nenv), dp))
      end if
   end subroutine context_from_libcint

   function basis_libcint(handle, aux, atm, natm, bas, nshell, env, nenv, cartesian) result(status)
      !! `cartesian` nonzero builds the basis in terco's own Cartesian AO
      !! space, as before. Zero tags it spherical instead: the SHELLS
      !! libcint's `bas` describes carry no cartesian/spherical bit of
      !! their own (that is a convention the caller picks, not something
      !! `atm`/`bas`/`env` encodes), so the Gaussians `from_libcint` builds
      !! are identical either way -- what changes is only `b%spherical`,
      !! the same tag `trc_basis_from_json` sets from a JSON file's
      !! "gto_spherical" function_type. Every downstream consumer
      !! (scf_space_nao, trc_run_scf, trc_set_guess, trc_density,
      !! trc_mo_coefficients) reads that tag, not how the basis was built,
      !! so a libcint caller lands in the same spherical AO space -- same
      !! ordering, same nao_sph -- as the JSON path (trc_sph.F90's
      !! transform is the only place either path's integrals are ever
      !! folded from Cartesian to spherical).
      type(c_ptr), value :: handle
      logical, intent(in) :: aux
      integer(c_int), value :: natm, nshell, nenv, cartesian
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_basis_t) :: b
      character(len=21), parameter :: who = "trc_set_basis_libcint", who_aux = "trc_set_aux_libcint  "
      character(len=21) :: whoer
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      whoer = merge(who_aux, who, aux)
      status = TRC_ERR_BADARG
      if (natm <= 0 .or. nshell <= 0 .or. nenv <= 0) then
         call refuse(cx, trim(whoer)//": natm, nshell and nenv must all be positive")
         return
      end if
      call context_from_libcint(cx, b, atm, natm, bas, nshell, env, nenv)
      b%spherical = (cartesian == 0)
      ! TRC_AO_NORM_UNIT and a spherical basis do not combine -- same
      ! refusal basis_json makes for the JSON path (see its own comment).
      if (b%spherical .and. cx%ao_norm_kind == TRC_AO_NORM_UNIT) then
         call refuse(cx, trim(whoer)//": cartesian = 0 is a spherical basis and this "// &
                     "context is set to TRC_AO_NORM_UNIT, which only applies to a Cartesian basis -- call "// &
                     "trc_set_ao_norm(TRC_AO_NORM_SHELL) first, or pass cartesian = 1")
         call b%release()
         return
      end if
      ! Orbital/auxiliary convention mismatch, both directions -- same
      ! policy as basis_json, except a libcint caller states `cartesian`
      ! explicitly rather than through a file tag, so a mismatch is always
      ! refused here rather than coerced: nothing silently reinterprets
      ! what the caller asked for (see this module's own header, "Neither
      ! direction silently falls back to one interpretation of the two").
      if (.not. aux) then
         if (cx%have_aux .and. b%spherical .neqv. cx%ab%b%spherical) then
            call refuse(cx, who//": orbital basis is "//trim(merge("spherical", "Cartesian", b%spherical))// &
                        ", but the auxiliary basis already set is "// &
                        trim(merge("spherical", "Cartesian", cx%ab%b%spherical))// &
                        " -- match cartesian flags, or load the orbital basis first")
            call b%release()
            return
         end if
         if (cx%have_jkfit .and. b%spherical .neqv. cx%jb%b%spherical) then
            call refuse(cx, who//": orbital basis is "//trim(merge("spherical", "Cartesian", b%spherical))// &
                        ", but the --df jkfit basis already set is "// &
                        trim(merge("spherical", "Cartesian", cx%jb%b%spherical))// &
                        " -- match cartesian flags, or load the orbital basis first")
            call b%release()
            return
         end if
      else
         if (cx%have_basis .and. b%spherical .neqv. cx%bb%b%spherical) then
            call refuse(cx, trim(who_aux)//": auxiliary basis is "//trim(merge("spherical", "Cartesian", b%spherical))// &
                        ", but the orbital basis already set is "// &
                        trim(merge("spherical", "Cartesian", cx%bb%b%spherical))// &
                        " -- match cartesian flags (pyscf's auxmol.cart = mol.cart rule)")
            call b%release()
            return
         end if
      end if
      call take_basis(cx, aux, b)
      status = TRC_OK
   end function basis_libcint

   function trc_set_basis_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) result(status) &
      bind(c, name="trc_set_basis_libcint")
      !! The orbital basis from libcint's packed `atm`/`bas`/`env`.
      !! `cartesian` nonzero builds terco's own Cartesian AO space; zero
      !! tags it spherical (libcint/pyscf order, m = -l..l except p as x,
      !! y, z) -- trc_nao and every SCF-facing entry then size and index
      !! at nao_sph, exactly as for a spherical trc_set_basis_json. A
      !! mismatch against an already-set auxiliary or jkfit basis's own
      !! convention is refused, not coerced.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm, nshell, nenv, cartesian
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: status
      status = basis_libcint(handle, .false., atm, natm, bas, nshell, env, nenv, cartesian)
   end function trc_set_basis_libcint

   function trc_set_aux_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) result(status) &
      bind(c, name="trc_set_aux_libcint")
      !! The auxiliary (RI-MP2 rifit) basis, same layout and `cartesian`
      !! convention as trc_set_basis_libcint; refused rather than coerced
      !! if it disagrees with the orbital basis already set.
      type(c_ptr), value :: handle
      integer(c_int), value :: natm, nshell, nenv, cartesian
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: status
      status = basis_libcint(handle, .true., atm, natm, bas, nshell, env, nenv, cartesian)
   end function trc_set_aux_libcint

   function basis_arrays(handle, aux, nshell, maxnp, l, nprim, exps, coefs, centres) result(status)
      !! Needs the molecule first: the nuclei come from it.
      type(c_ptr), value :: handle
      logical, intent(in) :: aux
      integer(c_int), value :: nshell, maxnp
      integer(c_int), intent(in) :: l(*), nprim(*)
      real(c_double), intent(in) :: exps(*), coefs(*), centres(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_basis_t) :: b
      integer :: i
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_basis_arrays: set the molecule first; the shells need its nuclei")
         return
      end if
      status = TRC_ERR_BADARG
      if (nshell <= 0 .or. maxnp <= 0) then
         call refuse(cx, "trc_set_basis_arrays: nshell and maxnp must be positive")
         return
      end if
      do i = 1, nshell
         if (nprim(i) < 1 .or. nprim(i) > maxnp .or. l(i) < 0) then
            call refuse(cx, "trc_set_basis_arrays: shell "//trim(itoa(i))//" has nprim "// &
                        trim(itoa(int(nprim(i))))//" (1..maxnp) or l "//trim(itoa(int(l(i))))//" (>= 0)")
            return
         end if
      end do
      call b%build(int(nshell), int(l(1:nshell)), int(nprim(1:nshell)), &
                   real(reshape(exps(1:maxnp*nshell), [int(maxnp), int(nshell)]), dp), &
                   real(reshape(coefs(1:maxnp*nshell), [int(maxnp), int(nshell)]), dp), &
                   real(reshape(centres(1:3*nshell), [3, int(nshell)]), dp), &
                   cx%natm, cx%at_z, cx%at_r, int(maxnp), at_elem=cx%elements)
      call take_basis(cx, aux, b)
      status = TRC_OK
   end function basis_arrays

   function trc_set_basis_arrays(handle, nshell, maxnp, l, nprim, exps, coefs, centres) result(status) &
      bind(c, name="trc_set_basis_arrays")
      type(c_ptr), value :: handle
      integer(c_int), value :: nshell, maxnp
      integer(c_int), intent(in) :: l(*), nprim(*)
      real(c_double), intent(in) :: exps(*), coefs(*), centres(*)
      integer(c_int) :: status
      status = basis_arrays(handle, .false., nshell, maxnp, l, nprim, exps, coefs, centres)
   end function trc_set_basis_arrays

   function trc_set_aux_arrays(handle, nshell, maxnp, l, nprim, exps, coefs, centres) result(status) &
      bind(c, name="trc_set_aux_arrays")
      type(c_ptr), value :: handle
      integer(c_int), value :: nshell, maxnp
      integer(c_int), intent(in) :: l(*), nprim(*)
      real(c_double), intent(in) :: exps(*), coefs(*), centres(*)
      integer(c_int) :: status
      status = basis_arrays(handle, .true., nshell, maxnp, l, nprim, exps, coefs, centres)
   end function trc_set_aux_arrays

   function basis_json(handle, aux, path) result(status)
      !! A MolSSI BSE JSON file, read for the molecule already set.
      type(c_ptr), value :: handle
      logical, intent(in) :: aux
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_basis_t) :: b
      type(error_t) :: err
      character(len=1024) :: fpath
      integer :: i
      logical :: force_cart_this
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_basis_json: set the molecule first; the file is read for its elements")
         return
      end if
      fpath = ""
      do i = 1, len(fpath)
         if (path(i) == c_null_char) exit
         fpath(i:i) = path(i)
      end do
      ! The orbital basis (aux = .false.) follows its own file's tag
      ! unless `force_cartesian` was set (trc_set_force_cartesian). The
      ! auxiliary/fitting basis follows the ORBITAL basis's convention
      ! instead of its own tag -- pyscf's own `auxmol.cart = mol.cart`
      ! rule (density fitting and RI-MP2 now support a spherical orbital
      ! basis; see trc_run_rimp2's own header and trc_scf_df.F90) -- so it
      ! is forced Cartesian whenever the orbital basis itself is (either
      ! `force_cartesian`, or an orbital basis already loaded and NOT
      ! spherical, e.g. a basis family -- 6-31G* -- with no spherical
      ! tag of its own). An auxiliary basis loaded before any orbital one
      ! falls back to `force_cartesian` alone, same as the orbital path.
      force_cart_this = cx%opts%force_cartesian
      if (aux .and. cx%have_basis) force_cart_this = force_cart_this .or. .not. cx%bb%b%spherical
      call trc_basis_from_json(trim(fpath), cx%natm, cx%elements, cx%at_r, b, err, &
                               force_cart=force_cart_this, overrides=trc_basis_json_overrides(charges=cx%at_z))
      if (err%has_error()) then
         cx%message = err%get_message()
         status = TRC_ERR_BADARG
         return
      end if
      ! TRC_AO_NORM_UNIT and a spherical basis do not combine (see
      ! trc_set_ao_norm's own comment on why); this is the other direction
      ! of that same refusal, for a basis loaded AFTER unit mode was
      ! already requested.
      if (b%spherical .and. cx%ao_norm_kind == TRC_AO_NORM_UNIT) then
         call refuse(cx, "trc_set_basis_json: "//trim(fpath)//" is a spherical basis "// &
                     "(its JSON tags a shell 'gto_spherical') and this context is set to "// &
                     "TRC_AO_NORM_UNIT, which only applies to a Cartesian basis -- call "// &
                     "trc_set_ao_norm(TRC_AO_NORM_SHELL) first, or set force_cartesian")
         call b%release()
         return
      end if
      ! The order-dependence this guards against: `force_cart_this` above
      ! only inherits the orbital basis's convention when the ORBITAL
      ! basis is already loaded (aux .and. cx%have_basis); an auxiliary
      ! or jkfit basis set BEFORE the orbital one has nothing to inherit
      ! from and falls back to force_cartesian alone (the comment above
      ! already says so). So when the ORBITAL basis itself is what is
      ! loading here, check it against whichever fitting basis is
      ! already in the context, in BOTH directions, rather than silently
      ! running RI-MP2/--df with an orbital/auxiliary mismatch pyscf's
      ! own auxmol.cart = mol.cart rule never allows.
      if (.not. aux) then
         if (cx%have_aux .and. b%spherical .neqv. cx%ab%b%spherical) then
            call refuse(cx, "trc_set_basis_json: "//trim(fpath)//" is "// &
                        trim(merge("spherical", "Cartesian", b%spherical))//", but the auxiliary "// &
                        "basis already set is "//trim(merge("spherical", "Cartesian", cx%ab%b%spherical))// &
                        " -- it was loaded before this orbital basis, so it could not inherit this "// &
                        "basis's own convention (pyscf's auxmol.cart = mol.cart rule); set "// &
                        "force_cartesian before either, or load the orbital basis first")
            call b%release()
            return
         end if
         if (cx%have_jkfit .and. b%spherical .neqv. cx%jb%b%spherical) then
            call refuse(cx, "trc_set_basis_json: "//trim(fpath)//" is "// &
                        trim(merge("spherical", "Cartesian", b%spherical))//", but the --df jkfit "// &
                        "basis already set is "//trim(merge("spherical", "Cartesian", cx%jb%b%spherical))// &
                        " -- it was loaded before this orbital basis, so it could not inherit this "// &
                        "basis's own convention (pyscf's auxmol.cart = mol.cart rule); set "// &
                        "force_cartesian before either, or load the orbital basis first")
            call b%release()
            return
         end if
      end if
      if (aux) cx%aux_path = fpath
      call take_basis(cx, aux, b)
      status = TRC_OK
   end function basis_json

   function trc_set_basis_json(handle, path) result(status) bind(c, name="trc_set_basis_json")
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
      status = basis_json(handle, .false., path)
   end function trc_set_basis_json

   function trc_set_aux_json(handle, path) result(status) bind(c, name="trc_set_aux_json")
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
      status = basis_json(handle, .true., path)
   end function trc_set_aux_json

   function trc_set_method(handle, functional, grid_level) result(status) bind(c, name="trc_set_method")
      !! An empty functional name is Hartree-Fock.
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: functional(*)
      integer(c_int), value :: grid_level
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_xc_functional_t) :: func
      integer :: i
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (grid_level < 0) then
         call refuse(cx, "trc_set_method: grid_level "//trim(itoa(int(grid_level)))// &
                     " -- terco builds its grids by level, 0 upward; it has no explicit point counts")
         return
      end if
      cx%opts%functional = ""
      do i = 1, len(cx%opts%functional)
         if (functional(i) == c_null_char) exit
         cx%opts%functional(i:i) = functional(i)
      end do
      if (len_trim(cx%opts%functional) > 0) then
         if (.not. xc_functional_by_name(cx%opts%functional, func)) then
            status = TRC_ERR_UNSUPPORTED
            call refuse(cx, "trc_set_method: unknown functional '"//trim(cx%opts%functional)// &
                        "'; terco has "//XC_AVAILABLE_FUNCTIONALS//", or an empty name for Hartree-Fock")
            cx%opts%functional = ""
            return
         end if
         if (len_trim(func%missing) > 0) then
            status = TRC_ERR_UNSUPPORTED
            call refuse(cx, "trc_set_method: functional '"//trim(cx%opts%functional)// &
                        "' is not available: "//trim(func%missing))
            cx%opts%functional = ""
            return
         end if
      end if
      cx%opts%grid_level = int(grid_level)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_method

   function trc_set_unrestricted(handle, on) result(status) bind(c, name="trc_set_unrestricted")
      !! Force UHF/UKS (two independent spin densities) even on a CLOSED
      !! shell (nalpha == nbeta), which otherwise runs restricted --
      !! `trc_scf_options_t%unrestricted`'s own comment, "forced; an open
      !! shell is unrestricted regardless": nalpha /= nbeta already gets
      !! nspin = 2 whether or not this is set, so this only changes
      !! anything on a singlet. `on` nonzero sets it; 0 (the default)
      !! restores restricted-when-closed-shell.
      type(c_ptr), value :: handle
      integer(c_int), value :: on
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%opts%unrestricted = (on /= 0)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_unrestricted

   function trc_set_nlc_grid(handle, grid_level) result(status) bind(c, name="trc_set_nlc_grid")
      !! Grid level for VV10's nonlocal correlation (default 1, coarser than
      !! the semilocal grid: VV10 is a double sum over its points). Only
      !! used by a functional with VV10, such as wB97X-V.
      type(c_ptr), value :: handle
      integer(c_int), value :: grid_level
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (grid_level < 0) then
         call refuse(cx, "trc_set_nlc_grid: grid_level "//trim(itoa(int(grid_level)))// &
                     " -- terco builds its grids by level, 0 upward")
         return
      end if
      cx%opts%nlc_grid_level = int(grid_level)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_nlc_grid

   !
   ! A polarisable continuum around the molecule for every SCF this context
   ! runs from here on (trc_pcm): `model` is "cpcm" or "iefpcm", or an empty
   ! string to turn it off again. `eps` is the solvent's dielectric constant
   ! and is required (water: 78.3553); there is no solvent-name table, so
   ! nothing is assumed from a name. `lebedev` is the points per atomic
   ! sphere and `radii_scale` the van der Waals radius to cavity radius
   ! factor; 0 for either takes the default (302 and 1.2, pyscf's). Additive:
   ! a context that never calls this runs exactly as before.
   !
   function trc_set_pcm(handle, model, eps, lebedev, radii_scale) result(status) bind(c, name="trc_set_pcm")
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: model(*)
      real(c_double), value :: eps, radii_scale
      integer(c_int), value :: lebedev
      integer(c_int) :: status
      type(context_box), pointer :: cx
      character(len=16) :: name
      integer :: i, kind, leb
      real(dp) :: scale
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      name = ""
      do i = 1, len(name)
         if (model(i) == c_null_char) exit
         name(i:i) = model(i)
      end do
      kind = pcm_model_from_name(name)
      status = TRC_ERR_UNSUPPORTED
      if (kind < 0) then
         call refuse(cx, "trc_set_pcm: model '"//trim(name)//"' -- terco solves 'cpcm' "// &
                     "(conductor-like) and 'iefpcm' (integral equation formalism), or '' for none; "// &
                     "they are different models, so neither stands in for another")
         return
      end if
      status = TRC_ERR_BADARG
      if (kind == PCM_MODEL_NONE) then
         cx%opts%pcm_model = ""
         call invalidate_runs(cx)
         status = TRC_OK
         return
      end if
      if (eps <= 1.0_c_double) then
         call refuse(cx, "trc_set_pcm: the dielectric constant must be greater than one "// &
                     "(water is 78.3553); there is no default, which would silently pick a solvent")
         return
      end if
      leb = int(lebedev)
      if (leb == 0) leb = PCM_DEFAULT_LEBEDEV
      if (.not. any(PCM_XI_ORDERS == leb)) then
         call refuse(cx, "trc_set_pcm: "//trim(itoa(leb))//" points per sphere has no fitted SWIG "// &
                     "exponent; use 0 (302) or one of 6, 14, 26, 38, 50, 86, 110, 146, 170, 194, "// &
                     "302, 350, 434, 590, 770, 974, ... 5810")
         return
      end if
      scale = real(radii_scale, dp)
      if (scale == 0.0_dp) scale = PCM_DEFAULT_RADII_SCALE
      if (scale < 0.0_dp) then
         call refuse(cx, "trc_set_pcm: radii_scale must be positive, or 0 for the default 1.2")
         return
      end if
      cx%opts%pcm_model = trim(name)
      cx%opts%pcm_eps = real(eps, dp)
      cx%opts%pcm_lebedev = leb
      cx%opts%pcm_radii_scale = scale
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_pcm

   function trc_set_pcm_radii(handle, natm, radii) result(status) bind(c, name="trc_set_pcm_radii")
      !! Override the PCM cavity radii table (trc_pcm_radii.F90) with
      !! explicit, FINAL (already-scaled) per-atom radii in bohr,
      !! trc_set_molecule's own atom order -- for a caller whose own PCM
      !! reference used a different van der Waals table or scale than
      !! terco's default (pyscf's modified Bondi, radii_scale 1.2): mqc's
      !! own convention is plain Bondi (hydrogen 1.20 Angstrom, not
      !! terco's 1.10) times 1.2, which the two backends otherwise
      !! disagree on by construction (water cpcm/PBE0: E_diel -0.011577
      !! here against -0.010200 there, 556 against 566 kept surface
      !! points) -- not a bug in either, a different cavity. This makes
      !! them agree by construction instead.
      !!
      !! natm = 0, or a NULL radii, clears the override, back to the
      !! table and radii_scale; natm must otherwise equal the molecule's
      !! own atom count. Invalidates the last SCF result, as trc_set_pcm
      !! does. Kept across trc_set_charge (the molecule, and so the
      !! cavity's atoms, are unchanged); cleared by trc_set_molecule,
      !! which can change the atom count and order the override was
      !! given in.
      !!
      !! A GHOST atom (trc_set_elements) would get no PCM sphere, so an
      !! entry here for one is IGNORED, not refused: the override is only
      !! ever read for an atom that still has a sphere to place. (PCM with
      !! ghosts is refused at trc_run_scf for now; see elements_refusal.)
      type(c_ptr), value :: handle
      integer(c_int), value :: natm
      type(c_ptr), value :: radii
      integer(c_int) :: status
      type(context_box), pointer :: cx
      real(c_double), pointer :: rp(:)
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      if (natm == 0 .or. .not. c_associated(radii)) then
         if (allocated(cx%opts%pcm_radii_override)) deallocate (cx%opts%pcm_radii_override)
         call invalidate_runs(cx)
         status = TRC_OK
         return
      end if
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_pcm_radii: set the molecule first")
         return
      end if
      status = TRC_ERR_BADARG
      if (int(natm) /= cx%natm) then
         call refuse(cx, "trc_set_pcm_radii: natm "//trim(itoa(int(natm)))// &
                     " does not match the molecule's "//trim(itoa(cx%natm))//" atoms")
         return
      end if
      if (allocated(cx%opts%pcm_radii_override)) deallocate (cx%opts%pcm_radii_override)
      call c_f_pointer(radii, rp, [int(natm)])
      allocate (cx%opts%pcm_radii_override(cx%natm))
      cx%opts%pcm_radii_override = real(rp, dp)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_pcm_radii

   function trc_pcm_energy(handle, e_diel, q_total, npoints) result(status) bind(c, name="trc_pcm_energy")
      !! The continuum's part of the last SCF: E_diel = q.v/2 (already in
      !! trc_energy's total), the total apparent surface charge, and the
      !! surface points the cavity kept. Zeros when no continuum was set.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: e_diel, q_total
      integer(c_int), intent(out) :: npoints
      integer(c_int) :: status
      type(context_box), pointer :: cx
      e_diel = 0.0_c_double
      q_total = 0.0_c_double
      npoints = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      e_diel = real(cx%res%e_pcm, c_double)
      q_total = real(cx%res%pcm_qtotal, c_double)
      npoints = int(cx%res%pcm_npoints, c_int)
      status = TRC_OK
   end function trc_pcm_energy

   function trc_set_convergence(handle, conv_energy, conv_diis, max_iter, ndiis) result(status) &
      bind(c, name="trc_set_convergence")
      type(c_ptr), value :: handle
      real(c_double), value :: conv_energy, conv_diis
      integer(c_int), value :: max_iter, ndiis
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (conv_energy <= 0.0_c_double) then
         call refuse(cx, "trc_set_convergence: conv_energy must be positive")
         return
      end if
      if (conv_diis <= 0.0_c_double) then
         call refuse(cx, "trc_set_convergence: conv_diis must be positive (the commutator norm gate)")
         return
      end if
      if (max_iter < 1) then
         call refuse(cx, "trc_set_convergence: max_iter must be at least 1, got "//trim(itoa(int(max_iter))))
         return
      end if
      if (ndiis < 1) then
         call refuse(cx, "trc_set_convergence: ndiis must be at least 1, got "//trim(itoa(int(ndiis))))
         return
      end if
      cx%opts%conv_energy = real(conv_energy, dp)
      cx%opts%conv_diis = real(conv_diis, dp)
      cx%opts%max_iter = int(max_iter)
      cx%opts%ndiis = int(ndiis)
      status = TRC_OK
   end function trc_set_convergence

   function trc_set_lindep(handle, threshold) result(status) bind(c, name="trc_set_lindep")
      !! Canonical orthogonalisation threshold: an eigenvalue of S below
      !! `threshold` is linear dependence and its direction is dropped from
      !! the SCF (and, for a positive threshold, from the minao guess's
      !! projection). Default 1e-8. A
      !! diffuse basis on a large molecule can need more: GPU4PySCF
      !! converges a def2-TZVPD nitrosamine anion only at 5e-6 (40 functions
      !! dropped). 0 means plain symmetric orthogonalisation in the SCF,
      !! which stops when S is not numerically positive definite; the minao
      !! guess's projection always keeps a floor, so at 0 it uses its own
      !! default (1e-8) rather than none. Takes effect at the next
      !! trc_run_scf.
      type(c_ptr), value :: handle
      real(c_double), value :: threshold
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      ! Written so a NaN fails it too: every comparison with NaN is false.
      if (.not. (threshold >= 0.0_c_double .and. threshold < 1.0_c_double)) then
         call refuse(cx, "trc_set_lindep: threshold must be in [0, 1)")
         return
      end if
      cx%opts%lindep_thresh = real(threshold, dp)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_lindep

   function trc_set_screening(handle, thresh) result(status) bind(c, name="trc_set_screening")
      type(c_ptr), value :: handle
      real(c_double), value :: thresh
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (thresh <= 0.0_c_double) then
         call refuse(cx, "trc_set_screening: the threshold must be positive; 1e-12 is the default")
         return
      end if
      cx%opts%eri_thresh = real(thresh, dp)
      call drop_eri(cx)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_screening

   !
   ! Incremental Fock builds: after the first iterations, G is accumulated
   ! from the density DIFFERENCE, which screens far harder as the SCF
   ! converges (trc_scf_options_t%incremental has the measurements and the
   ! safeguards: a 1e-4 tighter screen on each difference build, a start
   ! gated on the commutator, a full rebuild every 20). 0 off, the library's
   ! default; nonzero on, the terco CLI's default.
   !
   function trc_set_incremental(handle, on) result(status) bind(c, name="trc_set_incremental")
      type(c_ptr), value :: handle
      integer(c_int), value :: on
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%opts%incremental = on /= 0
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_incremental

   function trc_set_level_shift(handle, shift) result(status) bind(c, name="trc_set_level_shift")
      !! The virtual-orbital level shift, hartree; 0 turns it off.
      !! trc_scf_options_t%level_shift's own comment: F += mu (S - S D S
      !! / occ), applied ONLY WHILE DIIS (AND ADIIS) ARE OFF -- the
      !! damped start before either has engaged (trc_scf_driver.F90's own
      !! "level shift while DIIS and A/EDIIS are off" block). It is a
      !! CONVERGENCE AID for that early phase, not a persistent shift
      !! held through the whole SCF: once the commutator is small enough
      !! for CDIIS (or ADIIS takes over first), this stops being applied
      !! at all, however large `shift` is set. A caller wanting a shift
      !! that stays on throughout -- the sense a persistent level shift
      !! has in, for instance, a Fukui-ion deck's convergence recipe --
      !! is asking for something this entry does not do; nothing here
      !! changes that behaviour, only exposes the existing knob.
      type(c_ptr), value :: handle
      real(c_double), value :: shift
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (shift < 0.0_c_double) then
         call refuse(cx, "trc_set_level_shift: shift must be >= 0 (0 turns it off)")
         return
      end if
      cx%opts%level_shift = real(shift, dp)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_level_shift

   !
   ! Density-fitted J and K (`--df`). Needs a jkfit basis too
   ! (`trc_set_jkfit_json`, separate from `trc_set_aux_json`'s rifit set
   ! because one run can do DF-SCF and RI-MP2 with two different auxiliary
   ! bases) before `trc_run_scf` -- `trc_run_scf` reports that, not this,
   ! since nothing here can tell whether one was already set by the time
   ! this is called.
   !
   function trc_set_df(handle, on) result(status) bind(c, name="trc_set_df")
      type(c_ptr), value :: handle
      integer(c_int), value :: on
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%opts%df = on /= 0
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_df

   !
   ! The fitting (jkfit) basis for `--df`, a MolSSI BSE JSON file, read for
   ! the molecule already set -- the same file format `trc_set_aux_json`
   ! reads, into a SEPARATE slot (`cx%jb`). Follows the orbital basis's
   ! own convention, exactly `basis_json`'s rule for the rifit (`aux`)
   ! set: Cartesian when the orbital basis is (`force_cartesian`, or an
   ! orbital basis already loaded that is not itself spherical), its own
   ! JSON tag otherwise (def2-universal-jkfit is spherical, matching a
   ! spherical orbital basis -- `trc_scf_df.F90`'s `df_jk_setup` folds it,
   ! see that module's own header). Clears any fit already built on the
   ! previous jkfit basis.
   !
   function trc_set_jkfit_json(handle, path) result(status) bind(c, name="trc_set_jkfit_json")
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_basis_t) :: b
      type(error_t) :: err
      character(len=1024) :: fpath
      integer :: i
      logical :: force_cart_this
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_jkfit_json: set the molecule first; the file is read for its elements")
         return
      end if
      fpath = ""
      do i = 1, len(fpath)
         if (path(i) == c_null_char) exit
         fpath(i:i) = path(i)
      end do
      force_cart_this = cx%opts%force_cartesian
      if (cx%have_basis) force_cart_this = force_cart_this .or. .not. cx%bb%b%spherical
      call trc_basis_from_json(trim(fpath), cx%natm, cx%elements, cx%at_r, b, err, &
                               force_cart=force_cart_this, overrides=trc_basis_json_overrides(charges=cx%at_z))
      if (err%has_error()) then
         cx%message = err%get_message()
         status = TRC_ERR_BADARG
         return
      end if
      if (cx%have_jkfit) call cx%jb%b%release()
      cx%jb%b = b
      call cx%jb%b%to_device()
      cx%have_jkfit = .true.
      cx%jkfit_path = fpath
      call b%release()
      if (cx%df_fit_ok) call cx%dfjk%release()
      cx%df_fit_ok = .false.
      status = TRC_OK
   end function trc_set_jkfit_json

   !
   ! The same fitting basis as `trc_set_jkfit_json`, from libcint's packed
   ! `atm`/`bas`/`env` instead of a JSON file -- the same shape
   ! `trc_set_aux_libcint` takes. `cartesian` zero tags it spherical,
   ! exactly as `basis_libcint` above; refused rather than coerced if it
   ! disagrees with the orbital basis already set (pyscf's `auxmol.cart =
   ! mol.cart`), since there is no file tag here to silently override.
   !
   function trc_set_jkfit_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) result(status) &
      bind(c, name="trc_set_jkfit_libcint")
      !! The --df fitting basis from libcint's packed arrays, the same
      !! shape and `cartesian` convention as trc_set_aux_libcint; refused
      !! rather than coerced if it disagrees with the orbital basis
      !! already set (pyscf's auxmol.cart = mol.cart rule).
      type(c_ptr), value :: handle
      integer(c_int), value :: natm, nshell, nenv, cartesian
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_basis_t) :: b
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (natm <= 0 .or. nshell <= 0 .or. nenv <= 0) then
         call refuse(cx, "trc_set_jkfit_libcint: natm, nshell and nenv must all be positive")
         return
      end if
      call context_from_libcint(cx, b, atm, natm, bas, nshell, env, nenv)
      b%spherical = (cartesian == 0)
      if (cx%have_basis .and. b%spherical .neqv. cx%bb%b%spherical) then
         call refuse(cx, "trc_set_jkfit_libcint: jkfit basis is "// &
                     trim(merge("spherical", "Cartesian", b%spherical))//", but the orbital basis "// &
                     "already set is "//trim(merge("spherical", "Cartesian", cx%bb%b%spherical))// &
                     " -- match cartesian flags (pyscf's auxmol.cart = mol.cart rule)")
         call b%release()
         return
      end if
      if (cx%have_jkfit) call cx%jb%b%release()
      cx%jb%b = b
      call cx%jb%b%to_device()
      cx%have_jkfit = .true.
      cx%jkfit_path = ""
      call b%release()
      if (cx%df_fit_ok) call cx%dfjk%release()
      cx%df_fit_ok = .false.
      status = TRC_OK
   end function trc_set_jkfit_libcint

   !
   ! How finely the shell pairs are batched by their Schwarz bound: `res`
   ! buckets per decade. The pre-launch screen admits a pair of batches
   ! when their bounds multiply above the threshold, so a coarse bucket
   ! admits everything within its decade and leaves the rest to the test
   ! inside the kernel, where a rejected quartet has already cost a thread.
   ! Finer buckets catch more before the launch, at more and smaller
   ! batches. 1 is decades and is the default; 4 is the finest the bins are
   ! sized for. Rebuilds the ERI object, as a threshold change does.
   !
   function trc_set_batching(handle, res) result(status) bind(c, name="trc_set_batching")
      type(c_ptr), value :: handle
      integer(c_int), value :: res
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (res < 1 .or. res > 4) then
         call refuse(cx, "trc_set_batching: res is buckets per decade, 1 to 4; got "//trim(itoa(int(res))))
         return
      end if
      cx%opts%batch_res = int(res)
      call drop_eri(cx)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_batching

   function trc_set_force_cartesian(handle, flag) result(status) bind(c, name="trc_set_force_cartesian")
      !! Force the Cartesian AO space on every basis this context loads
      !! from here on, regardless of what a JSON file's function_type says
      !! (trc_scf_options_t%force_cartesian). `flag` non-zero sets it.
      !!
      !! NOT RETROACTIVE: a basis already loaded (trc_set_basis_json) kept
      !! whichever space it was built in at that call -- the same way
      !! trc_set_ao_norm's own comment describes for a TRC_GUESS_GIVEN
      !! density already copied. Call this BEFORE trc_set_basis_json/
      !! trc_set_aux_json for it to take effect on that basis; a caller
      !! wanting it to also apply to a basis already loaded has to load
      !! the file again.
      type(c_ptr), value :: handle
      integer(c_int), value :: flag
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%opts%force_cartesian = (flag /= 0)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_force_cartesian

   function trc_set_guess(handle, kind, dguess, nspin) result(status) bind(c, name="trc_set_guess")
      !! `kind` one of TRC_GUESS_*. For TRC_GUESS_GIVEN, `dguess` is
      !! (nao, nao, nspin) column-major and is COPIED here; the caller's array
      !! is free afterwards. The basis must be set first for that copy.
      !!
      !! CONTRAVARIANT under trc_set_ao_norm: a density given in
      !! TRC_AO_NORM_UNIT is converted to terco's own shell-normalised
      !! convention here, once, at the copy -- D_terco = N D_given N (see
      !! trc_ao_norm.F90's header; this is the "guess coming IN" direction,
      !! the covariant congruence run forward). Everything downstream of
      !! this copy, including trc_run_scf's own iteration, only ever sees
      !! terco's native convention.
      type(c_ptr), value :: handle
      integer(c_int), value :: kind, nspin
      type(c_ptr), value :: dguess
      integer(c_int) :: status
      type(context_box), pointer :: cx
      real(c_double), pointer :: dg(:, :, :)
      real(dp), allocatable :: nrm(:)
      integer :: n, s
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (kind < TRC_GUESS_CORE .or. kind > TRC_GUESS_MINAO) then
         call refuse(cx, "trc_set_guess: kind "//trim(itoa(int(kind))) &
                     //" -- 0 core, 1 GWH, 2 SAD, 3 given, 4 SADQ, 5 SAP, 6 projected, 7 minao")
         return
      end if
      if (allocated(cx%dguess)) deallocate (cx%dguess)
      if (kind == TRC_GUESS_GIVEN) then
         if (.not. c_associated(dguess)) then
            call refuse(cx, "trc_set_guess: kind 3 (given) needs a density pointer")
            return
         end if
         if (nspin < 1 .or. nspin > 2) then
            call refuse(cx, "trc_set_guess: nspin must be 1 or 2, got "//trim(itoa(int(nspin))))
            return
         end if
         status = TRC_ERR_STATE
         if (.not. cx%have_basis) then
            call refuse(cx, "trc_set_guess: set the basis before giving a density; its size comes from there")
            return
         end if
         ! The SCF's own space (nao_sph if this basis is spherical) -- a
         ! caller supplying a guess density sizes it from `trc_nao`, which
         ! now reports the same number (scf_space_nao).
         n = scf_space_nao(cx%bb)
         call c_f_pointer(dguess, dg, [n, n, int(nspin)])
         allocate (cx%dguess(n, n, nspin))
         cx%dguess = real(dg, dp)
         if (cx%bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
            allocate (nrm(n))
            call ao_norm_vector(cx%bb%b%nshell, cx%bb%b%sh_l, cx%bb%b%sh_ao, n, cx%bb%ao_norm_kind, nrm)
            do s = 1, int(nspin)
               call ao_scale_cov(nrm, cx%dguess(:, :, s))
            end do
            deallocate (nrm)
         end if
      end if
      cx%guess_kind = int(kind)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_guess

   function trc_set_ao_norm(handle, kind) result(status) bind(c, name="trc_set_ao_norm")
      !! How this context's AO-indexed C API entries normalise a Cartesian
      !! shell's components: TRC_AO_NORM_SHELL (0, the default -- terco's
      !! own convention, and libcint's) or TRC_AO_NORM_UNIT (1 -- every
      !! component individually normalised to 1, GAMESS's convention; see
      !! trc_ao_norm.F90's header for the closed form and which entries
      !! convert). Defaulting to SHELL means no caller that never calls
      !! this sees any difference from before this entry existed.
      !!
      !! Order-independent of trc_set_basis* and trc_set_aux*: settable before
      !! or after either, and takes effect on whichever basis is present now
      !! and on any set afterwards. NOT retroactive on an already-copied
      !! TRC_GUESS_GIVEN density, though -- that conversion happened once,
      !! at trc_set_guess's own call, under whatever mode was active then;
      !! call this before trc_set_guess if the two must agree.
      type(c_ptr), value :: handle
      integer(c_int), value :: kind
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (kind /= TRC_AO_NORM_SHELL .and. kind /= TRC_AO_NORM_UNIT) then
         call refuse(cx, "trc_set_ao_norm: kind must be 0 (shell) or 1 (unit)")
         return
      end if
      ! TRC_AO_NORM_UNIT rescales CARTESIAN components individually
      ! (GAMESS's convention); it says nothing about a SPHERICAL AO, which
      ! has no single Cartesian component to rescale -- a d_z2 is already a
      ! fixed linear combination of six of them, each carrying its own
      ! ratio. Combining the two would need a second, different
      ! normalisation convention for spherical AOs that nothing here
      ! defines, so this is refused rather than silently doing one of the
      ! two things it could plausibly mean. Cartesian bases (the default
      ! for every basis without a "gto_spherical" tag, or any basis loaded
      ! with force_cartesian) are unaffected.
      if (kind == TRC_AO_NORM_UNIT) then
         if ((cx%have_basis .and. cx%bb%b%spherical) .or. (cx%have_aux .and. cx%ab%b%spherical)) then
            call refuse(cx, "trc_set_ao_norm: TRC_AO_NORM_UNIT only applies to a Cartesian "// &
                        "basis; this context's basis is spherical (its JSON basis was "// &
                        "'gto_spherical' and force_cartesian was not set) -- set "// &
                        "force_cartesian, or load a Cartesian basis, before asking for UNIT")
            return
         end if
      end if
      cx%ao_norm_kind = int(kind)
      if (cx%have_basis) cx%bb%ao_norm_kind = cx%ao_norm_kind
      if (cx%have_aux) cx%ab%ao_norm_kind = cx%ao_norm_kind
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_ao_norm

   function trc_get_ao_norm(handle, kind) result(status) bind(c, name="trc_get_ao_norm")
      !! The AO normalisation this context is set to -- TRC_AO_NORM_SHELL
      !! unless trc_set_ao_norm was called with TRC_AO_NORM_UNIT.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: kind
      integer(c_int) :: status
      type(context_box), pointer :: cx
      kind = TRC_AO_NORM_SHELL
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      kind = int(cx%ao_norm_kind, c_int)
      status = TRC_OK
   end function trc_get_ao_norm

   function trc_set_comm(handle, fcomm) result(status) bind(c, name="trc_set_comm")
      !! `fcomm` is the Fortran handle of ANY intra-communicator the caller
      !! already has open and collective over -- MPI_Comm_c2f from C, or
      !! the integer a Fortran caller holds (mpi_f08's MPI_VAL, or the
      !! legacy `mpi` module's handle directly); -1 means one rank, no MPI.
      !! It need not be MPI_COMM_WORLD: GAMESS, for example, hands this
      !! DDI's compute communicator, whose data-server ranks are not part
      !! of it and never call trc_run_scf/trc_run_rimp2 at all. terco
      !! duplicates whatever `fcomm` names (MPI_Comm_dup) the first time a
      !! run needs it, so every rank of THAT communicator must reach the
      !! run together; the caller keeps ownership of its own communicator,
      !! since only terco's duplicate is ever freed.
      type(c_ptr), value :: handle
      integer(c_int), value :: fcomm
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%fcomm = int(fcomm)
      status = TRC_OK
   end function trc_set_comm

   function trc_set_verbose(handle, level) result(status) bind(c, name="trc_set_verbose")
      type(c_ptr), value :: handle
      integer(c_int), value :: level
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%opts%verbose = level > 0
      status = TRC_OK
   end function trc_set_verbose

   function trc_set_rimp2(handle, nfrozen, aux_block) result(status) bind(c, name="trc_set_rimp2")
      !! `nfrozen` -1 counts the core from the elements; `aux_block` <= 0 is
      !! the whole auxiliary basis at once.
      type(c_ptr), value :: handle
      integer(c_int), value :: nfrozen, aux_block
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%nfrozen = int(nfrozen)
      cx%aux_block = int(aux_block)
      cx%rimp2_ok = .false.
      status = TRC_OK
   end function trc_set_rimp2

   subroutine trc_build_comm(fcomm, comm, ok, status)
      !! THE ONE PLACE `fcomm` (a caller's Fortran communicator handle, as
      !! every C entry point below receives it) becomes terco's own comm_t.
      !! Every entry that takes or stores an `fcomm` -- capi_scf_mpi,
      !! capi_rimp2, and trc_run_scf/trc_run_rimp2/trc_context_eri through
      !! the context's cx%fcomm -- calls this rather than repeating the
      !! conversion.
      !!
      !! `fcomm` of -1 means one rank, no MPI: `ok` comes back .false. and
      !! `comm` is not built, so the caller takes its non-collective path.
      !! Any other value is the Fortran handle of an intra-communicator the
      !! caller already has open: MPI_Comm_c2f from C, the MPI_VAL of an
      !! mpi_f08 handle, or the plain integer the legacy `mpi` module uses.
      !! ANY such communicator is accepted -- there is no MPI_COMM_WORLD
      !! check here, deliberately, because GAMESS's DDI hands libraries its
      !! *compute* communicator, not MPI_COMM_WORLD: DDI's data-server
      !! ranks sit in MPI_COMM_WORLD but never call into terco, so waiting
      !! on MPI_COMM_WORLD from inside terco would hang on the first
      !! collective. A fragment method's per-fragment sub-communicator is
      !! the same story. This IS collective, over whatever communicator
      !! `fcomm` names: comm_from_handle duplicates it (MPI_Comm_dup), so
      !! every rank of that communicator must call this together.
      !!
      !! The duplicate is what makes ownership work: terco's `comm` is its
      !! own copy, so `comm%finalize()` after the run frees only that copy.
      !! The caller's own handle is never touched and stays open for
      !! whatever else it is used for (other DDI traffic, another library)
      !! during and after this call.
      !!
      !! A communicator of size 1 -- a one-rank fragment -- needs no special
      !! case: MPI_Comm_dup and the node split below are both defined, and
      !! well-behaved, on a single-rank communicator.
      !!
      !! Binds this rank's device by NODE-LOCAL rank, not by comm%rank():
      !! a DDI compute communicator or a fragment sub-communicator is not
      !! guaranteed to be node-ordered (rank 0 and rank 1 of it need not be
      !! on the same node, or even adjacent on one), so binding by comm%rank()
      !! directly could stack two ranks of one node on the same device
      !! while another device sat idle. Splitting by shared memory
      !! (comm%split(), MPI_COMM_TYPE_SHARED under the hood) finds the
      !! ranks that actually share a node; the node communicator is scratch
      !! and is freed here, before returning.
      integer(c_int), intent(in) :: fcomm
      type(comm_t), intent(out) :: comm
      logical, intent(out) :: ok
      integer(c_int), intent(out) :: status
      type(comm_t) :: node_comm
      integer :: dev
      ok = .false.
      status = TRC_OK
      if (fcomm == -1) return
      comm = comm_from_handle(int(fcomm))
      if (comm%is_null()) then
         status = TRC_ERR_BADARG
         return
      end if
      node_comm = comm%split()
      dev = trc_bind_device(node_comm%rank())
      call node_comm%finalize()
      ok = .true.
   end subroutine trc_build_comm

   function trc_run_scf(handle) result(status) bind(c, name="trc_run_scf")
      !! Needs the molecule and the basis. Builds the guess it was asked for,
      !! runs, keeps everything in the context.
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_scf_options_t) :: opts
      type(comm_t) :: comm
      type(error_t) :: err
      real(dp), allocatable :: dsad(:, :), dg(:, :, :), hsap(:, :)
      integer :: nelec, nunp, nalpha, nbeta, nspin, n, nc
      integer :: nocc(2)
      logical :: collective, talk, use_projected
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%scf_ok = .false.
      cx%rimp2_ok = .false.
      cx%grad_ok = .false.
      status = TRC_ERR_STATE
      if (.not. (cx%have_mol .and. cx%have_basis)) then
         cx%message = "trc_run_scf: set the molecule and the basis first"
         return
      end if
      status = TRC_ERR_BADARG
      if (len(elements_refusal(cx)) > 0) then
         cx%message = "trc_run_scf: "//elements_refusal(cx)
         return
      end if
      ! `n` is the SCF's own space (nao_sph if spherical); `nc` is always
      ! the Cartesian count, for the one place below (SAP) that still
      ! builds something in Cartesian space before handing it to the SCF.
      n = scf_space_nao(cx%bb)
      nc = cx%bb%b%nao
      nelec = nint(sum(cx%at_z)) - cx%charge
      nunp = cx%multiplicity - 1
      status = TRC_ERR_BADARG
      if (nelec <= 0 .or. nunp > nelec .or. mod(nelec - nunp, 2) /= 0) then
         cx%message = "trc_run_scf: charge and multiplicity do not fit the electron count "// &
                      "(must be a positive number of electrons)"
         return
      end if
      nalpha = (nelec + nunp)/2
      nbeta = nelec - nalpha
      if (nalpha > n) then
         cx%message = "trc_run_scf: more electrons than functions"
         return
      end if
      ! A shell above this build's four-centre kernels has no kernel: its
      ! integrals came out as NaN and the run died much later in the
      ! eigensolver. Checked even under --df, whose SAD guess still runs
      ! four-centre atomic SCFs; trc_max_angular_momentum gives the
      ! ceilings up front.
      if (.not. cx%opts%df .or. cx%guess_kind == TRC_GUESS_SAD .or. cx%guess_kind == TRC_GUESS_SADQ) then
         block
            character(len=:), allocatable :: refusal
            refusal = eri_lmax_refusal(cx%bb%b)
            if (len(refusal) > 0) then
               status = TRC_ERR_UNSUPPORTED
               cx%message = "trc_run_scf: "//refusal
               return
            end if
         end block
      end if
      call scf_spin_rule(nalpha, nbeta, cx%opts%unrestricted, cx%opts%frac_occ, cx%opts%nelec_frac, nspin, nocc)
      ! The communicator first, so the guess below binds the same device
      ! the SCF will and says what it did once, on rank 0, not once per
      ! rank. Every `return` from here on frees it.
      call trc_build_comm(cx%fcomm, comm, collective, status)
      if (status /= TRC_OK) return
      status = TRC_ERR_BADARG   ! what the guess's own failures below report
      talk = cx%opts%verbose
      if (collective) talk = talk .and. comm%rank() == 0
      ! Once per context; a later run only grows the workspace to n.
      call cx%la%init(n)
      ! The sizes first, before any guess runs: a SAD guess on a large
      ! system takes long enough that a wrong basis file should show now,
      ! not after it.
      if (talk) then
         call global_logger%info("   basis functions  : "//to_string(scf_space_nao(cx%bb))// &
                                 merge(" (spherical)", " (Cartesian)", cx%bb%b%spherical))
         if (cx%have_jkfit .and. cx%opts%df) then
            call global_logger%info("   auxiliary (JK)   : "//to_string(scf_space_nao(cx%jb))//" functions")
         end if
         if (cx%have_aux) then
            call global_logger%info("   auxiliary (MP2)  : "//to_string(scf_space_nao(cx%ab))//" functions")
         end if
      end if
      ! --- the guess, BEFORE the --df fit below -----------------------------
      ! SAD/SADQ's atomic sub-SCFs (trc_sad_build/trc_sadq_build) run the
      ! ordinary four-centre `eri%fock_resident` path regardless of
      ! `cx%opts%df` -- SAD never fits its own tiny aux basis, by design (see
      ! trc_sad.F90's own header). The FIRST four-centre kernel launch of the
      ! process -- here, the SAD atom's, however tiny its own basis is --
      ! pays a one-time cost that is never given back for the rest of the
      ! process: measured on a V100, ~4.7 GiB when it runs with ~24 GiB
      ! already committed (the order this file had before), ~8.9 GiB when it
      ! runs first against ~31 GiB free. It is NOT a fixed, problem-size-
      ! independent reservation -- its size tracks how much is free when it
      ! fires, which is exactly why paying it HERE, first, matters: every
      ! allocation after it (`cx%dfjk%setup`'s B, and its own free-memory-
      ! aware block sizing in `choose_setup_block`/`choose_qb`) then sees the
      ! device's true remaining budget instead of one this cost is about to
      ! shrink out from under it. Measured consequence at w40 cc-pVDZ --df
      ! (1000 functions, naux 10640, B alone ~19.8 GiB): with the fit BEFORE
      ! the guess (the old order), `choose_setup_block` sized its (mu nu|P)
      ! scratch against free memory that had not yet paid this cost, and the
      ! run died with CUDA_ERROR_OUT_OF_MEMORY inside that scratch's own
      ! allocation; with the fit AFTER (this order), the same run converges
      ! in 12 iterations with free memory flat at ~1.09 GiB from the first
      ! Fock build on. The two blocks set disjoint state (`dg`/`hsap`/
      ! `opts%guess` here, `cx%dfjk`/`cx%df_fit_ok` below) and neither reads
      ! the other's, so swapping their order changes nothing but the order
      ! memory gets claimed in.
      opts = cx%opts
      use_projected = .false.
      select case (cx%guess_kind)
      case (TRC_GUESS_CORE)
         opts%guess = "core"
      case (TRC_GUESS_GWH)
         opts%guess = "gwh"
      case (TRC_GUESS_SAD, TRC_GUESS_SADQ)
         ! SADQ is SAD with the molecule's charge put on its atoms, so it
         ! needs the electron count this run was asked for -- which is the
         ! one thing SAD never looks at. On a neutral molecule the two are
         ! the same density; see trc_sad.
         if (cx%guess_kind == TRC_GUESS_SADQ) then
            call trc_sadq_build(cx%bb%b, nelec, dsad, err, verbose=talk, la=cx%la)
         else
            call trc_sad_build(cx%bb%b, dsad, err, verbose=talk, la=cx%la)
         end if
         if (err%has_error()) then
            cx%message = err%get_message()
            if (collective) call comm%finalize()
            return
         end if
         allocate (dg(n, n, nspin))
         if (nspin == 1) then
            dg(:, :, 1) = dsad
         else
            dg(:, :, 1) = 0.5_dp*dsad
            dg(:, :, 2) = 0.5_dp*dsad
         end if
      case (TRC_GUESS_SAP)
         ! SAP hands back a one-electron operator, not a density, so it goes
         ! in through `hguess` below rather than `dg`. T and V_ne are what it
         ! needs and the context does not keep them, so they are built here;
         ! it is one pair list and one 1e pass against an SCF.
         !
         ! trc_sap_build itself stays entirely Cartesian (`hsap_c`, `nc`):
         ! it is built from the analytic one-electron matrices and a
         ! quadrature over free-atom potentials, neither of which this
         ! feature touches. The one transform this guess needs is here, at
         ! the boundary, the same as SAD's own (trc_sad.F90) -- except SAD
         ! produces its density already in the matching space because its
         ! atomic sub-SCFs run on a `subset` basis that carries `spherical`
         ! itself, where SAP's grid/one-electron machinery has no such
         ! basis-level hook to carry it through.
         block
            type(trc_pairlist_t) :: pl
            real(dp), allocatable :: s1(:, :), t1(:, :), v1(:, :), hsap_c(:, :)
            call pl%build(cx%bb%b, cx%opts%eri_thresh)
            call pl%to_device()
            allocate (s1(nc, nc), t1(nc, nc), v1(nc, nc), hsap_c(nc, nc))
            !$acc enter data create(s1, t1, v1)
            call trc_1e(cx%bb%b, pl, s1, t1, v1)
            !$acc update self(s1, t1, v1)
            !$acc exit data delete(s1, t1, v1)
            call trc_sap_build(cx%bb%b, pl, t1, v1, hsap_c, err, verbose=talk)
            call pl%release()
            deallocate (s1, t1, v1)
            if (.not. err%has_error()) then
               allocate (hsap(n, n))
               if (cx%bb%b%spherical) then
                  ! trc_sap_build already copied `hsap_c` back to the host
                  ! and dropped its own device copy, so it has to go back
                  ! up for the transform -- which, like every other matrix
                  ! this SCF touches, runs on the device.
                  !$acc enter data copyin(hsap_c) create(hsap)
                  call sph_cart_to_sph(cx%bb%b, hsap_c, hsap)
                  !$acc update self(hsap)
                  !$acc exit data delete(hsap_c, hsap)
               else
                  hsap = hsap_c
               end if
            end if
            deallocate (hsap_c)
         end block
         if (err%has_error()) then
            cx%message = err%get_message()
            if (collective) call comm%finalize()
            return
         end if
      case (TRC_GUESS_MINAO)
         ! No free-atom SCF and nothing cached (unlike SAD/SADQ above): this
         ! is a single cross-basis projection -- minao_guess_build, below,
         ! split out of this dispatch to keep it off trc_run_scf's own line
         ! count.
         call minao_guess_build(cx%bb%b, n, nspin, opts%lindep_thresh, opts%eri_thresh, dg, err)
         if (err%has_error()) then
            cx%message = err%get_message()
            if (collective) call comm%finalize()
            return
         end if
      case (TRC_GUESS_GIVEN)
         if (.not. allocated(cx%dguess)) then
            cx%message = "trc_run_scf: guess GIVEN but no density was given"
            status = TRC_ERR_STATE
            if (collective) call comm%finalize()
            return
         end if
         if (size(cx%dguess, 3) /= nspin) then
            cx%message = "trc_run_scf: the given density has the wrong spin count for this molecule"
            if (collective) call comm%finalize()
            return
         end if
         dg = cx%dguess
      case (TRC_GUESS_PROJECTED)
         ! `dg` stays unallocated here -- the new overlap this guess needs
         ! to project into only exists once `trc_scf_run`'s own
         ! `scf_context%setup` has built it, so the projection itself runs
         ! INSIDE that call (its own `cmo_project` argument), not here.
         ! This dispatch only validates that there IS a previous orbital
         ! set, at THIS run's own size, to project: `cx%bb%cmo` is stored
         ! by any run with at least one iteration (right after the
         ! trc_scf_run call below, gated on `iterations == 0`, not on
         ! `converged`), so "converged" below would overstate what is
         ! actually required.
         if (.not. allocated(cx%bb%cmo)) then
            cx%message = "trc_run_scf: guess projected but no previous SCF has run yet"
            status = TRC_ERR_STATE
            if (collective) call comm%finalize()
            return
         end if
         if (cx%opts%frac_occ) then
            cx%message = "trc_run_scf: guess projected does not support fractional occupation "// &
                         "(project_density's occupations are integer)"
            status = TRC_ERR_STATE
            if (collective) call comm%finalize()
            return
         end if
         if (size(cx%bb%cmo, 1) /= n .or. size(cx%bb%cmo, 3) /= nspin) then
            cx%message = "trc_run_scf: guess projected needs a previous SCF at the SAME size "// &
                         "and spin as this one (basis/molecule/unrestricted must not have changed)"
            status = TRC_ERR_STATE
            if (collective) call comm%finalize()
            return
         end if
         use_projected = .true.
      case default
         cx%message = "trc_run_scf: unknown guess kind"
         if (collective) call comm%finalize()
         return
      end select
      ! The guess just paid (or re-confirmed already-paid) the SAD atomic
      ! sub-SCF's own one-time CUDA local-memory reservation -- see
      ! trc_release_local_mem_reservation's own header. Nothing left in a
      ! --df run ever launches that four-centre class again (build_g skips
      ! eri%fock_resident entirely under --df), so give it back before `the
      ! fit` below sizes its own scratch against free memory that cost is
      ! about to vacate: that ordering is this branch's own fix for the w40
      ! OOM, and shrinking the reservation here is what makes the memory it
      ! freed actually available rather than merely counted as "free" by a
      ! query that cannot tell the difference.
      if (cx%opts%df) call trc_release_local_mem_reservation()
      ! --- --df: fit once per (molecule, basis, jkfit), carried across runs
      if (cx%opts%df .and. .not. cx%df_fit_ok) then
         if (.not. cx%have_jkfit) then
            cx%message = "trc_run_scf: --df needs a jkfit basis (trc_set_jkfit_json)"
            if (collective) call comm%finalize()
            return
         end if
         block
            type(trc_pairlist_t) :: pl
            type(error_t) :: dferr
            call pl%build(cx%bb%b, cx%opts%eri_thresh)
            call pl%to_device()
            call cx%dfjk%setup(cx%bb%b, cx%jb%b, pl, cx%la, dferr, nocc_hint=max(nalpha, 1))
            call pl%release()
            if (dferr%has_error()) then
               cx%message = dferr%get_message()
               if (collective) call comm%finalize()
               return
            end if
         end block
         cx%df_fit_ok = .true.
      end if
      if (use_projected) then
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, comm=comm, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri, cmo_project=cx%bb%cmo)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri, cmo_project=cx%bb%cmo)
         end if
      else if (allocated(dg)) then
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, dguess=dg, comm=comm, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, dguess=dg, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri)
         end if
      else if (allocated(hsap)) then
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, comm=comm, hguess=hsap, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, hguess=hsap, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri)
         end if
      else
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, comm=comm, la=cx%la, scf_dfjk=cx%dfjk, &
                             grid_cache=cx%kgrid, pair_eri_cache=cx%peri)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, la=cx%la, scf_dfjk=cx%dfjk, grid_cache=cx%kgrid, &
                             pair_eri_cache=cx%peri)
         end if
      end if
      if (collective) call comm%finalize()
      cx%message = cx%res%message
      if (len_trim(cx%res%message) > 0 .and. cx%res%iterations == 0) then
         status = TRC_ERR_UNSUPPORTED
         return
      end if
      ! the orbitals, where the old RI-MP2 entry also looks for them
      if (allocated(cx%bb%cmo)) deallocate (cx%bb%cmo, cx%bb%eps)
      allocate (cx%bb%cmo(n, n, nspin), cx%bb%eps(n, nspin))
      cx%bb%cmo = cx%res%cmo
      cx%bb%eps = cx%res%eps
      cx%bb%nalpha = nalpha
      cx%bb%nbeta = nbeta
      cx%nalpha = nalpha
      cx%nbeta = nbeta
      cx%scf_ok = .true.
      ! The settings `cx%kgrid` was just built at (empty/false if `opts%dft`
      ! was false), snapshotted for `trc_run_gradient`'s own staleness
      ! check -- see `context_box%kgrid`'s own header.
      cx%kgrid_level = opts%grid_level
      cx%kgrid_max_pts = opts%grid_max_pts
      cx%kgrid_rho_tol = opts%rho_tol
      cx%kgrid_nlc_level = opts%nlc_grid_level
      status = merge(TRC_OK, TRC_ERR_NOCONV, cx%res%converged)
   end function trc_run_scf

   !
   ! TRC_GUESS_MINAO's density, both spins -- split out of trc_run_scf's own
   ! dispatch purely to keep that procedure's own line count down; every
   ! other guess kind's block stays inline there.
   !
   subroutine minao_guess_build(b, n, nspin, lindep_thresh, eri_thresh, dg, err)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: n, nspin
      real(dp), intent(in) :: lindep_thresh, eri_thresh
      real(dp), allocatable, intent(out) :: dg(:, :, :)
      type(error_t), intent(inout) :: err
      real(dp), allocatable :: dminao(:, :)
      call trc_minao_build(b, dminao, err, threshold=lindep_thresh, eri_thresh=eri_thresh)
      if (err%has_error()) return
      allocate (dg(n, n, nspin))
      if (nspin == 1) then
         dg(:, :, 1) = dminao
      else
         dg(:, :, 1) = 0.5_dp*dminao
         dg(:, :, 2) = 0.5_dp*dminao
      end if
   end subroutine minao_guess_build

   !
   ! What trc_run_gradient refuses, factored out so trc_run_hessian can
   ! make the SAME check up front -- before spending any SCF time on a
   ! displacement loop that would fail at the first gradient anyway. Also
   ! resolves cx%opts%functional into `func`/`have_func`, since
   ! trc_run_gradient's own grid-reuse logic (its header) needs them too
   ! and must not look the name up twice under two different outcomes.
   ! `refusal` comes back zero-length when nothing here objects.
   !
   subroutine gradient_capability_check(cx, have_func, func, refusal, refusal_status)
      type(context_box), intent(in) :: cx
      logical, intent(out) :: have_func
      type(trc_xc_functional_t), intent(out) :: func
      character(len=:), allocatable, intent(out) :: refusal
      integer(c_int), intent(out) :: refusal_status
      !! TRC_ERR_UNSUPPORTED for everything below, EXCEPT the missing
      !! --df fit: that one is TRC_ERR_STATE (a setting the caller can
      !! still fix, by calling trc_run_scf again), matching
      !! trc_run_gradient's own status before this check was factored
      !! out -- test/check_grad_df_refusals.F90 checks that distinction.
      refusal = ""
      refusal_status = TRC_ERR_UNSUPPORTED
      have_func = len_trim(cx%opts%functional) > 0
      if (have_func) then
         ! `cx%opts%functional` only ever holds a name trc_set_method/
         ! scf_entry already validated (xc_functional_by_name succeeded,
         ! func%missing empty), so this lookup cannot fail -- it is
         ! re-done here only to read the fields below, not to re-validate.
         if (.not. xc_functional_by_name(cx%opts%functional, func)) then
            refusal = "functional '"//trim(cx%opts%functional)//"' failed to re-resolve -- this should be unreachable"
            return
         end if
      end if
      if (cx%rimp2_ok) then
         refusal = "RI-MP2/double-hybrid correlation gradient not implemented yet"
         return
      end if
      if (have_func) then
         if (func%family == XC_FAMILY_MGGA) then
            ! The energy runs; its tau term has no nuclear gradient in
            ! trc_xc_grad yet, which only does LDA and GGA.
            refusal = "meta-GGA ('"//trim(cx%opts%functional)//"') exchange-correlation gradient not implemented yet"
            return
         end if
      end if
      if (cx%opts%df) then
         ! The fitted gradient differentiates the SAME fit the SCF converged
         ! with (`cx%dfjk`'s B and metric Cholesky factor), so it needs that
         ! fit still in place, and inherits the DF-SCF's own scope:
         ! restricted, Cartesian, no range separation. trc_scf_settings
         ! and trc_scf_df already refuse the SCF outside it; checked again
         ! here so the message names the gradient -- and BEFORE the fit
         ! check, because a result taken from a file (trc_load_scf) reaches
         ! this point without ever having met the DF-SCF's own refusals.
         if (cx%res%nspin /= 1) then
            refusal = "--df gradient is restricted (closed-shell) only, as --df is"
            return
         end if
         if (cx%bb%b%spherical) then
            refusal = "--df gradient needs a Cartesian basis -- the DF derivative kernels are Cartesian only"
            return
         end if
         if (have_func) then
            if (func%exx_lr /= 0.0_dp) then
               refusal = "--df with a range-separated functional is refused (the DF-SCF fits no erf exchange)"
               return
            end if
         end if
         if (.not. cx%df_fit_ok) then
            refusal = "--df gradient needs the density fit the SCF ran with"
            refusal_status = TRC_ERR_STATE
            return
         end if
      end if
      if (len_trim(cx%opts%pcm_model) > 0) then
         refusal = "PCM gradient not implemented yet (no cavity/charge derivative)"
         return
      end if
   end subroutine gradient_capability_check

   !
   ! Move one atom's basis centre to `r` in a basis already built --
   ! EVERY shell whose sh_at (trc_api.F90's own exact-match shell->atom
   ! map, set by basis_build/basis_from_libcint regardless of which
   ! trc_set_basis_*/trc_set_aux_*/trc_set_jkfit_* path built it) equals
   ! `atom_idx`, plus the basis's own at_r(:,atom_idx) -- rather than
   ! rebuilding the basis from its JSON/libcint/array source. This is the
   ! cheap, source-independent route trc_run_hessian's displacement loop
   ! needs 6N times per Hessian: exponents and contraction coefficients
   ! never depend on an atom's position, only the centre each shell is
   ! built on does, so moving it is exact, not an approximation of a
   ! rebuild. A shell with sh_at == 0 (a ghost/probe function, sapcmp's
   ! far-field check -- trc_api.F90's own comment) is left alone: it is
   ! not ON any molecule atom, so no displacement of one moves it.
   !
   ! `!$acc update device` when this basis is device-resident, which on a
   ! GPU build it always is (take_basis calls to_device): the kernels must
   ! not read a device copy that still holds the OLD centre under -acc's
   ! separate memory.
   !
   subroutine move_basis_centre(b, atom_idx, r)
      type(trc_basis_t), intent(inout) :: b
      integer, intent(in) :: atom_idx
      real(dp), intent(in) :: r(3)
      integer :: i
      if (atom_idx < 1 .or. atom_idx > b%natm) return
      b%at_r(:, atom_idx) = r
      do i = 1, b%nshell
         if (b%sh_at(i) == atom_idx) b%sh_r(:, i) = r
      end do
      if (b%on_device) then
         !$acc update device(b%sh_r, b%at_r)
      end if
   end subroutine move_basis_centre

   !
   ! Move every atom of this context to `at_r` (3, natm) -- the orbital
   ! basis, and the auxiliary/jkfit bases if present, all follow -- then
   ! invalidate exactly as trc_set_molecule does for a geometry change:
   ! the last SCF/gradient/polarizability/Hessian result, the DFT grid
   ! cache (a NEW geometry's grid, not the old one, is what the next
   ! trc_run_scf must build), and the --df fit (its metric and
   ! three-centre tensor are geometry-dependent, same as
   ! trc_set_molecule's own comment on `df_fit_ok`). Used only by
   ! trc_run_hessian's displacement loop today; `cx%natm`/`cx%have_mol`
   ! are unchanged (same molecule, same atom count and order), so unlike
   ! trc_set_molecule this never touches `cx%at_z` or the PCM radii
   ! override.
   !
   subroutine hessian_set_geometry(cx, at_r)
      type(context_box), intent(inout) :: cx
      real(dp), intent(in) :: at_r(:, :)
      integer :: a
      cx%at_r = at_r
      do a = 1, cx%natm
         call move_basis_centre(cx%bb%b, a, at_r(:, a))
         if (cx%have_aux) call move_basis_centre(cx%ab%b, a, at_r(:, a))
         if (cx%have_jkfit) call move_basis_centre(cx%jb%b, a, at_r(:, a))
      end do
      call invalidate_runs(cx)
      if (cx%df_fit_ok) call cx%dfjk%release()
      cx%df_fit_ok = .false.
   end subroutine hessian_set_geometry

   function trc_set_geometry(handle, xyz) result(status) bind(c, name="trc_set_geometry")
      !! Move every atom to a NEW geometry, same molecule and atom order
      !! `trc_set_molecule` built (same count, same nuclear charges, same
      !! charge/multiplicity) -- the one call a geometry optimiser's step
      !! loop needs between gradients. Built on `hessian_set_geometry`
      !! (trc_run_hessian's own displacement primitive, PRIVATE until now):
      !! moves the orbital basis's centres in place, and the auxiliary/
      !! jkfit bases' if present, then invalidates exactly as that routine
      !! does -- the last SCF/gradient/polarizability/Hessian result, the
      !! DFT grid cache and the --df fit, all geometry-dependent.
      !!
      !! `trc_set_molecule` is NOT a substitute for this: it replaces
      !! `cx%at_z`/`cx%at_r` and invalidates the same things, but it never
      !! touches a basis already built (`cx%bb`/`cx%ab`/`cx%jb`) -- the
      !! shells stay at their OLD centres, which is correct for setting up
      !! a NEW molecule before any basis exists, and wrong for moving the
      !! atoms of one that already has one. A caller wanting the latter
      !! (an optimiser, or trc_run_hessian's own displacements) uses this.
      !!
      !! `xyz` is atom-major Bohr (x1 y1 z1 x2 ...), the same convention
      !! `trc_set_molecule` takes and `trc_gradient` returns, length
      !! 3*natm. A guess density/orbitals carried across the move (SAD
      !! stays SAD; TRC_GUESS_GIVEN's density and TRC_GUESS_PROJECTED's
      !! "this context's last converged orbitals" are both still valid at
      !! the new geometry -- that overlap-reprojection is exactly what
      !! TRC_GUESS_PROJECTED is for) is a caller choice, not this call's.
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: xyz(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      real(dp), allocatable :: new_r(:, :)
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_mol) then
         call refuse(cx, "trc_set_geometry: set the molecule first")
         return
      end if
      if (.not. cx%have_basis) then
         call refuse(cx, "trc_set_geometry: set the basis first -- this moves an existing "// &
                     "basis's centres, it does not build one")
         return
      end if
      allocate (new_r(3, cx%natm))
      new_r = reshape(real(xyz(1:3*cx%natm), dp), [3, cx%natm])
      call hessian_set_geometry(cx, new_r)
      deallocate (new_r)
      status = TRC_OK
   end function trc_set_geometry

   function trc_run_gradient(handle) result(status) bind(c, name="trc_run_gradient")
      !! The full analytic RHF/UHF nuclear gradient at the converged SCF
      !! this context already holds (trc_scf_grad%trc_scf_gradient: the
      !! one-electron terms, eri%gradient's two-electron part, and nuclear
      !! repulsion -- see that module's own header for the assembly order
      !! and the factor of 1/2). When the orbital basis is spherical,
      !! `trc_scf_gradient` itself folds D and W to the Cartesian space the
      !! derivative kernels need (T D T^T) before touching an integral;
      !! nothing here has to know which case it is in.
      !!
      !! Reuses `cx%la` (initialised once per context by `trc_run_scf`,
      !! released only at `trc_destroy`) rather than a `trc_linalg_t` of
      !! its own: a second handle here would need its own release before
      !! going out of scope (trc_linalg_t has no finalizer), for no reason
      !! -- the gradient's GEMMs are the same size `cx%la` was already
      !! sized for.
      !!
      !! NOT refused any more: a functional with VV10 nonlocal correlation
      !! (wB97X-V). `trc_scf_gradient` adds trc_vv10_grad%vv10_gradient_add's
      !! term on top of the semilocal one, on VV10's own grid
      !! (`cx%opts%nlc_grid_level`, the SAME level the SCF itself converged
      !! it on). wB97M-V is in terco's table too now (a meta-GGA with VV10
      !! and its own range-separated exchange), but the meta-GGA refusal
      !! below catches it before VV10 ever comes up: its energy runs, its
      !! gradient does not, same tau reason as every other meta-GGA.
      !!
      !! REFUSED, clearly, rather than attempted:
      !!   - meta-GGA: the energy runs, but trc_xc_grad has no tau term, so
      !!     the gradient is refused here with a message.
      !!   - RI-MP2 (or, once it exists, a double hybrid): no correlation
      !!     gradient exists; `trc_scf_gradient` is the SCF reference term
      !!     only.
      !!   - --df outside the DF-SCF's own scope (unrestricted, spherical,
      !!     range-separated -- the SCF refuses those first), or with the fit
      !!     gone. A --df gradient otherwise runs the FITTED J/K gradient
      !!     (trc_scf_grad_df.F90) over `cx%dfjk`, the fit the SCF used.
      !!   - PCM: no cavity/surface-charge derivative exists.
      !!
      !! A multi-rank context (`cx%fcomm /= -1`) is NOT refused: the
      !! gradient is not rank-split yet, so rank 0 computes it alone and
      !! broadcasts it, with a warning (see the comment at the call). Only
      !! rank 0's device rebuilds the full pair list, Schwarz screen and ERI
      !! structures for it, which is the memory a split SCF avoided; a real
      !! split is the follow-up. `eri%gradient` itself already supports one
      !! (`this%rank`/`this%nranks`, `sum_ranks`'s allreduce) once
      !! `trc_scf_gradient` passes its `eri%build` a `comm`.
      !!
      !! A plain LDA/GGA functional, a GLOBAL hybrid of one (b3lyp and
      !! friends), or a RANGE-SEPARATED hybrid (wb97x and friends: nonzero
      !! `func%omega`/`exx_lr`, `nlc_b`/`nlc_c` still zero) is NOT refused:
      !! `func` is looked up again here (same `xc_functional_by_name`
      !! trc_set_method already validated the name with) and handed to
      !! `trc_scf_gradient`, along with the SAME grid settings
      !! (`cx%opts%grid_level`, `grid_max_pts`, `rho_tol`) the SCF itself
      !! converged with -- the gradient's quadrature has to be the SCF's,
      !! or the two are not derivatives of the same energy.
      !! `trc_scf_gradient` itself builds the erf(omega r)/r long-range
      !! exchange term (a second eri%gradient pass, its own header) and
      !! trc_xc_grad%xc_gradient_add computes the functional's semilocal
      !! GGA part exactly as for any other GGA -- see that module's own
      !! header on why `omega`/`exx_lr` do not affect it.
      !!
      !! GRID REUSE. `cx%kgrid`, the SCF's own grid(s) kept alive since
      !! `trc_run_scf` (`context_box%kgrid`'s own header), is handed
      !! straight to `trc_scf_gradient` when it is built AND matches
      !! every setting the gradient's quadrature has to match the SCF's
      !! on (`reuse_grid`/`reuse_nlc`, below) -- the same list this
      !! routine's own header already promises (grid level, max points
      !! per batch, rho_tol, nlc grid level); geometry and basis are
      !! covered for free, because ANY change to either would have
      !! cleared `cx%scf_ok` (`invalidate_runs`) before this point could
      !! ever be reached. Otherwise `trc_scf_gradient` rebuilds exactly as
      !! it always has -- this is a pure performance path, never a
      !! correctness one: `grid`/`xg` are read-only there, so a mismatch
      !! caught here is the only thing standing between "reuse" and
      !! "quietly wrong".
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_xc_functional_t) :: func
      logical :: reuse_grid, reuse_nlc
      logical :: have_func
      character(len=:), allocatable :: refusal
      type(comm_t) :: comm
      logical :: collective, root
      real(dp), allocatable :: gbuf(:)
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%grad_ok = .false.
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) then
         cx%message = "trc_run_gradient: run the SCF first"
         return
      end if
      call gradient_capability_check(cx, have_func, func, refusal, status)
      if (len(refusal) > 0) then
         cx%message = "trc_run_gradient: "//refusal
         return
      end if
      ! GRID REUSE -- see this function's own header.
      reuse_grid = have_func .and. cx%kgrid%grid_built .and. &
                   cx%kgrid_level == cx%opts%grid_level .and. cx%kgrid_max_pts == cx%opts%grid_max_pts .and. &
                   cx%kgrid_rho_tol == cx%opts%rho_tol
      reuse_nlc = reuse_grid .and. cx%kgrid%nlc_built .and. (func%nlc_b /= 0.0_dp .or. func%nlc_c /= 0.0_dp) .and. &
                  cx%kgrid_nlc_level == cx%opts%nlc_grid_level
      ! More than one rank: the SCF ran split, but the gradient is not split
      ! yet, so rank 0 computes it alone on its own device and broadcasts
      ! it, and every rank's trc_gradient returns the same numbers. That is
      ! safe because every rank holds the whole converged SCF: the
      ! orbitals, the full XC grids (only the XC evaluation is split, by
      ! batch) and, for --df, the full fitted tensor (cx%dfjk%setup takes
      ! no communicator). It costs one GPU's time and memory for the
      ! gradient, hence the warning.
      call trc_build_comm(cx%fcomm, comm, collective, status)
      if (status /= TRC_OK) return
      root = .true.
      if (collective) then
         root = comm%rank() == 0
         if (root .and. comm%size() > 1) then
            call global_logger%warning("trc_run_gradient: the gradient is single GPU for now; rank 0 "// &
                                       "computes it and broadcasts it to the other ranks")
         end if
      end if
      if (allocated(cx%grad)) deallocate (cx%grad)
      allocate (cx%grad(3, cx%natm))
      if (root) then
         if (cx%opts%df) then
            ! The pair list at the threshold the fit's own was built with
            ! (trc_run_scf's `pl%build(cx%bb%b, cx%opts%eri_thresh)`), so the
            ! derivative runs over exactly the (mu nu|P) the energy used.
            if (have_func) then
               if (reuse_nlc) then
                  call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                        cx%bb%cmo, cx%bb%eps, cx%opts%eri_thresh, cx%grad, &
                                        func=func, grid_level=cx%opts%grid_level, &
                                        grid_max_pts=cx%opts%grid_max_pts, rho_tol=cx%opts%rho_tol, &
                                        nlc_grid_level=cx%opts%nlc_grid_level, dfjk=cx%dfjk, aux=cx%jb%b, &
                                        verbose=cx%opts%verbose, grid=cx%kgrid%grid, xg=cx%kgrid%xg, &
                                        grid_nlc=cx%kgrid%grid_nlc, xg_nlc=cx%kgrid%xg_nlc, pair_eri_cache=cx%peri)
               else if (reuse_grid) then
                  call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                        cx%bb%cmo, cx%bb%eps, cx%opts%eri_thresh, cx%grad, &
                                        func=func, grid_level=cx%opts%grid_level, &
                                        grid_max_pts=cx%opts%grid_max_pts, rho_tol=cx%opts%rho_tol, &
                                        nlc_grid_level=cx%opts%nlc_grid_level, dfjk=cx%dfjk, aux=cx%jb%b, &
                                        verbose=cx%opts%verbose, grid=cx%kgrid%grid, xg=cx%kgrid%xg, &
                                        pair_eri_cache=cx%peri)
               else
                  call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                        cx%bb%cmo, cx%bb%eps, cx%opts%eri_thresh, cx%grad, &
                                        func=func, grid_level=cx%opts%grid_level, &
                                        grid_max_pts=cx%opts%grid_max_pts, rho_tol=cx%opts%rho_tol, &
                                        nlc_grid_level=cx%opts%nlc_grid_level, dfjk=cx%dfjk, aux=cx%jb%b, &
                                        verbose=cx%opts%verbose, pair_eri_cache=cx%peri)
               end if
            else
               call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                     cx%bb%cmo, cx%bb%eps, cx%opts%eri_thresh, cx%grad, &
                                     dfjk=cx%dfjk, aux=cx%jb%b, verbose=cx%opts%verbose, pair_eri_cache=cx%peri)
            end if
         else if (have_func) then
            if (reuse_nlc) then
               call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                     cx%bb%cmo, cx%bb%eps, eri_thresh_env(cx%opts%eri_thresh), cx%grad, &
                                     func=func, grid_level=cx%opts%grid_level, &
                                     grid_max_pts=cx%opts%grid_max_pts, rho_tol=cx%opts%rho_tol, &
                                     nlc_grid_level=cx%opts%nlc_grid_level, verbose=cx%opts%verbose, &
                                     grid=cx%kgrid%grid, xg=cx%kgrid%xg, &
                                     grid_nlc=cx%kgrid%grid_nlc, xg_nlc=cx%kgrid%xg_nlc, pair_eri_cache=cx%peri)
            else if (reuse_grid) then
               call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                     cx%bb%cmo, cx%bb%eps, eri_thresh_env(cx%opts%eri_thresh), cx%grad, &
                                     func=func, grid_level=cx%opts%grid_level, &
                                     grid_max_pts=cx%opts%grid_max_pts, rho_tol=cx%opts%rho_tol, &
                                     nlc_grid_level=cx%opts%nlc_grid_level, verbose=cx%opts%verbose, &
                                     grid=cx%kgrid%grid, xg=cx%kgrid%xg, pair_eri_cache=cx%peri)
            else
               call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                     cx%bb%cmo, cx%bb%eps, eri_thresh_env(cx%opts%eri_thresh), cx%grad, &
                                     func=func, grid_level=cx%opts%grid_level, &
                                     grid_max_pts=cx%opts%grid_max_pts, rho_tol=cx%opts%rho_tol, &
                                     nlc_grid_level=cx%opts%nlc_grid_level, verbose=cx%opts%verbose, &
                                     pair_eri_cache=cx%peri)
            end if
         else
            call trc_scf_gradient(cx%bb%b, cx%la, cx%res%nspin, cx%nalpha, cx%nbeta, &
                                  cx%bb%cmo, cx%bb%eps, eri_thresh_env(cx%opts%eri_thresh), cx%grad, &
                                  verbose=cx%opts%verbose, pair_eri_cache=cx%peri)
         end if
      end if
      if (collective) then
         allocate (gbuf(3*cx%natm))
         if (root) gbuf = reshape(cx%grad, [3*cx%natm])
         call bcast(comm, gbuf, 3*cx%natm, 0)
         cx%grad = reshape(gbuf, [3, cx%natm])
         deallocate (gbuf)
         call comm%finalize()
      end if
      cx%grad_ok = .true.
      status = TRC_OK
   end function trc_run_gradient

   function trc_gradient(handle, grad) result(status) bind(c, name="trc_gradient")
      !! (3, natm), Eh/Bohr, column-major (x, y, z fastest, then atom, the
      !! same atom order `trc_set_molecule` was given) -- after a
      !! successful `trc_run_gradient`. NOT AO-indexed data (see
      !! trc_capi.F90's own "AO NORMALISATION" header): a nuclear gradient
      !! is a derivative of the SCF energy, which does not depend on which
      !! AO normalisation convention computed it, so `trc_set_ao_norm`
      !! plays no role here.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: grad(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%grad_ok) return
      grad(1:3*cx%natm) = real(reshape(cx%grad, [3*cx%natm]), c_double)
      status = TRC_OK
   end function trc_gradient

   function trc_set_hessian_step(handle, step) result(status) bind(c, name="trc_set_hessian_step")
      !! The Cartesian displacement trc_run_hessian's central differences
      !! use, in Bohr; must be positive. Default 0.005: its truncation
      !! error, about 1e-5 Eh/Bohr^2 on water/STO-3G, sits well above the
      !! SCF noise trc_run_hessian's tightened convergence leaves.
      !! Unchanged until this is called. Invalidates the last
      !! trc_run_hessian result (a Hessian at the OLD step is not a
      !! Hessian at the new one), not the SCF/gradient it is built from.
      type(c_ptr), value :: handle
      real(c_double), value :: step
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (step <= 0.0_c_double) then
         call refuse(cx, "trc_set_hessian_step: step must be positive, got "//trim(to_string(real(step, dp))))
         return
      end if
      cx%hess_step = real(step, dp)
      cx%hess_ok = .false.
      status = TRC_OK
   end function trc_set_hessian_step

   function trc_hessian_step(handle, step) result(status) bind(c, name="trc_hessian_step")
      !! The displacement trc_run_hessian would use (or last used) --
      !! trc_set_hessian_step's own default or override.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: step
      integer(c_int) :: status
      type(context_box), pointer :: cx
      step = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      step = real(cx%hess_step, c_double)
      status = TRC_OK
   end function trc_hessian_step

   function trc_run_hessian(handle) result(status) bind(c, name="trc_run_hessian")
      !! A finite-difference Cartesian Hessian from this context's own
      !! analytic gradient, central differences at trc_set_hessian_step's
      !! own h (default 0.005 Bohr) around the converged SCF this context
      !! already holds (trc_run_scf; trc_gradient's reference "run the
      !! SCF first" state requirement, same as trc_run_gradient):
      !!
      !!   H[k,j] = (g_j(x + h e_k) - g_j(x - h e_k)) / (2h)     Eh/Bohr^2
      !!
      !! k = 3*(atom-1)+{x,y,z}, 1-based, row = the displaced coordinate,
      !! then symmetrised, H <- (H + H^T)/2. 6N displaced SCF+gradient
      !! pairs plus one more to restore the reference state (below).
      !!
      !! REFUSED UP FRONT, before any displacement spends SCF time:
      !!   - no SCF yet, or one that did not converge (TRC_ERR_STATE);
      !!   - more than one rank (trc_set_comm): the loop below is already
      !!     written so a later split is "rank r takes k = r, r+P, ..."
      !!     plus a reduction over the gradients, but nothing here runs
      !!     that split yet (TRC_ERR_UNSUPPORTED);
      !!   - anything trc_run_gradient itself would refuse (meta-GGA,
      !!     RI-MP2/double-hybrid, PCM, --df outside its own scope --
      !!     gradient_capability_check, the SAME check trc_run_gradient
      !!     makes, so the two can never disagree about what runs);
      !!   - a non-positive step (trc_set_hessian_step already refuses
      !!     one, so this is the backstop, not the primary gate).
      !!
      !! SCF CONVERGENCE AT THE DISPLACED POINTS. FD noise is of order
      !! (gradient error)/(2h) = 100x the gradient error at h = 0.005, so
      !! the displaced SCFs need tighter convergence than a single point
      !! (water/STO-3G: HESSIAN_MAX_CONV_ENERGY/HESSIAN_MAX_CONV_DIIS
      !! below say how much). This context's OWN conv_energy/conv_diis (trc_set_
      !! convergence) are tightened to at least those floors for every
      !! displaced SCF -- never loosened, so a caller who already asked
      !! for tighter still gets it -- and restored to the caller's own
      !! values before this returns, success or failure.
      !!
      !! WARM START. Every displaced SCF starts from THIS run's own
      !! converged density (TRC_GUESS_GIVEN), not the previous
      !! displacement's -- no guess drift accumulates across 6N points.
      !! The reference density, this context's guess kind/density and
      !! conv_energy/conv_diis are all restored before returning.
      !!
      !! A displaced SCF that fails to converge, or a gradient that
      !! fails, stops the loop (TRC_ERR_NOCONV or whatever status that
      !! call returned) rather than filling the Hessian with whatever
      !! that point produced; the context is still restored to the
      !! reference geometry/density before returning, as on success.
      !!
      !! RESTORING THE REFERENCE STATE. After the loop (success or
      !! failure), the geometry goes back to the one this context held
      !! when trc_run_hessian was called, and the reference SCF is RE-RUN
      !! (warm-started from its own converged density, so it reconverges
      !! in very few iterations) rather than merely remembered -- cheap
      !! next to 6N displaced solves, and it leaves cx%res/cx%grad in
      !! whatever state an ordinary trc_run_scf leaves them (grad_ok
      !! false; call trc_run_gradient again for the reference gradient).
      !! This context's guess kind/density, as set before this call, are
      !! restored afterwards too, for the same reason the convergence
      !! settings are: trc_run_hessian's own internal warm start must not
      !! leak into what a caller sees after it returns.
      !!
      !! trc_hessian reads the result out; cx%hess is deallocated and
      !! hess_ok left false on any failure, so a caller cannot read a
      !! partial Hessian by mistake.
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_xc_functional_t) :: func
      logical :: have_func
      character(len=:), allocatable :: refusal
      real(dp), parameter :: HESSIAN_MAX_CONV_ENERGY = 1.0e-12_dp
      real(dp), parameter :: HESSIAN_MAX_CONV_DIIS = 1.0e-9_dp
      !! Measured on water/STO-3G against pyscf's analytic Hessian: at
      !! the CLI's own defaults (econv 1e-10, gconv/conv_diis 1e-5) the
      !! FD Hessian disagreed with pyscf's analytic one by far more than
      !! the h=0.005 truncation error alone predicts; tightening conv_diis
      !! (the DIIS commutator gate, trc_scf_options_t%conv_diis's own
      !! comment) to 1e-9 and conv_energy to 1e-12 brought the worst
      !! element back down to the truncation-error floor (~1e-5
      !! Eh/Bohr^2 at this h). conv_diis is the tighter requirement here:
      !! it gates the orthogonalised commutator |FDS-SDF|, which is what
      !! actually bounds the gradient's own residual error, not the
      !! energy (quadratic in the density error, and so a much blunter
      !! instrument for a quantity -- the gradient -- linear in it).
      integer :: natm, ndof, k, atom, comp, isign
      real(dp) :: h, sgn
      real(dp), allocatable :: ref_at_r(:, :), disp_r(:, :), d0(:, :, :)
      real(dp), allocatable :: gplus(:), gminus(:)
      real(dp) :: ref_conv_energy, ref_conv_diis, new_conv_energy, new_conv_diis
      integer :: ref_guess_kind
      real(dp), allocatable :: ref_dguess(:, :, :)
      logical :: failed
      character(len=256) :: fail_message
      integer(c_int) :: fail_status, st
      type(timer_type) :: wall_timer, k_timer
      integer :: niter_plus, niter_minus
      character(len=256) :: logline

      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)

      status = TRC_ERR_STATE
      if (.not. cx%scf_ok .or. .not. cx%res%converged) then
         cx%message = "trc_run_hessian: run a converged SCF first"
         return
      end if

      status = TRC_ERR_UNSUPPORTED
      if (cx%fcomm /= -1) then
         cx%message = "trc_run_hessian: more than one rank (trc_set_comm) is not supported yet -- "// &
                      "run this context single-rank; the displacement loop is already indexed so a "// &
                      "later split needs only a reduction over the gradients, not a redesign"
         return
      end if

      call gradient_capability_check(cx, have_func, func, refusal, status)
      if (len(refusal) > 0) then
         cx%message = "trc_run_hessian: "//refusal
         return
      end if

      status = TRC_ERR_BADARG
      h = cx%hess_step
      if (h <= 0.0_dp) then
         cx%message = "trc_run_hessian: the displacement step must be positive (trc_set_hessian_step)"
         return
      end if

      natm = cx%natm
      ndof = 3*natm
      allocate (ref_at_r(3, natm), disp_r(3, natm))
      ref_at_r = cx%at_r
      allocate (d0, source=cx%res%dmat)
      ref_guess_kind = cx%guess_kind
      if (ref_guess_kind == TRC_GUESS_GIVEN .and. allocated(cx%dguess)) allocate (ref_dguess, source=cx%dguess)

      ref_conv_energy = cx%opts%conv_energy
      ref_conv_diis = cx%opts%conv_diis
      new_conv_energy = min(cx%opts%conv_energy, HESSIAN_MAX_CONV_ENERGY)
      new_conv_diis = min(cx%opts%conv_diis, HESSIAN_MAX_CONV_DIIS)
      if (new_conv_energy < ref_conv_energy .or. new_conv_diis < ref_conv_diis) then
         write (logline, "(a,es9.2,a,es9.2,a)") "trc_run_hessian: tightening the displaced SCFs' convergence to "// &
            "dE ", new_conv_energy, ", |FDS-SDF| ", new_conv_diis, &
            " (FD noise is of order 100x the gradient error at h=0.005 Bohr)"
         call global_logger%info(trim(logline))
      end if
      cx%opts%conv_energy = new_conv_energy
      cx%opts%conv_diis = new_conv_diis

      if (allocated(cx%hess)) deallocate (cx%hess)
      allocate (cx%hess(ndof, ndof))
      cx%hess = 0.0_dp
      allocate (gplus(ndof), gminus(ndof))

      failed = .false.
      fail_status = TRC_OK
      fail_message = ""
      call wall_timer%start()

      k_loop: do k = 1, ndof
         atom = (k - 1)/3 + 1
         comp = mod(k - 1, 3) + 1
         call k_timer%start()
         do isign = 1, 2
            sgn = merge(1.0_dp, -1.0_dp, isign == 1)
            disp_r = ref_at_r
            disp_r(comp, atom) = ref_at_r(comp, atom) + sgn*h
            call hessian_set_geometry(cx, disp_r)
            cx%guess_kind = TRC_GUESS_GIVEN
            cx%dguess = d0
            st = trc_run_scf(handle)
            if (st /= TRC_OK) then
               failed = .true.
               fail_status = st
               fail_message = "trc_run_hessian: the SCF at displacement "//trim(to_string(k))//"/"// &
                              trim(to_string(ndof))//merge(" (+h)", " (-h)", isign == 1)// &
                              " did not converge: "//trim(cx%message)
               exit k_loop
            end if
            st = trc_run_gradient(handle)
            if (st /= TRC_OK) then
               failed = .true.
               fail_status = st
               fail_message = "trc_run_hessian: the gradient at displacement "//trim(to_string(k))//"/"// &
                              trim(to_string(ndof))//merge(" (+h)", " (-h)", isign == 1)// &
                              " failed: "//trim(cx%message)
               exit k_loop
            end if
            if (isign == 1) then
               niter_plus = cx%res%iterations
               gplus = reshape(cx%grad, [ndof])
            else
               niter_minus = cx%res%iterations
               gminus = reshape(cx%grad, [ndof])
            end if
         end do
         cx%hess(k, :) = (gplus - gminus)/(2.0_dp*h)
         ! One line per coordinate, not per SCF -- 6N of these at the
         ! INFO level is already as much as a Hessian run should print.
         write (logline, "(a,i0,a,i0,a,i0,a,i0,a,i0,a,f0.3,a)") "trc_run_hessian: displacement ", k, "/", &
            ndof, " (atom ", atom, ", xyz "//merge("x", merge("y", "z", comp == 2), comp == 1)// &
            "): SCF iterations +h/-h ", niter_plus, "/", niter_minus, &
            ", ", k_timer%get_elapsed_time(), " s"
         call global_logger%info(trim(logline))
      end do k_loop

      ! Restore the reference geometry and a consistent converged state
      ! regardless of success or failure -- see this function's own
      ! header.
      call hessian_set_geometry(cx, ref_at_r)
      cx%guess_kind = TRC_GUESS_GIVEN
      cx%dguess = d0
      st = trc_run_scf(handle)
      cx%opts%conv_energy = ref_conv_energy
      cx%opts%conv_diis = ref_conv_diis
      cx%guess_kind = ref_guess_kind
      if (ref_guess_kind == TRC_GUESS_GIVEN .and. allocated(ref_dguess)) then
         cx%dguess = ref_dguess
      else if (allocated(cx%dguess)) then
         deallocate (cx%dguess)
      end if

      if (failed) then
         cx%message = fail_message
         ! The caller must also learn that the context no longer holds a
         ! converged SCF, not only why the loop stopped.
         if (st /= TRC_OK) cx%message = trim(fail_message)//"; restoring the reference SCF then failed too"
         cx%hess_ok = .false.
         if (allocated(cx%hess)) deallocate (cx%hess)
         status = fail_status
         return
      end if
      if (st /= TRC_OK) then
         cx%message = "trc_run_hessian: every displacement succeeded, but restoring the reference SCF "// &
                      "did not reconverge: "//trim(cx%message)
         cx%hess_ok = .false.
         if (allocated(cx%hess)) deallocate (cx%hess)
         status = st
         return
      end if

      cx%hess = 0.5_dp*(cx%hess + transpose(cx%hess))
      cx%hess_ok = .true.
      write (logline, "(a,f0.3,a)") "trc_run_hessian: total time ", wall_timer%get_elapsed_time(), " s"
      call global_logger%info(trim(logline))
      status = TRC_OK
   end function trc_run_hessian

   function trc_hessian(handle, hess) result(status) bind(c, name="trc_hessian")
      !! (3*natm, 3*natm), Eh/Bohr^2, column-major, row/column k =
      !! 3*(atom-1)+{x,y,z} (trc_set_molecule's own atom order) -- after a
      !! successful trc_run_hessian. NOT AO-indexed data, same as
      !! trc_gradient: a Hessian is a second derivative of the SCF
      !! energy, independent of any AO normalisation convention.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: hess(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: ndof
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%hess_ok) return
      ndof = 3*cx%natm
      hess(1:ndof*ndof) = real(reshape(cx%hess, [ndof*ndof]), c_double)
      status = TRC_OK
   end function trc_hessian

   function trc_run_polarizability(handle) result(status) bind(c, name="trc_run_polarizability")
      !! The static dipole polarizability at the converged SCF this
      !! context already holds: CPHF (trc_cphf%trc_scf_polarizability)
      !! for RHF or UHF, Hartree-Fock only; CPKS
      !! (trc_cpks%trc_scf_polarizability_rks) when a functional is set --
      !! restricted only, trc_cpks has no unrestricted entry. Builds its
      !! own pair list and ERI engine at the context's own screening
      !! threshold, and (for CPKS) its own XC grid from `cx%res%grid_r`/
      !! `grid_w` -- the SAME quadrature `trc_run_scf` converged the
      !! density with, since a CPKS response built on a different grid is
      !! not a derivative of the energy that converged. All three are
      !! released before returning; `cx%la` is reused, as `trc_run_gradient`
      !! reuses it.
      !!
      !! REFUSED, as `trc_run_gradient` is and for the same reasons: a
      !! correlated (RI-MP2 or double-hybrid) reference, `--df`, PCM, or a
      !! multi-rank context. An unrestricted functional reference is
      !! refused here too (`trc_cpks` has nothing else to call). A
      !! meta-GGA or double-hybrid functional cannot reach this point at
      !! all -- `trc_set_method` already refused it.
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_xc_functional_t) :: func
      type(trc_xc_grid_t) :: xg
      type(trc_cphf_options_t) :: copts
      type(trc_cphf_result_t) :: cres
      logical :: have_func
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%polar_ok = .false.
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) then
         cx%message = "trc_run_polarizability: run the SCF first"
         return
      end if
      status = TRC_ERR_UNSUPPORTED
      if (cx%bb%b%spherical) then
         ! trc_cphf/trc_cpks know nothing of nao_sph: `cmo`/`eps` there are
         ! explicit-shape dummies at `b%nao`, the CARTESIAN count, while
         ! `cx%bb%cmo`/`cx%bb%eps` (what this entry would otherwise hand
         ! them) are sized at nao_sph for a spherical basis -- a silent
         ! shape mismatch, not merely an unimplemented case. Refused
         ! rather than risked.
         cx%message = "trc_run_polarizability: CPHF/CPKS are Cartesian-only; "// &
                      "set force_cartesian before loading the basis"
         return
      end if
      if (cx%rimp2_ok) then
         cx%message = "trc_run_polarizability: RI-MP2/double-hybrid response not implemented"
         return
      end if
      if (cx%opts%df) then
         cx%message = "trc_run_polarizability: --df response not implemented"
         return
      end if
      if (len_trim(cx%opts%pcm_model) > 0) then
         cx%message = "trc_run_polarizability: PCM response not implemented (no cavity/charge response)"
         return
      end if
      if (cx%fcomm /= -1) then
         cx%message = "trc_run_polarizability: the response path is not rank-split; run it on a single rank"
         return
      end if
      have_func = len_trim(cx%opts%functional) > 0
      if (have_func .and. cx%res%nspin /= 1) then
         cx%message = "trc_run_polarizability: CPKS (a functional is set) is restricted only"
         return
      end if
      ! `trc_cphf`/`trc_cpks` take a pair list and an ERI engine directly,
      ! not `cx%peri` (there is no response-side cache to hand them), so
      ! the build right below would otherwise sit ALONGSIDE the SCF's own
      ! engine, still resident in `cx%peri` since `trc_run_scf` -- doubling
      ! peak ERI device memory for no reason, since nothing reads the old
      ! one once the response starts. Released here instead of reused: a
      ! `trc_run_gradient` call after this one simply rebuilds its own,
      ! same cost as before `cx%peri` existed, which is cheap next to not
      ! fitting in device memory at all.
      call cx%peri%release()
      call pl%build(cx%bb%b, cx%opts%eri_thresh)
      call pl%to_device()
      call eri%build(cx%bb%b, cx%opts%eri_thresh)
      copts%conv_residual = 1.0e-8_dp
      copts%max_iter = 100
      copts%verbose = cx%opts%verbose
      if (have_func) then
         if (.not. xc_functional_by_name(cx%opts%functional, func)) then
            cx%message = "trc_run_polarizability: functional '"//trim(cx%opts%functional)// &
                         "' failed to re-resolve -- this should be unreachable"
            call pl%release()
            call eri%release()
            return
         end if
         call xg%build(size(cx%res%grid_w), cx%res%grid_r, cx%res%grid_w, cx%bb%b, max_pts=cx%opts%grid_max_pts)
         call xg%to_device()
         call trc_scf_polarizability_rks(cx%bb%b, pl, eri, xg, func, cx%nalpha, cx%bb%cmo(:, :, 1), &
                                         cx%bb%eps(:, 1), cx%res%dmat(:, :, 1), cx%alpha, copts, cres, cx%la)
         call xg%release()
      else
         call trc_scf_polarizability(cx%bb%b, pl, eri, cx%nalpha, cx%nbeta, cx%res, cx%alpha, copts, cres, cx%la)
      end if
      call pl%release()
      call eri%release()
      if (.not. cres%converged) then
         cx%message = "trc_run_polarizability: "//trim(cres%message)
         status = TRC_ERR_NOCONV
         return
      end if
      cx%polar_ok = .true.
      status = TRC_OK
   end function trc_run_polarizability

   function trc_polarizability(handle, alpha) result(status) bind(c, name="trc_polarizability")
      !! The static dipole polarizability tensor `trc_run_polarizability`
      !! computed, a.u., column-major (3, 3).
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: alpha(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%polar_ok) return
      alpha(1:9) = real(reshape(cx%alpha, [9]), c_double)
      status = TRC_OK
   end function trc_polarizability

   pure integer function core_orbitals(z) result(res)
      !! Frozen-core orbitals for one element, by period.
      real(dp), intent(in) :: z
      integer :: iz
      iz = nint(z)
      if (iz <= 2) then
         res = 0
      else if (iz <= 10) then
         res = 1
      else if (iz <= 18) then
         res = 5
      else if (iz <= 36) then
         res = 9
      else if (iz <= 54) then
         res = 18
      else
         res = 27
      end if
   end function core_orbitals

   function trc_run_rimp2(handle) result(status) bind(c, name="trc_run_rimp2")
      !! Needs a converged SCF in the context (restricted, or unrestricted --
      !! UHF or UKS, cx%res%nspin == 2) and an auxiliary basis. An
      !! unrestricted reference refuses more than one rank below (see
      !! trc_urimp2_run's own header for why the restricted ring is not a
      !! clean generalisation here). The orbital basis may be spherical:
      !! `cx%res%cmo`/`eps` (at whatever space the SCF ran in) are folded
      !! into the Cartesian, zero-padded space `trc_rimp2_run`/
      !! `trc_urimp2_run` themselves run in by `rimp2_mo_to_cart`
      !! (trc_rimp2_driver.F90's own header says why the pad is exact).
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_pairlist_t) :: pl
      type(trc_rimp2_result_t) :: res
      type(comm_t) :: comm
      integer :: nfrozen, nblk, i, nc
      logical :: collective
      real(dp), allocatable :: cmo_a(:, :), cmo_b(:, :), eps_a(:), eps_b(:)
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%rimp2_ok = .false.
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) then
         cx%message = "trc_run_rimp2: run the SCF first"
         return
      end if
      if (.not. cx%have_aux) then
         cx%message = "trc_run_rimp2: no auxiliary basis was set"
         return
      end if
      ! cfg%nspin is always 1 (nalpha == nbeta, restricted) or 2 (nalpha /=
      ! nbeta, UHF/UKS) -- trc_scf_driver's own rule -- so nothing is
      ! refused here on the reference's shape; ROHF does not exist in this
      ! build, so the spin check alone is enough.
      status = TRC_ERR_UNSUPPORTED
      ! A status, not the DF kernels' own error stop: a host passing an
      ! auxiliary basis (a *-rifit set carries f/g on heavy atoms) above
      ! what this build's density fitting was generated for gets a refusal
      ! it can report, rather than the whole process -- every MPI rank --
      ! stopping. trc_max_angular_momentum gives the ceilings up front.
      if (maxval(cx%bb%b%sh_l) > DF_LMAX) then
         cx%message = "trc_run_rimp2: the orbital basis has l = "// &
                      trim(itoa(int(maxval(cx%bb%b%sh_l))))//", above this build's density-fitting "// &
                      "ceiling "//trim(itoa(DF_LMAX))//" (TERCO_DF_LMAX)"
         return
      end if
      if (maxval(cx%ab%b%sh_l) > DF_LMAX_AUX) then
         cx%message = "trc_run_rimp2: the auxiliary basis has l = "// &
                      trim(itoa(int(maxval(cx%ab%b%sh_l))))//", above this build's density-fitting "// &
                      "auxiliary ceiling "//trim(itoa(DF_LMAX_AUX))//" (TERCO_DF_LMAX + 1)"
         return
      end if
      nfrozen = cx%nfrozen
      if (nfrozen < 0) then
         nfrozen = 0
         do i = 1, cx%natm
            nfrozen = nfrozen + core_orbitals(cx%at_z(i))
         end do
      end if
      status = TRC_ERR_BADARG
      if (nfrozen >= cx%nalpha .or. nfrozen >= cx%nbeta) then
         cx%message = "trc_run_rimp2: every occupied orbital is frozen in at least one spin channel"
         return
      end if
      nblk = cx%aux_block
      if (nblk <= 0) nblk = cx%ab%b%nao
      call trc_build_comm(cx%fcomm, comm, collective, status)
      if (status /= TRC_OK) return
      ! UHF/UKS RI-MP2 (cx%res%nspin == 2) refuses more than one rank --
      ! see trc_urimp2_run's own header for why the restricted ring is not
      ! a clean generalisation to the aa/bb/ab blocks an open-shell run
      ! needs. Refused here, before building the pair list, so nothing
      ! below runs for a request that cannot complete.
      if (cx%res%nspin == 2 .and. collective) then
         status = TRC_ERR_UNSUPPORTED
         cx%message = "trc_run_rimp2: unrestricted RI-MP2 runs on one rank only this round "// &
                      "-- the restricted ring distributes one block of B and lets it meet "// &
                      "itself, which does not generalise cleanly to the aa/bb/ab blocks an "// &
                      "open-shell reference needs (see trc_urimp2_run's own header)"
         call comm%finalize()
         return
      end if
      ! `cx%res%cmo`/`eps` are at the SCF's own space (`scf_space_nao`:
      ! `nao_sph` when `cx%bb%b%spherical`, the Cartesian `nao` otherwise)
      ! -- `trc_rimp2_run`/`trc_urimp2_run` themselves stay Cartesian-
      ! interface, so a spherical orbital basis is folded here, once, the
      ! same fold-and-pad `rimp2_mo_to_cart` is for (see its own header).
      ! A Cartesian basis copies straight through (nc == the SCF's own
      ! space already).
      nc = cx%bb%b%nao
      allocate (cmo_a(nc, nc), eps_a(nc))
      if (cx%res%nspin == 2) allocate (cmo_b(nc, nc), eps_b(nc))
      if (cx%bb%b%spherical) then
         call rimp2_mo_to_cart(cx%bb%b, cx%res%cmo(:, :, 1), cx%res%eps(:, 1), cmo_a, eps_a)
         if (cx%res%nspin == 2) call rimp2_mo_to_cart(cx%bb%b, cx%res%cmo(:, :, 2), cx%res%eps(:, 2), cmo_b, eps_b)
      else
         cmo_a = cx%res%cmo(:, :, 1)
         eps_a = cx%res%eps(:, 1)
         if (cx%res%nspin == 2) then
            cmo_b = cx%res%cmo(:, :, 2)
            eps_b = cx%res%eps(:, 2)
         end if
      end if
      call pl%build(cx%bb%b, cx%opts%eri_thresh)
      call pl%to_device()
      if (cx%res%nspin == 2) then
         ! Unrestricted: no DF chaining (--df is closed-shell only), one
         ! rank only (checked above).
         call trc_urimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cx%nbeta, cmo_a, cmo_b, eps_a, eps_b, res, &
                             nfrozen=nfrozen, aux_block=nblk, verbose=cx%opts%verbose, la=cx%la)
      ! Chained off --df's fit, ONLY when RI-MP2's own auxiliary basis is
      ! exactly the jkfit one (same JSON path), not merely the same shape:
      ! trc_rimp2_run checks dimensions, not that the basis is the same
      ! one, and a path match is this handle's own way to be sure of that.
      else if (cx%df_fit_ok .and. len_trim(cx%jkfit_path) > 0 .and. cx%jkfit_path == cx%aux_path) then
         if (collective) then
            call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cmo_a, eps_a, res, &
                               nfrozen=nfrozen, aux_block=nblk, comm=comm, verbose=cx%opts%verbose, la=cx%la, &
                               dfjk=cx%dfjk)
            call comm%finalize()
         else
            call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cmo_a, eps_a, res, &
                               nfrozen=nfrozen, aux_block=nblk, verbose=cx%opts%verbose, la=cx%la, dfjk=cx%dfjk)
         end if
      else if (collective) then
         call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cmo_a, eps_a, res, &
                            nfrozen=nfrozen, aux_block=nblk, comm=comm, verbose=cx%opts%verbose, la=cx%la)
         call comm%finalize()
      else
         call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cmo_a, eps_a, res, &
                            nfrozen=nfrozen, aux_block=nblk, verbose=cx%opts%verbose, la=cx%la)
      end if
      call pl%release()
      cx%e_os = res%e_os
      cx%e_ss = res%e_ss
      cx%rimp2_ok = .true.
      status = TRC_OK
   end function trc_run_rimp2

   function trc_save_scf(handle, path) result(status) bind(c, name="trc_save_scf")
      !! Write the context's SCF result -- orbitals, orbital energies,
      !! density, energies -- to `path` (NUL-terminated), with the molecule,
      !! basis and method it belongs to; trc_load_scf reads it back. One
      !! rank is enough: every rank holds the same result.
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      character(len=256) :: msg
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) then
         cx%message = "trc_save_scf: run the SCF first"
         return
      end if
      call trc_scf_save(c_path(path), cx%bb%b, cx%at_z, cx%at_r, cx%charge, cx%multiplicity, &
                        cx%nalpha, cx%nbeta, cx%opts%functional, cx%opts%pcm_model, cx%opts%pcm_eps, &
                        cx%res, msg)
      cx%message = msg
      status = merge(TRC_OK, TRC_ERR_BADARG, len_trim(msg) == 0)
   end function trc_save_scf

   function trc_load_scf(handle, path) result(status) bind(c, name="trc_load_scf")
      !! In place of trc_run_scf: take the SCF result from a file
      !! trc_save_scf wrote, so a correlated step (trc_run_rimp2) can follow
      !! without repeating the SCF. The molecule, basis, method and solvent
      !! must be set as for trc_run_scf, and must be the ones the file was
      !! made with; anything else is refused, naming the difference. Every
      !! rank reads the file itself. TRC_ERR_NOCONV if the saved SCF had not
      !! converged, as trc_run_scf would have said.
      type(c_ptr), value :: handle
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      character(len=256) :: msg
      integer :: nelec, nunp, nalpha, nbeta, nspin, n
      integer :: nocc(2)
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%scf_ok = .false.
      cx%rimp2_ok = .false.
      cx%grad_ok = .false.
      cx%polar_ok = .false.
      ! A loaded result has no grid of its own -- any grid `trc_run_scf`
      ! kept on this handle belongs to a DIFFERENT (possibly stale) run,
      ! so `trc_run_gradient` must rebuild (this function's own header,
      ! and trc_capi.F90's REUSE note on `context_box%kgrid`).
      call cx%kgrid%release()
      call cx%peri%release()
      cx%kgrid_level = -1
      cx%kgrid_max_pts = -1
      cx%kgrid_nlc_level = -1
      cx%kgrid_rho_tol = -1.0_dp
      status = TRC_ERR_STATE
      if (.not. (cx%have_mol .and. cx%have_basis)) then
         cx%message = "trc_load_scf: set the molecule and the basis first"
         return
      end if
      ! The same electron count and spin trc_run_scf would use.
      nelec = nint(sum(cx%at_z)) - cx%charge
      nunp = cx%multiplicity - 1
      nalpha = (nelec + nunp)/2
      nbeta = nelec - nalpha
      call scf_spin_rule(nalpha, nbeta, cx%opts%unrestricted, cx%opts%frac_occ, cx%opts%nelec_frac, nspin, nocc)
      status = TRC_ERR_BADARG
      call trc_scf_load(c_path(path), cx%bb%b, cx%at_z, cx%at_r, cx%charge, cx%multiplicity, nspin, &
                        cx%opts%functional, cx%opts%pcm_model, cx%opts%pcm_eps, cx%res, msg)
      cx%message = msg
      if (len_trim(msg) > 0) return
      n = scf_space_nao(cx%bb)
      ! The GEMM/eigensolver handles trc_run_scf would have made: a loaded
      ! result goes straight to trc_run_gradient or trc_run_rimp2, which
      ! reuse `cx%la` and would otherwise call cuBLAS with no handle.
      call cx%la%init(n)
      if (allocated(cx%bb%cmo)) deallocate (cx%bb%cmo, cx%bb%eps)
      allocate (cx%bb%cmo(n, n, nspin), cx%bb%eps(n, nspin))
      cx%bb%cmo = cx%res%cmo
      cx%bb%eps = cx%res%eps
      cx%bb%nalpha = nalpha
      cx%bb%nbeta = nbeta
      cx%nalpha = nalpha
      cx%nbeta = nbeta
      cx%scf_ok = .true.
      status = merge(TRC_OK, TRC_ERR_NOCONV, cx%res%converged)
   end function trc_load_scf

   function c_path(path) result(fpath)
      !! A NUL-terminated C path as a Fortran string.
      character(kind=c_char), intent(in) :: path(*)
      character(len=:), allocatable :: fpath
      integer :: i
      i = 0
      do while (path(i + 1) /= c_null_char)
         i = i + 1
      end do
      allocate (character(len=i) :: fpath)
      do i = 1, len(fpath)
         fpath(i:i) = path(i)
      end do
   end function c_path

   ! ---- read back ---------------------------------------------------------

   function trc_cart_order() result(order_id) bind(c, name="trc_cart_order")
      !! 0 (lexicographic) or 1 (gamess): the within-shell Cartesian order every
      !! AO-indexed array this API exchanges uses, fixed at compile time
      !! (TRC_CART_ORDER) and the same for every handle -- no handle argument,
      !! and it always succeeds.
      integer(c_int) :: order_id
      order_id = int(trc_cart_order_id(), c_int)
   end function trc_cart_order

   function trc_abi_version(major, minor) result(status) bind(c, name="trc_abi_version")
      !! This build's ABI version (TRC_ABI_MAJOR, TRC_ABI_MINOR above), so
      !! a host can compare it against the constants the header it
      !! compiled against carries before trusting any other entry point.
      !! Bump TRC_ABI_MINOR for an additive change (a new entry point, a
      !! new enum value appended at the end); bump TRC_ABI_MAJOR, and
      !! reset TRC_ABI_MINOR to 0, for anything that changes an existing
      !! signature or an existing enum value's meaning. No handle; always
      !! succeeds.
      integer(c_int), intent(out) :: major, minor
      integer(c_int) :: status
      major = TRC_ABI_MAJOR
      minor = TRC_ABI_MINOR
      status = TRC_OK
   end function trc_abi_version

   function trc_select_device(rank, dev) result(status) bind(c, name="trc_select_device")
      !! Bind this process to a GPU by its node-local rank: device
      !! mod(rank, number of devices), what trc_api's trc_bind_device does
      !! for terco's own MPI runs. A host code that runs one single-rank
      !! context per process and never calls trc_set_comm (one fragment per
      !! rank, say) calls this once, before the first trc_create, or every
      !! rank lands on device 0. `dev` is the device chosen, -1 on a build
      !! without the GPU or a node without one. Process-wide; no handle;
      !! always succeeds.
      integer(c_int), value :: rank
      integer(c_int), intent(out) :: dev
      integer(c_int) :: status
      dev = int(trc_bind_device(int(rank)), c_int)
      status = TRC_OK
   end function trc_select_device

   function trc_select_device_env(dev) result(status) bind(c, name="trc_select_device_env")
      !! trc_select_device by the node-local rank the launcher exports
      !! (OMPI_COMM_WORLD_LOCAL_RANK, PALS_LOCAL_RANKID, MPI_LOCALRANKID,
      !! MV2_COMM_WORLD_LOCAL_RANK, then SLURM_LOCALID), with no MPI needed.
      !! Call it BEFORE MPI_Init. A GPU-aware MPI records the current device
      !! there, and a process that switches devices later loses CUDA IPC
      !! between the node's ranks under Cray MPICH. trc_set_comm's own bind
      !! then picks the same device and leaves it alone. `dev` is the device
      !! chosen, -1 when no launcher variable is set (nothing is bound), on
      !! a build without the GPU, or on a node without one. Process-wide; no
      !! handle; always succeeds.
      integer(c_int), intent(out) :: dev
      integer(c_int) :: status
      dev = int(trc_bind_device_env(), c_int)
      status = TRC_OK
   end function trc_select_device_env

   function trc_physical_constant(id, value) result(status) bind(c, name="trc_physical_constant")
      !! One of trc_constants' own CODATA 2018 values, by the TRC_CONST_*
      !! selector above, so a host that keeps its own copy of these (a
      !! GAMESS-input front end converting Angstrom to Bohr, say) can read
      !! terco's instead of carrying a second literal that can drift from
      !! this library's. Process-wide; no handle; TRC_ERR_BADARG and
      !! `value` left at 0 for an `id` this build does not know (an older
      !! caller asking for a selector a newer header added, or a typo).
      !!
      !!   TRC_CONST_BOHR_TO_ANGSTROM    Bohr radius, Angstrom
      !!   TRC_CONST_HARTREE_TO_EV       Hartree energy, eV
      !!   TRC_CONST_HARTREE_TO_KCALMOL  Hartree energy, kcal/mol (the
      !!                                 thermochemical calorie)
      !!   TRC_CONST_AMU_TO_ME           atomic mass unit, electron masses
      !!   TRC_CONST_KB_SI               Boltzmann constant, J/K
      !!   TRC_CONST_NA_SI               Avogadro constant, mol^-1
      integer(c_int), value :: id
      real(c_double), intent(out) :: value
      integer(c_int) :: status
      value = 0.0_c_double
      status = TRC_ERR_BADARG
      select case (id)
      case (TRC_CONST_BOHR_TO_ANGSTROM)
         value = real(BOHR_TO_ANGSTROM, c_double)
      case (TRC_CONST_HARTREE_TO_EV)
         value = real(HARTREE_TO_EV, c_double)
      case (TRC_CONST_HARTREE_TO_KCALMOL)
         value = real(HARTREE_TO_KCALMOL, c_double)
      case (TRC_CONST_AMU_TO_ME)
         value = real(AMU_TO_ME, c_double)
      case (TRC_CONST_KB_SI)
         value = real(KB_SI, c_double)
      case (TRC_CONST_NA_SI)
         value = real(NA_SI, c_double)
      case default
         return
      end select
      status = TRC_OK
   end function trc_physical_constant

   function trc_functional_info(name, flags, exx, c_os, c_ss) result(status) bind(c, name="trc_functional_info")
      !! Look up a functional by name (case-insensitive, trc_xc_functional's
      !! own table) without a context: `flags` is a bitmask --
      !!
      !!   1  known      -- runnable today (trc_set_method would accept it)
      !!   2  hybrid     -- exx or exx_lr nonzero (some exact exchange)
      !!   4  range-separated -- omega nonzero
      !!   8  VV10       -- nlc_b nonzero (nonlocal correlation)
      !!   16 meta-GGA   -- depends on the kinetic-energy density tau; runs
      !!      (bit 1) unless it also needs the density Laplacian (r2SCAN-L),
      !!      which terco does not evaluate (bit 64)
      !!   32 double hybrid -- xc_is_double_hybrid (pt2_cos or pt2_css nonzero)
      !!   64 listed but not runnable -- trc_xc_functional_t%missing is set;
      !!      bit 1 is then NOT set, and `exx`/`c_os`/`c_ss`/the other bits
      !!      reflect only whatever the table's own entry set before
      !!      giving up (most "missing" entries set nothing at all, so
      !!      every other bit and output reads zero -- that is the table's
      !!      gap, not this entry's)
      !!
      !! `exx` is the fraction of ordinary (full-range) exact exchange;
      !! `c_os`/`c_ss` are the opposite-/same-spin RI-MP2 weights a double
      !! hybrid adds on top (both 0 for anything that is not one). An
      !! unknown name returns a non-OK status with flags = 0 and every
      !! output zeroed, rather than TRC_OK with nothing found.
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), intent(out) :: flags
      real(c_double), intent(out) :: exx, c_os, c_ss
      integer(c_int) :: status
      type(trc_xc_functional_t) :: f
      character(len=24) :: fname
      integer :: i
      flags = 0
      exx = 0.0_c_double
      c_os = 0.0_c_double
      c_ss = 0.0_c_double
      fname = ""
      do i = 1, len(fname)
         if (name(i) == c_null_char) exit
         fname(i:i) = name(i)
      end do
      status = TRC_ERR_UNSUPPORTED
      if (.not. xc_functional_by_name(trim(fname), f)) return
      status = TRC_OK
      if (len_trim(f%missing) > 0) then
         flags = ibset(flags, 6)                                    ! 64: listed, not runnable
         if (index(f%missing, "Laplacian") > 0) flags = ibset(flags, 4)  ! 16: meta-GGA
      else
         flags = ibset(flags, 0)                                    ! 1: known/runnable
         if (f%family == XC_FAMILY_MGGA) flags = ibset(flags, 4)    ! 16: meta-GGA
      end if
      if (f%exx /= 0.0_dp .or. f%exx_lr /= 0.0_dp) flags = ibset(flags, 1)  ! 2: hybrid
      if (f%omega /= 0.0_dp) flags = ibset(flags, 2)                       ! 4: range-separated
      if (f%nlc_b /= 0.0_dp) flags = ibset(flags, 3)                       ! 8: VV10
      if (xc_is_double_hybrid(f)) flags = ibset(flags, 5)                  ! 32: double hybrid
      exx = real(f%exx, c_double)
      c_os = real(f%pt2_cos, c_double)
      c_ss = real(f%pt2_css, c_double)
   end function trc_functional_info

   function trc_functional_list(buf, buflen, reqlen) result(status) bind(c, name="trc_functional_list")
      !! XC_AVAILABLE_FUNCTIONALS -- every name trc_set_method/
      !! trc_functional_info can look up that terco can actually run today
      !! -- space-separated and NUL-terminated into `buf` (capacity
      !! `buflen` bytes, NUL included). `reqlen` is ALWAYS set to the
      !! number of bytes that call needed, NUL included, whether or not it
      !! fit: a caller queries it with `buflen = 0` (`buf` is never
      !! touched then) and allocates exactly `reqlen` bytes for a second,
      !! filling call. Returns TRC_ERR_BADARG, leaving `buf` untouched, if
      !! `buflen` is positive but smaller than `reqlen`.
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int), value :: buflen
      integer(c_int), intent(out) :: reqlen
      integer(c_int) :: status
      integer :: i, k
      k = 0
      do i = 1, len_trim(XC_AVAILABLE_FUNCTIONALS)
         if (XC_AVAILABLE_FUNCTIONALS(i:i) == ",") cycle
         k = k + 1
      end do
      reqlen = int(k + 1, c_int)
      status = TRC_OK
      if (buflen <= 0) return
      if (buflen < reqlen) then
         status = TRC_ERR_BADARG
         return
      end if
      k = 0
      do i = 1, len_trim(XC_AVAILABLE_FUNCTIONALS)
         if (XC_AVAILABLE_FUNCTIONALS(i:i) == ",") cycle
         k = k + 1
         buf(k) = XC_AVAILABLE_FUNCTIONALS(i:i)
      end do
      buf(k + 1) = c_null_char
   end function trc_functional_list

   function trc_nao(handle, nao) result(status) bind(c, name="trc_nao")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: nao
      integer(c_int) :: status
      type(context_box), pointer :: cx
      nao = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_basis) return
      nao = int(scf_space_nao(cx%bb), c_int)
      status = TRC_OK
   end function trc_nao

   function trc_nmo(handle, nmo) result(status) bind(c, name="trc_nmo")
      !! Molecular-orbital count: `trc_nao`'s own AO count minus whatever
      !! canonical orthogonalisation dropped as linearly dependent at the
      !! last `trc_run_scf` (`trc_scf_result_t%n_lindep` -- 0 for every
      !! compact basis, nonzero only for a diffuse one close to linear
      !! dependency). Before the first `trc_run_scf` that count is
      !! unknown, so `nmo` is `trc_nao`'s own number.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: nmo
      integer(c_int) :: status
      type(context_box), pointer :: cx
      nmo = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_basis) return
      if (cx%scf_ok) then
         nmo = int(scf_space_nao(cx%bb) - cx%res%n_lindep, c_int)
      else
         nmo = int(scf_space_nao(cx%bb), c_int)
      end if
      status = TRC_OK
   end function trc_nmo

   function trc_mo_occupations(handle, occ) result(status) bind(c, name="trc_mo_occupations")
      !! Occupation numbers, (nmo, nspin) column-major (`trc_nmo`'s own
      !! size), for the last converged `trc_run_scf`: 2.0 per orbital up
      !! to `nalpha` and restricted (nspin 1), or 1.0 up to `nalpha`/
      !! `nbeta` on each spin's own column (nspin 2); 0 above. REFUSED
      !! (TRC_ERR_UNSUPPORTED) if this context's guess used fractional
      !! occupation: that per-orbital Aufbau split is private to
      !! `trc_scf_run` (its own `frac_occupations`) and not re-derivable
      !! from the stored result alone -- unreachable through this C API
      !! today, since nothing here sets fractional occupation, but
      !! refused rather than silently wrong for a Fortran-driver-built
      !! context that did.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: occ(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: n, nspin, k, i, s
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) then
         cx%message = "trc_mo_occupations: run the SCF first"
         return
      end if
      status = TRC_ERR_UNSUPPORTED
      if (cx%opts%frac_occ) then
         cx%message = "trc_mo_occupations: this context's guess used fractional occupation; "// &
                      "its per-orbital Aufbau split is not exposed here"
         return
      end if
      n = int(scf_space_nao(cx%bb)) - cx%res%n_lindep
      nspin = cx%res%nspin
      k = 0
      do s = 1, nspin
         do i = 1, n
            k = k + 1
            if (nspin == 1) then
               occ(k) = merge(2.0_c_double, 0.0_c_double, i <= cx%nalpha)
            else if (s == 1) then
               occ(k) = merge(1.0_c_double, 0.0_c_double, i <= cx%nalpha)
            else
               occ(k) = merge(1.0_c_double, 0.0_c_double, i <= cx%nbeta)
            end if
         end do
      end do
      status = TRC_OK
   end function trc_mo_occupations

   function trc_is_spherical(handle, flag) result(status) bind(c, name="trc_is_spherical")
      !! Non-zero if the orbital basis this context holds is running in
      !! the spherical AO space (nao_sph of trc_nao, rather than the
      !! Cartesian nao) -- the basis's own JSON function_type unless
      !! force_cartesian overrode it. For a caller checking which space
      !! `trc_nao`'s number, and every AO-indexed array it sizes, is in.
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: flag
      integer(c_int) :: status
      type(context_box), pointer :: cx
      flag = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_basis) return
      flag = merge(1_c_int, 0_c_int, cx%bb%b%spherical)
      status = TRC_OK
   end function trc_is_spherical

   function trc_nspin(handle, nspin) result(status) bind(c, name="trc_nspin")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: nspin
      integer(c_int) :: status
      type(context_box), pointer :: cx
      nspin = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      nspin = int(cx%res%nspin, c_int)
      status = TRC_OK
   end function trc_nspin

   function trc_energy(handle, energy) result(status) bind(c, name="trc_energy")
      !! The SCF total, plus the RI-MP2 correlation if that ran: the FULL
      !! E_os + E_ss for a plain post-SCF RI-MP2 (any functional, or
      !! Hartree-Fock), or c_os E_os + c_ss E_ss for a double hybrid --
      !! trc_rimp2_energy below always reports the raw, unscaled E_os/E_ss
      !! either way, so a caller after the PT2 part alone never sees this
      !! weighting. A caller asking for a double hybrid's total needs
      !! nothing beyond the ordinary trc_run_scf/trc_run_rimp2/trc_energy
      !! sequence any RI-MP2 run already uses -- the weighting is applied
      !! here, once, from the functional's own pt2_cos/pt2_css.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: energy
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_xc_functional_t) :: func
      energy = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      energy = real(cx%res%energy, c_double)
      if (cx%rimp2_ok) then
         if (len_trim(cx%opts%functional) > 0 .and. xc_functional_by_name(cx%opts%functional, func)) then
            if (xc_is_double_hybrid(func)) then
               energy = energy + real(func%pt2_cos*cx%e_os + func%pt2_css*cx%e_ss, c_double)
               status = TRC_OK
               return
            end if
         end if
         energy = energy + real(cx%e_os + cx%e_ss, c_double)
      end if
      status = TRC_OK
   end function trc_energy

   function trc_energy_parts(handle, e_nuc, e_one, e_two, e_xc) result(status) bind(c, name="trc_energy_parts")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: e_nuc, e_one, e_two, e_xc
      integer(c_int) :: status
      type(context_box), pointer :: cx
      e_nuc = 0.0_c_double
      e_one = 0.0_c_double
      e_two = 0.0_c_double
      e_xc = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      e_nuc = real(cx%res%e_nuc, c_double)
      e_one = real(cx%res%e_one, c_double)
      e_two = real(cx%res%e_two, c_double)
      e_xc = real(cx%res%e_xc, c_double)
      status = TRC_OK
   end function trc_energy_parts

   function trc_converged(handle, flag) result(status) bind(c, name="trc_converged")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: flag
      integer(c_int) :: status
      type(context_box), pointer :: cx
      flag = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      flag = merge(1_c_int, 0_c_int, cx%res%converged)
      status = TRC_OK
   end function trc_converged

   function trc_iterations(handle, niter) result(status) bind(c, name="trc_iterations")
      type(c_ptr), value :: handle
      integer(c_int), intent(out) :: niter
      integer(c_int) :: status
      type(context_box), pointer :: cx
      niter = 0
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      niter = int(cx%res%iterations, c_int)
      status = TRC_OK
   end function trc_iterations

   function trc_density(handle, dmat) result(status) bind(c, name="trc_density")
      !! (nao, nao, nspin), column-major; restricted holds the total.
      !!
      !! CONTRAVARIANT under trc_set_ao_norm: D_out = N^-1 D_terco N^-1 per
      !! spin slice when UNIT. Scaled on a COPY -- `cx%res%dmat` itself
      !! stays in terco's own convention, because trc_run_rimp2 and a
      !! second trc_run_scf both read it as such internally.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: dmat(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      real(dp), allocatable :: d(:, :, :), nrm(:)
      integer :: n, ns, s
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      n = scf_space_nao(cx%bb)
      ns = cx%res%nspin
      d = cx%res%dmat
      if (cx%bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         allocate (nrm(n))
         call ao_norm_vector(cx%bb%b%nshell, cx%bb%b%sh_l, cx%bb%b%sh_ao, n, cx%bb%ao_norm_kind, nrm)
         do s = 1, ns
            call ao_scale_contra(nrm, d(:, :, s))
         end do
         deallocate (nrm)
      end if
      dmat(1:n*n*ns) = real(reshape(d, [n*n*ns]), c_double)
      status = TRC_OK
   end function trc_density

   function trc_hirshfeld_charges(handle, charges) result(status) bind(c, name="trc_hirshfeld_charges")
      !! Hirshfeld atomic charges from the converged density, one per atom
      !! in `trc_set_molecule`'s own order -- `trc_hirshfeld.F90`'s
      !! stockholder partition (Hirshfeld, Theor. Chim. Acta 44, 129
      !! (1977)), with the free atoms H-Kr's fixed wB97X/def2-TZVP recipe.
      !! Works for RHF/UHF/RKS/UKS and with C-PCM/IEF-PCM (the charge
      !! analysis never touches the solvent operator, only the converged
      !! density). TRC_ERR_BADARG for an element outside H-Kr (terco has
      !! no effective-core-potential support to refuse separately);
      !! TRC_ERR_STATE when the free-atom data file cannot be read (set
      !! TRC_HIRSHFELD_DATA).
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: charges(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(error_t) :: err
      real(dp), allocatable :: q(:)
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return

      call state_hirshfeld_charges(cx, q, err)
      if (err%has_error()) then
         cx%message = err%get_message()
         status = hirshfeld_status(err)
         return
      end if
      charges(1:cx%bb%b%natm) = real(q, c_double)
      status = TRC_OK
   end function trc_hirshfeld_charges

   subroutine state_hirshfeld_charges(cx, q, error)
      !! The converged density THIS context currently holds, Hirshfeld-
      !! condensed -- trc_hirshfeld_charges' own conversion (the total
      !! density in the SCF's own space, then into the Cartesian AO space
      !! `trc_hirshfeld_populations` reads, the same transform
      !! `trc_scf_fock`'s own XC call goes through, sph_to_cart's `d_c` --
      !! see trc_sph.F90's own header), factored out here so
      !! trc_fukui_hirshfeld can call it at N, N+1 and N-1 without three
      !! copies of the same transform.
      type(context_box), intent(in) :: cx
      real(dp), allocatable, intent(out) :: q(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d_total_sph(:, :), d_total_cart(:, :)
      integer :: n, ns

      n = scf_space_nao(cx%bb)
      ns = cx%res%nspin
      allocate (d_total_sph(n, n))
      if (ns == 1) then
         d_total_sph = cx%res%dmat(:, :, 1)
      else
         d_total_sph = cx%res%dmat(:, :, 1) + cx%res%dmat(:, :, 2)
      end if
      if (cx%bb%b%spherical) then
         allocate (d_total_cart(cx%bb%b%nao, cx%bb%b%nao))
         call sph_sph_to_cart(cx%bb%b, d_total_sph, d_total_cart)
      else
         d_total_cart = d_total_sph
      end if

      call trc_hirshfeld_populations(cx%bb%b, d_total_cart, q, error)
   end subroutine state_hirshfeld_charges

   function trc_fukui_hirshfeld(handle, f_plus, f_minus, f_zero, q_n, q_np1, q_nm1) &
      result(status) bind(c, name="trc_fukui_hirshfeld")
      !! Hirshfeld-condensed Fukui indices by finite difference in the
      !! electron count, ported from metalquicha's mqc_czt_fukui
      !! (fukui_indices/fukui_from_states, "hirshfeld" scheme only): three
      !! SCFs on THIS handle -- the converged neutral already here, plus
      !! N+1 and N-1, both doublets --
      !!
      !!   f_A+ = q_A(N) - q_A(N+1)      f_A- = q_A(N-1) - q_A(N)
      !!
      !! condensed by trc_hirshfeld_populations, the SAME partition for all
      !! three states (mqc_czt_fukui.f90's own header: "the three
      !! calculations share one czt_molecule_t, so they see the same
      !! geometry and the same basis functions by construction" -- here,
      !! the same context_box, and trc_set_charge keeps the SAME basis and
      !! --df fit across all three, its own header).
      !!
      !! ONLY for a closed-shell neutral (multiplicity 1, at least two
      !! electrons): mqc's own fukui_indices refuses anything else,
      !! verbatim -- a radical's ions are not both doublets, and guessing
      !! what they should be is a caller decision this entry point does
      !! not make.
      !!
      !! Each ion runs UNSEEDED, from the handle's own guess kind (so
      !! `--guess minao` reaches the ions: a SAD guess stalls an NA-040-type
      !! nitrosamine anion that minao converges), except that SAD becomes
      !! SADQ, the charged-atom guess the CLI picks for any ion, and a guess
      !! that is a stored density (TRC_GUESS_GIVEN, TRC_GUESS_PROJECTED)
      !! belongs to the neutral, so the ions fall back to TRC_GUESS_GWH,
      !! mqc's own default when no neutral orbitals are offered. A caller wanting the
      !! neutral's orbitals to seed the ions instead (mqc's "neutral"
      !! mode) drives the three states by hand: trc_set_charge,
      !! trc_set_guess(TRC_GUESS_GIVEN, ...) and trc_run_scf, exactly what
      !! this function does internally but with its own density.
      !!
      !! Leaves the handle's charge and multiplicity as it found them:
      !! restored after the two ion SCFs, and the neutral SCF is RE-RUN
      !! (warm-started from its own converged density) before this
      !! returns, success or failure -- trc_run_hessian's own restore
      !! pattern (same file, below), because cx%res has to match
      !! cx%charge/cx%multiplicity again, not merely the two scalar
      !! fields.
      !!
      !! TRC_ERR_STATE: no converged SCF yet. TRC_ERR_BADARG: multiplicity
      !! is not 1, fewer than two electrons, an element outside
      !! Hirshfeld's H-Kr recipe, or the condensed indices failing their
      !! one-electron sum rule. TRC_ERR_NOCONV: an ion did not converge.
      !! With more than one rank (trc_set_comm) every rank must make this
      !! call, as for trc_run_scf: the three SCFs are collective, and the
      !! Hirshfeld analysis is replicated (every rank holds the same density
      !! and returns the same indices).
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: f_plus(*), f_minus(*), f_zero(*)
      real(c_double), intent(out) :: q_n(*), q_np1(*), q_nm1(*)
      integer(c_int) :: status

      type(context_box), pointer :: cx
      type(error_t) :: err
      real(dp), allocatable :: q_neutral(:), q_anion(:), q_cation(:)
      real(dp), allocatable :: fp(:), fm(:), fz(:)
      real(dp), allocatable :: d0(:, :, :), ref_dguess(:, :, :)
      integer :: ref_charge, ref_multiplicity, ref_guess_kind, nelec, ion_guess_kind
      integer(c_int) :: st, fail_status
      logical :: failed
      character(len=256) :: fail_message

      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)

      status = TRC_ERR_STATE
      if (.not. cx%scf_ok .or. .not. cx%res%converged) then
         cx%message = "trc_fukui_hirshfeld: run a converged SCF first"
         return
      end if

      status = TRC_ERR_BADARG
      if (cx%multiplicity /= 1) then
         cx%message = "trc_fukui_hirshfeld: Fukui indices need a closed-shell neutral "// &
                      "(multiplicity 1); this handle's multiplicity is "// &
                      trim(itoa(cx%multiplicity))//", and its ions would not both be doublets"
         return
      end if
      nelec = nint(sum(cx%at_z)) - cx%charge
      if (nelec < 2) then
         cx%message = "trc_fukui_hirshfeld: fewer than two electrons; no closed shell to "// &
                      "remove one from"
         return
      end if

      call state_hirshfeld_charges(cx, q_neutral, err)
      if (err%has_error()) then
         status = hirshfeld_status(err)
         cx%message = "trc_fukui_hirshfeld: the neutral's charges: "//err%get_message()
         return
      end if

      ! Save the reference state -- trc_run_hessian's own pattern, below.
      ref_charge = cx%charge
      ref_multiplicity = cx%multiplicity
      allocate (d0, source=cx%res%dmat)
      ref_guess_kind = cx%guess_kind
      if (ref_guess_kind == TRC_GUESS_GIVEN .and. allocated(cx%dguess)) then
         allocate (ref_dguess, source=cx%dguess)
      end if
      ! The ions start from the handle's own guess kind, with two changes:
      ! SAD becomes SADQ (the charged-atom guess, the CLI's own default for
      ! any ion), and a stored density (GIVEN, PROJECTED) belongs to the
      ! neutral, so the ions fall back to GWH, mqc's own unseeded default.
      select case (ref_guess_kind)
      case (TRC_GUESS_SAD)
         ion_guess_kind = TRC_GUESS_SADQ
      case (TRC_GUESS_GIVEN, TRC_GUESS_PROJECTED)
         ion_guess_kind = TRC_GUESS_GWH
      case default
         ion_guess_kind = ref_guess_kind
      end select

      failed = .false.
      fail_status = TRC_OK
      fail_message = ""

      ! The anion, N+1 electrons, a doublet, unseeded.
      st = trc_set_charge(handle, int(ref_charge - 1, c_int), 2_c_int)
      if (st == TRC_OK) then
         cx%guess_kind = ion_guess_kind
         st = trc_run_scf(handle)
      end if
      if (st /= TRC_OK .or. .not. cx%res%converged) then
         failed = .true.
         fail_status = merge(TRC_ERR_NOCONV, st, st == TRC_OK)
         fail_message = "trc_fukui_hirshfeld: the anion (N+1): "//trim(cx%message)
      else
         call state_hirshfeld_charges(cx, q_anion, err)
         if (err%has_error()) then
            failed = .true.
            fail_status = hirshfeld_status(err)
            fail_message = "trc_fukui_hirshfeld: the anion's charges: "//err%get_message()
         end if
      end if

      ! The cation, N-1 electrons, a doublet.
      if (.not. failed) then
         st = trc_set_charge(handle, int(ref_charge + 1, c_int), 2_c_int)
         if (st == TRC_OK) then
            cx%guess_kind = ion_guess_kind
            st = trc_run_scf(handle)
         end if
         if (st /= TRC_OK .or. .not. cx%res%converged) then
            failed = .true.
            fail_status = merge(TRC_ERR_NOCONV, st, st == TRC_OK)
            fail_message = "trc_fukui_hirshfeld: the cation (N-1): "//trim(cx%message)
         else
            call state_hirshfeld_charges(cx, q_cation, err)
            if (err%has_error()) then
               failed = .true.
               fail_status = hirshfeld_status(err)
               fail_message = "trc_fukui_hirshfeld: the cation's charges: "//err%get_message()
            end if
         end if
      end if

      if (.not. failed) then
         call trc_fukui_indices(q_neutral, q_anion, q_cation, fp, fm, fz, err)
         if (err%has_error()) then
            failed = .true.
            fail_status = TRC_ERR_BADARG
            fail_message = "trc_fukui_hirshfeld: "//err%get_message()
         end if
      end if

      ! Restore the reference state regardless of success or failure --
      ! trc_run_hessian's own pattern (same file): cx%res has to match
      ! cx%charge/cx%multiplicity again, not merely the two scalar fields.
      st = trc_set_charge(handle, int(ref_charge, c_int), int(ref_multiplicity, c_int))
      cx%guess_kind = TRC_GUESS_GIVEN
      cx%dguess = d0
      st = trc_run_scf(handle)
      cx%guess_kind = ref_guess_kind
      if (ref_guess_kind == TRC_GUESS_GIVEN .and. allocated(ref_dguess)) then
         cx%dguess = ref_dguess
      else if (allocated(cx%dguess)) then
         deallocate (cx%dguess)
      end if

      if (failed) then
         cx%message = fail_message
         if (st /= TRC_OK) then
            cx%message = trim(fail_message)//"; restoring the reference SCF then failed too"
         end if
         status = fail_status
         return
      end if
      if (st /= TRC_OK) then
         cx%message = "trc_fukui_hirshfeld: the ions converged, but restoring the reference "// &
                      "SCF did not reconverge: "//trim(cx%message)
         status = st
         return
      end if

      f_plus(1:size(fp)) = real(fp, c_double)
      f_minus(1:size(fm)) = real(fm, c_double)
      f_zero(1:size(fz)) = real(fz, c_double)
      q_n(1:size(q_neutral)) = real(q_neutral, c_double)
      q_np1(1:size(q_anion)) = real(q_anion, c_double)
      q_nm1(1:size(q_cation)) = real(q_cation, c_double)
      status = TRC_OK
   end function trc_fukui_hirshfeld

   pure integer(c_int) function hirshfeld_status(err) result(status)
      !! The C status for a failed Hirshfeld analysis: the data file is the
      !! environment's to fix (TRC_ERR_STATE), anything else is the
      !! molecule (an element outside H-Kr, TRC_ERR_BADARG).
      type(error_t), intent(in) :: err
      if (err%code == ERROR_IO) then
         status = TRC_ERR_STATE
      else
         status = TRC_ERR_BADARG
      end if
   end function hirshfeld_status

   function trc_mo_coefficients(handle, cmo) result(status) bind(c, name="trc_mo_coefficients")
      !! CONTRAVARIANT under trc_set_ao_norm, on the AO (row) index only:
      !! C_out = N^-1 C_terco per spin slice when UNIT -- see trc_density
      !! for why this scales a copy, not `cx%res%cmo` itself.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: cmo(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      real(dp), allocatable :: c(:, :, :), nrm(:)
      integer :: n, ns, s
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      n = scf_space_nao(cx%bb)
      ns = cx%res%nspin
      c = cx%res%cmo
      if (cx%bb%ao_norm_kind == TRC_AO_NORM_UNIT) then
         allocate (nrm(n))
         call ao_norm_vector(cx%bb%b%nshell, cx%bb%b%sh_l, cx%bb%b%sh_ao, n, cx%bb%ao_norm_kind, nrm)
         do s = 1, ns
            call ao_scale_contra_row(nrm, c(:, :, s))
         end do
         deallocate (nrm)
      end if
      cmo(1:n*n*ns) = real(reshape(c, [n*n*ns]), c_double)
      status = TRC_OK
   end function trc_mo_coefficients

   function trc_mo_energies(handle, eps) result(status) bind(c, name="trc_mo_energies")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: eps(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: n, ns
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      n = scf_space_nao(cx%bb)
      ns = cx%res%nspin
      eps(1:n*ns) = real(reshape(cx%res%eps, [n*ns]), c_double)
      status = TRC_OK
   end function trc_mo_energies

   function trc_rimp2_energy(handle, e_os, e_ss) result(status) bind(c, name="trc_rimp2_energy")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: e_os, e_ss
      integer(c_int) :: status
      type(context_box), pointer :: cx
      e_os = 0.0_c_double
      e_ss = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%rimp2_ok) return
      e_os = real(cx%e_os, c_double)
      e_ss = real(cx%e_ss, c_double)
      status = TRC_OK
   end function trc_rimp2_energy

   function trc_message(handle, buf, buflen) result(status) bind(c, name="trc_message")
      !! The last message, NUL-terminated into `buf` of `buflen` chars.
      type(c_ptr), value :: handle
      character(kind=c_char), intent(out) :: buf(*)
      integer(c_int), value :: buflen
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: i, n
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      n = min(len_trim(cx%message), int(buflen) - 1)
      do i = 1, n
         buf(i) = cx%message(i:i)
      end do
      if (buflen > 0) buf(n + 1) = c_null_char
      status = TRC_OK
   end function trc_message

   ! ---- borrowed handles for the Fock family ------------------------------

   function trc_context_basis(handle, basis) result(status) bind(c, name="trc_context_basis")
      !! The basis handle inside the context: BORROWED, never destroyed by
      !! the caller, valid until the context is destroyed or the basis reset.
      type(c_ptr), value :: handle
      type(c_ptr), intent(out) :: basis
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(basis_box), pointer :: bp
      basis = c_null_ptr
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_basis) return
      ! Through a pointer: a component of the box has no TARGET attribute of
      ! its own, and c_loc of it is what nvfortran warns about.
      bp => cx%bb
      basis = c_loc(bp)
      status = TRC_OK
   end function trc_context_basis

   function trc_context_eri(handle, eri) result(status) bind(c, name="trc_context_eri")
      !! The ERI context at the context's screening threshold, built on first
      !! use. BORROWED, as trc_context_basis.
      type(c_ptr), value :: handle
      type(c_ptr), intent(out) :: eri
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(eri_box), pointer :: ep
      type(comm_t) :: comm
      logical :: collective
      eri = c_null_ptr
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%have_basis) return
      if (.not. allocated(cx%eb)) then
         allocate (cx%eb)
         call trc_build_comm(cx%fcomm, comm, collective, status)
         if (status /= TRC_OK) then
            deallocate (cx%eb)
            return
         end if
         if (collective) then
            call cx%eb%e%build(cx%bb%b, cx%opts%eri_thresh, comm, batch_res=cx%opts%batch_res)
         else
            call cx%eb%e%build(cx%bb%b, cx%opts%eri_thresh, batch_res=cx%opts%batch_res)
         end if
      end if
      ep => cx%eb
      eri = c_loc(ep)
      status = TRC_OK
   end function trc_context_eri

end module trc_capi
