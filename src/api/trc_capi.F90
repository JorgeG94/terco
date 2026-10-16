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
! RI-MP2 (trc_run_rimp2) refuses a spherical orbital basis outright: its
! three-centre tensor is Cartesian and nothing converts the MO
! coefficients it is contracted against, which is the same reasoning as
! trc_compute_df2c/df3c's own exemption below, just enforced as an error
! instead of a silent exemption because RI-MP2 crosses actual MO data.
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
!     trc_iterations, trc_rimp2_energy, trc_nao, trc_nspin, trc_message,
!     trc_set_pcm, trc_pcm_energy
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
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, &
                        trc_1e, trc_df_2c, trc_df_3c, trc_multipoles, &
                        TRC_NMULT
   use trc_eri, only: trc_eri_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_settings, only: scf_spin_rule
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build, trc_sadq_build
   use trc_sap, only: trc_sap_build
   use trc_error, only: error_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_AVAILABLE_FUNCTIONALS
   use trc_tables, only: TRC_ERI_LMAX => LMAX
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use pic_mpi_lib, only: comm_t, comm_from_handle
   use trc_linalg, only: trc_linalg_t
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t
   use trc_scf_file, only: trc_scf_save, trc_scf_load
   use trc_cart_perm, only: trc_cart_order_id
   use trc_ao_norm, only: ao_norm_shell_ => TRC_AO_NORM_SHELL, ao_norm_unit_ => TRC_AO_NORM_UNIT, &
                          ao_norm_vector, ao_scale_cov, ao_scale_contra, ao_scale_contra_row
   use trc_sph, only: sph_cart_to_sph
   use trc_pcm, only: pcm_model_from_name, PCM_MODEL_NONE, PCM_DEFAULT_LEBEDEV, PCM_XI_ORDERS
   use trc_pcm_radii, only: PCM_DEFAULT_RADII_SCALE
   implicit none
   private

   public :: trc_basis_create, trc_basis_destroy, trc_basis_nao
   public :: trc_basis_create_libcint, trc_basis_maxl
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
   public :: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_libcint, trc_set_basis_arrays, &
             trc_set_basis_json, trc_set_aux_libcint, trc_set_aux_arrays, trc_set_aux_json, trc_set_method, trc_set_nlc_grid, &
             trc_set_convergence, trc_set_screening, trc_set_incremental, &
             trc_set_batching, trc_set_guess, trc_set_comm, trc_set_verbose, &
             trc_set_rimp2, trc_run_scf, trc_run_rimp2, trc_nao, trc_nspin, trc_energy, trc_energy_parts, &
             trc_converged, trc_iterations, trc_density, trc_mo_coefficients, trc_mo_energies, &
             trc_rimp2_energy, trc_message, trc_context_basis, trc_context_eri
   public :: trc_save_scf, trc_load_scf
   public :: trc_set_pcm, trc_pcm_energy
   public :: trc_set_ao_norm, trc_get_ao_norm
   public :: trc_set_force_cartesian, trc_is_spherical
   public :: trc_max_angular_momentum
   public :: capi_scf, capi_scf_mpi, capi_rimp2
   public :: trc_cart_order

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
                                        TRC_GUESS_SAP = 5
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
      logical :: have_mol = .false.
      ! bases; `bb` also keeps the last SCF's orbitals, as the old entries do
      type(basis_box) :: bb, ab
      logical :: have_basis = .false., have_aux = .false.
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
      call cx%la%release()
      if (cx%have_basis) call cx%bb%b%release()
      if (cx%have_aux) call cx%ab%b%release()
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

   subroutine invalidate_runs(cx)
      !! A setting the results depended on changed.
      type(context_box), intent(inout) :: cx
      cx%scf_ok = .false.
      cx%rimp2_ok = .false.
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
      allocate (cx%at_z(natm), cx%at_r(3, natm))
      do i = 1, natm
         cx%at_z(i) = real(z(i), dp)
         cx%at_r(:, i) = real(xyz(3*i - 2:3*i), dp)
      end do
      cx%natm = int(natm)
      cx%charge = int(charge)
      cx%multiplicity = int(multiplicity)
      cx%have_mol = .true.
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_molecule

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
   end subroutine take_basis

   function basis_libcint(handle, aux, atm, natm, bas, nshell, env, nenv, cartesian) result(status)
      type(c_ptr), value :: handle
      logical, intent(in) :: aux
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
         call refuse(cx, "trc_set_basis_libcint: natm, nshell and nenv must all be positive")
         return
      end if
      if (cartesian == 0) then
         status = TRC_ERR_UNSUPPORTED
         call refuse(cx, "trc_set_basis_libcint: terco builds every shell Cartesian; pass cartesian = 1 "// &
                     "and a Cartesian basis, a spherical one is not reinterpreted")
         return
      end if
      call b%from_libcint(int(atm(1:6*natm)), int(natm), int(bas(1:8*nshell)), int(nshell), real(env(1:nenv), dp))
      call take_basis(cx, aux, b)
      status = TRC_OK
   end function basis_libcint

   function trc_set_basis_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) result(status) &
      bind(c, name="trc_set_basis_libcint")
      type(c_ptr), value :: handle
      integer(c_int), value :: natm, nshell, nenv, cartesian
      integer(c_int), intent(in) :: atm(*), bas(*)
      real(c_double), intent(in) :: env(*)
      integer(c_int) :: status
      status = basis_libcint(handle, .false., atm, natm, bas, nshell, env, nenv, cartesian)
   end function trc_set_basis_libcint

   function trc_set_aux_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) result(status) &
      bind(c, name="trc_set_aux_libcint")
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
                   cx%natm, cx%at_z, cx%at_r, int(maxnp))
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
      ! An auxiliary/fitting basis stays Cartesian unconditionally --
      ! density fitting and RI-MP2 are untouched by this feature (the
      ! milestone's own scope), so there is no spherical path for `aux` to
      ! follow even if its own JSON says "gto_spherical" (def2's -jkfit
      ! sets often do). The orbital basis (aux = .false.) follows the file
      ! unless `force_cartesian` was set (trc_set_force_cartesian).
      call trc_basis_from_json(trim(fpath), cx%natm, nint(cx%at_z), cx%at_r, b, err, &
                               force_cart=(aux .or. cx%opts%force_cartesian))
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
      if (kind < TRC_GUESS_CORE .or. kind > TRC_GUESS_SAP) then
         call refuse(cx, "trc_set_guess: kind "//trim(itoa(int(kind))) &
                     //" -- 0 core, 1 GWH, 2 SAD, 3 given, 4 SADQ, 5 SAP")
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
      logical :: collective, talk
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%scf_ok = .false.
      cx%rimp2_ok = .false.
      status = TRC_ERR_STATE
      if (.not. (cx%have_mol .and. cx%have_basis)) then
         cx%message = "trc_run_scf: set the molecule and the basis first"
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
      if (nelec < 0 .or. nunp > nelec .or. mod(nelec - nunp, 2) /= 0) then
         cx%message = "trc_run_scf: charge and multiplicity do not fit the electron count"
         return
      end if
      nalpha = (nelec + nunp)/2
      nbeta = nelec - nalpha
      if (nalpha > n) then
         cx%message = "trc_run_scf: more electrons than functions"
         return
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
      opts = cx%opts
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
      case default
         cx%message = "trc_run_scf: unknown guess kind"
         if (collective) call comm%finalize()
         return
      end select
      if (allocated(dg)) then
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, dguess=dg, comm=comm, la=cx%la)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, dguess=dg, la=cx%la)
         end if
      else if (allocated(hsap)) then
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, comm=comm, hguess=hsap, la=cx%la)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, hguess=hsap, la=cx%la)
         end if
      else
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, comm=comm, la=cx%la)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, la=cx%la)
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
      status = merge(TRC_OK, TRC_ERR_NOCONV, cx%res%converged)
   end function trc_run_scf

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
      !! Needs a converged closed-shell SCF in the context and an auxiliary basis.
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_pairlist_t) :: pl
      type(trc_rimp2_result_t) :: res
      type(comm_t) :: comm
      integer :: nfrozen, nblk, i
      logical :: collective
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
      status = TRC_ERR_UNSUPPORTED
      if (cx%res%nspin /= 1 .or. cx%nalpha /= cx%nbeta) then
         cx%message = "trc_run_rimp2: restricted closed-shell reference only"
         return
      end if
      ! Density fitting and RI-MP2 are out of scope for the spherical AO
      ! space (this milestone's own boundary: trc_df_2c/trc_df_3c stay
      ! Cartesian-only, see this file's AO classification header). `cx%res
      ! %cmo` is AO-row-indexed in whatever space the SCF ran in, so an
      ! MO built from a spherical orbital basis does not line up with
      ! trc_rimp2_run's Cartesian three-centre tensor -- refused rather
      ! than silently contracting the wrong AO index against the right one.
      if (cx%bb%b%spherical) then
         cx%message = "trc_run_rimp2: the orbital basis is spherical; RI-MP2 only "// &
                      "supports a Cartesian basis today -- set force_cartesian"
         return
      end if
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
      if (nfrozen >= cx%nalpha) then
         cx%message = "trc_run_rimp2: every occupied orbital is frozen"
         return
      end if
      nblk = cx%aux_block
      if (nblk <= 0) nblk = cx%ab%b%nao
      call trc_build_comm(cx%fcomm, comm, collective, status)
      if (status /= TRC_OK) return
      call pl%build(cx%bb%b, cx%opts%eri_thresh)
      call pl%to_device()
      if (collective) then
         call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cx%res%cmo(:, :, 1), cx%res%eps(:, 1), res, &
                            nfrozen=nfrozen, aux_block=nblk, comm=comm, verbose=cx%opts%verbose, la=cx%la)
         call comm%finalize()
      else
         call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cx%res%cmo(:, :, 1), cx%res%eps(:, 1), res, &
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
      !! The SCF total, plus the RI-MP2 correlation if that ran.
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: energy
      integer(c_int) :: status
      type(context_box), pointer :: cx
      energy = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      energy = real(cx%res%energy, c_double)
      if (cx%rimp2_ok) energy = energy + real(cx%e_os + cx%e_ss, c_double)
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
