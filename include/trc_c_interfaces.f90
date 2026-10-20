!
! The authoritative Fortran declaration of terco's C ABI -- GENERATED,
! do not edit by hand.
!
! terco is nvfortran and its callers usually are not, so they reach it
! through `bind(c)`. tools/gen_capi_headers.py parses src/api/trc_capi.F90's
! own `bind(c, name=...)` procedures and public TRC_* parameters and
! writes this file and include/terco.h from the same pass, so a Fortran
! caller, a C caller and the library itself cannot drift from each other.
! Regenerate with that script; do not edit the output.
!
! `test/link_gfortran.F90` compiles this file with gfortran and asserts an
! SCF energy, so a drift anywhere in the chain is a wrong number in terco's
! own test suite before any host code sees it.
!
module trc_c_interfaces
   use, intrinsic :: iso_c_binding
   implicit none
   public

   integer(c_int), parameter :: TRC_OK = 0
   integer(c_int), parameter :: TRC_ERR_NULL = 1
   integer(c_int), parameter :: TRC_ERR_BADARG = 2
   integer(c_int), parameter :: TRC_ERR_UNSUPPORTED = 3
   integer(c_int), parameter :: TRC_ERR_NOCONV = 4
   integer(c_int), parameter :: TRC_ERR_STATE = 5
   integer(c_int), parameter :: TRC_GUESS_CORE = 0
   integer(c_int), parameter :: TRC_GUESS_GWH = 1
   integer(c_int), parameter :: TRC_GUESS_SAD = 2
   integer(c_int), parameter :: TRC_GUESS_GIVEN = 3
   integer(c_int), parameter :: TRC_GUESS_SADQ = 4
   integer(c_int), parameter :: TRC_GUESS_SAP = 5
   integer(c_int), parameter :: TRC_AO_NORM_SHELL = 0
   integer(c_int), parameter :: TRC_AO_NORM_UNIT = 1

   interface

      integer(c_int) function trc_basis_create(nshell, maxnp, l, nprim, exps, coefs, centres, natm, zatm, ratm, handle &
                                               ) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         integer(c_int), value :: nshell
         integer(c_int), value :: maxnp
         integer(c_int), intent(in) :: l(nshell)
         integer(c_int), intent(in) :: nprim(nshell)
         real(c_double), intent(in) :: exps(maxnp, nshell)
         real(c_double), intent(in) :: coefs(maxnp, nshell)
         real(c_double), intent(in) :: centres(3, nshell)
         integer(c_int), value :: natm
         real(c_double), intent(in) :: zatm(natm)
         real(c_double), intent(in) :: ratm(3, natm)
         type(c_ptr), intent(out) :: handle
      end function trc_basis_create

      integer(c_int) function trc_basis_destroy(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_basis_destroy

      integer(c_int) function trc_basis_nao(handle, nao) bind(c) result(res)
         !! Number of atomic orbitals, so a caller can size its own output buffers.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nao
      end function trc_basis_nao

      integer(c_int) function trc_pairs_create(basis, thresh, handle) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: basis
         real(c_double), value :: thresh
         type(c_ptr), intent(out) :: handle
      end function trc_pairs_create

      integer(c_int) function trc_pairs_destroy(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_pairs_destroy

      integer(c_int) function trc_pairs_count(handle, npair) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: npair
      end function trc_pairs_count

      integer(c_int) function trc_compute_1e(basis, pairs, smat, tmat, vmat) bind(c) result(res)
         !! COVARIANT under `basis`'s own ao_norm_kind (set through
         !! trc_set_ao_norm on a context that lent this handle via
         !! trc_context_basis; a raw trc_basis_create* handle is always
         !! TRC_AO_NORM_SHELL and this is then a no-op): S_out = N S N, same
         !! for T and V.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: basis
         type(c_ptr), value :: pairs
         real(c_double), intent(out) :: smat(*)
         real(c_double), intent(out) :: tmat(*)
         real(c_double), intent(out) :: vmat(*)
      end function trc_compute_1e

      integer(c_int) function trc_compute_df2c(aux, jmat) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: aux
         real(c_double), intent(out) :: jmat(*)
      end function trc_compute_df2c

      integer(c_int) function trc_compute_df3c(basis, pairs, aux, tens) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: basis
         type(c_ptr), value :: pairs
         type(c_ptr), value :: aux
         real(c_double), intent(out) :: tens(*)
      end function trc_compute_df3c

      integer(c_int) function trc_basis_create_libcint(atm, natm, bas, nshell, env, nenv, cartesian, handle) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nshell
         real(c_double), intent(in) :: env(*)
         integer(c_int), value :: nenv
         integer(c_int), value :: cartesian
         type(c_ptr), intent(out) :: handle
      end function trc_basis_create_libcint

      integer(c_int) function trc_basis_maxl(handle, maxl) bind(c) result(res)
         !! Highest angular momentum in the basis, so a caller can decide whether to
         !! use this library at all before it builds anything expensive. The
         !! four-centre kernels stop at d and the caller has to fall back rather than
         !! get a wrong answer.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: maxl
      end function trc_basis_maxl

      integer(c_int) function trc_max_angular_momentum(eri_lmax, df_lmax, df_lmax_aux) bind(c) result(res)
         !! The angular-momentum ceilings THIS build was compiled for, so a host
         !! can check a basis against them instead of hard-coding them: the
         !! four-centre ERIs / SCF (TERCO_LMAX: 2 by default, 3 for f), and the
         !! density-fitting orbital and auxiliary bases (TERCO_DF_LMAX; the
         !! auxiliary one is one higher). Compare with trc_basis_maxl.
         import :: c_int
         implicit none
         integer(c_int), intent(out) :: eri_lmax
         integer(c_int), intent(out) :: df_lmax
         integer(c_int), intent(out) :: df_lmax_aux
      end function trc_max_angular_momentum

      integer(c_int) function trc_multipole_count() bind(c) result(res)
         !! Number of Cartesian multipole components, so a caller can size its buffer
         !! without hardcoding 39 and drifting when the operator order changes.
         !! Named for what it counts rather than `trc_nmult`, which collides with the
         !! parameter of that name it returns.
         import :: c_int
         implicit none
      end function trc_multipole_count

      integer(c_int) function trc_compute_multipoles(basis, pairs, origin, mmat) bind(c) result(res)
         !! COVARIANT under `basis`'s own ao_norm_kind, one component at a
         !! time (see trc_compute_1e).
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: basis
         type(c_ptr), value :: pairs
         real(c_double), intent(in) :: origin(3)
         real(c_double), intent(out) :: mmat(*)
      end function trc_compute_multipoles

      integer(c_int) function trc_eri_create(basis, thresh, handle) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: basis
         real(c_double), value :: thresh
         type(c_ptr), intent(out) :: handle
      end function trc_eri_create

      integer(c_int) function trc_eri_destroy(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_eri_destroy

      integer(c_int) function trc_fock(eri, basis, dmat, gmat, jfac, kfac, dscreen, omega) bind(c) result(res)
         !! CONTRAVARIANT in (dmat), COVARIANT out (gmat), under `basis`'s own
         !! ao_norm_kind: this is GAMESS's own main caller (its Fock/J/K
         !! build, run_Fockbuild's role in hpcchem_driver.src), so it is wired
         !! rather than left at TRC_AO_NORM_SHELL like the entries the "not
         !! wired" list below still names. See fock_scale_in/fock_scale_out
         !! for why SHELL mode still costs nothing extra.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: eri
         type(c_ptr), value :: basis
         real(c_double), intent(in) :: dmat(*)
         real(c_double), intent(out) :: gmat(*)
         real(c_double), value :: jfac
         real(c_double), value :: kfac
         integer(c_int), value :: dscreen
         real(c_double), value :: omega
      end function trc_fock

      integer(c_int) function trc_fock_many(eri, basis, ndens, dmats, gmats, jfac, kfac, dscreen, omega) bind(c) result(res)
         !! CONTRAVARIANT in, COVARIANT out, under `basis`'s own ao_norm_kind
         !! -- see capi_fock; the density-index-fastest layout means the
         !! UNIT-mode scaling here is its own small routines
         !! (fock_many_scale_in/out) rather than a reuse of fock_scale_in/out.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: eri
         type(c_ptr), value :: basis
         integer(c_int), value :: ndens
         real(c_double), intent(in) :: dmats(*)
         real(c_double), intent(out) :: gmats(*)
         real(c_double), value :: jfac
         real(c_double), value :: kfac
         integer(c_int), value :: dscreen
         real(c_double), value :: omega
      end function trc_fock_many

      integer(c_int) function trc_fock_nosym(eri, basis, dmat, gmat, dscreen) bind(c) result(res)
         !! CONTRAVARIANT in, COVARIANT out, exactly as capi_fock -- the
         !! diagonal congruence ao_scale_cov/fock_scale_in apply does not
         !! care whether `dmat` is symmetric, so nothing here differs from
         !! the symmetric entry beyond which shared routines it reuses.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: eri
         type(c_ptr), value :: basis
         real(c_double), intent(in) :: dmat(*)
         real(c_double), intent(out) :: gmat(*)
         integer(c_int), value :: dscreen
      end function trc_fock_nosym

      integer(c_int) function trc_scf(basis, nalpha, nbeta, functional, grid_level, conv_energy, conv_density, max_iter, dguess, &
                                      verbose, energy, e_xc, dmat, eps, niter) bind(c) result(res)
         import :: c_char, c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: basis
         integer(c_int), value :: nalpha
         integer(c_int), value :: nbeta
         character(kind=c_char), intent(in) :: functional(*)
         integer(c_int), value :: grid_level
         real(c_double), value :: conv_energy
         real(c_double), value :: conv_density
         integer(c_int), value :: max_iter
         type(c_ptr), value :: dguess
         integer(c_int), value :: verbose
         real(c_double), intent(out) :: energy
         real(c_double), intent(out) :: e_xc
         real(c_double), intent(out) :: dmat(*)
         real(c_double), intent(out) :: eps(*)
         integer(c_int), intent(out) :: niter
      end function trc_scf

      integer(c_int) function trc_scf_mpi(fcomm, basis, nalpha, nbeta, functional, grid_level, conv_energy, conv_density, &
                                          max_iter, dguess, verbose, energy, e_xc, dmat, eps, niter) bind(c) result(res)
         import :: c_char, c_double, c_int, c_ptr
         implicit none
         integer(c_int), value :: fcomm
         type(c_ptr), value :: basis
         integer(c_int), value :: nalpha
         integer(c_int), value :: nbeta
         character(kind=c_char), intent(in) :: functional(*)
         integer(c_int), value :: grid_level
         real(c_double), value :: conv_energy
         real(c_double), value :: conv_density
         integer(c_int), value :: max_iter
         type(c_ptr), value :: dguess
         integer(c_int), value :: verbose
         real(c_double), intent(out) :: energy
         real(c_double), intent(out) :: e_xc
         real(c_double), intent(out) :: dmat(*)
         real(c_double), intent(out) :: eps(*)
         integer(c_int), intent(out) :: niter
      end function trc_scf_mpi

      integer(c_int) function trc_rimp2(fcomm, basis, aux, nfrozen, aux_block, e_os, e_ss) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         integer(c_int), value :: fcomm
         type(c_ptr), value :: basis
         type(c_ptr), value :: aux
         integer(c_int), value :: nfrozen
         integer(c_int), value :: aux_block
         real(c_double), intent(out) :: e_os
         real(c_double), intent(out) :: e_ss
      end function trc_rimp2

      integer(c_int) function trc_create(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(out) :: handle
      end function trc_create

      integer(c_int) function trc_destroy(handle) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_destroy

      integer(c_int) function trc_set_molecule(handle, natm, z, xyz, charge, multiplicity) bind(c) result(res)
         !! Atoms in Bohr, nuclear charges as doubles, atom-major coordinates.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natm
         real(c_double), intent(in) :: z(*)
         real(c_double), intent(in) :: xyz(*)
         integer(c_int), value :: charge
         integer(c_int), value :: multiplicity
      end function trc_set_molecule

      integer(c_int) function trc_set_basis_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) bind(c) result(res)
         !! The orbital basis from libcint's packed `atm`/`bas`/`env`.
         !! `cartesian` nonzero builds terco's own Cartesian AO space; zero
         !! tags it spherical (libcint/pyscf order, m = -l..l except p as x,
         !! y, z) -- trc_nao and every SCF-facing entry then size and index
         !! at nao_sph, exactly as for a spherical trc_set_basis_json. A
         !! mismatch against an already-set auxiliary or jkfit basis's own
         !! convention is refused, not coerced.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nshell
         real(c_double), intent(in) :: env(*)
         integer(c_int), value :: nenv
         integer(c_int), value :: cartesian
      end function trc_set_basis_libcint

      integer(c_int) function trc_set_aux_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) bind(c) result(res)
         !! The auxiliary (RI-MP2 rifit) basis, same layout and `cartesian`
         !! convention as trc_set_basis_libcint; refused rather than coerced
         !! if it disagrees with the orbital basis already set.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nshell
         real(c_double), intent(in) :: env(*)
         integer(c_int), value :: nenv
         integer(c_int), value :: cartesian
      end function trc_set_aux_libcint

      integer(c_int) function trc_set_basis_arrays(handle, nshell, maxnp, l, nprim, exps, coefs, centres) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: nshell
         integer(c_int), value :: maxnp
         integer(c_int), intent(in) :: l(*)
         integer(c_int), intent(in) :: nprim(*)
         real(c_double), intent(in) :: exps(*)
         real(c_double), intent(in) :: coefs(*)
         real(c_double), intent(in) :: centres(*)
      end function trc_set_basis_arrays

      integer(c_int) function trc_set_aux_arrays(handle, nshell, maxnp, l, nprim, exps, coefs, centres) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: nshell
         integer(c_int), value :: maxnp
         integer(c_int), intent(in) :: l(*)
         integer(c_int), intent(in) :: nprim(*)
         real(c_double), intent(in) :: exps(*)
         real(c_double), intent(in) :: coefs(*)
         real(c_double), intent(in) :: centres(*)
      end function trc_set_aux_arrays

      integer(c_int) function trc_set_basis_json(handle, path) bind(c) result(res)
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_set_basis_json

      integer(c_int) function trc_set_aux_json(handle, path) bind(c) result(res)
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_set_aux_json

      integer(c_int) function trc_set_method(handle, functional, grid_level) bind(c) result(res)
         !! An empty functional name is Hartree-Fock.
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: functional(*)
         integer(c_int), value :: grid_level
      end function trc_set_method

      integer(c_int) function trc_set_nlc_grid(handle, grid_level) bind(c) result(res)
         !! Grid level for VV10's nonlocal correlation (default 1, coarser than
         !! the semilocal grid: VV10 is a double sum over its points). Only
         !! used by a functional with VV10, such as wB97X-V.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: grid_level
      end function trc_set_nlc_grid

      integer(c_int) function trc_set_pcm(handle, model, eps, lebedev, radii_scale) bind(c) result(res)
         import :: c_char, c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: model(*)
         real(c_double), value :: eps
         integer(c_int), value :: lebedev
         real(c_double), value :: radii_scale
      end function trc_set_pcm

      integer(c_int) function trc_pcm_energy(handle, e_diel, q_total, npoints) bind(c) result(res)
         !! The continuum's part of the last SCF: E_diel = q.v/2 (already in
         !! trc_energy's total), the total apparent surface charge, and the
         !! surface points the cavity kept. Zeros when no continuum was set.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: e_diel
         real(c_double), intent(out) :: q_total
         integer(c_int), intent(out) :: npoints
      end function trc_pcm_energy

      integer(c_int) function trc_set_convergence(handle, conv_energy, conv_diis, max_iter, ndiis) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), value :: conv_energy
         real(c_double), value :: conv_diis
         integer(c_int), value :: max_iter
         integer(c_int), value :: ndiis
      end function trc_set_convergence

      integer(c_int) function trc_set_screening(handle, thresh) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), value :: thresh
      end function trc_set_screening

      integer(c_int) function trc_set_incremental(handle, on) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: on
      end function trc_set_incremental

      integer(c_int) function trc_set_df(handle, on) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: on
      end function trc_set_df

      integer(c_int) function trc_set_jkfit_json(handle, path) bind(c) result(res)
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_set_jkfit_json

      integer(c_int) function trc_set_jkfit_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) bind(c) result(res)
         !! The --df fitting basis from libcint's packed arrays, the same
         !! shape and `cartesian` convention as trc_set_aux_libcint; refused
         !! rather than coerced if it disagrees with the orbital basis
         !! already set (pyscf's auxmol.cart = mol.cart rule).
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value :: nshell
         real(c_double), intent(in) :: env(*)
         integer(c_int), value :: nenv
         integer(c_int), value :: cartesian
      end function trc_set_jkfit_libcint

      integer(c_int) function trc_set_batching(handle, res) bind(c) result(res_)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: res
      end function trc_set_batching

      integer(c_int) function trc_set_force_cartesian(handle, flag) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: flag
      end function trc_set_force_cartesian

      integer(c_int) function trc_set_guess(handle, kind, dguess, nspin) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: kind
         type(c_ptr), value :: dguess
         integer(c_int), value :: nspin
      end function trc_set_guess

      integer(c_int) function trc_set_ao_norm(handle, kind) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: kind
      end function trc_set_ao_norm

      integer(c_int) function trc_get_ao_norm(handle, kind) bind(c) result(res)
         !! The AO normalisation this context is set to -- TRC_AO_NORM_SHELL
         !! unless trc_set_ao_norm was called with TRC_AO_NORM_UNIT.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: kind
      end function trc_get_ao_norm

      integer(c_int) function trc_set_comm(handle, fcomm) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: fcomm
      end function trc_set_comm

      integer(c_int) function trc_set_verbose(handle, level) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: level
      end function trc_set_verbose

      integer(c_int) function trc_set_rimp2(handle, nfrozen, aux_block) bind(c) result(res)
         !! `nfrozen` -1 counts the core from the elements; `aux_block` <= 0 is
         !! the whole auxiliary basis at once.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: nfrozen
         integer(c_int), value :: aux_block
      end function trc_set_rimp2

      integer(c_int) function trc_run_scf(handle) bind(c) result(res)
         !! Needs the molecule and the basis. Builds the guess it was asked for,
         !! runs, keeps everything in the context.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_run_scf

      integer(c_int) function trc_run_gradient(handle) bind(c) result(res)
         !! ======================================================================
         !! !!! EXPERIMENTAL !!! Validated against pyscf and finite differences,
         !! but NOT performance-tuned: the --df gradient can cost far more than
         !! the energy, and only one rank is supported. Do not rely on its cost
         !! yet.
         !! ======================================================================
         !!
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
         !! it on). wB97M-V is meta-GGA and cannot reach this point at all (see
         !! below), so it is not a case this routine has to decide about.
         !!
         !! REFUSED, clearly, rather than attempted:
         !!   - meta-GGA or a double hybrid: `trc_xc_functional_t` has no
         !!     meta-GGA family and no MP2-fraction field on this branch, so
         !!     neither can reach this point at all -- `trc_set_method` already
         !!     refused it (func%missing) before a functional name naming one
         !!     could be stored in `cx%opts%functional`. Nothing to check here.
         !!   - RI-MP2 (or, once it exists, a double hybrid): no correlation
         !!     gradient exists; `trc_scf_gradient` is the SCF reference term
         !!     only.
         !!   - --df outside the DF-SCF's own scope (unrestricted, spherical,
         !!     range-separated -- the SCF refuses those first), or with the fit
         !!     gone. A --df gradient otherwise runs the FITTED J/K gradient
         !!     (trc_scf_grad_df.F90) over `cx%dfjk`, the fit the SCF used.
         !!   - PCM: no cavity/surface-charge derivative exists.
         !!   - a multi-rank context (`cx%fcomm /= -1`). `trc_scf_gradient`
         !!     takes no communicator: it rebuilds the pair list and ERI
         !!     engine from `cx%bb%b` alone, so on a collective run it would
         !!     rebuild the FULL (unsplit) pair list, Schwarz screen and ERI
         !!     structures on EVERY rank's own device -- not merely redundant
         !!     work, but re-materialising, on one GPU, exactly the memory a
         !!     multi-rank SCF split across ranks to avoid. `eri%gradient`
         !!     itself DOES support a rank split (`this%rank`/`this%nranks`,
         !!     `sum_ranks`'s allreduce) but only when its own `eri%build` is
         !!     given a `comm` -- a plumbing `trc_scf_gradient` does not do
         !!     today, found by reading trc_eri.F90's `eri_gradient` and
         !!     comparing it against `eri_build`'s own comm handling.
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_run_gradient

      integer(c_int) function trc_gradient(handle, grad) bind(c) result(res)
         !! (3, natm), Eh/Bohr, column-major (x, y, z fastest, then atom, the
         !! same atom order `trc_set_molecule` was given) -- after a
         !! successful `trc_run_gradient`. NOT AO-indexed data (see
         !! trc_capi.F90's own "AO NORMALISATION" header): a nuclear gradient
         !! is a derivative of the SCF energy, which does not depend on which
         !! AO normalisation convention computed it, so `trc_set_ao_norm`
         !! plays no role here.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: grad(*)
      end function trc_gradient

      integer(c_int) function trc_run_rimp2(handle) bind(c) result(res)
         !! Needs a converged SCF in the context (restricted, or unrestricted --
         !! UHF or UKS, cx%res%nspin == 2) and an auxiliary basis. An
         !! unrestricted reference refuses more than one rank below (see
         !! trc_urimp2_run's own header for why the restricted ring is not a
         !! clean generalisation here). The orbital basis may be spherical:
         !! `cx%res%cmo`/`eps` (at whatever space the SCF ran in) are folded
         !! into the Cartesian, zero-padded space `trc_rimp2_run`/
         !! `trc_urimp2_run` themselves run in by `rimp2_mo_to_cart`
         !! (trc_rimp2_driver.F90's own header says why the pad is exact).
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_run_rimp2

      integer(c_int) function trc_save_scf(handle, path) bind(c) result(res)
         !! Write the context's SCF result -- orbitals, orbital energies,
         !! density, energies -- to `path` (NUL-terminated), with the molecule,
         !! basis and method it belongs to; trc_load_scf reads it back. One
         !! rank is enough: every rank holds the same result.
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_save_scf

      integer(c_int) function trc_load_scf(handle, path) bind(c) result(res)
         !! In place of trc_run_scf: take the SCF result from a file
         !! trc_save_scf wrote, so a correlated step (trc_run_rimp2) can follow
         !! without repeating the SCF. The molecule, basis, method and solvent
         !! must be set as for trc_run_scf, and must be the ones the file was
         !! made with; anything else is refused, naming the difference. Every
         !! rank reads the file itself. TRC_ERR_NOCONV if the saved SCF had not
         !! converged, as trc_run_scf would have said.
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_load_scf

      integer(c_int) function trc_cart_order() bind(c) result(res)
         !! 0 (lexicographic) or 1 (gamess): the within-shell Cartesian order every
         !! AO-indexed array this API exchanges uses, fixed at compile time
         !! (TRC_CART_ORDER) and the same for every handle -- no handle argument,
         !! and it always succeeds.
         import :: c_int
         implicit none
      end function trc_cart_order

      integer(c_int) function trc_nao(handle, nao) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nao
      end function trc_nao

      integer(c_int) function trc_is_spherical(handle, flag) bind(c) result(res)
         !! Non-zero if the orbital basis this context holds is running in
         !! the spherical AO space (nao_sph of trc_nao, rather than the
         !! Cartesian nao) -- the basis's own JSON function_type unless
         !! force_cartesian overrode it. For a caller checking which space
         !! `trc_nao`'s number, and every AO-indexed array it sizes, is in.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: flag
      end function trc_is_spherical

      integer(c_int) function trc_nspin(handle, nspin) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nspin
      end function trc_nspin

      integer(c_int) function trc_energy(handle, energy) bind(c) result(res)
         !! The SCF total, plus the RI-MP2 correlation if that ran: the FULL
         !! E_os + E_ss for a plain post-SCF RI-MP2 (any functional, or
         !! Hartree-Fock), or c_os E_os + c_ss E_ss for a double hybrid --
         !! trc_rimp2_energy below always reports the raw, unscaled E_os/E_ss
         !! either way, so a caller after the PT2 part alone never sees this
         !! weighting. A caller asking for a double hybrid's total needs
         !! nothing beyond the ordinary trc_run_scf/trc_run_rimp2/trc_energy
         !! sequence any RI-MP2 run already uses -- the weighting is applied
         !! here, once, from the functional's own pt2_cos/pt2_css.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: energy
      end function trc_energy

      integer(c_int) function trc_energy_parts(handle, e_nuc, e_one, e_two, e_xc) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: e_nuc
         real(c_double), intent(out) :: e_one
         real(c_double), intent(out) :: e_two
         real(c_double), intent(out) :: e_xc
      end function trc_energy_parts

      integer(c_int) function trc_converged(handle, flag) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: flag
      end function trc_converged

      integer(c_int) function trc_iterations(handle, niter) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: niter
      end function trc_iterations

      integer(c_int) function trc_density(handle, dmat) bind(c) result(res)
         !! (nao, nao, nspin), column-major; restricted holds the total.
         !!
         !! CONTRAVARIANT under trc_set_ao_norm: D_out = N^-1 D_terco N^-1 per
         !! spin slice when UNIT. Scaled on a COPY -- `cx%res%dmat` itself
         !! stays in terco's own convention, because trc_run_rimp2 and a
         !! second trc_run_scf both read it as such internally.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: dmat(*)
      end function trc_density

      integer(c_int) function trc_mo_coefficients(handle, cmo) bind(c) result(res)
         !! CONTRAVARIANT under trc_set_ao_norm, on the AO (row) index only:
         !! C_out = N^-1 C_terco per spin slice when UNIT -- see trc_density
         !! for why this scales a copy, not `cx%res%cmo` itself.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: cmo(*)
      end function trc_mo_coefficients

      integer(c_int) function trc_mo_energies(handle, eps) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: eps(*)
      end function trc_mo_energies

      integer(c_int) function trc_rimp2_energy(handle, e_os, e_ss) bind(c) result(res)
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: e_os
         real(c_double), intent(out) :: e_ss
      end function trc_rimp2_energy

      integer(c_int) function trc_message(handle, buf, buflen) bind(c) result(res)
         !! The last message, NUL-terminated into `buf` of `buflen` chars.
         import :: c_char, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         character(kind=c_char), intent(out) :: buf(*)
         integer(c_int), value :: buflen
      end function trc_message

      integer(c_int) function trc_context_basis(handle, basis) bind(c) result(res)
         !! The basis handle inside the context: BORROWED, never destroyed by
         !! the caller, valid until the context is destroyed or the basis reset.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), intent(out) :: basis
      end function trc_context_basis

      integer(c_int) function trc_context_eri(handle, eri) bind(c) result(res)
         !! The ERI context at the context's screening threshold, built on first
         !! use. BORROWED, as trc_context_basis.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), intent(out) :: eri
      end function trc_context_eri

   end interface

end module trc_c_interfaces
