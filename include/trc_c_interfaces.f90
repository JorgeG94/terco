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
   integer(c_int), parameter :: TRC_GUESS_PROJECTED = 6
   integer(c_int), parameter :: TRC_AO_NORM_SHELL = 0
   integer(c_int), parameter :: TRC_AO_NORM_UNIT = 1
   integer(c_int), parameter :: TRC_ABI_MAJOR = 1
   integer(c_int), parameter :: TRC_ABI_MINOR = 6
   integer(c_int), parameter :: TRC_CONST_BOHR_TO_ANGSTROM = 0
   integer(c_int), parameter :: TRC_CONST_HARTREE_TO_EV = 1
   integer(c_int), parameter :: TRC_CONST_HARTREE_TO_KCALMOL = 2
   integer(c_int), parameter :: TRC_CONST_AMU_TO_ME = 3
   integer(c_int), parameter :: TRC_CONST_KB_SI = 4
   integer(c_int), parameter :: TRC_CONST_NA_SI = 5

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

      integer(c_int) function trc_basis_shell_count(handle, nshell) bind(c) result(res)
         !! Number of shells, so a caller can size trc_basis_shell_info's
         !! buffers. NOT split by decontraction (src/basis/trc_decontract.F90):
         !! that split lives inside trc_eri_t, built from a copy of this
         !! basis, and never touches the trc_basis_t this handle itself
         !! holds -- this is always the basis as the caller (or
         !! trc_basis_from_json) built it.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nshell
      end function trc_basis_shell_count

      integer(c_int) function trc_basis_shell_max_primitives(handle, maxnp) bind(c) result(res)
         !! The second dimension trc_basis_shell_primitives' exps/coefs
         !! buffers need (MAXNP*NSHELL each): the widest shell's own
         !! primitive count, box%b%maxnp -- the same bound
         !! trc_basis_create's own `maxnp` argument sets on the way in.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: maxnp
      end function trc_basis_shell_max_primitives

      integer(c_int) function trc_basis_shell_info(handle, atom, l, nprim) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: atom(*)
         integer(c_int), intent(out) :: l(*)
         integer(c_int), intent(out) :: nprim(*)
      end function trc_basis_shell_info

      integer(c_int) function trc_basis_shell_primitives(handle, exps, coefs) bind(c) result(res)
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
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: exps(*)
         real(c_double), intent(out) :: coefs(*)
      end function trc_basis_shell_primitives

      integer(c_int) function trc_basis_is_spherical(handle, spherical) bind(c) result(res)
         !! Non-zero if this basis handle's shells are spherical rather than
         !! Cartesian -- box%b%spherical, set by trc_basis_from_json from a
         !! JSON basis's own function_type, the same field trc_is_spherical
         !! reads off a context's current basis (see that entry's own
         !! comment for what it changes downstream). This is the basis-
         !! handle-level answer, for a caller that holds the handle directly
         !! rather than the context it came from -- trc_context_basis hands
         !! out exactly that handle.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: spherical
      end function trc_basis_is_spherical

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
         !! The four-centre engine for `basis` at screening threshold `thresh`.
         !! TRC_ERR_UNSUPPORTED, and no handle, when the basis has a shell
         !! above this build's four-centre ceiling (trc_max_angular_momentum).
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

      integer(c_int) function trc_set_charge(handle, charge, multiplicity) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: charge
         integer(c_int), value :: multiplicity
      end function trc_set_charge

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

      integer(c_int) function trc_set_unrestricted(handle, on) bind(c) result(res)
         !! Force UHF/UKS (two independent spin densities) even on a CLOSED
         !! shell (nalpha == nbeta), which otherwise runs restricted --
         !! `trc_scf_options_t%unrestricted`'s own comment, "forced; an open
         !! shell is unrestricted regardless": nalpha /= nbeta already gets
         !! nspin = 2 whether or not this is set, so this only changes
         !! anything on a singlet. `on` nonzero sets it; 0 (the default)
         !! restores restricted-when-closed-shell.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: on
      end function trc_set_unrestricted

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

      integer(c_int) function trc_set_pcm_radii(handle, natm, radii) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: natm
         type(c_ptr), value :: radii
      end function trc_set_pcm_radii

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

      integer(c_int) function trc_set_level_shift(handle, shift) bind(c) result(res)
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
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), value :: shift
      end function trc_set_level_shift

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

      integer(c_int) function trc_set_geometry(handle, xyz) bind(c) result(res)
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
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(in) :: xyz(*)
      end function trc_set_geometry

      integer(c_int) function trc_run_gradient(handle) bind(c) result(res)
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
         !! it on). wB97M-V is not in terco's table, and a meta-GGA is refused
         !! below, so it is not a case this routine has to decide about.
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

      integer(c_int) function trc_set_hessian_step(handle, step) bind(c) result(res)
         !! The Cartesian displacement trc_run_hessian's central differences
         !! use, in Bohr; must be positive. Default 0.005: its truncation
         !! error, about 1e-5 Eh/Bohr^2 on water/STO-3G, sits well above the
         !! SCF noise trc_run_hessian's tightened convergence leaves.
         !! Unchanged until this is called. Invalidates the last
         !! trc_run_hessian result (a Hessian at the OLD step is not a
         !! Hessian at the new one), not the SCF/gradient it is built from.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), value :: step
      end function trc_set_hessian_step

      integer(c_int) function trc_hessian_step(handle, step) bind(c) result(res)
         !! The displacement trc_run_hessian would use (or last used) --
         !! trc_set_hessian_step's own default or override.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: step
      end function trc_hessian_step

      integer(c_int) function trc_run_hessian(handle) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_run_hessian

      integer(c_int) function trc_hessian(handle, hess) bind(c) result(res)
         !! (3*natm, 3*natm), Eh/Bohr^2, column-major, row/column k =
         !! 3*(atom-1)+{x,y,z} (trc_set_molecule's own atom order) -- after a
         !! successful trc_run_hessian. NOT AO-indexed data, same as
         !! trc_gradient: a Hessian is a second derivative of the SCF
         !! energy, independent of any AO normalisation convention.
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: hess(*)
      end function trc_hessian

      integer(c_int) function trc_run_polarizability(handle) bind(c) result(res)
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
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
      end function trc_run_polarizability

      integer(c_int) function trc_polarizability(handle, alpha) bind(c) result(res)
         !! The static dipole polarizability tensor `trc_run_polarizability`
         !! computed, a.u., column-major (3, 3).
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: alpha(*)
      end function trc_polarizability

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

      integer(c_int) function trc_abi_version(major, minor) bind(c) result(res)
         !! This build's ABI version (TRC_ABI_MAJOR, TRC_ABI_MINOR above), so
         !! a host can compare it against the constants the header it
         !! compiled against carries before trusting any other entry point.
         !! Bump TRC_ABI_MINOR for an additive change (a new entry point, a
         !! new enum value appended at the end); bump TRC_ABI_MAJOR, and
         !! reset TRC_ABI_MINOR to 0, for anything that changes an existing
         !! signature or an existing enum value's meaning. No handle; always
         !! succeeds.
         import :: c_int
         implicit none
         integer(c_int), intent(out) :: major
         integer(c_int), intent(out) :: minor
      end function trc_abi_version

      integer(c_int) function trc_select_device(rank, dev) bind(c) result(res)
         !! Bind this process to a GPU by its node-local rank: device
         !! mod(rank, number of devices), what trc_api's trc_bind_device does
         !! for terco's own MPI runs. A host code that runs one single-rank
         !! context per process and never calls trc_set_comm (one fragment per
         !! rank, say) calls this once, before the first trc_create, or every
         !! rank lands on device 0. `dev` is the device chosen, -1 on a build
         !! without the GPU or a node without one. Process-wide; no handle;
         !! always succeeds.
         import :: c_int
         implicit none
         integer(c_int), value :: rank
         integer(c_int), intent(out) :: dev
      end function trc_select_device

      integer(c_int) function trc_select_device_env(dev) bind(c) result(res)
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
         import :: c_int
         implicit none
         integer(c_int), intent(out) :: dev
      end function trc_select_device_env

      integer(c_int) function trc_physical_constant(id, value) bind(c) result(res)
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
         import :: c_double, c_int
         implicit none
         integer(c_int), value :: id
         real(c_double), intent(out) :: value
      end function trc_physical_constant

      integer(c_int) function trc_functional_info(name, flags, exx, c_os, c_ss) bind(c) result(res)
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
         import :: c_char, c_double, c_int
         implicit none
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), intent(out) :: flags
         real(c_double), intent(out) :: exx
         real(c_double), intent(out) :: c_os
         real(c_double), intent(out) :: c_ss
      end function trc_functional_info

      integer(c_int) function trc_functional_list(buf, buflen, reqlen) bind(c) result(res)
         !! XC_AVAILABLE_FUNCTIONALS -- every name trc_set_method/
         !! trc_functional_info can look up that terco can actually run today
         !! -- space-separated and NUL-terminated into `buf` (capacity
         !! `buflen` bytes, NUL included). `reqlen` is ALWAYS set to the
         !! number of bytes that call needed, NUL included, whether or not it
         !! fit: a caller queries it with `buflen = 0` (`buf` is never
         !! touched then) and allocates exactly `reqlen` bytes for a second,
         !! filling call. Returns TRC_ERR_BADARG, leaving `buf` untouched, if
         !! `buflen` is positive but smaller than `reqlen`.
         import :: c_char, c_int
         implicit none
         character(kind=c_char), intent(out) :: buf(*)
         integer(c_int), value :: buflen
         integer(c_int), intent(out) :: reqlen
      end function trc_functional_list

      integer(c_int) function trc_nao(handle, nao) bind(c) result(res)
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nao
      end function trc_nao

      integer(c_int) function trc_nmo(handle, nmo) bind(c) result(res)
         !! Molecular-orbital count: `trc_nao`'s own AO count minus whatever
         !! canonical orthogonalisation dropped as linearly dependent at the
         !! last `trc_run_scf` (`trc_scf_result_t%n_lindep` -- 0 for every
         !! compact basis, nonzero only for a diffuse one close to linear
         !! dependency). Before the first `trc_run_scf` that count is
         !! unknown, so `nmo` is `trc_nao`'s own number.
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nmo
      end function trc_nmo

      integer(c_int) function trc_mo_occupations(handle, occ) bind(c) result(res)
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
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: occ(*)
      end function trc_mo_occupations

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
