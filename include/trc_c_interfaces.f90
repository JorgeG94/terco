!
! The authoritative Fortran declaration of terco's C ABI.
!
! WHY THIS FILE EXISTS
! --------------------
! terco is nvfortran and its callers usually are not, so they reach it through
! `bind(c)`. A Fortran caller then needs an `interface` block for every entry
! point -- and the obvious thing is for each caller to write its own. That is a
! copy of a contract, kept in step by hand, and there is nothing to catch it
! drifting: adding an argument on the library side and forgetting one copy of
! the block still COMPILES and still LINKS, because a C symbol carries no
! signature. It fails at run time, by reading whatever was in the register.
!
! So the declaration lives here, once, in terco's own repository, and every
! consumer compiles this file rather than restating it. This is plain
! `iso_c_binding` -- interface bodies and parameters, no executable code -- so
! any Fortran compiler can build it, which is the point: gfortran, ifx or
! nvfortran all get the same declarations from the same source.
!
! `test/link_gfortran.F90` compiles this file too, with gfortran, and asserts
! an SCF energy. That is what makes the arrangement self-checking: change
! `trc_capi.F90` without changing this file and the energy comes out wrong,
! in terco's own test suite, before any host code sees it.
!
! HOW TO USE IT
! -------------
!     use trc_c_interfaces
!     rc = trc_basis_create_libcint(atm, natm, bas, nbas, env, nenv, 1, h)
!     if (rc /= TRC_OK) ...
!
! ARRAY LAYOUT, WHICH THE DECLARATIONS CANNOT STATE
! --------------------------------------------------
! Everything is Fortran order -- first index fastest. Two cases are worth
! naming because a caller would otherwise guess wrong:
!
!   * `exps` and `coefs` in `trc_basis_create` are `(maxnp, nshell)`. A C
!     caller writing `exps[shell][prim]` has them transposed.
!   * `dmats` and `gmats` in `trc_fock_many` are `(ndens, nao, nao)` -- the
!     DENSITY INDEX IS FIRST. That is not arbitrary: with it last, the N
!     atomic updates for one (mu,nu) sit nao^2 apart and the batch stops
!     paying past N = 4.
!
module trc_c_interfaces
   use, intrinsic :: iso_c_binding
   implicit none
   public

   ! Zero is success; everything else is a reason. A library that stops the
   ! process on bad input is unusable inside someone else's SCF.
   integer(c_int), parameter :: TRC_OK              = 0
   integer(c_int), parameter :: TRC_ERR_NULL        = 1
   integer(c_int), parameter :: TRC_ERR_BADARG      = 2
   integer(c_int), parameter :: TRC_ERR_UNSUPPORTED = 3
   integer(c_int), parameter :: TRC_ERR_NOCONV      = 4   !! the SCF ran out of iterations
   integer(c_int), parameter :: TRC_ERR_STATE = 5   !! called before what it needs was set
   integer(c_int), parameter :: TRC_GUESS_CORE = 0, TRC_GUESS_GWH = 1, TRC_GUESS_SAD = 2, &
                                TRC_GUESS_GIVEN = 3

   interface

      ! --- basis ------------------------------------------------------------

      integer(c_int) function trc_basis_create(nshell, maxnp, l, nprim, exps, &
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
            coefs, centres, natm, zatm, ratm, handle) bind(c)
         !! Build a basis from flat per-shell arrays.
         import :: c_int, c_double, c_ptr
         integer(c_int), value :: nshell, maxnp, natm
         integer(c_int), intent(in) :: l(*), nprim(*)
         real(c_double), intent(in) :: exps(*), coefs(*), centres(*)
         real(c_double), intent(in) :: zatm(*), ratm(*)
         type(c_ptr), intent(out) :: handle
      end function trc_basis_create

      integer(c_int) function trc_basis_create_libcint(atm, natm, bas, nshell, &
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
            env, nenv, cartesian, handle) bind(c)
         !! Build a basis from libcint's packed `atm`/`bas`/`env`.
         !!
         !! CARTESIAN ONLY -- pass `cartesian = 1`. A spherical basis is refused
         !! with TRC_ERR_UNSUPPORTED rather than reinterpreted, because
         !! reinterpreting it gives wrong numbers instead of an error.
         import :: c_int, c_double, c_ptr
         integer(c_int), value :: natm, nshell, nenv, cartesian
         integer(c_int), intent(in) :: atm(*), bas(*)
         real(c_double), intent(in) :: env(*)
         type(c_ptr), intent(out) :: handle
      end function trc_basis_create_libcint

      integer(c_int) function trc_basis_nao(handle, nao) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nao
      end function trc_basis_nao

      integer(c_int) function trc_basis_maxl(handle, maxl) bind(c)
         !! Highest angular momentum present, so a caller can decide to fall back
         !! before building anything expensive. The four-centre kernels stop at d.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: maxl
      end function trc_basis_maxl

      integer(c_int) function trc_basis_destroy(handle) bind(c)
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_basis_destroy

      ! --- shell pairs ------------------------------------------------------

      integer(c_int) function trc_pairs_create(basis, thresh, handle) bind(c)
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: basis
         real(c_double), value :: thresh
         type(c_ptr), intent(out) :: handle
      end function trc_pairs_create

      integer(c_int) function trc_pairs_count(handle, npair) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: npair
      end function trc_pairs_count

      integer(c_int) function trc_pairs_destroy(handle) bind(c)
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_pairs_destroy

      ! --- one electron -----------------------------------------------------

      integer(c_int) function trc_compute_1e(basis, pairs, smat, tmat, vmat) &
            bind(c)
         !! Overlap, kinetic and nuclear attraction, `nao*nao` each.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: basis, pairs
         real(c_double), intent(out) :: smat(*), tmat(*), vmat(*)
      end function trc_compute_1e

      integer(c_int) function trc_multipole_count() bind(c)
         !! Number of Cartesian multipole components, so a caller sizes its
         !! buffer from the library rather than hardcoding it.
         import :: c_int
      end function trc_multipole_count

      integer(c_int) function trc_compute_multipoles(basis, pairs, origin, &
            mmat) bind(c)
         !! Multipoles through the octupole about `origin`, `nao*nao*count`,
         !! component index slowest.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: basis, pairs
         real(c_double), intent(in) :: origin(*)
         real(c_double), intent(out) :: mmat(*)
      end function trc_compute_multipoles

      ! --- four-centre ------------------------------------------------------

      integer(c_int) function trc_eri_create(basis, thresh, handle) bind(c)
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
         !! The per-geometry ERI object: Schwarz bounds, binned pair list and
         !! primitive-pair tables. Expensive to build and meant to be reused by
         !! every Fock build; one per SCF iteration defeats the purpose.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: basis
         real(c_double), value :: thresh
         type(c_ptr), intent(out) :: handle
      end function trc_eri_create

      integer(c_int) function trc_eri_destroy(handle) bind(c)
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_eri_destroy

      integer(c_int) function trc_fock(eri, basis, dmat, gmat, jfac, kfac, &
            dscreen) bind(c)
         !! `g = jfac*J - kfac*K/2` for one density.
         !!
         !! `dscreen` nonzero weights the Schwarz bound by the density, which is
         !! what an SCF wants. A COUPLED-PERTURBED SOLVE MUST PASS ZERO: its
         !! density is a trial vector driven towards zero, so a density-keyed
         !! screen tightens as the solve proceeds and the operator stops being
         !! the same linear map between matvecs.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: eri, basis
         real(c_double), intent(in) :: dmat(*)
         real(c_double), intent(out) :: gmat(*)
         real(c_double), value :: jfac, kfac
         integer(c_int), value :: dscreen
      end function trc_fock

      integer(c_int) function trc_scf(basis, nalpha, nbeta, functional, grid_level, &
            conv_energy, conv_density, max_iter, dguess, verbose, energy, e_xc, dmat, eps, &
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
            niter) bind(c)
         !! The whole SCF, HF or Kohn-Sham, R or U. `functional` NUL-terminated
         !! (empty for HF); `dguess` NULL or a (nao,nao,nspin) density; `verbose`
         !! non-zero prints one line per iteration; `dmat` (nao,nao,nspin) and
         !! `eps` (nao,nspin) out, nspin = 2 iff nalpha /= nbeta.
         !! TRC_ERR_NOCONV means the last iterate is in the outputs.
         import :: c_int, c_double, c_ptr, c_char
         type(c_ptr), value :: basis
         integer(c_int), value :: nalpha, nbeta, grid_level, max_iter, verbose
         character(kind=c_char), intent(in) :: functional(*)
         real(c_double), value :: conv_energy, conv_density
         type(c_ptr), value :: dguess
         real(c_double), intent(out) :: energy, e_xc
         real(c_double), intent(out) :: dmat(*), eps(*)
         integer(c_int), intent(out) :: niter
      end function trc_scf

      integer(c_int) function trc_scf_mpi(fcomm, basis, nalpha, nbeta, functional, grid_level, &
            conv_energy, conv_density, max_iter, dguess, verbose, energy, e_xc, dmat, eps, &
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
            niter) bind(c)
         !! trc_scf run collectively by every rank of the communicator whose
         !! Fortran handle is `fcomm` (MPI_Comm_c2f from C); only the world
         !! communicator is accepted, anything else is TRC_ERR_UNSUPPORTED, and
         !! -1 runs on one rank as trc_scf does. Every rank binds its own device
         !! and returns the identical result. Without MPI in the build any
         !! handle but -1 means the single rank there is.
         import :: c_int, c_double, c_ptr, c_char
         integer(c_int), value :: fcomm
         type(c_ptr), value :: basis
         integer(c_int), value :: nalpha, nbeta, grid_level, max_iter, verbose
         character(kind=c_char), intent(in) :: functional(*)
         real(c_double), value :: conv_energy, conv_density
         type(c_ptr), value :: dguess
         real(c_double), intent(out) :: energy, e_xc
         real(c_double), intent(out) :: dmat(*), eps(*)
         integer(c_int), intent(out) :: niter
      end function trc_scf_mpi

      integer(c_int) function trc_rimp2(fcomm, basis, aux, nfrozen, aux_block, e_os, e_ss) bind(c)
         !! DEPRECATED: the context surface below (trc_create ...) replaces this;
         !! kept one release for callers already written to it.
         !! RI-MP2 on the orbitals of the last trc_scf/trc_scf_mpi on `basis`
         !! (restricted only), with `aux` the auxiliary basis, `nfrozen` frozen
         !! doubly occupied orbitals, `aux_block` the depth of the three-index
         !! tensor held at once (0: all). E_os and E_ss separately; MP2 is the
         !! sum. `fcomm` -1 runs on one rank; the world communicator's Fortran
         !! handle splits the occupied orbitals over its ranks.
         import :: c_int, c_double, c_ptr
         integer(c_int), value :: fcomm
         type(c_ptr), value :: basis, aux
         integer(c_int), value :: nfrozen, aux_block
         real(c_double), intent(out) :: e_os, e_ss
      end function trc_rimp2

      integer(c_int) function trc_fock_many(eri, basis, ndens, dmats, gmats, &
            jfac, kfac, dscreen) bind(c)
         !! N densities against one pass over the integrals.
         !! `(ndens, nao, nao)`, density index FIRST -- see the header.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: eri, basis
         integer(c_int), value :: ndens, dscreen
         real(c_double), intent(in) :: dmats(*)
         real(c_double), intent(out) :: gmats(*)
         real(c_double), value :: jfac, kfac
      end function trc_fock_many

      integer(c_int) function trc_fock_nosym(eri, basis, dmat, gmat, dscreen) &
            bind(c)
         !! The same build for a density that is NOT symmetric.
         !!
         !! A separate name rather than a flag on `trc_fock`: the folded
         !! eight-fold digestion assumes D(mu,nu) = D(nu,mu), and feeding it an
         !! antisymmetric response density gives a plausible wrong answer with no
         !! diagnostic. Two names cannot be confused by accident.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: eri, basis
         real(c_double), intent(in) :: dmat(*)
         real(c_double), intent(out) :: gmat(*)
         integer(c_int), value :: dscreen
      end function trc_fock_nosym

      ! --- density fitting --------------------------------------------------

      integer(c_int) function trc_compute_df2c(aux, jmat) bind(c)
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: aux
         real(c_double), intent(out) :: jmat(*)
      end function trc_compute_df2c

      integer(c_int) function trc_compute_df3c(basis, pairs, aux, tens) bind(c)
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: basis, pairs, aux
         real(c_double), intent(out) :: tens(*)
      end function trc_compute_df3c


      ! ==================================================================
      ! THE CONTEXT. One handle for a whole calculation: create, set up in
      ! any order, run, read back, destroy. The three-handle entries above
      ! stay for one release and are deprecated; new callers use these.
      ! Statuses: 0 TRC_OK, 1 NULL, 2 BADARG, 3 UNSUPPORTED, 4 NOCONV,
      ! 5 STATE (called before what it needs was set).
      ! ==================================================================

      integer(c_int) function trc_create(handle) bind(c)
         !! A context with defaults in every setting: HF, SAD guess, one rank.
         import :: c_int, c_ptr
         type(c_ptr), intent(out) :: handle
      end function trc_create

      integer(c_int) function trc_destroy(handle) bind(c)
         !! Release everything the context holds, device data included.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_destroy

      integer(c_int) function trc_set_molecule(handle, natm, z, xyz, charge, multiplicity) bind(c)
         !! Nuclear charges as doubles, coordinates (3, natm) in Bohr.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: natm, charge, multiplicity
         real(c_double), intent(in) :: z(*), xyz(*)
      end function trc_set_molecule

      integer(c_int) function trc_set_basis_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) bind(c)
         !! The orbital basis from libcint's packed arrays. Cartesian only.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: natm, nshell, nenv, cartesian
         integer(c_int), intent(in) :: atm(*), bas(*)
         real(c_double), intent(in) :: env(*)
      end function trc_set_basis_libcint

      integer(c_int) function trc_set_aux_libcint(handle, atm, natm, bas, nshell, env, nenv, cartesian) bind(c)
         !! The auxiliary basis for RI-MP2, same layout.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: natm, nshell, nenv, cartesian
         integer(c_int), intent(in) :: atm(*), bas(*)
         real(c_double), intent(in) :: env(*)
      end function trc_set_aux_libcint

      integer(c_int) function trc_set_basis_arrays(handle, nshell, maxnp, l, nprim, exps, coefs, centres) bind(c)
         !! From flat shell arrays, `exps`/`coefs` (maxnp, nshell) column-major,
         !! `centres` (3, nshell). Needs the molecule set first.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: nshell, maxnp
         integer(c_int), intent(in) :: l(*), nprim(*)
         real(c_double), intent(in) :: exps(*), coefs(*), centres(*)
      end function trc_set_basis_arrays

      integer(c_int) function trc_set_aux_arrays(handle, nshell, maxnp, l, nprim, exps, coefs, centres) bind(c)
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: nshell, maxnp
         integer(c_int), intent(in) :: l(*), nprim(*)
         real(c_double), intent(in) :: exps(*), coefs(*), centres(*)
      end function trc_set_aux_arrays

      integer(c_int) function trc_set_basis_json(handle, path) bind(c)
         !! From a MolSSI BSE JSON file, NUL-terminated path. Needs the molecule first.
         import :: c_int, c_ptr, c_char
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_set_basis_json

      integer(c_int) function trc_set_aux_json(handle, path) bind(c)
         import :: c_int, c_ptr, c_char
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: path(*)
      end function trc_set_aux_json

      integer(c_int) function trc_set_method(handle, functional, grid_level) bind(c)
         !! NUL-terminated functional name; empty is Hartree-Fock.
         import :: c_int, c_ptr, c_char
         type(c_ptr), value :: handle
         character(kind=c_char), intent(in) :: functional(*)
         integer(c_int), value :: grid_level
      end function trc_set_method

      integer(c_int) function trc_set_convergence(handle, conv_energy, conv_diis, max_iter, ndiis) bind(c)
         !! The SCF stops when the energy change is under conv_energy AND the
         !! commutator norm |X^T(FDS-SDF)X| is under conv_diis.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), value :: conv_energy, conv_diis
         integer(c_int), value :: max_iter, ndiis
      end function trc_set_convergence

      integer(c_int) function trc_set_screening(handle, thresh) bind(c)
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), value :: thresh
      end function trc_set_screening

      integer(c_int) function trc_set_guess(handle, kind, dguess, nspin) bind(c)
         !! `kind`: 0 core, 1 GWH, 2 SAD (the default; built here, one atomic
         !! SCF per element), 3 a density given in `dguess`, (nao, nao, nspin)
         !! column-major, copied at this call. Set the basis before giving one.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: kind, nspin
         type(c_ptr), value :: dguess
      end function trc_set_guess

      integer(c_int) function trc_set_comm(handle, fcomm) bind(c)
         !! A Fortran communicator handle to split the run over; -1 is one
         !! rank. MPI_COMM_WORLD only, and every rank must make the same calls.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: fcomm
      end function trc_set_comm

      integer(c_int) function trc_set_verbose(handle, level) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: level
      end function trc_set_verbose

      integer(c_int) function trc_set_rimp2(handle, nfrozen, aux_block) bind(c)
         !! `nfrozen` -1 counts the core from the elements; `aux_block` <= 0
         !! takes the auxiliary basis in one block.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), value :: nfrozen, aux_block
      end function trc_set_rimp2

      integer(c_int) function trc_run_scf(handle) bind(c)
         !! TRC_ERR_NOCONV when the iterations ran out; the last iteration's
         !! numbers are still readable.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_run_scf

      integer(c_int) function trc_run_rimp2(handle) bind(c)
         !! After a converged closed-shell SCF, with an auxiliary basis set.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_run_rimp2

      integer(c_int) function trc_nao(handle, nao) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nao
      end function trc_nao

      integer(c_int) function trc_nspin(handle, nspin) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: nspin
      end function trc_nspin

      integer(c_int) function trc_energy(handle, energy) bind(c)
         !! The SCF total, plus the RI-MP2 correlation once that has run.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: energy
      end function trc_energy

      integer(c_int) function trc_energy_parts(handle, e_nuc, e_one, e_two, e_xc) bind(c)
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: e_nuc, e_one, e_two, e_xc
      end function trc_energy_parts

      integer(c_int) function trc_converged(handle, flag) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: flag
      end function trc_converged

      integer(c_int) function trc_iterations(handle, niter) bind(c)
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: niter
      end function trc_iterations

      integer(c_int) function trc_density(handle, dmat) bind(c)
         !! (nao, nao, nspin) column-major; restricted holds the total.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: dmat(*)
      end function trc_density

      integer(c_int) function trc_mo_coefficients(handle, cmo) bind(c)
         !! (nao, nao, nspin) column-major, columns are orbitals.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: cmo(*)
      end function trc_mo_coefficients

      integer(c_int) function trc_mo_energies(handle, eps) bind(c)
         !! (nao, nspin).
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: eps(*)
      end function trc_mo_energies

      integer(c_int) function trc_rimp2_energy(handle, e_os, e_ss) bind(c)
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: e_os, e_ss
      end function trc_rimp2_energy

      integer(c_int) function trc_message(handle, buf, buflen) bind(c)
         !! The last message, NUL-terminated into `buf`.
         import :: c_int, c_ptr, c_char
         type(c_ptr), value :: handle
         character(kind=c_char), intent(out) :: buf(*)
         integer(c_int), value :: buflen
      end function trc_message

      integer(c_int) function trc_context_basis(handle, basis) bind(c)
         !! The basis handle inside the context, for trc_compute_* and the
         !! Fock family. BORROWED: never pass it to trc_basis_destroy.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         type(c_ptr), intent(out) :: basis
      end function trc_context_basis

      integer(c_int) function trc_context_eri(handle, eri) bind(c)
         !! The ERI context at the context's screening threshold, built on
         !! first use. BORROWED: never pass it to trc_eri_destroy.
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
         type(c_ptr), intent(out) :: eri
      end function trc_context_eri

   end interface

end module trc_c_interfaces
