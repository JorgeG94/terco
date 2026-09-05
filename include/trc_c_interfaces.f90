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

   interface

      ! --- basis ------------------------------------------------------------

      integer(c_int) function trc_basis_create(nshell, maxnp, l, nprim, exps, &
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
         import :: c_int, c_ptr
         type(c_ptr), value :: handle
      end function trc_basis_destroy

      ! --- shell pairs ------------------------------------------------------

      integer(c_int) function trc_pairs_create(basis, thresh, handle) bind(c)
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
         !! The per-geometry ERI object: Schwarz bounds, binned pair list and
         !! primitive-pair tables. Expensive to build and meant to be reused by
         !! every Fock build; one per SCF iteration defeats the purpose.
         import :: c_int, c_double, c_ptr
         type(c_ptr), value :: basis
         real(c_double), value :: thresh
         type(c_ptr), intent(out) :: handle
      end function trc_eri_create

      integer(c_int) function trc_eri_destroy(handle) bind(c)
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

   end interface

end module trc_c_interfaces
