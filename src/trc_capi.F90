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
module trc_capi
   use, intrinsic :: iso_c_binding
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, &
                        trc_1e, trc_df_2c, trc_df_3c, trc_multipoles, &
                        TRC_NMULT
   use trc_eri, only: trc_eri_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build
   use trc_error, only: error_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use pic_mpi_lib, only: comm_t, comm_world
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t
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
             trc_set_basis_json, trc_set_aux_libcint, trc_set_aux_arrays, trc_set_aux_json, trc_set_method, &
             trc_set_convergence, trc_set_screening, trc_set_guess, trc_set_comm, trc_set_verbose, &
             trc_set_rimp2, trc_run_scf, trc_run_rimp2, trc_nao, trc_nspin, trc_energy, trc_energy_parts, &
             trc_converged, trc_iterations, trc_density, trc_mo_coefficients, trc_mo_energies, &
             trc_rimp2_energy, trc_message, trc_context_basis, trc_context_eri
   public :: capi_scf

   ! Status codes. Zero is success; everything else is a reason.
   integer(c_int), parameter, public :: TRC_OK            = 0
   integer(c_int), parameter, public :: TRC_ERR_NULL      = 1
   integer(c_int), parameter, public :: TRC_ERR_BADARG    = 2
   integer(c_int), parameter, public :: TRC_ERR_UNSUPPORTED = 3
   integer(c_int), parameter, public :: TRC_ERR_NOCONV = 4   !! the SCF ran out of iterations
   integer(c_int), parameter, public :: TRC_ERR_STATE = 5    !! called before what it needs was set
   ! Guess kinds for trc_set_guess.
   integer(c_int), parameter, public :: TRC_GUESS_CORE = 0, TRC_GUESS_GWH = 1, TRC_GUESS_SAD = 2, &
                                        TRC_GUESS_GIVEN = 3

   ! Wrappers so a bare derived type can be pointed at from C.
   type :: basis_box
      type(trc_basis_t) :: b
      ! The orbitals of the last trc_scf on this basis, kept for a
      ! correlated step that follows: trc_scf hands the caller a density and
      ! eigenvalues, not orbitals, and its signature is not changing under
      ! callers already written to it.
      real(dp), allocatable :: cmo(:, :, :), eps(:, :)
      integer :: nalpha = -1, nbeta = -1
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
      ! the ERI context, only if something borrows it
      type(eri_box), allocatable :: eb
      ! results
      type(trc_scf_result_t) :: res
      integer :: nalpha = 0, nbeta = 0
      logical :: scf_ok = .false., rimp2_ok = .false.
      real(dp) :: e_os = 0.0_dp, e_ss = 0.0_dp
      character(len=256) :: message = ""
   end type context_box

contains

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

   !> Number of atomic orbitals, so a caller can size its own output buffers.
   function trc_basis_nao(handle, nao) result(status) &
      bind(c, name="trc_basis_nao")
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
      type(c_ptr), value :: basis, pairs
      real(c_double), intent(out) :: smat(*), tmat(*), vmat(*)
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(pairs_box), pointer :: pb
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :)
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
      tens(1:int(n, 8)*n*m) = reshape(t, [n*n*m])
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

   !> Highest angular momentum in the basis, so a caller can decide whether to
   !> use this library at all before it builds anything expensive. The
   !> four-centre kernels stop at d and the caller has to fall back rather than
   !> get a wrong answer.
   function trc_basis_maxl(handle, maxl) result(status) &
      bind(c, name="trc_basis_maxl")
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

   !> Number of Cartesian multipole components, so a caller can size its buffer
   !> without hardcoding 39 and drifting when the operator order changes.
   !> Named for what it counts rather than `trc_nmult`, which collides with the
   !> parameter of that name it returns.
   function trc_multipole_count() result(n) &
      bind(c, name="trc_multipole_count")
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
      type(c_ptr), value :: basis, pairs
      real(c_double), intent(in)  :: origin(3)
      real(c_double), intent(out) :: mmat(*)
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(pairs_box), pointer :: pb
      real(dp), allocatable :: m(:, :, :)
      integer :: n
      integer(kind=8) :: ntot

      status = TRC_ERR_NULL
      if (.not. c_associated(basis) .or. .not. c_associated(pairs)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(pairs, pb)
      n = bb%b%nao
      ntot = int(n, 8)*int(n, 8)*int(TRC_NMULT, 8)
      allocate (m(n, n, TRC_NMULT))
      !$acc enter data create(m)
      call trc_multipoles(bb%b, pb%p, real(origin, dp), m)
      !$acc update self(m)
      !$acc exit data delete(m)
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
   function capi_fock(eri, basis, dmat, gmat, jfac, kfac, dscreen) &
      result(status) bind(c, name="trc_fock")
      type(c_ptr), value :: eri, basis
      real(c_double), intent(in)  :: dmat(*)
      real(c_double), intent(out) :: gmat(*)
      real(c_double), value :: jfac, kfac
      !> Nonzero to weight the Schwarz bound by the density (an SCF), zero for
      !> the plain Schwarz screen (a coupled-perturbed solve, whose trial
      !> vector the solver drives to zero -- a density-keyed screen tightens
      !> as it does and the operator stops being a fixed linear map).
      !> An int rather than a logical: Fortran's logical has no guaranteed C
      !> representation and c_bool is one byte against C's four-byte int.
      integer(c_int), value :: dscreen
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      integer :: n

      status = TRC_ERR_NULL
      if (.not. c_associated(eri) .or. .not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(eri, eb)
      n = bb%b%nao
      call fock_shaped(bb%b, eb%e, n, dmat, gmat, real(jfac, dp), &
                       real(kfac, dp), dscreen /= 0)
      status = TRC_OK
   end function capi_fock

   !
   ! NO COPY. `dmat` and `gmat` arrive as assumed-size, and handing an
   ! assumed-size actual to an explicit-shape `(n, n)` dummy is sequence
   ! association -- the callee addresses the caller's memory directly.
   !
   ! The first version reshaped both into allocatable locals, which at 2773
   ! functions is four passes over 61 MB per Fock build for nothing. Reshaping
   ! is the obvious way to write this and the obvious way is wrong: the bytes
   ! are already contiguous and already in the right order, so all the copy
   ! bought was a second name for them.
   !
   subroutine fock_shaped(b, e, n, d, g, jfac, kfac, dscreen)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: e
      integer, intent(in) :: n
      real(dp), intent(in)  :: d(n, n)
      real(dp), intent(out) :: g(n, n)
      real(dp), intent(in) :: jfac, kfac
      logical, intent(in) :: dscreen
      call e%fock(b, d, g, k_scale=kfac, j_scale=jfac, density_screen=dscreen)
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
   function trc_fock_many(eri, basis, ndens, dmats, gmats, jfac, kfac, dscreen) &
      result(status) bind(c, name="trc_fock_many")
      type(c_ptr), value :: eri, basis
      integer(c_int), value :: ndens
      real(c_double), intent(in)  :: dmats(*)
      real(c_double), intent(out) :: gmats(*)
      real(c_double), value :: jfac, kfac
      integer(c_int), value :: dscreen
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      integer :: n

      status = TRC_ERR_NULL
      if (.not. c_associated(eri) .or. .not. c_associated(basis)) return
      if (ndens < 1) then
         status = TRC_ERR_BADARG
         return
      end if
      call c_f_pointer(basis, bb)
      call c_f_pointer(eri, eb)
      n = bb%b%nao
      call fock_many_shaped(bb%b, eb%e, int(ndens), n, dmats, gmats, &
                            real(jfac, dp), real(kfac, dp), dscreen /= 0)
      status = TRC_OK
   end function trc_fock_many

   !> As `fock_shaped`, and for the same reason -- the batch is the case where
   !> the copy would have hurt most, being N times the size.
   subroutine fock_many_shaped(b, e, nd, n, d, g, jfac, kfac, dscreen)
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: e
      integer, intent(in) :: nd, n
      real(dp), intent(in)  :: d(nd, n, n)
      real(dp), intent(out) :: g(nd, n, n)
      real(dp), intent(in) :: jfac, kfac
      logical, intent(in) :: dscreen
      call e%fock_many(b, nd, d, g, k_scale=kfac, j_scale=jfac, &
                       density_screen=dscreen)
   end subroutine fock_many_shaped

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
      type(c_ptr), value :: eri, basis
      real(c_double), intent(in)  :: dmat(*)
      real(c_double), intent(out) :: gmat(*)
      integer(c_int), value :: dscreen
      integer(c_int) :: status
      type(basis_box), pointer :: bb
      type(eri_box), pointer :: eb
      integer :: n

      status = TRC_ERR_NULL
      if (.not. c_associated(eri) .or. .not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(eri, eb)
      n = bb%b%nao
      call fock_nosym_shaped(bb%b, eb%e, n, dmat, gmat, dscreen /= 0)
      status = TRC_OK
   end function trc_fock_nosym

   !> As `fock_shaped`.
   subroutine fock_nosym_shaped(b, e, n, d, g, dscreen)
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
   ! together. Each rank binds the device numbered by its rank, modulo the
   ! devices present, before anything else -- the caller does not.
   !
   ! Only the world communicator is accepted for now: pic-mpi builds its
   ! communicator object by duplicating MPI_COMM_WORLD, and terco has no
   ! way to wrap an arbitrary handle. Any other handle is
   ! TRC_ERR_UNSUPPORTED, except -1, which is trc_scf on one rank in any
   ! build. In a build without MPI every other handle is the single rank.
   !
   function capi_scf_mpi(fcomm, basis, nalpha, nbeta, functional, grid_level, conv_energy, &
                         conv_density, max_iter, dguess, verbose, energy, e_xc, dmat, eps, &
                         niter) result(status) bind(c, name="trc_scf_mpi")
#ifdef TERCO_HAVE_MPI
      use mpi_f08, only: MPI_COMM_WORLD
#endif
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
      integer :: dev
      if (fcomm == -1) then
         ! One rank, no communicator: trc_scf.
         status = scf_entry(basis, nalpha, nbeta, functional, grid_level, conv_energy, conv_density, &
                            max_iter, dguess, verbose, energy, e_xc, dmat, eps, niter)
         return
      end if
#ifdef TERCO_HAVE_MPI
      if (fcomm /= MPI_COMM_WORLD%mpi_val) then
         status = TRC_ERR_UNSUPPORTED
         energy = 0.0_c_double; e_xc = 0.0_c_double; niter = 0
         return
      end if
#endif
      comm = comm_world()
      dev = trc_bind_device(comm%rank())
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
   ! weights on them. With `fcomm` the world communicator's Fortran handle
   ! the occupied orbitals are split over its ranks and every rank gets the
   ! same numbers; pass any value in a build without MPI, or run on one
   ! rank with fcomm = -1 in any build.
   !
   function capi_rimp2(fcomm, basis, aux, nfrozen, aux_block, e_os, e_ss) &
      result(status) bind(c, name="trc_rimp2")
#ifdef TERCO_HAVE_MPI
      use mpi_f08, only: MPI_COMM_WORLD
#endif
      integer(c_int), value :: fcomm
      type(c_ptr), value :: basis, aux
      integer(c_int), value :: nfrozen, aux_block
      real(c_double), intent(out) :: e_os, e_ss
      integer(c_int) :: status
      type(basis_box), pointer :: bb, ab
      type(trc_pairlist_t) :: pl
      type(trc_rimp2_result_t) :: res
      type(comm_t) :: comm
      integer :: dev, nblk
      logical :: collective

      e_os = 0.0_c_double; e_ss = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(basis) .or. .not. c_associated(aux)) return
      call c_f_pointer(basis, bb)
      call c_f_pointer(aux, ab)
      status = TRC_ERR_BADARG
      if (.not. allocated(bb%cmo)) return           ! no SCF on this basis yet
      if (bb%nalpha /= bb%nbeta) return             ! restricted reference only
      if (nfrozen < 0 .or. nfrozen >= bb%nalpha) return
      collective = fcomm /= -1
#ifdef TERCO_HAVE_MPI
      if (collective .and. fcomm /= MPI_COMM_WORLD%mpi_val) then
         status = TRC_ERR_UNSUPPORTED
         return
      end if
#endif
      nblk = int(aux_block)
      if (nblk <= 0) nblk = ab%b%nao
      if (.not. ab%b%on_device) call ab%b%to_device()
      call pl%build(bb%b, 1.0e-12_dp)
      call pl%to_device()
      if (collective) then
         comm = comm_world()
         dev = trc_bind_device(comm%rank())
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
      integer :: n, nspin, i, s, k

      status = TRC_ERR_NULL
      energy = 0.0_c_double; e_xc = 0.0_c_double; niter = 0
      if (.not. c_associated(basis)) return
      call c_f_pointer(basis, bb)
      n = bb%b%nao
      status = TRC_ERR_BADARG
      if (nalpha < 0 .or. nbeta < 0 .or. nalpha + nbeta > 2*n) return
      if (grid_level < 0 .or. max_iter < 1) return
      nspin = merge(2, 1, nalpha /= nbeta)

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
         call trc_scf_run(bb%b, int(nalpha), int(nbeta), opts, res, dguess=real(dg, dp), comm=comm)
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
      k = 0
      do s = 1, nspin
         do i = 1, n*n
            k = k + 1
            dmat(k) = real(res%dmat(mod(i - 1, n) + 1, (i - 1)/n + 1, s), c_double)
         end do
      end do
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
      bb%nalpha = int(nalpha); bb%nbeta = int(nbeta)
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
      write (t, '(i0)') i
   end function itoa

   subroutine invalidate_runs(cx)
      !! A setting the results depended on changed.
      type(context_box), intent(inout) :: cx
      cx%scf_ok = .false.; cx%rimp2_ok = .false.
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
      cx%natm = int(natm); cx%charge = int(charge); cx%multiplicity = int(multiplicity)
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
         cx%have_aux = .true.
         cx%rimp2_ok = .false.
         call b%release()
         return
      else
         if (cx%have_basis) call cx%bb%b%release()
         cx%bb%b = b
         call cx%bb%b%to_device()
         if (allocated(cx%bb%cmo)) deallocate (cx%bb%cmo, cx%bb%eps)
         cx%have_basis = .true.
         call drop_eri(cx)
      end if
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
      call trc_basis_from_json(trim(fpath), cx%natm, nint(cx%at_z), cx%at_r, b, err)
      if (err%has_error()) then
         cx%message = err%get_message()
         status = TRC_ERR_BADARG
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
                        "'; terco has lda_x (slater), svwn, svwn_rpa, pbe, blyp, pbe0, b3lyp, b3lyp5, "// &
                        "or an empty name for Hartree-Fock")
            cx%opts%functional = ""
            return
         end if
      end if
      cx%opts%grid_level = int(grid_level)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_method

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
         call refuse(cx, "trc_set_convergence: conv_energy must be positive"); return
      end if
      if (conv_diis <= 0.0_c_double) then
         call refuse(cx, "trc_set_convergence: conv_diis must be positive (the commutator norm gate)"); return
      end if
      if (max_iter < 1) then
         call refuse(cx, "trc_set_convergence: max_iter must be at least 1, got "//trim(itoa(int(max_iter)))); return
      end if
      if (ndiis < 1) then
         call refuse(cx, "trc_set_convergence: ndiis must be at least 1, got "//trim(itoa(int(ndiis)))); return
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

   function trc_set_guess(handle, kind, dguess, nspin) result(status) bind(c, name="trc_set_guess")
      !! `kind` one of TRC_GUESS_*. For TRC_GUESS_GIVEN, `dguess` is
      !! (nao, nao, nspin) column-major and is COPIED here; the caller's array
      !! is free afterwards. The basis must be set first for that copy.
      type(c_ptr), value :: handle
      integer(c_int), value :: kind, nspin
      type(c_ptr), value :: dguess
      integer(c_int) :: status
      type(context_box), pointer :: cx
      real(c_double), pointer :: dg(:, :, :)
      integer :: n
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_BADARG
      if (kind < TRC_GUESS_CORE .or. kind > TRC_GUESS_GIVEN) then
         call refuse(cx, "trc_set_guess: kind "//trim(itoa(int(kind)))//" -- 0 core, 1 GWH, 2 SAD, 3 given")
         return
      end if
      if (allocated(cx%dguess)) deallocate (cx%dguess)
      if (kind == TRC_GUESS_GIVEN) then
         if (.not. c_associated(dguess)) then
            call refuse(cx, "trc_set_guess: kind 3 (given) needs a density pointer"); return
         end if
         if (nspin < 1 .or. nspin > 2) then
            call refuse(cx, "trc_set_guess: nspin must be 1 or 2, got "//trim(itoa(int(nspin)))); return
         end if
         status = TRC_ERR_STATE
         if (.not. cx%have_basis) then
            call refuse(cx, "trc_set_guess: set the basis before giving a density; its size comes from there")
            return
         end if
         n = cx%bb%b%nao
         call c_f_pointer(dguess, dg, [n, n, int(nspin)])
         allocate (cx%dguess(n, n, nspin))
         cx%dguess = real(dg, dp)
      end if
      cx%guess_kind = int(kind)
      call invalidate_runs(cx)
      status = TRC_OK
   end function trc_set_guess

   function trc_set_comm(handle, fcomm) result(status) bind(c, name="trc_set_comm")
      !! A Fortran communicator handle; -1 means one rank. MPI_COMM_WORLD only.
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
      cx%nfrozen = int(nfrozen); cx%aux_block = int(aux_block)
      cx%rimp2_ok = .false.
      status = TRC_OK
   end function trc_set_rimp2

   subroutine world_comm(cx, comm, ok, status)
      !! The communicator a run splits over, if the context names one.
#ifdef TERCO_HAVE_MPI
      use mpi_f08, only: MPI_COMM_WORLD
#endif
      type(context_box), intent(in) :: cx
      type(comm_t), intent(out) :: comm
      logical, intent(out) :: ok
      integer(c_int), intent(out) :: status
      integer :: dev
      ok = .false.; status = TRC_OK
      if (cx%fcomm == -1) return
#ifdef TERCO_HAVE_MPI
      if (cx%fcomm /= MPI_COMM_WORLD%mpi_val) then
         status = TRC_ERR_UNSUPPORTED
         return
      end if
#endif
      comm = comm_world()
      dev = trc_bind_device(comm%rank())
      ok = .true.
   end subroutine world_comm

   function trc_run_scf(handle) result(status) bind(c, name="trc_run_scf")
      !! Needs the molecule and the basis. Builds the guess it was asked for,
      !! runs, keeps everything in the context.
      type(c_ptr), value :: handle
      integer(c_int) :: status
      type(context_box), pointer :: cx
      type(trc_scf_options_t) :: opts
      type(comm_t) :: comm
      type(error_t) :: err
      real(dp), allocatable :: dsad(:, :), dg(:, :, :)
      integer :: nelec, nunp, nalpha, nbeta, nspin, n
      logical :: collective
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      cx%scf_ok = .false.; cx%rimp2_ok = .false.
      status = TRC_ERR_STATE
      if (.not. (cx%have_mol .and. cx%have_basis)) then
         cx%message = "trc_run_scf: set the molecule and the basis first"
         return
      end if
      n = cx%bb%b%nao
      nelec = nint(sum(cx%at_z)) - cx%charge
      nunp = cx%multiplicity - 1
      status = TRC_ERR_BADARG
      if (nelec < 0 .or. nunp > nelec .or. mod(nelec - nunp, 2) /= 0) then
         cx%message = "trc_run_scf: charge and multiplicity do not fit the electron count"
         return
      end if
      nalpha = (nelec + nunp)/2; nbeta = nelec - nalpha
      if (nalpha > n) then
         cx%message = "trc_run_scf: more electrons than functions"
         return
      end if
      nspin = merge(2, 1, nalpha /= nbeta .or. cx%opts%unrestricted)
      opts = cx%opts
      select case (cx%guess_kind)
      case (TRC_GUESS_CORE)
         opts%guess = "core"
      case (TRC_GUESS_GWH)
         opts%guess = "gwh"
      case (TRC_GUESS_SAD)
         call trc_sad_build(cx%bb%b, dsad, err, verbose=cx%opts%verbose)
         if (err%has_error()) then
            cx%message = err%get_message()
            return
         end if
         allocate (dg(n, n, nspin))
         if (nspin == 1) then
            dg(:, :, 1) = dsad
         else
            dg(:, :, 1) = 0.5_dp*dsad; dg(:, :, 2) = 0.5_dp*dsad
         end if
      case (TRC_GUESS_GIVEN)
         if (.not. allocated(cx%dguess)) then
            cx%message = "trc_run_scf: guess GIVEN but no density was given"
            status = TRC_ERR_STATE
            return
         end if
         if (size(cx%dguess, 3) /= nspin) then
            cx%message = "trc_run_scf: the given density has the wrong spin count for this molecule"
            return
         end if
         dg = cx%dguess
      end select
      call world_comm(cx, comm, collective, status)
      if (status /= TRC_OK) return
      if (allocated(dg)) then
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, dguess=dg, comm=comm)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, dguess=dg)
         end if
      else
         if (collective) then
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res, comm=comm)
         else
            call trc_scf_run(cx%bb%b, nalpha, nbeta, opts, cx%res)
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
      cx%bb%cmo = cx%res%cmo; cx%bb%eps = cx%res%eps
      cx%bb%nalpha = nalpha; cx%bb%nbeta = nbeta
      cx%nalpha = nalpha; cx%nbeta = nbeta
      cx%scf_ok = .true.
      status = merge(TRC_OK, TRC_ERR_NOCONV, cx%res%converged)
   end function trc_run_scf

   pure integer function core_orbitals(z)
      !! Frozen-core orbitals for one element, by period.
      real(dp), intent(in) :: z
      integer :: iz
      iz = nint(z)
      if (iz <= 2) then
         core_orbitals = 0
      else if (iz <= 10) then
         core_orbitals = 1
      else if (iz <= 18) then
         core_orbitals = 5
      else if (iz <= 36) then
         core_orbitals = 9
      else if (iz <= 54) then
         core_orbitals = 18
      else
         core_orbitals = 27
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
      call world_comm(cx, comm, collective, status)
      if (status /= TRC_OK) return
      call pl%build(cx%bb%b, cx%opts%eri_thresh)
      call pl%to_device()
      if (collective) then
         call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cx%res%cmo(:, :, 1), cx%res%eps(:, 1), res, &
                            nfrozen=nfrozen, aux_block=nblk, comm=comm, verbose=cx%opts%verbose)
         call comm%finalize()
      else
         call trc_rimp2_run(cx%bb%b, cx%ab%b, pl, cx%nalpha, cx%res%cmo(:, :, 1), cx%res%eps(:, 1), res, &
                            nfrozen=nfrozen, aux_block=nblk, verbose=cx%opts%verbose)
      end if
      call pl%release()
      cx%e_os = res%e_os; cx%e_ss = res%e_ss
      cx%rimp2_ok = .true.
      status = TRC_OK
   end function trc_run_rimp2

   ! ---- read back ---------------------------------------------------------

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
      nao = int(cx%bb%b%nao, c_int)
      status = TRC_OK
   end function trc_nao

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
      e_nuc = 0.0_c_double; e_one = 0.0_c_double; e_two = 0.0_c_double; e_xc = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      e_nuc = real(cx%res%e_nuc, c_double); e_one = real(cx%res%e_one, c_double)
      e_two = real(cx%res%e_two, c_double); e_xc = real(cx%res%e_xc, c_double)
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
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: dmat(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: n, ns
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      n = cx%bb%b%nao; ns = cx%res%nspin
      dmat(1:n*n*ns) = real(reshape(cx%res%dmat, [n*n*ns]), c_double)
      status = TRC_OK
   end function trc_density

   function trc_mo_coefficients(handle, cmo) result(status) bind(c, name="trc_mo_coefficients")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: cmo(*)
      integer(c_int) :: status
      type(context_box), pointer :: cx
      integer :: n, ns
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%scf_ok) return
      n = cx%bb%b%nao; ns = cx%res%nspin
      cmo(1:n*n*ns) = real(reshape(cx%res%cmo, [n*n*ns]), c_double)
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
      n = cx%bb%b%nao; ns = cx%res%nspin
      eps(1:n*ns) = real(reshape(cx%res%eps, [n*ns]), c_double)
      status = TRC_OK
   end function trc_mo_energies

   function trc_rimp2_energy(handle, e_os, e_ss) result(status) bind(c, name="trc_rimp2_energy")
      type(c_ptr), value :: handle
      real(c_double), intent(out) :: e_os, e_ss
      integer(c_int) :: status
      type(context_box), pointer :: cx
      e_os = 0.0_c_double; e_ss = 0.0_c_double
      status = TRC_ERR_NULL
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, cx)
      status = TRC_ERR_STATE
      if (.not. cx%rimp2_ok) return
      e_os = real(cx%e_os, c_double); e_ss = real(cx%e_ss, c_double)
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
         call world_comm(cx, comm, collective, status)
         if (status /= TRC_OK) then
            deallocate (cx%eb)
            return
         end if
         if (collective) then
            call cx%eb%e%build(cx%bb%b, cx%opts%eri_thresh, comm)
         else
            call cx%eb%e%build(cx%bb%b, cx%opts%eri_thresh)
         end if
      end if
      ep => cx%eb
      eri = c_loc(ep)
      status = TRC_OK
   end function trc_context_eri

end module trc_capi
