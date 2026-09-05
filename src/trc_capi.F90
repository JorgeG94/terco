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
   public :: capi_scf

   ! Status codes. Zero is success; everything else is a reason.
   integer(c_int), parameter, public :: TRC_OK            = 0
   integer(c_int), parameter, public :: TRC_ERR_NULL      = 1
   integer(c_int), parameter, public :: TRC_ERR_BADARG    = 2
   integer(c_int), parameter, public :: TRC_ERR_UNSUPPORTED = 3
   integer(c_int), parameter, public :: TRC_ERR_NOCONV = 4   !! the SCF ran out of iterations

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

end module trc_capi
