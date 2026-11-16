module trc_tddft
   !!
   !! Excitation energies of a closed-shell RHF/RKS reference in the
   !! Tamm-Dancoff approximation (TDA; CIS for Hartree-Fock): the lowest
   !! `nroots` singlet roots of the matrix A by a block Davidson iteration, with
   !! the transition dipoles and oscillator strengths.
   !!
   !! THE EQUATION (singlet, restricted)
   !! ----------------------------------
   !! With X_ia the amplitude of the i -> a excitation (i occupied, a virtual)
   !! and T = C_occ X C_vir^T its AO transition density (NON-symmetric):
   !!
   !!     [A X]_ia = (e_a - e_i) X_ia
   !!                + [C_occ^T (2 J[T] - c_x K[T] + 2 f_xc[T]) C_vir]_ia
   !!
   !! J[T]_mn = (mn|ls) T_ls and K[T]_mn = (ml|ns) T_ls, the pyscf get_jk(hermi=0)
   !! convention; c_x is the functional's exact-exchange fraction (`func%exx`,
   !! 1 for Hartree-Fock, where f_xc is absent). The exchange of a
   !! non-symmetric density is not the exchange of its symmetric part, which is
   !! why this uses `eri%fock_nosym_many`: it returns j J - k K/2, so the call
   !! passes j_scale = 2 and k_scale = 2 c_x. f_xc[T] is
   !! `trc_xc_fxc_rks_apply_many` (it symmetrises T itself, the only part the
   !! density response sees), against a cache built once from the ground-state
   !! density. A is real symmetric and, for a stable ground state, positive
   !! definite; the eigenvalues are the excitation energies in hartree.
   !!
   !! THE SOLVER
   !! ----------
   !! Every new trial vector of an iteration goes through ONE `fock_nosym_many`
   !! call and ONE `fxc_apply_many` call, i.e. one pass over the integrals per
   !! iteration whatever the number of roots. `density_screen` is false: a
   !! density-keyed screen on trial vectors would change the operator between
   !! matvecs (see `eri%fock_nosym_many`).
   !!
   !! Guesses are the unit vectors of the lowest orbital-energy differences (a
   !! tie at the cut is kept whole). Each
   !! iteration diagonalises the subspace matrix, forms the residuals of the
   !! lowest `nroots` Ritz pairs and, for every unconverged one, the
   !! diagonally preconditioned correction -R/(d - w). The corrections are
   !! orthonormalised against the subspace and each other (classical
   !! Gram-Schmidt twice) and a correction that is numerically inside the span
   !! is dropped. When the subspace would exceed `max_space` it collapses to
   !! the current Ritz vectors. All the roots iterate together as a block, so a
   !! degenerate multiplet is converged as one: nothing here assumes the roots
   !! are non-degenerate (the amplitudes within a degenerate multiplet are an
   !! arbitrary orthonormal rotation of each other).
   !!
   !! MISSED ROOTS. Davidson follows the lowest `nroots` RITZ values, and a
   !! preconditioned correction of a Ritz vector stays in the symmetry of
   !! that vector. A root of a symmetry the guesses represent poorly has a
   !! Ritz value far above its eigenvalue, so with few guesses it can sit
   !! above the tracked window and never be iterated: measured on formamide
   !! HF (Cs, A'/A''), nroots = 2 with 2*nroots guesses returned the 1st and
   !! 3rd roots and skipped the 2nd (0.3464 Eh, which has 92% of its weight on
   !! the four lowest differences, so it was not orthogonal to the guesses).
   !! The Ritz estimates of that other sector stayed at 0.389 Eh (the guesses
   !! describe it badly), above the 3rd root's 0.348, so it sat outside the
   !! tracked window and was never corrected, even with 10 guesses. Three
   !! measures: (1) a BUFFER of max(2, nroots/2) extra Ritz vectors is iterated
   !! along with the wanted ones, to a looser residual (BUFFER_TOL), so such a
   !! sector's Ritz value descends into the window and takes its place;
   !! (2) generous guesses, a tie at the cut kept whole; (3) one probe with the
   !! next-lowest differences before convergence is accepted. None is a proof:
   !! no Davidson variant guarantees the lowest roots without symmetry
   !! information, and the tests sweep nroots against pyscf for that reason.
   !!
   !! A breakdown (no usable new direction while a root is unconverged), a NaN
   !! or Inf anywhere, a non-positive excitation energy or running out of
   !! iterations is an error (`err`, `res%converged = .false.`), never a
   !! silently returned approximate answer. A NaN residual reads as not
   !! converged (`.not. (r < tol)`), as in trc_cpks.
   !!
   !! WHAT THIS DOES NOT DO
   !! ---------------------
   !! The full RPA (Casida) problem, triplets, UHF/UKS references, meta-GGA,
   !! VV10, range-separated and double hybrids: `trc_tddft_refuse_reason` names
   !! each refusal, and the solver reports it as an error before any work.
   !!
   use pic_types, only: dp, int64, default_int
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_sorting, only: sort_index
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_multipoles, TRC_NMULT
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, XC_FAMILY_LDA, XC_FAMILY_GGA
   use trc_xc, only: trc_xc_fxc_cache_t, trc_xc_fxc_rks_apply_many
   implicit none
   private

   public :: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_op_t
   public :: trc_tddft_solve_rks, trc_tddft_refuse_reason

   real(dp), parameter :: OSC_PREFACTOR = 2.0_dp/3.0_dp
   !! The 2/3 of the oscillator strength f = (2/3) w |mu|^2, atomic units.
   real(dp), parameter :: BUFFER_TOL = 1.0e-3_dp
   !! Residual the buffer Ritz vectors (those past `nroots`) are iterated to; at
   !! least 100 conv_tol. Their Ritz values are then good to ~1e-6.
   real(dp), parameter :: TIE_TOL = 1.0e-4_dp
   !! Orbital-energy differences closer than this (hartree) count as tied
   !! when the guesses are cut.
   real(dp), parameter :: PRECOND_FLOOR = 1.0e-4_dp
   !! The preconditioner's denominator (d - w) is kept at least this far from
   !! zero, sign preserved.
   real(dp), parameter :: DROP_REL = 1.0e-6_dp
   !! A correction whose norm falls below this fraction of its norm before
   !! orthogonalisation lies inside the subspace and is dropped.

   type :: trc_tddft_options_t
      integer :: nroots = 3
      !! Lowest singlet roots wanted.
      logical :: tda = .true.
      !! Only the Tamm-Dancoff approximation is implemented; false is refused.
      logical :: singlet = .true.
      !! Only singlets are implemented; false is refused.
      real(dp) :: conv_tol = 1.0e-6_dp
      !! Converged when ||A X - w X|| (Frobenius, X normalised) is below this
      !! for every root. The excitation energy error is of order
      !! conv_tol^2/gap, the amplitude error of order conv_tol/gap.
      integer :: max_space = 0
      !! Largest subspace before the collapse; 0 takes 12*nroots. Raised to
      !! 3*(nroots + buffer) if smaller, and capped at the dimension of the problem.
      integer :: max_iter = 50
      !! Most Davidson iterations (matvec passes).
      integer :: nguess = 0
      !! Unit-vector guesses; 0 takes max(2*nroots, nroots + 4), capped at the
      !! dimension and extended over a tie at the cut. A probe with the
      !! next-lowest differences runs before convergence is accepted.
      logical :: verbose = .false.
   end type trc_tddft_options_t

   type :: trc_tddft_result_t
      logical :: converged = .false.
      integer :: nroots = 0
      integer :: iterations = 0
      !! Davidson iterations taken.
      integer :: nmatvec = 0
      !! Trial vectors multiplied by A in total.
      integer :: space = 0
      !! Subspace dimension at the end.
      real(dp), allocatable :: omega(:)
      !! (nroots) excitation energies, hartree, ascending.
      real(dp), allocatable :: xamp(:, :, :)
      !! (nocc, nvir, nroots) orthonormal TDA amplitudes X_ia.
      real(dp), allocatable :: tdip(:, :)
      !! (3, nroots) transition dipoles <0|r|n>, atomic units
      !! (sqrt(2) sum_ia X_ia <i|r|a>).
      real(dp), allocatable :: osc(:)
      !! (nroots) oscillator strengths (2/3) w |mu|^2.
      real(dp), allocatable :: resid(:)
      !! (nroots) the residual norms ||A X - w X|| at the end, from the
      !! solver's own bookkeeping; a caller who wants an independent figure
      !! multiplies X by `trc_tddft_op_t%apply`.
   end type trc_tddft_result_t

   type :: trc_tddft_op_t
      !! The matvec X -> A X of the module header, for any number of vectors
      !! in one integral pass. Holds the orbitals and the XC cache; the
      !! orbitals are device-resident for the life of the object (mapped on
      !! their own, never the type), the trial densities per call.
      integer :: nao = 0, nocc = 0, nvir = 0
      logical :: ready = .false.
      logical :: dft = .false.
      real(dp) :: cx = 1.0_dp
      !! Exact-exchange fraction applied to K.
      integer(int64) :: budget = 0
      real(dp), allocatable :: cocc(:, :), cvir(:, :), de(:, :)
      !! (nao, nocc), (nao, nvir), (nocc, nvir) orbital-energy differences.
      type(trc_xc_functional_t) :: func
      type(trc_xc_fxc_cache_t) :: cache
   contains
      procedure :: setup => tddft_op_setup
      procedure :: apply => tddft_op_apply
      procedure :: transition_dipoles => tddft_op_dipoles
      procedure :: release => tddft_op_release
   end type trc_tddft_op_t

contains

   pure function trc_tddft_refuse_reason(opts, func, nspin) result(msg)
      !!
      !! Why this calculation cannot run through trc_tddft, or an empty string if
      !! it can. `func` absent is Hartree-Fock (CIS); `nspin` is the reference's
      !! spin count (1 restricted; default 1).
      !!
      type(trc_tddft_options_t), intent(in) :: opts
      type(trc_xc_functional_t), intent(in), optional :: func
      integer, intent(in), optional :: nspin
      character(len=256) :: msg
      msg = ""
      if (present(nspin)) then
         if (nspin /= 1) then
            msg = "trc_tddft: only a closed-shell restricted reference is implemented (UHF/UKS is refused)"
            return
         end if
      end if
      if (.not. opts%singlet) then
         msg = "trc_tddft: only singlet excitations are implemented (triplets are refused)"
      else if (.not. opts%tda) then
         msg = "trc_tddft: only the Tamm-Dancoff approximation is implemented (tda = .false., the full RPA, " &
               //"is refused)"
      else if (opts%nroots < 1) then
         msg = "trc_tddft: nroots must be at least 1"
      else if (.not. present(func)) then
         msg = ""
      else if (len_trim(func%missing) > 0) then
         msg = "trc_tddft: functional '"//trim(func%name)//"' is not available: "//trim(func%missing)
      else if (func%omega /= 0.0_dp .or. func%exx_lr /= 0.0_dp) then
         msg = "trc_tddft: range-separated hybrids ('"//trim(func%name)//"') are refused -- the response needs " &
               //"a second, short-range-screened exchange term this module does not build"
      else if (func%nlc_b /= 0.0_dp .or. func%nlc_c /= 0.0_dp) then
         msg = "trc_tddft: VV10 nonlocal correlation ('"//trim(func%name)//"') has no response term here"
      else if (func%pt2_cos /= 0.0_dp .or. func%pt2_css /= 0.0_dp) then
         msg = "trc_tddft: double hybrids ('"//trim(func%name)//"') are refused -- the PT2 part of the " &
               //"excitation energy is not implemented"
      else if (func%family /= XC_FAMILY_LDA .and. func%family /= XC_FAMILY_GGA) then
         msg = "trc_tddft: only LDA and GGA kernels are implemented ('"//trim(func%name)//"' is a meta-GGA)"
      end if
   end function trc_tddft_refuse_reason

   subroutine tddft_op_setup(this, b, la, nocc, cmo, eps, dmat, err, func, xg, budget)
      !!
      !! Keep the orbitals (`cmo` (nao, nmo), occupied columns first, host) on
      !! the device and, for a functional, build the f_xc cache from the
      !! ground-state density `dmat` (nao, nao, host) on `xg`, which must
      !! be the grid the SCF converged on and already on the device. `func`
      !! absent: Hartree-Fock. `err` is set, and nothing is mapped, if the
      !! functional is refused or `xg` is missing.
      !!
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(:, :), eps(:), dmat(:, :)
      type(error_t), intent(out) :: err
      type(trc_xc_functional_t), intent(in), optional :: func
      type(trc_xc_grid_t), intent(in), optional :: xg
      integer(int64), intent(in), optional :: budget
      type(trc_tddft_options_t) :: probe
      character(len=256) :: reason
      integer :: n, nvir, i, a

      call this%release()
      if (present(func)) then
         reason = trc_tddft_refuse_reason(probe, func)
      else
         reason = ""
      end if
      if (len_trim(reason) > 0) then
         call err%set(ERROR_VALIDATION, trim(reason))
         return
      end if
      if (present(func) .and. .not. present(xg)) then
         call err%set(ERROR_VALIDATION, "trc_tddft: a functional needs the SCF's xc grid")
         return
      end if
      n = size(cmo, 1)
      nvir = size(cmo, 2) - nocc
      if (nocc < 1 .or. nvir < 1 .or. n /= b%nao .or. size(eps) < size(cmo, 2)) then
         call err%set(ERROR_VALIDATION, "trc_tddft: needs at least one occupied and one virtual orbital, " &
                      //"and orbitals of the basis's size")
         return
      end if
      this%nao = n
      this%nocc = nocc
      this%nvir = nvir
      if (present(budget)) this%budget = budget
      allocate (this%cocc(n, nocc), this%cvir(n, nvir), this%de(nocc, nvir))
      this%cocc = cmo(:, 1:nocc)
      this%cvir = cmo(:, nocc + 1:nocc + nvir)
      do a = 1, nvir
         do i = 1, nocc
            this%de(i, a) = eps(nocc + a) - eps(i)
         end do
      end do
      this%cx = 1.0_dp
      this%dft = present(func)
      if (this%dft) then
         this%func = func
         this%cx = func%exx
      end if
      !$acc enter data copyin(this%cocc, this%cvir, this%de)
      if (this%dft) call this%cache%setup(b, xg, func, dmat, la)
      this%ready = .true.
   end subroutine tddft_op_setup

   subroutine tddft_add_diag(nov, nvec, de, x, ax)
      !! ax += (e_a - e_i) x for every vector, over the amplitudes flattened
      !! to one index (the callers pass their (nocc, nvir, nvec) arrays by
      !! sequence association): a loop over small run-time extents is one
      !! index here, not a multi-index do concurrent.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: de(nov), x(nov, nvec)
      real(dp), intent(inout) :: ax(nov, nvec)
      integer :: idx, k
      do concurrent(k=1:nvec, idx=1:nov)
         ax(idx, k) = ax(idx, k) + de(idx)*x(idx, k)
      end do
   end subroutine tddft_add_diag

   subroutine tddft_op_apply(this, b, eri, la, nvec, x, ax, xg)
      !!
      !! ax(:, :, k) := A x(:, :, k) for k = 1..nvec, x and ax HOST arrays
      !! (nocc, nvir, nvec): ONE `fock_nosym_many` and ONE
      !! `fxc_apply_many` call whatever nvec is. Scratch is three
      !! (nao, nao, nvec) arrays (two for Hartree-Fock). `xg` is needed
      !! exactly when `setup` was given a functional.
      !!
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nvec
      real(dp), intent(in) :: x(this%nocc, this%nvir, nvec)
      real(dp), intent(out) :: ax(this%nocc, this%nvir, nvec)
      type(trc_xc_grid_t), intent(in), optional :: xg
      real(dp), allocatable :: tmat(:, :, :), gmat(:, :, :), vxc(:, :, :), half(:, :), tmp(:, :)
      integer :: n, nocc, nvir, k

      if (.not. this%ready) error stop "trc_tddft: operator used before setup"
      if (this%dft .and. .not. present(xg)) error stop "trc_tddft: a functional's matvec needs the xc grid"
      if (nvec < 1) return
      n = this%nao
      nocc = this%nocc
      nvir = this%nvir
      allocate (tmat(n, n, nvec), gmat(n, n, nvec), half(n, nvir), tmp(nocc, n))
      if (this%dft) then
         allocate (vxc(n, n, nvec))
      else
         allocate (vxc(1, 1, 1))
      end if
      !$acc enter data copyin(x) create(ax, tmat, gmat, vxc, half, tmp)
      do k = 1, nvec
         call la%gemm_mat(this%cocc, x(:, :, k), half)
         call la%gemm_mat(half, this%cvir, tmat(:, :, k), transb="T")
      end do
      !$acc update self(tmat)
      call eri%fock_nosym_many(b, nvec, tmat, gmat, k_scale=2.0_dp*this%cx, j_scale=2.0_dp, &
                               density_screen=.false.)
      !$acc update device(gmat)
      if (this%dft) call trc_xc_fxc_rks_apply_many(b, xg, this%cache, nvec, tmat, vxc, this%budget, la)
      do k = 1, nvec
         call la%gemm_mat(this%cocc, gmat(:, :, k), tmp, transa="T")
         call la%gemm_mat(tmp, this%cvir, ax(:, :, k))
         if (this%dft) then
            call la%gemm_mat(this%cocc, vxc(:, :, k), tmp, transa="T")
            call la%gemm_mat(tmp, this%cvir, ax(:, :, k), alpha=2.0_dp, beta=1.0_dp)
         end if
      end do
      call tddft_add_diag(nocc*nvir, nvec, this%de, x, ax)
      !$acc update self(ax)
      !$acc exit data delete(x, ax, tmat, gmat, vxc, half, tmp)
      deallocate (tmat, gmat, vxc, half, tmp)
   end subroutine tddft_op_apply

   subroutine tddft_op_dipoles(this, b, pl, la, nvec, x, tdip)
      !! tdip(:, k) = sqrt(2) sum_ia x(i, a, k) <i|r|a>, the singlet transition
      !! dipole of an excitation with unit-norm amplitudes, origin at zero
      !! (immaterial: <i|a> = 0). x host (nocc, nvir, nvec).
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nvec
      real(dp), intent(in) :: x(this%nocc, this%nvir, nvec)
      real(dp), intent(out) :: tdip(3, nvec)
      real(dp), allocatable :: mmat(:, :, :), dip(:, :, :), dmo(:, :, :), half(:, :)
      integer :: n, c, k

      if (.not. this%ready) error stop "trc_tddft: operator used before setup"
      n = this%nao
      ! All TRC_NMULT components are evaluated for the 3 dipoles: terco has no
      ! dipole-only multipole routine (trc_cphf/trc_cpks do the same), and this
      ! runs once per solve.
      allocate (mmat(n, n, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)
      allocate (dip(n, n, 3), dmo(this%nocc, this%nvir, 3), half(n, this%nvir))
      dip = mmat(:, :, 1:3)
      deallocate (mmat)
      !$acc enter data copyin(dip) create(dmo, half)
      do c = 1, 3
         call la%gemm_mat(dip(:, :, c), this%cvir, half)
         call la%gemm_mat(this%cocc, half, dmo(:, :, c), transa="T")
      end do
      !$acc update self(dmo)
      !$acc exit data delete(dip, dmo, half)
      do k = 1, nvec
         do c = 1, 3
            tdip(c, k) = sqrt(2.0_dp)*sum(x(:, :, k)*dmo(:, :, c))
         end do
      end do
      deallocate (dip, dmo, half)
   end subroutine tddft_op_dipoles

   subroutine tddft_op_release(this)
      class(trc_tddft_op_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cocc, this%cvir, this%de)
      end if
      call this%cache%release()
      if (allocated(this%cocc)) deallocate (this%cocc)
      if (allocated(this%cvir)) deallocate (this%cvir)
      if (allocated(this%de)) deallocate (this%de)
      this%ready = .false.
      this%dft = .false.
      this%nao = 0
      this%nocc = 0
      this%nvir = 0
      this%budget = 0
   end subroutine tddft_op_release

   subroutine tddft_fail(res, err, code, msg)
      type(trc_tddft_result_t), intent(inout) :: res
      type(error_t), intent(inout) :: err
      integer, intent(in) :: code
      character(len=*), intent(in) :: msg
      res%converged = .false.
      call err%set(code, msg)
   end subroutine tddft_fail

   subroutine tddft_orthonormalise(nov, m, mmax, vs, nc, tc, nnew)
      !! Orthonormalise the columns of tc(:, 1:nc) against vs(:, 1:m) and one
      !! another (classical Gram-Schmidt, twice), keeping those that stay
      !! above DROP_REL of their starting norm. The survivors are compacted into
      !! the first nnew columns of tc and appended to vs after column m. At most
      !! mmax - m are kept: vs has mmax columns, and when mmax is the problem
      !! dimension roundoff must not let a dependent vector past it.
      integer, intent(in) :: nov, m, mmax, nc
      real(dp), intent(inout) :: vs(nov, *), tc(nov, nc)
      integer, intent(out) :: nnew
      real(dp) :: n0, n1
      real(dp) :: proj(max(m + nc, 1))
      integer :: j, p, q, mm
      nnew = 0
      mm = m
      do j = 1, nc
         n0 = sqrt(sum(tc(:, j)**2))
         if (mm >= mmax) exit
         if (.not. (n0 > 0.0_dp)) cycle
         do p = 1, 2
            do q = 1, mm
               proj(q) = sum(vs(:, q)*tc(:, j))
            end do
            do q = 1, mm
               tc(:, j) = tc(:, j) - proj(q)*vs(:, q)
            end do
         end do
         n1 = sqrt(sum(tc(:, j)**2))
         if (.not. (n1 > DROP_REL*n0)) cycle
         mm = mm + 1
         nnew = nnew + 1
         vs(:, mm) = tc(:, j)/n1
      end do
      do j = 1, nnew
         tc(:, j) = vs(:, m + j)
      end do
   end subroutine tddft_orthonormalise

   subroutine tddft_probe(nov, nt, nprobe, nguess, mmax, ord, xr, axr, tp, m, vs, avs, nnew)
      !! Append the unit vectors of the nprobe orbital-energy differences after
      !! the nguess guessed ones (those the guesses never held), collapsing
      !! onto the nt Ritz vectors xr/axr first if there is no room. nnew is the
      !! number kept (0: all lay inside the space), m is advanced by it.
      integer, intent(in) :: nov, nt, nprobe, nguess, mmax
      integer, intent(in) :: ord(nov)
      real(dp), intent(in) :: xr(nov, nt), axr(nov, nt)
      real(dp), intent(inout) :: tp(nov, nprobe), vs(nov, *), avs(nov, *)
      integer, intent(inout) :: m
      integer, intent(out) :: nnew
      integer :: j
      tp = 0.0_dp
      do j = 1, nprobe
         tp(ord(nguess + j), j) = 1.0_dp
      end do
      if (m + nprobe > mmax) then
         vs(:, 1:nt) = xr
         avs(:, 1:nt) = axr
         m = nt
      end if
      call tddft_orthonormalise(nov, m, mmax, vs, nprobe, tp, nnew)
      m = m + nnew
   end subroutine tddft_probe

   subroutine trc_tddft_solve_rks(b, pl, eri, la, nocc, cmo, eps, dmat, opts, res, err, func, xg, budget)
      !!
      !! The lowest `opts%nroots` singlet TDA roots of a closed-shell reference
      !! (module header). `cmo` (nao, nmo) and `eps` (nmo) are the SCF's
      !! orbitals and energies, `dmat` (nao, nao) its total density, all host as
      !! `trc_scf_run` returns them; `func` and `xg` (the grid the SCF
      !! converged on, built and on the device) are given together for a
      !! Kohn-Sham reference and absent for Hartree-Fock. `err` carries a
      !! refusal, a breakdown, a NaN or a failure to converge; `res%converged` is
      !! then false and the other fields are not to be used.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(:, :), eps(:), dmat(:, :)
      type(trc_tddft_options_t), intent(in) :: opts
      type(trc_tddft_result_t), intent(out) :: res
      type(error_t), intent(out) :: err
      type(trc_xc_functional_t), intent(in), optional :: func
      type(trc_xc_grid_t), intent(in), optional :: xg
      integer(int64), intent(in), optional :: budget

      type(trc_tddft_op_t) :: op
      real(dp), allocatable :: vs(:, :), avs(:, :), hsub(:, :), csub(:, :), wsub(:)
      real(dp), allocatable :: xr(:, :), axr(:, :), tc(:, :), rnorm(:), dgn(:), dsort(:), tp(:, :)
      real(dp) :: den
      integer, allocatable :: ord(:)
      integer :: nov, nvir, nr, nt, mmax, nguess, m, nnew, nunc, nmain, i, a, j, k, ioff, nprobe
      integer(default_int) :: linfo
      character(len=256) :: reason
      character(len=1024) :: logline
      logical :: done, probed
      real(dp) :: buf_tol

      res%converged = .false.
      res%nroots = opts%nroots
      reason = trc_tddft_refuse_reason(opts, func)
      if (len_trim(reason) > 0) then
         call tddft_fail(res, err, ERROR_VALIDATION, trim(reason))
         return
      end if
      if (present(func) .neqv. present(xg)) then
         call tddft_fail(res, err, ERROR_VALIDATION, "trc_tddft: a functional and its xc grid are given together")
         return
      end if
      if (present(func)) then
         call op%setup(b, la, nocc, cmo, eps, dmat, err, func=func, xg=xg, budget=budget)
      else
         call op%setup(b, la, nocc, cmo, eps, dmat, err, budget=budget)
      end if
      if (err%has_error()) then
         res%converged = .false.
         call op%release()
         return
      end if

      nvir = op%nvir
      nov = nocc*nvir
      nr = opts%nroots
      if (nr > nov) then
         call tddft_fail(res, err, ERROR_VALIDATION, "trc_tddft: nroots = "//to_string(nr)// &
                         " exceeds the "//to_string(nov)//" occupied-virtual pairs")
         call op%release()
         return
      end if
      mmax = opts%max_space
      if (mmax <= 0) mmax = 12*nr
      ! The block iterated is the nroots wanted plus a buffer of further Ritz
      ! vectors (see the header, MISSED ROOTS).
      nt = min(nr + max(2, nr/2), nov)
      buf_tol = max(BUFFER_TOL, 100.0_dp*opts%conv_tol)
      mmax = min(max(mmax, 3*nt), nov)
      nguess = opts%nguess
      if (nguess <= 0) nguess = max(2*nr, nr + 4)
      nguess = min(max(nguess, nt), nov, mmax)

      allocate (vs(nov, mmax), avs(nov, mmax), hsub(mmax, mmax), csub(mmax, mmax), wsub(mmax))
      allocate (xr(nov, nt), axr(nov, nt), tc(nov, nt), rnorm(nt), dgn(nov), dsort(nov))
      allocate (ord(nov))

      ! Diagonal (e_a - e_i) in the flat (i fastest) ordering of the amplitudes;
      ! the guesses are its lowest entries.
      do a = 1, nvir
         do i = 1, nocc
            dgn(i + (a - 1)*nocc) = op%de(i, a)
         end do
      end do
      ! sort_index sorts its array in place; the preconditioner needs dgn as is.
      dsort = dgn
      call sort_index(dsort, ord)
      ! Do not cut a cluster of (near-)degenerate differences in two: the
      ! partner left out is a direction the solve may never generate.
      do while (nguess < min(nov, mmax))
         if (dsort(nguess + 1) - dsort(nguess) >= TIE_TOL) exit
         nguess = nguess + 1
      end do
      ! The probe must fit beside the Ritz vectors it may collapse onto.
      nprobe = max(min(nguess, nov - nguess, mmax - nt), 0)
      allocate (tp(nov, max(nprobe, 1)))
      probed = .false.
      vs(:, 1:nguess) = 0.0_dp
      do j = 1, nguess
         vs(ord(j), j) = 1.0_dp
      end do
      m = nguess
      nnew = nguess
      res%iterations = 0
      res%nmatvec = 0
      done = .false.

      do
         ! A on the vectors added since the last pass: one integral pass.
         if (nnew > 0) then
            call op%apply(b, eri, la, nnew, vs(:, m - nnew + 1:m), avs(:, m - nnew + 1:m), xg)
            res%nmatvec = res%nmatvec + nnew
            if (.not. all(ieee_is_finite(avs(:, m - nnew + 1:m)))) then
               call tddft_fail(res, err, ERROR_GENERIC, "trc_tddft: the matvec returned NaN or Inf")
               exit
            end if
         end if
         res%iterations = res%iterations + 1

         ! Subspace problem.
         call pic_gemm(vs(:, 1:m), avs(:, 1:m), hsub(1:m, 1:m), "T", "N", 1.0_dp, 0.0_dp)
         do j = 1, m
            do k = 1, j - 1
               hsub(k, j) = 0.5_dp*(hsub(k, j) + hsub(j, k))
            end do
         end do
         csub(1:m, 1:m) = hsub(1:m, 1:m)
         call pic_syev(csub(1:m, 1:m), wsub(1:m), "V", "U", linfo)
         if (linfo /= 0 .or. .not. all(ieee_is_finite(wsub(1:m)))) then
            call tddft_fail(res, err, ERROR_GENERIC, "trc_tddft: the subspace eigenproblem broke down (NaN or " &
                            //"LAPACK failure)")
            exit
         end if

         ! Ritz vectors, A times them, and the residuals of the lowest nt.
         call pic_gemm(vs(:, 1:m), csub(1:m, 1:nt), xr, "N", "N", 1.0_dp, 0.0_dp)
         call pic_gemm(avs(:, 1:m), csub(1:m, 1:nt), axr, "N", "N", 1.0_dp, 0.0_dp)
         nunc = 0
         nmain = 0
         do k = 1, nt
            tc(:, k) = axr(:, k) - wsub(k)*xr(:, k)
            rnorm(k) = sqrt(sum(tc(:, k)**2))
            if (.not. (rnorm(k) < merge(opts%conv_tol, buf_tol, k <= nr))) then
               nunc = nunc + 1
               if (k <= nr) nmain = nmain + 1
            end if
         end do
         if (opts%verbose) then
            write (logline, "(a,i0,a,i0,a,i0,a,es9.2,a,f14.8)") "   tddft it ", res%iterations, ": space ", m, &
               ", unconverged ", nunc, ", max resid ", maxval(rnorm(1:nr)), ", w1 ", wsub(1)
            call global_logger%info(trim(logline))
         end if
         if (.not. all(ieee_is_finite(rnorm))) then
            call tddft_fail(res, err, ERROR_GENERIC, "trc_tddft: a residual is NaN or Inf")
            exit
         end if
         if (nunc == 0) then
            ! One probe before accepting: unit vectors of the next-lowest
            ! differences, none of which the guesses held. A root whose
            ! character the guesses (and so the corrections, which stay in
            ! their span) never reached shows up as a Ritz value below the
            ! converged ones and the iteration resumes. One more matvec pass.
            if (.not. probed .and. nprobe > 0) then
               probed = .true.
               call tddft_probe(nov, nt, nprobe, nguess, mmax, ord, xr, axr, tp, m, vs, avs, nnew)
               if (nnew > 0) cycle
            end if
            done = .true.
            exit
         end if
         if (res%iterations >= opts%max_iter) then
            call tddft_fail(res, err, ERROR_GENERIC, "trc_tddft: not converged in "//to_string(opts%max_iter)// &
                            " iterations, residual "//to_string(maxval(rnorm(1:nr)), "es9.2")//" (conv_tol "// &
                            to_string(opts%conv_tol, "es9.2")//")")
            exit
         end if

         ! Preconditioned corrections of the unconverged roots, packed first.
         ioff = 0
         do k = 1, nt
            if (rnorm(k) < merge(opts%conv_tol, buf_tol, k <= nr)) cycle
            ioff = ioff + 1
            do j = 1, nov
               den = dgn(j) - wsub(k)
               if (abs(den) < PRECOND_FLOOR) den = sign(PRECOND_FLOOR, den)
               tc(j, ioff) = -(axr(j, k) - wsub(k)*xr(j, k))/den
            end do
         end do

         ! Collapse onto the current Ritz vectors when there is no room.
         if (m + nunc > mmax) then
            vs(:, 1:nt) = xr
            avs(:, 1:nt) = axr
            m = nt
         end if
         call tddft_orthonormalise(nov, m, mmax, vs, nunc, tc, nnew)
         if (nnew == 0 .and. nmain == 0) then
            ! Only buffer vectors were still unconverged and their corrections
            ! lie inside the space: nothing more to learn from them. The probe
            ! still runs once, as on the ordinary converged path.
            if (.not. probed .and. nprobe > 0) then
               probed = .true.
               call tddft_probe(nov, nt, nprobe, nguess, mmax, ord, xr, axr, tp, m, vs, avs, nnew)
               if (nnew > 0) cycle
            end if
            done = .true.
            exit
         end if
         if (nnew == 0) then
            call tddft_fail(res, err, ERROR_GENERIC, "trc_tddft: breakdown, no usable new trial vector while "// &
                            to_string(nunc)//" root(s) are unconverged (residual "// &
                            to_string(maxval(rnorm(1:nr)), "es9.2")//")")
            exit
         end if
         m = m + nnew
      end do

      if (done) then
         if (wsub(1) <= 0.0_dp) then
            call tddft_fail(res, err, ERROR_GENERIC, "trc_tddft: a non-positive excitation energy ("// &
                            to_string(wsub(1), "es10.2")//" Eh): the reference is unstable (or a triplet-like " &
                            //"instability)")
         else
            res%space = m
            allocate (res%omega(nr), res%xamp(nocc, nvir, nr), res%tdip(3, nr), res%osc(nr), res%resid(nr))
            res%omega = wsub(1:nr)
            res%resid = rnorm(1:nr)
            do k = 1, nr
               xr(:, k) = xr(:, k)/sqrt(sum(xr(:, k)**2))
               res%xamp(:, :, k) = reshape(xr(:, k), [nocc, nvir])
            end do
            call op%transition_dipoles(b, pl, la, nr, res%xamp, res%tdip)
            do k = 1, nr
               res%osc(k) = OSC_PREFACTOR*res%omega(k)*sum(res%tdip(:, k)**2)
            end do
            res%converged = .true.
         end if
      end if
      call op%release()
      deallocate (vs, avs, hsub, csub, wsub, xr, axr, tc, rnorm, dgn, dsort, tp, ord)
   end subroutine trc_tddft_solve_rks

end module trc_tddft
