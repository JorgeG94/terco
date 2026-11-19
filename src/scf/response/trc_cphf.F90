module trc_cphf
   !!
   !! Coupled-perturbed Hartree-Fock: the orbital-response linear system this
   !! milestone solves, and the static dipole polarizability built on top of it
   !! as the validating property. Intended foundation for RI-MP2's Z-vector
   !! equation and, later, a coupled-perturbed Kohn-Sham variant for
   !! double-hybrid gradients -- neither is implemented here, this module only
   !! solves the ground-state orbital Hessian system and the one property that
   !! needs nothing but it.
   !!
   !! THE EQUATION (RHF)
   !! -------------------
   !! For a real, symmetric one-electron perturbation (a static field has no
   !! AO-basis response, so this is the whole right-hand side), the response
   !! orbital-rotation amplitudes U_ai solve
   !!
   !!     (e_a - e_i) U_ai + [A U]_ai = -B_ai
   !!
   !! with the orbital Hessian action, built from one resident Fock-like build
   !! rather than an AO integral derivative:
   !!
   !!     D_U = C_vir U C_occ^T + (C_vir U C_occ^T)^T          (symmetric, AO)
   !!     [A U]_ai = 2 [C_vir^T G[D_U] C_occ]_ai,   G[D_U] = J[D_U] - K[D_U]/2
   !!
   !! i.e. exactly eri%fock_resident(b, D_U, G_U, j_scale=1, k_scale=1,
   !! density_screen=.false.) -- the trial density is not the SCF's physical
   !! one, so it must not drive the Schwarz screen (trc_eri.F90's own warning
   !! on both `fock_resident` and `fock`). The derivation (4(ai|jb) -
   !! (ij|ab) - (ib|ja), contracted against the full permutation symmetry of
   !! real two-electron integrals) is the standard RHF orbital-Hessian
   !! result -- Pople, Krishnan, Schlegel and Binkley, Int. J. Quantum Chem.
   !! 16, 225 (1979), or Helgaker/Jorgensen/Olsen, Molecular Electronic-
   !! Structure Theory, the CPHF section of ch. 10 -- not rederived here;
   !! see trc_eri.F90's header for why `fock_resident` is the entry point
   !! this was always meant to reach ("the coupled-perturbed equations run
   !! a hundred of them").
   !!
   !! THE EQUATION (UHF)
   !! -------------------
   !! Each spin channel's Fock response only picks up a J contribution from
   !! the OTHER spin's trial density; its own K is spin-local:
   !!
   !!     D_U^s = C_vir,s U^s C_occ,s^T + transpose,   D_U^tot = D_U^a + D_U^b
   !!     [A U]_ai^s = [C_vir,s^T (J[D_U^tot] - K[D_U^s]) C_occ,s]_ai
   !!
   !! J[D_U^tot] is one fock_resident call (j_scale=1, k_scale=0) shared by
   !! both spins; K[D_U^s] is fock_resident(j_scale=0, k_scale=2) per spin --
   !! the same "k_scale=2 on the undoubled density" trick trc_scf_fock.F90's
   !! `scf_fock_build_g` uses for the ground-state UKS exchange term. No extra factor
   !! of 2 here: RHF's 2 comes from folding two identical spin channels into
   !! one solve, which UHF does not do.
   !!
   !! THE SOLVER
   !! ----------
   !! Preconditioned conjugate gradient, the diagonal (e_a - e_i) as the
   !! preconditioner (SPD for a stable ground state, so CG is the right
   !! choice over a general Krylov method), with an iteration count reported
   !! per column. The right-hand sides are solved TOGETHER: one CG loop over
   !! every column that has not converged, each with its own alpha, beta and
   !! residual (so each column's iterates, and its iteration count, are those
   !! of a solve of that column alone), whose matvecs are built in ONE pass
   !! over the integrals per iteration (`eri%fock_resident_many_lent`,
   !! density_screen off) rather than one pass per column. A column drops out
   !! of the batch as it converges, so the last iterations are cheap. The
   !! Fock pass costs less for several densities than for the same number of
   !! passes (the quartet is evaluated once, the digestion is per density),
   !! which is what the sharing buys; the Z-vector of every root of an
   !! excited-state gradient is the intended customer.
   !! `trc_cphf_options_t%max_batch` caps the columns in flight (memory: about
   !! 4 nao^2 doubles each for HF, 5 with the XC response, the integrals'
   !! working basis counted for nao: more for a decontracted basis) and 1
   !! reverts to one column per pass. `trc_cphf_cg_t` is the shared CG state
   !! machine: a solver builds the matvec for `cg%act(1:cg%na)` and calls
   !! `advance`.
   !! The density-fitted solves (`trc_cphf_solve_rhf_fitted`) and the UHF
   !! solve are still one column at a time: the MO-basis fitted response has
   !! no multi-vector form, and UHF's two-density matvec is not batched here.
   !!
   !! RESIDENCY. `cmo`, the right-hand sides and the solution are assumed
   !! already on the device (the same contract eri%fock_resident itself
   !! documents); only `eps` (nao-size) and the CG's own scalars (dot
   !! products, read back by trc_linalg%dot_2d, exactly as DIIS's history
   !! does) cross the bus. Every per-iteration matrix -- D_U, G_U, the GEMM
   !! scratch, the CG vectors themselves -- stays resident for the life of
   !! the solve.
   !!
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_multipoles, TRC_NMULT
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_settings, only: trc_scf_result_t
   use trc_scf_df, only: df_jk_t
   use trc_scf_df_response, only: df_rotation_t
   implicit none
   private

   public :: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_cg_t
   public :: trc_cphf_solve_rhf, trc_cphf_solve_uhf, trc_cphf_solve_rhf_fitted
   public :: trc_scf_polarizability_rhf, trc_scf_polarizability_uhf, trc_scf_polarizability
   ! Shared with trc_cpks (the CPKS extension of this solver): the same CG
   ! loop, the same elementwise `do concurrent` helpers -- see this
   ! module's header on why they are split out as explicit-shape-dummy
   ! subroutines rather than inlined.
   public :: cphf_zero_matrix, cphf_copy_matrix, cphf_neg_matrix, cphf_add_matrix, cphf_axpy_matrix, &
             cphf_precond_matrix, cphf_combine_p, cphf_add_diag

   type :: trc_cphf_options_t
      real(dp) :: conv_residual = 1.0e-8_dp
      !! Converged when the Frobenius norm of the residual (e_a-e_i) U +
      !! [A U] + B drops below this, for every right-hand side.
      integer :: max_iteration = 100
      logical :: verbose = .false.
      integer :: max_batch = 8
      !! Right-hand sides whose matvecs share one pass over the integrals
      !! (RHF and RKS exact solves). The scratch is that many Fock-sized
      !! matrices (about 4 nao^2 doubles each for HF, 5 with the XC response),
      !! so the cap keeps it modest for a large basis, where the extra
      !! columns buy little: the integral evaluation is shared, but each
      !! density's digestion is not. Eight is about 250 MB for HF and 320 MB
      !! for KS at nao = 1000, more for a decontracted basis (the Fock pass's
      !! scratch is sized by the working basis). Right-hand sides beyond it
      !! are solved in further groups of at most this many. 1 solves one
      !! column at a time, the previous behaviour; a value below 1 is read as 1.
   end type trc_cphf_options_t

   type :: trc_cphf_cg_t
      !! Preconditioned CG over a group of right-hand sides at once, one
      !! independent recurrence per column: the same arithmetic, column by
      !! column, as a loop that solves them one after another, with the
      !! matvec left to the caller so that every column still iterating can
      !! be fed to one batched Fock build. Usage: `init`, then per group
      !! `begin`; while `na > 0`, fill `ap(:, :, act(k))` with A p for
      !! k = 1..na (`p` holds each column's search direction) and `advance`;
      !! then `report`; finally `release`. Device arrays (r, z, p, ap,
      !! diagv) are mapped on their own, never the type.
      integer :: nvirt = 0, nocc = 0, nb = 0, nc = 0
      !! The columns the state holds, and those of the current group.
      logical :: started = .false.
      !! Set by `begin`; `advance` before it is an error.
      integer :: na = 0
      !! Columns still iterating; their group-local indices are act(1:na).
      integer :: c0 = 0
      !! First column (of the caller's right-hand sides) of the current group.
      integer :: maxit = 0
      real(dp) :: tol = 0.0_dp
      integer, allocatable :: act(:), it(:)
      logical, allocatable :: broke(:)
      real(dp), allocatable :: rz(:), resid(:)
      real(dp), allocatable :: r(:, :, :), z(:, :, :), p(:, :, :), ap(:, :, :), diagv(:, :)
      !! r, z, p, ap are (nvirt, nocc, nb); diagv (nvirt, nocc).
   contains
      procedure :: init => cg_init
      procedure :: begin => cg_begin
      procedure :: advance => cg_advance
      procedure :: report => cg_report
      procedure :: release => cg_release
   end type trc_cphf_cg_t

   type :: trc_cphf_result_t
      logical :: converged = .true.
      integer, allocatable :: iterations(:)   !! CG iterations, one per right-hand side
      character(len=256) :: message = ""
   end type trc_cphf_result_t

contains

   subroutine trc_cphf_solve_rhf(b, eri, lp, nao, nocc, nvirt, cmo, eps, ncomp, bvec, u, opts, res)
      !!
      !! RHF CPHF: `ncomp` right-hand sides, solved together (this module's
      !! own header: one CG recurrence per column, the matvecs of the columns
      !! still iterating built in one batched Fock pass, in groups of at most
      !! `opts%max_batch`). `cmo` (nao, nao), `bvec` (nvirt, nocc, ncomp) and `u`
      !! (same shape, overwritten with the solution) must already be
      !! device-resident; `eps` (nao) is host, read once to build the diagonal
      !! preconditioner. `res%iterations(c)` is column c's own count.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, nocc, nvirt, ncomp
      real(dp), intent(in) :: cmo(nao, nao)
      !! Occupied columns 1:nocc, virtual nocc+1:nao, ascending orbital
      !! energy -- trc_scf_driver's own convention.
      real(dp), intent(in) :: eps(nao)
      real(dp), intent(in) :: bvec(nvirt, nocc, ncomp)
      real(dp), intent(inout) :: u(nvirt, nocc, ncomp)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res

      real(dp), allocatable :: d_u(:, :, :), d_ut(:, :), g_u(:, :, :)
      real(dp), allocatable :: temporary_no(:, :), temporary_no2(:, :), work(:)
      type(trc_cphf_cg_t) :: cg
      type(trc_fock_options_t) :: fopts
      integer :: nbmax, c0, nc

      nbmax = max(1, min(opts%max_batch, ncomp))
      allocate (d_u(nao, nao, nbmax), d_ut(nao, nao), g_u(nao, nao, nbmax))
      allocate (temporary_no(nao, nocc), temporary_no2(nao, nocc))
      !$acc enter data create(d_u, d_ut, g_u, temporary_no, temporary_no2)
      if (nbmax > 1) then
         allocate (work(eri%fock_work_size(nbmax)))
         !$acc enter data create(work)
      end if
      ! density_screen off: a trial vector must not drive the screen
      ! (this module's header).
      fopts = trc_fock_options_t(k_scale=1.0_dp, j_scale=1.0_dp, density_screen=.false.)
      call cg%init(nvirt, nocc, nbmax, eps, opts)

      allocate (res%iterations(ncomp))
      res%iterations = 0
      res%converged = .true.
      res%message = ""

      c0 = 1
      do while (c0 <= ncomp)
         nc = min(nbmax, ncomp - c0 + 1)
         call cg%begin(lp, ncomp, bvec, u, c0, nc)
         do while (cg%na > 0)
            call matvec()
            call cg%advance(lp, ncomp, u)
         end do
         call cg%report(res, "cphf", opts%verbose)
         c0 = c0 + nc
      end do

      call cg%release()
      if (nbmax > 1) then
         !$acc exit data delete(work)
         deallocate (work)
      end if
      !$acc exit data delete(d_u, d_ut, g_u, temporary_no, temporary_no2)

   contains

      subroutine matvec()
         !! [A p]_ai = 2 [C_vir^T G[D_p] C_occ]_ai + (e_a-e_i) p_ai for every
         !! active column, D_p the symmetrised AO density this p rotates into
         !! (this module's header): all the D_p, one Fock pass, all the A p.
         integer :: k, kk
         do k = 1, cg%na
            kk = cg%act(k)
            call lp%gemm_mat(cmo(:, nocc + 1:nao), cg%p(:, :, kk), temporary_no)
            call lp%gemm_mat(temporary_no, cmo(:, 1:nocc), d_u(:, :, k), transb="T")
            call lp%transpose(nao, nao, d_u(:, :, k), nao, d_ut, nao)
            call cphf_add_matrix(nao, nao, d_u(:, :, k), d_ut, d_u(:, :, k))
         end do
         if (nbmax == 1) then
            call eri%fock_resident(b, d_u(:, :, 1), g_u(:, :, 1), opts=fopts)
         else
            call eri%fock_resident_many_lent(b, cg%na, d_u, g_u, fopts, work)
         end if
         do k = 1, cg%na
            kk = cg%act(k)
            call lp%gemm_mat(g_u(:, :, k), cmo(:, 1:nocc), temporary_no2)
            call lp%gemm_mat(cmo(:, nocc + 1:nao), temporary_no2, cg%ap(:, :, kk), transa="T", alpha=2.0_dp)
            call cphf_add_diag(nvirt, nocc, cg%diagv, cg%p(:, :, kk), cg%ap(:, :, kk))
         end do
      end subroutine matvec

   end subroutine trc_cphf_solve_rhf

   subroutine cg_init(this, nvirt, nocc, nb, eps, opts)
      !! Size the state for groups of up to `nb` columns, and map it. `eps`
      !! (host, indexed over all orbitals, occupied first) gives the diagonal
      !! preconditioner (e_a - e_i); `opts` the tolerance and iteration cap.
      class(trc_cphf_cg_t), intent(inout) :: this
      integer, intent(in) :: nvirt, nocc, nb
      real(dp), intent(in) :: eps(:)
      type(trc_cphf_options_t), intent(in) :: opts
      integer :: i, a
      call this%release()
      this%nvirt = nvirt
      this%nocc = nocc
      this%nb = nb
      this%maxit = opts%max_iteration
      this%tol = opts%conv_residual
      allocate (this%act(nb), this%it(nb), this%broke(nb), this%rz(nb), this%resid(nb))
      allocate (this%r(nvirt, nocc, nb), this%z(nvirt, nocc, nb), this%p(nvirt, nocc, nb), &
                this%ap(nvirt, nocc, nb), this%diagv(nvirt, nocc))
      do i = 1, nocc
         do a = 1, nvirt
            this%diagv(a, i) = eps(nocc + a) - eps(i)
         end do
      end do
      this%it = 0
      this%broke = .false.
      this%na = 0
      !$acc enter data create(this%r, this%z, this%p, this%ap)
      !$acc enter data copyin(this%diagv)
   end subroutine cg_init

   subroutine cg_begin(this, lp, ncomp, bvec, u, c0, nc)
      !! Start columns c0 .. c0+nc-1 of `bvec` (nvirt, nocc, ncomp), device
      !! resident, as the group's columns 1..nc: u = 0, r = -b, z = r/diag,
      !! p = z, and the loop-entry test of each column (`.not. (resid < tol)`,
      !! so a NaN right-hand side starts iterating and is caught by the
      !! breakdown guard of `advance` instead of passing as converged).
      class(trc_cphf_cg_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: ncomp, c0, nc
      real(dp), intent(in) :: bvec(this%nvirt, this%nocc, ncomp)
      real(dp), intent(inout) :: u(this%nvirt, this%nocc, ncomp)
      integer :: k, nv, no
      nv = this%nvirt
      no = this%nocc
      if (nc > this%nb) error stop "trc_cphf: a CG group larger than the state was sized for"
      this%c0 = c0
      this%na = 0
      this%started = .true.
      do k = 1, nc
         call cphf_neg_matrix(nv, no, bvec(:, :, c0 + k - 1), this%r(:, :, k))
         call cphf_zero_matrix(nv, no, u(:, :, c0 + k - 1))
         call cphf_precond_matrix(nv, no, this%diagv, this%r(:, :, k), this%z(:, :, k))
         call cphf_copy_matrix(nv, no, this%z(:, :, k), this%p(:, :, k))
         this%rz(k) = lp%dot_2d(this%r(:, :, k), this%z(:, :, k))
         this%resid(k) = sqrt(lp%dot_2d(this%r(:, :, k), this%r(:, :, k)))
         this%it(k) = 0
         this%broke(k) = .false.
         if (this%it(k) < this%maxit .and. .not. (this%resid(k) < this%tol)) then
            this%na = this%na + 1
            this%act(this%na) = k
         end if
      end do
      this%nc = nc
   end subroutine cg_begin

   subroutine cg_advance(this, lp, ncomp, u)
      !! One CG step for every active column, given `ap` = A p for each
      !! (filled by the caller): that column's own alpha = rz/(p.Ap), its
      !! own beta, residual and iteration count. A column leaves the active
      !! set when its residual drops below tol (converged), at the iteration
      !! cap, or on a breakdown (p.Ap <= 0, or NaN -- NaN > 0 is false): then
      !! it stops BEFORE the update, never falling through to an Inf/NaN alpha
      !! that the residual test could not catch either (trc_cphf_solve_rhf's
      !! old comment, which still holds).
      class(trc_cphf_cg_t), intent(inout) :: this
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: ncomp
      real(dp), intent(inout) :: u(this%nvirt, this%nocc, ncomp)
      integer :: j, k, nn, nv, no, col
      real(dp) :: pap, alpha_cg, beta_cg, rz_old
      if (.not. this%started) error stop "trc_cphf: CG advance before begin"
      nv = this%nvirt
      no = this%nocc
      nn = 0
      do j = 1, this%na
         k = this%act(j)
         col = this%c0 + k - 1
         pap = lp%dot_2d(this%p(:, :, k), this%ap(:, :, k))
         if (.not. (pap > 0.0_dp)) then
            this%broke(k) = .true.
            cycle
         end if
         alpha_cg = this%rz(k)/pap
         call cphf_axpy_matrix(nv, no, alpha_cg, this%p(:, :, k), u(:, :, col))
         call cphf_axpy_matrix(nv, no, -alpha_cg, this%ap(:, :, k), this%r(:, :, k))
         this%resid(k) = sqrt(lp%dot_2d(this%r(:, :, k), this%r(:, :, k)))
         this%it(k) = this%it(k) + 1
         if (this%resid(k) < this%tol) cycle
         call cphf_precond_matrix(nv, no, this%diagv, this%r(:, :, k), this%z(:, :, k))
         rz_old = this%rz(k)
         this%rz(k) = lp%dot_2d(this%r(:, :, k), this%z(:, :, k))
         beta_cg = this%rz(k)/rz_old
         call cphf_combine_p(nv, no, this%z(:, :, k), beta_cg, this%p(:, :, k))
         if (this%it(k) < this%maxit) then
            nn = nn + 1
            this%act(nn) = k
         end if
      end do
      this%na = nn
   end subroutine cg_advance

   subroutine cg_report(this, res, label, verbose)
      !! Copy the finished group's iteration counts into `res%iterations`
      !! and fold its outcome into res%converged and res%message, column by
      !! column in order (so the message is that of the last column that
      !! failed, as when the columns were solved one after another).
      !! `label` ("cphf", "cpks") names the solver in the message and log.
      class(trc_cphf_cg_t), intent(in) :: this
      type(trc_cphf_result_t), intent(inout) :: res
      character(len=*), intent(in) :: label
      logical, intent(in) :: verbose
      integer :: k
      character(len=1024) :: logline
      do k = 1, this%nc
         res%iterations(this%c0 + k - 1) = this%it(k)
         if (this%broke(k)) then
            res%converged = .false.
            res%message = "trc_"//label//": CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (this%resid(k) < this%tol)) then
            res%converged = .false.
            res%message = "trc_"//label//": a right-hand side did not converge"
         end if
         if (verbose) then
            write (logline, "(a,a,i0,a,i0,a,es9.2)") "   ", label//" rhs ", this%c0 + k - 1, ": ", this%it(k), &
               " iterations, residual ", this%resid(k)
            call global_logger%info(trim(logline))
         end if
      end do
   end subroutine cg_report

   subroutine cg_release(this)
      class(trc_cphf_cg_t), intent(inout) :: this
      if (allocated(this%r)) then
         !$acc exit data delete(this%diagv)
         !$acc exit data delete(this%r, this%z, this%p, this%ap)
      end if
      if (allocated(this%act)) deallocate (this%act, this%it, this%broke, this%rz, this%resid)
      if (allocated(this%r)) deallocate (this%r, this%z, this%p, this%ap, this%diagv)
      this%na = 0
      this%nc = 0
      this%started = .false.
   end subroutine cg_release

   subroutine trc_cphf_solve_rhf_fitted(resp, dfjk, lp, bvec, u, opts, res)
      !!
      !! RHF CPHF with the DENSITY-FITTED orbital Hessian: the same equation
      !! and the same preconditioned CG as trc_cphf_solve_rhf, but
      !! [A U] - (e_a - e_i) U = 2 C_v^T (J_fit - K_fit/2)[D_U] C_o from `resp`
      !! (trc_scf_df_response's df_rotation_t, over `dfjk`'s fit and the
      !! orbitals `resp` was set up with: the MO-basis fitted orbital Hessian,
      !! or its AO form when the MO blocks do not fit). This is the response
      !! operator of a --df SCF: its orbitals are stationary for the fitted
      !! energy, so a Z-vector or a property over that reference has to use it,
      !! not the exact one. No screening is involved at all.
      !!
      !! Dimensions come from `resp` and the arrays: bvec/u (nvirt, nocc,
      !! ncomp), nocc and nvirt `resp`'s own. bvec and u device-resident. Same
      !! convergence test, breakdown and NaN handling as trc_cphf_solve_rhf.
      !!
      type(df_rotation_t), intent(inout) :: resp
      type(df_jk_t), intent(in) :: dfjk
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in) :: bvec(:, :, :)
      real(dp), intent(inout) :: u(:, :, :)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res

      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :), diagv(:, :)
      real(dp) :: rz, rz_old, pap, alpha_cg, beta_cg, resid
      integer :: nocc, nvirt, ncomp, it, comp, i, a
      logical :: broke

      nvirt = resp%nvir
      nocc = resp%nocc
      ncomp = size(bvec, 3)
      if (size(bvec, 1) /= nvirt .or. size(bvec, 2) /= nocc) error stop &
         "trc_cphf_solve_rhf_fitted: the right-hand side does not match the response's orbitals"
      allocate (r(nvirt, nocc), z(nvirt, nocc), p(nvirt, nocc), ap(nvirt, nocc), diagv(nvirt, nocc))
      do i = 1, nocc
         do a = 1, nvirt
            diagv(a, i) = resp%eps(nocc + a) - resp%eps(i)
         end do
      end do
      !$acc enter data create(r, z, p, ap)
      !$acc enter data copyin(diagv)

      allocate (res%iterations(ncomp))
      res%converged = .true.
      res%message = ""
      do comp = 1, ncomp
         call cphf_neg_matrix(nvirt, nocc, bvec(:, :, comp), r)
         call cphf_zero_matrix(nvirt, nocc, u(:, :, comp))
         call cphf_precond_matrix(nvirt, nocc, diagv, r, z)
         call cphf_copy_matrix(nvirt, nocc, z, p)
         rz = lp%dot_2d(r, z)
         resid = sqrt(lp%dot_2d(r, r))
         it = 0
         broke = .false.
         ! `.not. (resid < conv)`: a NaN residual must read as not converged
         ! (trc_cphf_solve_rhf's own comment).
         do while (it < opts%max_iteration .and. .not. (resid < opts%conv_residual))
            ! [A p] = 2 C_v^T G_fit[D_p] C_o + (e_a - e_i) p
            call resp%apply(dfjk, lp, p, ap)
            call cphf_add_diag(nvirt, nocc, diagv, p, ap)
            pap = lp%dot_2d(p, ap)
            if (.not. (pap > 0.0_dp)) then
               broke = .true.
               exit
            end if
            alpha_cg = rz/pap
            call cphf_axpy_matrix(nvirt, nocc, alpha_cg, p, u(:, :, comp))
            call cphf_axpy_matrix(nvirt, nocc, -alpha_cg, ap, r)
            resid = sqrt(lp%dot_2d(r, r))
            it = it + 1
            if (resid < opts%conv_residual) exit
            call cphf_precond_matrix(nvirt, nocc, diagv, r, z)
            rz_old = rz
            rz = lp%dot_2d(r, z)
            beta_cg = rz/rz_old
            call cphf_combine_p(nvirt, nocc, z, beta_cg, p)
         end do
         res%iterations(comp) = it
         if (broke) then
            res%converged = .false.
            res%message = "trc_cphf: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < opts%conv_residual)) then
            res%converged = .false.
            res%message = "trc_cphf: a right-hand side did not converge"
         end if
         if (opts%verbose) then
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,i0,a,es9.2)") &
                  "   cphf (fitted) rhs ", comp, ": ", it, " iterations, residual ", resid
               call global_logger%info(trim(logline))
            end block
         end if
      end do

      !$acc exit data delete(diagv)
      !$acc exit data delete(r, z, p, ap)
   end subroutine trc_cphf_solve_rhf_fitted

   subroutine trc_cphf_solve_uhf(b, eri, lp, nao, noa, nva, nob, nvb, cmo_a, cmo_b, eps_a, eps_b, &
                                 ncomp, bvec_a, bvec_b, u_a, u_b, opts, res)
      !!
      !! UHF CPHF: the two spin channels solved together (one combined CG over
      !! the stacked [U^a; U^b] vector), because the matvec couples them
      !! through J -- see this module's own header. `cmo_a`/`cmo_b`, the
      !! right-hand sides and the solutions must already be device-resident;
      !! `eps_a`/`eps_b` are host.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nao, noa, nva, nob, nvb, ncomp
      real(dp), intent(in) :: cmo_a(nao, nao), cmo_b(nao, nao)
      real(dp), intent(in) :: eps_a(nao), eps_b(nao)
      real(dp), intent(in) :: bvec_a(nva, noa, ncomp), bvec_b(nvb, nob, ncomp)
      real(dp), intent(inout) :: u_a(nva, noa, ncomp), u_b(nvb, nob, ncomp)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res

      real(dp), allocatable :: d_ua(:, :), d_ub(:, :), d_t(:, :), d_temporary(:, :)
      real(dp), allocatable :: g_j(:, :), g_ka(:, :), g_kb(:, :), g_eff_a(:, :), g_eff_b(:, :)
      real(dp), allocatable :: temporary_a(:, :), temporary_a2(:, :), temporary_b(:, :), temporary_b2(:, :)
      real(dp), allocatable :: ra(:, :), za(:, :), pa(:, :), apa(:, :), diag_a(:, :)
      real(dp), allocatable :: rb(:, :), zb(:, :), pb(:, :), apb(:, :), diag_b(:, :)
      real(dp) :: rz, rz_old, pap, alpha_cg, beta_cg, resid
      integer :: it, comp, i, a
      logical :: broke

      allocate (d_ua(nao, nao), d_ub(nao, nao), d_t(nao, nao), d_temporary(nao, nao))
      allocate (g_j(nao, nao), g_ka(nao, nao), g_kb(nao, nao), g_eff_a(nao, nao), g_eff_b(nao, nao))
      allocate (temporary_a(nao, noa), temporary_a2(nao, noa), temporary_b(nao, nob), temporary_b2(nao, nob))
      allocate (ra(nva, noa), za(nva, noa), pa(nva, noa), apa(nva, noa), diag_a(nva, noa))
      allocate (rb(nvb, nob), zb(nvb, nob), pb(nvb, nob), apb(nvb, nob), diag_b(nvb, nob))
      do i = 1, noa
         do a = 1, nva
            diag_a(a, i) = eps_a(noa + a) - eps_a(i)
         end do
      end do
      do i = 1, nob
         do a = 1, nvb
            diag_b(a, i) = eps_b(nob + a) - eps_b(i)
         end do
      end do

      !$acc enter data create(d_ua, d_ub, d_t, d_temporary, g_j, g_ka, g_kb, g_eff_a, g_eff_b, &
      !$acc                   temporary_a, temporary_a2, temporary_b, temporary_b2, ra, za, pa, apa, rb, zb, pb, apb)
      !$acc enter data copyin(diag_a, diag_b)

      allocate (res%iterations(ncomp))
      res%converged = .true.
      res%message = ""

      do comp = 1, ncomp
         call cphf_neg_matrix(nva, noa, bvec_a(:, :, comp), ra)
         call cphf_neg_matrix(nvb, nob, bvec_b(:, :, comp), rb)
         call cphf_zero_matrix(nva, noa, u_a(:, :, comp))
         call cphf_zero_matrix(nvb, nob, u_b(:, :, comp))
         call cphf_precond_matrix(nva, noa, diag_a, ra, za)
         call cphf_precond_matrix(nvb, nob, diag_b, rb, zb)
         call cphf_copy_matrix(nva, noa, za, pa)
         call cphf_copy_matrix(nvb, nob, zb, pb)
         rz = lp%dot_2d(ra, za) + lp%dot_2d(rb, zb)
         resid = sqrt(lp%dot_2d(ra, ra) + lp%dot_2d(rb, rb))
         it = 0
         broke = .false.
         ! See trc_cphf_solve_rhf's own comment on `.not. (resid < conv)`
         ! vs `resid >= conv`: a NaN resid must read as "not converged",
         ! not fall through both comparisons as false.
         do while (it < opts%max_iteration .and. .not. (resid < opts%conv_residual))
            call matvec(pa, pb, apa, apb)
            pap = lp%dot_2d(pa, apa) + lp%dot_2d(pb, apb)
            ! See trc_cphf_solve_rhf's own comment on this guard: a
            ! breakdown (pap <= 0, or NaN) must stop here, never fall
            ! through to an Inf/NaN alpha_cg that the residual check
            ! below cannot catch.
            if (.not. (pap > 0.0_dp)) then
               broke = .true.
               exit
            end if
            alpha_cg = rz/pap
            call cphf_axpy_matrix(nva, noa, alpha_cg, pa, u_a(:, :, comp))
            call cphf_axpy_matrix(nvb, nob, alpha_cg, pb, u_b(:, :, comp))
            call cphf_axpy_matrix(nva, noa, -alpha_cg, apa, ra)
            call cphf_axpy_matrix(nvb, nob, -alpha_cg, apb, rb)
            resid = sqrt(lp%dot_2d(ra, ra) + lp%dot_2d(rb, rb))
            it = it + 1
            if (resid < opts%conv_residual) exit
            call cphf_precond_matrix(nva, noa, diag_a, ra, za)
            call cphf_precond_matrix(nvb, nob, diag_b, rb, zb)
            rz_old = rz
            rz = lp%dot_2d(ra, za) + lp%dot_2d(rb, zb)
            beta_cg = rz/rz_old
            call cphf_combine_p(nva, noa, za, beta_cg, pa)
            call cphf_combine_p(nvb, nob, zb, beta_cg, pb)
         end do
         res%iterations(comp) = it
         if (broke) then
            res%converged = .false.
            res%message = "trc_cphf: CG breakdown (p^T A p <= 0, or NaN) on a right-hand side"
         else if (.not. (resid < opts%conv_residual)) then
            res%converged = .false.
            res%message = "trc_cphf: a right-hand side did not converge"
         end if
         if (opts%verbose) then
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,i0,a,es9.2)") &
                  "   cphf (uhf) rhs ", comp, ": ", it, " iterations, residual ", resid
               call global_logger%info(trim(logline))
            end block
         end if
      end do

      !$acc exit data delete(diag_a, diag_b)
      !$acc exit data delete(d_ua, d_ub, d_t, d_temporary, g_j, g_ka, g_kb, g_eff_a, g_eff_b, &
      !$acc                  temporary_a, temporary_a2, temporary_b, temporary_b2, ra, za, pa, apa, rb, zb, pb, apb)

   contains

      subroutine matvec(pina, pinb, apouta, apoutb)
         real(dp), intent(in) :: pina(nva, noa), pinb(nvb, nob)
         real(dp), intent(inout) :: apouta(nva, noa), apoutb(nvb, nob)
         ! D_U^alpha
         call lp%gemm_mat(cmo_a(:, noa + 1:nao), pina, temporary_a)
         call lp%gemm_mat(temporary_a, cmo_a(:, 1:noa), d_ua, transb="T")
         call lp%transpose(nao, nao, d_ua, nao, d_temporary, nao)
         call cphf_add_matrix(nao, nao, d_ua, d_temporary, d_ua)
         ! D_U^beta
         call lp%gemm_mat(cmo_b(:, nob + 1:nao), pinb, temporary_b)
         call lp%gemm_mat(temporary_b, cmo_b(:, 1:nob), d_ub, transb="T")
         call lp%transpose(nao, nao, d_ub, nao, d_temporary, nao)
         call cphf_add_matrix(nao, nao, d_ub, d_temporary, d_ub)
         call cphf_add_matrix(nao, nao, d_ua, d_ub, d_t)
         call eri%fock_resident(b, d_t, g_j, &
                                opts=trc_fock_options_t(j_scale=1.0_dp, k_scale=0.0_dp, density_screen=.false.))
         call eri%fock_resident(b, d_ua, g_ka, &
                                opts=trc_fock_options_t(j_scale=0.0_dp, k_scale=2.0_dp, density_screen=.false.))
         call eri%fock_resident(b, d_ub, g_kb, &
                                opts=trc_fock_options_t(j_scale=0.0_dp, k_scale=2.0_dp, density_screen=.false.))
         call cphf_add_matrix(nao, nao, g_j, g_ka, g_eff_a)
         call cphf_add_matrix(nao, nao, g_j, g_kb, g_eff_b)
         call lp%gemm_mat(g_eff_a, cmo_a(:, 1:noa), temporary_a2)
         call lp%gemm_mat(cmo_a(:, noa + 1:nao), temporary_a2, apouta, transa="T")
         call lp%gemm_mat(g_eff_b, cmo_b(:, 1:nob), temporary_b2)
         call lp%gemm_mat(cmo_b(:, nob + 1:nao), temporary_b2, apoutb, transa="T")
         call cphf_add_diag(nva, noa, diag_a, pina, apouta)
         call cphf_add_diag(nvb, nob, diag_b, pinb, apoutb)
      end subroutine matvec

   end subroutine trc_cphf_solve_uhf

   recursive subroutine trc_scf_polarizability_rhf(b, pl, eri, nocc, cmo, eps, alpha, opts, res, la)
      !!
      !! Static dipole polarizability, RHF: alpha_xy = -2 [sum_ai mu^x_ai U^y_ai
      !! + sum_ai mu^y_ai U^x_ai] (symmetrised -- see this routine's own
      !! comment), with U solving (e_a-e_i) U^y + [A U^y] = -mu^y (the "B" of
      !! this module's equation is the dipole integral itself, occ-virt MO
      !! block -- a static field adds no AO-basis response, so there is no
      !! other term).
      !! `cmo`/`eps` are `res%cmo(:,:,1)`/`res%eps(:,1)`, HOST, exactly as
      !! trc_scf_run hands them back. `pl` and `eri` are the caller's, already
      !! built on the SAME basis `b` (not necessarily the ones trc_scf_run used
      !! internally -- that engine is released when the SCF returns).
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(b%nao, b%nao), eps(b%nao)
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional :: la
      !! The caller's handles; absent, made and released for this call.

      real(dp), allocatable :: mmat(:, :, :), temporary(:, :), bvec(:, :, :), u(:, :, :)
      integer :: nao, nvirt, x, y

      nao = b%nao
      nvirt = nao - nocc
      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call la_local%init(nao)
            call trc_scf_polarizability_rhf(b, pl, eri, nocc, cmo, eps, alpha, opts, res, la_local)
            call la_local%release()
         end block
         return
      end if

      allocate (mmat(nao, nao, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)

      allocate (temporary(nao, nocc), bvec(nvirt, nocc, 3), u(nvirt, nocc, 3))
      !$acc enter data copyin(cmo, mmat) create(temporary, bvec, u)
      do x = 1, 3
         call la%gemm_mat(mmat(:, :, x), cmo(:, 1:nocc), temporary)
         call la%gemm_mat(cmo(:, nocc + 1:nao), temporary, bvec(:, :, x), transa="T")
      end do

      call trc_cphf_solve_rhf(b, eri, la, nao, nocc, nvirt, cmo, eps, 3, bvec, u, opts, res)

      ! alpha_xy = -2 [sum_ai mu^x_ai U^y_ai + sum_ai mu^y_ai U^x_ai], the
      ! symmetrised form pyscf's rhf.polarizability computes (e2 + e2^T):
      ! a not-fully-converged CG residual breaks mu^x.U^y = mu^y.U^x, so
      ! the symmetrised sum is needed to cancel that asymmetry, not just
      ! a cosmetic average.
      do y = 1, 3
         do x = 1, 3
            alpha(x, y) = -2.0_dp*(la%dot_2d(bvec(:, :, x), u(:, :, y)) + la%dot_2d(bvec(:, :, y), u(:, :, x)))
         end do
      end do

      !$acc exit data delete(cmo, mmat, temporary, bvec, u)
      deallocate (mmat, temporary, bvec, u)
   end subroutine trc_scf_polarizability_rhf

   recursive subroutine trc_scf_polarizability_uhf(b, pl, eri, noa, nob, cmo_a, cmo_b, eps_a, eps_b, alpha, opts, &
                                                   res, la)
      !!
      !! Static dipole polarizability, UHF: alpha_xy = -sum_s sum_ai mu^x_ai,s
      !! U^y_ai,s -- the sum over both spin channels already carries RHF's
      !! factor of 2 (one identical contribution per spin there; two distinct
      !! ones, each counted once, here).
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: noa, nob
      real(dp), intent(in) :: cmo_a(b%nao, b%nao), cmo_b(b%nao, b%nao), eps_a(b%nao), eps_b(b%nao)
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional :: la
      !! The caller's handles; absent, made and released for this call.

      real(dp), allocatable :: mmat(:, :, :), temporary_a(:, :), temporary_b(:, :)
      real(dp), allocatable :: bvec_a(:, :, :), bvec_b(:, :, :), u_a(:, :, :), u_b(:, :, :)
      integer :: nao, nva, nvb, x, y

      nao = b%nao
      nva = nao - noa
      nvb = nao - nob
      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_local
            call la_local%init(nao)
            call trc_scf_polarizability_uhf(b, pl, eri, noa, nob, cmo_a, cmo_b, eps_a, eps_b, alpha, opts, res, &
                                            la_local)
            call la_local%release()
         end block
         return
      end if

      allocate (mmat(nao, nao, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)

      allocate (temporary_a(nao, noa), temporary_b(nao, nob))
      allocate (bvec_a(nva, noa, 3), bvec_b(nvb, nob, 3), u_a(nva, noa, 3), u_b(nvb, nob, 3))
      !$acc enter data copyin(cmo_a, cmo_b, mmat) create(temporary_a, temporary_b, bvec_a, bvec_b, u_a, u_b)
      do x = 1, 3
         call la%gemm_mat(mmat(:, :, x), cmo_a(:, 1:noa), temporary_a)
         call la%gemm_mat(cmo_a(:, noa + 1:nao), temporary_a, bvec_a(:, :, x), transa="T")
         call la%gemm_mat(mmat(:, :, x), cmo_b(:, 1:nob), temporary_b)
         call la%gemm_mat(cmo_b(:, nob + 1:nao), temporary_b, bvec_b(:, :, x), transa="T")
      end do

      call trc_cphf_solve_uhf(b, eri, la, nao, noa, nva, nob, nvb, cmo_a, cmo_b, eps_a, eps_b, &
                              3, bvec_a, bvec_b, u_a, u_b, opts, res)

      ! alpha_xy = -[sum_s sum_ai mu^x_ai,s U^y_ai,s + (x<->y)] -- pyscf's
      ! uhf.polarizability e2 + e2^T; see trc_scf_polarizability_rhf's own
      ! comment on why the symmetrised form, not the bare one, is right.
      do y = 1, 3
         do x = 1, 3
            alpha(x, y) = -(la%dot_2d(bvec_a(:, :, x), u_a(:, :, y)) + la%dot_2d(bvec_b(:, :, x), u_b(:, :, y)) + &
                            la%dot_2d(bvec_a(:, :, y), u_a(:, :, x)) + la%dot_2d(bvec_b(:, :, y), u_b(:, :, x)))
         end do
      end do

      !$acc exit data delete(cmo_a, cmo_b, mmat, temporary_a, temporary_b, bvec_a, bvec_b, u_a, u_b)
      deallocate (mmat, temporary_a, temporary_b, bvec_a, bvec_b, u_a, u_b)
   end subroutine trc_scf_polarizability_uhf

   subroutine trc_scf_polarizability(b, pl, eri, nalpha, nbeta, res_scf, alpha, opts, res, la)
      !!
      !! Dispatch on `res%nspin`: the one entry point a caller holding a
      !! converged `trc_scf_result_t` needs. `nalpha`/`nbeta` are the same
      !! electron counts the SCF run itself took -- trc_scf_result_t does not
      !! carry them (trc_scf_driver.F90's own `cfg%nocc` is private to that
      !! run), so the caller passes them again, exactly as scf_rimp2.F90 does
      !! for `trc_rimp2_run`'s `nocc`.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_result_t), intent(in) :: res_scf
      real(dp), intent(out) :: alpha(3, 3)
      type(trc_cphf_options_t), intent(in) :: opts
      type(trc_cphf_result_t), intent(out) :: res
      type(trc_linalg_t), intent(inout), optional :: la

      if (res_scf%nspin == 1) then
         call trc_scf_polarizability_rhf(b, pl, eri, nalpha, res_scf%cmo(:, :, 1), res_scf%eps(:, 1), &
                                         alpha, opts, res, la)
      else
         call trc_scf_polarizability_uhf(b, pl, eri, nalpha, nbeta, res_scf%cmo(:, :, 1), res_scf%cmo(:, :, 2), &
                                         res_scf%eps(:, 1), res_scf%eps(:, 2), alpha, opts, res, la)
      end if
   end subroutine trc_scf_polarizability

   subroutine cphf_zero_matrix(m, n, a)
      !! ----------------------------------------------------------------------
      !! Explicit-shape elementwise helpers shared by the RHF/UHF CG loops
      !! above: never touching anything but their own dummies, so each is a
      !! legal `do concurrent` body on its own -- see trc_scf_fock.F90's own
      !! header for why that split matters on the GPU build.
      !! ----------------------------------------------------------------------
      integer, intent(in) :: m, n
      real(dp), intent(inout) :: a(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         a(i, j) = 0.0_dp
      end do
   end subroutine cphf_zero_matrix

   subroutine cphf_copy_matrix(m, n, a, c)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = a(i, j)
      end do
   end subroutine cphf_copy_matrix

   subroutine cphf_neg_matrix(m, n, a, c)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = -a(i, j)
      end do
   end subroutine cphf_neg_matrix

   subroutine cphf_add_matrix(m, n, a, bb, c)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: a(m, n), bb(m, n)
      real(dp), intent(inout) :: c(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         c(i, j) = a(i, j) + bb(i, j)
      end do
   end subroutine cphf_add_matrix

   subroutine cphf_axpy_matrix(m, n, alpha, x, y)
      !! y := y + alpha*x
      integer, intent(in) :: m, n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: x(m, n)
      real(dp), intent(inout) :: y(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         y(i, j) = y(i, j) + alpha*x(i, j)
      end do
   end subroutine cphf_axpy_matrix

   subroutine cphf_precond_matrix(m, n, diagv, r, z)
      !! z := r / diagv, the diagonal (e_a-e_i) preconditioner.
      integer, intent(in) :: m, n
      real(dp), intent(in) :: diagv(m, n), r(m, n)
      real(dp), intent(inout) :: z(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         z(i, j) = r(i, j)/diagv(i, j)
      end do
   end subroutine cphf_precond_matrix

   subroutine cphf_combine_p(m, n, z, beta, p)
      !! p := z + beta*p
      integer, intent(in) :: m, n
      real(dp), intent(in) :: beta
      real(dp), intent(in) :: z(m, n)
      real(dp), intent(inout) :: p(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         p(i, j) = z(i, j) + beta*p(i, j)
      end do
   end subroutine cphf_combine_p

   subroutine cphf_add_diag(m, n, diagv, pin, apout)
      !! apout := apout + diagv * pin (elementwise)
      integer, intent(in) :: m, n
      real(dp), intent(in) :: diagv(m, n), pin(m, n)
      real(dp), intent(inout) :: apout(m, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:m)
         apout(i, j) = apout(i, j) + diagv(i, j)*pin(i, j)
      end do
   end subroutine cphf_add_diag

end module trc_cphf
