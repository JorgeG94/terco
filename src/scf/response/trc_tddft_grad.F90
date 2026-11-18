module trc_tddft_grad
   !!
   !! Analytic nuclear gradients of TDA (CIS for Hartree-Fock) singlet excited
   !! states of a closed-shell RHF/RKS reference, for several roots at once.
   !! The method is the Z-vector (Furche-Ahlrichs) gradient, term by term that
   !! of pyscf's grad/tdrks.py `grad_elec` (singlet, Y = 0), and this module
   !! cites it where a line follows it.
   !!
   !! NORMALISATION. terco's X (trc_tddft) has unit norm, pyscf's x half of it
   !! (sum x^2 = 1/2). With T = C_occ X C_vir^T (the AO transition density,
   !! NON-symmetric) and R = T^T/sqrt(2) = pyscf's `dmxpy` = `dmxmy`:
   !!
   !!     S = R + R^T = (T + T^T)/sqrt(2)      pyscf dm[1], symmetric
   !!     A = R - R^T = (T^T - T)/sqrt(2)      pyscf dm[2], antisymmetric
   !!     P_d = C_occ (-X X^T) C_occ^T + C_vir (X^T X) C_vir^T
   !!                                          pyscf `dmzoo`, the unrelaxed
   !!                                          difference density
   !!
   !! and the operator G[D] = 2 J[D] - c_x K[D] (+ 2 f_xc[D]) is the one the
   !! excitation-energy matvec uses.
   !!
   !! THE STEPS (one pass over all roots each)
   !!  1. G[P_d], G[S] and G[A] for every root through ONE
   !!     `eri%fock_nosym_many` call; f_xc[P_d] and f_xc[S] through one
   !!     `trc_xc_fxc_rks_apply_many`; the k_xc term k1ao = E_xc'''[S, S, .]
   !!     through one `trc_xc_kxc_rks_apply`. Together with the MO blocks of
   !!     G[S] and G[A] these are the right-hand sides w_ai of the Z-vector
   !!     equation (pyscf `wvo`).
   !!  2. ONE CPKS (CPHF for Hartree-Fock) solve with ncomp = nroots, for the
   !!     Z vectors of all roots. trc_cpks_solve_rks and trc_cphf_solve_rhf run
   !!     the right-hand sides one after another, each CG with its own
   !!     integral passes: only the right-hand sides' setup is batched.
   !!  3. The relaxed difference density P = P_d + sym(Z) and the
   !!     energy-weighted density W of the excitation (pyscf `im0 + zeta dm1`,
   !!     without its ground-state part, which is the ground state's W).
   !!  4. The gradient of root k is the ground-state gradient plus
   !!       h^x . P - S^x . W                          (trc_1e_gradient)
   !!       + d/dR [ P.G_hf[C C^T] + 1/2 S.G_hf[S] - c_x/2 A.K[A] ]
   !!                                                  (eri%gradient_pairs)
   !!       + the XC terms (trc_xc_grad_response_many with P and 2 sym(T)),
   !!     where the two-electron derivative of every root and the ground state
   !!     come from the same sweep of the integral derivatives.
   !!
   !! THE TWO-ELECTRON PAIRS. pyscf contracts the bra derivative integrals with
   !! the ground density (veff1[0] with P, veff1[1] with C C^T), with S
   !! (veff1[2], twice) and with the antisymmetric A (veff1[3]). Those are the
   !! four-centre derivatives of three bilinear forms at fixed densities, with
   !! B(U, V) = U . (2 J - c_x K)[V]:
   !!     B(P, o) with o = C_occ C_occ^T       the relaxed density against the
   !!                                          ground one; a cross pair counts
   !!                                          both orderings, so (coulomb 1,
   !!                                          exchange c_x/2)
   !!     B(S, S)/2                            (S, S, 1, c_x/2)
   !!     -c_x A.K[A]/2                        (A, A, 0, -c_x/2): the exchange of a
   !!                                          non-symmetric X is not that of its
   !!                                          symmetric part, and A is where the
   !!                                          rest lives. The sign: an
   !!                                          antisymmetric pair's exchange is
   !!                                          taken with the transposed density,
   !!                                          K[A^T] = -K[A] (gradient_pairs)
   !! and the ground state's energy is the pair (o, o, 2, c_x). Each root is a
   !! column of the sweep; `roots_per_sweep` (default 2, at most 4) caps how many roots share one pass.
   !!
   !! GRID. The XC terms are Pulay only unless `grid_response` is set (pyscf's
   !! tdrks has no grid response either: compare with it off). The default is
   !! therefore NOT the exact derivative of terco's energy with a grid that
   !! moves with the nuclei (the SCF's grids do): it misses the grid's
   !! point-motion and partition-weight terms, which check_tddft_grad's
   !! finite-difference case measures. Set `grid_response` for the exact
   !! gradient; the grid must then come from `build_from_grid` with `partition`
   !! set to its own, exactly as trc_xc_grad_response_many requires.
   !!
   !! SCOPE. Restricted, singlet, TDA, LDA/GGA and hybrid functionals with an
   !! exchange fraction c_x, Cartesian basis, no density fitting; the refusals are
   !! trc_tddft_refuse_reason's, and a spherical basis. A degenerate root has no
   !! unique gradient (any rotation inside the multiplet is an equally good X);
   !! a warning is logged when two computed roots lie within 1e-6 hartree. The
   !! dense algebra (AO/MO transforms, the Z-vector's right-hand side) runs on the
   !! host.
   !!
   use pic_types, only: dp, int32, int64
   use pic_blas_interfaces, only: pic_gemm
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_1e_gradient
   use trc_eri, only: trc_eri_t, trc_fock_options_t, trc_grad_pair_t
   use trc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc, only: trc_xc_fxc_cache_t, trc_xc_fxc_rks_apply_many, trc_xc_kxc_rks_apply
   use trc_xc_grad, only: xc_grad_add_on_grid, xc_grad_options, xc_grad_step_names, xc_grad_options_t
   use trc_xc_grad_response, only: trc_xc_grad_response_many, trc_xc_grad_response_options_t
   use trc_dft_partition, only: dft_partition_options_t
   use trc_step_table, only: step_table_t
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_cphf_solve_rhf
   use trc_cpks, only: trc_cpks_solve_rks
   use trc_scf_grad, only: scf_grad_nuclear_repulsion
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_reference_t, trc_tddft_refuse_reason
   implicit none
   private

   public :: trc_tddft_grad_options_t, trc_tddft_gradient_rks, trc_tddft_grad_refuse_reason
   ! The building blocks the RPA gradient (trc_tddft_rpa_grad) shares with this module.
   public :: tddft_grad_state_t, tddft_grad_gemm, tddft_grad_basis_change, tddft_grad_mo_block, tddft_grad_symmetrize, &
             tddft_grad_fock, tddft_grad_fxc, tddft_grad_zvector, tddft_grad_chunk, tddft_grad_xc, &
             tddft_grad_one_electron
   public :: SQRT_HALF

   integer, parameter :: MAX_ROOTS_PER_SWEEP = 4
   integer, parameter :: DEFAULT_ROOTS_PER_SWEEP = 2
   real(dp), parameter :: DEGENERACY_TOL = 1.0e-6_dp
   real(dp), parameter :: SQRT_HALF = 0.70710678118654752440_dp

   type :: trc_tddft_grad_options_t
      !! The settings of `trc_tddft_gradient_rks`, and what it returns.
      integer :: roots_per_sweep = 2
      !! Roots sharing one pass over the integral derivatives (and one pool of
      !! densities in memory), clamped to 1..MAX_ROOTS_PER_SWEEP; zero or less
      !! takes the default 2. Each extra column of a pass costs about one
      !! single-column sweep (the fold is per column and the wider kernels
      !! spill), so batching saves only integral setup and screening and the
      !! default is small. A tuning knob, not a correctness one: the gradients
      !! do not depend on it.
      real(dp) :: cpks_tol = 1.0e-9_dp
      !! The Z-vector's residual norm convergence (trc_cphf_options_t%conv_residual).
      integer :: cpks_max_iteration = 100
      !! CG iterations per root.
      logical :: grid_response = .false.
      !! Differentiate the grid's points and weights too (the XC terms'
      !! `grid_response`); needs `partition`. Off matches pyscf's tdrks.
      type(dft_partition_options_t) :: partition
      !! The partition `xg`'s grid was built with, read only with `grid_response`.
      logical :: verbose = .false.
      real(dp), allocatable :: grad(:, :, :)
      !! Out: (3, natm, nroots), hartree/bohr, the full gradient of every root's
      !! total energy (nuclear repulsion and ground state included).
      type(error_t) :: err
      !! Out: a refusal, or a Z-vector or sweep failure; `grad` is then not allocated.
      integer :: nsweep = 0
      !! Out: the integral-derivative sweeps made (one per chunk of roots).
      integer, allocatable :: cpks_iterations(:)
      !! Out: CG iterations each root's Z vector took.
   end type trc_tddft_grad_options_t

   type :: tddft_grad_state_t
      !! Everything the steps hand each other. AO matrices are (nao, nao, nroots).
      integer :: n = 0, nocc = 0, nvir = 0, nroots = 0
      real(dp) :: cx = 1.0_dp
      logical :: dft = .false.
      logical :: has_budget = .false.
      integer(int64) :: budget = 0
      real(dp), allocatable :: cocc(:, :), cvir(:, :), eps(:)
      real(dp), allocatable :: xamp(:, :, :)     !! (nocc, nvir, nroots) unit-norm amplitudes
      real(dp), allocatable :: tmat(:, :, :)     !! T = C_occ X C_vir^T
      real(dp), allocatable :: pdiff(:, :, :)    !! P_d
      real(dp), allocatable :: veff0(:, :, :)    !! pyscf veff0doo = G[P_d] + 2 f_xc[P_d] + 2 k1ao
      real(dp), allocatable :: voo(:, :, :), vvo(:, :, :), vvv(:, :, :)
      !! MO blocks of G[S] + 2 f_xc[S] (pyscf veff0mop): (nocc,nocc), (nvir,nocc), (nvir,nvir).
      real(dp), allocatable :: aoo(:, :, :), avo(:, :, :), avv(:, :, :)
      !! MO blocks of G[A] (pyscf veff0mom).
      real(dp), allocatable :: wvo(:, :, :), zvec(:, :, :)
      !! (nvir, nocc, nroots): the Z-vector's right-hand side and solution.
      real(dp), allocatable :: prel(:, :, :)     !! sym(P_d + Z)
      real(dp), allocatable :: wexc(:, :, :)     !! the excitation's energy-weighted density, symmetric
   end type tddft_grad_state_t

contains

   pure function trc_tddft_grad_refuse_reason(ref, td, nao_cart) result(msg)
      !! Why this gradient cannot run, or an empty string: what trc_tddft
      !! refuses for the functional, an unconverged or empty result, a basis
      !! that is not Cartesian-sized, or orbitals that are not (nao, nao).
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_result_t), intent(in) :: td
      integer, intent(in) :: nao_cart
      character(len=256) :: msg
      type(trc_tddft_options_t) :: probe
      msg = ""
      if (ref%dft) then
         msg = trc_tddft_refuse_reason(probe, ref%func)
      else
         msg = trc_tddft_refuse_reason(probe)
      end if
      if (len_trim(msg) > 0) return
      if (.not. td%converged .or. .not. allocated(td%xamp)) then
         msg = "trc_tddft_grad: needs the amplitudes of a converged TDA solve"
      else if (size(ref%cmo, 1) /= nao_cart .or. size(ref%cmo, 2) /= nao_cart) then
         msg = "trc_tddft_grad: the orbitals must be (nao, nao) in the Cartesian basis the integrals use " &
               //"(spherical bases and linearly dependent sets are not implemented)"
      else if (size(td%xamp, 1) /= ref%nocc .or. size(td%xamp, 2) /= nao_cart - ref%nocc) then
         msg = "trc_tddft_grad: the amplitudes do not match the reference's occupied and virtual counts"
      end if
   end function trc_tddft_grad_refuse_reason

   subroutine tddft_grad_gemm(a, b, c, transa, transb, alpha, beta)
      !! c = alpha op(a) op(b) + beta c on the host.
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)
      character(len=1), intent(in), optional :: transa, transb
      real(dp), intent(in), optional :: alpha, beta
      character(len=1) :: ta, tb
      real(dp) :: al, be
      ta = "N"
      tb = "N"
      al = 1.0_dp
      be = 0.0_dp
      if (present(transa)) ta = transa
      if (present(transb)) tb = transb
      if (present(alpha)) al = alpha
      if (present(beta)) be = beta
      call pic_gemm(a, b, c, ta, tb, al, be)
   end subroutine tddft_grad_gemm

   subroutine tddft_grad_basis_change(left, dmo, right, out)
      !! out = left dmo right^T, the AO matrix of an MO-basis block.
      real(dp), intent(in) :: left(:, :), dmo(:, :), right(:, :)
      real(dp), intent(out) :: out(:, :)
      real(dp), allocatable :: half(:, :)
      allocate (half(size(left, 1), size(dmo, 2)))
      call tddft_grad_gemm(left, dmo, half)
      call tddft_grad_gemm(half, right, out, transb="T")
   end subroutine tddft_grad_basis_change

   subroutine tddft_grad_mo_block(left, ao, right, out)
      !! out = left^T ao right, an MO block of an AO matrix.
      real(dp), intent(in) :: left(:, :), ao(:, :), right(:, :)
      real(dp), intent(out) :: out(:, :)
      real(dp), allocatable :: half(:, :)
      allocate (half(size(ao, 1), size(right, 2)))
      call tddft_grad_gemm(ao, right, half)
      call tddft_grad_gemm(left, half, out, transa="T")
   end subroutine tddft_grad_mo_block

   subroutine tddft_grad_symmetrize(a)
      real(dp), intent(inout) :: a(:, :)
      a = 0.5_dp*(a + transpose(a))
   end subroutine tddft_grad_symmetrize

   subroutine tddft_grad_init(st, ref, td)
      !! The orbitals, T and the unrelaxed difference density of every root.
      type(tddft_grad_state_t), intent(out) :: st
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_result_t), intent(in) :: td
      real(dp), allocatable :: doo(:, :), dvv(:, :), half(:, :), xm(:, :)
      integer :: k

      st%n = size(ref%cmo, 1)
      st%nocc = ref%nocc
      st%nvir = st%n - ref%nocc
      st%nroots = size(td%xamp, 3)
      st%dft = ref%dft
      st%cx = 1.0_dp
      if (ref%dft) st%cx = ref%func%exx
      st%has_budget = ref%has_budget
      st%budget = ref%budget
      st%cocc = ref%cmo(:, 1:st%nocc)
      st%cvir = ref%cmo(:, st%nocc + 1:st%n)
      st%eps = ref%eps(1:st%n)
      st%xamp = td%xamp
      allocate (st%tmat(st%n, st%n, st%nroots), st%pdiff(st%n, st%n, st%nroots))
      allocate (doo(st%nocc, st%nocc), dvv(st%nvir, st%nvir), half(st%n, st%n))
      do k = 1, st%nroots
         xm = st%xamp(:, :, k)
         call tddft_grad_basis_change(st%cocc, xm, st%cvir, st%tmat(:, :, k))
         ! P_d: doo = -X X^T, dvv = X^T X (pyscf doo, dvv with x = X/sqrt(2)).
         call tddft_grad_gemm(xm, xm, doo, transb="T", alpha=-1.0_dp)
         call tddft_grad_gemm(xm, xm, dvv, transa="T")
         call tddft_grad_basis_change(st%cocc, doo, st%cocc, st%pdiff(:, :, k))
         call tddft_grad_basis_change(st%cvir, dvv, st%cvir, half)
         st%pdiff(:, :, k) = st%pdiff(:, :, k) + half
         call tddft_grad_symmetrize(st%pdiff(:, :, k))
      end do
   end subroutine tddft_grad_init

   subroutine tddft_grad_fock(b, eri, ndens, dmats, gmats, jscale, kscale)
      !! gmats = jscale J - kscale K/2 of every (possibly non-symmetric)
      !! density, one pass.
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      integer, intent(in) :: ndens
      real(dp), intent(in) :: dmats(:, :, :)
      real(dp), intent(out) :: gmats(:, :, :)
      real(dp), intent(in) :: jscale, kscale
      call eri%fock_nosym_many(b, ndens, dmats, gmats, &
                               trc_fock_options_t(j_scale=jscale, k_scale=kscale, density_screen=.false.))
   end subroutine tddft_grad_fock

   subroutine tddft_grad_fxc(st, b, xg, cache, la, dmats, vresp)
      !! f_xc[D_k] of every density, one call; `dmats` and `vresp` are (nao, nao, nvec).
      type(tddft_grad_state_t), intent(in) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_xc_fxc_cache_t), intent(in) :: cache
      type(trc_linalg_t), intent(inout) :: la
      real(dp), intent(in) :: dmats(:, :, :)
      real(dp), intent(out) :: vresp(:, :, :)
      integer :: nvec
      nvec = size(dmats, 3)
      !$acc enter data copyin(dmats) create(vresp)
      if (st%has_budget) then
         call trc_xc_fxc_rks_apply_many(b, xg, cache, nvec, dmats, vresp, st%budget, la)
      else
         call trc_xc_fxc_rks_apply_many(b, xg, cache, nvec, dmats, vresp, la=la)
      end if
      !$acc update self(vresp)
      !$acc exit data delete(dmats, vresp)
   end subroutine tddft_grad_fxc

   subroutine tddft_grad_rhs(st, b, eri, la, ref, xg, err)
      !! Steps 1: the right-hand side of the Z-vector equation for every root
      !! (pyscf `wvo`), with the MO blocks it and the energy-weighted
      !! density need. `err` is set on a refused kernel.
      type(tddft_grad_state_t), intent(inout) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in), optional :: xg
      type(error_t), intent(out) :: err

      type(trc_xc_fxc_cache_t) :: cache
      real(dp), allocatable :: dm3(:, :, :), g3(:, :, :), dmx(:, :, :), fx(:, :, :), dmk(:, :, :), vk(:, :, :)
      real(dp), allocatable :: veff_s(:, :), veff_a(:, :), xpy(:, :), half(:, :)
      integer :: n, nocc, nvir, nr, k
      n = st%n
      nocc = st%nocc
      nvir = st%nvir
      nr = st%nroots
      allocate (dm3(n, n, 3*nr), g3(n, n, 3*nr), veff_s(n, n), veff_a(n, n), xpy(nvir, nocc))
      allocate (st%veff0(n, n, nr), st%voo(nocc, nocc, nr), st%vvo(nvir, nocc, nr), st%vvv(nvir, nvir, nr))
      allocate (st%aoo(nocc, nocc, nr), st%avo(nvir, nocc, nr), st%avv(nvir, nvir, nr), st%wvo(nvir, nocc, nr))
      do k = 1, nr
         dm3(:, :, 3*k - 2) = st%pdiff(:, :, k)
         dm3(:, :, 3*k - 1) = SQRT_HALF*(st%tmat(:, :, k) + transpose(st%tmat(:, :, k)))
         dm3(:, :, 3*k) = SQRT_HALF*(transpose(st%tmat(:, :, k)) - st%tmat(:, :, k))
      end do
      ! 2 J[D] - c_x K[D] of P_d, S and A (pyscf vj*2 - vk*hyb): A's J is zero,
      ! its K is the antisymmetric veff0mom.
      call tddft_grad_fock(b, eri, 3*nr, dm3, g3, 2.0_dp, 2.0_dp*st%cx)

      allocate (fx(n, n, 2*nr), dmx(n, n, 2*nr), dmk(n, n, 1 + 2*nr), vk(n, n, nr))
      fx = 0.0_dp
      vk = 0.0_dp
      if (st%dft) then
         call cache%setup(b, xg, ref%func, ref%dmat, la)
         do k = 1, nr
            dmx(:, :, 2*k - 1) = dm3(:, :, 3*k - 2)
            dmx(:, :, 2*k) = dm3(:, :, 3*k - 1)
         end do
         call tddft_grad_fxc(st, b, xg, cache, la, dmx, fx)
         call cache%release()
         ! k1ao = E_xc'''[S, S, .] (pyscf k1ao[0]); packed as D, A_k = S_k, B_k = S_k.
         dmk(:, :, 1) = ref%dmat
         do k = 1, nr
            dmk(:, :, 1 + k) = dm3(:, :, 3*k - 1)
            dmk(:, :, 1 + nr + k) = dm3(:, :, 3*k - 1)
         end do
         if (st%has_budget) then
            call trc_xc_kxc_rks_apply(b, xg, ref%func, dmk, vk, err, st%budget, la)
         else
            call trc_xc_kxc_rks_apply(b, xg, ref%func, dmk, vk, err, la=la)
         end if
         if (err%has_error()) return
      end if

      allocate (half(nvir, n))
      do k = 1, nr
         st%veff0(:, :, k) = g3(:, :, 3*k - 2) + 2.0_dp*fx(:, :, 2*k - 1) + 2.0_dp*vk(:, :, k)
         veff_s = g3(:, :, 3*k - 1) + 2.0_dp*fx(:, :, 2*k)
         veff_a = g3(:, :, 3*k)
         call tddft_grad_mo_block(st%cocc, veff_s, st%cocc, st%voo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_s, st%cocc, st%vvo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_s, st%cvir, st%vvv(:, :, k))
         call tddft_grad_mo_block(st%cocc, veff_a, st%cocc, st%aoo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_a, st%cocc, st%avo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_a, st%cvir, st%avv(:, :, k))
         ! x = xpy = xmy = X^T/sqrt(2), (nvir, nocc).
         xpy = SQRT_HALF*transpose(st%xamp(:, :, k))
         ! wvo = 2 C_v^T veff0doo C_o
         call tddft_grad_mo_block(st%cvir, st%veff0(:, :, k), st%cocc, st%wvo(:, :, k))
         st%wvo(:, :, k) = 2.0_dp*st%wvo(:, :, k)
         ! - 2 xpy Voo^T + 2 Vvv^T xpy   (and the same with A's blocks: xmy = xpy)
         call tddft_grad_gemm(xpy, st%voo(:, :, k), st%wvo(:, :, k), transb="T", alpha=-2.0_dp, beta=1.0_dp)
         call tddft_grad_gemm(st%vvv(:, :, k), xpy, st%wvo(:, :, k), transa="T", alpha=2.0_dp, beta=1.0_dp)
         call tddft_grad_gemm(xpy, st%aoo(:, :, k), st%wvo(:, :, k), transb="T", alpha=-2.0_dp, beta=1.0_dp)
         call tddft_grad_gemm(st%avv(:, :, k), xpy, st%wvo(:, :, k), transa="T", alpha=2.0_dp, beta=1.0_dp)
      end do
   end subroutine tddft_grad_rhs

   subroutine tddft_grad_zvector(st, b, eri, la, ref, xg, opts)
      !! Step 2: one CPKS (CPHF) solve with ncomp = nroots for the Z vectors,
      !! (e_a - e_i) z + A z = -wvo (pyscf's cphf.solve, same sign).
      type(tddft_grad_state_t), intent(inout) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in), optional :: xg
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      type(trc_cphf_options_t) :: copts
      type(trc_cphf_result_t) :: cres
      real(dp), allocatable :: cmo(:, :), bvec(:, :, :), u(:, :, :)

      copts%conv_residual = opts%cpks_tol
      copts%max_iteration = opts%cpks_max_iteration
      copts%verbose = opts%verbose
      cmo = ref%cmo
      bvec = st%wvo
      allocate (u(st%nvir, st%nocc, st%nroots))
      u = 0.0_dp
      ! The solvers want the orbitals, the right-hand sides and the solution
      ! resident; eps and the density are read from the host.
      !$acc enter data copyin(cmo, bvec, u)
      if (st%dft) then
         if (st%has_budget) then
            call trc_cpks_solve_rks(b, eri, la, xg, ref%func, st%n, st%nocc, st%nvir, cmo, st%eps, ref%dmat, &
                                    st%nroots, bvec, u, copts, cres, st%budget)
         else
            call trc_cpks_solve_rks(b, eri, la, xg, ref%func, st%n, st%nocc, st%nvir, cmo, st%eps, ref%dmat, &
                                    st%nroots, bvec, u, copts, cres)
         end if
      else
         call trc_cphf_solve_rhf(b, eri, la, st%n, st%nocc, st%nvir, cmo, st%eps, st%nroots, bvec, u, copts, cres)
      end if
      !$acc update self(u)
      !$acc exit data delete(cmo, bvec, u)
      st%zvec = u
      opts%cpks_iterations = cres%iterations
      if (.not. cres%converged) call opts%err%set(ERROR_GENERIC, "trc_tddft_grad: the Z-vector solve failed: " &
                                                  //trim(cres%message))
   end subroutine tddft_grad_zvector

   subroutine tddft_grad_relaxed(st, b, eri, la, ref, xg, err)
      !! Step 3: P = P_d + sym(Z) and the excitation's energy-weighted
      !! density (pyscf im0 + zeta*dm1 without the ground state's 2 I).
      type(tddft_grad_state_t), intent(inout) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in), optional :: xg
      type(error_t), intent(out) :: err

      type(trc_xc_fxc_cache_t) :: cache
      real(dp), allocatable :: dmz(:, :, :), gz(:, :, :), fz(:, :, :), xpy(:, :), xk(:, :)
      real(dp), allocatable :: woo(:, :), wvv(:, :), wvo(:, :), doo(:, :), dvv(:, :), ao(:, :), veffz(:, :)
      integer :: n, nocc, nvir, nr, k, i, j
      n = st%n
      nocc = st%nocc
      nvir = st%nvir
      nr = st%nroots
      allocate (dmz(n, n, nr), gz(n, n, nr), fz(n, n, nr), veffz(n, n), ao(n, n))
      allocate (woo(nocc, nocc), wvv(nvir, nvir), wvo(nvir, nocc), doo(nocc, nocc), dvv(nvir, nvir), xpy(nvir, nocc))
      allocate (st%prel(n, n, nr), st%wexc(n, n, nr))
      do k = 1, nr
         call tddft_grad_basis_change(st%cvir, st%zvec(:, :, k), st%cocc, ao)
         dmz(:, :, k) = ao + transpose(ao)
      end do
      ! veff = vresp(Z + Z^T) = (J - c_x K/2 + f_xc)[Z + Z^T]
      call tddft_grad_fock(b, eri, nr, dmz, gz, 1.0_dp, st%cx)
      fz = 0.0_dp
      if (st%dft) then
         call cache%setup(b, xg, ref%func, ref%dmat, la)
         call tddft_grad_fxc(st, b, xg, cache, la, dmz, fz)
         call cache%release()
      end if
      do k = 1, nr
         veffz = st%veff0(:, :, k) + gz(:, :, k) + fz(:, :, k)
         xk = st%xamp(:, :, k)
         xpy = SQRT_HALF*transpose(xk)
         ! oo: C_o^T (veff0doo + veff) C_o + V_vo^T x + A_vo^T x
         call tddft_grad_mo_block(st%cocc, veffz, st%cocc, woo)
         call tddft_grad_gemm(st%vvo(:, :, k), xpy, woo, transa="T", beta=1.0_dp)
         call tddft_grad_gemm(st%avo(:, :, k), xpy, woo, transa="T", beta=1.0_dp)
         ! vv: x V_vo^T + x A_vo^T
         call tddft_grad_gemm(xpy, st%vvo(:, :, k), wvv, transb="T")
         call tddft_grad_gemm(xpy, st%avo(:, :, k), wvv, transb="T", beta=1.0_dp)
         ! vo: 2 x V_oo^T + 2 x A_oo^T
         call tddft_grad_gemm(xpy, st%voo(:, :, k), wvo, transb="T", alpha=2.0_dp)
         call tddft_grad_gemm(xpy, st%aoo(:, :, k), wvo, transb="T", alpha=2.0_dp, beta=1.0_dp)
         ! + zeta * dm1: dm1 = (-X X^T, X^T X, Z), zeta_pq = (e_p + e_q)/2, zeta_ai = e_i
         call tddft_grad_gemm(xk, xk, doo, transb="T", alpha=-1.0_dp)
         call tddft_grad_gemm(xk, xk, dvv, transa="T")
         do j = 1, nocc
            do i = 1, nocc
               woo(i, j) = woo(i, j) + 0.5_dp*(st%eps(i) + st%eps(j))*doo(i, j)
            end do
         end do
         do j = 1, nvir
            do i = 1, nvir
               wvv(i, j) = wvv(i, j) + 0.5_dp*(st%eps(nocc + i) + st%eps(nocc + j))*dvv(i, j)
            end do
         end do
         do j = 1, nocc
            do i = 1, nvir
               wvo(i, j) = wvo(i, j) + st%eps(j)*st%zvec(i, j, k)
            end do
         end do
         call tddft_grad_basis_change(st%cocc, woo, st%cocc, st%wexc(:, :, k))
         call tddft_grad_basis_change(st%cvir, wvv, st%cvir, ao)
         st%wexc(:, :, k) = st%wexc(:, :, k) + ao
         call tddft_grad_basis_change(st%cvir, wvo, st%cocc, ao)
         st%wexc(:, :, k) = st%wexc(:, :, k) + ao
         call tddft_grad_symmetrize(st%wexc(:, :, k))
         st%prel(:, :, k) = st%pdiff(:, :, k) + 0.5_dp*dmz(:, :, k)
      end do
      err = error_t()
   end subroutine tddft_grad_relaxed

   subroutine tddft_grad_chunk(st, opts, c0, nc)
      !! The roots c0 .. c0 + nc - 1 of the sweep that starts at root c0.
      type(tddft_grad_state_t), intent(in) :: st
      type(trc_tddft_grad_options_t), intent(in) :: opts
      integer, intent(in) :: c0
      integer, intent(out) :: nc
      integer :: size_chunk
      size_chunk = DEFAULT_ROOTS_PER_SWEEP
      if (opts%roots_per_sweep > 0) size_chunk = min(opts%roots_per_sweep, MAX_ROOTS_PER_SWEEP)
      size_chunk = min(size_chunk, st%nroots)
      nc = min(size_chunk, st%nroots - c0 + 1)
   end subroutine tddft_grad_chunk

   subroutine tddft_grad_two_electron(st, b, eri, opts, gground)
      !! The two-electron derivative terms of every root, and of the ground
      !! state (`gground`), from one `gradient_pairs` sweep per chunk of roots.
      type(tddft_grad_state_t), intent(in) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      real(dp), intent(out) :: gground(3, b%natm)

      type(trc_grad_pair_t), allocatable :: pairs(:)
      real(dp), allocatable :: pool(:, :, :), gp(:, :, :), occ0(:, :)
      integer :: c0, nc, ncol, j, k, np, ip
      logical :: first

      allocate (occ0(st%n, st%n))
      call tddft_grad_gemm(st%cocc, st%cocc, occ0, transb="T")
      gground = 0.0_dp
      c0 = 1
      do while (c0 <= st%nroots)
         call tddft_grad_chunk(st, opts, c0, nc)
         first = c0 == 1
         ncol = nc + merge(1, 0, first)
         allocate (pool(st%n, st%n, 1 + 3*nc), pairs(3*nc + merge(1, 0, first)), gp(3, b%natm, ncol))
         pool(:, :, 1) = occ0
         np = 0
         do j = 1, nc
            k = c0 + j - 1
            pool(:, :, 3*j - 1) = st%prel(:, :, k)
            pool(:, :, 3*j) = SQRT_HALF*(st%tmat(:, :, k) + transpose(st%tmat(:, :, k)))
            pool(:, :, 3*j + 1) = SQRT_HALF*(transpose(st%tmat(:, :, k)) - st%tmat(:, :, k))
            ! P against the ground density (cross: both orderings), S and A on their own.
            pairs(np + 1) = trc_grad_pair_t(2_int32 + 3_int32*int(j - 1, int32), 1_int32, int(j, int32), 1.0_dp, 0.5_dp*st%cx)
            pairs(np + 2) = trc_grad_pair_t(3_int32 + 3_int32*int(j - 1, int32), 3_int32 + 3_int32*int(j - 1, int32), &
                                            int(j, int32), 1.0_dp, 0.5_dp*st%cx)
            ! gradient_pairs differentiates the exchange of an antisymmetric
            ! pair with the transposed right density, K[A^T] = -K[A], so the
            ! literal -c_x/2 A.K[A] is the weight +c_x/2 turned negative.
            pairs(np + 3) = trc_grad_pair_t(4_int32 + 3_int32*int(j - 1, int32), 4_int32 + 3_int32*int(j - 1, int32), &
                                            int(j, int32), 0.0_dp, -0.5_dp*st%cx)
            np = np + 3
         end do
         if (first) then
            ! E_2 = o (2 J - c_x K) o of the ground state, o = C_occ C_occ^T
            ip = np + 1
            pairs(ip) = trc_grad_pair_t(1_int32, 1_int32, int(ncol, int32), 2.0_dp, st%cx)
         end if
         call eri%gradient_pairs(b, pool, pairs, ncol, gp, opts%err)
         if (opts%err%has_error()) return
         opts%nsweep = opts%nsweep + 1
         opts%grad(:, :, c0:c0 + nc - 1) = opts%grad(:, :, c0:c0 + nc - 1) + gp(:, :, 1:nc)
         if (first) gground = gp(:, :, ncol)
         deallocate (pool, pairs, gp)
         c0 = c0 + nc
      end do
   end subroutine tddft_grad_two_electron

   subroutine tddft_grad_xc(st, b, la, ref, xg, opts, gground)
      !! The XC terms: the ground state's (added to `gground`) and, per chunk
      !! of roots, trc_xc_grad_response_many with P (the relaxed difference
      !! density) and 2 sym(T) = sqrt(2) S, for which 1/2 E''[T, T] is the
      !! response energy 2 E''[sym(X), sym(X)] the excitation energy carries.
      type(tddft_grad_state_t), intent(in) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in) :: xg
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      real(dp), intent(inout) :: gground(3, b%natm)

      type(trc_xc_grad_response_options_t) :: xo
      type(xc_grad_options_t) :: xopts
      type(step_table_t) :: xst
      real(dp), allocatable :: density(:, :, :), g0(:, :), columns(:, :, :, :), gx(:, :, :)
      integer :: c0, nc, j, k

      xopts = xc_grad_options(grid_response=opts%grid_response, partition=opts%partition)
      call xst%init("xc grad", xc_grad_step_names())
      allocate (density(st%n, st%n, 1), g0(3, b%natm))
      density(:, :, 1) = ref%dmat
      g0 = 0.0_dp
      !$acc enter data copyin(density, g0)
      call xc_grad_add_on_grid(b, la, ref%func, density, g0, xopts, xst, xg)
      !$acc update self(g0)
      !$acc exit data delete(density, g0)
      gground = gground + g0

      xo%gradient = xopts
      c0 = 1
      do while (c0 <= st%nroots)
         call tddft_grad_chunk(st, opts, c0, nc)
         allocate (columns(st%n, st%n, 2, nc), gx(3, b%natm, nc))
         do j = 1, nc
            k = c0 + j - 1
            columns(:, :, 1, j) = st%prel(:, :, k)
            columns(:, :, 2, j) = st%tmat(:, :, k) + transpose(st%tmat(:, :, k))
         end do
         gx = 0.0_dp
         !$acc enter data copyin(gx)
         call trc_xc_grad_response_many(b, la, ref%func, xg, ref%dmat, columns, gx, xo)
         !$acc update self(gx)
         !$acc exit data delete(gx)
         if (xo%err%has_error()) then
            opts%err = xo%err
            return
         end if
         opts%grad(:, :, c0:c0 + nc - 1) = opts%grad(:, :, c0:c0 + nc - 1) + gx
         deallocate (columns, gx)
         c0 = c0 + nc
      end do
   end subroutine tddft_grad_xc

   subroutine tddft_grad_one_electron(st, b, pl, ref, opts)
      !! h^x . (C C^T 2 + P) - S^x . (W_ground + W): the ground state's
      !! and every root's one-electron and overlap terms, one call per root.
      type(tddft_grad_state_t), intent(in) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      real(dp), allocatable :: w0(:, :), wk(:, :), dk(:, :), g1(:, :), half(:, :)
      integer :: k, i

      allocate (w0(st%n, st%n), wk(st%n, st%n), dk(st%n, st%n), g1(3, b%natm), half(st%n, st%nocc))
      ! W_ground = 2 sum_i e_i c_i c_i^T
      do i = 1, st%nocc
         half(:, i) = 2.0_dp*st%eps(i)*st%cocc(:, i)
      end do
      call tddft_grad_gemm(half, st%cocc, w0, transb="T")
      do k = 1, st%nroots
         wk = w0 + st%wexc(:, :, k)
         dk = ref%dmat + st%prel(:, :, k)
         g1 = 0.0_dp
         !$acc enter data copyin(wk, dk, g1)
         call trc_1e_gradient(b, pl, wk, dk, g1)
         !$acc update self(g1)
         !$acc exit data delete(wk, dk, g1)
         opts%grad(:, :, k) = opts%grad(:, :, k) + g1
      end do
   end subroutine tddft_grad_one_electron

   subroutine trc_tddft_gradient_rks(b, pl, eri, la, ref, td, opts, xg)
      !!
      !! The nuclear gradient of the total energy (ground state plus
      !! excitation energy) of every root of `td`, a converged
      !! `trc_tddft_solve_rks` result for the closed-shell reference `ref`:
      !! `opts%grad(3, natm, nroots)`, hartree per bohr. `b` and `pl` are on the
      !! device, `eri` is built for `b` (the engine the solve used does),
      !! `xg` is the grid the SCF converged on for Kohn-Sham and absent for
      !! Hartree-Fock. A refusal (trc_tddft_grad_refuse_reason, an unconverged
      !! `td`, a failed Z-vector solve) is in `opts%err`, with `opts%grad`
      !! unallocated. This module's header has the method.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_result_t), intent(in) :: td
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      type(trc_xc_grid_t), intent(in), optional :: xg

      type(tddft_grad_state_t) :: st
      real(dp), allocatable :: gground(:, :), gnuc(:, :)
      character(len=256) :: reason
      integer :: k

      opts%err = error_t()
      opts%nsweep = 0
      if (allocated(opts%grad)) deallocate (opts%grad)
      if (allocated(opts%cpks_iterations)) deallocate (opts%cpks_iterations)
      reason = trc_tddft_grad_refuse_reason(ref, td, b%nao)
      if (len_trim(reason) == 0 .and. b%spherical) reason = "trc_tddft_grad: a spherical basis is not implemented"
      if (len_trim(reason) == 0 .and. ref%dft .and. .not. present(xg)) then
         reason = "trc_tddft_grad: a functional needs the SCF's xc grid"
      end if
      if (len_trim(reason) > 0) then
         call opts%err%set(ERROR_VALIDATION, trim(reason))
         return
      end if

      ! The gradient of a root inside a degenerate multiplet is an arbitrary
      ! rotation of its partners' (any rotation of X is as good a solution).
      if (allocated(td%omega)) then
         do k = 1, size(td%omega) - 1
            if (abs(td%omega(k + 1) - td%omega(k)) < DEGENERACY_TOL) then
               call global_logger%warning("trc_tddft_grad: roots "//to_string(k)//" and "//to_string(k + 1)// &
                                          " are degenerate (|d omega| < 1e-6): their gradients are not unique")
            end if
         end do
      end if

      call tddft_grad_init(st, ref, td)
      call tddft_grad_rhs(st, b, eri, la, ref, xg, opts%err)
      if (opts%err%has_error()) return
      call tddft_grad_zvector(st, b, eri, la, ref, xg, opts)
      if (opts%err%has_error()) return
      call tddft_grad_relaxed(st, b, eri, la, ref, xg, opts%err)
      if (opts%err%has_error()) return

      allocate (opts%grad(3, b%natm, st%nroots), gground(3, b%natm), gnuc(3, b%natm))
      opts%grad = 0.0_dp
      call tddft_grad_two_electron(st, b, eri, opts, gground)
      if (opts%err%has_error()) then
         deallocate (opts%grad)
         return
      end if
      if (st%dft) then
         call tddft_grad_xc(st, b, la, ref, xg, opts, gground)
         if (opts%err%has_error()) then
            deallocate (opts%grad)
            return
         end if
      end if
      call tddft_grad_one_electron(st, b, pl, ref, opts)
      gnuc = 0.0_dp
      call scf_grad_nuclear_repulsion(b, gnuc)
      do k = 1, st%nroots
         opts%grad(:, :, k) = opts%grad(:, :, k) + gground + gnuc
      end do
      if (opts%verbose) then
         call global_logger%info("tddft grad: "//to_string(st%nroots)//" roots, "//to_string(opts%nsweep)// &
                                 " integral-derivative sweep(s)")
      end if
   end subroutine trc_tddft_gradient_rks

end module trc_tddft_grad
