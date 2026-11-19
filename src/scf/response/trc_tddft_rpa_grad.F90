module trc_tddft_rpa_grad
   !!
   !! Analytic nuclear gradients of full-RPA (TDHF / Casida TDDFT) singlet
   !! excited states of a closed-shell RHF/RKS reference, for several roots at
   !! once: the Furche-Ahlrichs Z-vector gradient with Y /= 0, term by term
   !! that of pyscf's grad/tdrhf.py and grad/tdrks.py `grad_elec`.
   !!
   !! It is trc_tddft_grad (read its header first) with the (X, Y) pair in
   !! place of the TDA amplitude. What changes, with terco's X.X - Y.Y = 1
   !! (pyscf's x, y are these over sqrt(2)):
   !!
   !!     xpy = (X + Y)^T/sqrt(2),  xmy = (X - Y)^T/sqrt(2)   pyscf `xpy`, `xmy`
   !!     T+ = C_occ (X + Y) C_vir^T,  T- = C_occ (X - Y) C_vir^T
   !!     S = (T+ + T+^T)/sqrt(2)     symmetric, from X + Y   (pyscf dm[1])
   !!     A = (T-^T - T-)/sqrt(2)     antisymmetric, from X - Y (pyscf dm[2])
   !!     P_d = C_occ (-(X X^T + Y Y^T)) C_occ^T + C_vir (X^T X + Y^T Y) C_vir^T
   !!                                 the unrelaxed difference density
   !!
   !! and then
   !!  - the Z-vector's right-hand side takes the occupied-occupied and
   !!    virtual-virtual blocks of G[S] with xpy and those of G[A] with xmy;
   !!  - so do the energy-weighted density's blocks (pyscf `im0`);
   !!  - the XC response (f_xc, k_xc, and the grid-point derivatives) sees
   !!    only X + Y, through S;
   !!  - the two-electron derivative is the same three bilinear forms, with
   !!    S and A now built from different matrices.
   !! With Y = 0 every line reduces to trc_tddft_grad's, and a TDA result
   !! (whose `yamp` is zero) gives its numbers.
   !!
   !! The state, the Fock, f_xc, Z-vector, XC-response and one-electron steps
   !! are trc_tddft_grad's, shared as they stand; the right-hand side, the
   !! relaxed density and the pairs of the derivative sweep are here.
   !!
   !! SCOPE. As trc_tddft_grad: restricted, singlet, LDA/GGA and global
   !! hybrids, Cartesian basis, no density fitting.
   !!
   use pic_types, only: dp, int32
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_eri, only: trc_eri_t, trc_grad_pair_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc, only: trc_xc_fxc_cache_t, trc_xc_kxc_rks_apply
   use trc_scf_grad, only: scf_grad_nuclear_repulsion
   use trc_tddft, only: trc_tddft_result_t, trc_tddft_reference_t
   use trc_tddft_grad, only: trc_tddft_grad_options_t, trc_tddft_grad_refuse_reason, tddft_grad_state_t, tddft_grad_gemm, &
                             tddft_grad_basis_change, tddft_grad_mo_block, tddft_grad_symmetrize, tddft_grad_fock, &
                             tddft_grad_fock_triples, &
                             tddft_grad_fxc, tddft_grad_zvector, tddft_grad_chunk, tddft_grad_xc, &
                             tddft_grad_one_electron, SQRT_HALF
   implicit none
   private

   public :: trc_tddft_rpa_gradient_rks

   real(dp), parameter :: DEGENERACY_TOL = 1.0e-6_dp
   !! Two adjacent requested roots closer than this (hartree) are warned about as degenerate.

   type, extends(tddft_grad_state_t) :: rpa_grad_state_t
      !! trc_tddft_grad's state, with what the RPA adds. `xamp` is X, `tmat` is T+.
      real(dp), allocatable :: yamp(:, :, :)     !! (nocc, nvir, nroots) Y
      real(dp), allocatable :: tminus(:, :, :)   !! T- = C_occ (X - Y) C_vir^T
      real(dp), allocatable :: xpy(:, :, :)      !! (nvir, nocc, nroots) (X + Y)^T/sqrt(2), pyscf's xpy
      real(dp), allocatable :: xmy(:, :, :)      !! (nvir, nocc, nroots) (X - Y)^T/sqrt(2), pyscf's xmy
   end type rpa_grad_state_t

contains

   subroutine rpa_grad_init(st, ref, td)
      !! The orbitals, T+, T- and the unrelaxed difference density of every root.
      type(rpa_grad_state_t), intent(out) :: st
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_result_t), intent(in) :: td
      real(dp), allocatable :: doo(:, :), dvv(:, :), half(:, :), xk(:, :), yk(:, :), plus(:, :), minus(:, :)
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
      st%yamp = td%yamp
      allocate (st%tmat(st%n, st%n, st%nroots), st%tminus(st%n, st%n, st%nroots), st%pdiff(st%n, st%n, st%nroots))
      allocate (st%xpy(st%nvir, st%nocc, st%nroots), st%xmy(st%nvir, st%nocc, st%nroots))
      allocate (doo(st%nocc, st%nocc), dvv(st%nvir, st%nvir), half(st%n, st%n))
      do k = 1, st%nroots
         xk = st%xamp(:, :, k)
         yk = st%yamp(:, :, k)
         plus = xk + yk
         minus = xk - yk
         call tddft_grad_basis_change(st%cocc, plus, st%cvir, st%tmat(:, :, k))
         call tddft_grad_basis_change(st%cocc, minus, st%cvir, st%tminus(:, :, k))
         st%xpy(:, :, k) = SQRT_HALF*transpose(plus)
         st%xmy(:, :, k) = SQRT_HALF*transpose(minus)
         call rpa_grad_diff_blocks(xk, yk, doo, dvv)
         call tddft_grad_basis_change(st%cocc, doo, st%cocc, st%pdiff(:, :, k))
         call tddft_grad_basis_change(st%cvir, dvv, st%cvir, half)
         st%pdiff(:, :, k) = st%pdiff(:, :, k) + half
         call tddft_grad_symmetrize(st%pdiff(:, :, k))
      end do
   end subroutine rpa_grad_init

   subroutine rpa_grad_diff_blocks(xk, yk, doo, dvv)
      !! The occupied and virtual blocks of the unrelaxed difference density in
      !! the MO basis (pyscf doo, dvv): -(X X^T + Y Y^T) and X^T X + Y^T Y.
      real(dp), intent(in) :: xk(:, :), yk(:, :)
      real(dp), intent(out) :: doo(:, :), dvv(:, :)
      call tddft_grad_gemm(xk, xk, doo, transb="T", alpha=-1.0_dp)
      call tddft_grad_gemm(yk, yk, doo, transb="T", alpha=-1.0_dp, beta=1.0_dp)
      call tddft_grad_gemm(xk, xk, dvv, transa="T")
      call tddft_grad_gemm(yk, yk, dvv, transa="T", beta=1.0_dp)
   end subroutine rpa_grad_diff_blocks

   subroutine rpa_grad_densities(st, dm3)
      !! The three AO densities of every root the first Fock pass contracts:
      !! P_d, S (from T+) and A (from T-), as columns 3k-2, 3k-1, 3k.
      type(rpa_grad_state_t), intent(in) :: st
      real(dp), intent(out) :: dm3(:, :, :)
      integer :: k
      do k = 1, st%nroots
         dm3(:, :, 3*k - 2) = st%pdiff(:, :, k)
         dm3(:, :, 3*k - 1) = SQRT_HALF*(st%tmat(:, :, k) + transpose(st%tmat(:, :, k)))
         dm3(:, :, 3*k) = SQRT_HALF*(transpose(st%tminus(:, :, k)) - st%tminus(:, :, k))
      end do
   end subroutine rpa_grad_densities

   subroutine rpa_grad_xc_terms(st, b, la, ref, xg, dm3, xcmat, err)
      !! f_xc[P_d], f_xc[S] and k1ao = E_xc'''[S, S, .] of every root (pyscf
      !! f1oo, f1vo, k1ao), from the densities of `rpa_grad_densities`, packed
      !! in `xcmat` (n, n, 3 nroots): columns 2k-1, 2k are f_xc[P_d], f_xc[S] of
      !! root k and 2 nroots + k its k1ao. The XC kernels see X + Y only.
      type(rpa_grad_state_t), intent(in) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: dm3(:, :, :)
      real(dp), intent(out) :: xcmat(:, :, :)
      type(error_t), intent(out) :: err

      type(trc_xc_fxc_cache_t) :: cache
      real(dp), allocatable :: dmx(:, :, :), dmk(:, :, :), fx(:, :, :), vk(:, :, :)
      integer :: k, nr
      nr = st%nroots
      allocate (dmx(st%n, st%n, 2*nr), dmk(st%n, st%n, 1 + 2*nr), fx(st%n, st%n, 2*nr), vk(st%n, st%n, nr))
      call cache%setup(b, xg, ref%func, ref%dmat, la)
      do k = 1, nr
         dmx(:, :, 2*k - 1) = dm3(:, :, 3*k - 2)
         dmx(:, :, 2*k) = dm3(:, :, 3*k - 1)
      end do
      call tddft_grad_fxc(st%tddft_grad_state_t, b, xg, cache, la, dmx, fx)
      call cache%release()
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
      xcmat(:, :, 1:2*nr) = fx
      xcmat(:, :, 2*nr + 1:3*nr) = vk
   end subroutine rpa_grad_xc_terms

   subroutine rpa_grad_rhs(st, b, eri, la, ref, xg, err)
      !! Step 1: the right-hand side of the Z-vector equation of every root
      !! (pyscf `wvo`), with the MO blocks it and the energy-weighted density
      !! need. `err` is set on a refused kernel.
      type(rpa_grad_state_t), intent(inout) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in), optional :: xg
      type(error_t), intent(out) :: err

      real(dp), allocatable :: dm3(:, :, :), g3(:, :, :), xcmat(:, :, :), veff_s(:, :), veff_a(:, :)
      integer :: n, nocc, nvir, nr, k
      n = st%n
      nocc = st%nocc
      nvir = st%nvir
      nr = st%nroots
      allocate (dm3(n, n, 3*nr), g3(n, n, 3*nr), veff_s(n, n), veff_a(n, n))
      allocate (st%veff0(n, n, nr), st%voo(nocc, nocc, nr), st%vvo(nvir, nocc, nr), st%vvv(nvir, nvir, nr))
      allocate (st%aoo(nocc, nocc, nr), st%avo(nvir, nocc, nr), st%avv(nvir, nvir, nr), st%wvo(nvir, nocc, nr))
      call rpa_grad_densities(st, dm3)
      ! 2 J[D] - c_x K[D] of P_d, S and A (pyscf vj*2 - vk*hyb). A is antisymmetric, so the triples
      ! helper builds it through its own exchange-only batch: the symmetric folded build
      ! (tddft_grad_fock) is not valid for it.
      call tddft_grad_fock_triples(b, eri, nr, dm3, g3, 2.0_dp, 2.0_dp*st%cx)
      allocate (xcmat(n, n, 3*nr))
      xcmat = 0.0_dp
      if (st%dft) then
         call rpa_grad_xc_terms(st, b, la, ref, xg, dm3, xcmat, err)
         if (err%has_error()) return
      end if
      do k = 1, nr
         st%veff0(:, :, k) = g3(:, :, 3*k - 2) + 2.0_dp*xcmat(:, :, 2*k - 1) + 2.0_dp*xcmat(:, :, 2*nr + k)
         veff_s = g3(:, :, 3*k - 1) + 2.0_dp*xcmat(:, :, 2*k)
         veff_a = g3(:, :, 3*k)
         call tddft_grad_mo_block(st%cocc, veff_s, st%cocc, st%voo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_s, st%cocc, st%vvo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_s, st%cvir, st%vvv(:, :, k))
         call tddft_grad_mo_block(st%cocc, veff_a, st%cocc, st%aoo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_a, st%cocc, st%avo(:, :, k))
         call tddft_grad_mo_block(st%cvir, veff_a, st%cvir, st%avv(:, :, k))
         ! wvo = 2 C_v^T veff0doo C_o - 2 xpy Voo^T + 2 Vvv^T xpy
         !                            - 2 xmy Aoo^T + 2 Avv^T xmy
         call tddft_grad_mo_block(st%cvir, st%veff0(:, :, k), st%cocc, st%wvo(:, :, k))
         st%wvo(:, :, k) = 2.0_dp*st%wvo(:, :, k)
         call tddft_grad_gemm(st%xpy(:, :, k), st%voo(:, :, k), st%wvo(:, :, k), transb="T", alpha=-2.0_dp, beta=1.0_dp)
         call tddft_grad_gemm(st%vvv(:, :, k), st%xpy(:, :, k), st%wvo(:, :, k), transa="T", alpha=2.0_dp, beta=1.0_dp)
         call tddft_grad_gemm(st%xmy(:, :, k), st%aoo(:, :, k), st%wvo(:, :, k), transb="T", alpha=-2.0_dp, beta=1.0_dp)
         call tddft_grad_gemm(st%avv(:, :, k), st%xmy(:, :, k), st%wvo(:, :, k), transa="T", alpha=2.0_dp, beta=1.0_dp)
      end do
   end subroutine rpa_grad_rhs

   subroutine rpa_grad_weighted(st, k, veffz, wexc)
      !! The excitation's energy-weighted density of root k (pyscf `im0` plus
      !! `zeta dm1`, without the ground state's part): `veffz` is
      !! veff0doo + veff of the relaxed Z response, AO.
      type(rpa_grad_state_t), intent(in) :: st
      integer, intent(in) :: k
      real(dp), intent(in) :: veffz(:, :)
      real(dp), intent(out) :: wexc(:, :)
      real(dp), allocatable :: woo(:, :), wvv(:, :), wvo(:, :), doo(:, :), dvv(:, :), ao(:, :)
      integer :: nocc, nvir, i, j
      nocc = st%nocc
      nvir = st%nvir
      allocate (woo(nocc, nocc), wvv(nvir, nvir), wvo(nvir, nocc), doo(nocc, nocc), dvv(nvir, nvir), ao(st%n, st%n))
      ! oo: C_o^T (veff0doo + veff) C_o + V_vo^T xpy + A_vo^T xmy
      call tddft_grad_mo_block(st%cocc, veffz, st%cocc, woo)
      call tddft_grad_gemm(st%vvo(:, :, k), st%xpy(:, :, k), woo, transa="T", beta=1.0_dp)
      call tddft_grad_gemm(st%avo(:, :, k), st%xmy(:, :, k), woo, transa="T", beta=1.0_dp)
      ! vv: xpy V_vo^T + xmy A_vo^T
      call tddft_grad_gemm(st%xpy(:, :, k), st%vvo(:, :, k), wvv, transb="T")
      call tddft_grad_gemm(st%xmy(:, :, k), st%avo(:, :, k), wvv, transb="T", beta=1.0_dp)
      ! vo: 2 xpy V_oo^T + 2 xmy A_oo^T
      call tddft_grad_gemm(st%xpy(:, :, k), st%voo(:, :, k), wvo, transb="T", alpha=2.0_dp)
      call tddft_grad_gemm(st%xmy(:, :, k), st%aoo(:, :, k), wvo, transb="T", alpha=2.0_dp, beta=1.0_dp)
      ! + zeta * dm1: dm1 = (doo, dvv, Z), zeta_pq = (e_p + e_q)/2, zeta_ai = e_i
      call rpa_grad_diff_blocks(st%xamp(:, :, k), st%yamp(:, :, k), doo, dvv)
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
      call tddft_grad_basis_change(st%cocc, woo, st%cocc, wexc)
      call tddft_grad_basis_change(st%cvir, wvv, st%cvir, ao)
      wexc = wexc + ao
      call tddft_grad_basis_change(st%cvir, wvo, st%cocc, ao)
      wexc = wexc + ao
      call tddft_grad_symmetrize(wexc)
   end subroutine rpa_grad_weighted

   subroutine rpa_grad_relaxed(st, b, eri, la, ref, xg, err)
      !! Step 3: P = P_d + sym(Z) and the excitation's energy-weighted density.
      type(rpa_grad_state_t), intent(inout) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_xc_grid_t), intent(in), optional :: xg
      type(error_t), intent(out) :: err

      type(trc_xc_fxc_cache_t) :: cache
      real(dp), allocatable :: dmz(:, :, :), gz(:, :, :), fz(:, :, :), ao(:, :)
      integer :: n, nr, k
      n = st%n
      nr = st%nroots
      allocate (dmz(n, n, nr), gz(n, n, nr), fz(n, n, nr), ao(n, n))
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
         call tddft_grad_fxc(st%tddft_grad_state_t, b, xg, cache, la, dmz, fz)
         call cache%release()
      end if
      do k = 1, nr
         call rpa_grad_weighted(st, k, st%veff0(:, :, k) + gz(:, :, k) + fz(:, :, k), st%wexc(:, :, k))
         st%prel(:, :, k) = st%pdiff(:, :, k) + 0.5_dp*dmz(:, :, k)
      end do
      err = error_t()
   end subroutine rpa_grad_relaxed

   subroutine rpa_grad_two_electron(st, b, eri, opts, gground)
      !! The two-electron derivative terms of every root, and of the ground
      !! state (`gground`), from one `gradient_pairs` sweep per chunk of roots:
      !! trc_tddft_grad's pairs, with A built from T-. See that module's header
      !! for the pairs and the sign of the antisymmetric one.
      type(rpa_grad_state_t), intent(in) :: st
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      real(dp), intent(out) :: gground(3, b%natm)

      type(trc_grad_pair_t), allocatable :: pairs(:)
      real(dp), allocatable :: pool(:, :, :), gp(:, :, :), occ0(:, :)
      integer :: c0, nc, ncol, j, k, np
      integer(int32) :: pcol, scol, acol, col
      logical :: first

      allocate (occ0(st%n, st%n))
      call tddft_grad_gemm(st%cocc, st%cocc, occ0, transb="T")
      gground = 0.0_dp
      c0 = 1
      do while (c0 <= st%nroots)
         call tddft_grad_chunk(st%tddft_grad_state_t, opts, c0, nc)
         first = c0 == 1
         ncol = nc + merge(1, 0, first)
         allocate (pool(st%n, st%n, 1 + 3*nc), pairs(3*nc + merge(1, 0, first)), gp(3, b%natm, ncol))
         pool(:, :, 1) = occ0
         np = 0
         do j = 1, nc
            k = c0 + j - 1
            col = int(j, int32)
            pool(:, :, 3*j - 1) = st%prel(:, :, k)
            pool(:, :, 3*j) = SQRT_HALF*(st%tmat(:, :, k) + transpose(st%tmat(:, :, k)))
            pool(:, :, 3*j + 1) = SQRT_HALF*(transpose(st%tminus(:, :, k)) - st%tminus(:, :, k))
            pcol = 3_int32*col - 1_int32
            scol = pcol + 1_int32
            acol = pcol + 2_int32
            ! P against the ground density (cross: both orderings), S and A on their own.
            pairs(np + 1) = trc_grad_pair_t(pcol, 1_int32, col, 1.0_dp, 0.5_dp*st%cx)
            pairs(np + 2) = trc_grad_pair_t(scol, scol, col, 1.0_dp, 0.5_dp*st%cx)
            ! The antisymmetric pair's exchange is differentiated with the
            ! transposed right density, K[A^T] = -K[A]: the literal
            ! -c_x/2 A.K[A] is the weight +c_x/2 turned negative.
            pairs(np + 3) = trc_grad_pair_t(acol, acol, col, 0.0_dp, -0.5_dp*st%cx)
            np = np + 3
         end do
         if (first) then
            ! E_2 = o (2 J - c_x K) o of the ground state, o = C_occ C_occ^T
            pairs(np + 1) = trc_grad_pair_t(1_int32, 1_int32, int(ncol, int32), 2.0_dp, st%cx)
         end if
         call eri%gradient_pairs(b, pool, pairs, ncol, gp, opts%err)
         if (opts%err%has_error()) return
         opts%nsweep = opts%nsweep + 1
         opts%grad(:, :, c0:c0 + nc - 1) = opts%grad(:, :, c0:c0 + nc - 1) + gp(:, :, 1:nc)
         if (first) gground = gp(:, :, ncol)
         deallocate (pool, pairs, gp)
         c0 = c0 + nc
      end do
   end subroutine rpa_grad_two_electron

   function rpa_grad_refuse_reason(ref, td, b) result(msg)
      !! Why this gradient cannot run, or an empty string: trc_tddft_grad's
      !! reasons, plus a result without Y amplitudes of the shape of X.
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_result_t), intent(in) :: td
      type(trc_basis_t), intent(in) :: b
      character(len=256) :: msg
      msg = ""
      if (.not. td%converged .or. .not. allocated(td%xamp)) then
         msg = "trc_tddft_rpa_grad: needs the amplitudes of a converged TDA or RPA solve"
         return
      end if
      msg = trc_tddft_grad_refuse_reason(ref, td, b%nao)
      if (len_trim(msg) > 0) return
      if (.not. allocated(td%yamp)) then
         msg = "trc_tddft_rpa_grad: the result has no Y amplitudes (yamp); an RPA solve fills them"
      else if (any(shape(td%yamp) /= shape(td%xamp))) then
         msg = "trc_tddft_rpa_grad: yamp and xamp differ in shape"
      else if (b%spherical) then
         msg = "trc_tddft_rpa_grad: a spherical basis is not implemented"
      end if
   end function rpa_grad_refuse_reason

   subroutine trc_tddft_rpa_gradient_rks(b, pl, eri, la, ref, td, opts, xg)
      !!
      !! The nuclear gradient of the total energy (ground state plus
      !! excitation energy) of every root of `td`, a converged
      !! `trc_tddft_solve_rks` result (`tda = .false.`, so with Y; a TDA
      !! result has Y = 0 and gives the numbers of `trc_tddft_gradient_rks`)
      !! for the closed-shell reference `ref`: `opts%grad(3, natm, nroots)`,
      !! hartree per bohr. The arguments and `opts` are those of
      !! `trc_tddft_gradient_rks`. A refusal (an unconverged `td`, a result
      !! without `yamp`, what trc_tddft_grad_refuse_reason says, a failed
      !! Z-vector solve) is in `opts%err`, with `opts%grad` unallocated.
      !! This module's header has the method.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_result_t), intent(in) :: td
      type(trc_tddft_grad_options_t), intent(inout) :: opts
      type(trc_xc_grid_t), intent(in), optional :: xg

      type(rpa_grad_state_t) :: st
      real(dp), allocatable :: gground(:, :), gnuc(:, :)
      character(len=256) :: reason
      integer :: k

      opts%err = error_t()
      opts%nsweep = 0
      if (allocated(opts%grad)) deallocate (opts%grad)
      if (allocated(opts%cpks_iterations)) deallocate (opts%cpks_iterations)
      reason = rpa_grad_refuse_reason(ref, td, b)
      if (len_trim(reason) == 0 .and. ref%dft .and. .not. present(xg)) then
         reason = "trc_tddft_rpa_grad: a functional needs the SCF's xc grid"
      end if
      if (len_trim(reason) > 0) then
         call opts%err%set(ERROR_VALIDATION, trim(reason))
         return
      end if

      ! The gradient of a root inside a degenerate multiplet is an arbitrary
      ! rotation of its partners' (any rotation of (X, Y) is as good a solution).
      ! Only adjacent roots of the requested set are compared; a partner of
      ! the multiplet beyond the last requested root is not seen.
      if (allocated(td%omega)) then
         do k = 1, size(td%omega) - 1
            if (abs(td%omega(k + 1) - td%omega(k)) < DEGENERACY_TOL) then
               call global_logger%warning("trc_tddft_rpa_grad: roots "//to_string(k)//" and "//to_string(k + 1)// &
                                          " are degenerate (|d omega| < 1e-6): their gradients are not unique")
            end if
         end do
      end if

      call rpa_grad_init(st, ref, td)
      call rpa_grad_rhs(st, b, eri, la, ref, xg, opts%err)
      if (opts%err%has_error()) return
      call tddft_grad_zvector(st%tddft_grad_state_t, b, eri, la, ref, xg, opts)
      if (opts%err%has_error()) return
      call rpa_grad_relaxed(st, b, eri, la, ref, xg, opts%err)
      if (opts%err%has_error()) return

      allocate (opts%grad(3, b%natm, st%nroots), gground(3, b%natm), gnuc(3, b%natm))
      opts%grad = 0.0_dp
      call rpa_grad_two_electron(st, b, eri, opts, gground)
      if (opts%err%has_error()) then
         deallocate (opts%grad)
         return
      end if
      if (st%dft) then
         call tddft_grad_xc(st%tddft_grad_state_t, b, la, ref, xg, opts, gground)
         if (opts%err%has_error()) then
            deallocate (opts%grad)
            return
         end if
      end if
      call tddft_grad_one_electron(st%tddft_grad_state_t, b, pl, ref, opts)
      gnuc = 0.0_dp
      call scf_grad_nuclear_repulsion(b, gnuc)
      do k = 1, st%nroots
         opts%grad(:, :, k) = opts%grad(:, :, k) + gground + gnuc
      end do
      if (opts%verbose) then
         call global_logger%info("tddft rpa grad: "//to_string(st%nroots)//" roots, "//to_string(opts%nsweep)// &
                                 " integral-derivative sweep(s)")
      end if
   end subroutine trc_tddft_rpa_gradient_rks

end module trc_tddft_rpa_grad
