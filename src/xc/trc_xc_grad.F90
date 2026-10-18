!
! The exchange-correlation contribution to the nuclear gradient, LDA and
! GGA, restricted and unrestricted. Added into an existing `grad(3,natm)`,
! the same convention trc_1e_gradient and nuclear_repulsion_gradient use
! (trc_scf_grad.F90).
!
! THE FORMULA. For every atom A, the only thing that moves is the basis
! functions centred on it (the grid and its Becke weights are held fixed
! -- see GRID RESPONSE below), so differentiating
!
!     E_xc = sum_g w_g f(rho_g, sigma_g)
!
! with rho_g = sum_uv D_uv chi_u(g) chi_v(g) and sigma_g = |grad rho_g|^2
! gives, restricted,
!
!     dE_xc/dR_{A,k} = -2 sum_{u in A} sum_v D_uv sum_g w_g
!                      [ vrho_g dchi_u/dx_k(g) chi_v(g)
!                        + 2 vsigma_g ( chi_v(g) (grad_rho_g . H_u(g))_k
!                                       + dchi_u/dx_k(g) (grad_rho_g . grad_chi_v(g)) ) ]
!
! where H_u is chi_u's Hessian. LDA is the same formula with vsigma = 0 (a
! pure-LDA kernel composition already returns vsigma = 0 from
! xc_eval_point, so nothing here is skipped for correctness -- only for
! cost: LDA never collocates a Hessian). Unrestricted replaces the single
! vrho/vsigma pair with the five spin-resolved derivatives
! (va, vb, v_aa, v_ab, v_bb) and a density matrix per channel; see
! xc_grad_point_gga_body for the channel-by-channel expansion. Checked,
! term for term, against metalquicha's independently-validated
! accumulate_channel/accumulate_gga_channel (mqc_czt_gradient.f90).
!
! GRID RESPONSE -- NOT IMPLEMENTED, BY DEFAULT. Two more terms exist and
! are left out here: the grid POINTS are atom-centred and travel with
! their owning atom, and the Becke PARTITION WEIGHT of every point depends
! on every nuclear position. Both are variational-grid artefacts (a grid
! matched to no finite order quietly leaks net force) rather than part of
! the Pulay term above, and both shrink as the grid is refined. Measured
! on terco's own water/6-31G* (check_xc_grad_fd.F90, comparing this
! module's analytic gradient to a finite difference of terco's own
! converged energy, at TWO grid levels) rather than assumed: see that
! test's output for the level-3-vs-level-5 numbers this header does not
! hardcode. Omitting the term is the DEFAULT and is not silent: a caller
! that needs the physical gradient to the last digit (a tight geometry
! optimisation, a numerical Hessian built from this gradient) should use a
! finer grid instead, the same remedy metalquicha's own fixed_grid flag
! documents on its two-density counterpart.
!
! BATCHING -- deliberately SIMPLER than the energy path's. trc_xc.F90's
! xc_drive pads every batch's tile to a common shape and groups batches of
! equal padded shape into one strided-batched cuBLAS call, because it runs
! every SCF iteration and the grouping is what makes that affordable. This
! module runs once, after convergence, so it processes one spatial batch
! at a time -- collocate, gather, a handful of exact-sized GEMMs, scatter
! -- with no padding and no cross-batch grouping. Still cuBLAS (`trc_linalg`),
! still batched (one pass per `trc_xc_grid_t` batch, same screening), just
! not padding-grouped; see this branch's own report for the measured cost
! against the two-electron gradient. Grouping batches the way the energy
! path does is the natural next step if this ever dominates a gradient's
! wall time.
!
! CARTESIAN ONLY. XC is Cartesian throughout terco (the SCF folds to
! spherical); `b` and `dens` here are exactly what trc_scf_grad.F90 hands
! in, with no spherical folding of its own.
!
! GPU: CONFIRMED WRONG, NOT MERELY UNTESTED. check_full_grad_fd (the
! existing one/two-electron gradient) and scf_rks (the existing XC
! ENERGY path, same atomic-scatter idiom as this module's own
! xc_grad_scatter) both pass on this exact GPU build, which rules out the
! basis/pairlist/eri/trc_linalg/XC-energy machinery this module reuses.
! What is NEW here -- and the only thing left to suspect -- is this
! module's own per-batch device lifetime: unlike xc_drive, which
! allocates its device tiles ONCE outside the batch loop and reuses them
! (see BATCHING above), xc_grad_one_batch allocates (`!$acc enter data
! create`) and frees (`!$acc exit data delete`) a fresh set of device
! buffers on EVERY batch. check_xc_grad_fixedgrid_fd and
! check_xc_grad_vs_pyscf both fail the same way on GPU 2, deterministically
! (repeated runs give the identical worst|diff|, not a flaky range), with
! the error growing with how much per-batch work there is (UKS GGA worse
! than RKS GGA worse than RKS LDA) -- consistent with a stale or
! overlapping device allocation across batches rather than a sign/index
! mistake, which CPU's bitwise-identical formula code rules out anyway.
! NOT fixed this round: restructuring the per-batch allocation to the
! energy path's hoist-and-reuse convention is a real change to this
! module's data flow, not a one-line fix, and is left for the GPU
! follow-up this header's BATCHING section already calls out as the next
! step. Use the CPU build for this module until then.
!
module trc_xc_grad
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_eval_point, xc_eval_point_polar, &
                                 XC_FAMILY_LDA, XC_FAMILY_GGA, XC_MAX_KERNELS
   use trc_collocation, only: shell_collocate, shell_collocate2, NCART_MAX, NHESS
   use trc_linalg, only: trc_linalg_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid, DEFAULT_GRID_LEVEL
   use trc_error, only: error_t
   implicit none
   private

   public :: xc_gradient_add, XC_GRAD_RHO_TOL

   real(dp), parameter :: XC_GRAD_RHO_TOL = 1.0e-12_dp
   !! Same role and value as trc_xc.F90's XC_RHO_TOL: a batch whose bound on
   !! rho is below this contributes nothing and is skipped before any
   !! device work for it.

contains

   !
   ! The one entry point. Builds its own DFT quadrature at `b`'s geometry
   ! (trc_scf_run keeps neither its grid nor its functional alive past
   ! convergence, the same reason trc_scf_gradient rebuilds the pair list
   ! and the ERI engine), refuses what it cannot do, and ADDS the result
   ! into `grad`.
   !
   subroutine xc_gradient_add(b, lp, func, nspin, dens, grad, grid_level, grid_max_pts, rho_tol)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      integer, intent(in) :: nspin
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      !! The converged density, per spin channel. Present on both host and
      !! device: the GEMM pipeline below reads the device copy, the
      !! per-batch screen reads the host one (trc_scf_gradient's own
      !! convention -- see its header -- keeps both valid at this point).
      real(dp), intent(inout) :: grad(3, b%natm)
      !! ADDED into; device-resident for the whole call, exactly as
      !! trc_1e_gradient and nuclear_repulsion_gradient require of their
      !! own `grad` argument.
      integer, intent(in), optional :: grid_level, grid_max_pts
      real(dp), intent(in), optional :: rho_tol

      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(error_t) :: err
      logical :: is_gga
      integer :: glevel, gmaxpts, ish, i
      real(dp) :: rtol
      integer, allocatable :: ao_atom(:)

      ! --- refuse what this module cannot do, with a reason -------------
      !
      ! Meta-GGA and double hybrids are not merely unimplemented here: on
      ! this branch `trc_xc_functional_t` has no meta-GGA family and no
      ! MP2-fraction field at all (see its own header), so neither can
      ! reach a converged SCF result in the first place -- there is no
      ! flag left to check for a double hybrid. The two checks below are
      ! for what CAN reach here: a range-separated hybrid (wb97x and
      ! friends) and VV10, both of which run in the energy/Fock path today.
      !
      if (func%family /= XC_FAMILY_LDA .and. func%family /= XC_FAMILY_GGA) then
         error stop "trc_xc_grad: only LDA and GGA exchange-correlation gradients "// &
            "are implemented (meta-GGA is not representable in trc_xc_functional_t "// &
            "on this branch, so this should be unreachable)"
      end if
      if (func%omega /= 0.0_dp .or. func%exx_lr /= 0.0_dp) then
         error stop "trc_xc_grad: range-separated (erf(omega r)/r) exchange has no "// &
            "gradient kernel; refusing rather than returning one missing that term"
      end if
      if (func%nlc_b /= 0.0_dp .or. func%nlc_c /= 0.0_dp) then
         error stop "trc_xc_grad: VV10 nonlocal correlation has no gradient kernel; "// &
            "refusing rather than returning one missing that term"
      end if

      is_gga = func%family == XC_FAMILY_GGA

      glevel = DEFAULT_GRID_LEVEL
      if (present(grid_level)) glevel = grid_level
      gmaxpts = 512
      if (present(grid_max_pts)) gmaxpts = grid_max_pts
      rtol = XC_GRAD_RHO_TOL
      if (present(rho_tol)) rtol = rho_tol

      call build_dft_grid(b%at_r, nint(b%at_z), grid, err, level=glevel)
      if (err%has_error()) error stop "trc_xc_grad: "//trim(err%get_message())
      call xg%build(grid%n_points, grid%coords, grid%weights, b, max_pts=gmaxpts)
      call xg%to_device()

      if (xg%npts > 0) then
         ! AO -> atom, built once from the shell -> atom map `basis_build`
         ! already produced (sh_at), the same map the two-electron and
         ! one-electron gradients use.
         allocate (ao_atom(b%nao))
         do ish = 1, b%nshell
            do i = 0, (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2 - 1
               ao_atom(b%sh_ao(ish) + i) = b%sh_at(ish)
            end do
         end do
         !$acc enter data copyin(ao_atom)
         call xc_grad_drive(b, lp, is_gga, nspin, func%nk, func%kid, func%coef, dens, xg, rtol, &
                            ao_atom, grad)
         !$acc exit data delete(ao_atom)
         deallocate (ao_atom)
      end if

      call xg%release()
      call grid%destroy()
   end subroutine xc_gradient_add

   !
   ! One pass over every batch of the grid: screen, then collocate + gather
   ! + GEMM + scatter for the ones that survive. See this module's header
   ! for why this is per-batch rather than the energy path's padded,
   ! shape-grouped chunks.
   !
   subroutine xc_grad_drive(b, lp, is_gga, nspin, nk, kid, coef, dens, xg, rtol, ao_atom, grad)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      logical, intent(in) :: is_gga
      integer, intent(in) :: nspin, nk
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: rtol
      integer, intent(in) :: ao_atom(b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)

      integer :: ib, npb, nloc

      do ib = 1, xg%nbatch
         npb = xg%b_off(ib + 1) - xg%b_off(ib)
         nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
         if (npb == 0 .or. nloc == 0) cycle
         if (.not. xc_grad_batch_active(ib, xg, b%nao, nspin, dens, rtol)) cycle
         call xc_grad_one_batch(b, lp, is_gga, nspin, nk, kid, coef, dens, xg, ib, npb, nloc, &
                                ao_atom, grad)
      end do
   end subroutine xc_grad_drive

   !
   ! Host-side density screen, bounding rho on a batch through the same
   ! per-function maxima (`xg%b_amax`) trc_xc.F90's xc_screen_body uses --
   ! plain Fortran, not a device kernel, because `dens`'s host copy already
   ! holds the right numbers (trc_scf_gradient's `!$acc update self(dens)`,
   ! before this is ever called) and nloc is small enough that this costs
   ! nothing next to the GEMMs a surviving batch pays for.
   !
   logical function xc_grad_batch_active(ib, xg, nao, nspin, dens, rtol) result(active)
      integer, intent(in) :: ib, nao, nspin
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: dens(nao, nao, nspin)
      real(dp), intent(in) :: rtol
      real(dp) :: bound, su
      integer :: s, u, v, au, av
      bound = 0.0_dp
      do s = 1, nspin
         do u = xg%b_aooff(ib), xg%b_aooff(ib + 1) - 1
            au = xg%b_ao(u)
            su = 0.0_dp
            do v = xg%b_aooff(ib), xg%b_aooff(ib + 1) - 1
               av = xg%b_ao(v)
               su = su + abs(dens(av, au, s))*xg%b_amax(v)
            end do
            bound = bound + su*xg%b_amax(u)
         end do
      end do
      active = bound > rtol
   end function xc_grad_batch_active

   !
   ! One batch: collocate (order 1, or order 1+2 for a GGA), gather D onto
   ! the local basis, X = chi D, the per-point functional and force
   ! tiles, the force GEMMs, and the atom scatter. Every array here is
   ! sized to this batch exactly -- see the module header on why there is
   ! no padding to match.
   !
   subroutine xc_grad_one_batch(b, lp, is_gga, nspin, nk, kid, coef, dens, xg, ib, npb, nloc, &
                                ao_atom, grad)
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      logical, intent(in) :: is_gga
      integer, intent(in) :: nspin, nk
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: dens(b%nao, b%nao, nspin)
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: ib, npb, nloc
      integer, intent(in) :: ao_atom(b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)

      real(dp), allocatable :: chi(:, :), gchi(:, :, :), hchi(:, :, :)
      real(dp), allocatable :: dloc(:, :, :), xmat(:, :, :)
      real(dp), allocatable :: pall(:, :, :), q1all(:, :, :), q2all(:, :, :), rall(:, :, :)
      real(dp), allocatable :: tall(:, :, :)
      integer :: s

      allocate (chi(npb, nloc), gchi(npb, nloc, 3))
      allocate (dloc(nloc, nloc, nspin), xmat(npb, nloc, nspin))
      allocate (pall(npb, 3*nloc, nspin), tall(3*nloc, nloc, nspin))
      !$acc enter data create(chi, gchi, dloc, xmat, pall, tall)

      if (is_gga) then
         allocate (hchi(npb, nloc, NHESS))
         allocate (q1all(npb, 3*nloc, nspin), q2all(npb, 3*nloc, nspin), rall(npb, nloc, nspin))
         !$acc enter data create(hchi, q1all, q2all, rall)
         call xc_grad_collocate2(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                                 xg%npts, xg%r, xg%nbatch, xg%b_off, xg%b_shoff, size(xg%b_sh), &
                                 xg%b_sh, xg%b_shao, ib, npb, nloc, chi, gchi, hchi)
      else
         call xc_grad_collocate1(b%nshell, b%maxnp, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, &
                                 xg%npts, xg%r, xg%nbatch, xg%b_off, xg%b_shoff, size(xg%b_sh), &
                                 xg%b_sh, xg%b_shao, ib, npb, nloc, chi, gchi)
      end if

      call xc_grad_gather(b%nao, size(xg%b_ao), xg%b_ao, xg%b_aooff(ib), nloc, nspin, dens, dloc)

      do s = 1, nspin
         call lp%gemm_mat(chi, dloc(:, :, s), xmat(:, :, s))
      end do

      if (is_gga) then
         call xc_grad_point_gga(npb, xg%npts, xg%w, xg%b_off(ib), nloc, nspin, nk, kid, coef, &
                                chi, gchi, hchi, xmat, pall, q1all, q2all, rall)
      else
         call xc_grad_point_lda(npb, xg%npts, xg%w, xg%b_off(ib), nloc, nspin, nk, kid, coef, &
                                chi, gchi, xmat, pall)
      end if

      do s = 1, nspin
         call lp%gemm_mat(pall(:, :, s), chi, tall(:, :, s), transa="T", beta=0.0_dp)
      end do
      if (is_gga) then
         do s = 1, nspin
            call lp%gemm_mat(q1all(:, :, s), rall(:, :, s), tall(:, :, s), transa="T", beta=1.0_dp)
         end do
         if (nspin == 2) then
            call lp%gemm_mat(q2all(:, :, 1), rall(:, :, 2), tall(:, :, 1), transa="T", beta=1.0_dp)
            call lp%gemm_mat(q2all(:, :, 2), rall(:, :, 1), tall(:, :, 2), transa="T", beta=1.0_dp)
         end if
      end if

      call xc_grad_scatter(ib, xg%nbatch, xg%b_aooff, size(xg%b_ao), xg%b_ao, ao_atom, b%nao, &
                           nloc, nspin, dloc, tall, b%natm, grad)

      !$acc exit data delete(chi, gchi, dloc, xmat, pall, tall)
      deallocate (chi, gchi, dloc, xmat, pall, tall)
      if (is_gga) then
         !$acc exit data delete(hchi, q1all, q2all, rall)
         deallocate (hchi, q1all, q2all, rall)
      end if
   end subroutine xc_grad_one_batch

   !
   ! Collocation, order 1: one thread per (local shell, local point) of
   ! this one batch. Same mapping trc_xc.F90's xc_collocate_body uses
   ! (`b_shao` is the shell's first LOCAL function, 0-based), just over one
   ! batch's own (ks, gl) pair instead of a flattened multi-batch work list.
   !
   subroutine xc_grad_collocate1(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                 nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, chi, gchi)
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3)
      integer :: nsh_b, ks, gl

      nsh_b = b_shoff(ib + 1) - b_shoff(ib)
      do concurrent(ks=1:nsh_b, gl=1:npb)
         call xc_grad_colloc1_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                   nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                   ks, gl, chi, gchi)
      end do
   end subroutine xc_grad_collocate1

   pure subroutine xc_grad_colloc1_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                        nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                        ks, gl, chi, gchi)
      !$acc routine seq
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc, ks, gl
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3)
      integer :: ish, iao, l, nc, m, gg
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX)

      ish = b_sh(b_shoff(ib) + ks - 1)
      iao = b_shao(b_shoff(ib) + ks - 1)
      gg = b_off(ib) + gl - 1
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call shell_collocate(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs)
      nc = (l + 1)*(l + 2)/2
      do m = 1, nc
         chi(gl, iao + m) = cs(m)
         gchi(gl, iao + m, 1) = gs(1, m)
         gchi(gl, iao + m, 2) = gs(2, m)
         gchi(gl, iao + m, 3) = gs(3, m)
      end do
   end subroutine xc_grad_colloc1_body

   !
   ! As `xc_grad_collocate1`, order 1 and 2 (the GGA Hessian term).
   !
   subroutine xc_grad_collocate2(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                 nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                 chi, gchi, hchi)
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      integer :: nsh_b, ks, gl

      nsh_b = b_shoff(ib + 1) - b_shoff(ib)
      do concurrent(ks=1:nsh_b, gl=1:npb)
         call xc_grad_colloc2_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                   nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                   ks, gl, chi, gchi, hchi)
      end do
   end subroutine xc_grad_collocate2

   pure subroutine xc_grad_colloc2_body(nshell, maxnp, sh_l, sh_np, sh_e, sh_c, sh_r, npts, r, &
                                        nbatch, b_off, b_shoff, nbsh, b_sh, b_shao, ib, npb, nloc, &
                                        ks, gl, chi, gchi, hchi)
      !$acc routine seq
      integer, intent(in) :: nshell, maxnp, npts, nbatch, nbsh, ib, npb, nloc, ks, gl
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell), r(3, npts)
      integer, intent(in) :: b_off(nbatch + 1), b_shoff(nbatch + 1), b_sh(nbsh), b_shao(nbsh)
      real(dp), intent(inout) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      integer :: ish, iao, l, nc, m, gg
      real(dp) :: d(3), cs(NCART_MAX), gs(3, NCART_MAX), hs(NHESS, NCART_MAX)

      ish = b_sh(b_shoff(ib) + ks - 1)
      iao = b_shao(b_shoff(ib) + ks - 1)
      gg = b_off(ib) + gl - 1
      l = sh_l(ish)
      d(1) = r(1, gg) - sh_r(1, ish)
      d(2) = r(2, gg) - sh_r(2, ish)
      d(3) = r(3, gg) - sh_r(3, ish)
      call shell_collocate2(l, sh_np(ish), sh_e(1:sh_np(ish), ish), sh_c(1:sh_np(ish), ish), d, cs, gs, hs)
      nc = (l + 1)*(l + 2)/2
      do m = 1, nc
         chi(gl, iao + m) = cs(m)
         gchi(gl, iao + m, 1) = gs(1, m)
         gchi(gl, iao + m, 2) = gs(2, m)
         gchi(gl, iao + m, 3) = gs(3, m)
         hchi(gl, iao + m, 1) = hs(1, m)
         hchi(gl, iao + m, 2) = hs(2, m)
         hchi(gl, iao + m, 3) = hs(3, m)
         hchi(gl, iao + m, 4) = hs(4, m)
         hchi(gl, iao + m, 5) = hs(5, m)
         hchi(gl, iao + m, 6) = hs(6, m)
      end do
   end subroutine xc_grad_colloc2_body

   !
   ! D sliced onto this batch's local basis, every (mu, nu, spin) of it --
   ! no padding, no screening within the batch, unlike xc_gather_body: this
   ! module's tiles are exact-sized, so there is no padded element to zero.
   !
   subroutine xc_grad_gather(nao, nbao, b_ao, aooff0, nloc, nspin, dens, dloc)
      integer, intent(in) :: nao, nbao, aooff0, nloc, nspin
      integer, intent(in) :: b_ao(nbao)
      real(dp), intent(in) :: dens(nao, nao, nspin)
      real(dp), intent(inout) :: dloc(nloc, nloc, nspin)
      integer :: mu, nu, s
      do concurrent(mu=1:nloc, nu=1:nloc, s=1:nspin)
         dloc(mu, nu, s) = dens(b_ao(aooff0 + mu - 1), b_ao(aooff0 + nu - 1), s)
      end do
   end subroutine xc_grad_gather

   !
   ! LDA: one thread per point. Builds rho, grad rho (needed even at LDA,
   ! since the moving function's OWN gradient enters the Pulay term), the
   ! potential, and the "P" tile `pall(g, (k-1)*nloc+mu, s)` the force
   ! GEMM `pall^T @ chi` contracts into T_k(mu,nu). Restricted and
   ! unrestricted in the same body, as trc_xc.F90's own xc_point_body does.
   !
   subroutine xc_grad_point_lda(npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, xmat, pall)
      integer, intent(in) :: npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin)
      integer :: g
      do concurrent(g=1:npb)
         call xc_grad_point_lda_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, xmat, pall)
      end do
   end subroutine xc_grad_point_lda

   pure subroutine xc_grad_point_lda_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, &
                                          xmat, pall)
      !$acc routine seq
      integer, intent(in) :: g, npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin)
      real(dp) :: rho(2), grho(3, 2)
      real(dp) :: sigma, eps, vrho, vsigma, wg
      real(dp) :: va, vb, s_aa, s_ab, s_bb, v_aa, v_ab, v_bb
      integer :: s, u, mu, k

      wg = w(p0 + g - 1)
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            rho(s) = rho(s) + xmat(g, u, s)*chi(g, u)
            grho(1, s) = grho(1, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 3)
         end do
      end do

      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point(nk, kid, coef, rho(1), sigma, eps, vrho, vsigma)
         do mu = 1, nloc
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*vrho*gchi(g, mu, k)
            end do
         end do
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, &
                                  eps, va, vb, v_aa, v_ab, v_bb)
         do mu = 1, nloc
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*va*gchi(g, mu, k)
               pall(g, (k - 1)*nloc + mu, 2) = wg*vb*gchi(g, mu, k)
            end do
         end do
      end if
   end subroutine xc_grad_point_lda_body

   !
   ! GGA: as the LDA point kernel, plus the Hessian-dependent "P" piece and
   ! the "Q"/"R" pair the second force GEMM (`q^T @ r`) contracts. R_nu(g)
   ! = grad_rho(g) . grad_chi_nu(g) is the ONE quantity shared by every mu
   ! on this point -- built once per nu, not recomputed per (mu, nu) pair.
   !
   ! Unrestricted carries two cross terms the restricted branch does not:
   ! channel a's P and Q1/Q2 read v_ab against channel b's Hessian/R and
   ! vice versa (Ca_k = Pa_k phi_nu + Qa1_k Ra_nu + Qa2_k Rb_nu, this
   ! module's header), so `q2all(:,:,1)` pairs with `rall(:,:,2)` and
   ! `q2all(:,:,2)` with `rall(:,:,1)` in the caller's GEMM step.
   !
   subroutine xc_grad_point_gga(npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, hchi, xmat, &
                                pall, q1all, q2all, rall)
      integer, intent(in) :: npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      real(dp), intent(in) :: xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin), q1all(npb, 3*nloc, nspin)
      real(dp), intent(inout) :: q2all(npb, 3*nloc, nspin), rall(npb, nloc, nspin)
      integer :: g
      do concurrent(g=1:npb)
         call xc_grad_point_gga_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, hchi, &
                                     xmat, pall, q1all, q2all, rall)
      end do
   end subroutine xc_grad_point_gga

   pure subroutine xc_grad_point_gga_body(g, npb, npts, w, p0, nloc, nspin, nk, kid, coef, chi, gchi, &
                                          hchi, xmat, pall, q1all, q2all, rall)
      !$acc routine seq
      integer, intent(in) :: g, npb, npts, p0, nloc, nspin, nk
      real(dp), intent(in) :: w(npts)
      integer, intent(in) :: kid(XC_MAX_KERNELS)
      real(dp), intent(in) :: coef(XC_MAX_KERNELS)
      real(dp), intent(in) :: chi(npb, nloc), gchi(npb, nloc, 3), hchi(npb, nloc, NHESS)
      real(dp), intent(in) :: xmat(npb, nloc, nspin)
      real(dp), intent(inout) :: pall(npb, 3*nloc, nspin), q1all(npb, 3*nloc, nspin)
      real(dp), intent(inout) :: q2all(npb, 3*nloc, nspin), rall(npb, nloc, nspin)
      real(dp) :: rho(2), grho(3, 2)
      real(dp) :: sigma, eps, vrho, vsigma, wg
      real(dp) :: va, vb, s_aa, s_ab, s_bb, v_aa, v_ab, v_bb
      real(dp) :: hdota(3), hdotb(3)
      integer :: s, u, mu, k

      wg = w(p0 + g - 1)
      rho = 0.0_dp
      grho = 0.0_dp
      do s = 1, nspin
         do u = 1, nloc
            rho(s) = rho(s) + xmat(g, u, s)*chi(g, u)
            grho(1, s) = grho(1, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 1)
            grho(2, s) = grho(2, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 2)
            grho(3, s) = grho(3, s) + 2.0_dp*xmat(g, u, s)*gchi(g, u, 3)
         end do
      end do

      if (nspin == 1) then
         sigma = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         call xc_eval_point(nk, kid, coef, rho(1), sigma, eps, vrho, vsigma)
         do mu = 1, nloc
            ! (grad_rho . H_mu)_k, hchi packed xx, xy, xz, yy, yz, zz.
            hdota(1) = grho(1, 1)*hchi(g, mu, 1) + grho(2, 1)*hchi(g, mu, 2) + grho(3, 1)*hchi(g, mu, 3)
            hdota(2) = grho(1, 1)*hchi(g, mu, 2) + grho(2, 1)*hchi(g, mu, 4) + grho(3, 1)*hchi(g, mu, 5)
            hdota(3) = grho(1, 1)*hchi(g, mu, 3) + grho(2, 1)*hchi(g, mu, 5) + grho(3, 1)*hchi(g, mu, 6)
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*(vrho*gchi(g, mu, k) + 2.0_dp*vsigma*hdota(k))
               q1all(g, (k - 1)*nloc + mu, 1) = wg*2.0_dp*vsigma*gchi(g, mu, k)
            end do
         end do
         do u = 1, nloc
            rall(g, u, 1) = grho(1, 1)*gchi(g, u, 1) + grho(2, 1)*gchi(g, u, 2) + grho(3, 1)*gchi(g, u, 3)
         end do
      else
         s_aa = grho(1, 1)**2 + grho(2, 1)**2 + grho(3, 1)**2
         s_bb = grho(1, 2)**2 + grho(2, 2)**2 + grho(3, 2)**2
         s_ab = grho(1, 1)*grho(1, 2) + grho(2, 1)*grho(2, 2) + grho(3, 1)*grho(3, 2)
         call xc_eval_point_polar(nk, kid, coef, rho(1), rho(2), s_aa, s_ab, s_bb, &
                                  eps, va, vb, v_aa, v_ab, v_bb)
         do mu = 1, nloc
            hdota(1) = grho(1, 1)*hchi(g, mu, 1) + grho(2, 1)*hchi(g, mu, 2) + grho(3, 1)*hchi(g, mu, 3)
            hdota(2) = grho(1, 1)*hchi(g, mu, 2) + grho(2, 1)*hchi(g, mu, 4) + grho(3, 1)*hchi(g, mu, 5)
            hdota(3) = grho(1, 1)*hchi(g, mu, 3) + grho(2, 1)*hchi(g, mu, 5) + grho(3, 1)*hchi(g, mu, 6)
            hdotb(1) = grho(1, 2)*hchi(g, mu, 1) + grho(2, 2)*hchi(g, mu, 2) + grho(3, 2)*hchi(g, mu, 3)
            hdotb(2) = grho(1, 2)*hchi(g, mu, 2) + grho(2, 2)*hchi(g, mu, 4) + grho(3, 2)*hchi(g, mu, 5)
            hdotb(3) = grho(1, 2)*hchi(g, mu, 3) + grho(2, 2)*hchi(g, mu, 5) + grho(3, 2)*hchi(g, mu, 6)
            do k = 1, 3
               pall(g, (k - 1)*nloc + mu, 1) = wg*(va*gchi(g, mu, k) + 2.0_dp*v_aa*hdota(k) + v_ab*hdotb(k))
               pall(g, (k - 1)*nloc + mu, 2) = wg*(vb*gchi(g, mu, k) + 2.0_dp*v_bb*hdotb(k) + v_ab*hdota(k))
               q1all(g, (k - 1)*nloc + mu, 1) = wg*2.0_dp*v_aa*gchi(g, mu, k)   !! x rall(:,:,1)
               q2all(g, (k - 1)*nloc + mu, 1) = wg*v_ab*gchi(g, mu, k)         !! x rall(:,:,2)
               q1all(g, (k - 1)*nloc + mu, 2) = wg*2.0_dp*v_bb*gchi(g, mu, k)   !! x rall(:,:,2)
               q2all(g, (k - 1)*nloc + mu, 2) = wg*v_ab*gchi(g, mu, k)         !! x rall(:,:,1)
            end do
         end do
         do u = 1, nloc
            rall(g, u, 1) = grho(1, 1)*gchi(g, u, 1) + grho(2, 1)*gchi(g, u, 2) + grho(3, 1)*gchi(g, u, 3)
            rall(g, u, 2) = grho(1, 2)*gchi(g, u, 1) + grho(2, 2)*gchi(g, u, 2) + grho(3, 2)*gchi(g, u, 3)
         end do
      end if
   end subroutine xc_grad_point_gga_body

   !
   ! T_k(mu, nu) = tall((k-1)*nloc+mu, nu, s), contracted with D and summed
   ! over spin, lands as -2 * (that sum) on the atom owning local function
   ! mu -- one thread per (mu, k) of this batch, atomic because two batches
   ! (or two points of the atom-scatter loop below) can own the same atom.
   ! Same #ifdef pattern as trc_xc.F90's xc_scatter, for the same reason:
   ! a host compiler that takes `do concurrent` at its word reorders the
   ! read-modify-writes where the device's atomic makes them correct.
   !
   subroutine xc_grad_scatter(ib, nbatch, b_aooff, nbao, b_ao, ao_atom, nao, nloc, nspin, &
                              dloc, tall, natm, grad)
      integer, intent(in) :: ib, nbatch, nbao, nao, nloc, nspin, natm
      integer, intent(in) :: b_aooff(nbatch + 1), b_ao(nbao), ao_atom(nao)
      real(dp), intent(in) :: dloc(nloc, nloc, nspin), tall(3*nloc, nloc, nspin)
      real(dp), intent(inout) :: grad(3, natm)
      integer :: mu, k, at
      real(dp) :: fval
#ifdef _OPENACC
      do concurrent(mu=1:nloc, k=1:3) local(fval, at)
         call xc_grad_scatter_body(mu, k, nloc, nspin, dloc, tall, fval)
         at = ao_atom(b_ao(b_aooff(ib) + mu - 1))
         !$acc atomic update
         grad(k, at) = grad(k, at) + fval
      end do
#else
      do mu = 1, nloc
         do k = 1, 3
            call xc_grad_scatter_body(mu, k, nloc, nspin, dloc, tall, fval)
            at = ao_atom(b_ao(b_aooff(ib) + mu - 1))
            grad(k, at) = grad(k, at) + fval
         end do
      end do
#endif
   end subroutine xc_grad_scatter

   pure subroutine xc_grad_scatter_body(mu, k, nloc, nspin, dloc, tall, fval)
      !$acc routine seq
      integer, intent(in) :: mu, k, nloc, nspin
      real(dp), intent(in) :: dloc(nloc, nloc, nspin), tall(3*nloc, nloc, nspin)
      real(dp), intent(out) :: fval
      integer :: nu, s
      real(dp) :: acc
      acc = 0.0_dp
      do s = 1, nspin
         do nu = 1, nloc
            acc = acc + dloc(mu, nu, s)*tall((k - 1)*nloc + mu, nu, s)
         end do
      end do
      fval = -2.0_dp*acc
   end subroutine xc_grad_scatter_body

end module trc_xc_grad
