module trc_xc_grad_response
   !!
   !! The exchange-correlation part of the excited-state (TDDFT) nuclear
   !! gradient, for many roots in one sweep of the grid: for each column I,
   !! the nuclear derivative at FIXED AO matrices of
   !!
   !!     G_I = Tr(P_I V_xc[D]) + 1/2 E_xc''[T_I, T_I]
   !!
   !! with D the reference density, P_I a relaxed or difference density, T_I
   !! a transition density, and, restricted, f(rho, sigma) the functional of
   !! the TOTAL density rho:
   !!
   !!     Tr(P V_xc[D]) = sum_g w_g [ f_rho rho_P + 2 f_sigma grad rho . grad rho_P ]
   !!     E_xc''[T, T]  = sum_g w_g [ f_rr rho_T^2 + 2 f_rs rho_T sigma'
   !!                                 + f_ss sigma'^2 + 2 f_s |grad rho_T|^2 ],
   !!                     sigma' = 2 grad rho . grad rho_T,
   !!
   !! rho_X = sum_uv sym(X)_uv chi_u chi_v. E'' is the quantity
   !! trc_xc_fxc_rks contracts: 1/2 E''[T, T] = 1/2 Tr(T trc_xc_fxc_rks[T]).
   !! The TDDFT gradient of Furche and Ahlrichs (J. Chem. Phys. 117, 7433
   !! (2002)) needs exactly these two pieces; the CALLER owns the spin factor,
   !! which is quadratic in T: pass the T for which Tr(T f_xc[T]) is the
   !! response energy the same solve used (the T its trc_xc_fxc_rks_apply_many
   !! call was given). Only the symmetric part of P and T matters (V_xc and f_xc
   !! are symmetric operators); both are symmetrised here.
   !!
   !! WHAT MOVES. The same as xc_grad_potential_add and trc_xc_grad: the basis
   !! functions (Pulay), with `options%grid_response` the grid points with
   !! their owner and the Becke weights (trc_xc_grad's header, GRID RESPONSE),
   !! the grid's own response being the integrand times dP/dR. The matrices
   !! D, P_I and T_I do not move. Collected by which density's basis
   !! derivative multiplies what, per column (a drho + b . d grad rho, with
   !! kernels at D):
   !!
   !!   on P:  a = f_rho,                        b = 2 f_sigma grad rho
   !!   on T:  a = f_rr rho_T + f_rs sigma',     b = 2 (f_rs rho_T + f_ss sigma') grad rho
   !!                                                + 2 f_s grad rho_T
   !!   on D:  the P term's  a = f_rr rho_P + 2 f_rs grad rho . grad rho_P,
   !!                        b = 2 f_rs rho_P grad rho + 4 f_ss (grad rho . grad rho_P) grad rho
   !!                            + 2 f_s grad rho_P
   !!          plus the T term's k_xc part: the derivative of
   !!          Psi = 1/2 (a_T rho_T + b_T . grad rho_T) with respect to rho and
   !!          grad rho at fixed T, by central difference of the analytic f_xc
   !!          (trc_xc_kxc_rks_apply's step rule; trc_xc_grad_coef_rks).
   !!
   !! The T channel's coefficients are the T-derivative of Psi, which is f_xc[T]
   !! itself, and the P and D ones those of xc_grad_potential_add, so with T = 0
   !! a column is exactly xc_grad_potential_add's.
   !!
   !! ONE SWEEP. The force sweep (xc_grad_drive_channels) collocates chi, its
   !! gradient and Hessian once per chunk of batches for every column in the
   !! sweep, and forms X = chi D_s per density slot, R and RD per channel; the
   !! coefficient sweep (trc_xc_grad_coef_rks) is one more pass, without
   !! Hessians. The partition's response is one dft_partition_becke_grad_contract
   !! per column (the Becke derivatives are recomputed for each).
   !!
   !! SCRATCH. With npts points, tile = the largest padded batch tile
   !! (npp x nlp) and nc columns in a sweep, the device holds
   !! 13 npts nc doubles of coefficients and energies plus
   !! (11 + 8 nc) tiles of at most tcap doubles each, plus the pair tiles.
   !! `budget` (doubles; default 32 XC_GRAD_TILE_CAP) bounds the sum: nc is the
   !! most columns that fit with the largest single tile, and tcap shrinks
   !! as nc grows (at least one batch per chunk). More columns than fit are
   !! processed in further sweeps. trc_xc_grad_response_plan returns both.
   !!
   !! LDA and GGA, restricted, Cartesian AO. The refusals are
   !! trc_xc_response_refuse_reason's. Validated by test/check_xc_grad_response.F90.
   !!
   use pic_types, only: dp, int64
   use trc_basis_types, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc, only: trc_xc_grad_coef_rks, trc_xc_grad_coef_t, trc_xc_response_refuse_reason, XC_NP_PAD, XC_NL_PAD
   use trc_xc_grad, only: xc_grad_drive_channels, xc_grad_options_t, xc_grad_step_names, xc_grad_batch_active, &
                          XC_GRAD_TILE_CAP
   use trc_dft_partition, only: dft_partition_becke_grad_contract
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: trc_xc_grad_response_many, trc_xc_grad_response_plan

contains

   subroutine trc_xc_grad_response_plan(b, xg, ncol, budget, nc, tcap)
      !! The sweep shape trc_xc_grad_response_many uses for `ncol` columns at
      !! `budget` doubles (absent: 32 XC_GRAD_TILE_CAP): `nc` columns per sweep
      !! (>= 1) and `tcap`, the doubles one tile of a chunk may hold.
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: ncol
      integer(int64), intent(in), optional :: budget
      integer, intent(out) :: nc
      integer(int64), intent(out) :: tcap

      integer(int64) :: cap, tmin, fixed, per_col, npb, nloc
      integer :: ib

      cap = 32_int64*XC_GRAD_TILE_CAP
      if (present(budget)) cap = max(budget, 1_int64)
      tmin = 1
      do ib = 1, xg%nbatch
         npb = ((int(xg%b_off(ib + 1) - xg%b_off(ib), int64) + XC_NP_PAD - 1)/XC_NP_PAD)*XC_NP_PAD
         nloc = ((int(xg%b_aooff(ib + 1) - xg%b_aooff(ib), int64) + XC_NL_PAD - 1)/XC_NL_PAD)*XC_NL_PAD
         tmin = max(tmin, npb*nloc)
      end do
      fixed = 11_int64*tmin
      per_col = 8_int64*tmin + 13_int64*int(max(xg%npts, 1), int64)
      nc = int(max(1_int64, min(int(max(ncol, 1), int64), (cap - fixed)/per_col)))
      tcap = (cap - 13_int64*int(max(xg%npts, 1), int64)*nc)/(11_int64 + 8_int64*nc)
      tcap = max(tmin, min(tcap, XC_GRAD_TILE_CAP))
      tcap = max(tcap, tmin)
   end subroutine trc_xc_grad_response_plan

   subroutine trc_xc_grad_response_many(b, lp, func, xg, dmat, pmats, tmats, ncol, grad, options, err, budget, nsweep, no_kxc)
      !!
      !! grad(:, :, I) += d/dR [ Tr(P_I V_xc[D]) + 1/2 E_xc''[T_I, T_I] ],
      !! I = 1..ncol, D, P_I, T_I held fixed (this module's header has the
      !! terms and what moves).
      !!
      !! `dmat` (nao, nao) is the reference density, `pmats` and `tmats`
      !! (nao, nao, ncol) the columns' densities, all HOST arrays in the
      !! Cartesian AO space; P and T are symmetrised on entry. `xg` is the
      !! reference's grid, device-resident and built with its atoms and
      !! quadrature weights (`build_from_grid`); `grad` (3, natm, ncol) is
      !! device-resident (as xc_grad_potential_add's) and ADDED into.
      !! `options`: density screen, grid response and the grid's partition, as
      !! xc_grad_add_on_grid's. A batch the screen drops on D carries no term,
      !! for T as well. `err` carries trc_xc_response_refuse_reason and
      !! `grad` is then untouched. `nsweep`, if given, is the number of grid
      !! sweeps made (columns are split by `budget`, trc_xc_grad_response_plan).
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: ncol
      real(dp), intent(in) :: dmat(b%nao, b%nao), pmats(b%nao, b%nao, ncol), tmats(b%nao, b%nao, ncol)
      real(dp), intent(inout) :: grad(3, b%natm, ncol)
      type(xc_grad_options_t), intent(in) :: options
      type(error_t), intent(out) :: err
      integer(int64), intent(in), optional :: budget
      integer, intent(out), optional :: nsweep
      logical, intent(in), optional :: no_kxc
      !! Diagnostic: omit the k_xc part (the negative control of the tests).

      type(trc_xc_grad_coef_t) :: pts
      type(step_table_t) :: timing
      character(len=256) :: why
      real(dp), allocatable :: dpack(:, :, :), dsym1(:, :, :), gsw(:, :, :), ghost(:, :, :), grad_extra(:, :), pcoef(:)
      integer, allocatable :: active(:), at_elem(:)
      integer(int64) :: tcap
      integer :: n, ib, nc, c0, nci, k, icol, nsw, npts, nslot

      if (present(nsweep)) nsweep = 0
      why = trc_xc_response_refuse_reason(func, "trc_xc_grad_response")
      if (len_trim(why) > 0) then
         call err%set(ERROR_VALIDATION, trim(why))
         return
      end if
      if (xg%npts == 0 .or. ncol < 1) return
      n = b%nao
      npts = xg%npts
      call trc_xc_grad_response_plan(b, xg, ncol, budget, nc, tcap)

      ! The reference's batch screen (xc_grad_potential_add's): host, on D.
      allocate (dsym1(n, n, 1), active(xg%nbatch))
      dsym1(:, :, 1) = 0.5_dp*(dmat + transpose(dmat))
      do ib = 1, xg%nbatch
         active(ib) = 0
         if (xg%b_off(ib + 1) > xg%b_off(ib) .and. xg%b_aooff(ib + 1) > xg%b_aooff(ib)) then
            if (xc_grad_batch_active(ib, xg, n, 1, dsym1, options%rho_tol)) active(ib) = 1
         end if
      end do
      call timing%init("xc response grad", xc_grad_step_names())
      allocate (ghost(3, b%natm, ncol), source=0.0_dp)
      allocate (at_elem(b%natm), source=b%at_elem)
      if (options%grid_response) allocate (grad_extra(3, b%natm), pcoef(npts))

      nsw = 0
      do c0 = 1, ncol, nc
         nci = min(nc, ncol - c0 + 1)
         nslot = 1 + 2*nci
         allocate (dpack(n, n, nslot), gsw(3, b%natm, nci))
         dpack(:, :, 1) = dsym1(:, :, 1)
         do k = 1, nci
            dpack(:, :, 1 + k) = 0.5_dp*(pmats(:, :, c0 + k - 1) + transpose(pmats(:, :, c0 + k - 1)))
            dpack(:, :, 1 + nci + k) = 0.5_dp*(tmats(:, :, c0 + k - 1) + transpose(tmats(:, :, c0 + k - 1)))
         end do
         gsw = 0.0_dp
         !$acc enter data copyin(dpack, gsw)
         call trc_xc_grad_coef_rks(b, xg, func, dpack, nci, active, tcap, pts, lp, no_kxc)
         call xc_grad_drive_channels(b, lp, xg, dpack, nci, pts%cf, active, tcap, options%grid_response, gsw, timing)
         nsw = nsw + 1
         !$acc update self(gsw)
         ghost(:, :, c0:c0 + nci - 1) = ghost(:, :, c0:c0 + nci - 1) + gsw

         if (options%grid_response) then
            ! The partition's response, both pieces (trc_xc_grad's header, B
            ! and C), per column from the integrand times the quadrature
            ! weight alone.
            do k = 1, nci
               icol = c0 + k - 1
               grad_extra = 0.0_dp
               pcoef = xg%qw
               !$acc enter data copyin(pcoef)
               call xc_grad_response_scale(npts, nci, k, pts%energy, pcoef)
               call dft_partition_becke_grad_contract(xg%r, b%at_r, at_elem, xg%atom, pcoef, options%partition, &
                                                      grad_extra, err)
               !$acc exit data delete(pcoef)
               if (err%has_error()) exit
               ghost(:, :, icol) = ghost(:, :, icol) + grad_extra
            end do
         end if
         !$acc exit data delete(dpack, gsw)
         call pts%release()
         deallocate (dpack, gsw)
         if (err%has_error()) return
      end do

      !$acc update self(grad)
      grad = grad + ghost
      !$acc update device(grad)
      if (present(nsweep)) nsweep = nsw
   end subroutine trc_xc_grad_response_many

   subroutine xc_grad_response_scale(npts, ncol, icol, energy, pcoef)
      !! pcoef(g) := energy(g, icol) pcoef(g).
      integer, intent(in) :: npts, ncol, icol
      real(dp), intent(in) :: energy(npts, ncol)
      real(dp), intent(inout) :: pcoef(npts)
      integer :: g
      do concurrent(g=1:npts)
         pcoef(g) = energy(g, icol)*pcoef(g)
      end do
   end subroutine xc_grad_response_scale

end module trc_xc_grad_response
