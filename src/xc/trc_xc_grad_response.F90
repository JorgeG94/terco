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
   !! SCREENING. None by default. trc_xc_fxc_rks and trc_xc_kxc_rks_apply, which
   !! define the T term and drive the response solve, are unscreened, and the
   !! gradient differentiates the same operator, so every batch carrying
   !! points and functions is kept. `options%screen_tol` > 0 applies
   !! trc_xc_grad's density screen on D at that threshold instead (as
   !! xc_grad_potential_add does at 1e-12); `options%gradient%rho_tol` is NOT read. Measured on water/6-31G*
   !! (test/check_xc_grad_response.F90, SVWN): a 1e-12
   !! screen changes a column by 1e-14, a 1e-4 screen by about 1e-10, and a
   !! 1e-1 screen by 2e-2, which the test asserts it can tell apart.
   !!
   !! SCRATCH. With npts points, nao orbitals and nc columns in a sweep, the
   !! device holds, in doubles:
   !!   13 npts nc            coefficients and energies (trc_xc_grad_coef_t)
   !!   nao^2 (1 + 2 nc)      the symmetrised density slots
   !!   (11 + 8 nc) T         chi, grad chi, Hessian, X per slot, R and RD per
   !!                         channel, T the largest chunk's tile
   !!   (1 + 2 nc) Pr         the gathered pair tiles, Pr the largest chunk's
   !!                         sum of nlp^2
   !!   3 nmu 32 nc           the Pulay partials (nmu the largest chunk's
   !!                         functions; 32 lanes on the device, 1 on the host)
   !! `budget` (doubles; default 32 XC_GRAD_TILE_CAP) bounds that sum:
   !! trc_xc_grad_response_plan takes the most columns per sweep for which a
   !! tile cap of at least one batch brings the sum under it, shrinking the
   !! tile cap (more chunks of batches per sweep) before the column count.
   !! More columns than fit are processed in further sweeps. A budget too
   !! small for even one column at one batch per chunk is exceeded, with a
   !! warning from the plan (pic_logger), never silently. The plan counts
   !! every batch as kept, so it is an upper bound.
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
                          XC_GRAD_TILE_CAP, XC_GRAD_LANES
   use trc_dft_partition, only: dft_partition_becke_grad_contract
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: trc_xc_grad_response_many, trc_xc_grad_response_plan, trc_xc_grad_response_options_t
   public :: xc_grad_response_core, xc_grad_response_hooks_t
   ! xc_grad_response_core: trc_xc_grad_response_many plus test hooks; not API.

   type :: trc_xc_grad_response_options_t
      !! The settings of trc_xc_grad_response_many, and its status on return.
      type(xc_grad_options_t) :: gradient
      !! Grid response and the grid's partition, as xc_grad_add_on_grid's. Its
      !! density screen is not read: `screen_tol`.
      integer(int64) :: budget = 0_int64
      !! Doubles bounding a sweep's scratch (this module's header, SCRATCH);
      !! zero is 32 XC_GRAD_TILE_CAP.
      real(dp) :: screen_tol = 0.0_dp
      !! Density screen on D, as xc_grad_batch_active; zero: none.
      integer :: nsweep = 0
      !! Out: the grid sweeps made. Reset to 0 on entry.
      type(error_t) :: err
      !! Out: trc_xc_response_refuse_reason's message, or the partition's error.
      !! Cleared on entry, so one options object can be reused after a refusal.
   end type trc_xc_grad_response_options_t

   type, extends(trc_xc_grad_response_options_t) :: xc_grad_response_hooks_t
      !! The options of xc_grad_response_core: the test hooks on top.
      logical :: no_kxc = .false.
      !! Omit the k_xc part of the D channel (the negative control).
      integer(int64) :: tile_cap = 0_int64
      !! Replace the plan's tile cap (> 0), so two calls with different column
      !! counts share one chunking.
   end type xc_grad_response_hooks_t

contains

   subroutine xc_grad_response_chunks(xg, npp, nlp, tcap, nchunk, tile, pairs, nmu)
      !! The chunks xc_grad_drive_channels makes at tile cap `tcap` with every
      !! batch kept: their number, and the largest chunk's tile, pair-tile and
      !! function counts.
      type(trc_xc_grid_t), intent(in) :: xg
      integer(int64), intent(in) :: tcap
      integer(int64), intent(in) :: npp(xg%nbatch), nlp(xg%nbatch)
      integer, intent(out) :: nchunk
      integer(int64), intent(out) :: tile, pairs, nmu
      integer(int64) :: need, pr, mu
      integer :: ib, b0

      nchunk = 0
      tile = 1
      pairs = 1
      nmu = 1
      b0 = 1
      do while (b0 <= xg%nbatch)
         need = 0
         pr = 0
         mu = 0
         ib = b0
         do while (ib <= xg%nbatch)
            if (need + npp(ib)*nlp(ib) > tcap .and. ib > b0) exit
            need = need + npp(ib)*nlp(ib)
            pr = pr + nlp(ib)*nlp(ib)
            mu = mu + (xg%b_aooff(ib + 1) - xg%b_aooff(ib))
            ib = ib + 1
         end do
         nchunk = nchunk + 1
         tile = max(tile, need)
         pairs = max(pairs, pr)
         nmu = max(nmu, mu)
         b0 = ib
      end do
   end subroutine xc_grad_response_chunks

   subroutine trc_xc_grad_response_plan(b, xg, ncol, budget, nc, tcap, nchunk, need)
      !! The sweep shape trc_xc_grad_response_many uses for `ncol` columns at
      !! `budget` doubles (absent: 32 XC_GRAD_TILE_CAP; this module's header,
      !! SCRATCH): `nc` columns per sweep (>= 1), `tcap`, the doubles one tile
      !! of a chunk may hold, `nchunk` the chunks of batches a sweep makes at
      !! that cap, and `need` the doubles a sweep then holds. Warns when `need`
      !! exceeds the budget even for one column at one batch per chunk.
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      integer, intent(in) :: ncol
      integer(int64), intent(in), optional :: budget
      integer, intent(out) :: nc
      integer(int64), intent(out) :: tcap
      integer, intent(out), optional :: nchunk
      integer(int64), intent(out), optional :: need

      integer(int64), allocatable :: npp(:), nlp(:)
      integer(int64) :: cap, tmin, tc, tile, pairs, nmu, tot, base, nao2
      integer :: ib, nct, nch, ntry
      logical :: fits

      cap = 32_int64*XC_GRAD_TILE_CAP
      if (present(budget)) cap = max(budget, 1_int64)
      allocate (npp(max(xg%nbatch, 1)), nlp(max(xg%nbatch, 1)))
      tmin = 1
      do ib = 1, xg%nbatch
         npp(ib) = ((int(xg%b_off(ib + 1) - xg%b_off(ib), int64) + XC_NP_PAD - 1)/XC_NP_PAD)*XC_NP_PAD
         nlp(ib) = ((int(xg%b_aooff(ib + 1) - xg%b_aooff(ib), int64) + XC_NL_PAD - 1)/XC_NL_PAD)*XC_NL_PAD
         tmin = max(tmin, npp(ib)*nlp(ib))
      end do
      nao2 = int(b%nao, int64)**2
      fits = .false.
      nc = 1
      tcap = tmin
      nch = 0
      tot = 0
      do ntry = max(ncol, 1), 1, -1
         nct = ntry
         base = 13_int64*int(max(xg%npts, 1), int64)*nct + nao2*(1_int64 + 2_int64*nct)
         tc = (cap - base)/(11_int64 + 8_int64*nct + 1_int64 + 2_int64*nct)
         tc = max(tmin, min(tc, XC_GRAD_TILE_CAP))
         do
            call xc_grad_response_chunks(xg, npp, nlp, tc, nch, tile, pairs, nmu)
            tot = base + (11_int64 + 8_int64*nct)*tile + (1_int64 + 2_int64*nct)*pairs &
                  + 3_int64*nmu*XC_GRAD_LANES*nct
            if (tot <= cap) then
               fits = .true.
               exit
            end if
            if (tc <= tmin) exit
            tc = max(tmin, (tc*3_int64)/4_int64)
         end do
         if (fits .or. ntry == 1) then
            nc = nct
            tcap = tc
            exit
         end if
      end do
      if (.not. fits) then
         call global_logger%warning("trc_xc_grad_response_plan: budget of "//to_string(cap)//" doubles is too small " &
                                    //"for one column at one batch per chunk ("//to_string(tot)//" needed); exceeding it")
      end if
      if (present(nchunk)) nchunk = nch
      if (present(need)) need = tot
   end subroutine trc_xc_grad_response_plan

   subroutine trc_xc_grad_response_many(b, lp, func, xg, dmat, columns, grad, options)
      !!
      !! grad(:, :, I) += d/dR [ Tr(P_I V_xc[D]) + 1/2 E_xc''[T_I, T_I] ],
      !! I = 1..ncol, D, P_I, T_I held fixed (this module's header has the
      !! terms, what moves, the screening and the scratch).
      !!
      !! `dmat` (nao, nao) is the reference density and `columns`
      !! (nao, nao, 2, ncol) the columns' P (slice 1) and T (slice 2), all HOST
      !! arrays in the Cartesian AO space; P and T are symmetrised on entry.
      !! `xg` is the reference's grid, device-resident and built with its atoms
      !! and quadrature weights (`build_from_grid`); `grad` (3, natm, ncol) is
      !! device-resident (as xc_grad_potential_add's) and ADDED into.
      !! `options` holds the settings and, on return, `err`
      !! (trc_xc_response_refuse_reason, and `grad` is then untouched) and
      !! `nsweep`, the grid sweeps made (columns are split by the budget,
      !! trc_xc_grad_response_plan).
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: dmat(:, :)
      real(dp), intent(in) :: columns(:, :, :, :)
      real(dp), intent(inout) :: grad(:, :, :)
      type(trc_xc_grad_response_options_t), intent(inout) :: options
      type(xc_grad_response_hooks_t) :: hooks

      hooks%trc_xc_grad_response_options_t = options
      call xc_grad_response_core(b, lp, func, xg, dmat, columns, grad, hooks)
      options = hooks%trc_xc_grad_response_options_t
   end subroutine trc_xc_grad_response_many

   subroutine xc_grad_response_core(b, lp, func, xg, dmat, columns, grad, options)
      !! trc_xc_grad_response_many with the test hooks of
      !! xc_grad_response_hooks_t, which are not API.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: dmat(:, :)
      real(dp), intent(in) :: columns(:, :, :, :)
      real(dp), intent(inout) :: grad(:, :, :)
      type(xc_grad_response_hooks_t), intent(inout) :: options

      type(trc_xc_grad_coef_t) :: pts
      type(step_table_t) :: timing
      character(len=256) :: why
      real(dp), allocatable :: dpack(:, :, :), dsym1(:, :, :), gsw(:, :, :), ghost(:, :, :), grad_extra(:, :), pcoef(:)
      integer, allocatable :: active(:), at_elem(:)
      integer(int64) :: tcap
      integer :: n, ib, nc, c0, nci, k, icol, nsw, npts, nslot, ncol
      real(dp) :: tol
      integer(int64) :: budget

      ncol = size(columns, 4)
      options%nsweep = 0
      call options%err%clear()
      why = trc_xc_response_refuse_reason(func, "trc_xc_grad_response")
      if (len_trim(why) > 0) then
         call options%err%set(ERROR_VALIDATION, trim(why))
         return
      end if
      if (xg%npts == 0 .or. ncol < 1) return
      n = b%nao
      npts = xg%npts
      budget = options%budget
      if (budget <= 0) budget = 32_int64*XC_GRAD_TILE_CAP
      call trc_xc_grad_response_plan(b, xg, ncol, budget, nc, tcap)
      if (options%tile_cap > 0) tcap = options%tile_cap
      tol = options%screen_tol

      allocate (dsym1(n, n, 1), active(xg%nbatch))
      dsym1(:, :, 1) = 0.5_dp*(dmat + transpose(dmat))
      do ib = 1, xg%nbatch
         active(ib) = 0
         if (xg%b_off(ib + 1) > xg%b_off(ib) .and. xg%b_aooff(ib + 1) > xg%b_aooff(ib)) then
            active(ib) = 1
            if (tol > 0.0_dp) then
               if (.not. xc_grad_batch_active(ib, xg, n, 1, dsym1, tol)) active(ib) = 0
            end if
         end if
      end do
      call timing%init("xc response grad", xc_grad_step_names())
      allocate (ghost(3, b%natm, ncol), source=0.0_dp)
      allocate (at_elem(b%natm), source=b%at_elem)
      if (options%gradient%grid_response) allocate (grad_extra(3, b%natm), pcoef(npts))

      nsw = 0
      do c0 = 1, ncol, nc
         nci = min(nc, ncol - c0 + 1)
         nslot = 1 + 2*nci
         allocate (dpack(n, n, nslot), gsw(3, b%natm, nci))
         dpack(:, :, 1) = dsym1(:, :, 1)
         do k = 1, nci
            dpack(:, :, 1 + k) = 0.5_dp*(columns(:, :, 1, c0 + k - 1) + transpose(columns(:, :, 1, c0 + k - 1)))
            dpack(:, :, 1 + nci + k) = 0.5_dp*(columns(:, :, 2, c0 + k - 1) + transpose(columns(:, :, 2, c0 + k - 1)))
         end do
         gsw = 0.0_dp
         !$acc enter data copyin(dpack, gsw)
         pts%active = active
         pts%tcap = tcap
         pts%no_kxc = options%no_kxc
         call trc_xc_grad_coef_rks(b, xg, func, dpack, pts, lp)
         call xc_grad_drive_channels(b, lp, xg, dpack, pts, options%gradient%grid_response, gsw, timing)
         nsw = nsw + 1
         !$acc update self(gsw)
         ghost(:, :, c0:c0 + nci - 1) = ghost(:, :, c0:c0 + nci - 1) + gsw

         if (options%gradient%grid_response) then
            ! The partition's response, both pieces (trc_xc_grad's header, B
            ! and C), per column from the integrand times the quadrature
            ! weight alone.
            do k = 1, nci
               icol = c0 + k - 1
               grad_extra = 0.0_dp
               pcoef = xg%qw
               !$acc enter data copyin(pcoef)
               call xc_grad_response_scale(npts, nci, k, pts%energy, pcoef)
               call dft_partition_becke_grad_contract(xg%r, b%at_r, at_elem, xg%atom, pcoef, options%gradient%partition, &
                                                      grad_extra, options%err)
               !$acc exit data delete(pcoef)
               if (options%err%has_error()) exit
               ghost(:, :, icol) = ghost(:, :, icol) + grad_extra
            end do
         end if
         !$acc exit data delete(dpack, gsw)
         call pts%release()
         deallocate (dpack, gsw)
         if (options%err%has_error()) return
      end do

      !$acc update self(grad)
      grad = grad + ghost
      !$acc update device(grad)
      options%nsweep = nsw
   end subroutine xc_grad_response_core

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
