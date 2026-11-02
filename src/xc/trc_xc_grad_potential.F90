!
! d/dR Tr(P V_xc[D]) with both densities held fixed: the exchange-
! correlation term a double hybrid's gradient needs and no gradient before
! it did. Its perturbative part is not stationary in the orbitals, so the
! Z-vector leaves the relaxed correlation density P contracted with the
! derivative of the REFERENCE operator, and a Kohn-Sham reference's
! operator contains V_xc[D]. Omitting it leaves a gradient of plausible
! size, short by a piece the order of the Coulomb one beside it
! (metalquicha's xc_potential_gradient, mqc_czt_gradient.f90, which this
! ports: same four channels, same coefficients).
!
! THE TERM. P is relaxed and symmetric but indefinite and traceless, and
! appears only linearly, so this is the derivative of a linear form, not of
! an energy. Restricted, f(rho, sigma) the functional at D:
!
!     G(R) = sum_g w_g [ v_rho rho_P + 2 v_sigma grad rho . grad rho_P ]
!
! Every object in it moves with the nuclei the way trc_xc_grad's header
! lays out for E_xc -- the basis functions (Pulay), the points with their
! owning atom, the Becke weights -- and v_rho, v_sigma move with D's own
! basis functions too, which is where f_xc enters. Collected by WHICH
! density's basis derivative multiplies what (the channels):
!
!   d rho_D       : f_rr rho_P + 2 f_rs (grad rho . grad rho_P)
!   d grad rho_D  : 2 f_rs rho_P grad rho + 4 f_ss (grad rho . grad rho_P) grad rho
!                   + 2 v_sigma grad rho_P
!   d rho_P       : v_rho
!   d grad rho_P  : 2 v_sigma grad rho
!
! Each pair is exactly what trc_xc_grad's Pulay and point-motion kernels
! contract for an energy, whose coefficients are v_rho on d rho and 2
! v_sigma grad rho on d grad rho -- so each channel is ONE xc_grad_drive
! call in its external-potential mode, with D or P as the density and
! (v_rho, v_sigma, grad rho) set so that those two coefficients are the
! channel's: the first pair with v_sigma = 1/2 and "grad rho" the vector
! coefficient itself, the second with the potential at D. The point-motion
! term separates the same way, being the same per-(point, function)
! products summed over functions instead of points. The partition-weight
! response needs only the integrand per point (`energy` below, times the
! quadrature weight qw), as for E_xc.
!
! The per-point coefficients come from trc_xc's trc_xc_kernel_points_rks
! (first derivatives from xc_eval_point, the energy's own; second from
! xc_eval_point_fxc, the CPKS kernel's), on the batches the reference
! potential was built over: the same density screen as trc_xc_grad, on D.
!
! LDA and GGA, restricted, Cartesian AO (both densities in the Cartesian
! space the grid's collocation runs in). Validated by the identity that
! defines it: E_xc[D + t P] = E_xc[D] + t Tr(P V_xc[D]) + O(t^2) at every
! geometry, so this is the t-derivative of trc_xc_grad's own gradient at
! D + t P (test/check_dh_grad.F90 differences it).
!
module trc_xc_grad_potential
   use pic_types, only: dp
   use trc_basis_types, only: trc_basis_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, XC_FAMILY_LDA, XC_FAMILY_GGA
   use trc_linalg, only: trc_linalg_t
   use trc_xc, only: trc_xc_kernel_points_rks, trc_xc_kernel_points_t
   use trc_xc_grad, only: xc_grad_drive, xc_grad_options_t, xc_grad_point_io_t, xc_grad_step_names, &
                          xc_grad_batch_active
   use trc_dft_partition, only: dft_partition_becke_grad_contract
   use trc_step_table, only: step_table_t
   use trc_error, only: error_t
   implicit none
   private

   public :: xc_grad_potential_add

contains

   subroutine xc_grad_potential_add(b, lp, func, xg, density, pmat, grad, options)
      !! grad += d/dR Tr(P V_xc[D]), D = `density` and P = `pmat` held fixed
      !! (this module's header). Both (b%nao, b%nao) HOST arrays in the
      !! Cartesian AO space; `xg` the reference's grid, device-resident and
      !! built with its atoms and quadrature weights (`build_from_grid`);
      !! `grad` device-resident, as xc_grad_add_on_grid's. `options` as
      !! xc_grad_add_on_grid's: screen, grid response, the grid's partition.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: lp
      type(trc_xc_functional_t), intent(in) :: func
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: density(:, :), pmat(:, :)
      real(dp), intent(inout) :: grad(3, b%natm)
      type(xc_grad_options_t), intent(in) :: options

      type(trc_xc_kernel_points_t) :: pts
      type(xc_grad_point_io_t) :: pio
      type(trc_xc_functional_t) :: func_ext
      type(step_table_t) :: timing
      type(error_t) :: err
      real(dp), allocatable :: dens3(:, :, :), pmat3(:, :, :), grad_extra(:, :), pcoef(:)
      integer, allocatable :: active(:), at_elem(:)
      integer :: n, ib, npts

      if (func%family /= XC_FAMILY_LDA .and. func%family /= XC_FAMILY_GGA) then
         error stop "trc_xc_grad_potential: only LDA and GGA functionals have this term here"
      end if
      if (xg%npts == 0) return
      n = b%nao
      npts = xg%npts
      allocate (dens3(n, n, 1), pmat3(n, n, 1), active(xg%nbatch))
      dens3(:, :, 1) = density
      pmat3(:, :, 1) = pmat
      ! The reference potential's own batches: trc_xc_grad's density screen,
      ! on D (a batch it drops carries no V_xc[D], so nothing here either).
      do ib = 1, xg%nbatch
         active(ib) = 0
         if (xg%b_off(ib + 1) > xg%b_off(ib) .and. xg%b_aooff(ib + 1) > xg%b_aooff(ib)) then
            if (xc_grad_batch_active(ib, xg, n, 1, dens3, options%rho_tol)) active(ib) = 1
         end if
      end do
      call trc_xc_kernel_points_rks(b, xg, func, density, pmat, active, pts, lp)

      ! The external-potential mode reads the functional's family only; an
      ! LDA's GGA coefficients are all zero, so running it as a GGA is exact.
      func_ext = func
      func_ext%family = XC_FAMILY_GGA
      call timing%init("xc potential grad", xc_grad_step_names())
      pio%options = options
      allocate (pio%active(xg%nbatch), source=active)
      allocate (pio%vrho(npts), pio%vsigma(npts), pio%grho(3, npts))
      !$acc enter data copyin(dens3, pmat3) create(pio%vrho, pio%vsigma, pio%grho)

      ! d rho_D, d grad rho_D: the second derivatives' channel, on D.
      call copy_channel(npts, pts%coef_d, pts%vec_d, pio%vrho, pio%grho)
      call fill_value(npts, 0.5_dp, pio%vsigma)
      call xc_grad_drive(b, lp, func_ext, dens3, xg, grad, timing, pio)
      ! d rho_P, d grad rho_P: the potential at D, on P.
      call copy_channel(npts, pts%vrho, pts%grho, pio%vrho, pio%grho)
      call copy_values(npts, pts%vsigma, pio%vsigma)
      call xc_grad_drive(b, lp, func_ext, pmat3, xg, grad, timing, pio)

      if (options%grid_response) then
         ! The partition's response, both pieces (trc_xc_grad's header, B and
         ! C), from the integrand times the quadrature weight alone.
         allocate (at_elem(b%natm), source=b%at_elem)
         allocate (grad_extra(3, b%natm), pcoef(npts))
         grad_extra = 0.0_dp
         pcoef = xg%qw
         !$acc enter data copyin(pcoef)
         call scale_values(npts, pts%energy, pcoef)
         call dft_partition_becke_grad_contract(xg%r, b%at_r, at_elem, xg%atom, pcoef, options%partition, &
                                                grad_extra, err)
         if (err%has_error()) error stop "trc_xc_grad_potential: "//trim(err%get_message())
         !$acc exit data delete(pcoef)
         !$acc update self(grad)
         grad = grad + grad_extra
         !$acc update device(grad)
      end if
      !$acc exit data delete(dens3, pmat3, pio%vrho, pio%vsigma, pio%grho)
      call pts%release()
   end subroutine xc_grad_potential_add

   ! ----------------------------------------------------------------------
   ! Elementwise device copies: the point arrays stay where they were made.
   ! ----------------------------------------------------------------------

   subroutine copy_channel(npts, scalar, vector, scalar_out, vector_out)
      !! scalar_out := scalar, vector_out := vector, (npts) and (3, npts).
      integer, intent(in) :: npts
      real(dp), intent(in) :: scalar(npts), vector(3, npts)
      real(dp), intent(inout) :: scalar_out(npts), vector_out(3, npts)
      integer :: g
      do concurrent(g=1:npts)
         scalar_out(g) = scalar(g)
         vector_out(1, g) = vector(1, g)
         vector_out(2, g) = vector(2, g)
         vector_out(3, g) = vector(3, g)
      end do
   end subroutine copy_channel

   subroutine copy_values(npts, x, y)
      !! y := x.
      integer, intent(in) :: npts
      real(dp), intent(in) :: x(npts)
      real(dp), intent(inout) :: y(npts)
      integer :: g
      do concurrent(g=1:npts)
         y(g) = x(g)
      end do
   end subroutine copy_values

   subroutine fill_value(npts, value, y)
      !! y := value.
      integer, intent(in) :: npts
      real(dp), intent(in) :: value
      real(dp), intent(inout) :: y(npts)
      integer :: g
      do concurrent(g=1:npts)
         y(g) = value
      end do
   end subroutine fill_value

   subroutine scale_values(npts, x, y)
      !! y := x y, elementwise.
      integer, intent(in) :: npts
      real(dp), intent(in) :: x(npts)
      real(dp), intent(inout) :: y(npts)
      integer :: g
      do concurrent(g=1:npts)
         y(g) = x(g)*y(g)
      end do
   end subroutine scale_values

end module trc_xc_grad_potential
