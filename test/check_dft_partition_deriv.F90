!
! dft_partition_becke_derivatives (trc_dft_partition.F90): the O(n_atoms^2)-
! per-point production routine against a hand-kept COPY of the
! O(n_atoms^3)-per-point implementation it replaced -- the one
! trc_xc_grad.F90's grid-response term was validated against (FD at two
! grid levels, translational invariance, and pyscf's own grid_response=True
! to ~1e-13/1e-11 -- see that module's header) before the rewrite. Water
! (3 atoms) and formamide (6 atoms), both PARTITION_BECKE/ADJUST_TREUTLER
! (terco's own defaults) and PARTITION_STRATMANN/ADJUST_NONE (the other
! combination trc_dft_partition.F90 offers), level 3, every grid point --
! not a handful of synthetic ones, since the point of this check is that
! the rewrite did not change the answer anywhere on a real quadrature.
!
! The reference below is NOT independently re-derived: it is the routine
! this file's own git history had until the O(n_atoms^2) rewrite, kept
! here verbatim (down to variable names) so this test is a diff against
! what used to ship, not a second attempt at the same algebra.
!
module check_dft_partition_deriv_ref
   use pic_types, only: dp
   use trc_dft_partition, only: PARTITION_BECKE
   implicit none
   private
   public :: reference_derivatives

contains

   ! Re-exposes what trc_dft_partition.F90 keeps private, duplicated here
   ! rather than exported from the module itself: production code has no
   ! business calling the O(n_atoms^3) form it was rewritten away from.
   subroutine reference_derivatives(points, atom_coords, shift, inv_distance, owner, scheme, dweights)
      real(dp), intent(in) :: points(:, :), atom_coords(:, :)
      real(dp), intent(in) :: shift(:, :), inv_distance(:, :)
      integer, intent(in) :: owner(:), scheme
      real(dp), intent(out) :: dweights(:, :, :)

      real(dp), allocatable :: atom_r(:), unit_vec(:, :)
      real(dp), allocatable :: cell(:), dcell(:, :, :)
      real(dp) :: dtotal(3), rij(3), ds_i(3), ds_j(3), dmu_i(3), dmu_j(3)
      real(dp) :: mu, nu, s, total, dnu_dmu, ds_dnu, invd
      integer :: n_atoms, n_points, i, j, k, a

      n_atoms = size(atom_coords, 2)
      n_points = size(points, 2)
      dweights = 0.0_dp
      if (n_atoms == 1) return

      allocate (atom_r(n_atoms), unit_vec(3, n_atoms))
      allocate (cell(n_atoms), dcell(3, n_atoms, n_atoms))

      do k = 1, n_points
         do i = 1, n_atoms
            atom_r(i) = norm2(points(:, k) - atom_coords(:, i))
            if (atom_r(i) > 0.0_dp) then
               unit_vec(:, i) = (points(:, k) - atom_coords(:, i))/atom_r(i)
            else
               unit_vec(:, i) = 0.0_dp
            end if
         end do

         cell = 1.0_dp
         dcell = 0.0_dp

         do i = 1, n_atoms
            do j = i + 1, n_atoms
               invd = inv_distance(i, j)
               rij = atom_coords(:, i) - atom_coords(:, j)

               mu = (atom_r(i) - atom_r(j))*invd
               nu = mu + shift(i, j)*(1.0_dp - mu*mu)

               if (scheme == PARTITION_BECKE) then
                  s = becke_cutoff_ref(nu)
                  ds_dnu = becke_cutoff_derivative_ref(nu)
               else
                  s = stratmann_cutoff_ref(nu)
                  ds_dnu = stratmann_cutoff_derivative_ref(nu)
               end if

               dnu_dmu = 1.0_dp - 2.0_dp*shift(i, j)*mu

               dmu_i = -unit_vec(:, i)*invd - (atom_r(i) - atom_r(j))*invd**3*rij
               dmu_j = unit_vec(:, j)*invd + (atom_r(i) - atom_r(j))*invd**3*rij

               ds_i = ds_dnu*dnu_dmu*dmu_i
               ds_j = ds_dnu*dnu_dmu*dmu_j

               do a = 1, n_atoms
                  dcell(:, a, i) = dcell(:, a, i)*s
                  dcell(:, a, j) = dcell(:, a, j)*(1.0_dp - s)
               end do
               dcell(:, i, i) = dcell(:, i, i) + cell(i)*ds_i
               dcell(:, j, i) = dcell(:, j, i) + cell(i)*ds_j
               dcell(:, i, j) = dcell(:, i, j) - cell(j)*ds_i
               dcell(:, j, j) = dcell(:, j, j) - cell(j)*ds_j

               cell(i) = cell(i)*s
               cell(j) = cell(j)*(1.0_dp - s)
            end do
         end do

         total = sum(cell)
         if (total <= 0.0_dp) cycle

         do a = 1, n_atoms
            dtotal = 0.0_dp
            do i = 1, n_atoms
               dtotal = dtotal + dcell(:, a, i)
            end do
            dweights(:, a, k) = (dcell(:, a, owner(k))*total &
                                 - cell(owner(k))*dtotal)/(total*total)
         end do
      end do
   end subroutine reference_derivatives

   pure function becke_cutoff_ref(nu) result(s)
      real(dp), intent(in) :: nu
      real(dp) :: s, f
      f = nu
      f = 0.5_dp*f*(3.0_dp - f*f)
      f = 0.5_dp*f*(3.0_dp - f*f)
      f = 0.5_dp*f*(3.0_dp - f*f)
      s = 0.5_dp*(1.0_dp - f)
   end function becke_cutoff_ref

   pure function becke_cutoff_derivative_ref(nu) result(ds)
      real(dp), intent(in) :: nu
      real(dp) :: ds, f0, f1, f2
      f0 = nu
      f1 = 0.5_dp*f0*(3.0_dp - f0*f0)
      f2 = 0.5_dp*f1*(3.0_dp - f1*f1)
      ds = -0.5_dp*(1.5_dp*(1.0_dp - f2*f2))*(1.5_dp*(1.0_dp - f1*f1))*(1.5_dp*(1.0_dp - f0*f0))
   end function becke_cutoff_derivative_ref

   pure function stratmann_cutoff_ref(mu) result(s)
      real(dp), intent(in) :: mu
      real(dp) :: s, z, z2
      real(dp), parameter :: A = 0.64_dp
      if (mu <= -A) then
         s = 1.0_dp
      else if (mu >= A) then
         s = 0.0_dp
      else
         z = mu/A
         z2 = z*z
         z = (1.0_dp/16.0_dp)*(z*(35.0_dp + z2*(-35.0_dp + z2*(21.0_dp - 5.0_dp*z2))))
         s = 0.5_dp*(1.0_dp - z)
      end if
   end function stratmann_cutoff_ref

   pure function stratmann_cutoff_derivative_ref(mu) result(ds)
      real(dp), intent(in) :: mu
      real(dp) :: ds, z, z2
      real(dp), parameter :: A = 0.64_dp
      if (mu <= -A .or. mu >= A) then
         ds = 0.0_dp
      else
         z = mu/A
         z2 = z*z
         ds = -0.5_dp*(1.0_dp/16.0_dp)*(35.0_dp + z2*(-105.0_dp + z2*(105.0_dp - 35.0_dp*z2)))/A
      end if
   end function stratmann_cutoff_derivative_ref

end module check_dft_partition_deriv_ref

program check_dft_partition_deriv
   use pic_types, only: dp
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_dft_grid, only: dft_grid_t, dft_grid_options_t, dft_grid_build
   use pic_types, only: int64
   use trc_dft_partition, only: dft_partition_becke_derivatives, dft_partition_becke_grad_contract, &
                                dft_partition_options_t, PARTITION_BECKE, PARTITION_STRATMANN, &
                                ADJUST_TREUTLER, ADJUST_NONE
   use check_dft_partition_deriv_ref, only: reference_derivatives
   implicit none

   logical :: any_fail
   any_fail = .false.

   call one_case('water.xyz', PARTITION_BECKE, ADJUST_TREUTLER, 'water, Becke/Treutler', any_fail)
   call one_case('water.xyz', PARTITION_STRATMANN, ADJUST_NONE, 'water, Stratmann/none', any_fail)
   call one_case('formamide.xyz', PARTITION_BECKE, ADJUST_TREUTLER, 'formamide, Becke/Treutler', any_fail)
   call one_case('formamide.xyz', PARTITION_STRATMANN, ADJUST_NONE, 'formamide, Stratmann/none', any_fail)
   ! A water cluster, where most of a point's active atoms have negligible
   ! cells: the case dft_partition_becke_grad_contract's cell screening is for.
   call large_case('w5.xyz', 'w5, Becke/Treutler, level 1', any_fail)

   if (any_fail) then
      print '(a)', 'check_dft_partition_deriv: RESULT: FAIL'; stop 1
   end if
   print '(a)', 'check_dft_partition_deriv: RESULT: PASS'

contains

   subroutine one_case(xyzfile, scheme, adjust, label, any_fail)
      character(*), intent(in) :: xyzfile, label
      integer, intent(in) :: scheme, adjust
      logical, intent(inout) :: any_fail

      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      real(dp), allocatable :: shift(:, :), inv_distance(:, :)
      real(dp), allocatable :: dw_new(:, :, :), dw_ref(:, :, :)
      real(dp) :: worst
      integer :: i, j
      type(dft_grid_options_t) :: gopts
      type(dft_partition_options_t) :: unscreened

      call read_xyz(xyzfile, nat, at_z, r)
      gopts%partition%scheme = scheme
      gopts%partition%adjust = adjust
      call dft_grid_build(r, at_z, grid, err, level=3, options=gopts)
      if (err%has_error()) then
         print '(a)', 'check_dft_partition_deriv: '//err%get_message(); stop 1
      end if

      allocate (dw_new(3, nat, grid%n_points), dw_ref(3, nat, grid%n_points))
      ! nu_max=2 disables the production routine's own active-set
      ! screening (trc_dft_partition.F90's own header: screened and
      ! unscreened Becke/Treutler are NOT expected to agree, by design --
      ! this check is isolating the O(n_atoms^2) rewrite's algebra from
      ! that approximation, not re-measuring it). `reference_derivatives`
      ! below has no screening of its own at all, so this is the
      ! apples-to-apples comparison.
      unscreened = grid%partition
      unscreened%nu_max = 2.0_dp
      call dft_partition_becke_derivatives(grid%coords, r, at_z, grid%atom, unscreened, dw_new, err)
      if (err%has_error()) then
         print '(a)', 'check_dft_partition_deriv: '//err%get_message(); stop 1
      end if

      ! The reference needs `shift`/`inv_distance` exactly as the
      ! production routine builds them -- re-derived here from the same
      ! public `size_adjustment`-equivalent inputs is not available
      ! (private), so this reconstructs them from the same formula
      ! trc_dft_partition.F90's own `size_adjustment` uses, which is
      ! public knowledge from that module's header (Becke's appendix,
      ! Treutler's sqrt variant) and not duplicated blindly: a mismatch
      ! here would show up as a wrong answer below, not a silent pass.
      allocate (shift(nat, nat), inv_distance(nat, nat))
      call shift_and_distance(nat, r, at_z, adjust, shift, inv_distance)
      call reference_derivatives(grid%coords, r, shift, inv_distance, grid%atom, scheme, dw_ref)

      worst = maxval(abs(dw_new - dw_ref))
      print '(a,a,a,es10.2,a,i0)', '  ', label, ': worst|new-old| ', worst, '  over ', grid%n_points
      if (worst > 1.0e-13_dp) any_fail = .true.

      ! dft_partition_becke_grad_contract against the dense derivatives
      ! contracted the way the XC gradient used to: with the production
      ! screening and with none, a per-point factor that is zero on every
      ! seventh point (the skip path), and a budget small enough to force
      ! many chunks.
      call contract_case(grid, nat, r, at_z, scheme, adjust, label//', screened', -1.0_dp, any_fail)
      call contract_case(grid, nat, r, at_z, scheme, adjust, label//', unscreened', 2.0_dp, any_fail)

      call grid%destroy()
      deallocate (at_z, r, dw_new, dw_ref, shift, inv_distance)
   end subroutine one_case

   subroutine large_case(xyzfile, label, any_fail)
      character(*), intent(in) :: xyzfile, label
      logical, intent(inout) :: any_fail
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      type(dft_grid_t) :: grid
      type(error_t) :: err
      type(dft_grid_options_t) :: gopts

      call read_xyz(xyzfile, nat, at_z, r)
      gopts%partition%scheme = PARTITION_BECKE
      gopts%partition%adjust = ADJUST_TREUTLER
      call dft_grid_build(r, at_z, grid, err, level=1, options=gopts)
      if (err%has_error()) then
         print '(a)', 'check_dft_partition_deriv: '//err%get_message(); stop 1
      end if
      call contract_case(grid, nat, r, at_z, PARTITION_BECKE, ADJUST_TREUTLER, label, -1.0_dp, any_fail)
      call grid%destroy()
   end subroutine large_case

   subroutine contract_case(grid, nat, r, at_z, scheme, adjust, label, numax, any_fail)
      type(dft_grid_t), intent(in) :: grid
      integer, intent(in) :: nat, scheme, adjust
      real(dp), intent(in) :: r(3, nat), numax
      integer, intent(in) :: at_z(nat)
      character(*), intent(in) :: label
      logical, intent(inout) :: any_fail
      real(dp), allocatable :: dw(:, :, :), coef(:), g_ref(:, :), g_new(:, :)
      type(error_t) :: err
      real(dp) :: worst, scale
      integer :: k, a, own
      type(dft_partition_options_t) :: popts

      allocate (dw(3, nat, grid%n_points), coef(grid%n_points), g_ref(3, nat), g_new(3, nat))
      do k = 1, grid%n_points
         coef(k) = grid%quad_weights(k)*(1.0_dp + 0.3_dp*sin(real(k, dp)))
         if (mod(k, 7) == 0) coef(k) = 0.0_dp
      end do
      popts%scheme = scheme
      popts%adjust = adjust
      if (numax > 0.0_dp) popts%nu_max = numax
      call dft_partition_becke_derivatives(grid%coords, r, at_z, grid%atom, popts, dw, err)
      if (err%has_error()) then
         print '(a)', 'check_dft_partition_deriv: '//err%get_message(); stop 1
      end if
      g_ref = 0.0_dp
      do k = 1, grid%n_points
         own = grid%atom(k)
         do a = 1, nat
            g_ref(:, a) = g_ref(:, a) + coef(k)*dw(:, a, k)
            g_ref(:, own) = g_ref(:, own) - coef(k)*dw(:, a, k)
         end do
      end do
      g_new = 1.0_dp   ! ADDED into: start from a non-zero gradient
      popts%budget = int(4*nat*37, int64)
      call dft_partition_becke_grad_contract(grid%coords, r, at_z, grid%atom, coef, popts, g_new, err)
      if (err%has_error()) then
         print '(a)', 'check_dft_partition_deriv: '//err%get_message(); stop 1
      end if
      g_new = g_new - 1.0_dp
      scale = max(1.0_dp, maxval(abs(g_ref)))
      worst = maxval(abs(g_new - g_ref))
      print '(a,a,a,es10.2,a,es10.2)', '  ', label, ': contracted worst|new-dense| ', worst, &
         '  of max ', maxval(abs(g_ref))
      if (worst > 1.0e-13_dp*scale .or. maxval(abs(g_ref)) == 0.0_dp) any_fail = .true.
      deallocate (dw, coef, g_ref, g_new)
   end subroutine contract_case

   ! Bragg radii and the two adjustment rules, duplicated from
   ! trc_dft_partition.F90%size_adjustment (private there) -- see
   ! one_case's own comment on why a mismatch cannot hide.
   subroutine shift_and_distance(n_atoms, atom_coords, atomic_numbers, adjust, shift, inv_distance)
      use trc_dft_radial, only: dft_radial_bragg_radius
      integer, intent(in) :: n_atoms, atomic_numbers(n_atoms), adjust
      real(dp), intent(in) :: atom_coords(3, n_atoms)
      real(dp), intent(out) :: shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms)
      real(dp), parameter :: TINY_R = 1.0e-200_dp, MAX_ADJ = 0.5_dp
      real(dp) :: radius(n_atoms), chi, a
      integer :: i, j

      select case (adjust)
      case (ADJUST_TREUTLER)
         do i = 1, n_atoms
            radius(i) = sqrt(dft_radial_bragg_radius(atomic_numbers(i))) + TINY_R
         end do
      case default
         radius = 1.0_dp
      end select

      shift = 0.0_dp
      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i == j) cycle
            inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
            if (adjust == ADJUST_TREUTLER) then
               chi = radius(i)/radius(j)
               a = 0.25_dp*(1.0_dp/chi - chi)
               shift(i, j) = max(-MAX_ADJ, min(MAX_ADJ, a))
            end if
         end do
      end do
   end subroutine shift_and_distance

end program check_dft_partition_deriv
