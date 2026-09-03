!
! Where the XC time goes, on a real molecule in a real basis.
!
!   bench_xc <xyz> [level] [reps] [max_pts]
!
! Builds 6-31G on the molecule, a grid at the given level with the default
! pruning, a fixed positive density, and times `trc_xc_rks` for PBE:
! the batching, then `reps` evaluations split into the point kernel
! (collocation, density, functional) and the pair kernel (V_xc). A timing
! is not a pass, so this is built and not registered as a test.
!
program bench_xc
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: trc_xc_rks
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   character(len=256) :: xyzfile, arg
   integer :: level, reps, max_pts, nat, nsh, maxnp, nao, i, j, r
   integer, allocatable :: z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), zat(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: c(:, :), dmat(:, :), vxc(:, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
   real(dp) :: trds
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
   type(error_t) :: err
   real(dp) :: exc, nelec, seed, t0, t1, tp, tq, sp, sq, tmin
   integer :: nsk
   integer(kind=8) :: lcg
   logical :: ok

   level = 3; reps = 5; max_pts = 512
   if (command_argument_count() < 1) then
      print '(a)', 'usage: bench_xc <xyz> [level] [reps] [max_pts]'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read (arg, *) level
   end if
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read (arg, *) reps
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg); read (arg, *) max_pts
   end if

   call read_xyz(trim(xyzfile), nat, z, at_r)
   allocate (zat(nat)); zat = real(z, dp)
   call build_631g(nat, z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call bas%to_device()
   nao = bas%nao

   call build_dft_grid(at_r, z, grid, err, level=level)
   if (err%has_error()) then
      print '(a)', 'bench_xc: grid failed: '//err%get_message(); stop 1
   end if
   t0 = wall()
   call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=max_pts)
   call xg%to_device()
   t1 = wall()
   print '(a,i0,a,i0,a,i0,a,i0)', 'bench_xc: ', nat, ' atoms, ', nao, ' AOs, ', grid%n_points, &
      ' points at level ', level
   print '(a,i0,a,i0,a,f6.1,a,f8.3,a)', '  batches ', xg%nbatch, ', largest local basis ', xg%max_nloc, &
      ', mean ', real(xg%b_aooff(xg%nbatch + 1) - 1, dp)/real(xg%nbatch, dp), &
      '   (batching ', t1 - t0, ' s)'

   allocate (c(nao, 8), dmat(nao, nao), vxc(nao, nao))
   lcg = 61_8
   do j = 1, 8
      do i = 1, nao
         lcg = mod(25214903917_8*lcg + 11_8, 281474976710656_8)
         seed = real(lcg, dp)/281474976710656.0_dp
         c(i, j) = 0.3_dp*(seed - 0.5_dp)
      end do
   end do
   dmat = 2.0_dp*matmul(c, transpose(c))
   ok = xc_functional_by_name('pbe', func)
   ! The control: sum_g w rho must be Tr(DS) to the grid's accuracy.
   call pl%build(bas, 1.0e-30_dp); call pl%to_device()
   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
   call trc_1e(bas, pl, smat, tmat, vmat)
   trds = sum(dmat*smat)

   sp = 0.0_dp; sq = 0.0_dp; tmin = huge(1.0_dp)
   do r = 1, reps + 1
      t0 = wall()
      call trc_xc_rks(bas, xg, func, dmat, vxc, exc, nelec, t_points=tp, t_pairs=tq, n_skipped=nsk)
      t1 = wall()
      if (r == 1) cycle   ! warm-up: first launches, allocations
      sp = sp + tp; sq = sq + tq; tmin = min(tmin, t1 - t0)
   end do
   print '(a,f9.4,a,f9.4,a,f9.4,a,i0,a,i0,a)', '  per call: points ', sp/reps, ' s   pairs ', sq/reps, &
      ' s   total (best) ', tmin, ' s   (', nsk, ' of ', xg%nbatch, ' batches screened out)'
   print '(a,f12.8,a,f14.10,a,f14.10)', '  E_xc(PBE) ', exc, '   N ', nelec, '   Tr(DS) ', trds

   call xg%release(); call grid%destroy(); call bas%release()

contains

   function wall() result(t)
      real(dp) :: t
      integer(kind=8) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

end program bench_xc
