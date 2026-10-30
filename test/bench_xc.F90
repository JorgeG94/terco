!
! Where the XC time goes, on a real molecule in a real basis.
!
!   bench_xc <xyz> [level] [reps] [max_pts] [basis.json] [functional]
!
! Builds 6-31G on the molecule by default, or the MolSSI BSE basis named by
! the optional 5th argument (e.g. basis_sets/cc-pvdz.json on sih4.xyz, to put
! a d shell through `collocation_shell` -- 6-31G alone has none, s/p only) --
! a grid at the given level with the default pruning, a fixed positive
! density, and times `trc_xc_rks` for PBE (or the functional named by the
! optional 6th argument, r2scan say, for a meta-GGA): the batching, then `reps`
! evaluations split into the point kernel (collocation, density, functional)
! and the pair kernel (V_xc). A timing is not a pass, so this is built and
! not registered as a test.
!
program bench_xc
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: trc_xc_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   implicit none

   character(len=256) :: xyzfile, arg, bas_json, fname
   integer :: level, reps, max_pts, nat, nao, i, j, r
   integer, allocatable :: z(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp), allocatable :: c(:, :), dmat(:, :), vxc(:, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
   real(dp) :: trds
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(trc_xc_functional_t) :: func
   type(error_t) :: err
   real(dp) :: exc, nelec, seed, t1, tp, tq, sp, sq, tmin
   type(timer_type) :: tm
   integer :: nsk
   integer(int64) :: lcg
   logical :: ok

   level = 3; reps = 5; max_pts = 512; bas_json = ''; fname = 'pbe'
   if (command_argument_count() < 1) then
      print '(a)', 'usage: bench_xc <xyz> [level] [reps] [max_pts] [basis.json] [functional]'
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
   if (command_argument_count() >= 5) call get_command_argument(5, bas_json)
   if (command_argument_count() >= 6) call get_command_argument(6, fname)

   call read_xyz(trim(xyzfile), nat, z, at_r)
   if (len_trim(bas_json) > 0) then
      call trc_basis_from_json(trim(bas_json), nat, z, at_r, bas, err)
      if (err%has_error()) then
         print '(a)', 'bench_xc: '//err%get_message(); stop 1
      end if
   else
      call trc_basis_from_json('basis_sets/6-31g.json', nat, z, at_r, bas, err)
      if (err%has_error()) then
         print '(a)', 'bench_xc: '//err%get_message(); stop 1
      end if
   end if
   call bas%to_device()
   nao = bas%nao

   call tm%start()
   call dft_grid_build(at_r, z, grid, err, level=level)
   t1 = tm%get_elapsed_time()
   if (err%has_error()) then
      print '(a)', 'bench_xc: grid failed: '//err%get_message(); stop 1
   end if
   print '(a,f8.3,a)', '  grid build (Lebedev, radial, partition) ', t1, ' s'
   call tm%start()
   call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=max_pts)
   call xg%to_device()
   t1 = tm%get_elapsed_time()
   print '(a,i0,a,i0,a,i0,a,i0)', 'bench_xc: ', nat, ' atoms, ', nao, ' AOs, ', grid%n_points, &
      ' points at level ', level
   print '(a,i0,a,i0,a,f6.1,a,f8.3,a)', '  batches ', xg%nbatch, ', largest local basis ', xg%max_nloc, &
      ', mean ', real(xg%b_aooff(xg%nbatch + 1) - 1, dp)/real(xg%nbatch, dp), &
      '   (batching ', t1, ' s)'

   allocate (c(nao, 8), dmat(nao, nao), vxc(nao, nao))
   lcg = 61_int64
   do j = 1, 8
      do i = 1, nao
         lcg = mod(48271_int64*lcg, 2147483647_int64)
         seed = real(lcg, dp)/2147483647.0_dp
         c(i, j) = 0.3_dp*(seed - 0.5_dp)
      end do
   end do
   dmat = 2.0_dp*matmul(c, transpose(c))
   ok = xc_functional_by_name(trim(fname), func)
   if (.not. ok .or. len_trim(func%missing) > 0) then
      print '(a)', 'bench_xc: cannot run functional '//trim(fname)
      stop 1
   end if
   ! The control: sum_g w rho must be Tr(DS) to the grid's accuracy.
   call pl%build(bas, 1.0e-30_dp); call pl%to_device()
   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
   call trc_1e(bas, pl, smat, tmat, vmat)
   trds = sum(dmat*smat)

   sp = 0.0_dp; sq = 0.0_dp; tmin = huge(1.0_dp)
   do r = 1, reps + 1
      call tm%start()
      call trc_xc_rks(bas, xg, func, dmat, vxc, exc, nelec, t_points=tp, t_pairs=tq, n_skipped=nsk)
      t1 = tm%get_elapsed_time()
      if (r == 1) cycle   ! warm-up: first launches, allocations
      sp = sp + tp; sq = sq + tq; tmin = min(tmin, t1)
   end do
   print '(a,f9.4,a,f9.4,a,f9.4,a,i0,a,i0,a)', '  per call: points ', sp/reps, ' s   pairs ', sq/reps, &
      ' s   total (best) ', tmin, ' s   (', nsk, ' of ', xg%nbatch, ' batches screened out)'
   print '(a,f16.8,a,f14.8,a,f14.8)', '  E_xc('//trim(fname)//') ', exc, '   N ', nelec, '   Tr(DS) ', trds

   call xg%release(); call grid%destroy(); call bas%release()

contains


end program bench_xc
