!
! Where VV10's time goes, on a real molecule in a real basis.
!
!   bench_vv10 <xyz> [level] [basis.json]
!
! Builds 6-31G on the molecule by default, a VV10 grid at the given level,
! a fixed positive density, and times `trc_vv10_rks` once (its own header
! explains why this is O(N_points^2) and not meant to be timed over many
! reps the way bench_xc times the semilocal pipeline). A timing is not a
! pass, so this is built and not registered as a test.
!
program bench_vv10
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_vv10, only: trc_vv10_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   implicit none

   character(len=256) :: xyzfile, arg, bas_json
   integer :: level, nat, nao, i, j
   integer, allocatable :: z(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp), allocatable :: c(:, :), dmat(:, :), vmat(:, :)
   type(trc_basis_t) :: bas
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(error_t) :: err
   real(dp) :: e_nlc, seed, t0, t1
   integer(int64) :: lcg
   integer, parameter :: nocc = 100

   level = 1; bas_json = ''
   if (command_argument_count() < 1) then
      print '(a)', 'usage: bench_vv10 <xyz> [level] [basis.json]'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read (arg, *) level
   end if
   if (command_argument_count() >= 3) call get_command_argument(3, bas_json)

   call read_xyz(trim(xyzfile), nat, z, at_r)
   if (len_trim(bas_json) > 0) then
      call trc_basis_from_json(trim(bas_json), nat, z, at_r, bas, err, force_cart=.true.)
   else
      call trc_basis_from_json('basis_sets/6-31g.json', nat, z, at_r, bas, err, force_cart=.true.)
   end if
   if (err%has_error()) then
      print '(a)', 'bench_vv10: '//err%get_message(); stop 1
   end if
   call bas%to_device()
   nao = bas%nao

   t0 = wall()
   call build_dft_grid(at_r, z, grid, err, level=level)
   t1 = wall()
   if (err%has_error()) then
      print '(a)', 'bench_vv10: grid failed: '//err%get_message(); stop 1
   end if
   print '(a,f8.3,a)', '  grid build ', t1 - t0, ' s'
   t0 = wall()
   call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=512)
   call xg%to_device()
   t1 = wall()
   print '(a,i0,a,i0,a,i0,a,i0)', 'bench_vv10: ', nat, ' atoms, ', nao, ' AOs, ', grid%n_points, &
      ' points at level ', level
   print '(a,i0,a,f8.3,a)', '  batching ', xg%nbatch, ' s ', t1 - t0, ' s'

   allocate (c(nao, min(nocc, nao)), dmat(nao, nao), vmat(nao, nao))
   lcg = 7_int64
   do j = 1, min(nocc, nao)
      do i = 1, nao
         lcg = mod(48271_int64*lcg, 2147483647_int64)
         seed = real(lcg, dp)/2147483647.0_dp
         c(i, j) = 0.1_dp*(seed - 0.5_dp)
      end do
   end do
   dmat = 2.0_dp*matmul(c, transpose(c))
   vmat = 0.0_dp

   !$acc enter data copyin(dmat, vmat)
   t0 = wall()
   call trc_vv10_rks(bas, xg, dmat, 6.0_dp, 0.01_dp, e_nlc, vmat)
   t1 = wall()
   !$acc exit data delete(dmat, vmat)

   print '(a,i0,a,f20.14,a,f10.3,a)', '  VV10 at ', grid%n_points, ' points: E_nlc = ', e_nlc, &
      '  time = ', t1 - t0, ' s'

   call xg%release(); call grid%destroy(); call bas%release()

contains

   function wall() result(t)
      real(dp) :: t
      integer(int64) :: cc, rate
      call system_clock(cc, rate)
      t = real(cc, dp)/real(rate, dp)
   end function wall

end program bench_vv10
