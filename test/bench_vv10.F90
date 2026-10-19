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
! Under `mpiexec -np N`, one rank per device (`trc_bind_device`, the same
! call app/main.F90 makes), `trc_vv10_rks` is handed the world comm and
! pass 2 splits by rank -- see vv10_drive's own header (MULTI-GPU). Every
! rank prints its own E_nlc and wall time, so a scaling run (1 vs N
! ranks) can diff them: E_nlc must agree to the SCF floor regardless of
! rank count, since the allreduce after pass 2 reconstructs the same
! arrays every point's contribution was always going to add up to.
!
! NO internal `wall()` function here, unlike bench_xc/bench_sph: a
! `contains`-ed procedure in this program's own scope hit a known
! nvfortran + HPC-X MPI defect at more than one rank (silently wrong
! MPI_Allreduce results -- a scalar sum of 1.0 and 2.0 across two ranks
! came back 0.0, not 3.0, with nothing else in this file touched), the
! same class of bug PICKUP.md already records for a `contains` procedure
! in an MPI main program ("any `contains` procedure in the main program
! zeroed G at >1 rank"). `system_clock` is inlined at each call site
! instead, exactly as test/scf_mpi.F90 already does.
!
program bench_vv10
   use pic_types, only: dp, int64
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_error, only: error_t
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_vv10, only: trc_vv10_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   implicit none

   character(len=256) :: xyzfile, arg, bas_json
   integer :: level, nat, nao, i, j, dev
   integer, allocatable :: z(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp), allocatable :: c(:, :), dmat(:, :), vmat(:, :)
   type(trc_basis_t) :: bas
   type(dft_grid_t) :: grid
   type(trc_xc_grid_t) :: xg
   type(error_t) :: err
   type(comm_t) :: comm
   real(dp) :: e_nlc, seed, t0, t1
   integer(int64) :: lcg, cc, rate
   integer, parameter :: nocc = 100

   call pic_mpi_init()
   comm = comm_world()
   dev = trc_bind_device(comm%rank())

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

   call system_clock(cc, rate); t0 = real(cc, dp)/real(rate, dp)
   call build_dft_grid(at_r, z, grid, err, level=level)
   call system_clock(cc, rate); t1 = real(cc, dp)/real(rate, dp)
   if (err%has_error()) then
      print '(a)', 'bench_vv10: grid failed: '//err%get_message(); stop 1
   end if
   print '(a,f8.3,a)', '  grid build ', t1 - t0, ' s'
   call system_clock(cc, rate); t0 = real(cc, dp)/real(rate, dp)
   call xg%build(grid%n_points, grid%coords, grid%weights, bas, max_pts=512)
   call xg%to_device()
   call system_clock(cc, rate); t1 = real(cc, dp)/real(rate, dp)
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
   call system_clock(cc, rate); t0 = real(cc, dp)/real(rate, dp)
   call trc_vv10_rks(bas, xg, dmat, 6.0_dp, 0.01_dp, e_nlc, vmat, comm=comm)
   call system_clock(cc, rate); t1 = real(cc, dp)/real(rate, dp)
   !$acc exit data delete(dmat, vmat)

   print '(a,i0,a,i0,a,i0,a,i0,a,f20.14,a,f10.3,a)', '  rank ', comm%rank(), ' of ', comm%size(), &
      ' (device ', dev, ') VV10 at ', grid%n_points, ' points: E_nlc = ', e_nlc, &
      '  time = ', t1 - t0, ' s'

   call xg%release(); call grid%destroy(); call bas%release()
   call pic_mpi_finalize()

end program bench_vv10
