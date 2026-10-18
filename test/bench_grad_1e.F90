!
! Timing only: the one-electron gradient (trc_1e_gradient) next to
! eri%gradient's own two-electron part, on the same system and basis, GPU 0.
! Not a correctness check -- see check_1e_grad_vs_pyscf.F90 and
! check_full_grad_fd.F90 for that -- just where the time goes. The density
! is a fixed synthetic one (same closed form as the gradient correctness
! checks): a timing run has no SCF to converge, only shell pairs and
! primitives to walk, and that is independent of which density is plugged
! in.
!
!   bench_grad_1e <xyz> <basis.json> [thresh]
!
program bench_grad_1e
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e_gradient, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none

   character(len=256) :: xyzfile, basis_json, arg
   real(dp) :: thresh
   integer :: nat, nao, i, j, dev
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: at_r(:, :)
   type(error_t) :: err
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t) :: eri
   real(dp), allocatable :: wmat(:, :), dmat(:, :), grad(:, :)
   integer(int64) :: c0, c1, rate
   real(dp) :: t_1e, t_2e

   thresh = 1.0e-10_dp
   if (command_argument_count() < 2) then
      print '(a)', 'usage: bench_grad_1e <xyz> <basis.json> [thresh]'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read (arg, *) thresh
   end if

   dev = trc_bind_device(0)
   call read_xyz(trim(xyzfile), nat, at_z, at_r)
   call trc_basis_from_json(trim(basis_json), nat, at_z, at_r, bas, err)
   if (err%has_error()) then
      print '(a)', 'bench_grad_1e: '//err%get_message(); stop 1
   end if
   nao = bas%nao
   print '(a,i0,a,i0,a,i0)', 'bench_grad_1e: natm=', nat, '  nao=', nao, '  dev=', dev

   call bas%to_device()
   call pl%build(bas, thresh)
   call pl%to_device()
   call eri%build(bas, thresh)

   allocate (wmat(nao, nao), dmat(nao, nao), grad(3, nat))
   do j = 1, nao
      do i = 1, nao
         wmat(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j, dp))
         dmat(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j + 7, dp))
      end do
   end do
   wmat = 0.5_dp*(wmat + transpose(wmat))
   dmat = 0.5_dp*(dmat + transpose(dmat))

   ! warm-up (device init, first-touch, JIT-free by now but still: discard)
   grad = 0.0_dp
   call eri%gradient(bas, reshape(dmat, [nao, nao, 1]), grad)
   grad = 0.0_dp
   !$acc enter data copyin(wmat, dmat, grad)
   call trc_1e_gradient(bas, pl, wmat, dmat, grad)
   !$acc exit data delete(wmat, dmat, grad)

   call system_clock(c0, rate)
   grad = 0.0_dp
   call eri%gradient(bas, reshape(dmat, [nao, nao, 1]), grad)
   call system_clock(c1)
   t_2e = real(c1 - c0, dp)/real(rate, dp)

   call system_clock(c0, rate)
   grad = 0.0_dp
   !$acc enter data copyin(wmat, dmat, grad)
   call trc_1e_gradient(bas, pl, wmat, dmat, grad)
   !$acc update self(grad)
   !$acc exit data delete(wmat, dmat, grad)
   call system_clock(c1)
   t_1e = real(c1 - c0, dp)/real(rate, dp)

   print '(a,f10.4,a)', 'eri%gradient (2e):     ', t_2e, ' s'
   print '(a,f10.4,a)', 'trc_1e_gradient (1e):  ', t_1e, ' s'
   print '(a,f10.4)',   '1e / 2e ratio:         ', t_1e/max(t_2e, 1.0e-12_dp)

   call eri%release()
   call pl%release()
   call bas%release()
end program bench_grad_1e
