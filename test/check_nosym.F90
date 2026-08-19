!
! The folded and enumerated digestions, on a SYMMETRIC density.
!
! They must agree for ANY symmetric density -- that is the whole premise of
! having two of them. `fock` folds the eight-fold permutational symmetry into
! six atomic updates and assumes D(mu,nu) = D(nu,mu); `fock_nosym` enumerates
! the permutations and assumes nothing. Where the assumption holds the two are
! the same operator, so a disagreement here is a bug in one of them and a
! disagreement on an ASYMMETRIC density means nothing until this passes.
!
! Phase 3 recorded this agreeing to 1.24e-14 and then nothing re-ran it. This
! is that control, as a standing test.
!
program check_nosym
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_fock, only: trc_eri_t
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   character(len=256) :: xyzfile
   integer :: nat, nsh, maxnp, nao, i, j
   integer,  allocatable :: at_z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), zat(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: d(:, :), gs(:, :), gn(:, :)
   type(trc_basis_t)    :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t)      :: eri
   real(dp) :: worst, thresh

   thresh = 1.0e-10_dp
   if (command_argument_count() < 1) then
      print '(a)', 'usage: check_nosym <xyz>'; stop 1
   end if
   call get_command_argument(1, xyzfile)

   call read_xyz(xyzfile, nat, at_z, at_r)
   call build_631g(nat, at_z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   allocate (zat(nat)); zat = real(at_z, dp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call pl%build(bas, thresh)
   nao = bas%nao
   call bas%to_device(); call pl%to_device()
   call eri%build(bas, thresh)

   ! Any symmetric density will do; a deterministic one that is not sparse
   ! exercises more of the digestion than a converged one would.
   allocate (d(nao, nao), gs(nao, nao), gn(nao, nao))
   do j = 1, nao
      do i = 1, nao
         d(i, j) = 0.3_dp*cos(real(2*i + 3*j, dp)) + 0.1_dp*real(min(i, j), dp)
      end do
   end do
   d = 0.5_dp*(d + transpose(d))

   call eri%fock(bas, d, gs)
   call eri%fock_nosym(bas, d, gn)
   worst = maxval(abs(gs - gn))

   print '(a,i0,a,i0)', '  ', nao, ' functions, shells ', nsh
   print '(a,es12.4)',  '  asymmetry of the density ', maxval(abs(d - transpose(d)))
   print '(a,es24.16)', '  sum folded     ', sum(gs)
   print '(a,es24.16)', '  sum enumerated ', sum(gn)
   print '(a,f18.10)',  '  ratio          ', sum(gn)/sum(gs)
   print '(a,es12.4)',  '  worst |folded - enumerated| ', worst
   if (worst < 1.0e-10_dp) then
      print '(a)', '  RESULT: PASS'
   else
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if

   call eri%release(); call pl%release(); call bas%release()
end program check_nosym
