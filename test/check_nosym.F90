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
! IT IS ALSO THE ORACLE FOR THE HYBRID DISPATCH, which is why it takes a basis
! ceiling. `fock` routes each class to an unrolled kernel or, when none exists,
! to the general path over the same segments; `fock_nosym` sends every class
! through the shared kernel. So running this with f functions tests both halves
! of the routing at once: a dropped f class shows up as a disagreement, and so
! does an s..d class digested twice.
!
! Without the ceiling argument it is 6-31G, which stops at d and therefore says
! nothing about any of that.
!
program check_nosym
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_eri, only: trc_eri_t
   use trc_tables, only: LMAX
   use trc_test_basis, only: read_xyz, build_631g, build_aux
   implicit none

   character(len=256) :: xyzfile, arg
   integer :: nat, nsh, maxnp, nao, i, j, auxl
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
      print '(a)', 'usage: check_nosym <xyz> [aux-lmax]'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   ! An uncontracted l = 0..auxl set, because 6-31G stops at d and the point of
   ! the second argument is to reach past it.
   auxl = -1
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read (arg, *) auxl
   end if

   call read_xyz(xyzfile, nat, at_z, at_r)
   if (auxl >= 0) then
      call build_aux(nat, at_z, at_r, auxl, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   else
      call build_631g(nat, at_z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   end if
   if (maxval(sh_l) > LMAX) then
      print '(a,i0,a,i0,a)', '  SKIP: basis needs l = ', maxval(sh_l), &
         ' and this build has LMAX = ', LMAX, '.'
      stop 0
   end if

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

   print '(a,i0,a,i0,a,i0)', '  ', nao, ' functions, shells ', nsh, &
      ', max l ', maxval(sh_l)
   print '(a,es12.4)',  '  asymmetry of the density ', maxval(abs(d - transpose(d)))
   print '(a,es24.16)', '  sum folded     ', sum(gs)
   print '(a,es24.16)', '  sum enumerated ', sum(gn)
   print '(a,f18.10)',  '  ratio          ', sum(gn)/sum(gs)
   print '(a,es12.4)',  '  largest |G|                 ', maxval(abs(gn))
   print '(a,es12.4)',  '  worst |folded - enumerated| ', worst
   print '(a,es12.4)',  '  relative to largest |G|     ', worst/maxval(abs(gn))
   !
   ! RELATIVE, because the tolerance has to mean the same thing at both sizes
   ! it is run at. 6-31G on water reaches |G| ~ 20; the uncontracted l = 0..3
   ! set reaches ~1e3 over 140 functions, and holding that to the same ABSOLUTE
   ! 1e-10 demands 1e-13 relative -- below what summing that many terms in
   ! double precision returns. Scaling by the largest element asks both for the
   ! same number of correct digits.
   !
   if (worst < 1.0e-11_dp*max(1.0_dp, maxval(abs(gn)))) then
      print '(a)', '  RESULT: PASS'
   else
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if

   call eri%release(); call pl%release(); call bas%release()
end program check_nosym
