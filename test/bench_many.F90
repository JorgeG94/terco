!
! Batched multi-density Fock: correctness against the loop, and where it stops
! paying.
!
! WHY THIS MATTERS MORE THAN THE SINGLE-BUILD NUMBER
! --------------------------------------------------
! makeEFP is not one Fock build, it is about a hundred: the coupled-perturbed
! equations for the dynamic polarizabilities need nine perturbations times
! twelve imaginary frequencies, each a Fock build on a different response
! density. Contracting one quartet against every density in hand is the
! difference between one integral pass and N of them.
!
! WHAT TO EXPECT
! --------------
! Cost is `eval + N*digest`, so the ceiling is one over the digestion
! fraction. terco's digestion measured 24% (Gly30) to 37% (RNA3), which puts
! the ceiling between 2.7x and 4.2x -- approached slowly, since N has to be
! large before `N*digest` dominates `eval`.
!
! The GPU constraint that mqc's CPU version does not have is register
! pressure. This benchmark exists to find out whether it bites: the density
! loop is OUTSIDE the block accumulators precisely so the batch does not add
! registers, and that claim should be checked rather than believed.
!
!   bench_many <xyz> <density-file> [nmax] [thresh]
!
program bench_many
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_eri, only: trc_eri_t
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   character(len=256) :: xyzfile, dmfile, arg
   integer  :: nmax
   real(dp) :: thresh

   integer :: nat, nsh, maxnp, nao, i, j, d, n, r
   integer,  allocatable :: at_z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), zat(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: d0(:, :), dmats(:, :, :), gmany(:, :, :), gone(:, :)
   type(trc_basis_t)    :: bas
   type(trc_pairlist_t) :: pl
   type(trc_eri_t)      :: eri
   real(dp) :: t0, t1, t_one, t_many, worst

   nmax = 16; thresh = 1.0e-10_dp
   if (command_argument_count() < 2) then
      print '(a)', 'usage: bench_many <xyz> <density-file> [nmax] [thresh]'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, dmfile)
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read (arg, *) nmax
   end if
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg); read (arg, *) thresh
   end if

   call read_xyz(xyzfile, nat, at_z, at_r)
   call build_631g(nat, at_z, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   allocate (zat(nat)); zat = real(at_z, dp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, zat, at_r, maxnp)
   call pl%build(bas, thresh)
   nao = bas%nao

   allocate (d0(nao, nao))
   call read_density(dmfile, nao, d0)

   call bas%to_device()
   call pl%to_device()
   call eri%build(bas, thresh)

   print '(a,a,a,i0,a,i0)', '  ', trim(xyzfile), '   nao ', nao, '   shells ', nsh
   print '(a,es9.2)', '  threshold ', thresh
   print '(a)', ''
   print '(a)', '     N   t(batched)   t(N separate)   speedup   worst |diff|'

   allocate (gone(nao, nao))
   n = 1
   do while (n <= nmax)
      allocate (dmats(n, nao, nao), gmany(n, nao, nao))
      !
      ! Perturbed densities, deterministic and all different. Scaled down so
      ! they stay physically plausible in magnitude; the point is that no two
      ! sets are the same, so a batched kernel that silently reused one
      ! density would be caught.
      !
      do d = 1, n
         do j = 1, nao
            do i = 1, nao
               dmats(d, i, j) = d0(i, j)*(1.0_dp + 0.01_dp*real(d, dp)) &
                                + 1.0e-3_dp*cos(real(i + 2*j + 7*d, dp))
            end do
         end do
         ! the kernel contracts a symmetric density
         dmats(d, :, :) = 0.5_dp*(dmats(d, :, :) + transpose(dmats(d, :, :)))
      end do

      call eri%fock_many(bas, n, dmats, gmany)      ! warm
      call tick(t0)
      call eri%fock_many(bas, n, dmats, gmany)
      call tick(t1)
      t_many = t1 - t0

      call tick(t0)
      do d = 1, n
         call eri%fock(bas, dmats(d, :, :), gone)
      end do
      call tick(t1)
      t_one = t1 - t0

      ! and the batch must agree with the loop, set by set
      worst = 0.0_dp
      do d = 1, n
         call eri%fock(bas, dmats(d, :, :), gone)
         do j = 1, nao
            do i = 1, nao
               worst = max(worst, abs(gmany(d, i, j) - gone(i, j)))
            end do
         end do
      end do

      print '(i6,f13.4,f16.4,f10.2,a,es15.3)', n, t_many, t_one, &
         t_one/max(t_many, 1.0e-12_dp), 'x', worst
      deallocate (dmats, gmany)
      n = n*2
   end do

   call eri%release(); call pl%release(); call bas%release()

contains

   subroutine tick(t)
      real(dp), intent(out) :: t
      integer(kind=8) :: c, rate
      call system_clock(c, rate)
      t = real(c, dp)/real(rate, dp)
   end subroutine tick

   subroutine read_density(fn, nao, d)
      character(len=*), intent(in) :: fn
      integer, intent(in) :: nao
      real(dp), intent(out) :: d(nao, nao)
      integer :: u, ios
      integer(kind=8) :: nchk
      open (newunit=u, file=fn, access='stream', form='unformatted', &
            status='old', iostat=ios)
      if (ios /= 0) then
         print '(a,a)', '  cannot open ', trim(fn); stop 1
      end if
      read (u) nchk
      if (int(nchk) /= nao) then
         print '(a,i0,a,i0)', '  density nao ', nchk, ' /= ', nao; stop 1
      end if
      read (u) d
      close (u)
   end subroutine read_density

end program bench_many
