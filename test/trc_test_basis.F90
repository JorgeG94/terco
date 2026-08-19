!
! Basis construction for the test programs.
!
! The library does NOT ship a basis set -- a caller supplies exponents and
! coefficients. That is the right contract, but it means every test program
! needs its own basis, and three of them wanting the same 6-31G is duplication
! rather than independence. This module is that shared code, and it lives in
! `test/` rather than `src/` precisely because it is not part of the library.
!
module trc_test_basis
   use trc_boys, only: dp
   implicit none
   public

   real(dp), parameter :: ANG = 1.8897261254578281_dp

contains

   subroutine read_xyz(fn, nat, z, r)
      character(len=*), intent(in) :: fn
      integer, intent(out) :: nat
      integer,  allocatable, intent(out) :: z(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      integer :: u, i, ios
      character(len=8) :: sym
      character(len=256) :: line
      real(dp) :: x, y, zz
      open (newunit=u, file=fn, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print '(a,a)', '  cannot open ', trim(fn); stop 1
      end if
      read (u, *) nat
      read (u, '(a)') line
      allocate (z(nat), r(3, nat))
      do i = 1, nat
         read (u, *) sym, x, y, zz
         r(:, i) = [x, y, zz]*ANG
         select case (trim(adjustl(sym)))
         case ('H', 'h'); z(i) = 1
         case ('C', 'c'); z(i) = 6
         case ('N', 'n'); z(i) = 7
         case ('O', 'o'); z(i) = 8
         case default
            print '(a,a)', '  unsupported element ', trim(sym); stop 1
         end select
      end do
      close (u)
   end subroutine read_xyz

   !> 6-31G for H, C, N, O. Same values as bench_binwater, from the BSE JSON.
   subroutine build_631g(nat, z, r, nsh, l, np, e, c, sr, mnp)
      integer,  intent(in)  :: nat, z(nat)
      real(dp), intent(in)  :: r(3, nat)
      integer,  intent(out) :: nsh, mnp
      integer,  allocatable, intent(out) :: l(:), np(:)
      real(dp), allocatable, intent(out) :: e(:, :), c(:, :), sr(:, :)

      integer :: i, k
      real(dp) :: e6(6), c6(6), e3(3), cs3(3), cp3(3), e1(1)

      nsh = 0
      do i = 1, nat
         if (z(i) > 2) then
            nsh = nsh + 5
         else
            nsh = nsh + 2
         end if
      end do
      mnp = 6
      allocate (l(nsh), np(nsh), e(mnp, nsh), c(mnp, nsh), sr(3, nsh))
      e = 0.0_dp; c = 0.0_dp
      k = 0
      do i = 1, nat
         select case (z(i))
         case (1)
            call put(k, 0, 3, [18.7311370_dp, 2.8253944_dp, 0.6401217_dp], &
                     [0.0334946_dp, 0.2347270_dp, 0.8137573_dp], r(:, i), &
                     l, np, e, c, sr, mnp)
            call put(k, 0, 1, [0.1612778_dp], [1.0_dp], r(:, i), &
                     l, np, e, c, sr, mnp)
            cycle
         case (6)
            e6 = [3047.5248800_dp, 457.3695180_dp, 103.9486850_dp, &
                  29.2101553_dp, 9.2866630_dp, 3.1639270_dp]
            c6 = [0.0018347_dp, 0.0140373_dp, 0.0688426_dp, &
                  0.2321844_dp, 0.4679413_dp, 0.3623120_dp]
            e3 = [7.8682723_dp, 1.8812885_dp, 0.5442493_dp]
            cs3 = [-0.1193324_dp, -0.1608542_dp, 1.1434564_dp]
            cp3 = [0.0689991_dp, 0.3164240_dp, 0.7443083_dp]
            e1 = [0.1687145_dp]
         case (7)
            e6 = [4173.5114600_dp, 627.4579110_dp, 142.9020930_dp, &
                  40.2343293_dp, 12.8202129_dp, 4.3904370_dp]
            c6 = [0.0018348_dp, 0.0139946_dp, 0.0685866_dp, &
                  0.2322409_dp, 0.4690699_dp, 0.3604552_dp]
            e3 = [11.6263619_dp, 2.7162798_dp, 0.7722184_dp]
            cs3 = [-0.1149612_dp, -0.1691175_dp, 1.1458519_dp]
            cp3 = [0.0675797_dp, 0.3239073_dp, 0.7408951_dp]
            e1 = [0.2120315_dp]
         case (8)
            e6 = [5484.6716600_dp, 825.2349460_dp, 188.0469580_dp, &
                  52.9645000_dp, 16.8975704_dp, 5.7996353_dp]
            c6 = [0.0018311_dp, 0.0139502_dp, 0.0684451_dp, &
                  0.2327143_dp, 0.4701929_dp, 0.3585209_dp]
            e3 = [15.5396162_dp, 3.5999336_dp, 1.0137618_dp]
            cs3 = [-0.1107775_dp, -0.1480263_dp, 1.1307670_dp]
            cp3 = [0.0708743_dp, 0.3397528_dp, 0.7271586_dp]
            e1 = [0.2700058_dp]
         case default
            print '(a,i0)', '  no 6-31G for Z = ', z(i); stop 1
         end select
         call put(k, 0, 6, e6, c6, r(:, i), l, np, e, c, sr, mnp)
         call put(k, 0, 3, e3, cs3, r(:, i), l, np, e, c, sr, mnp)
         call put(k, 1, 3, e3, cp3, r(:, i), l, np, e, c, sr, mnp)
         call put(k, 0, 1, e1, [1.0_dp], r(:, i), l, np, e, c, sr, mnp)
         call put(k, 1, 1, e1, [1.0_dp], r(:, i), l, np, e, c, sr, mnp)
      end do
   end subroutine build_631g

   !
   ! Synthetic auxiliary basis: uncontracted shells per atom, l = 0..lmax,
   ! exponents on a geometric ladder. Heavier atoms get one more shell per l,
   ! which is roughly how real fitting sets scale.
   !
   subroutine build_aux(nat, z, r, lmax, nsh, l, np, e, c, sr, mnp)
      integer,  intent(in)  :: nat, z(nat), lmax
      real(dp), intent(in)  :: r(3, nat)
      integer,  intent(out) :: nsh, mnp
      integer,  allocatable, intent(out) :: l(:), np(:)
      real(dp), allocatable, intent(out) :: e(:, :), c(:, :), sr(:, :)

      integer :: i, ll, m, k, per
      real(dp) :: ex

      nsh = 0
      do i = 1, nat
         per = merge(3, 2, z(i) > 2)
         nsh = nsh + per*(lmax + 1)
      end do
      mnp = 1
      allocate (l(nsh), np(nsh), e(mnp, nsh), c(mnp, nsh), sr(3, nsh))
      k = 0
      do i = 1, nat
         per = merge(3, 2, z(i) > 2)
         do ll = 0, lmax
            do m = 0, per - 1
               ex = 0.4_dp*(3.5_dp**m)*(1.0_dp + 0.6_dp*real(ll, dp))
               if (z(i) > 2) ex = ex*2.0_dp
               call put(k, ll, 1, [ex], [1.0_dp], r(:, i), l, np, e, c, sr, mnp)
            end do
         end do
      end do
   end subroutine build_aux

   subroutine put(k, ll, n, ex, co, rr, l, np, e, c, sr, mnp)
      integer,  intent(inout) :: k
      integer,  intent(in)    :: ll, n, mnp
      real(dp), intent(in)    :: ex(:), co(:), rr(3)
      integer,  intent(inout) :: l(:), np(:)
      ! explicit shape, not `e(mnp, :)`: Fortran forbids mixing a fixed
      ! extent with an assumed one in the same dummy
      real(dp), intent(inout) :: e(mnp, size(l)), c(mnp, size(l)), sr(3, size(l))
      integer :: q
      k = k + 1
      l(k) = ll; np(k) = n; sr(:, k) = rr
      do q = 1, n
         e(q, k) = ex(q)
         c(q, k) = co(q)*gto_norm(ll, ex(q))
      end do
   end subroutine put

   !> libcint's primitive normalisation, so the caller's coefficients mean
   !> what a basis-set file means.
   pure real(dp) function gto_norm(l, a)
      integer,  intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: t
      t = 2.0_dp**(2*l + 3)*fact(l + 1)*(2.0_dp*a)**(l + 1.5_dp) &
          /(fact(2*l + 2)*sqrt(3.14159265358979323846_dp))
      gto_norm = sqrt(t)
   end function gto_norm

   pure real(dp) function fact(n)
      integer, intent(in) :: n
      integer :: i
      fact = 1.0_dp
      do i = 2, n
         fact = fact*real(i, dp)
      end do
   end function fact

end module trc_test_basis
