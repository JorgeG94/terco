!
! One-electron reference generator: overlap, kinetic, nuclear attraction.
!
! Built with GFORTRAN against libfint, for the reason gen_reference.F90 already
! documents: libfint cannot be compiled by nvfortran (real128 is -1 there, and
! all 68 of its real(qp) declarations are hard errors), so the oracle cannot be
! linked into an nvfortran test harness.  It runs once here and writes a
! plain-text dump the nvfortran checks read back.
!
! Same three systems as the ERI reference, so a failure can be compared against
! the two-electron result on identical geometry and basis:
!
!   1  H2, two uncontracted s shells
!   2  two centres carrying s, p and d, plus a contracted s
!   3  s and p only, contracted -- the 6-31G regime
!
! FORMAT
! ------
!   SYSTEM <n> / NATM / NBAS / NENV / ATM / BAS / ENV     (as gen_reference)
!   BLOCK <OVLP|KIN|NUC> <i> <j> <di> <dj>
!   <di*dj doubles, %25.17e, libcint ordering>
!   END
!
! %25.17e round-trips IEEE double exactly.
!
program gen_reference1e
   use cint_const, only: dp
   use cint_bas, only: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                       BAS_SLOTS, CHARGE_OF, PTR_COORD, NUC_MOD_OF, ATM_SLOTS, &
                       POINT_NUC, cint_cgto_cart, cint_gto_norm
   use cint_workspace, only: cint_ws
   use cint_1e, only: int1e_ovlp_cart, int1e_nuc_cart
   use cint_gen_intor1, only: int1e_kin_cart
   implicit none

   integer, parameter :: PTR_ENV_START = 20
   integer, parameter :: MAXENV = 8000, MAXBAS = 64, MAXATM = 8

   integer, allocatable :: atm(:), bas(:)
   real(dp), allocatable :: env(:)
   integer :: natm, nbas, nenv
   type(cint_ws) :: ws
   integer :: u, isys

   open (newunit=u, file='reference1e.dat', status='replace', action='write')

   do isys = 1, 3
      call build_system(isys, atm, bas, env, natm, nbas, nenv)
      call dump_system(u, isys, atm, bas, env, natm, nbas, nenv, ws)
      deallocate (atm, bas, env)
   end do

   write (u, '(a)') 'END'
   close (u)
   print '(a)', 'wrote reference1e.dat'

contains

   !
   ! Identical to gen_reference.F90's systems.  Duplicated rather than shared
   ! so the two oracles cannot drift apart silently through one edit.
   !
   subroutine build_system(which, atm, bas, env, natm, nbas, nenv)
      integer, intent(in) :: which
      integer, allocatable, intent(out) :: atm(:), bas(:)
      real(dp), allocatable, intent(out) :: env(:)
      integer, intent(out) :: natm, nbas, nenv

      integer :: off, ib, ia, l, nprim
      real(dp) :: e(4), c(4)

      allocate (atm(0:ATM_SLOTS*MAXATM - 1), bas(0:BAS_SLOTS*MAXBAS - 1))
      allocate (env(0:MAXENV - 1))
      atm = 0; bas = 0; env = 0.0_dp
      off = PTR_ENV_START

      select case (which)

      case (1)
         natm = 2; nbas = 2
         do ia = 0, 1
            atm(ATM_SLOTS*ia + CHARGE_OF) = 1
            atm(ATM_SLOTS*ia + PTR_COORD) = off
            atm(ATM_SLOTS*ia + NUC_MOD_OF) = POINT_NUC
            env(off) = 0.0_dp; env(off + 1) = 0.0_dp
            env(off + 2) = 1.4_dp*real(ia, dp)
            off = off + 3
         end do
         do ib = 0, 1
            call put_shell(bas, env, off, ib, ib, 0, [1.0_dp], [1.0_dp], 1)
         end do

      case (2)
         natm = 2; nbas = 8
         do ia = 0, 1
            atm(ATM_SLOTS*ia + CHARGE_OF) = 4
            atm(ATM_SLOTS*ia + PTR_COORD) = off
            atm(ATM_SLOTS*ia + NUC_MOD_OF) = POINT_NUC
            env(off) = 0.0_dp; env(off + 1) = 0.0_dp
            env(off + 2) = 1.8_dp*real(ia, dp)
            off = off + 3
         end do
         ib = 0
         do ia = 0, 1
            do l = 0, 2
               e(1) = 0.8_dp + 0.35_dp*real(l, dp) + 0.11_dp*real(ia, dp)
               call put_shell(bas, env, off, ib, ia, l, e(1:1), [1.0_dp], 1)
               ib = ib + 1
            end do
            nprim = 3
            e(1:3) = [6.0_dp, 1.3_dp, 0.42_dp]
            c(1:3) = [0.15_dp, 0.53_dp, 0.44_dp]
            call put_shell(bas, env, off, ib, ia, 0, e(1:3), c(1:3), nprim)
            ib = ib + 1
         end do

      case (3)
         natm = 2; nbas = 6
         do ia = 0, 1
            atm(ATM_SLOTS*ia + CHARGE_OF) = 4
            atm(ATM_SLOTS*ia + PTR_COORD) = off
            atm(ATM_SLOTS*ia + NUC_MOD_OF) = POINT_NUC
            env(off) = 0.0_dp; env(off + 1) = 0.0_dp
            env(off + 2) = 1.7_dp*real(ia, dp)
            off = off + 3
         end do
         ib = 0
         do ia = 0, 1
            e(1) = 0.9_dp + 0.2_dp*real(ia, dp)
            call put_shell(bas, env, off, ib, ia, 0, e(1:1), [1.0_dp], 1)
            ib = ib + 1
            e(1) = 0.55_dp + 0.13_dp*real(ia, dp)
            call put_shell(bas, env, off, ib, ia, 1, e(1:1), [1.0_dp], 1)
            ib = ib + 1
            nprim = 3
            e(1:3) = [5.2_dp, 1.1_dp, 0.38_dp]
            c(1:3) = [0.17_dp, 0.51_dp, 0.46_dp]
            call put_shell(bas, env, off, ib, ia, 0, e(1:3), c(1:3), nprim)
            ib = ib + 1
         end do

      end select

      nenv = off
   end subroutine build_system

   subroutine put_shell(bas, env, off, ib, iatom, l, e, c, nprim)
      integer, intent(inout) :: bas(0:)
      real(dp), intent(inout) :: env(0:)
      integer, intent(inout) :: off
      integer, intent(in) :: ib, iatom, l, nprim
      real(dp), intent(in) :: e(:), c(:)
      integer :: k

      bas(BAS_SLOTS*ib + ATOM_OF) = iatom
      bas(BAS_SLOTS*ib + ANG_OF) = l
      bas(BAS_SLOTS*ib + NPRIM_OF) = nprim
      bas(BAS_SLOTS*ib + NCTR_OF) = 1
      bas(BAS_SLOTS*ib + PTR_EXP) = off
      do k = 1, nprim
         env(off + k - 1) = e(k)
      end do
      off = off + nprim
      bas(BAS_SLOTS*ib + PTR_COEFF) = off
      do k = 1, nprim
         env(off + k - 1) = c(k)*cint_gto_norm(l, e(k))
      end do
      off = off + nprim
   end subroutine put_shell

   subroutine dump_system(u, isys, atm, bas, env, natm, nbas, nenv, ws)
      integer, intent(in) :: u, isys, natm, nbas, nenv
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      type(cint_ws), intent(inout) :: ws

      integer :: i, di

      write (u, '(a,i0)') 'SYSTEM ', isys
      write (u, '(a,i0)') 'NATM ', natm
      write (u, '(a,i0)') 'NBAS ', nbas
      write (u, '(a,i0)') 'NENV ', nenv
      write (u, '(a)') 'ATM'
      write (u, '(i0)') (atm(i), i=0, ATM_SLOTS*natm - 1)
      write (u, '(a)') 'BAS'
      write (u, '(i0)') (bas(i), i=0, BAS_SLOTS*nbas - 1)
      write (u, '(a)') 'ENV'
      write (u, '(es25.17)') (env(i), i=0, nenv - 1)

      call dump_kind(u, 'OVLP', 1, atm, bas, env, natm, nbas, ws)
      call dump_kind(u, 'KIN ', 2, atm, bas, env, natm, nbas, ws)
      call dump_kind(u, 'NUC ', 3, atm, bas, env, natm, nbas, ws)
      di = 0
   end subroutine dump_system

   subroutine dump_kind(u, name, which, atm, bas, env, natm, nbas, ws)
      integer, intent(in) :: u, which, natm, nbas
      character(len=*), intent(in) :: name
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      type(cint_ws), intent(inout) :: ws

      integer :: i, j, di, dj, n, m, shls(0:1), dims(0:1)
      real(dp), allocatable :: buf(:)
      logical :: hv

      do i = 0, nbas - 1
      do j = 0, nbas - 1
         di = cint_cgto_cart(i, bas); dj = cint_cgto_cart(j, bas)
         n = di*dj
         allocate (buf(n)); buf = 0.0_dp
         shls(0) = i; shls(1) = j
         dims(0) = di; dims(1) = dj
         select case (which)
         case (1)
            hv = int1e_ovlp_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (2)
            hv = int1e_kin_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (3)
            hv = int1e_nuc_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         end select
         if (hv) then
            write (u, '(a,a,4(1x,i0))') 'BLOCK ', trim(name), i, j, di, dj
            write (u, '(es25.17)') (buf(m), m=1, n)
         end if
         deallocate (buf)
      end do
      end do
   end subroutine dump_kind

end program gen_reference1e
