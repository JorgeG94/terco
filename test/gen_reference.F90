!
! Reference-integral generator.  Built with GFORTRAN, linked against libfint.
!
! WHY A SEPARATE PROGRAM
! ---------------------
! libfint cannot be built with nvfortran: `iso_fortran_env`'s real128 is -1
! there, so all 68 of its `real(qp)` declarations are hard errors (libfint's
! own CI keeps nvidia-hpc in the matrix as an expected failure and says so).
! The quad kind is the Rys fallback and is not droppable.
!
! So the oracle cannot be linked into an nvfortran test harness.  Instead this
! program runs once under gfortran, writes a plain-text dump of a system
! definition and every integral libfint says it has, and the nvfortran checks
! read that dump back.  Two consequences, both good: the tests need no second
! compiler at run time, and the reference data is a reviewable artefact rather
! than a number that appears from a library call.
!
! FORMAT
! ------
! Free-form text, one token per line group, deliberately dull so a human can
! read it and a two-line reader can parse it:
!
!   SYSTEM <name>
!   NATM <n> / NBAS <n> / NENV <n>
!   ATM  <ATM_SLOTS*natm integers>
!   BAS  <BAS_SLOTS*nbas integers>
!   ENV  <nenv doubles, %25.17e>
!   QUARTET <i> <j> <k> <l> <di> <dj> <dk> <dl>
!   <di*dj*dk*dl doubles, %25.17e, libcint's own ordering>
!   END
!
! %25.17e is round-trip exact for IEEE double, so the dump loses nothing.
!
program gen_reference
   use cint_const, only: dp
   use cint_bas, only: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                       BAS_SLOTS, CHARGE_OF, PTR_COORD, NUC_MOD_OF, ATM_SLOTS, &
                       POINT_NUC, cint_cgto_cart, cint_gto_norm
   use cint_workspace, only: cint_ws
   use cint_2e, only: int2e_cart
   implicit none

   integer, parameter :: PTR_ENV_START = 20
   integer, parameter :: MAXENV = 8000, MAXBAS = 64, MAXATM = 8

   integer, allocatable :: atm(:), bas(:)
   real(dp), allocatable :: env(:)
   integer :: natm, nbas, nenv
   type(cint_ws) :: ws
   integer :: u, isys

   open (newunit=u, file='reference.dat', status='replace', action='write')

   do isys = 1, 3
      call build_system(isys, atm, bas, env, natm, nbas, nenv)
      call dump_system(u, isys, atm, bas, env, natm, nbas, nenv, ws)
      deallocate (atm, bas, env)
   end do

   write (u, '(a)') 'END'
   close (u)
   print '(a)', 'wrote reference.dat'

contains

   !
   ! Two systems.  The first is the smallest thing that can be wrong; the
   ! second reaches d functions, which is where the register wall of the plan's
   ! section 0.2 actually bites and therefore the case worth being sure about.
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
         ! H2, two uncontracted s shells.  r = 1.4 bohr.
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
         ! Two centres carrying s, p and d, uncontracted, plus one contracted
         ! s so the primitive loop is exercised rather than trivially K=1.
         ! Total L across a quartet reaches 8, which is the Boys m_max the
         ! plan's section 3.1 says a through-d target needs.
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
            ! a genuinely contracted s shell
            nprim = 3
            e(1:3) = [6.0_dp, 1.3_dp, 0.42_dp]
            c(1:3) = [0.15_dp, 0.53_dp, 0.44_dp]
            call put_shell(bas, env, off, ib, ia, 0, e(1:3), c(1:3), nprim)
            ib = ib + 1
         end do

      case (3)
         ! s and p only, contracted -- the 6-31G regime.  Exists so an
         ! LMAX=1 build (which cannot represent system 2's d shells) still has
         ! something non-trivial to validate against.
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

   !
   ! One contracted shell.  libcint wants the coefficients in env already
   ! multiplied by the primitive normalisation, which is what cint_gto_norm
   ! returns -- passing raw coefficients is the classic way to get answers that
   ! are wrong by a smooth factor and look plausible.
   !
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

      integer :: i, j, k, l, di, dj, dk, dl, nq, shls(0:3), dims(0:3), m
      real(dp), allocatable :: buf(:)
      logical :: hv

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

      do i = 0, nbas - 1
      do j = 0, nbas - 1
      do k = 0, nbas - 1
      do l = 0, nbas - 1
         di = cint_cgto_cart(i, bas); dj = cint_cgto_cart(j, bas)
         dk = cint_cgto_cart(k, bas); dl = cint_cgto_cart(l, bas)
         nq = di*dj*dk*dl
         allocate (buf(nq)); buf = 0.0_dp
         shls(0) = i; shls(1) = j; shls(2) = k; shls(3) = l
         dims(0) = di; dims(1) = dj; dims(2) = dk; dims(3) = dl
         hv = int2e_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         if (hv) then
            write (u, '(a,8(1x,i0))') 'QUARTET', i, j, k, l, di, dj, dk, dl
            write (u, '(es25.17)') (buf(m), m=1, nq)
         end if
         deallocate (buf)
      end do
      end do
      end do
      end do
   end subroutine dump_system

end program gen_reference
