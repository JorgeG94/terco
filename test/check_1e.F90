!
! One-electron kernels against libfint: overlap, kinetic, nuclear attraction.
!
! Reads reference1e.dat (written once by gen_reference1e under gfortran, since
! libfint cannot be built with nvfortran) and compares every shell-pair block
! for all three reference systems -- including system 2, which carries d
! shells.
!
! NORMALISATION
! -------------
! libcint's `env` already holds contraction coefficients multiplied by the
! primitive norm, which `put_shell` did on the reference side. What it does NOT
! hold is `common_fac_sp`, the per-shell factor libcint applies separately:
! 1/sqrt(4 pi) for s, sqrt(3/(4 pi)) for p, 1 beyond. The four-centre path was
! wrong by exactly this factor once already; it is applied here per shell, so
! twice per block.
!
! Tolerance is mixed, absolute plus relative, for the reason the ERI check
! documents: dividing a 1e-16 difference by a 1e-12 reference manufactures
! failures out of exact results.
!
! Checks are grouped one per (system, kind) -- kind is OVLP/KIN/NUC -- each
! reporting that system's worst diff, rather than one check per matrix
! element: the file carries every shell-pair block for three systems, which
! would otherwise be thousands of individual check() calls.
!
module check_1e_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_1e_kernels, only: one_e_dispatch, ONE_E_LMAX
   use trc_cart_perm, only: cart_perm_b2l
   implicit none
   private
   public :: collect_1e_tests

   integer, parameter :: ATM_SLOTS = 6, BAS_SLOTS = 8
   integer, parameter :: ATOM_OF = 0, ANG_OF = 1, NPRIM_OF = 2
   integer, parameter :: PTR_EXP = 5, PTR_COEFF = 6, PTR_COORD = 1

contains

   subroutine collect_1e_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overlap/kinetic/nuclear kernels match libfint, worst diff per system/kind", &
                               test_vs_libfint) &
                  ]
   end subroutine collect_1e_tests

   subroutine test_vs_libfint(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: atm(:), bas(:)
      real(dp), allocatable :: env(:)
      integer :: natm, nbas, nenv, u, isys
      character(len=64) :: t

      open (newunit=u, file='reference1e.dat', status='old', action='read')
      do
         read (u, '(a)', end=99) t
         if (t(1:3) == 'END') exit
         if (t(1:6) /= 'SYSTEM') cycle
         read (t(8:), *) isys
         call read_system(u, atm, bas, env, natm, nbas, nenv)
         call check_system(u, isys, atm, bas, env, natm, nbas, error)
         deallocate (atm, bas, env)
         if (allocated(error)) then
            close (u)
            return
         end if
      end do
99    continue
      close (u)
   end subroutine test_vs_libfint

   subroutine read_system(u, atm, bas, env, natm, nbas, nenv)
      integer, intent(in) :: u
      integer, allocatable, intent(out) :: atm(:), bas(:)
      real(dp), allocatable, intent(out) :: env(:)
      integer, intent(out) :: natm, nbas, nenv
      integer :: i
      character(len=64) :: t

      read (u, '(a)') t; read (t(6:), *) natm
      read (u, '(a)') t; read (t(6:), *) nbas
      read (u, '(a)') t; read (t(6:), *) nenv
      allocate (atm(0:ATM_SLOTS*natm - 1), bas(0:BAS_SLOTS*nbas - 1), env(0:nenv - 1))
      read (u, '(a)') t; read (u, *) (atm(i), i=0, ATM_SLOTS*natm - 1)
      read (u, '(a)') t; read (u, *) (bas(i), i=0, BAS_SLOTS*nbas - 1)
      read (u, '(a)') t; read (u, *) (env(i), i=0, nenv - 1)
   end subroutine read_system

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   pure real(dp) function common_fac_sp(l)
      !! libcint's per-shell factor, applied outside `env`.
      integer, intent(in) :: l
      select case (l)
      case (0); common_fac_sp = 0.282094791773878143_dp
      case (1); common_fac_sp = 0.488602511902919921_dp
      case default; common_fac_sp = 1.0_dp
      end select
   end function common_fac_sp

   subroutine check_system(u, isys, atm, bas, env, natm, nbas, error)
      integer, intent(in) :: u, isys, natm, nbas
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      type(error_type), allocatable, intent(out) :: error

      character(len=80) :: line
      character(len=8)  :: kind
      integer :: i, j, di, dj, n, m, ml, k, la, lb, npa, npb, ic
      real(dp), allocatable :: ref(:), sout(:), tout(:), vout(:)
      real(dp), allocatable :: zatm(:), ratm(:, :)
      real(dp) :: ra(3), rb(3), d, tol, fac
      real(dp) :: worst(3)
      integer :: nbad(3)
      integer :: nblk

      allocate (zatm(natm), ratm(3, natm))
      do ic = 1, natm
         zatm(ic) = real(atm(ATM_SLOTS*(ic - 1) + 0), dp)
         ratm(:, ic) = env(atm(ATM_SLOTS*(ic - 1) + PTR_COORD): &
                           atm(ATM_SLOTS*(ic - 1) + PTR_COORD) + 2)
      end do

      worst = 0.0_dp
      nbad = 0
      nblk = 0
      do
         read (u, '(a)', end=98) line
         if (line(1:3) == 'END') then
            backspace (u); exit
         end if
         if (line(1:6) == 'SYSTEM') then
            backspace (u); exit
         end if
         if (line(1:5) /= 'BLOCK') cycle
         read (line(7:), *) kind, i, j, di, dj
         n = di*dj
         allocate (ref(n))
         read (u, *) (ref(m), m=1, n)

         la = bas(BAS_SLOTS*i + ANG_OF)
         lb = bas(BAS_SLOTS*j + ANG_OF)
         if (max(la, lb) > ONE_E_LMAX) then
            deallocate (ref); cycle
         end if
         npa = bas(BAS_SLOTS*i + NPRIM_OF)
         npb = bas(BAS_SLOTS*j + NPRIM_OF)
         ra = env(atm(ATM_SLOTS*bas(BAS_SLOTS*i + 0) + PTR_COORD): &
                  atm(ATM_SLOTS*bas(BAS_SLOTS*i + 0) + PTR_COORD) + 2)
         rb = env(atm(ATM_SLOTS*bas(BAS_SLOTS*j + 0) + PTR_COORD): &
                  atm(ATM_SLOTS*bas(BAS_SLOTS*j + 0) + PTR_COORD) + 2)

         allocate (sout(n), tout(n), vout(n))
         call one_e_dispatch(la, lb, npa, npb, &
              env(bas(BAS_SLOTS*i + PTR_EXP):bas(BAS_SLOTS*i + PTR_EXP) + npa - 1), &
              env(bas(BAS_SLOTS*i + PTR_COEFF):bas(BAS_SLOTS*i + PTR_COEFF) + npa - 1), &
              env(bas(BAS_SLOTS*j + PTR_EXP):bas(BAS_SLOTS*j + PTR_EXP) + npb - 1), &
              env(bas(BAS_SLOTS*j + PTR_COEFF):bas(BAS_SLOTS*j + PTR_COEFF) + npb - 1), &
              ra, rb, natm, zatm, ratm, sout, tout, vout)

         fac = common_fac_sp(la)*common_fac_sp(lb)
         do m = 1, n
            ! sout/tout/vout(m) is this BUILD's own flat index
            ! (1 + ia_build + di*ib_build); ref(m) is libfint's, LEXICOGRAPHIC
            ! (1 + ia_lex + di*ib_lex). Decode m's build-order (ia,ib),
            ! convert each to lexicographic, recombine, and read `ref` at
            ! THAT index -- the two flat spaces agree only at the default
            ! order, where cart_perm_b2l is the identity.
            ml = 1 + cart_perm_b2l(la, mod(m - 1, di)) &
                   + di*cart_perm_b2l(lb, (m - 1)/di)
            select case (trim(kind))
            case ('OVLP'); d = abs(fac*sout(m) - ref(ml)); k = 1
            case ('KIN');  d = abs(fac*tout(m) - ref(ml)); k = 2
            case ('NUC');  d = abs(fac*vout(m) - ref(ml)); k = 3
            case default;  d = 0.0_dp; k = 1
            end select
            tol = 1.0e-12_dp + 1.0e-11_dp*abs(ref(ml))
            worst(k) = max(worst(k), d)
            if (d > tol) nbad(k) = nbad(k) + 1
         end do
         nblk = nblk + 1
         deallocate (ref, sout, tout, vout)
      end do
98    continue
      print '(a,i0,a,i0,a)', '  system ', isys, ': ', nblk, ' blocks checked'
      deallocate (zatm, ratm)

      call check(error, nbad(1) == 0, &
                 "system "//to_string(isys)//" OVLP: worst diff "//to_string(worst(1), "es10.2")// &
                 " exceeds mixed tol ("//to_string(nbad(1))//" entries)")
      if (allocated(error)) return
      call check(error, nbad(2) == 0, &
                 "system "//to_string(isys)//" KIN: worst diff "//to_string(worst(2), "es10.2")// &
                 " exceeds mixed tol ("//to_string(nbad(2))//" entries)")
      if (allocated(error)) return
      call check(error, nbad(3) == 0, &
                 "system "//to_string(isys)//" NUC: worst diff "//to_string(worst(3), "es10.2")// &
                 " exceeds mixed tol ("//to_string(nbad(3))//" entries)")
   end subroutine check_system

end module check_1e_suite

program check_1e
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_boys, only: boys_init
   use check_1e_suite, only: collect_1e_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   call boys_init()
   stat = 0
   testsuites = [new_testsuite("check_1e", collect_1e_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)

   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_1e
