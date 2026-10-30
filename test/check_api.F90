!
! The batched API against libfint, end to end.
!
! check_1e and check_df pin the KERNELS: they call one class at a time with
! hand-assembled arguments. This pins everything between a caller and those
! kernels -- the libcint-layout import, the AO offsets, device residency, and
! the `do concurrent` drivers -- by building the basis object the way a
! consumer would and comparing whole matrices against the same references.
!
! It is a separate check for the reason check_binfock is separate from
! check_batch: if the container mislabelled a shell, kernel-level tests would
! still pass, because they never use the container.
!
module check_api_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_df_3c, &
                      trc_multipoles, TRC_NMULT
   use trc_cart_perm, only: cart_perm_l2b
   implicit none
   private
   public :: collect_api_tests

   integer, parameter :: ATM_SLOTS = 6, BAS_SLOTS = 8
   integer, parameter :: ANG_OF = 1

contains

   subroutine collect_api_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("one-electron matrices (OVLP/KIN/NUC) match libfint, every system", &
                               test_1e_vs_libfint), &
                  new_unittest("three-centre DF tensor matches libfint, every system", test_3c_vs_libfint) &
                  ]
   end subroutine collect_api_tests

   subroutine test_1e_vs_libfint(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: atm(:), bas(:)
      real(dp), allocatable :: env(:)
      integer :: natm, nbas, nenv, u, isys, nbad
      real(dp) :: w1(3)
      character(len=80) :: t

      nbad = 0; w1 = 0.0_dp

      open (newunit=u, file='reference1e.dat', status='old', action='read')
      do
         read (u, '(a)', end=91) t
         if (t(1:3) == 'END') exit
         if (t(1:6) /= 'SYSTEM') cycle
         read (t(8:), *) isys
         call read_system(u, atm, bas, env, natm, nbas, nenv)
         call check_1e_system(u, isys, atm, bas, env, natm, nbas, w1, nbad)
         deallocate (atm, bas, env)
      end do
91    continue
      close (u)

      call check(error, nbad == 0, &
                 "1e API vs libfint: worst overlap "//to_string(w1(1), "es10.2")// &
                 ", kinetic "//to_string(w1(2), "es10.2")//", nuclear "//to_string(w1(3), "es10.2")// &
                 "; "//to_string(nbad)//" block entries outside tol")
   end subroutine test_1e_vs_libfint

   subroutine test_3c_vs_libfint(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: atm(:), bas(:)
      real(dp), allocatable :: env(:)
      integer :: natm, nbas, nenv, u, isys, nbad
      real(dp) :: w3
      character(len=80) :: t

      nbad = 0; w3 = 0.0_dp

      open (newunit=u, file='reference_df.dat', status='old', action='read')
      do
         read (u, '(a)', end=92) t
         if (t(1:3) == 'END') exit
         if (t(1:6) /= 'SYSTEM') cycle
         read (t(8:), *) isys
         call read_system(u, atm, bas, env, natm, nbas, nenv)
         call check_3c_system(u, isys, atm, bas, env, natm, nbas, w3, nbad)
         deallocate (atm, bas, env)
      end do
92    continue
      close (u)

      call check(error, nbad == 0, &
                 "3-centre API vs libfint: worst diff "//to_string(w3, "es10.2")// &
                 "; "//to_string(nbad)//" tensor entries outside tol")
   end subroutine test_3c_vs_libfint

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

   subroutine ao_offsets(bas, nbas, off, nao)
      !! AO offset of each shell, computed the way the container does.
      integer, intent(in)  :: bas(0:), nbas
      integer, intent(out) :: off(nbas), nao
      integer :: i, l
      nao = 0
      do i = 1, nbas
         l = bas(BAS_SLOTS*(i - 1) + ANG_OF)
         off(i) = nao + 1
         nao = nao + (l + 1)*(l + 2)/2
      end do
   end subroutine ao_offsets

   subroutine check_1e_system(u, isys, atm, bas, env, natm, nbas, w1, nbad)
      integer, intent(in) :: u, isys, natm, nbas
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      real(dp), intent(inout) :: w1(3)
      integer, intent(inout) :: nbad

      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: pl
      character(len=80) :: line
      character(len=8)  :: kind
      integer :: i, j, di, dj, n, m, a, c, mu, nu, k, nao, nkept
      integer, allocatable :: off(:)
      real(dp), allocatable :: ref(:), smat(:, :), tmat(:, :), vmat(:, :)
      real(dp) :: d, tol

      allocate (off(nbas))
      call ao_offsets(bas, nbas, off, nao)

      call b%from_libcint(atm, natm, bas, nbas, env)
      ! A correctness check must screen nothing, so the threshold is far below
      ! anything the tolerance could notice.
      call pl%build(b, 1.0e-30_dp)
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      !$acc enter data create(smat, tmat, vmat)
      call b%to_device()
      call pl%to_device()
      call trc_1e(b, pl, smat, tmat, vmat)
      !$acc update self(smat, tmat, vmat)
      !$acc exit data delete(smat, tmat, vmat)
      ! Multipoles about the origin, dumped for pyscf to check. libfint has no
      ! plain multipole entry point, so the oracle is pyscf's int1e_r,
      ! int1e_rr and int1e_rrr -- an outside code, which is the point.
      if (isys == 2) then
         block
            real(dp), allocatable :: mm(:, :, :)
            integer :: um
            allocate (mm(nao, nao, TRC_NMULT))
            !$acc enter data create(mm)
            call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mm)
            !$acc update self(mm)
            !$acc exit data delete(mm)
            open (newunit=um, file='multipole_probe.bin', access='stream', &
                  form='unformatted', status='replace')
            write (um) int(nao, int64), int(TRC_NMULT, int64)
            write (um) mm
            ! S alongside, as the control: if pyscf's rebuild of this system
            ! disagrees on the OVERLAP, the fault is in the test system and
            ! not in the multipole kernel.
            write (um) smat
            close (um)
            deallocate (mm)
         end block
      end if

      nkept = pl%npair        ! before release, which zeroes it
      call pl%release()
      call b%release()

      do
         read (u, '(a)', end=98) line
         if (line(1:3) == 'END' .or. line(1:6) == 'SYSTEM') then
            backspace (u); exit
         end if
         if (line(1:5) /= 'BLOCK') cycle
         read (line(7:), *) kind, i, j, di, dj
         n = di*dj
         allocate (ref(n))
         read (u, *) (ref(m), m=1, n)
         k = 0
         ! a, c walk `ref` in the order it was WRITTEN -- libfint's own,
         ! LEXICOGRAPHIC -- so each is a lexicographic local offset;
         ! cart_perm_l2b maps it onto this build's own local offset before
         ! it becomes part of an index into smat/tmat/vmat, which come from
         ! THIS build's trc_1e and so are already build-order.
         do c = 0, dj - 1
            nu = off(j + 1) + cart_perm_l2b(bas(BAS_SLOTS*j + ANG_OF), c)
            do a = 0, di - 1
               mu = off(i + 1) + cart_perm_l2b(bas(BAS_SLOTS*i + ANG_OF), a)
               k = k + 1
               select case (trim(kind))
               case ('OVLP'); d = abs(smat(mu, nu) - ref(k)); m = 1
               case ('KIN');  d = abs(tmat(mu, nu) - ref(k)); m = 2
               case ('NUC');  d = abs(vmat(mu, nu) - ref(k)); m = 3
               case default;  d = 0.0_dp; m = 1
               end select
               tol = 1.0e-12_dp + 1.0e-11_dp*abs(ref(k))
               w1(m) = max(w1(m), d)
               if (d > tol) nbad = nbad + 1
            end do
         end do
         deallocate (ref)
      end do
98    continue
      print '(a,i0,a,i0,a,i0,a,i0,a)', '  1e system ', isys, ' (nao ', nao, &
         ', ', nkept, ' of ', nbas*(nbas + 1)/2, ' pairs kept) checked'
      deallocate (off, smat, tmat, vmat)
   end subroutine check_1e_system

   !
   ! The three-index tensor, with the ORBITAL and AUXILIARY bases both taken
   ! to be the whole reference basis. Physically odd -- an auxiliary basis is
   ! not normally the orbital one -- but it exercises every (la,lb|lc) the
   ! generator emits, which is what a correctness check is for.
   !
   subroutine check_3c_system(u, isys, atm, bas, env, natm, nbas, w3, nbad)
      integer, intent(in) :: u, isys, natm, nbas
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      real(dp), intent(inout) :: w3
      integer, intent(inout) :: nbad

      type(trc_basis_t) :: b, aux
      type(trc_pairlist_t) :: pl
      character(len=100) :: line
      integer :: i, j, kk, di, dj, dk, n, m, a, c, e, mu, nu, pp, idx, nao
      integer, allocatable :: off(:)
      real(dp), allocatable :: ref(:), tens(:, :, :)
      real(dp) :: d, tol

      allocate (off(nbas))
      call ao_offsets(bas, nbas, off, nao)

      call b%from_libcint(atm, natm, bas, nbas, env)
      call aux%from_libcint(atm, natm, bas, nbas, env)
      call pl%build(b, 1.0e-30_dp)
      allocate (tens(nao, nao, nao))
      !$acc enter data create(tens)
      call b%to_device()
      call aux%to_device()
      call pl%to_device()
      call trc_df_3c(b, pl, aux, tens)
      !$acc update self(tens)
      !$acc exit data delete(tens)
      call pl%release()
      call b%release()
      call aux%release()

      do
         read (u, '(a)', end=97) line
         if (line(1:3) == 'END' .or. line(1:6) == 'SYSTEM') then
            backspace (u); exit
         end if
         if (line(1:7) == 'BLOCK2C') then
            read (line(8:), *) i, j, di, dj
            n = di*dj
            allocate (ref(n)); read (u, *) (ref(m), m=1, n)
            deallocate (ref)
            cycle
         end if
         if (line(1:7) /= 'BLOCK3C') cycle
         read (line(8:), *) i, j, kk, di, dj, dk
         n = di*dj*dk
         allocate (ref(n))
         read (u, *) (ref(m), m=1, n)
         idx = 0
         ! e, c, a: same lexicographic-local-offset story as check_1e_system.
         do e = 0, dk - 1
            pp = off(kk + 1) + cart_perm_l2b(bas(BAS_SLOTS*kk + ANG_OF), e)
            do c = 0, dj - 1
               nu = off(j + 1) + cart_perm_l2b(bas(BAS_SLOTS*j + ANG_OF), c)
               do a = 0, di - 1
                  mu = off(i + 1) + cart_perm_l2b(bas(BAS_SLOTS*i + ANG_OF), a)
                  idx = idx + 1
                  d = abs(tens(mu, nu, pp) - ref(idx))
                  tol = 1.0e-12_dp + 1.0e-11_dp*abs(ref(idx))
                  w3 = max(w3, d)
                  if (d > tol) nbad = nbad + 1
               end do
            end do
         end do
         deallocate (ref)
      end do
97    continue
      print '(a,i0,a)', '  3c system ', isys, ' checked'
      deallocate (off, tens)
   end subroutine check_3c_system

end module check_api_suite

program check_api
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_api_suite, only: collect_api_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_api", collect_api_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)

   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_api
