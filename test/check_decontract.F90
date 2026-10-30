!
! The decontraction transform, against the one-electron integrals.
!
! `trc_decontract` splits a general contraction into single-primitive shells
! and carries the contraction as C. The claim is that C is exact, so for ANY
! one-electron operator over the split basis, C^T O_p C is the same matrix the
! contracted basis gives directly. Overlap and kinetic energy are two
! independent operators over the same functions, and `trc_1e` computes both
! without going anywhere near the four-centre path.
!
! This replaces a test that compared terco's blocked kernel against its scalar
! one. That comparison was between two paths that shared the same
! representation, so it stayed green through two real bugs in this transform:
! alpha derived before `basis_build` applied common_fac_sp, which left every s
! and p shell short by that factor, and a two-dimensional `do concurrent`
! whose `local` accumulator was not private. Both show up here immediately,
! because there is an external answer to be wrong against.
!
! Registered twice (ctest terco/check_decontract on water.xyz/cc-pvdz.json,
! terco/check_decontract_si on sih4.xyz/cc-pvdz.json -- twelve silicon s
! primitives under four contraction columns, the case the whole exercise
! exists for): one binary, argv picked by ctest, so argv is read in the
! program unit before the testsuite runs and the single resulting case is
! what test_decontract checks.
!
module check_decontract_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json
   use trc_decontract, only: decon_t, decontract_basis, decontract_fold_host
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_decontract_tests, g_xyzfile, g_basfile

   real(dp), parameter :: TOL = 1.0e-12_dp

   character(len=256) :: g_xyzfile = 'water.xyz'
   character(len=256) :: g_basfile = 'basis_sets/cc-pvdz.json'
      !! Set from argv by the program unit before run_testsuite is called;
      !! test_decontract reads these, never argv itself.

contains

   subroutine collect_decontract_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("C^T S_p C = S and C^T T_p C = T", test_decontract) &
                  ]
   end subroutine collect_decontract_tests

   subroutine test_decontract(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b, pb
      type(trc_pairlist_t) :: p, pp
      type(decon_t) :: m
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp), allocatable :: s(:, :), t(:, :), v(:, :), sp(:, :), tp(:, :), vp(:, :), f(:, :)
      integer :: natm, i, j
      real(dp) :: ws, wt, scl

      call read_xyz(trim(g_xyzfile), natm, zint, at_r)
      call trc_basis_from_json(trim(g_basfile), natm, zint, at_r, b, err)
      call check(error,.not. err%has_error(), "check_decontract: basis: "//err%get_message())
      if (allocated(error)) return

      call decontract_basis(b, pb, m)
      print '(a,l2,a,i0,a,i0,a,i0,a,i0)', '  active', m%active, '   nao ', b%nao, ' -> ', &
         pb%nao, '   shells ', b%nshell, ' -> ', pb%nshell
      if (.not. m%active) then
         ! A segmented basis is left alone, and that is the correct answer.
         call check(error, .true., "decontraction left the basis alone (nothing to split)")
         return
      end if

      call b%to_device(); call pb%to_device()
      call p%build(b, 1.0e-14_dp); call pp%build(pb, 1.0e-14_dp)
      allocate (s(b%nao, b%nao), t(b%nao, b%nao), v(b%nao, b%nao), f(b%nao, b%nao))
      allocate (sp(pb%nao, pb%nao), tp(pb%nao, pb%nao), vp(pb%nao, pb%nao))
      call trc_1e(b, p, s, t, v)
      call trc_1e(pb, pp, sp, tp, vp)

      call decontract_fold_host(m, sp, f)
      ws = 0.0_dp; scl = 0.0_dp
      do j = 1, b%nao
         do i = 1, b%nao
            ws = max(ws, abs(f(i, j) - s(i, j))); scl = max(scl, abs(s(i, j)))
         end do
      end do
      call decontract_fold_host(m, tp, f)
      wt = 0.0_dp
      do j = 1, b%nao
         do i = 1, b%nao
            wt = max(wt, abs(f(i, j) - t(i, j)))
         end do
      end do

      print '(a,es12.4,a,es12.4)', '  worst |C^T S_p C - S| ', ws, '   scale ', scl
      print '(a,es12.4)', '  worst |C^T T_p C - T| ', wt

      call check(error, ws, 0.0_dp, thr=TOL, message="worst |C^T S_p C - S| exceeds TOL")
      if (allocated(error)) return
      call check(error, wt, 0.0_dp, thr=TOL, message="worst |C^T T_p C - T| exceeds TOL")
   end subroutine test_decontract

end module check_decontract_suite

program check_decontract
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_decontract_suite, only: collect_decontract_tests, g_xyzfile, g_basfile
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   ! argv is read here, before the testsuite runs: ctest hands this binary
   ! water.xyz/cc-pvdz.json (terco/check_decontract) or
   ! sih4.xyz/cc-pvdz.json (terco/check_decontract_si), and test_decontract
   ! only ever sees the pair chosen for this one invocation.
   if (command_argument_count() >= 1) call get_command_argument(1, g_xyzfile)
   if (command_argument_count() >= 2) call get_command_argument(2, g_basfile)

   stat = 0
   testsuites = [new_testsuite("check_decontract", collect_decontract_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_decontract
