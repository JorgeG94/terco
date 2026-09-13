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
program check_decontract
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_decontract, only: decon_t, decontract_basis, decon_fold_host
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   type(trc_basis_t) :: b, pb
   type(trc_pairlist_t) :: p, pp
   type(decon_t) :: m
   type(error_t) :: err
   character(len=256) :: xyzfile, basfile
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp), allocatable :: s(:, :), t(:, :), v(:, :), sp(:, :), tp(:, :), vp(:, :), f(:, :)
   integer :: natm, i, j, dev
   real(dp) :: ws, wt, scl
   real(dp), parameter :: TOL = 1.0e-12_dp

   dev = trc_bind_device(0)
   xyzfile = 'water.xyz'; basfile = 'basis_sets/cc-pvdz.json'
   if (command_argument_count() >= 1) call get_command_argument(1, xyzfile)
   if (command_argument_count() >= 2) call get_command_argument(2, basfile)
   call read_xyz(trim(xyzfile), natm, zint, at_r)
   call trc_basis_from_json(trim(basfile), natm, zint, at_r, b, err)
   if (err%has_error()) error stop 'check_decontract: basis'

   call decontract_basis(b, pb, m)
   print '(a,l2,a,i0,a,i0,a,i0,a,i0)', '  active', m%active, '   nao ', b%nao, ' -> ', &
      pb%nao, '   shells ', b%nshell, ' -> ', pb%nshell
   if (.not. m%active) then
      ! A segmented basis is left alone, and that is the correct answer.
      print '(a)', '  RESULT: PASS (nothing to split)'
      stop 0
   end if

   call b%to_device(); call pb%to_device()
   call p%build(b, 1.0e-14_dp); call pp%build(pb, 1.0e-14_dp)
   allocate (s(b%nao, b%nao), t(b%nao, b%nao), v(b%nao, b%nao), f(b%nao, b%nao))
   allocate (sp(pb%nao, pb%nao), tp(pb%nao, pb%nao), vp(pb%nao, pb%nao))
   call trc_1e(b, p, s, t, v)
   call trc_1e(pb, pp, sp, tp, vp)

   call decon_fold_host(m, sp, f)
   ws = 0.0_dp; scl = 0.0_dp
   do j = 1, b%nao
      do i = 1, b%nao
         ws = max(ws, abs(f(i, j) - s(i, j))); scl = max(scl, abs(s(i, j)))
      end do
   end do
   call decon_fold_host(m, tp, f)
   wt = 0.0_dp
   do j = 1, b%nao
      do i = 1, b%nao
         wt = max(wt, abs(f(i, j) - t(i, j)))
      end do
   end do

   print '(a,es12.4,a,es12.4)', '  worst |C^T S_p C - S| ', ws, '   scale ', scl
   print '(a,es12.4)', '  worst |C^T T_p C - T| ', wt
   if (ws > TOL .or. wt > TOL) then
      print '(a)', '  RESULT: FAIL'
      error stop 1
   end if
   print '(a)', '  RESULT: PASS'
end program check_decontract
