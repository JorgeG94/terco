!
! Rys against HGP, class by class, on the SAME density and screening.
!
! With TRC_RYS_CLASSES=all, every class the build generated a kernel for
! routes through Rys instead of
! HGP (trc_rys_select%rys_route); everything else -- including classes this
! build never generated a Rys kernel for -- still goes through HGP. Building
! the SAME Fock matrix twice, once per setting, and comparing is a stronger
! check than trusting either kernel alone: an integral that is wrong the
! same way in both schemes would pass a check against itself but not this
! one, because HGP and Rys compute it by unrelated recurrences.
!
! Two systems:
!   - water/6-31G*: light, segmented, exercises the d-shell classes without
!     general contraction.
!   - water/cc-pVDZ: the oxygen s and p shells are generally contracted, so
!     this exercises trc_decontract's split alongside the same comparison.
!
! TRC_RYS_CLASSES is read ONCE PER FOCK BUILD (trc_rys_select), so toggling
! it between the two eri%fock calls below -- rather than, say, building two
! separate eri_t objects -- is exactly the intended use: same bins, same
! density, same screening, only the env var differs.
!
program check_rys_vs_hgp
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   interface
      function c_setenv(name, value, overwrite) bind(c, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
      end function c_setenv
   end interface

   real(dp), parameter :: ATOL = 1.0e-12_dp
   integer :: nat, i, j
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: at_r(:, :)
   real(dp) :: worst
   integer :: nbad
   character(len=256) :: basis_json

   worst = 0.0_dp; nbad = 0

   ! One basis per process, not two trc_eri_t builds in one -- the OpenACC
   ! present-table bookkeeping underneath trc_api's device-resident types
   ! is not validated for a second build+teardown in the same process (a
   ! property of that bookkeeping, not of either integral scheme), so two
   ! ctest invocations cover the two systems (water/6-31G*, water/cc-pVDZ).
   if (command_argument_count() < 1) then
      print '(a)', 'usage: check_rys_vs_hgp <basis.json>'; stop 1
   end if
   call get_command_argument(1, basis_json)
   call run_system(trim(basis_json), trim(basis_json))

   print '(a,es10.2)', '  worst |G_rys - G_hgp| : ', worst
   print '(a,i0)',     '  outside tol          : ', nbad
   if (nbad > 0) then
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if
   print '(a)', '  RESULT: PASS'

contains

   subroutine set_rys(v)
      character(*), intent(in) :: v
      integer(c_int) :: r
      r = c_setenv("TRC_RYS_CLASSES"//c_null_char, trim(v)//c_null_char, 1_c_int)
      if (r /= 0) then
         print '(a)', '  setenv(TRC_RYS_CLASSES) failed'; stop 1
      end if
   end subroutine set_rys

   subroutine run_system(basis_json, tag)
      character(*), intent(in) :: basis_json, tag
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      real(dp), allocatable :: dmat(:, :), g_hgp(:, :), g_rys(:, :)
      real(dp) :: thresh, local_worst
      integer :: nao, ii, jj, nbad0

      call read_xyz('water.xyz', nat, at_z, at_r)
      call trc_basis_from_json(basis_json, nat, at_z, at_r, bas, err)
      if (err%has_error()) then
         print '(a)', 'check_rys_vs_hgp: '//err%get_message(); stop 1
      end if

      thresh = 1.0e-13_dp
      nao = bas%nao
      call bas%to_device()
      call eri%build(bas, thresh)

      allocate (dmat(nao, nao), g_hgp(nao, nao), g_rys(nao, nao))
      ! Deterministic, dense, symmetric -- same density for both builds, as
      ! the gate requires, and not sparse enough to hide an indexing slip.
      do j = 1, nao
         do i = 1, nao
            dmat(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j, dp))
         end do
      end do
      dmat = 0.5_dp*(dmat + transpose(dmat))

      call set_rys('none')
      call eri%fock(bas, dmat, g_hgp, density_screen=.false.)
      call set_rys('all')
      call eri%fock(bas, dmat, g_rys, density_screen=.false.)
      call set_rys('none')

      local_worst = 0.0_dp; nbad0 = 0
      do jj = 1, nao
         do ii = 1, nao
            local_worst = max(local_worst, abs(g_hgp(ii, jj) - g_rys(ii, jj)))
            if (abs(g_hgp(ii, jj) - g_rys(ii, jj)) > ATOL) nbad0 = nbad0 + 1
         end do
      end do
      print '(a,a,a,i0,a,es10.2,a,i0)', '  ', tag, ': nao ', nao, &
         '  worst diff ', local_worst, '  outside tol ', nbad0
      worst = max(worst, local_worst)
      nbad = nbad + nbad0

      deallocate (dmat, g_hgp, g_rys, at_z, at_r)
   end subroutine run_system

end program check_rys_vs_hgp
