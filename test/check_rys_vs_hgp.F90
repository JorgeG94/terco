!
! Rys against HGP, class by class, on the SAME density and screening.
!
! With TRC_RYS_CLASSES=all, every class the build generated a kernel for
! routes through Rys instead of
! HGP (trc_rys_select%rys_select_route); everything else -- including classes this
! build never generated a Rys kernel for -- still goes through HGP. Building
! the SAME Fock matrix twice, once per setting, and comparing is a stronger
! check than trusting either kernel alone: an integral that is wrong the
! same way in both schemes would pass a check against itself but not this
! one, because HGP and Rys compute it by unrelated recurrences.
!
! Two systems, each its own ctest invocation (command-line argument, see
! the module header below for why they cannot share a process):
!   - water/6-31G*: light, segmented, exercises the d-shell classes without
!     general contraction.
!   - water/cc-pVDZ: the oxygen s and p shells are generally contracted, so
!     this exercises trc_decontract's split alongside the same comparison.
!
! The default routing (TRC_RYS_CLASSES unset: HGP, except the d classes
! gen_rys.py's FOCK_RYS_DEFAULT sends to the chunked Rys kernels) is a
! third build, compared against HGP the same way.
!
! TRC_RYS_CLASSES is read ONCE PER FOCK BUILD (trc_rys_select), so toggling
! it between the eri%fock calls below -- rather than, say, building two
! separate eri_t objects -- is exactly the intended use: same bins, same
! density, same screening, only the env var differs.
!
module check_rys_vs_hgp_suite
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_rys_vs_hgp, basis_json_arg

   interface
      function c_setenv(name, value, overwrite) bind(c, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
      end function c_setenv
   end interface

   real(dp), parameter :: ATOL = 1.0e-12_dp

   character(len=256) :: basis_json_arg = ''
   !! Set by the program from argv(1) before run_testsuite: ONE basis per
   !! process, not two trc_eri_t builds in one -- the OpenACC present-table
   !! bookkeeping underneath trc_api's device-resident types is not
   !! validated for a second build+teardown in the same process, so two
   !! ctest invocations (different ARGS) cover the two systems.

contains

   subroutine collect_rys_vs_hgp(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Rys vs HGP on water/<basis>", test_vs_hgp) &
                  ]
   end subroutine collect_rys_vs_hgp

   subroutine set_rys(v)
      character(*), intent(in) :: v
      integer(c_int) :: r
      r = c_setenv("TRC_RYS_CLASSES"//c_null_char, trim(v)//c_null_char, 1_c_int)
      if (r /= 0) error stop 'setenv(TRC_RYS_CLASSES) failed'
   end subroutine set_rys

   subroutine test_vs_hgp(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      real(dp), allocatable :: dmat(:, :), g_hgp(:, :), g_rys(:, :), g_def(:, :)
      real(dp) :: thresh, worst_rys, worst_def
      integer :: nat, nao, i, j
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: at_r(:, :)

      call read_xyz('water.xyz', nat, at_z, at_r)
      call trc_basis_from_json(trim(basis_json_arg), nat, at_z, at_r, bas, err)
      call check(error, .not. err%has_error(), "basis build: "//err%get_message())
      if (allocated(error)) return

      thresh = 1.0e-13_dp
      nao = bas%nao
      call bas%to_device()
      call eri%build(bas, thresh)

      allocate (dmat(nao, nao), g_hgp(nao, nao), g_rys(nao, nao), g_def(nao, nao))
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
      call set_rys('default')
      call eri%fock(bas, dmat, g_def, density_screen=.false.)
      call set_rys('none')

      worst_rys = maxval(abs(g_hgp - g_rys))
      worst_def = maxval(abs(g_hgp - g_def))
      print '(a,i0,a,es10.2,a,es10.2)', '  nao ', nao, &
         '  worst |G_rys - G_hgp| ', worst_rys, '  worst |G_def - G_hgp| ', worst_def

      call eri%release(); call bas%release()
      deallocate (dmat, g_hgp, g_rys, g_def, at_z, at_r)

      call check(error, worst_rys <= ATOL, &
                 "TRC_RYS_CLASSES=all: worst |G_rys - G_hgp| = "//to_string(worst_rys, "es10.2")// &
                 " exceeds "//to_string(ATOL, "es10.2"))
      if (allocated(error)) return
      call check(error, worst_def <= ATOL, &
                 "TRC_RYS_CLASSES=default: worst |G_def - G_hgp| = "//to_string(worst_def, "es10.2")// &
                 " exceeds "//to_string(ATOL, "es10.2"))
   end subroutine test_vs_hgp

end module check_rys_vs_hgp_suite

program check_rys_vs_hgp
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rys_vs_hgp_suite, only: collect_rys_vs_hgp, basis_json_arg
   implicit none
   integer :: stat

   if (command_argument_count() < 1) then
      write (error_unit, '(a)') 'usage: check_rys_vs_hgp <basis.json>'
      error stop 1
   end if
   call get_command_argument(1, basis_json_arg)

   stat = 0
   call run_testsuite(collect_rys_vs_hgp, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rys_vs_hgp
