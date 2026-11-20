module check_ecp_screen_suite
   !!
   !! The 3-centre screen of the ECP matrix (trc_ecp.F90's header). Two
   !! properties, on systems where the bound is tightest (s-s pairs on one
   !! heavy atom: Xe, SnH4) and where most tasks are screened (CH3I units
   !! tens of bohr apart):
   !!
   !!   * the bound is a bound with margin: on an unscreened run, every task's
   !!     largest element is at most MAX_RATIO times its bound (a ratio
   !!     above 1 would be a violation), and no task has a zero bound with a
   !!     nonzero element;
   !!   * the screened and the unscreened matrix agree to SCREEN_TOL, and the
   !!     screen really dropped tasks where it should (not vacuous). A bound
   !!     that dropped a contributing block changes a matrix element by far
   !!     more than SCREEN_TOL: the CH3I units sit 9 bohr apart, close enough
   !!     that the iodine channels of one reach the shells of the other.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_basis_types, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_json, only: trc_ecp_from_json
   use trc_ecp, only: trc_ecp_matrix_cart, trc_ecp_screen_stats_t
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_ecp_screen_tests

   real(dp), parameter :: SCREEN_TOL = 1.0e-15_dp
   real(dp), parameter :: MAX_RATIO = 0.25_dp
   !! Largest allowed (largest element)/bound. Below 1 is the property; 0.25
   !! is the margin the bound has here with 4x to spare (measured 0.053-0.056
   !! on these systems): dropping the 4*pi (ratio 0.67) or keeping only the
   !! coarsest level (0.37) still leaves a bound above every element on them,
   !! so only a margin check, not `<= 1`, catches those two.
   character(len=*), parameter :: DEF2_SVP = 'basis_sets/def2-svp.json'

contains

   subroutine collect_ecp_screen_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("3 CH3I, 9 and 60 bohr apart", test_three_units), &
                   new_unittest("Xe: s-s pairs on the ECP atom", test_xe), &
                   new_unittest("SnH4: one ECP atom, four partners", test_snh4), &
                   new_unittest("I2 bonded: nothing to drop", test_i2)]
   end subroutine collect_ecp_screen_tests

   subroutine ch3i_units(nunit, offsets, nat, elements, r)
      !! `nunit` CH3I (the geometry of check_ecp_matrix's ch3i case, bohr),
      !! unit u shifted along x by offsets(u).
      integer, intent(in) :: nunit
      real(dp), intent(in) :: offsets(nunit)
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: elements(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      real(dp), parameter :: UNIT_R(3, 5) = reshape([ &
                                                    0.0_dp, 0.0_dp, 0.0_dp, &
                                                    0.0_dp, 0.0_dp, 4.04_dp, &
                                                    1.9647_dp, 0.0_dp, -0.6195_dp, &
                                                    -0.98235_dp, 1.70148_dp, -0.6195_dp, &
                                                    -0.98235_dp, -1.70148_dp, -0.6195_dp], [3, 5])
      integer, parameter :: UNIT_Z(5) = [6, 53, 1, 1, 1]
      integer :: u, k

      nat = 5*nunit
      allocate (elements(nat), r(3, nat))
      do u = 1, nunit
         do k = 1, 5
            elements(5*(u - 1) + k) = UNIT_Z(k)
            r(:, 5*(u - 1) + k) = UNIT_R(:, k) + [offsets(u), 0.0_dp, 0.0_dp]
         end do
      end do
   end subroutine ch3i_units

   subroutine screened_vs_unscreened(error, label, basis_file, nat, elements, r, drops)
      !! `drops`: the screen must keep some tasks and drop some (.true.), or
      !! keep every one (.false.).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label, basis_file
      integer, intent(in) :: nat, elements(:)
      real(dp), intent(in) :: r(:, :)
      logical, intent(in) :: drops
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      type(trc_ecp_screen_stats_t) :: on, off
      real(dp), allocatable :: u_on(:, :), u_off(:, :)
      real(dp) :: worst

      call trc_basis_from_json(basis_file, nat, elements, r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      call check(error,.not. err%has_error(), label//": basis must build: "//err%get_message())
      if (allocated(error)) return
      call trc_ecp_from_json(basis_file, nat, elements, ecp, err)
      call check(error,.not. err%has_error(), label//": ECP must read: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         return
      end if

      call trc_ecp_matrix_cart(bas, ecp, u_on, err, screen=.true., stats=on)
      call check(error,.not. err%has_error(), label//": screened build: "//err%get_message())
      if (.not. allocated(error)) then
         call trc_ecp_matrix_cart(bas, ecp, u_off, err, screen=.false., stats=off)
         call check(error,.not. err%has_error(), label//": unscreened build: "//err%get_message())
      end if
      if (.not. allocated(error)) then
         call check(error, off%nkept == off%ntask .and. on%ntask == off%ntask, &
                    label//": the unscreened build must keep every task")
      end if
      if (.not. allocated(error)) then
         call check(error, maxval(abs(u_off)) > 1.0e-3_dp, label//": the matrix must not be trivially zero")
      end if
      if (.not. allocated(error)) then
         call check(error, off%nzero_bound == 0, label//": "//to_string(off%nzero_bound)// &
                    " tasks have a zero bound and a nonzero element")
      end if
      if (.not. allocated(error)) then
         call check(error, off%worst_ratio <= MAX_RATIO, label//": largest element/bound = "// &
                    to_string(off%worst_ratio, "es12.4")//", the bound must stay at least "// &
                    to_string(1.0_dp/MAX_RATIO, "f4.1")//"x above every element")
      end if
      if (.not. allocated(error)) then
         worst = maxval(abs(u_on - u_off))
         call check(error, worst <= SCREEN_TOL, label//": max|screened - unscreened| = "// &
                    to_string(worst, "es12.4")//" with "//to_string(on%nkept)//" of "//to_string(on%ntask)//" tasks kept")
      end if
      if (.not. allocated(error)) then
         if (drops) then
            call check(error, on%nkept > 0 .and. on%nkept < on%ntask, label// &
                       ": the screen must drop some tasks and keep some, kept "// &
                       to_string(on%nkept)//" of "//to_string(on%ntask))
         else
            call check(error, on%nkept == on%ntask, label//": the screen must keep every task, kept "// &
                       to_string(on%nkept)//" of "//to_string(on%ntask))
         end if
      end if
      call bas%release()
      call ecp%release()
   end subroutine screened_vs_unscreened

   subroutine test_three_units(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nat
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :)

      call ch3i_units(3, [0.0_dp, 9.0_dp, 60.0_dp], nat, elements, r)
      call screened_vs_unscreened(error, "3 CH3I", DEF2_SVP, nat, elements, r, .true.)
   end subroutine test_three_units

   subroutine test_xe(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:, :)

      allocate (r(3, 1))
      r(:, 1) = 0.0_dp
      call screened_vs_unscreened(error, "Xe", DEF2_SVP, 1, [54], r, .false.)
   end subroutine test_xe

   subroutine test_snh4(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:, :)

      allocate (r(3, 5))
      r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      r(:, 2) = [1.8652_dp, 1.8652_dp, 1.8652_dp]
      r(:, 3) = [1.8652_dp, -1.8652_dp, -1.8652_dp]
      r(:, 4) = [-1.8652_dp, 1.8652_dp, -1.8652_dp]
      r(:, 5) = [-1.8652_dp, -1.8652_dp, 1.8652_dp]
      call screened_vs_unscreened(error, "SnH4", DEF2_SVP, 5, [50, 1, 1, 1, 1], r, .false.)
   end subroutine test_snh4

   subroutine test_i2(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:, :)

      allocate (r(3, 2))
      r(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      r(:, 2) = [0.0_dp, 0.0_dp, 5.04_dp]
      call screened_vs_unscreened(error, "I2", DEF2_SVP, 2, [53, 53], r, .false.)
   end subroutine test_i2

end module check_ecp_screen_suite

program check_ecp_screen
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_ecp_screen_suite, only: collect_ecp_screen_tests
   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(collect_ecp_screen_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ecp_screen
