module check_fock_many_screen_suite
   !!
   !! `eri%fock_resident_many_lent` with density_screen = .false. (what a
   !! batched coupled-perturbed solve passes) against `eri%fock_resident` on
   !! each density alone, at a deliberately loose threshold so that a density
   !! keyed screen is visible: the two densities below are an O(1) one
   !! supported on the first half of the basis functions and a 1e-8 one over
   !! all of them. The screening bound is the maximum over the densities per
   !! shell pair, so with the screen on the quartets the tiny density needs
   !! away from the big one's support are dropped and its G is wrong by O(1)
   !! relative (a coupled-perturbed trial vector shrinks the same way as it
   !! converges); with it off every density sees the plain Schwarz bound, as
   !! when built alone. Both branches of the routine are covered: the shared
   !! bins pass, and the one-density-at-a-time loop it takes when MD J is
   !! requested (j_engine), which has to hand the flag on.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_eri, only: trc_eri_t, trc_fock_options_t, TRC_J_ENGINE_MDJ, TRC_J_ENGINE_AUTO
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_fock_many_screen

   real(dp), parameter :: LOOSE_THRESH = 1.0e-3_dp
   integer, parameter :: NDENS = 3

contains

   subroutine collect_fock_many_screen(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("shared bins pass, density screen off", test_bins), &
                  new_unittest("one density at a time (MD J requested), density screen off", test_md_loop) &
                  ]
   end subroutine collect_fock_many_screen

   subroutine densities(n, d)
      !! D1: O(1), nonzero only on the first half of the functions; D2: 1e-8
      !! everywhere; D3: 1e-3 on the second half.
      integer, intent(in) :: n
      real(dp), intent(out) :: d(n, n, NDENS)
      integer :: i, j
      d = 0.0_dp
      do j = 1, n
         do i = 1, n
            if (i <= n/2 .and. j <= n/2) d(i, j, 1) = cos(0.37_dp*real(i + 2*j, dp)) + cos(0.37_dp*real(j + 2*i, dp))
            d(i, j, 2) = 1.0e-8_dp*(sin(0.21_dp*real(i*j, dp)) + sin(0.21_dp*real(i + j, dp)) + 1.5_dp)
            if (i > n/2 .and. j > n/2) d(i, j, 3) = 1.0e-3_dp*(cos(0.5_dp*real(i, dp)*real(j, dp)) &
                                                              + cos(0.5_dp*real(i + j, dp)))
         end do
      end do
   end subroutine densities

   function rel_diff(a, b) result(w)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp) :: w
      w = maxval(abs(a - b))/max(maxval(abs(b)), tiny(1.0_dp))
   end function rel_diff

   subroutine compare(error, opts_base, label)
      !! many_lent with the screen off against fock_resident per density
      !! (screen off), and the screen-on build to show the threshold bites.
      type(error_type), allocatable, intent(out) :: error
      type(trc_fock_options_t), intent(in) :: opts_base
      character(len=*), intent(in) :: label
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(error_t) :: err
      type(trc_fock_options_t) :: o_off, o_on
      integer :: natm, dev, n, k
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), d(:, :, :), g_many(:, :, :), g_on(:, :, :), g_one(:, :, :), work(:)
      real(dp) :: w, w_on

      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/6-31gs.json', natm, zint, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call eri%build(bas, LOOSE_THRESH)
      n = bas%nao
      allocate (d(n, n, NDENS), g_many(n, n, NDENS), g_on(n, n, NDENS), g_one(n, n, NDENS))
      allocate (work(eri%fock_work_size(NDENS)))
      call densities(n, d)
      o_off = opts_base
      o_off%density_screen = .false.
      o_on = opts_base
      o_on%density_screen = .true.

      !$acc enter data copyin(d) create(g_many, g_on, g_one, work)
      call eri%fock_resident_many_lent(bas, NDENS, d, g_many, o_off, work)
      call eri%fock_resident_many_lent(bas, NDENS, d, g_on, o_on, work)
      do k = 1, NDENS
         call eri%fock_resident(bas, d(:, :, k), g_one(:, :, k), opts=o_off)
      end do
      !$acc update self(g_many, g_on, g_one)
      !$acc exit data delete(d, g_many, g_on, g_one, work)

      ! The screen-on build must visibly differ for the tiny density, or this
      ! threshold proves nothing.
      w_on = rel_diff(g_on(:, :, 2), g_one(:, :, 2))
      call check(error, w_on > 1.0e-3_dp, label//": the density screen should matter at this threshold, got "// &
                 to_string(w_on, "es10.2"))
      if (allocated(error)) then
         call eri%release(); call bas%release()
         return
      end if
      do k = 1, NDENS
         w = rel_diff(g_many(:, :, k), g_one(:, :, k))
         call check(error, w < 1.0e-10_dp, label//": density "//to_string(k)//" many vs alone, relative diff "// &
                    to_string(w, "es10.2")//" (screen on: "//to_string(w_on, "es10.2")//" for density 2)")
         if (allocated(error)) exit
      end do
      call eri%release()
      call bas%release()
   end subroutine compare

   subroutine test_bins(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_fock_options_t) :: o
      o = trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, j_engine=TRC_J_ENGINE_AUTO)
      call compare(error, o, "bins")
   end subroutine test_bins

   subroutine test_md_loop(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_fock_options_t) :: o
      o = trc_fock_options_t(j_scale=1.0_dp, k_scale=1.0_dp, j_engine=TRC_J_ENGINE_MDJ)
      call compare(error, o, "MD J loop")
   end subroutine test_md_loop

end module check_fock_many_screen_suite

program check_fock_many_screen
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_fock_many_screen_suite, only: collect_fock_many_screen
   implicit none
   integer :: stat

   stat = 0
   call run_testsuite(collect_fock_many_screen, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_fock_many_screen
