module check_exchange_screen_suite
   !! The density screens of an exchange-only bin launch (jfac = 0, which is
   !! what the bins run whenever the MD J engine forms J).
   !!
   !! Two properties, on a chain (40 waters in a line, STO-3G) with a
   !! density that decays with the distance between the two AO centres, as
   !! a real insulator's does:
   !!
   !! - the kernels' density screen counts only the blocks whose factor is
   !!   non-zero, so an exchange-only launch keeps strictly fewer quartets
   !!   than a Coulomb + exchange one, and its K still matches the one from
   !!   the full build (G_JK - G_J) to well inside the screening threshold;
   !! - the exchange tiles (trc_bin_tiles) prune the launch plan before
   !!   launch, and drop nothing the kernels would keep: the same quartets
   !!   survive, quartet for quartet, and the same matrix comes back.
   !!
   !! Both counts come from fock_resident's count_survivors, the kernels'
   !! own screen replayed on the plan the build launches.
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_bin_tiles, only: bin_tiles_set_enabled
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_exchange_screen

   real(dp), parameter :: THRESH = 1.0e-10_dp
   real(dp), parameter :: DECAY = 1.0_dp
   !! the model density's decay, per bohr between the two AO centres
   real(dp), parameter :: TOL_SAME = 1.0e-12_dp
   !! tiles on against off, relative to max|K|: summation order only
   real(dp), parameter :: TOL_SCREEN = 1.0e-8_dp
   !! K-only against G_JK - G_J, absolute: what the two screens may differ
   !! by is a sum of quartets each below THRESH

contains

   subroutine collect_exchange_screen(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("exchange-only screen keeps fewer quartets", test_exchange_screen), &
                  new_unittest("tiles prune the launch and drop nothing", test_tiles) &
                  ]
   end subroutine collect_exchange_screen

   subroutine setup(error, bas, d)
      !! The water chain in STO-3G and the decaying model density.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(out) :: bas
      real(dp), allocatable, intent(out) :: d(:, :)
      type(error_t) :: err
      integer :: nat, nao, si, sj, i, j, ni, nj
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: dist

      call read_xyz("water_chain40.xyz", nat, at_z, r)
      call trc_basis_from_json("basis_sets/sto-3g.json", nat, at_z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao))
      do sj = 1, bas%nshell
         nj = (bas%sh_l(sj) + 1)*(bas%sh_l(sj) + 2)/2
         do si = 1, bas%nshell
            ni = (bas%sh_l(si) + 1)*(bas%sh_l(si) + 2)/2
            dist = norm2(bas%sh_r(:, si) - bas%sh_r(:, sj))
            do j = bas%sh_ao(sj), bas%sh_ao(sj) + nj - 1
               do i = bas%sh_ao(si), bas%sh_ao(si) + ni - 1
                  d(i, j) = 0.4_dp*exp(-DECAY*dist)*(1.0_dp + 0.1_dp*cos(real(i + j, dp)))
               end do
            end do
         end do
      end do
      call bas%to_device()
   end subroutine setup

   subroutine build(eri, bas, nao, d, g, jfac, kfac, nkept, nwork)
      !! One resident build of G at these scalings (fock_resident's
      !! conventions), with its kept and launched quartet counts.
      type(trc_eri_t), intent(inout) :: eri
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: nao
      real(dp), intent(in) :: d(nao, nao)
      real(dp), intent(out) :: g(nao, nao)
      real(dp), intent(in) :: jfac, kfac
      integer(int64), intent(out) :: nkept, nwork
      !$acc enter data copyin(d) create(g)
      call eri%fock_resident(bas, d, g, opts=trc_fock_options_t(k_scale=kfac, j_scale=jfac, count_survivors=.true.))
      !$acc update self(g)
      !$acc exit data delete(d, g)
      nkept = eri%nkept
      nwork = eri%nwork
   end subroutine build

   subroutine test_exchange_screen(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      real(dp), allocatable :: d(:, :), gk(:, :), gjk(:, :), gj(:, :)
      integer(int64) :: kept_k, kept_jk, nw
      character(len=160) :: msg
      real(dp) :: worst

      call setup(error, bas, d)
      if (allocated(error)) return
      allocate (gk, gjk, gj, mold=d)
      call bin_tiles_set_enabled(.false.)
      call eri%build(bas, THRESH)
      call build(eri, bas, bas%nao, d, gj, 1.0_dp, 0.0_dp, kept_jk, nw)
      call build(eri, bas, bas%nao, d, gjk, 1.0_dp, 1.0_dp, kept_jk, nw)
      call build(eri, bas, bas%nao, d, gk, 0.0_dp, 1.0_dp, kept_k, nw)
      call bin_tiles_set_enabled(.true.)
      call eri%release()
      call bas%release()

      write (msg, "(a,i0,a,i0)") "kept quartets: exchange-only ", kept_k, ", Coulomb + exchange ", kept_jk
      print "(4x,a)", trim(msg)
      call check(error, kept_k < kept_jk, trim(msg))
      if (allocated(error)) return
      worst = maxval(abs(gk - (gjk - gj)))
      write (msg, "(a,es9.2,a,es9.2)") "max |G_K - (G_JK - G_J)| = ", worst, ", max|G_K| = ", maxval(abs(gk))
      print "(4x,a)", trim(msg)
      call check(error, worst <= TOL_SCREEN, trim(msg))
   end subroutine test_exchange_screen

   subroutine test_tiles(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      real(dp), allocatable :: d(:, :), gon(:, :), goff(:, :)
      integer(int64) :: kept_on, kept_off, work_on, work_off
      character(len=200) :: msg
      real(dp) :: worst

      call setup(error, bas, d)
      if (allocated(error)) return
      allocate (gon, goff, mold=d)
      call eri%build(bas, THRESH)
      call bin_tiles_set_enabled(.false.)
      call build(eri, bas, bas%nao, d, goff, 0.0_dp, 1.0_dp, kept_off, work_off)
      call bin_tiles_set_enabled(.true.)
      call build(eri, bas, bas%nao, d, gon, 0.0_dp, 1.0_dp, kept_on, work_on)
      call eri%release()
      call bas%release()

      write (msg, "(a,i0,a,i0,a,i0,a,i0)") "launched ", work_off, " -> ", work_on, &
         " quartets, kept ", kept_off, " -> ", kept_on
      print "(4x,a)", trim(msg)
      ! The plan must actually be pruned on this chain, or the rest proves
      ! nothing.
      call check(error, real(work_on, dp) < 0.9_dp*real(work_off, dp), "no pruning: "//trim(msg))
      if (allocated(error)) return
      call check(error, kept_on == kept_off, "a kept quartet was dropped: "//trim(msg))
      if (allocated(error)) return
      worst = maxval(abs(gon - goff))/maxval(abs(goff))
      write (msg, "(a,es9.2)") "max |G_tiles - G_whole|/max|G| = ", worst
      print "(4x,a)", trim(msg)
      call check(error, worst <= TOL_SAME, trim(msg))
   end subroutine test_tiles

end module check_exchange_screen_suite

program check_exchange_screen
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_exchange_screen_suite, only: collect_exchange_screen
   implicit none
   integer :: stat
   stat = 0
   ! Serial: the cases share the tiles switch and device-resident engines.
   call run_testsuite(collect_exchange_screen, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_exchange_screen
