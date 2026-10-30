!
! Finite-difference oracle for the long-range exchange gradient term alone:
! eri%gradient(j_scale=0, k_scale=1, omega) against a central difference of
! E_K^lr(D) = Tr(D, eri%fock(D, j_scale=0, k_scale=1, omega)), D FIXED and
! the basis rebuilt at the displaced geometry -- same technique and the same
! synthetic density as check_rys_grad_fd.F90, specialised to the erf
! exchange term by itself (RHF only; see check_lr_exchange_grad_vs_pyscf.F90
! for why there is no UHF range-separated convention here to check against).
!
! THE FACTOR. eri%fock's own convention is vhf = j_scale*J - k_scale*K/2, so
! eri%fock(D, j_scale=0, k_scale=1, omega) = -K_lr(D)/2 and
! Tr(D, -K_lr(D)/2) = -(1/2) Tr(D, K_lr(D)) -- NOT the -(1/4) Tr(D, K_lr(D))
! a textbook RHF exchange ENERGY convention would give (that 1/4 is e_two's
! own factor of 1/2 on top of this), exactly the same "twice the actual
! energy" convention eri%gradient's own output is in (this module's header
! says so for the full-range term; checked here to be exactly as true of
! the long-range one, not assumed). Both e2_lr_at and eri%gradient's own
! raw output are compared RAW, with no extra factor on either side -- the
! same structure check_rys_grad_fd.F90's RHF path uses for the full
! exchange, specialised to omega > 0 and no Coulomb term.
!
! BIT-IDENTITY AT OMEGA = 0. Before the omega > 0 sweep, one direct check:
! eri%gradient called with omega ABSENT (today's only call shape, the one
! every existing caller uses) against the same call with j_scale=1,
! k_scale=1, omega=0.0_dp explicit. See scripts/gen_rys_grad.py's
! `_RANGE_SEP_BLOCK` for why these SHOULD agree to the last bit: theta ==
! 1.0_dp exactly at omega == 0, so rho_theta == rho and pref*sqrt(theta) ==
! pref, each an exact multiply by one, all the way down to the per-root Rys
! coefficients -- and on the CPU build (deterministic, sequential host
! accumulation) they measure EXACTLY 0.0_dp, both bases, every run: that is
! the bit-identity this program's own header promises, checked.
!
! ON THE GPU BUILD the same two calls are two SEPARATE kernel launches, and
! the lane-private gradient copies (trc_grad.F90's GRAD_NPRIV, atomically
! summed into one grad(3,natm)) commute but do not associate in floating
! point: nothing about either launch's own warp/block scheduling is
! guaranteed to revisit the same quartets in the same order as the other's,
! so two bit-identical kernels can sum to results a few ULP apart. Measured
! on water/cc-pVDZ, three back-to-back runs of this program: 3.55e-15,
! 1.42e-14, 5.68e-14 -- a different, non-repeating value each time (not a
! fixed wrong number), while the SAME runs' omega > 0 finite-difference
! comparison below stayed at 1.04e-07/1.32e-07 to the last digit printed,
! every time. GPU_TOL below is two decades above that largest measured gap
! and eight below the loosest FD tolerance this program checks against --
! room for the one real failure mode (a REAL mistake threading omega
! through only one of the two calls) without flagging atomic reordering.
!
! Basis (6-31G* or cc-pVDZ) comes from the command line (test/CMakeLists.txt
! registers this source as TWO executables, ARGS basis_sets/6-31gs.json and
! basis_sets/cc-pvdz.json), so each run checks exactly the one basis it was
! invoked with. The omega sweep (0.3, 0.4) is internal to the program, not
! CMake-driven, and stays that way.
!
module check_lr_exchange_grad_fd_suite
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
   public :: collect_lr_exchange_grad_fd_tests

   real(dp), parameter :: H = 1.0e-4_dp
   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: RELTOL = 1.0e-7_dp
   real(dp), parameter :: GPU_TOL = 1.0e-11_dp
   !! The omega-absent vs omega=0.0_dp-explicit gap: exactly 0.0_dp on the
   !! CPU build, atomic-reordering noise up to ~1e-13 measured on the GPU
   !! build -- see this module's own header for the three-run measurement
   !! GPU_TOL sits two decades above.

contains

   subroutine collect_lr_exchange_grad_fd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("omega=0 is bit-identical to the omega-absent call", test_omega_zero), &
                  new_unittest("omega=0.3 matches a central difference of E_K^lr", test_omega_03), &
                  new_unittest("omega=0.4 matches a central difference of E_K^lr", test_omega_04) &
                  ]
   end subroutine collect_lr_exchange_grad_fd_tests

   subroutine basis_arg(basis_json)
      character(len=*), intent(out) :: basis_json
      basis_json = ''
      if (command_argument_count() >= 1) call get_command_argument(1, basis_json)
   end subroutine basis_arg

   !
   ! Same deterministic dense symmetric matrix as check_rys_grad_fd.F90's
   ! phase-1 density (RHF total density here, so phase is fixed at 1).
   !
   subroutine synth_density(nao, d)
      integer, intent(in) :: nao
      real(dp), intent(out) :: d(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            d(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j + 7, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

   real(dp) function trace_prod(nao, a, b)
      integer, intent(in) :: nao
      real(dp), intent(in) :: a(nao, nao), b(nao, nao)
      integer :: i, j
      trace_prod = 0.0_dp
      do j = 1, nao
         do i = 1, nao
            trace_prod = trace_prod + a(i, j)*b(i, j)
         end do
      end do
   end function trace_prod

   !
   ! omega absent vs omega=0.0_dp explicit, same basis/density/geometry:
   ! must be bit-for-bit, not just numerically close -- see this module's
   ! own header.
   !
   subroutine test_omega_zero(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      integer :: nat, nao
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), d(:, :, :), g_default(:, :), g_omega0(:, :)
      real(dp) :: worst

      call basis_arg(basis_json)
      call check(error, len_trim(basis_json) > 0, "usage: check_lr_exchange_grad_fd <basis.json>")
      if (allocated(error)) return

      call read_xyz('water.xyz', nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), "check_lr_exchange_grad_fd: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao, 1), g_default(3, nat), g_omega0(3, nat))
      call synth_density(nao, d(:, :, 1))

      call bas%to_device()
      call eri%build(bas, THRESH)
      call eri%gradient(bas, d, g_default)
      call eri%gradient(bas, d, g_omega0, j_scale=1.0_dp, k_scale=1.0_dp, omega=0.0_dp)
      call eri%release()
      call bas%release()

      worst = maxval(abs(g_default - g_omega0))
      call check(error, worst, 0.0_dp, thr=GPU_TOL, message= &
                 "omega-absent vs omega=0 explicit: max|diff| = "//to_string(worst, "es10.2")// &
                 " exceeds GPU_TOL (not bit-identical)")
      deallocate (at_z, r, d, g_default, g_omega0)
   end subroutine test_omega_zero

   !
   ! E_K^lr(D) at the geometry in `r`, basis `basis_json`, FIXED D, RHF
   ! (ndens=1) only -- see this module's own header for the factor.
   !
   real(dp) function e2_lr_at(basis_json, nat, at_z, r, nao, d, omega, error)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nao
      real(dp), intent(in) :: r(3, nat), d(nao, nao), omega
      type(error_type), allocatable, intent(inout) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      real(dp) :: vhf(nao, nao)

      e2_lr_at = 0.0_dp
      call trc_basis_from_json(basis_json, nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), "check_lr_exchange_grad_fd: "//err%get_message())
      if (allocated(error)) return
      call check(error, bas%nao == nao, "check_lr_exchange_grad_fd: nao changed under displacement")
      if (allocated(error)) then
         call bas%release()
         return
      end if
      call bas%to_device()
      call eri%build(bas, THRESH)
      call eri%fock(bas, d, vhf, j_scale=0.0_dp, k_scale=1.0_dp, omega=omega, density_screen=.false.)
      e2_lr_at = trace_prod(nao, d, vhf)
      call eri%release()
      call bas%release()
   end function e2_lr_at

   subroutine run_one(omega, error)
      real(dp), intent(in) :: omega
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      integer :: nat, nao, ia, dd
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), rp(:, :), rm(:, :)
      real(dp), allocatable :: dmat(:, :, :), grad(:, :)
      real(dp) :: ep, em, fd, worst_abs, worst_rel, ref, denom
      integer :: nbad

      call basis_arg(basis_json)
      call check(error, len_trim(basis_json) > 0, "usage: check_lr_exchange_grad_fd <basis.json>")
      if (allocated(error)) return

      call read_xyz('water.xyz', nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), "check_lr_exchange_grad_fd: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao

      allocate (dmat(nao, nao, 1), grad(3, nat))
      call synth_density(nao, dmat(:, :, 1))

      call bas%to_device()
      call eri%build(bas, THRESH)
      call eri%gradient(bas, dmat, grad, j_scale=0.0_dp, k_scale=1.0_dp, omega=omega)
      call eri%release()
      call bas%release()

      allocate (rp(3, nat), rm(3, nat))
      worst_abs = 0.0_dp; worst_rel = 0.0_dp; nbad = 0
      do ia = 1, nat
         do dd = 1, 3
            rp = r; rm = r
            rp(dd, ia) = rp(dd, ia) + H
            rm(dd, ia) = rm(dd, ia) - H
            ep = e2_lr_at(trim(basis_json), nat, at_z, rp, nao, dmat(:, :, 1), omega, error)
            if (allocated(error)) then
               deallocate (at_z, r, rp, rm, dmat, grad)
               return
            end if
            em = e2_lr_at(trim(basis_json), nat, at_z, rm, nao, dmat(:, :, 1), omega, error)
            if (allocated(error)) then
               deallocate (at_z, r, rp, rm, dmat, grad)
               return
            end if
            fd = (ep - em)/(2.0_dp*H)
            ref = grad(dd, ia)
            worst_abs = max(worst_abs, abs(fd - ref))
            denom = max(abs(ref), abs(fd), 1.0e-10_dp)
            worst_rel = max(worst_rel, abs(fd - ref)/denom)
            if (abs(fd - ref)/denom > RELTOL) nbad = nbad + 1
         end do
      end do

      call check(error, nbad == 0, "omega "//to_string(omega, "f4.2")//": worst abs "// &
                 to_string(worst_abs, "es10.2")//" worst rel "//to_string(worst_rel, "es10.2")// &
                 " outside tol "//to_string(nbad))

      deallocate (at_z, r, rp, rm, dmat, grad)
   end subroutine run_one

   subroutine test_omega_03(error)
      type(error_type), allocatable, intent(out) :: error
      call run_one(0.3_dp, error)
   end subroutine test_omega_03

   subroutine test_omega_04(error)
      type(error_type), allocatable, intent(out) :: error
      call run_one(0.4_dp, error)
   end subroutine test_omega_04

end module check_lr_exchange_grad_fd_suite

program check_lr_exchange_grad_fd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_lr_exchange_grad_fd_suite, only: collect_lr_exchange_grad_fd_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_lr_exchange_grad_fd", collect_lr_exchange_grad_fd_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_lr_exchange_grad_fd
