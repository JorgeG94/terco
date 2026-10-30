!
! Finite-difference oracle for eri%gradient: the two-electron part of the
! nuclear gradient against a central difference of E2 = Tr(D * vhf(D)),
! with D FIXED and the basis rebuilt at the displaced geometry so the basis
! functions move with their atom.
!
! D is a deterministic, dense, symmetric synthetic density -- same
! construction check_rys_vs_hgp.F90 uses -- not a converged SCF density.
! This is a kernel check (is d(ab|cd)/dR right?), not an SCF check, and a
! fixed synthetic D exercises every six-block path exactly as a converged
! one would while costing nothing to construct and needing no separate SCF
! driver here.
!
! RHF: E2 = Tr(D * eri%fock(D)), one Fock build per displacement.
! UHF (Da /= Db): split by linearity --
!     E2 = Tr(Dt*J(Dt)) - Tr(Da*Ka(Da)) - Tr(Db*Kb(Db))
!   with NO 1/2 anywhere (checked against the Da = Db reduction to RHF --
!   eri%gradient(Da, Da) must equal eri%gradient(2*Da) exactly, since UHF
!   with identical spin densities IS the RHF solution -- not taken from a
!   textbook formula's placement of the 1/2 on faith). eri%fock's k_scale/
!   j_scale isolate J and K:
!     J(Dt)  = eri%fock(Dt, k_scale=0)
!     -Ka(Da) = eri%fock(Da, j_scale=0, k_scale=2)   (k_scale=2 cancels
!                                                      eri%fock's own J-K/2)
!     -Kb(Db) = eri%fock(Db, j_scale=0, k_scale=2)
!   so E2 = Tr(Dt, J(Dt)) + Tr(Da, -Ka(Da)) + Tr(Db, -Kb(Db)).
!
! Translational invariance: sum over atoms of eri%gradient, every component,
! against zero -- a property of the ANALYTIC gradient alone, not of the
! finite difference.
!
! The basis (6-31G* or cc-pVDZ) comes from the command line (test/CMakeLists.txt
! registers this source as TWO executables, ARGS basis_sets/6-31gs.json and
! basis_sets/cc-pvdz.json, each also setting TRC_RYS_CLASSES=none so the
! HGP path is what gets checked), so each run checks exactly the one basis
! it was invoked with.
!
module check_rys_grad_fd_suite
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
   public :: collect_rys_grad_fd_tests

   real(dp), parameter :: H = 1.0e-4_dp
   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: RELTOL = 1.0e-7_dp
   real(dp), parameter :: INVTOL = 1.0e-10_dp

contains

   subroutine collect_rys_grad_fd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF gradient matches FD and is translation-invariant", test_rhf), &
                  new_unittest("UHF gradient matches FD and is translation-invariant", test_uhf) &
                  ]
   end subroutine collect_rys_grad_fd_tests

   subroutine basis_arg(basis_json, error)
      character(len=*), intent(out) :: basis_json
      type(error_type), allocatable, intent(out) :: error
      basis_json = ''
      if (command_argument_count() >= 1) call get_command_argument(1, basis_json)
      call check(error, len_trim(basis_json) > 0, "usage: check_rys_grad_fd <basis.json>")
   end subroutine basis_arg

   !
   ! A deterministic dense symmetric matrix, same closed form as
   ! check_rys_vs_hgp -- `phase` separates Da from Db so UHF's Da /= Db.
   !
   subroutine synth_density(nao, phase, d)
      integer, intent(in) :: nao, phase
      real(dp), intent(out) :: d(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            d(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j + 7*phase, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

   !
   ! Tr(a^T b) over the full square matrix -- both symmetric here, so this
   ! is Tr(a b).
   !
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
   ! E2 at the geometry in `r`, basis `basis_json`, for a FIXED set of
   ! densities. ndens = size(dmat, 3): 1 is RHF (dmat(:,:,1) = Dtot), 2 is
   ! UHF (dmat(:,:,1:2) = Da, Db).
   !
   real(dp) function e2_at(basis_json, nat, at_z, r, nao, ndens, dmat, error)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nao, ndens
      real(dp), intent(in) :: r(3, nat), dmat(nao, nao, ndens)
      type(error_type), allocatable, intent(inout) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      real(dp) :: vhf(nao, nao), dtot(nao, nao)

      e2_at = 0.0_dp
      call trc_basis_from_json(basis_json, nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), "check_rys_grad_fd: "//err%get_message())
      if (allocated(error)) return
      call check(error, bas%nao == nao, "check_rys_grad_fd: nao changed under displacement")
      if (allocated(error)) then
         call bas%release()
         return
      end if
      call bas%to_device()
      call eri%build(bas, THRESH)

      if (ndens == 1) then
         call eri%fock(bas, dmat(:, :, 1), vhf, density_screen=.false.)
         e2_at = trace_prod(nao, dmat(:, :, 1), vhf)
      else
         ! UHF: E2 = Tr(Dt, J(Dt)) - Tr(Da, Ka(Da)) - Tr(Db, Kb(Db)). No 1/2
         ! anywhere -- checked against the gradient this is meant to match
         ! (eri%gradient(Da, Da) == eri%gradient(2Da) exactly, the Da = Db
         ! reduction to RHF), not against a textbook formula taken on faith.
         ! eri%fock returns J - K/2, so k_scale=2 cancels that /2 and leaves
         ! -K(D) directly.
         dtot = dmat(:, :, 1) + dmat(:, :, 2)
         call eri%fock(bas, dtot, vhf, k_scale=0.0_dp, density_screen=.false.)
         e2_at = trace_prod(nao, dtot, vhf)
         call eri%fock(bas, dmat(:, :, 1), vhf, j_scale=0.0_dp, k_scale=2.0_dp, &
                      density_screen=.false.)
         e2_at = e2_at + trace_prod(nao, dmat(:, :, 1), vhf)
         call eri%fock(bas, dmat(:, :, 2), vhf, j_scale=0.0_dp, k_scale=2.0_dp, &
                      density_screen=.false.)
         e2_at = e2_at + trace_prod(nao, dmat(:, :, 2), vhf)
      end if

      call eri%release()
      call bas%release()
   end function e2_at

   !
   ! Runs the FD sweep once for `tag`/`ndens` and checks BOTH properties
   ! this program has always checked: the FD/analytic agreement, then (if
   ! that passed) the translational-invariance residual.
   !
   subroutine run_one(tag, ndens, error)
      character(*), intent(in) :: tag
      integer, intent(in) :: ndens
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      integer :: nat, nao, ia, d, id
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), rp(:, :), rm(:, :)
      real(dp), allocatable :: dmat(:, :, :), grad(:, :), inv(:)
      real(dp) :: ep, em, fd, worst_abs, worst_rel, ref, denom
      integer :: nbad

      call basis_arg(basis_json, error)
      if (allocated(error)) return

      call read_xyz('water.xyz', nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), "check_rys_grad_fd: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao

      allocate (dmat(nao, nao, ndens), grad(3, nat))
      do id = 1, ndens
         call synth_density(nao, id, dmat(:, :, id))
      end do

      call bas%to_device()
      call eri%build(bas, THRESH)
      call eri%gradient(bas, dmat, grad)

      call eri%release()
      call bas%release()

      allocate (rp(3, nat), rm(3, nat), inv(3))
      worst_abs = 0.0_dp; worst_rel = 0.0_dp; nbad = 0
      do ia = 1, nat
         do d = 1, 3
            rp = r; rm = r
            rp(d, ia) = rp(d, ia) + H
            rm(d, ia) = rm(d, ia) - H
            ep = e2_at(trim(basis_json), nat, at_z, rp, nao, ndens, dmat, error)
            if (allocated(error)) then
               deallocate (at_z, r, rp, rm, inv, dmat, grad)
               return
            end if
            em = e2_at(trim(basis_json), nat, at_z, rm, nao, ndens, dmat, error)
            if (allocated(error)) then
               deallocate (at_z, r, rp, rm, inv, dmat, grad)
               return
            end if
            fd = (ep - em)/(2.0_dp*H)
            ref = grad(d, ia)
            worst_abs = max(worst_abs, abs(fd - ref))
            denom = max(abs(ref), abs(fd), 1.0e-10_dp)
            worst_rel = max(worst_rel, abs(fd - ref)/denom)
            if (abs(fd - ref)/denom > RELTOL) nbad = nbad + 1
         end do
      end do

      inv = 0.0_dp
      do ia = 1, nat
         inv = inv + grad(:, ia)
      end do

      call check(error, nbad == 0, tag//": worst abs "//to_string(worst_abs, "es10.2")// &
                 " worst rel "//to_string(worst_rel, "es10.2")//" outside tol "//to_string(nbad))
      if (allocated(error)) then
         deallocate (at_z, r, rp, rm, inv, dmat, grad)
         return
      end if
      call check(error, maxval(abs(inv)) <= INVTOL, tag// &
                 ": translational-invariance residual "//to_string(maxval(abs(inv)), "es10.2")// &
                 " exceeds INVTOL")

      deallocate (at_z, r, rp, rm, inv, dmat, grad)
   end subroutine run_one

   subroutine test_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_one('RHF', 1, error)
   end subroutine test_rhf

   subroutine test_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_one('UHF', 2, error)
   end subroutine test_uhf

end module check_rys_grad_fd_suite

program check_rys_grad_fd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_rys_grad_fd_suite, only: collect_rys_grad_fd_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_rys_grad_fd", collect_rys_grad_fd_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rys_grad_fd
