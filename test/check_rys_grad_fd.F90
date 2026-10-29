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
program check_rys_grad_fd
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   real(dp), parameter :: H = 1.0e-4_dp
   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: RELTOL = 1.0e-7_dp
   real(dp), parameter :: INVTOL = 1.0e-10_dp

   character(len=256) :: basis_json
   logical :: any_fail

   any_fail = .false.
   if (command_argument_count() < 1) then
      print '(a)', 'usage: check_rys_grad_fd <basis.json>'; stop 1
   end if
   call get_command_argument(1, basis_json)

   call run_rhf(trim(basis_json))
   call run_uhf(trim(basis_json))

   if (any_fail) then
      print '(a)', '  RESULT: FAIL'
      stop 1
   end if
   print '(a)', '  RESULT: PASS'

contains

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
   real(dp) function e2_at(basis_json, nat, at_z, r, nao, ndens, dmat)
      character(*), intent(in) :: basis_json
      integer, intent(in) :: nat, at_z(nat), nao, ndens
      real(dp), intent(in) :: r(3, nat), dmat(nao, nao, ndens)
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      real(dp) :: vhf(nao, nao), dtot(nao, nao)

      ! overrides%prune = .false.: this check compares the analytic gradient against a
      ! finite difference of the SAME path's own energy, at a RELTOL margin
      ! tight enough that a near-zero gradient component (water's symmetry
      ! plane) is sensitive to which exact rounding the kernels take --
      ! unrelated to correctness, but a reason to keep this comparison on
      ! the arithmetic path it was tuned against rather than
      ! trc_basis_json's prune_block reduction.
      call trc_basis_from_json(basis_json, nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (err%has_error()) then
         print '(a)', 'check_rys_grad_fd: '//err%get_message(); stop 1
      end if
      if (bas%nao /= nao) then
         print '(a)', 'check_rys_grad_fd: nao changed under displacement'; stop 1
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

   subroutine run_one(basis_json, tag, ndens)
      character(*), intent(in) :: basis_json, tag
      integer, intent(in) :: ndens
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      integer :: nat, nao, ia, d, id
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), rp(:, :), rm(:, :)
      real(dp), allocatable :: dmat(:, :, :), grad(:, :), inv(:)
      real(dp) :: ep, em, fd, worst_abs, worst_rel, ref, denom
      integer :: nbad

      call read_xyz('water.xyz', nat, at_z, r)
      ! overrides%prune = .false.: same reasoning as e2_at above, and it must agree
      ! with it -- e2_at rebuilds the basis at each displaced geometry.
      call trc_basis_from_json(basis_json, nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (err%has_error()) then
         print '(a)', 'check_rys_grad_fd: '//err%get_message(); stop 1
      end if
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
            ep = e2_at(basis_json, nat, at_z, rp, nao, ndens, dmat)
            em = e2_at(basis_json, nat, at_z, rm, nao, ndens, dmat)
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

      print '(a,a,a,i0,a,es10.2,a,es10.2,a,i0)', '  ', tag, ': nat ', nat, &
         '  worst abs ', worst_abs, '  worst rel ', worst_rel, '  outside tol ', nbad
      print '(a,3es10.2)', '    translational-invariance residual: ', inv

      if (nbad > 0) any_fail = .true.
      if (maxval(abs(inv)) > INVTOL) any_fail = .true.

      deallocate (at_z, r, rp, rm, inv, dmat, grad)
   end subroutine run_one

   subroutine run_rhf(basis_json)
      character(*), intent(in) :: basis_json
      call run_one(basis_json, 'RHF', 1)
   end subroutine run_rhf

   subroutine run_uhf(basis_json)
      character(*), intent(in) :: basis_json
      call run_one(basis_json, 'UHF', 2)
   end subroutine run_uhf

end program check_rys_grad_fd
