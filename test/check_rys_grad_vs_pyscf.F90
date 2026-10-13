!
! eri%gradient against pyscf's int2e_ip1, contracted with the same Gamma,
! RHF and UHF, on water. Writes the dump named by its second argument;
! rys_grad_ref.py reads it, builds the same geometry and the same BSE JSON
! basis in pyscf (cart=True, so the AO order already matches this build's
! lexicographic one -- no cart_perm needed), and does the contraction itself rather than
! trusting pyscf's own Gradients class, which assumes a converged SCF
! density (the energy-weighted density term) that does not apply here: D
! is the same fixed, synthetic, non-SCF density check_rys_grad_fd.F90
! uses, for the same reason -- this is a kernel check on an arbitrary
! density, not an SCF check. Also writes terco's overlap matrix, for
! rys_grad_ref.py's per-AO normalisation check/fix (see its own comment).
!
program check_rys_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: int64
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   real(dp), parameter :: THRESH = 1.0e-14_dp
   character(len=256) :: basis_json, dump
   type(trc_basis_t) :: bas
   type(trc_eri_t)   :: eri
   type(error_t) :: err
   integer :: nat, nao, i, j, id, unit
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :)
   real(dp), allocatable :: drhf(:, :, :), duhf(:, :, :)
   real(dp), allocatable :: grad_rhf(:, :), grad_uhf(:, :)
   type(trc_pairlist_t) :: pl
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)

   if (command_argument_count() < 2) then
      print '(a)', 'usage: check_rys_grad_vs_pyscf <basis.json> <dump.bin>'; stop 1
   end if
   call get_command_argument(1, basis_json)
   call get_command_argument(2, dump)

   call read_xyz('water.xyz', nat, at_z, r)
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
   if (err%has_error()) then
      print '(a)', 'check_rys_grad_vs_pyscf: '//err%get_message(); stop 1
   end if
   nao = bas%nao

   allocate (drhf(nao, nao, 1), duhf(nao, nao, 2))
   allocate (grad_rhf(3, nat), grad_uhf(3, nat))
   call synth_density(nao, 1, drhf(:, :, 1))
   call synth_density(nao, 1, duhf(:, :, 1))
   call synth_density(nao, 2, duhf(:, :, 2))

   call bas%to_device()
   call eri%build(bas, THRESH)
   call eri%gradient(bas, drhf, grad_rhf)
   call eri%gradient(bas, duhf, grad_uhf)
   call eri%release()

   ! The overlap, for rys_grad_ref.py's AO-normalisation check (a fixed
   ! density is NOT invariant to a per-AO rescaling the way an SCF energy
   ! is, and pyscf's cart=True d functions are not all unit-normalised the
   ! way this build's are -- off-diagonal xy/xz/yz carry a different factor
   ! from xx/yy/zz). S_terco = N S_pyscf N with n_i =
   ! sqrt(S_terco(i,i)/S_pyscf(i,i)) is the check; D_pyscf = N D_terco N is
   ! the fix, if the check holds.
   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
   call pl%build(bas, THRESH)
   call pl%to_device()
   !$acc enter data create(smat, tmat, vmat)
   call trc_1e(bas, pl, smat, tmat, vmat)
   !$acc update self(smat)
   !$acc exit data delete(smat, tmat, vmat)
   call pl%release()
   call bas%release()

   open (newunit=unit, file=trim(dump), access='stream', &
         form='unformatted', status='replace')
   write (unit) int(nat, kind=8), int(nao, kind=8)
   write (unit) int(at_z, kind=8)
   write (unit) r
   write (unit) drhf(:, :, 1)
   write (unit) duhf(:, :, 1), duhf(:, :, 2)
   write (unit) grad_rhf
   write (unit) grad_uhf
   write (unit) smat
   close (unit)
   print '(a,i0,a,i0,a)', 'check_rys_grad_vs_pyscf: wrote '//trim(dump)//' (', &
      nat, ' atoms, ', nao, ' AOs) for rys_grad_ref.py'

   deallocate (at_z, r, drhf, duhf, grad_rhf, grad_uhf)

contains

   !
   ! Same deterministic, dense, symmetric construction as
   ! check_rys_grad_fd.F90, so both checks exercise the same density.
   !
   subroutine synth_density(n, phase, d)
      integer, intent(in) :: n, phase
      real(dp), intent(out) :: d(n, n)
      integer :: ii, jj
      do jj = 1, n
         do ii = 1, n
            d(ii, jj) = 0.4_dp + 0.1_dp*cos(real(3*ii + 5*jj + 7*phase, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

end program check_rys_grad_vs_pyscf
