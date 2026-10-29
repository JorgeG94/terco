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
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_overrides_t
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   implicit none

   real(dp), parameter :: THRESH = 1.0e-14_dp
   character(len=256) :: basis_json, dump, case_name
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
   case_name = 'ccpvdz'
   if (command_argument_count() >= 3) call get_command_argument(3, case_name)

   call read_xyz('water.xyz', nat, at_z, r)
   ! overrides%prune = .false.: the pyscf reference is built from this file's literal
   ! columns (bse_json_basis), same reasoning as check_rys_fock_vs_pyscf.F90.
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err, &
                               overrides=trc_basis_overrides_t(prune=.false.))
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
   write (unit) int(nat, int64), int(nao, int64)
   write (unit) int(at_z, int64)
   write (unit) r
   write (unit) drhf(:, :, 1)
   write (unit) duhf(:, :, 1), duhf(:, :, 2)
   write (unit) grad_rhf
   write (unit) grad_uhf
   write (unit) smat
   close (unit)
   print '(a,i0,a,i0,a)', 'check_rys_grad_vs_pyscf: wrote '//trim(dump)//' (', &
      nat, ' atoms, ', nao, ' AOs) for the pyscf-reference comparison'

   ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
   block
      integer :: u, nat_ref, nao_ref, nfail
      real(dp), allocatable :: s_ref(:, :), grad_rhf_ref(:, :), grad_uhf_ref(:, :)
      real(dp) :: ds, worst_rhf, worst_uhf
      call ref_require('rys_grad_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      if (nat_ref /= nat .or. nao_ref /= nao) then
         print '(a)', 'check_rys_grad_vs_pyscf: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (s_ref(nao, nao), grad_rhf_ref(3, nat), grad_uhf_ref(3, nat))
      call ref_read_matrix(u, 'S', nao, nao, s_ref)
      call ref_read_matrix(u, 'GRAD_RHF', 3, nat, grad_rhf_ref)
      call ref_read_matrix(u, 'GRAD_UHF', 3, nat, grad_uhf_ref)
      call ref_close(u)
      nfail = 0
      ds = maxval(abs(smat - s_ref))
      call ref_report(nfail, ds < 1.0e-9_dp, 'overlap check vs pyscf, max|diff|', ds)
      worst_rhf = maxval(abs(grad_rhf - grad_rhf_ref))
      worst_uhf = maxval(abs(grad_uhf - grad_uhf_ref))
      call ref_report(nfail, worst_rhf < 1.0e-9_dp, 'RHF gradient vs pyscf, worst|diff|', worst_rhf)
      call ref_report(nfail, worst_uhf < 1.0e-9_dp, 'UHF gradient vs pyscf, worst|diff|', worst_uhf)
      if (nfail > 0) then
         print '(a)', 'check_rys_grad_vs_pyscf: FAILED against the pyscf reference'; stop 1
      end if
      print '(a)', 'check_rys_grad_vs_pyscf: RESULT: PASS'
   end block

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
