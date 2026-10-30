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
module check_rys_grad_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix
   implicit none
   private
   public :: collect_rys_grad_vs_pyscf, basis_json_arg, dump_arg, case_name_arg

   real(dp), parameter :: THRESH = 1.0e-14_dp

   character(len=256) :: basis_json_arg = ''
   character(len=256) :: dump_arg = ''
   character(len=256) :: case_name_arg = 'ccpvdz'
   !! Set by the program from argv before run_testsuite: ONE basis per
   !! process, same reasoning as check_rys_vs_hgp's own module-level args.

contains

   subroutine collect_rys_grad_vs_pyscf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF/UHF gradient vs pyscf", test_vs_pyscf) &
                  ]
   end subroutine collect_rys_grad_vs_pyscf

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

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      integer :: nat, nao, unit
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: drhf(:, :, :), duhf(:, :, :)
      real(dp), allocatable :: grad_rhf(:, :), grad_uhf(:, :)
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)

      call read_xyz('water.xyz', nat, at_z, r)
      call trc_basis_from_json(trim(basis_json_arg), nat, at_z, r, bas, err)
      call check(error, .not. err%has_error(), "basis build: "//err%get_message())
      if (allocated(error)) return
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

      open (newunit=unit, file=trim(dump_arg), access='stream', &
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
      print '(a,i0,a,i0,a)', 'check_rys_grad_vs_pyscf: wrote '//trim(dump_arg)//' (', &
         nat, ' atoms, ', nao, ' AOs) for the pyscf-reference comparison'

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      block
         integer :: u, nat_ref, nao_ref
         real(dp), allocatable :: s_ref(:, :), grad_rhf_ref(:, :), grad_uhf_ref(:, :)
         real(dp) :: ds, worst_rhf, worst_uhf
         call ref_require('rys_grad_'//trim(case_name_arg)//'.ref', u)
         call ref_read_int(u, 'NAT', nat_ref)
         call ref_read_int(u, 'NAO', nao_ref)
         call check(error, nat_ref == nat .and. nao_ref == nao, &
                    "pyscf reference shape does not match this build")
         if (allocated(error)) then
            call ref_close(u)
            deallocate (at_z, r, drhf, duhf, grad_rhf, grad_uhf, smat, tmat, vmat)
            return
         end if
         allocate (s_ref(nao, nao), grad_rhf_ref(3, nat), grad_uhf_ref(3, nat))
         call ref_read_matrix(u, 'S', nao, nao, s_ref)
         call ref_read_matrix(u, 'GRAD_RHF', 3, nat, grad_rhf_ref)
         call ref_read_matrix(u, 'GRAD_UHF', 3, nat, grad_uhf_ref)
         call ref_close(u)

         ds = maxval(abs(smat - s_ref))
         worst_rhf = maxval(abs(grad_rhf - grad_rhf_ref))
         worst_uhf = maxval(abs(grad_uhf - grad_uhf_ref))

         call check(error, ds < 1.0e-9_dp, &
                    "overlap check vs pyscf, max|diff| = "//to_string(ds, "es10.2"))
         if (.not. allocated(error)) then
            call check(error, worst_rhf < 1.0e-9_dp, &
                       "RHF gradient vs pyscf, worst|diff| = "//to_string(worst_rhf, "es10.2"))
         end if
         if (.not. allocated(error)) then
            call check(error, worst_uhf < 1.0e-9_dp, &
                       "UHF gradient vs pyscf, worst|diff| = "//to_string(worst_uhf, "es10.2"))
         end if
         deallocate (s_ref, grad_rhf_ref, grad_uhf_ref)
      end block

      deallocate (at_z, r, drhf, duhf, grad_rhf, grad_uhf, smat, tmat, vmat)
   end subroutine test_vs_pyscf

end module check_rys_grad_vs_pyscf_suite

program check_rys_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rys_grad_vs_pyscf_suite, only: collect_rys_grad_vs_pyscf, basis_json_arg, dump_arg, case_name_arg
   implicit none
   integer :: stat

   if (command_argument_count() < 2) then
      write (error_unit, '(a)') 'usage: check_rys_grad_vs_pyscf <basis.json> <dump.bin> [case_name]'
      error stop 1
   end if
   call get_command_argument(1, basis_json_arg)
   call get_command_argument(2, dump_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, case_name_arg)

   stat = 0
   call run_testsuite(collect_rys_grad_vs_pyscf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rys_grad_vs_pyscf
