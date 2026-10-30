!
! The long-range exchange gradient term alone -- eri%gradient(j_scale=0,
! k_scale=1, omega) -- against pyscf's own range-separated int2e_ip1
! (mol.with_range_coulomb(omega)), contracted by hand exactly the way
! check_rys_grad_vs_pyscf.F90 contracts the full-range int2e_ip1: a fixed,
! deterministic, non-SCF density (pyscf's Gradients class assumes a
! converged SCF density, which does not apply here), RHF only -- a
! range-separated UHF gradient is out of scope on this branch (feat/grad-1e
! has no exx-scaled full-range term either, and mqc's own reference
! implementation, backends/cenzontle/derivatives/mqc_czt_gradient.f90,
! refuses UHF + range-separated for the same reason: no validated
! convention to check against).
!
! TWO omega values, both real functionals already in this tree
! (src/xc/trc_xc_functional.F90): 0.3 (wB97X-V/wB97X) and 0.4 (wB97). k_scale
! is 1 here regardless -- the kernel is linear in kfac by construction (see
! scripts/gen_perclass.py's own prescreen comment), so checking it at 1
! checks it at any exx_lr once this passes.
!
! eri%gradient's own convention (this module's header comment, same as
! check_rys_grad_vs_pyscf.F90): TWICE d(e_two)/dR, no halving -- so the
! pyscf side below is contracted with the identical "*4.0" (not "*2.0")
! factor check_rys_grad_vs_pyscf.F90's grad2e_rhf uses, and the comparison
! is against that raw number on both sides.
!
! Basis (6-31G* or cc-pVDZ) and dump path come from the command line
! (test/CMakeLists.txt's `foreach(b 631gs ccpvdz)`), so each run checks
! exactly the one case it was invoked with.
!
module check_lr_exchange_grad_vs_pyscf_suite
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
   public :: collect_lr_exchange_grad_vs_pyscf_tests

   integer, parameter :: NOMEGA = 2
   real(dp), parameter :: OMEGA_LIST(NOMEGA) = [0.3_dp, 0.4_dp]
   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: TOL = 1.0e-9_dp

contains

   subroutine collect_lr_exchange_grad_vs_pyscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overlap and long-range exchange gradient vs pyscf", test_case) &
                  ]
   end subroutine collect_lr_exchange_grad_vs_pyscf_tests

   !
   ! Same deterministic, dense, symmetric construction as
   ! check_rys_grad_vs_pyscf.F90/check_rys_grad_fd.F90, phase 1.
   !
   subroutine synth_density(n, d)
      integer, intent(in) :: n
      real(dp), intent(out) :: d(n, n)
      integer :: ii, jj
      do jj = 1, n
         do ii = 1, n
            d(ii, jj) = 0.4_dp + 0.1_dp*cos(real(3*ii + 5*jj + 7, dp))
         end do
      end do
      d = 0.5_dp*(d + transpose(d))
   end subroutine synth_density

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json, dump, case_name
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(error_t) :: err
      integer :: nat, nao, io
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: drhf(:, :, :)
      real(dp), allocatable :: grad_lr(:, :, :)
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)

      call check(error, command_argument_count() >= 2, &
                 "usage: check_lr_exchange_grad_vs_pyscf <basis.json> <dump.bin>")
      if (allocated(error)) return
      call get_command_argument(1, basis_json)
      call get_command_argument(2, dump)
      case_name = 'ccpvdz'
      if (command_argument_count() >= 3) call get_command_argument(3, case_name)

      call read_xyz('water.xyz', nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), "check_lr_exchange_grad_vs_pyscf: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao

      allocate (drhf(nao, nao, 1), grad_lr(3, nat, NOMEGA))
      call synth_density(nao, drhf(:, :, 1))

      call bas%to_device()
      call eri%build(bas, THRESH)
      do io = 1, NOMEGA
         call eri%gradient(bas, drhf, grad_lr(:, :, io), &
                           j_scale=0.0_dp, k_scale=1.0_dp, omega=OMEGA_LIST(io))
      end do
      call eri%release()

      ! The overlap, for the reference generator's AO-normalisation check --
      ! same reasoning as check_rys_grad_vs_pyscf.F90's own smat dump.
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      call pl%build(bas, THRESH)
      call pl%to_device()
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(smat)
      !$acc exit data delete(smat, tmat, vmat)
      call pl%release()
      call bas%release()

      open (newunit=io, file=trim(dump), access='stream', &
            form='unformatted', status='replace')
      write (io) int(nat, int64), int(nao, int64)
      write (io) int(at_z, int64)
      write (io) r
      write (io) drhf(:, :, 1)
      write (io) grad_lr(:, :, 1)
      write (io) grad_lr(:, :, 2)
      write (io) smat
      close (io)

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      block
         integer :: u, nat_ref, nao_ref
         real(dp), allocatable :: s_ref(:, :), grad_ref(:, :, :)
         real(dp) :: ds, worst
         call ref_require('lr_exchange_grad_'//trim(case_name)//'.ref', u)
         call ref_read_int(u, 'NAT', nat_ref)
         call ref_read_int(u, 'NAO', nao_ref)
         call check(error, nat_ref == nat .and. nao_ref == nao, &
                    "check_lr_exchange_grad_vs_pyscf: pyscf reference shape does not match this build")
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
         allocate (s_ref(nao, nao), grad_ref(3, nat, NOMEGA))
         call ref_read_matrix(u, 'S', nao, nao, s_ref)
         call ref_read_matrix(u, 'GRAD_LR_OM1', 3, nat, grad_ref(:, :, 1))
         call ref_read_matrix(u, 'GRAD_LR_OM2', 3, nat, grad_ref(:, :, 2))
         call ref_close(u)

         ds = maxval(abs(smat - s_ref))
         call check(error, ds, 0.0_dp, thr=TOL, message= &
                    "overlap check vs pyscf, max|diff| = "//to_string(ds, "es12.4"))
         if (allocated(error)) then
            deallocate (s_ref, grad_ref)
            return
         end if

         do io = 1, NOMEGA
            worst = maxval(abs(grad_lr(:, :, io) - grad_ref(:, :, io)))
            call check(error, worst, 0.0_dp, thr=TOL, message= &
                       "omega "//to_string(OMEGA_LIST(io), "f4.2")// &
                       " LR exchange gradient vs pyscf, worst|diff| = "//to_string(worst, "es12.4"))
            if (allocated(error)) then
               deallocate (s_ref, grad_ref)
               return
            end if
         end do
         deallocate (s_ref, grad_ref)
      end block

      deallocate (at_z, r, drhf, grad_lr, smat, tmat, vmat)
   end subroutine test_case

end module check_lr_exchange_grad_vs_pyscf_suite

program check_lr_exchange_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_lr_exchange_grad_vs_pyscf_suite, only: collect_lr_exchange_grad_vs_pyscf_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_lr_exchange_grad_vs_pyscf", collect_lr_exchange_grad_vs_pyscf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_lr_exchange_grad_vs_pyscf
