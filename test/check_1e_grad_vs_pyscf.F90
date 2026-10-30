!
! trc_1e_gradient against pyscf, each one-electron gradient term SEPARATELY,
! at a FIXED synthetic W (energy-weighted density) and D (density) -- same
! kernel-check reasoning as check_rys_grad_vs_pyscf.F90: this is "is d(S|T|V)
! /dR right", not an SCF check, so the density is deterministic and
! non-variational.
!
! SEPARATING THE THREE TERMS. trc_1e_gradient always contracts -W:dS/dR +
! D:dT/dR + D:dV/dR together (see its own header: no derivative matrix is
! ever materialised, so there is nothing to pull a single term out of after
! the call). Three calls isolate them instead, by zeroing whichever input
! makes the unwanted terms vanish identically:
!
!   dS   : W = synthetic, D = 0            -- kills T and V (both need D)
!   dT   : W = 0, D = synthetic, Z = 0      -- kills S (needs W) and V (needs Z)
!   dT+dV: W = 0, D = synthetic, Z = real   -- dV = (dT+dV) - dT
!
! Z = 0 is `bas%at_z` overwritten after `trc_basis_from_json` builds the
! basis from the REAL elements (so the basis functions themselves are
! unchanged) and restored before the third call -- `!$acc update
! device(bas%at_z)` both times, since trc_1e_gradient reads it off the
! device copy.
!
! THE OPERATOR TERM. dV above is the Hellmann-Feynman term ALREADY SUMMED
! with the basis-function derivative (trc_1e_gradient's ga/gb from the bra/
! ket shift, gc from every nucleus, added together by the driver before
! this program ever sees them) -- there is no separate hook to pull the
! Hellmann-Feynman piece out on its own, so the reference below is built
! the same way: pyscf's int1e_ipnuc (basis derivative) plus its own
! `with_rinv_at_nucleus` Hellmann-Feynman loop (test/pyscf_refs/generate.py's
! gen_1e_grad), summed per atom before the comparison.
!
! One ctest case per executable run: the molecule/basis/dump-path/case-name
! come from the command line (test/CMakeLists.txt's `foreach(case
! water_631gs water_ccpvdz formamide_ccpvdz ...)`), so this program checks
! exactly the one case it was invoked with, the way it always has.
!
module check_1e_grad_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_1e_gradient
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix
   implicit none
   private
   public :: collect_1e_grad_vs_pyscf_tests

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: TOL = 1.0e-9_dp

contains

   subroutine collect_1e_grad_vs_pyscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overlap and dS/dT/dV terms vs pyscf", test_case) &
                  ]
   end subroutine collect_1e_grad_vs_pyscf_tests

   !
   ! Same deterministic, dense, symmetric construction as
   ! check_rys_grad_fd.F90/check_rys_grad_vs_pyscf.F90, so every gradient
   ! check in this tree exercises the same kind of density.
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

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=256) :: basis_json, xyzfile, dump, case_name
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(error_t) :: err
      integer :: nat, nao, unit
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), zsave(:)
      real(dp), allocatable :: wmat(:, :), dmat(:, :), zero(:, :)
      real(dp), allocatable :: gs(:, :), gt(:, :), gv(:, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)

      call check(error, command_argument_count() >= 3, &
                 "usage: check_1e_grad_vs_pyscf <xyz> <basis.json> <dump.bin> [case_name]")
      if (allocated(error)) return
      call get_command_argument(1, xyzfile)
      call get_command_argument(2, basis_json)
      call get_command_argument(3, dump)
      case_name = 'default'
      if (command_argument_count() >= 4) call get_command_argument(4, case_name)

      call read_xyz(trim(xyzfile), nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      call check(error,.not. err%has_error(), trim(case_name)//" must build: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      call bas%to_device()
      call pl%build(bas, THRESH)
      call pl%to_device()

      allocate (wmat(nao, nao), dmat(nao, nao), zero(nao, nao))
      call synth_density(nao, 1, wmat)
      call synth_density(nao, 2, dmat)
      zero = 0.0_dp

      allocate (gs(3, nat), gt(3, nat), gv(3, nat), zsave(nat))
      gs = 0.0_dp; gt = 0.0_dp; gv = 0.0_dp

      !$acc enter data copyin(wmat, dmat, zero)

      ! --- dS: W only ---
      !$acc enter data copyin(gs)
      call trc_1e_gradient(bas, pl, wmat, zero, gs)
      !$acc update self(gs)
      !$acc exit data delete(gs)

      ! --- dT: D only, Z = 0 ---
      zsave = bas%at_z
      bas%at_z = 0.0_dp
      !$acc update device(bas%at_z)
      !$acc enter data copyin(gt)
      call trc_1e_gradient(bas, pl, zero, dmat, gt)
      !$acc update self(gt)
      !$acc exit data delete(gt)
      bas%at_z = zsave
      !$acc update device(bas%at_z)

      ! --- dT + dV: D only, Z real; dV = this - dT ---
      !$acc enter data copyin(gv)
      call trc_1e_gradient(bas, pl, zero, dmat, gv)
      !$acc update self(gv)
      !$acc exit data delete(gv)
      gv = gv - gt

      !$acc exit data delete(wmat, dmat, zero)

      ! the overlap, for rys_grad_ref.py-style AO normalisation (see
      ! check_rys_grad_vs_pyscf.F90's own comment -- unchanged here)
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
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
      write (unit) wmat
      write (unit) dmat
      write (unit) gs
      write (unit) gt
      write (unit) gv
      write (unit) smat
      close (unit)

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      block
         integer :: u, nat_ref, nao_ref
         real(dp), allocatable :: s_ref(:, :), gs_ref(:, :), gt_ref(:, :), gv_ref(:, :)
         real(dp) :: ds, worst_s, worst_t, worst_v
         call ref_require('onee_grad_'//trim(case_name)//'.ref', u)
         call ref_read_int(u, 'NAT', nat_ref)
         call ref_read_int(u, 'NAO', nao_ref)
         call check(error, nat_ref == nat .and. nao_ref == nao, &
                    trim(case_name)//": pyscf reference shape does not match this build")
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
         allocate (s_ref(nao, nao), gs_ref(3, nat), gt_ref(3, nat), gv_ref(3, nat))
         call ref_read_matrix(u, 'S', nao, nao, s_ref)
         call ref_read_matrix(u, 'GRAD_S', 3, nat, gs_ref)
         call ref_read_matrix(u, 'GRAD_T', 3, nat, gt_ref)
         call ref_read_matrix(u, 'GRAD_V', 3, nat, gv_ref)
         call ref_close(u)

         ds = maxval(abs(smat - s_ref))
         call check(error, ds, 0.0_dp, thr=TOL, message=trim(case_name)// &
                    ": overlap check vs pyscf, max|diff| = "//to_string(ds, "es12.4"))
         if (allocated(error)) then
            deallocate (s_ref, gs_ref, gt_ref, gv_ref)
            return
         end if

         worst_s = maxval(abs(gs - gs_ref))
         call check(error, worst_s, 0.0_dp, thr=TOL, message=trim(case_name)// &
                    ": dS gradient vs pyscf, worst|diff| = "//to_string(worst_s, "es12.4"))
         if (allocated(error)) then
            deallocate (s_ref, gs_ref, gt_ref, gv_ref)
            return
         end if

         worst_t = maxval(abs(gt - gt_ref))
         call check(error, worst_t, 0.0_dp, thr=TOL, message=trim(case_name)// &
                    ": dT gradient vs pyscf, worst|diff| = "//to_string(worst_t, "es12.4"))
         if (allocated(error)) then
            deallocate (s_ref, gs_ref, gt_ref, gv_ref)
            return
         end if

         worst_v = maxval(abs(gv - gv_ref))
         call check(error, worst_v, 0.0_dp, thr=TOL, message=trim(case_name)// &
                    ": dV gradient vs pyscf, worst|diff| = "//to_string(worst_v, "es12.4"))
         deallocate (s_ref, gs_ref, gt_ref, gv_ref)
      end block

      deallocate (at_z, r, wmat, dmat, zero, gs, gt, gv, zsave, smat, tmat, vmat)
   end subroutine test_case

end module check_1e_grad_vs_pyscf_suite

program check_1e_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_1e_grad_vs_pyscf_suite, only: collect_1e_grad_vs_pyscf_tests
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [new_testsuite("check_1e_grad_vs_pyscf", collect_1e_grad_vs_pyscf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_1e_grad_vs_pyscf
