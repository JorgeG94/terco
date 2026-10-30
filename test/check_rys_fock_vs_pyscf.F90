!
! eri%fock (the two-electron Fock build, J - K/2 for the total density, see
! src/eri/trc_eri.F90's own header for that convention) against pyscf's raw
! int2e, on a fixed, deterministic, non-SCF density -- the first class to
! exercise an f shell (l=3) in a four-centre test.
!
! Modelled on check_rys_vs_hgp.F90 (the fixed-density Fock comparison this
! reuses almost verbatim) and check_rys_grad_vs_pyscf.F90 (the dump/pyscf-
! reference split, the overlap check for AO normalisation, and the AO-order
! permutation -- pyscf's cart=True Mole.build() does not keep a JSON's
! literal shell order, see rys_fock_ref.py's own header). Writes the dump
! named by its third argument; rys_fock_ref.py reads it, rebuilds the same
! geometry and basis in pyscf, and does the comparison.
!
! SKIPS CLEANLY when the basis needs l = 3 and this build's LMAX does not
! reach it (TRC_LMAX < 3), exactly as check_nosym.F90 does for check_nosym_f
! -- `stop 0` after printing SKIP, so ctest sees a pass rather than a red
! build that never had the kernels to run. Registered at every ceiling for
! the same reason check_nosym_f is: a build that CAN reach f always runs it,
! and nothing here depends on how that ceiling gets raised. This self-skip
! happens in run_case, BEFORE run_testsuite is ever called, so the stop-0/
! "SKIP" text ctest's SKIP_REGULAR_EXPRESSION matches is unchanged.
!
! At LMAX=2 against cc-pVDZ or 6-31G* (no f on H/C/N/O), this proves the
! OTHER half of the plumbing that does not need l=3 at all: the JSON basis
! build, the AO permutation, and the per-AO normalisation this comparison
! shares with every other scheme working around the same pyscf cart=True
! quirk.
!
module check_rys_fock_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_tables, only: LMAX
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   use golden_io, only: golden_begin, golden_fock, golden_end
   implicit none
   private
   public :: collect_rys_fock_tests, run_case

   real(dp), parameter :: THRESH = 1.0e-14_dp

   real(dp) :: g_ds, g_worst, g_tol
      !! overlap and Fock max|diff| against pyscf, and the tolerance the
      !! command line selected, filled in by run_case.

contains

   subroutine collect_rys_fock_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overlap and Fock build match pyscf", test_vs_pyscf) &
                  ]
   end subroutine collect_rys_fock_tests

   !
   ! Everything the old program did: argument parsing, the LMAX self-skip,
   ! the fixed-density Fock build (plus its golden dump, left exactly where
   ! it was relative to the build), the AO-normalisation overlap, writing
   ! the probe dump, and reading the pyscf reference back into the two
   ! max|diff| values test_vs_pyscf checks.
   !
   subroutine run_case(skip, message)
      logical, intent(out) :: skip
      character(len=*), intent(out) :: message

      character(len=256) :: xyzfile, basis_json, dump, tol_arg, case_name
      type(trc_basis_t) :: bas
      type(trc_eri_t)   :: eri
      type(trc_pairlist_t) :: pl
      type(error_t) :: err
      integer :: nat, nao, unit
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: d(:, :), g(:, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      integer :: u, nat_ref, nao_ref, nfail
      real(dp), allocatable :: s_ref(:, :), g_ref(:, :)

      skip = .false.
      message = ''

      if (command_argument_count() < 3) then
         print '(a)', 'usage: check_rys_fock_vs_pyscf <xyz> <basis.json> <dump.bin>'; stop 1
      end if
      call get_command_argument(1, xyzfile)
      call get_command_argument(2, basis_json)
      call get_command_argument(3, dump)
      g_tol = 1.0e-9_dp
      if (command_argument_count() >= 4) then
         call get_command_argument(4, tol_arg)
         read (tol_arg, *) g_tol
      end if
      ! The case name, for the pyscf reference's own filename
      ! (test/pyscf_refs/rys_fock_<case>.ref) -- derived from the dump's own
      ! basename, which test/CMakeLists.txt always names
      ! rys_fock_probe_<case>.bin, rather than a fifth CLI argument that would
      ! just repeat it.
      block
         integer :: islash, idot
         islash = index(trim(dump), '/', back=.true.)
         idot = index(trim(dump), '.bin', back=.true.)
         case_name = dump(islash + 1 + len('rys_fock_probe_'):idot - 1)
      end block
      call golden_begin('rys_fock_'//trim(case_name))

      call read_xyz(trim(xyzfile), nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      if (err%has_error()) then
         print '(a)', 'check_rys_fock_vs_pyscf: '//err%get_message(); stop 1
      end if

      ! See this program's own header: a basis that needs l = 3 is not an
      ! error on a build whose kernels stop at d, it is nothing to test yet.
      if (maxval(bas%sh_l) > LMAX) then
         write (message, '(a,i0,a,i0,a)') '  SKIP: basis needs l = ', maxval(bas%sh_l), &
            ' and this build has LMAX = ', LMAX, '.'
         skip = .true.; return
      end if

      nao = bas%nao
      allocate (d(nao, nao), g(nao, nao))
      call synth_density(nao, d)

      call bas%to_device()
      call eri%build(bas, THRESH)
      call eri%fock(bas, d, g, density_screen=.false.)
      call eri%release()
      call golden_fock('G', d, g)
      call golden_end()

      ! The overlap, for rys_fock_ref.py's per-AO normalisation check -- same
      ! reasoning as check_rys_grad_vs_pyscf.F90's own smat dump: a fixed
      ! density (and the G it produces) is not invariant to a per-AO rescaling
      ! the way an SCF energy is, and pyscf's cart=True d/f functions are not
      ! all unit-normalised the way this build's are.
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
      write (unit) d
      write (unit) g
      write (unit) smat
      close (unit)
      print '(a,i0,a,i0,a)', 'check_rys_fock_vs_pyscf: wrote '//trim(dump)//' (', &
         nat, ' atoms, ', nao, ' AOs) for the pyscf-reference comparison'

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      call ref_require('rys_fock_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      if (nat_ref /= nat .or. nao_ref /= nao) then
         print '(a)', 'check_rys_fock_vs_pyscf: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (s_ref(nao, nao), g_ref(nao, nao))
      call ref_read_matrix(u, 'S', nao, nao, s_ref)
      call ref_read_matrix(u, 'G', nao, nao, g_ref)
      call ref_close(u)
      nfail = 0
      g_ds = maxval(abs(smat - s_ref))
      call ref_report(nfail, g_ds < g_tol, 'overlap check vs pyscf, max|diff|', g_ds)
      g_worst = maxval(abs(g - g_ref))
      call ref_report(nfail, g_worst < g_tol, 'G vs pyscf, worst|diff|', g_worst)

      deallocate (at_z, r, d, g, smat, tmat, vmat, s_ref, g_ref)
   end subroutine run_case

   !
   ! Same deterministic, dense, symmetric construction as check_nosym.F90
   ! and check_rys_vs_hgp.F90 -- not sparse, not SCF-converged, so the
   ! comparison exercises the whole digestion rather than whatever a
   ! converged density happens to leave non-negligible.
   !
   subroutine synth_density(n, dd)
      integer, intent(in) :: n
      real(dp), intent(out) :: dd(n, n)
      integer :: ii, jj
      do jj = 1, n
         do ii = 1, n
            dd(ii, jj) = 0.4_dp + 0.1_dp*cos(real(3*ii + 5*jj, dp))
         end do
      end do
      dd = 0.5_dp*(dd + transpose(dd))
   end subroutine synth_density

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, g_ds < g_tol, "overlap check vs pyscf, max|diff| = "// &
                 to_string(g_ds, "es12.4")//" exceeds "//to_string(g_tol, "es12.4"))
      if (allocated(error)) return
      call check(error, g_worst < g_tol, "G vs pyscf, worst|diff| = "// &
                 to_string(g_worst, "es12.4")//" exceeds "//to_string(g_tol, "es12.4"))
   end subroutine test_vs_pyscf

end module check_rys_fock_vs_pyscf_suite

program check_rys_fock_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_rys_fock_vs_pyscf_suite, only: collect_rys_fock_tests, run_case
   implicit none
   integer :: stat
   logical :: skip
   character(len=256) :: message
   type(testsuite_type), allocatable :: testsuites(:)

   call run_case(skip, message)
   if (skip) then
      print '(a)', trim(message)
      stop 0
   end if

   stat = 0
   testsuites = [new_testsuite("check_rys_fock_vs_pyscf", collect_rys_fock_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_rys_fock_vs_pyscf
