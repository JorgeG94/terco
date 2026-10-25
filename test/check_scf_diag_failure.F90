!
! A malformed SCF guess density -- one entry replaced by a NaN, simulating
! a warm-start density read past the end of an undersized buffer (a
! per-image buffer sized before the first image's own converged nao was
! known, app/trc_dlfind_bridge.F90's own reset_image_state) fed to
! trc_set_guess(TRC_GUESS_GIVEN) -- must come back through trc_run_scf's
! own status/message, not `error stop` inside src/core/trc_linalg.F90's
! eigensolver.
!
! Root cause, read directly off the source (not guessed from the symptom):
! src/scf/trc_scf_driver.F90's per-iteration Fock diagonalisation
! (orbitals_and_density, F' = X^T F X then linalg_provider%syev(fp,
! eps_s)) is reached from a GIVEN density on the very first iteration --
! `dguess` skips the guess-construction path (guess_from) entirely and
! goes straight into the main loop -- so a malformed GIVEN density is the
! most direct, deterministic way to reach exactly the call site a
! real warm-start crash would hit.
!
! A signalling NaN is used to seed the one bad entry (ieee_signaling_nan,
! below) only because it is the straightforward intrinsic for "write a
! deliberately invalid real" -- dsyev reports a nonzero info for a NaN
! input whether it started quiet or signalling, since ordinary arithmetic
! on a signalling NaN quiets it before it ever reaches the eigensolver
! (test/check_linalg_mat.F90%test_syev_failure's own header has the
! measurement).
!
! A test-drive suite.
!
module check_scf_diag_failure_suite
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_null_ptr, c_null_char
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
   use pic_types, only: dp
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                       trc_set_convergence, trc_set_guess, trc_run_scf, trc_nao, trc_converged, &
                       trc_iterations, trc_message, TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_GIVEN, TRC_GUESS_SAD
   use trc_test_basis, only: read_xyz
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_scf_diag_failure

   integer, parameter :: MSG_LEN = 1024

contains

   subroutine collect_scf_diag_failure(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("a NaN guess density comes back through trc_run_scf's own status, not error stop", &
                      test_diag_failure) &
         ]
   end subroutine collect_scf_diag_failure

   subroutine test_diag_failure(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: st, nao, conv, niter
      integer :: nat, nfail
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      real(c_double), allocatable, target :: dbuf(:)
      character(kind=c_char) :: msg_buf(MSG_LEN)

      nfail = 0
      call read_xyz('water.xyz', nat, z, r)

      st = trc_create(h)
      st = trc_set_molecule(h, int(nat, c_int), real(z, c_double), reshape(real(r, c_double), [3*nat]), &
                            0_c_int, 1_c_int)
      st = trc_set_basis_json(h, c_path('basis_sets/sto-3g.json'))
      st = trc_set_method(h, c_path(''), 3_c_int)
      st = trc_set_convergence(h, 1.0e-10_c_double, 1.0e-6_c_double, 100_c_int, 10_c_int)

      ! A normal, converged SCF first -- so the density this test corrupts
      ! is otherwise exactly what a real warm start would have handed back.
      st = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
      st = trc_run_scf(h)
      call check(nfail, st == TRC_OK, 'the baseline SAD SCF converges', real(st, dp))

      st = trc_nao(h, nao)
      call check(nfail, st == TRC_OK, 'trc_nao', real(st, dp))
      allocate (dbuf(int(nao)*int(nao)))
      dbuf = 0.0_c_double
      ! A plausible-looking density: occupied-ish diagonal, nothing exotic
      ! -- the ONE thing wrong with it is a single signalling-NaN
      ! off-diagonal entry (symmetric, as a density must be).
      block
         integer :: i
         do i = 1, int(nao)
            dbuf((i - 1)*int(nao) + i) = 1.0_c_double
         end do
      end block
      block
         real(c_double) :: snan
         integer :: i, j, n
         n = int(nao)
         i = min(2, n); j = min(5, n)
         if (i == j) j = min(i + 1, n)
         snan = ieee_value(1.0_c_double, ieee_signaling_nan)
         dbuf((j - 1)*n + i) = snan
         dbuf((i - 1)*n + j) = snan
      end block

      st = trc_set_guess(h, TRC_GUESS_GIVEN, c_loc(dbuf(1)), 1_c_int)
      call check(nfail, st == TRC_OK, 'trc_set_guess(GIVEN) itself accepts the buffer', real(st, dp))

      st = trc_run_scf(h)
      print '(a,i0)', '  trc_run_scf status = ', st
      call check(nfail, st == TRC_ERR_NOCONV, &
                 'a signalling-NaN guess density comes back as TRC_ERR_NOCONV, not a crash', real(st, dp))

      st = trc_converged(h, conv)
      call check(nfail, conv == 0, 'trc_converged reports false', 0.0_dp)

      st = trc_iterations(h, niter)
      print '(a,i0)', '  iterations = ', niter
      call check(nfail, niter >= 1, &
                 'at least one iteration ran before the failure was reported (this is TRC_ERR_NOCONV, '// &
                 'not the "never started" TRC_ERR_UNSUPPORTED path)', real(niter, dp))

      st = trc_message(h, msg_buf, int(size(msg_buf), c_int))
      block
         character(len=:), allocatable :: msg
         integer :: k
         character(len=MSG_LEN) :: tmp
         tmp = ''
         do k = 1, MSG_LEN
            if (msg_buf(k) == c_null_char) exit
            tmp(k:k) = msg_buf(k)
         end do
         msg = trim(tmp)
         print '(a,a)', '  message: ', msg
         call check(nfail, index(msg, 'diagonalis') > 0 .or. index(msg, 'syev') > 0, &
                    'trc_message names the diagonalisation, not a bare non-convergence', 0.0_dp)
      end block

      ! A SECOND, ordinary SAD run on the SAME handle must still work --
      ! the failed run must not have left the context in a state that
      ! poisons whatever runs next (trc_set_geometry's own
      ! invalidate_runs convention, extended to this failure path).
      st = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
      st = trc_run_scf(h)
      call check(nfail, st == TRC_OK, 'an ordinary SCF on the SAME handle recovers after the failure', &
                 real(st, dp))

      st = trc_destroy(h)

      if (nfail > 0) call test_failed(error, "check_scf_diag_failure: "//"check(s) failed (see the output above)")
   end subroutine test_diag_failure

   pure function c_path(s) result(cs)
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function c_path

   subroutine check(nf, cond, what, val)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in) :: val
      if (cond) then
         print '(a,a)', '  PASS  ', what
      else
         print '(a,a,a,es12.4)', '  FAIL  ', what, '   value=', val
         nf = nf + 1
      end if
   end subroutine check

end module check_scf_diag_failure_suite

program check_scf_diag_failure
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_scf_diag_failure_suite, only: collect_scf_diag_failure
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_diag_failure, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_diag_failure
