!
! A malformed SCF guess density -- one entry replaced by a signalling NaN,
! simulating the warm-start density app/trc_dlfind_bridge.F90's own
! --opt-neb per-image buffer fed to trc_set_guess(TRC_GUESS_GIVEN) once it
! was read PAST the end of an undersized buffer (found and fixed on
! feat/opt-neb: the array was allocated at size 1 instead of nao*nao*nspin,
! so reading it back fed trc_set_guess whatever memory happened to follow
! it) -- must come back through trc_run_scf's own status/message, not
! `error stop` inside src/core/trc_linalg.F90's eigensolver.
!
! Root cause, read directly off the source (not guessed from the symptom):
! src/scf/trc_scf_driver.F90's per-iteration Fock diagonalisation
! (orbitals_and_density, F' = X^T F X then linalg_provider%syev(fp,
! eps_s)) is reached from a GIVEN density on the very first iteration --
! `dguess` skips the guess-construction path (guess_from) entirely and
! goes straight into the main loop -- so a malformed GIVEN density is the
! most direct, deterministic way to reach exactly the call site the
! original crash came from. A QUIET NaN (the IEEE default) propagates
! through LAPACK's dsyev with info=0 on this build (measured separately,
! test/check_linalg_mat.F90%test_syev_failure's own header); a SIGNALLING
! one does not, which is why this uses ieee_signaling_nan specifically.
!
program check_scf_diag_failure
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_loc, c_null_ptr, c_null_char
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan
   use pic_types, only: dp
   use trc_capi, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_method, &
                       trc_set_convergence, trc_set_guess, trc_run_scf, trc_nao, trc_converged, &
                       trc_iterations, trc_message, TRC_OK, TRC_ERR_NOCONV, TRC_GUESS_GIVEN, TRC_GUESS_SAD
   use trc_test_basis, only: read_xyz
   implicit none

   integer, parameter :: MSG_LEN = 1024
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

   ! A normal, converged SCF first -- so the density this test corrupts is
   ! otherwise exactly what a real warm start would have handed back.
   st = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
   st = trc_run_scf(h)
   call report(nfail, st == TRC_OK, 'check_scf_diag_failure: the baseline SAD SCF converges', real(st, dp))

   st = trc_nao(h, nao)
   call report(nfail, st == TRC_OK, 'check_scf_diag_failure: trc_nao', real(st, dp))
   allocate (dbuf(int(nao)*int(nao)))
   dbuf = 0.0_c_double
   ! A plausible-looking density: occupied-ish diagonal, nothing exotic --
   ! the ONE thing wrong with it is a single signalling-NaN off-diagonal
   ! entry (symmetric, as a density must be).
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
   call report(nfail, st == TRC_OK, 'check_scf_diag_failure: trc_set_guess(GIVEN) itself accepts the buffer', &
              real(st, dp))

   st = trc_run_scf(h)
   print '(a,i0)', '  trc_run_scf status = ', st
   call report(nfail, st == TRC_ERR_NOCONV, &
              'check_scf_diag_failure: a signalling-NaN guess density comes back as TRC_ERR_NOCONV, '// &
              'not a crash', real(st, dp))

   st = trc_converged(h, conv)
   call report(nfail, conv == 0, 'check_scf_diag_failure: trc_converged reports false')

   st = trc_iterations(h, niter)
   print '(a,i0)', '  iterations = ', niter
   call report(nfail, niter >= 1, &
              'check_scf_diag_failure: at least one iteration ran before the failure was reported '// &
              '(this is TRC_ERR_NOCONV, not the "never started" TRC_ERR_UNSUPPORTED path)', real(niter, dp))

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
      call report(nfail, index(msg, 'diagonalis') > 0 .or. index(msg, 'syev') > 0, &
                 'check_scf_diag_failure: trc_message names the diagonalisation, not a bare non-convergence')
   end block

   ! A SECOND, ordinary SAD run on the SAME handle must still work --
   ! the failed run must not have left the context in a state that
   ! poisons whatever runs next (trc_set_geometry's own invalidate_runs
   ! convention, extended to this failure path).
   st = trc_set_guess(h, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
   st = trc_run_scf(h)
   call report(nfail, st == TRC_OK, &
              'check_scf_diag_failure: an ordinary SCF on the SAME handle recovers after the failure', &
              real(st, dp))

   st = trc_destroy(h)

   if (nfail > 0) then
      print '(a,i0,a)', 'check_scf_diag_failure: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_scf_diag_failure: RESULT: PASS'

contains

   pure function c_path(s) result(cs)
      character(len=*), intent(in) :: s
      character(kind=c_char, len=len_trim(s) + 1) :: cs
      cs = trim(s)//c_null_char
   end function c_path

   subroutine report(nf, cond, what, val)
      integer, intent(inout) :: nf
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      real(dp), intent(in), optional :: val
      if (cond) then
         print '(a,a)', '  PASS  ', what
      else
         if (present(val)) then
            print '(a,a,a,es12.4)', '  FAIL  ', what, '   value=', val
         else
            print '(a,a)', '  FAIL  ', what
         end if
         nf = nf + 1
      end if
   end subroutine report

end program check_scf_diag_failure
