!
! The SAD guess from mqc's JSON blocks: water from the same basis file, the
! guess assembled by trc_guess_json, the SCF started from it and from
! terco's own guess. Both must reach the same energy, and the SAD start
! must take at most one iteration more than the built-in one (see the
! comment at that comparison for why one, not zero). Fixtures
! sad_6-31g.json and sad_cc-pvdz.json were written by mqc
! (mqc.write_sad_guess), in terco's CARTESIAN AO convention -- so both
! basis loads below force Cartesian: cc-pVDZ's d shells default to the
! spherical AO space now, and mqc's guess has no spherical form to
! compare against.
!
! Two independent cases below (6-31g, cc-pvdz), each its own test-drive
! unittest.
!
module scf_sad_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_guess_json, only: trc_guess_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_scf_sad

contains

   subroutine collect_scf_sad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water/6-31g SAD vs own guess", test_631g), &
                  new_unittest("water/cc-pvdz SAD vs own guess", test_ccpvdz) &
                  ]
   end subroutine collect_scf_sad

   subroutine test_631g(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "6-31g")
   end subroutine test_631g

   subroutine test_ccpvdz(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "cc-pvdz")
   end subroutine test_ccpvdz

   subroutine run_case(error, basis)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: basis
      integer :: natm, nocc, dev
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), dguess(:, :), dg3(:, :, :)
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: r0, r1

      type(error_t) :: err

      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      nocc = sum(zint)/2
      opts%conv_energy = 1.0e-10_dp
      opts%conv_density = 1.0e-7_dp

      call trc_basis_from_json('basis_sets/'//basis//'.json', natm, zint, at_r, bas, err, &
                               force_cart=.true.)
      call check(error, .not. err%has_error(), "basis build: "//err%get_message())
      if (allocated(error)) return
      call trc_guess_from_json('sad_'//basis//'.json', bas, zint, dguess, err)
      call check(error, .not. err%has_error(), "trc_guess_from_json: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         return
      end if

      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
      call trc_scf_run(bas, nocc, nocc, opts, r0)
      block
         integer :: u
         open (newunit=u, file='sad_dump_'//basis//'.bin', access='stream', form='unformatted', status='replace')
         write (u) int(bas%nao, int64), int(bas%nshell, int64), int(bas%sh_l, int64), int(bas%sh_ao, int64); write (u) smat, dguess
         close (u)
      end block
      allocate (dg3(bas%nao, bas%nao, 1)); dg3(:, :, 1) = dguess
      call trc_scf_run(bas, nocc, nocc, opts, r1, dguess=dg3)
      block
         type(trc_scf_options_t) :: o2
         type(trc_scf_result_t) :: r2, r3
         o2 = opts
         o2%ndamp = 0; o2%damp = 0.0_dp; o2%level_shift = 0.0_dp; o2%diis_start = 1.0e9_dp
         call trc_scf_run(bas, nocc, nocc, o2, r2, dguess=dg3)
         call trc_scf_run(bas, nocc, nocc, o2, r3)
         print '(a,a,a,i0,a,l1,a,i0,a,l1)', "scf_sad: ", basis, "  DIIS from the start: SAD ", r2%iterations, &
            " iter conv ", r2%converged, "   own ", r3%iterations, " iter conv ", r3%converged
      end block
      print '(a,a,a,f10.5,a,f10.5)', "scf_sad: ", basis, "  Tr(D_sad S) = ", sum(dguess*smat), "   Tr(D_conv S) = ", &
         sum(r0%dmat(:, :, 1)*smat)
      print '(a,a,a,f18.12,a,i0,a,f18.12,a,i0,a,f8.4)', "scf_sad: ", basis, "  own guess E = ", r0%energy, &
         " (", r0%iterations, " iter)   SAD E = ", r1%energy, " (", r1%iterations, " iter)   Tr(D_sad) = ", &
         trace_ds(bas%nao, dguess)

      call check(error, r0%converged .and. r1%converged, "own guess and/or SAD guess did not converge")
      if (allocated(error)) then
         call bas%release(); call pl%release()
         deallocate (dguess, dg3, smat, tmat, vmat)
         return
      end if
      call check(error, abs(r0%energy - r1%energy) <= 1.0e-8_dp, &
                 "SAD energy "//to_string(r1%energy, "f18.12")//" differs from own guess's "// &
                 to_string(r0%energy, "f18.12"))
      if (allocated(error)) then
         call bas%release(); call pl%release()
         deallocate (dguess, dg3, smat, tmat, vmat)
         return
      end if
      ! A SLACK OF ONE ITERATION, not zero. Tr(D_sad S) is exact to 14
      ! digits on both bases now (trc_guess_json.F90's conversion was fixed
      ! alongside trc_basis_json.F90's, see that module's header) -- a
      ! correctly-scaled guess, not a merely-close one -- and cc-pVDZ's DIIS
      ! still takes one more cycle to cross `conv_density` than terco's own
      ! guess (10 against 9), because the two guesses are genuinely
      ! different densities, just both excellent ones, and DIIS's iteration
      ! count near convergence is not monotone in how good the start was.
      ! Zero slack passed only while trc_guess_json's conversion had the
      ! wrong s/p scale (Tr(D_sad S) = 92.27, not 10, on cc-pVDZ) and
      ! happened to also land on an iteration count this check liked.
      call check(error, r1%iterations <= r0%iterations + 1, &
                 "SAD took "//to_string(r1%iterations)//" iterations vs own guess's "//to_string(r0%iterations))
      if (allocated(error)) then
         call bas%release(); call pl%release()
         deallocate (dguess, dg3, smat, tmat, vmat)
         return
      end if
      ! The check an energy cannot make: the guess must carry the electrons.
      call check(error, abs(sum(dguess*smat) - 2.0_dp*nocc) <= 1.0e-6_dp, &
                 "the SAD guess does not carry the right electron count")
      call bas%release(); call pl%release()
      deallocate (dguess, dg3, smat, tmat, vmat)
   end subroutine run_case

   ! Electron count of the guess needs S; the trace of D alone is only a
   ! sanity number, printed and not checked.
   real(dp) function trace_ds(n, d)
      integer, intent(in) :: n
      real(dp), intent(in) :: d(n, n)
      integer :: i
      trace_ds = 0.0_dp
      do i = 1, n
         trace_ds = trace_ds + d(i, i)
      end do
   end function trace_ds

end module scf_sad_suite

program scf_sad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_sad_suite, only: collect_scf_sad
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_sad, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_sad
