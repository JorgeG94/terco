!
! The SAD guess from mqc's JSON blocks: water from the same basis file, the
! guess assembled by trc_guess_json, the SCF started from it and from
! terco's own guess. Both must reach the same energy, and the SAD start
! must take no more iterations than the built-in one. Fixtures sad_6-31g.json
! and sad_cc-pvdz.json were written by mqc (mqc.write_sad_guess).
!
program scf_sad
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_guess_json, only: trc_guess_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   character(len=*), parameter :: bases(2) = [character(len=24) :: "6-31g", "cc-pvdz"]
   integer :: natm, ib, nocc, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), dguess(:, :), dg3(:, :, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: r0, r1
   type(error_t) :: err
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nocc = sum(zint)/2
   opts%conv_energy = 1.0e-10_dp
   opts%conv_density = 1.0e-7_dp
   ok = .true.
   do ib = 1, 2
      call trc_basis_from_json('basis_sets/'//trim(bases(ib))//'.json', natm, zint, at_r, bas, err)
      if (.not. err%has_error()) call trc_guess_from_json('sad_'//trim(bases(ib))//'.json', bas, zint, dguess, err)
      if (err%has_error()) then
         print '(a)', "scf_sad: "//err%get_message()
         stop 1
      end if
      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
      call trc_scf_run(bas, nocc, nocc, opts, r0)
      block
         integer :: u
         open (newunit=u, file='sad_dump_'//trim(bases(ib))//'.bin', access='stream', form='unformatted', status='replace')
         write (u) int(bas%nao, 8), int(bas%nshell, 8), int(bas%sh_l, 8), int(bas%sh_ao, 8); write (u) smat, dguess
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
         print '(a,a8,a,i0,a,l1,a,i0,a,l1)', "scf_sad: ", bases(ib), "  DIIS from the start: SAD ", r2%iterations, &
            " iter conv ", r2%converged, "   own ", r3%iterations, " iter conv ", r3%converged
      end block
      print '(a,a8,a,f10.5,a,f10.5)', "scf_sad: ", bases(ib), "  Tr(D_sad S) = ", sum(dguess*smat), "   Tr(D_conv S) = ", &
         sum(r0%dmat(:, :, 1)*smat)
      print '(a,a8,a,f18.12,a,i0,a,f18.12,a,i0,a,f8.4)', "scf_sad: ", bases(ib), "  own guess E = ", r0%energy, &
         " (", r0%iterations, " iter)   SAD E = ", r1%energy, " (", r1%iterations, " iter)   Tr(D_sad) = ", &
         trace_ds(bas%nao, dguess)
      if (.not. (r0%converged .and. r1%converged)) ok = .false.
      if (abs(r0%energy - r1%energy) > 1.0e-8_dp) ok = .false.
      if (r1%iterations > r0%iterations) ok = .false.
      call bas%release()
      deallocate (dguess, dg3, smat, tmat, vmat); call pl%release()
   end do
   print '(a)', merge("scf_sad: PASS", "scf_sad: FAIL", ok)
   if (.not. ok) stop 1

contains

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

end program scf_sad
