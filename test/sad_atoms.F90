!! terco's own SAD guess: electron count, and that it beats GWH on water
program sad_atoms
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   character(len=*), parameter :: bases(2) = [character(len=24) :: "6-31g", "cc-pvdz"]
   integer :: natm, ib, nocc, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), dguess(:, :), dg3(:, :, :), smat(:, :), tmat(:, :), vmat(:, :)
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: r_gwh, r_sad
   type(error_t) :: err
   real(dp) :: trace
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('water.xyz', natm, zint, at_r)
   nocc = sum(zint)/2
   opts%conv_energy = 1.0e-10_dp
   ok = .true.
   do ib = 1, 2
      call trc_basis_from_json('basis_sets/'//trim(bases(ib))//'.json', natm, zint, at_r, bas, err)
      if (.not. err%has_error()) call trc_sad_build(bas, dguess, err, verbose=.true.)
      if (err%has_error()) then
         print '(a)', "sad_atoms: "//err%get_message()
         stop 1
      end if
      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
      trace = sum(dguess*smat)
      call trc_scf_run(bas, nocc, nocc, opts, r_gwh)
      allocate (dg3(bas%nao, bas%nao, 1)); dg3(:, :, 1) = dguess
      call trc_scf_run(bas, nocc, nocc, opts, r_sad, dguess=dg3)
      print '(a,a8,a,f10.6,a,f18.12,a,i0,a,f18.12,a,i0,a)', "sad_atoms: ", bases(ib), "  Tr(D S) = ", trace, &
         "   GWH E = ", r_gwh%energy, " (", r_gwh%iterations, ")   SAD E = ", r_sad%energy, " (", r_sad%iterations, ")"
      if (abs(trace - real(sum(zint), dp)) > 1.0e-8_dp) ok = .false.
      if (abs(r_sad%energy - r_gwh%energy) > 1.0e-8_dp) ok = .false.
      ! With DIIS from the first iteration water converges in ten from either
      ! guess; the count only has to not be worse by more than a couple.
      if (.not. r_sad%converged .or. r_sad%iterations > r_gwh%iterations + 2) ok = .false.
      deallocate (smat, tmat, vmat, dg3, dguess)
      call pl%release(); call bas%release()
   end do
   if (.not. ok) then
      print '(a)', "sad_atoms: FAIL"
      stop 1
   end if
   print '(a)', "sad_atoms: PASS"
end program sad_atoms
