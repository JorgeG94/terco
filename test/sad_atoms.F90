module sad_atoms_suite
   !! terco's own SAD guess: electron count, and that it matches GWH on water.
   !!
   !! FORCED CARTESIAN: `smat`/`tmat`/`vmat` below come from `trc_1e` directly,
   !! which is always Cartesian-sized (b%nao), and are compared against
   !! `dguess` and used to build `dg3`. A spherical cc-pVDZ's `dguess` (from
   !! trc_sad_build) would come back nao_sph-sized instead, a shape mismatch
   !! against the Cartesian matrices this test builds by hand -- this test is
   !! about the SAD algorithm (electron count, matching GWH), not about which
   !! AO space it runs in, so it forces Cartesian rather than growing a
   !! parallel spherical path for a one-off check.
   !!
   !! Two independent cases below (6-31g, cc-pvdz), each its own test-drive
   !! unittest.
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sad_atoms

contains

   subroutine collect_sad_atoms(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water/6-31g SAD vs GWH", test_631g), &
                  new_unittest("water/cc-pvdz SAD vs GWH", test_ccpvdz) &
                  ]
   end subroutine collect_sad_atoms

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
      real(dp), allocatable :: at_r(:, :), dguess(:, :), dg3(:, :, :), smat(:, :), tmat(:, :), vmat(:, :)
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: r_gwh, r_sad
      type(error_t) :: err
      real(dp) :: trace

      dev = trc_bind_device(0)
      call read_xyz('water.xyz', natm, zint, at_r)
      nocc = sum(zint)/2
      opts%conv_energy = 1.0e-10_dp

      call trc_basis_from_json('basis_sets/'//basis//'.json', natm, zint, at_r, bas, err, &
                               force_cart=.true.)
      call check(error, .not. err%has_error(), "basis build: "//err%get_message())
      if (allocated(error)) return

      call trc_sad_build(bas, dguess, err, verbose=.true.)
      call check(error, .not. err%has_error(), "trc_sad_build: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         return
      end if

      call bas%to_device()
      call pl%build(bas, 1.0e-12_dp); call pl%to_device()
      allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
      call trc_1e(bas, pl, smat, tmat, vmat)
      trace = sum(dguess*smat)
      call trc_scf_run(bas, nocc, nocc, opts, r_gwh)
      allocate (dg3(bas%nao, bas%nao, 1)); dg3(:, :, 1) = dguess
      call trc_scf_run(bas, nocc, nocc, opts, r_sad, dguess=dg3)
      print '(a,a,a,f10.6,a,f18.12,a,i0,a,f18.12,a,i0,a)', "sad_atoms: ", basis, "  Tr(D S) = ", trace, &
         "   GWH E = ", r_gwh%energy, " (", r_gwh%iterations, ")   SAD E = ", r_sad%energy, " (", r_sad%iterations, ")"

      deallocate (smat, tmat, vmat)
      call pl%release(); call bas%release()

      call check(error, abs(trace - real(sum(zint), dp)) <= 1.0e-8_dp, &
                 basis//": Tr(D S) = "//to_string(trace, "f18.12")//", want "//to_string(real(sum(zint), dp), "f18.12"))
      if (allocated(error)) then
         deallocate (dg3, dguess)
         return
      end if
      call check(error, abs(r_sad%energy - r_gwh%energy) <= 1.0e-8_dp, &
                 basis//": SAD energy "//to_string(r_sad%energy, "f18.12")// &
                 " differs from GWH's "//to_string(r_gwh%energy, "f18.12"))
      if (allocated(error)) then
         deallocate (dg3, dguess)
         return
      end if
      ! With DIIS from the first iteration water converges in ten from either
      ! guess; the count only has to not be worse by more than a couple.
      call check(error, r_sad%converged .and. r_sad%iterations <= r_gwh%iterations + 2, &
                 basis//": SAD took "//to_string(r_sad%iterations)// &
                 " iterations vs GWH's "//to_string(r_gwh%iterations))
      deallocate (dg3, dguess)
   end subroutine run_case

end module sad_atoms_suite

program sad_atoms
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use sad_atoms_suite, only: collect_sad_atoms
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sad_atoms, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program sad_atoms
