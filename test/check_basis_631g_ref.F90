!
! build_631g against the JSON path, the ONE place in the tree that is
! allowed to call build_631g for its own basis (see that subroutine's own
! comment in trc_test_basis.F90 for why). build_631g is a second, literal
! transcription of the same BSE 6-31G numbers basis_sets/6-31g.json holds,
! typed by hand to 7 significant figures; trc_basis_from_json reads the
! JSON file's own full-precision strings. Agreement between the two, on
! water through a twenty-water cluster, is the external check that the
! JSON path (one shell per coefficient column, each normalised and then
! divided by common_fac_sp so trc_basis_t%build's own fold lands on the
! right number -- see trc_basis_json.F90's header) builds the basis
! build_631g was always known to build correctly.
!
! THE NUMBER THIS CHECK WOULD HAVE CAUGHT: before trc_basis_json's
! `normalise` was corrected to divide by common_fac_sp, an s or p shell's
! coefficients came out short by exactly that factor, invisible on one
! water (the energies agreed to 4e-8 Eh) and breaking the SCF outright on
! twenty (off by hundreds of Hartree) -- see trc_basis_json.F90's header
! for the mechanism. A per-shell coefficient scale is not something an
! energy check at ONE geometry reliably catches (the two paths can still
! agree to several digits at small sizes even when wrong), which is why
! this check runs the comparison at more than one cluster size rather than
! trusting a single geometry's agreement.
!
module check_basis_631g_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_test_basis, only: read_xyz, build_631g
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_basis_631g_tests

   real(dp), parameter :: TOL_C = 2.0e-4_dp
      !! Per-primitive coefficient agreement. Does not grow with cluster size
      !! (it is a per-shell comparison, not an energy), so one tolerance
      !! serves every geometry below.

contains

   subroutine collect_basis_631g_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water.xyz: build_631g vs the JSON path", test_water), &
                  new_unittest("w3.xyz: build_631g vs the JSON path", test_w3), &
                  new_unittest("w5.xyz: build_631g vs the JSON path", test_w5) &
                  ]
   end subroutine collect_basis_631g_tests

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_one("water.xyz", error)
   end subroutine test_water

   subroutine test_w3(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_one("w3.xyz", error)
   end subroutine test_w3

   subroutine test_w5(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_one("w5.xyz", error)
   end subroutine test_w5

   subroutine compare_one(xyzfile, error)
      character(len=*), intent(in) :: xyzfile
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b1, b2
      type(error_t) :: err
      integer, allocatable :: zint(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: at_r(:, :), at_z(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer :: natm, nsh, maxnp, i, k, nmismatch
      real(dp) :: maxdc, maxde

      call read_xyz(xyzfile, natm, zint, at_r)
      allocate (at_z(natm)); at_z = real(zint, dp)
      call build_631g(natm, zint, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      call b1%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
      call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, b2, err)
      call check(error,.not. err%has_error(), xyzfile//": "//err%get_message())
      if (allocated(error)) return

      call check(error, b1%nshell == b2%nshell .and. b1%nao == b2%nao, &
                 xyzfile//": shell/AO count mismatch")
      if (allocated(error)) then
         call b1%release(); call b2%release()
         return
      end if

      nmismatch = 0
      maxdc = 0.0_dp; maxde = 0.0_dp
      do i = 1, b1%nshell
         if (b1%sh_l(i) /= b2%sh_l(i) .or. b1%sh_np(i) /= b2%sh_np(i)) then
            nmismatch = nmismatch + 1
            cycle
         end if
         do k = 1, b1%sh_np(i)
            maxde = max(maxde, abs(b1%sh_e(k, i) - b2%sh_e(k, i)))
            maxdc = max(maxdc, abs(b1%sh_c(k, i) - b2%sh_c(k, i)))
         end do
      end do
      print '(a,es10.2,a,es10.2)', "check_basis_631g_ref: "//xyzfile//"  max|de| ", maxde, &
         "   max|dc| ", maxdc

      call check(error, nmismatch == 0, xyzfile//": "//to_string(nmismatch)// &
                 " shell(s) with a shape mismatch")
      if (allocated(error)) then
         call b1%release(); call b2%release()
         return
      end if

      ! de is the hand-copied literal's own precision (BSE exponents given
      ! to 7 figures): a few times 1e-8 is the literal, not a bug. dc
      ! follows from de through the normalisation, so TOL_C is generous
      ! enough to pass the literal's rounding and tight enough that the
      ! common_fac_sp regression above (a fractional-constant-sized error,
      ! 0.28 or 0.49 relative) fails it outright.
      call check(error, maxde < 1.0e-6_dp, xyzfile//": exponents differ, max|de| = "// &
                 to_string(maxde, "es10.2"))
      if (allocated(error)) then
         call b1%release(); call b2%release()
         return
      end if
      call check(error, maxdc < TOL_C, xyzfile//": coefficients differ, max|dc| = "// &
                 to_string(maxdc, "es10.2"))
      call b1%release(); call b2%release()
   end subroutine compare_one

end module check_basis_631g_suite

program check_basis_631g_ref
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use trc_api, only: trc_bind_device
   use check_basis_631g_suite, only: collect_basis_631g_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, dev

   dev = trc_bind_device(0)
   stat = 0
   testsuites = [new_testsuite("check_basis_631g_ref", collect_basis_631g_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_basis_631g_ref
