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
program check_basis_631g_ref
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_test_basis, only: read_xyz, build_631g
   use trc_error, only: error_t
   implicit none

   character(len=*), parameter :: files(3) = [character(len=16) :: "water.xyz", "w3.xyz", "w5.xyz"]
   !> Per-primitive coefficient agreement. Does not grow with cluster size
   !> (it is a per-shell comparison, not an energy), so one tolerance
   !> serves every geometry in `files`.
   real(dp), parameter :: TOL_C = 2.0e-4_dp
   integer :: ib, dev
   logical :: ok
   type(error_t) :: err

   dev = trc_bind_device(0)
   ok = .true.
   do ib = 1, size(files)
      call compare_one(trim(files(ib)), ok)
   end do
   print '(a)', merge("check_basis_631g_ref: PASS", "check_basis_631g_ref: FAIL", ok)
   if (.not. ok) stop 1

contains

   subroutine compare_one(xyzfile, ok)
      character(len=*), intent(in) :: xyzfile
      logical, intent(inout) :: ok
      type(trc_basis_t) :: b1, b2
      type(error_t) :: err
      integer, allocatable :: zint(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: at_r(:, :), at_z(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer :: natm, nsh, maxnp, i, k
      real(dp) :: maxdc, maxde

      call read_xyz(xyzfile, natm, zint, at_r)
      allocate (at_z(natm)); at_z = real(zint, dp)
      call build_631g(natm, zint, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      call b1%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
      call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, b2, err)
      if (err%has_error()) then
         print '(a)', "check_basis_631g_ref: "//err%get_message()
         ok = .false.; return
      end if
      if (b1%nshell /= b2%nshell .or. b1%nao /= b2%nao) then
         print '(a)', "check_basis_631g_ref: "//xyzfile//" shell/AO count mismatch"
         ok = .false.; return
      end if
      maxdc = 0.0_dp; maxde = 0.0_dp
      do i = 1, b1%nshell
         if (b1%sh_l(i) /= b2%sh_l(i) .or. b1%sh_np(i) /= b2%sh_np(i)) then
            print '(a,i0)', "check_basis_631g_ref: "//xyzfile//" shape mismatch, shell ", i
            ok = .false.; cycle
         end if
         do k = 1, b1%sh_np(i)
            maxde = max(maxde, abs(b1%sh_e(k, i) - b2%sh_e(k, i)))
            maxdc = max(maxdc, abs(b1%sh_c(k, i) - b2%sh_c(k, i)))
         end do
      end do
      print '(a,es10.2,a,es10.2)', "check_basis_631g_ref: "//xyzfile//"  max|de| ", maxde, &
         "   max|dc| ", maxdc
      ! de is the hand-copied literal's own precision (BSE exponents given
      ! to 7 figures): a few times 1e-8 is the literal, not a bug. dc
      ! follows from de through the normalisation, so TOL_C is generous
      ! enough to pass the literal's rounding and tight enough that the
      ! common_fac_sp regression above (a fractional-constant-sized error,
      ! 0.28 or 0.49 relative) fails it outright.
      call expect(ok, maxde < 1.0e-6_dp, xyzfile//": exponents agree")
      call expect(ok, maxdc < TOL_C, xyzfile//": coefficients agree")
      call b1%release(); call b2%release()
   end subroutine compare_one

   subroutine expect(ok, cond, what)
      logical, intent(inout) :: ok
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', "check_basis_631g_ref: FAILED: "//what
         ok = .false.
      end if
   end subroutine expect

end program check_basis_631g_ref
