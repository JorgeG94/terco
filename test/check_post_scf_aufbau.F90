!
! Post-SCF methods that divide by eps_a - eps_i (RI-MP2, the CPHF/CPKS
! response) refuse a reference whose occupation is not aufbau, as a
! test-drive suite.
!
! 1. scf_settings_is_aufbau on synthetic orbital energies: aufbau, a swapped
!    HOMO/LUMO, an open shell with only one spin swapped, degenerate levels
!    within the tolerance, levels beyond `nmo` (the lifted null states of a
!    linearly dependent basis) ignored, and nothing to compare.
!
! 2. A MOM-held doubly excited water (the converged HOMO and LUMO swapped,
!    given as the guess with mom on, as check_scf_stabilise.F90 does) is
!    refused by trc_run_rimp2 and trc_run_polarizability with a message that
!    says why. The same handle without the swap still runs both.
!
! 3. The aufbau references still give the anchored numbers: RI-MP2 on water
!    cc-pVDZ (E_os, E_ss of rimp2_ccpvdz.ref, 1e-9 as scf_rimp2.F90 holds it)
!    and the CPHF polarizability on water 6-31G* (cphf_water_6-31gs.ref, the
!    tolerance check_capi_polarizability.F90 holds).
!
module check_post_scf_aufbau_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_loc, c_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_scf_settings, only: scf_settings_is_aufbau, AUFBAU_TOL
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_force_cartesian, trc_set_method, trc_set_convergence, trc_set_rimp2, &
                               trc_set_scf_stabilise, trc_set_guess, trc_run_scf, trc_run_rimp2, trc_rimp2_energy, &
                               trc_run_polarizability, trc_polarizability, trc_energy, trc_nao, trc_density, &
                               trc_mo_coefficients, trc_message, TRC_OK, TRC_ERR_UNSUPPORTED, TRC_GUESS_GIVEN
   implicit none
   private
   public :: collect_post_scf_aufbau

   real(dp), parameter :: MIN_EXCITATION = 0.3_dp
   !! The swapped determinant must lie at least this far above the ground
   !! state (check_scf_stabilise.F90's number) or it is not a different state.

contains

   subroutine collect_post_scf_aufbau(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("is_aufbau: aufbau and swapped, closed shell", test_helper_closed), &
                  new_unittest("is_aufbau: an open shell is aufbau only if both spins are", test_helper_open), &
                  new_unittest("is_aufbau: the tolerance and the nmo truncation", test_helper_tolerance_nmo), &
                  new_unittest("RI-MP2: an aufbau water gives pyscf's E_os and E_ss", test_rimp2_aufbau), &
                  new_unittest("RI-MP2: a mom-held excited water is refused", test_rimp2_refused), &
                  new_unittest("CPHF: an aufbau water gives pyscf's polarizability", test_polar_aufbau), &
                  new_unittest("CPHF: a mom-held excited water is refused", test_polar_refused) &
                  ]
   end subroutine collect_post_scf_aufbau

   subroutine test_helper_closed(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(6, 1)
      integer, parameter :: NOCC(2) = [3, 0]
      eps(:, 1) = [-2.0_dp, -1.0_dp, -0.5_dp, 0.1_dp, 0.4_dp, 0.9_dp]
      call check(error, scf_settings_is_aufbau(6, 1, 6, NOCC, eps), "an aufbau ladder must be aufbau")
      if (allocated(error)) return
      ! Occupied columns 1:3 hold the swapped determinant: the HOMO column is
      ! now the old LUMO, above the first virtual.
      eps(:, 1) = [-2.0_dp, -1.0_dp, 0.4_dp, 0.1_dp, -0.5_dp, 0.9_dp]
      call check(error, .not. scf_settings_is_aufbau(6, 1, 6, NOCC, eps), "a swapped HOMO/LUMO must not be aufbau")
   end subroutine test_helper_closed

   subroutine test_helper_open(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(5, 2)
      integer, parameter :: NOCC(2) = [3, 2]
      ! The beta LUMO below the alpha HOMO is a radical, not a swap.
      eps(:, 1) = [-2.0_dp, -1.0_dp, -0.2_dp, 0.3_dp, 0.6_dp]
      eps(:, 2) = [-2.1_dp, -1.1_dp, -0.4_dp, 0.2_dp, 0.5_dp]
      call check(error, scf_settings_is_aufbau(5, 2, 5, NOCC, eps), "an aufbau open shell must be aufbau")
      if (allocated(error)) return
      eps(:, 2) = [-2.1_dp, 0.2_dp, -0.4_dp, -1.1_dp, 0.5_dp]
      call check(error, .not. scf_settings_is_aufbau(5, 2, 5, NOCC, eps), "a swapped beta spin alone must not be aufbau")
      if (allocated(error)) return
      eps(:, 2) = [-2.1_dp, -1.1_dp, -0.4_dp, 0.2_dp, 0.5_dp]
      eps(:, 1) = [-2.0_dp, 0.3_dp, -0.2_dp, -1.0_dp, 0.6_dp]
      call check(error, .not. scf_settings_is_aufbau(5, 2, 5, NOCC, eps), "a swapped alpha spin alone must not be aufbau")
   end subroutine test_helper_open

   subroutine test_helper_tolerance_nmo(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: eps(6, 1)
      integer, parameter :: NOCC(2) = [3, 0]
      ! Degenerate frontier levels filled in either order by SCF noise.
      eps(:, 1) = [-2.0_dp, -1.0_dp, -0.5_dp + 0.5_dp*AUFBAU_TOL, -0.5_dp, 0.4_dp, 0.9_dp]
      call check(error, scf_settings_is_aufbau(6, 1, 6, NOCC, eps), "a split below the tolerance must pass")
      if (allocated(error)) return
      eps(:, 1) = [-2.0_dp, -1.0_dp, -0.5_dp + 100.0_dp*AUFBAU_TOL, -0.5_dp, 0.4_dp, 0.9_dp]
      call check(error, .not. scf_settings_is_aufbau(6, 1, 6, NOCC, eps), "a split well above the tolerance must fail")
      if (allocated(error)) return
      ! nmo = 5 of 6: the last level is a lifted null state, far below the
      ! occupied ones here on purpose; it must not count as a virtual.
      eps(:, 1) = [-2.0_dp, -1.0_dp, -0.5_dp, 0.1_dp, 0.4_dp, -50.0_dp]
      call check(error, scf_settings_is_aufbau(6, 1, 5, NOCC, eps), "levels beyond nmo must be ignored")
      if (allocated(error)) return
      ! No virtual at all (nocc >= nmo) and no electrons: nothing to compare.
      call check(error, scf_settings_is_aufbau(6, 1, 3, NOCC, eps), "no virtual orbital: aufbau")
      if (allocated(error)) return
      call check(error, scf_settings_is_aufbau(6, 1, 6, [0, 0], eps), "no occupied orbital: aufbau")
   end subroutine test_helper_tolerance_nmo

   subroutine make_handle(h, basis, aux, mom, dguess, error)
      !! A Cartesian water handle (RI-MP2 and the response are Cartesian only),
      !! with `aux` when it is not empty, and, when `mom` is set, the density
      !! `dguess` given as the starting density with mom on.
      type(c_ptr), intent(out) :: h
      character(len=*), intent(in) :: basis, aux
      logical, intent(in) :: mom
      real(dp), intent(in), target, optional :: dguess(:, :)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm
      integer(c_int) :: rc
      h = c_null_ptr
      call read_xyz('water.xyz', natm, z, at_r)
      rc = trc_create(h)
      call check(error, rc == TRC_OK, "trc_create")
      if (allocated(error)) return
      rc = trc_set_molecule(h, int(natm, c_int), real(z, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_force_cartesian(h, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, 'basis_sets/'//basis//'.json'//c_null_char)
      if (rc == TRC_OK .and. len(aux) > 0) rc = trc_set_aux_json(h, 'basis_sets/'//aux//'.json'//c_null_char)
      if (rc == TRC_OK .and. len(aux) > 0) rc = trc_set_rimp2(h, 0_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 100_c_int, 10_c_int)
      if (rc == TRC_OK .and. mom) rc = trc_set_scf_stabilise(h, 1_c_int, 0.0_c_double)
      if (rc == TRC_OK .and. present(dguess)) rc = trc_set_guess(h, TRC_GUESS_GIVEN, c_loc(dguess(1, 1)), 1_c_int)
      call check(error, rc == TRC_OK, "handle setup (status "//to_string(int(rc))//")")
      if (allocated(error)) then
         rc = trc_destroy(h)
         h = c_null_ptr
      end if
   end subroutine make_handle

   subroutine swapped_density(basis, aux, dex, e0, error)
      !! The converged ground-state water in `basis` and the density of its
      !! determinant with the HOMO and LUMO swapped.
      character(len=*), intent(in) :: basis, aux
      real(dp), allocatable, intent(out) :: dex(:, :)
      real(dp), intent(out) :: e0
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: cmo(:, :)
      integer, allocatable :: cols(:)
      integer(c_int) :: rc, nao
      real(c_double) :: e
      integer :: k
      e0 = 0.0_dp
      call make_handle(h, basis, aux, .false., error=error)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "ground state trc_run_scf")
      if (.not. allocated(error)) then
         rc = trc_nao(h, nao)
         allocate (cmo(nao, nao), dex(nao, nao))
         rc = trc_mo_coefficients(h, cmo)
         rc = trc_energy(h, e)
         e0 = e
         cols = [(k, k=1, 4), 6]
         dex = 2.0_dp*matmul(cmo(:, cols), transpose(cmo(:, cols)))
      end if
      rc = trc_destroy(h)
   end subroutine swapped_density

   function message_of(h) result(text)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char) :: buf(1024)
      integer(c_int) :: rc
      integer :: i
      buf = c_null_char
      rc = trc_message(h, buf, 1024_c_int)
      text = ""
      do i = 1, 1024
         if (buf(i) == c_null_char) exit
         text = text//buf(i)
      end do
   end function message_of

   subroutine excited_handle(h, basis, aux, error)
      !! A converged MOM-held excited water handle in `basis`.
      type(c_ptr), intent(out) :: h
      character(len=*), intent(in) :: basis, aux
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable, target :: dex(:, :)
      real(dp) :: e0
      real(c_double) :: e
      integer(c_int) :: rc
      call swapped_density(basis, aux, dex, e0, error)
      if (allocated(error)) return
      call make_handle(h, basis, aux, .true., dex, error)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "mom trc_run_scf from the swapped density (status "//to_string(int(rc))//")")
      if (.not. allocated(error)) then
         rc = trc_energy(h, e)
         call check(error, real(e, dp) > e0 + MIN_EXCITATION, &
                    "the mom run did not hold the excited occupation: E - E0 = "//to_string(real(e, dp) - e0, "f12.8"))
      end if
      if (allocated(error)) then
         rc = trc_destroy(h)
         h = c_null_ptr
      end if
   end subroutine excited_handle

   subroutine test_rimp2_aufbau(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: e_os, e_ss
      real(dp) :: e_os_ref, e_ss_ref, e_scf_ref
      real(c_double) :: e_scf
      integer :: u
      call make_handle(h, 'cc-pvdz', 'cc-pvdz-rifit', .false., error=error)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_energy(h, e_scf)
      if (rc == TRC_OK) rc = trc_run_rimp2(h)
      call check(error, rc == TRC_OK, "SCF and trc_run_rimp2 on an aufbau reference (status "//to_string(int(rc))//")")
      if (.not. allocated(error)) then
         rc = trc_rimp2_energy(h, e_os, e_ss)
         call ref_require('rimp2_ccpvdz.ref', u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_real(u, 'E_OS', e_os_ref)
         call ref_read_real(u, 'E_SS', e_ss_ref)
         call ref_close(u)
         call check(error, abs(real(e_scf, dp) - e_scf_ref) < 1.0e-9_dp, "E_RHF vs pyscf, |dE| = "// &
                    to_string(abs(real(e_scf, dp) - e_scf_ref), "es10.2"))
      end if
      if (.not. allocated(error)) then
         call check(error, abs(real(e_os, dp) - e_os_ref) < 1.0e-9_dp, "E_os vs pyscf, |dE| = "// &
                    to_string(abs(real(e_os, dp) - e_os_ref), "es10.2"))
      end if
      if (.not. allocated(error)) call check(error, abs(real(e_ss, dp) - e_ss_ref) < 1.0e-9_dp, "E_ss vs pyscf, |dE| = "// &
                                             to_string(abs(real(e_ss, dp) - e_ss_ref), "es10.2"))
      rc = trc_destroy(h)
   end subroutine test_rimp2_aufbau

   subroutine test_rimp2_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      call excited_handle(h, 'cc-pvdz', 'cc-pvdz-rifit', error)
      if (allocated(error)) return
      rc = trc_run_rimp2(h)
      call check(error, rc == TRC_ERR_UNSUPPORTED, "RI-MP2 on a non-aufbau reference must be refused (status "// &
                 to_string(int(rc))//")")
      if (.not. allocated(error)) call check(error, index(message_of(h), "needs an aufbau reference") > 0, &
                                             "the refusal should say why: "//message_of(h))
      rc = trc_destroy(h)
   end subroutine test_rimp2_refused

   subroutine test_polar_aufbau(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      real(c_double) :: alpha(3, 3)
      real(dp) :: alpha_ref(3, 3), e_scf_ref
      integer :: u
      call make_handle(h, '6-31gs', '', .false., error=error)
      if (allocated(error)) return
      rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_run_polarizability(h)
      if (rc == TRC_OK) rc = trc_polarizability(h, alpha)
      call check(error, rc == TRC_OK, "SCF and trc_run_polarizability on an aufbau reference (status "// &
                 to_string(int(rc))//")")
      if (.not. allocated(error)) then
         call ref_require('cphf_water_6-31gs.ref', u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
         call ref_close(u)
         call check(error, all(abs(real(alpha, dp) - alpha_ref) < 1.0e-6_dp), "alpha vs pyscf, worst |dalpha| = "// &
                    to_string(maxval(abs(real(alpha, dp) - alpha_ref)), "es10.2"))
      end if
      rc = trc_destroy(h)
   end subroutine test_polar_aufbau

   subroutine test_polar_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      call excited_handle(h, '6-31gs', '', error)
      if (allocated(error)) return
      rc = trc_run_polarizability(h)
      call check(error, rc == TRC_ERR_UNSUPPORTED, "the response on a non-aufbau reference must be refused (status "// &
                 to_string(int(rc))//")")
      if (.not. allocated(error)) call check(error, index(message_of(h), "needs an aufbau reference") > 0, &
                                             "the refusal should say why: "//message_of(h))
      rc = trc_destroy(h)
   end subroutine test_polar_refused

end module check_post_scf_aufbau_suite

program check_post_scf_aufbau
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_post_scf_aufbau_suite, only: collect_post_scf_aufbau
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_post_scf_aufbau, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_post_scf_aufbau
