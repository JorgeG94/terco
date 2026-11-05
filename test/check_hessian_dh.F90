!
! The finite-difference Hessian (trc_run_hessian) over a correlated
! reference: a double hybrid, E = E_KS + c_os E_os + c_ss E_ss with the
! RI-MP2 on the KS orbitals. trc_run_hessian differences the analytic
! gradient trc_run_gradient returns, and with an RI-MP2 result in the
! context it reruns trc_run_rimp2 at every displaced geometry, so the
! Hessian is the double hybrid's. Before that it refused the reference
! outright ("each displaced point would need its own trc_run_rimp2").
!
! The references are terco's own double-hybrid ENERGY, never its gradient:
!
! A. Energy/gradient consistency at a non-stationary geometry (water.xyz is
!    far from the B2PLYP/def2-SVP minimum, |g| ~ 1e-2): a central
!    difference of trc_energy's total along one Cartesian direction against
!    the analytic component. This is the pairing --opt hands DL-FIND at
!    every step (trc_dlfind_bridge's trc_opt_options_t%rimp2), so it is
!    checked on its own. trc_energy's double-hybrid total is anchored to
!    pyscf by test/scf_dh.F90.
!
! B. The Hessian at the B2PLYP/def2-SVP minimum against second differences
!    of that same energy, (E(x+h) - 2 E(x) + E(x-h))/h^2 on every diagonal
!    element: no gradient code in the reference at all. Plus symmetry, the
!    reference state (and its correlated energy) restored afterwards, and a
!    sane harmonic analysis: six zero modes, three real ones, no NaN. The
!    same at the RHF + RI-MP2/def2-SVP minimum (plain --rimp2 --hessian,
!    E_RHF + E_os + E_ss).
!
! C. What the RI-MP2 gradient refuses, trc_run_hessian still refuses, with
!    the gradient's message: an unrestricted RI-MP2 reference and RI-MP2
!    over PCM.
!
! Water, def2-SVP (spherical) with def2-SVP-RIFIT, frozen 1, grid level 3.
!
! A test-drive suite.
!
module check_hessian_dh_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use, intrinsic :: iso_fortran_env, only: output_unit
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_constants, only: trc_atomic_mass, BOHR_TO_ANGSTROM
   use trc_vib, only: trc_vibrational_analysis
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_method, trc_set_convergence, trc_set_rimp2, trc_set_pcm, trc_set_geometry, &
                               trc_run_scf, trc_run_rimp2, trc_run_gradient, trc_gradient, trc_energy, &
                               trc_run_hessian, trc_hessian, trc_set_hessian_step, trc_message, TRC_OK, &
                               TRC_ERR_UNSUPPORTED
   implicit none
   private
   public :: collect_hessian_dh

   real(dp), parameter :: FD_STEP = 1.0e-3_dp
   !! Bohr, the central difference of A.
   real(dp), parameter :: FD_GRAD_TOL = 5.0e-7_dp
   !! A: measured 6.2e-8 (O z) and 1.4e-7 (H y) on the gfortran 15 host
   !! build -- the h^2 truncation of a plain central difference (check_dh_
   !! grad.F90 measures up to 2.8e-7 for cc-pVDZ at this h).
   real(dp), parameter :: HESSIAN_STEP = 5.0e-3_dp
   !! Bohr: trc_run_hessian's own default, and the second difference's.
   real(dp), parameter :: HESSIAN_DIAG_TOL = 1.0e-5_dp
   !! B: measured worst 3.7e-6 (b2plyp) and 3.8e-6 (RI-MP2) Eh/Bohr^2, both
   !! H(3,3). Both sides carry an O(h^2)
   !! truncation at h = 0.005 (check_fd_hessian_vs_pyscf.F90 measures
   !! 1.4e-5 for the gradient side alone on water/HF).
   real(dp), parameter :: ENERGY_TOL = 1.0e-9_dp
   real(dp), parameter :: MIN_ANGSTROM(3, 3) = reshape([ &
                                               0.0_dp, 0.0_dp, 0.126501316433_dp, &
                                               0.0_dp, 0.754065492362_dp, -0.473800658216_dp, &
                                               0.0_dp, -0.754065492362_dp, -0.473800658216_dp], [3, 3])
   !! The B2PLYP/def2-SVP minimum, Angstrom: terco --opt (this setup), max
   !! |g| 1.1e-5 Eh/Bohr there.
   real(dp), parameter :: MIN_RIMP2_ANGSTROM(3, 3) = reshape([ &
                                                     0.0_dp, 0.0_dp, 0.128286001881_dp, &
                                                     0.0_dp, 0.750231168460_dp, -0.474693000940_dp, &
                                                     0.0_dp, -0.750231168460_dp, -0.474693000940_dp], [3, 3])
   !! The RHF + RI-MP2/def2-SVP minimum, Angstrom: terco --rimp2 --opt
   !! --opt-gtol 1e-6, max |g| 3.2e-8 Eh/Bohr there.
   real(dp), parameter :: Z_WATER(3) = [8.0_dp, 1.0_dp, 1.0_dp]

contains

   subroutine collect_hessian_dh(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("b2plyp: a central difference of the total energy matches the analytic gradient", &
                               test_fd_gradient), &
                  new_unittest("b2plyp: the FD Hessian matches second differences of the total energy", &
                               test_hessian), &
                  new_unittest("RHF + RI-MP2: the FD Hessian matches second differences of the total energy", &
                               test_hessian_rimp2), &
                  new_unittest("the Hessian over an unrestricted RI-MP2 reference is still refused", &
                               test_refuses_unrestricted), &
                  new_unittest("the Hessian over RI-MP2 + PCM is still refused", test_refuses_pcm) &
                  ]
   end subroutine collect_hessian_dh

   subroutine water_context(error, h, functional, xyz, charge, multiplicity, pcm)
      !! A handle with the molecule, basis and method set, no SCF yet.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(out) :: h
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: xyz(3, 3)
      integer, intent(in) :: charge, multiplicity
      logical, intent(in) :: pcm
      integer(c_int) :: rc
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, 3_c_int, real(Z_WATER, c_double), reshape(real(xyz, c_double), [9]), &
                                              int(charge, c_int), int(multiplicity, c_int))
      if (rc == TRC_OK) rc = trc_set_basis_json(h, "basis_sets/def2-svp.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_aux_json(h, "basis_sets/def2-svp-rifit.json"//c_null_char)
      if (rc == TRC_OK) rc = trc_set_method(h, functional//c_null_char, 3_c_int)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 8_c_int)
      if (rc == TRC_OK) rc = trc_set_rimp2(h, 1_c_int, 0_c_int)
      if (rc == TRC_OK .and. pcm) rc = trc_set_pcm(h, "cpcm"//c_null_char, 78.3553_c_double, 0_c_int, 0.0_c_double)
      call check(error, rc == TRC_OK, "context setup failed: "//message_of(h))
   end subroutine water_context

   subroutine correlated_energy(error, h, xyz, energy)
      !! trc_energy's total (SCF + correlation) at `xyz`, Bohr.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(in) :: h
      real(dp), intent(in) :: xyz(3, 3)
      real(dp), intent(out) :: energy
      integer(c_int) :: rc
      real(c_double) :: e
      energy = 0.0_dp
      rc = trc_set_geometry(h, reshape(real(xyz, c_double), [9]))
      if (rc == TRC_OK) rc = trc_run_scf(h)
      if (rc == TRC_OK) rc = trc_run_rimp2(h)
      if (rc == TRC_OK) rc = trc_energy(h, e)
      call check(error, rc == TRC_OK, "SCF + RI-MP2 energy failed: "//message_of(h))
      energy = real(e, dp)
   end subroutine correlated_energy

   subroutine test_fd_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer :: natm, kk, k, atom, comp
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: x0(3, 3), xd(3, 3), eplus, eminus, e0, fd, worst
      real(c_double) :: gbuf(9)
      integer(c_int) :: rc
      integer, parameter :: COORDS(2) = [3, 5]
      !! O z and H(1) y: the two largest components at this geometry.
      call read_xyz("water.xyz", natm, z, r)
      x0 = r
      call water_context(error, h, "b2plyp", x0, 0, 1, .false.)
      if (allocated(error)) return
      call correlated_energy(error, h, x0, e0)
      if (allocated(error)) return
      rc = trc_run_gradient(h)
      if (rc == TRC_OK) rc = trc_gradient(h, gbuf)
      call check(error, rc == TRC_OK, "trc_run_gradient: "//message_of(h))
      if (allocated(error)) return
      call check(error, maxval(abs(gbuf)) > 1.0e-3_c_double, "the start geometry should be far from stationary")
      if (allocated(error)) return
      worst = 0.0_dp
      do kk = 1, size(COORDS)
         k = COORDS(kk)
         atom = (k - 1)/3 + 1
         comp = mod(k - 1, 3) + 1
         xd = x0
         xd(comp, atom) = x0(comp, atom) + FD_STEP
         call correlated_energy(error, h, xd, eplus)
         if (allocated(error)) return
         xd(comp, atom) = x0(comp, atom) - FD_STEP
         call correlated_energy(error, h, xd, eminus)
         if (allocated(error)) return
         fd = (eplus - eminus)/(2.0_dp*FD_STEP)
         write (output_unit, "(a,i0,a,es22.14,a,es22.14,a,es10.3)") "  coordinate ", k, ": analytic ", gbuf(k), &
            "  FD ", fd, "  diff ", abs(fd - gbuf(k))
         worst = max(worst, abs(fd - real(gbuf(k), dp)))
      end do
      call check(error, worst <= FD_GRAD_TOL, "the total energy's central difference disagrees with the analytic "// &
                 "double-hybrid gradient", more="worst "//real_text(worst))
      rc = trc_destroy(h)
   end subroutine test_fd_gradient

   subroutine test_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      call hessian_case(error, "b2plyp", MIN_ANGSTROM)
   end subroutine test_hessian

   subroutine test_hessian_rimp2(error)
      !! Plain --rimp2 --hessian: an RHF reference, E_RHF + E_os + E_ss.
      type(error_type), allocatable, intent(out) :: error
      call hessian_case(error, "", MIN_RIMP2_ANGSTROM)
   end subroutine test_hessian_rimp2

   subroutine hessian_case(error, functional, minimum_angstrom)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: functional
      real(dp), intent(in) :: minimum_angstrom(3, 3)
      type(c_ptr) :: h
      real(dp) :: x0(3, 3), xd(3, 3), e0, e_after, eplus, eminus, fd, worst, asym
      real(dp), allocatable :: hess(:, :), mass(:), freq(:), redmass(:), force_const(:), displacement(:, :, :)
      real(c_double) :: hbuf(81), e
      integer(c_int) :: rc
      integer :: k, atom, comp, n_tr, i
      logical :: is_linear
      type(error_t) :: verr
      x0 = minimum_angstrom/BOHR_TO_ANGSTROM
      call water_context(error, h, functional, x0, 0, 1, .false.)
      if (allocated(error)) return
      call correlated_energy(error, h, x0, e0)
      if (allocated(error)) return
      rc = trc_set_hessian_step(h, real(HESSIAN_STEP, c_double))
      if (rc == TRC_OK) rc = trc_run_hessian(h)
      call check(error, rc == TRC_OK, "trc_run_hessian over a correlated reference: "//message_of(h))
      if (allocated(error)) return
      rc = trc_hessian(h, hbuf)
      call check(error, rc == TRC_OK, "trc_hessian after a successful trc_run_hessian")
      if (allocated(error)) return
      hess = reshape(real(hbuf, dp), [9, 9])

      ! The reference state comes back with its RI-MP2 result, so trc_energy
      ! is the correlated total again, not the SCF energy alone.
      rc = trc_energy(h, e)
      e_after = real(e, dp)
      call check(error, rc == TRC_OK .and. abs(e_after - e0) <= ENERGY_TOL, &
                 "trc_energy after trc_run_hessian is not the reference's correlated total", &
                 more="before "//real_text(e0)//", after "//real_text(e_after))
      if (allocated(error)) return

      asym = maxval(abs(hess - transpose(hess)))
      call check(error, asym <= 1.0e-12_dp, "the Hessian is not symmetric", more=real_text(asym))
      if (allocated(error)) return

      worst = 0.0_dp
      do k = 1, 9
         atom = (k - 1)/3 + 1
         comp = mod(k - 1, 3) + 1
         xd = x0
         xd(comp, atom) = x0(comp, atom) + HESSIAN_STEP
         call correlated_energy(error, h, xd, eplus)
         if (allocated(error)) return
         xd(comp, atom) = x0(comp, atom) - HESSIAN_STEP
         call correlated_energy(error, h, xd, eminus)
         if (allocated(error)) return
         fd = (eplus - 2.0_dp*e0 + eminus)/HESSIAN_STEP**2
         write (output_unit, "(a,i0,a,es20.12,a,es20.12,a,es10.3)") "  H(", k, ",k): FD-gradient ", hess(k, k), &
            "  energy 2nd difference ", fd, "  diff ", abs(fd - hess(k, k))
         worst = max(worst, abs(fd - hess(k, k)))
      end do
      call check(error, worst <= HESSIAN_DIAG_TOL, "the Hessian's diagonal disagrees with second differences "// &
                 "of the correlated energy", more="worst "//real_text(worst))
      if (allocated(error)) return

      allocate (mass(3))
      do i = 1, 3
         mass(i) = trc_atomic_mass(nint(Z_WATER(i)), verr)
      end do
      call trc_vibrational_analysis(mass, x0, hess, freq, redmass, force_const, displacement, is_linear, n_tr, verr)
      call check(error, .not. verr%has_error(), "trc_vibrational_analysis failed")
      if (allocated(error)) return
      write (output_unit, "(a,3f12.2)") "  frequencies (cm-1): ", freq
      call check(error, n_tr == 6 .and. .not. is_linear .and. size(freq) == 3, "expected 3N-6 = 3 vibrational modes")
      if (allocated(error)) return
      call check(error, all(freq == freq) .and. all(freq > 1000.0_dp) .and. all(freq < 4500.0_dp), &
                 "a frequency is NaN, imaginary or out of range for water at its minimum")
      rc = trc_destroy(h)
   end subroutine hessian_case

   subroutine test_refuses_unrestricted(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_hessian_refusal(error, "", 1, 2, .false., "restricted (closed-shell RHF reference) only")
   end subroutine test_refuses_unrestricted

   subroutine test_refuses_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      call expect_hessian_refusal(error, "", 0, 1, .true., "no PCM response in the Z-vector")
   end subroutine test_refuses_pcm

   subroutine expect_hessian_refusal(error, functional, charge, multiplicity, pcm, expected)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: functional, expected
      integer, intent(in) :: charge, multiplicity
      logical, intent(in) :: pcm
      type(c_ptr) :: h
      real(dp) :: x0(3, 3), e
      integer(c_int) :: rc
      x0 = MIN_ANGSTROM/BOHR_TO_ANGSTROM
      call water_context(error, h, functional, x0, charge, multiplicity, pcm)
      if (allocated(error)) return
      call correlated_energy(error, h, x0, e)
      if (allocated(error)) return
      rc = trc_run_hessian(h)
      call check(error, rc == TRC_ERR_UNSUPPORTED .and. index(message_of(h), expected) > 0, &
                 "trc_run_hessian should refuse with '"//expected//"'", more="got: "//message_of(h))
      rc = trc_destroy(h)
   end subroutine expect_hessian_refusal

   function message_of(h) result(text)
      type(c_ptr), intent(in) :: h
      character(len=:), allocatable :: text
      character(kind=c_char, len=1024) :: buffer
      integer(c_int) :: rc
      integer :: k
      buffer = ""
      rc = trc_message(h, buffer, 1024_c_int)
      k = index(buffer, c_null_char)
      if (k == 0) k = len_trim(buffer) + 1
      text = buffer(1:k - 1)
   end function message_of

   function real_text(x) result(text)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      write (buffer, "(es22.14)") x
      text = trim(adjustl(buffer))
   end function real_text

end module check_hessian_dh_suite

program check_hessian_dh
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_hessian_dh_suite, only: collect_hessian_dh
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_hessian_dh, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_hessian_dh
