!
! The analytic RI-MP2 nuclear gradient (src/scf/trc_rimp2_grad.F90, exact-
! ERI RHF reference + fitted correlation) against two references:
!
! A. metalquicha's RI-MP2 gradient references, on the same molecule, basis,
!    auxiliary set and frozen core, read off its validation suite:
!    metalquicha validation/validation_tests_cpu.json (commit f132dbb877),
!    entries "RI-MP2 gradient H2O cc-pvdz/cc-pvdz-rifit frozen 0", "... H2O
!    6-31g/cc-pvdz-rifit frozen 1" and "... NH3 def2-svp/def2-svp-rifit
!    frozen 0", geometries validation/inputs/sample_inputs/w1.xyz and
!    nh3.xyz (copied here as water_mqc.xyz and nh3_mqc.xyz). Provenance of
!    those numbers: tools/cpu_validation/gen_cpu_validation.py's
!    mqc_ri_mp2_gradient, i.e. tools/cpu_validation/ri_mp2_gradient.py, an
!    analytic RI-MP2 gradient in numpy over PySCF 2.14 integrals (exact RHF
!    reference, exact Z-vector operator, the same BSE JSON files read by
!    PySCF), itself validated against finite differences of PySCF's own
!    RHF + DF-MP2 energy (PySCF implements no RI-MP2 gradient:
!    DFMP2.nuc_grad_method raises NotImplementedError). The orbital and
!    auxiliary JSON files are byte-identical to metalquicha's basis_sets/
!    (def2-svp-rifit.json copied from there); spherical throughout, as
!    PySCF built them (`mol.cart = False`, auxmol inheriting it). The total
!    energy is checked against metalquicha's as well, the proof that the
!    two codes differentiate the same number.
!
! B. Central finite differences (step 1e-4 Bohr) of terco's OWN RHF +
!    RI-MP2 energy (trc_scf_run + trc_rimp2_run, the pair anchored to pyscf
!    by rimp2*.ref and sph_corr.ref), on cases metalquicha has none for:
!    Cartesian bases, formamide, frozen core on a molecule with three
!    cores, and a spherical orbital basis with its padded virtuals.
!
! Every case also checks the correlation energy the gradient's amplitudes
! give against trc_rimp2_run's (same B, same pairs) and translational
! invariance -- necessary, not sufficient.
!
! A test-drive suite.
!
module check_rimp2_grad_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t, trc_rimp2_options_t, rimp2_driver_mo_to_cart
   use trc_rimp2_grad, only: trc_rimp2_gradient, trc_rimp2_grad_options_t, trc_rimp2_grad_result_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_rimp2_grad

   real(dp), parameter :: MQC_GRAD_TOL = 5.0e-9_dp
   !! Analytic against analytic, two codes, same basis files and settings.
   real(dp), parameter :: MQC_ENERGY_TOL = 1.0e-8_dp
   real(dp), parameter :: FD_STEP = 1.0e-4_dp
   real(dp), parameter :: FD_TOL = 5.0e-8_dp
   real(dp), parameter :: ECORR_TOL = 1.0e-11_dp
   !! The gradient's own amplitudes against trc_rimp2_run's pair energies.
   real(dp), parameter :: TRANS_TOL = 1.0e-10_dp

   type :: case_t
      character(len=64) :: name = ""
      character(len=64) :: xyz = ""
      character(len=64) :: basis = ""
      character(len=64) :: aux = ""
      integer :: nfrozen = 0
      logical :: cart = .false.
   end type case_t

contains

   subroutine collect_rimp2_grad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water cc-pVDZ/cc-pVDZ-RIFIT (sph) vs metalquicha", test_mqc_water_ccpvdz), &
                  new_unittest("water 6-31G/cc-pVDZ-RIFIT, frozen 1, vs metalquicha", test_mqc_water_631g_fc1), &
                  new_unittest("NH3 def2-SVP/def2-SVP-RIFIT (sph) vs metalquicha", test_mqc_nh3_def2svp), &
                  new_unittest("formamide cc-pVDZ (Cartesian) vs finite differences", test_fd_formamide_cart), &
                  new_unittest("formamide cc-pVDZ (sph), frozen 3, vs finite differences", test_fd_formamide_sph_fc3), &
                  new_unittest("water cc-pVDZ (Cartesian), frozen 1, vs finite differences", test_fd_water_cart_fc1) &
                  ]
   end subroutine collect_rimp2_grad

   ! ---- A. metalquicha -----------------------------------------------------

   subroutine test_mqc_water_ccpvdz(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cc-pVDZ fc0"
      c%xyz = "water_mqc.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      call against_mqc(error, c, -76.230257623889_dp, reshape([ &
                       -0.0_dp, 2.31e-10_dp, -0.015717146972_dp, &
                       0.0_dp, 0.002913907061_dp, 0.007858573573_dp, &
                       -0.0_dp, -0.002913907293_dp, 0.0078585734_dp], [3, 3]))
   end subroutine test_mqc_water_ccpvdz

   subroutine test_mqc_water_631g_fc1(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water 6-31G fc1"
      c%xyz = "water_mqc.xyz"
      c%basis = "basis_sets/6-31g.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%nfrozen = 1
      call against_mqc(error, c, -76.112532751837_dp, reshape([ &
                       0.0_dp, 2.3e-10_dp, -0.011977548613_dp, &
                       0.0_dp, -0.016321577203_dp, 0.005988774392_dp, &
                       -0.0_dp, 0.016321576973_dp, 0.005988774222_dp], [3, 3]))
   end subroutine test_mqc_water_631g_fc1

   subroutine test_mqc_nh3_def2svp(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "NH3 def2-SVP fc0"
      c%xyz = "nh3_mqc.xyz"
      c%basis = "basis_sets/def2-svp.json"
      c%aux = "basis_sets/def2-svp-rifit.json"
      call against_mqc(error, c, -56.338085238337_dp, reshape([ &
                       -0.0_dp, -0.0_dp, 0.008308503754_dp, &
                       -0.002395148441_dp, 0.0_dp, -0.002769501251_dp, &
                       0.001197574221_dp, -0.002074259396_dp, -0.002769501251_dp, &
                       0.001197574221_dp, 0.002074259396_dp, -0.002769501251_dp], [3, 4]))
   end subroutine test_mqc_nh3_def2svp

   subroutine against_mqc(error, c, e_ref, g_ref)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      real(dp), intent(in) :: e_ref, g_ref(:, :)
      real(dp), allocatable :: grad(:, :), r(:, :)
      integer, allocatable :: z(:)
      integer :: nat
      real(dp) :: e_tot, worst

      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, grad, e_tot)
      if (allocated(error)) return
      worst = maxval(abs(grad - g_ref))
      print "(a,a,es10.2,a,es10.2,a,es10.2)", "  ", trim(c%name)//": |E - E_mqc| ", abs(e_tot - e_ref), &
         "  max |g - g_mqc| ", worst, "  |sum g| ", maxval(abs(sum(grad, dim=2)))
      call check(error, abs(e_tot - e_ref) < MQC_ENERGY_TOL, trim(c%name)// &
                 ": the RHF + RI-MP2 energy differs from metalquicha's -- not the same energy")
      if (allocated(error)) return
      call check(error, worst < MQC_GRAD_TOL, trim(c%name)//": gradient differs from metalquicha's")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(grad, dim=2))) < TRANS_TOL, trim(c%name)//": not translationally invariant")
   end subroutine against_mqc

   ! ---- B. finite differences of terco's own energy -------------------------

   subroutine test_fd_formamide_cart(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "formamide cart fc0"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%cart = .true.
      ! O x, N y, the carbonyl-side H x: three atoms, both in-plane axes.
      call against_fd(error, c, reshape([1, 2, 2, 3, 1, 4], [2, 3]))
   end subroutine test_fd_formamide_cart

   subroutine test_fd_formamide_sph_fc3(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "formamide sph fc3"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%nfrozen = 3
      ! C y, N x, an amide H y.
      call against_fd(error, c, reshape([2, 1, 1, 3, 2, 5], [2, 3]))
   end subroutine test_fd_formamide_sph_fc3

   subroutine test_fd_water_cart_fc1(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cart fc1"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%nfrozen = 1
      c%cart = .true.
      ! O z, H y, H z.
      call against_fd(error, c, reshape([3, 1, 2, 2, 3, 3], [2, 3]))
   end subroutine test_fd_water_cart_fc1

   subroutine against_fd(error, c, comps)
      !! `comps(:, k)` = (component, atom) of the k-th differenced coordinate.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: comps(:, :)
      real(dp), allocatable :: grad(:, :), r(:, :), rd(:, :)
      integer, allocatable :: z(:)
      integer :: nat, k, ic, ia
      real(dp) :: e_tot, e_plus, e_minus, fd, worst

      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, grad, e_tot)
      if (allocated(error)) return
      worst = 0.0_dp
      do k = 1, size(comps, 2)
         ic = comps(1, k)
         ia = comps(2, k)
         rd = r
         rd(ic, ia) = r(ic, ia) + FD_STEP
         call total_energy(error, c, nat, z, rd, e_plus)
         if (allocated(error)) return
         rd(ic, ia) = r(ic, ia) - FD_STEP
         call total_energy(error, c, nat, z, rd, e_minus)
         if (allocated(error)) return
         fd = (e_plus - e_minus)/(2.0_dp*FD_STEP)
         print "(a,a,i0,a,i0,a,f16.10,a,f16.10,a,es9.2)", "  ", trim(c%name)//": atom ", ia, " comp ", ic, &
            "  analytic ", grad(ic, ia), "  fd ", fd, "  diff ", abs(grad(ic, ia) - fd)
         worst = max(worst, abs(grad(ic, ia) - fd))
      end do
      call check(error, worst < FD_TOL, trim(c%name)//": analytic gradient differs from finite differences")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(grad, dim=2))) < TRANS_TOL, trim(c%name)//": not translationally invariant")
   end subroutine against_fd

   ! ---- shared ------------------------------------------------------------

   subroutine build_bases(error, c, nat, z, r, bas, aux)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      type(trc_basis_t), intent(out) :: bas, aux
      type(error_t) :: err
      ! Unpruned: the AO space metalquicha's and PySCF's readers build from
      ! the same file (pruning is span-preserving, so the energy would not
      ! move, but the comparison should not depend on that).
      call trc_basis_from_json(trim(c%basis), nat, z, r, bas, err, force_cart=c%cart, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (.not. err%has_error()) call trc_basis_from_json(trim(c%aux), nat, z, r, aux, err, force_cart=c%cart, &
                                                           overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), trim(c%name)//": basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call aux%to_device()
   end subroutine build_bases

   subroutine converged_scf(error, c, bas, nocc, res)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      type(trc_basis_t), intent(inout) :: bas
      integer, intent(in) :: nocc
      type(trc_scf_result_t), intent(out) :: res
      type(trc_scf_options_t) :: opts
      opts%conv_energy = 1.0e-13_dp
      opts%conv_density = 1.0e-11_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      call trc_scf_run(bas, nocc, nocc, opts, res)
      call check(error, res%converged, trim(c%name)//": the SCF did not converge")
      if (allocated(error)) return
      call check(error, res%n_lindep == 0, trim(c%name)//": functions dropped; the gradient needs them all")
   end subroutine converged_scf

   subroutine total_energy(error, c, nat, z, r, e_tot)
      !! E_RHF + E_RI-MP2 at geometry r, the way the energy path computes it.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      real(dp), intent(out) :: e_tot
      type(trc_basis_t) :: bas, aux
      type(trc_scf_result_t) :: res
      real(dp) :: e_corr
      integer :: nocc

      e_tot = 0.0_dp
      nocc = sum(z)/2
      call build_bases(error, c, nat, z, r, bas, aux)
      if (allocated(error)) return
      call converged_scf(error, c, bas, nocc, res)
      if (allocated(error)) return
      call rimp2_energy(error, c, bas, aux, nocc, res, e_corr)
      if (allocated(error)) return
      e_tot = res%energy + e_corr
      call aux%release()
      call bas%release()
   end subroutine total_energy

   subroutine rimp2_energy(error, c, bas, aux, nocc, res, e_corr)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      type(trc_basis_t), intent(in) :: bas, aux
      integer, intent(in) :: nocc
      type(trc_scf_result_t), intent(in) :: res
      real(dp), intent(out) :: e_corr
      type(trc_pairlist_t) :: pl
      type(trc_rimp2_result_t) :: mp
      type(trc_rimp2_options_t) :: ropts
      real(dp), allocatable :: cmo_c(:, :), eps_c(:)
      allocate (cmo_c(bas%nao, bas%nao), eps_c(bas%nao))
      if (bas%spherical) then
         call rimp2_driver_mo_to_cart(bas, res%cmo(:, :, 1), res%eps(:, 1), cmo_c, eps_c)
      else
         cmo_c = res%cmo(:, :, 1)
         eps_c = res%eps(:, 1)
      end if
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      ropts%nfrozen = c%nfrozen
      call trc_rimp2_run(bas, aux, pl, nocc, cmo_c, eps_c, mp, opts=ropts)
      call pl%release()
      e_corr = mp%e_corr
      call check(error, len_trim(mp%message) == 0, trim(c%name)//": RI-MP2: "//trim(mp%message))
   end subroutine rimp2_energy

   subroutine analytic_gradient(error, c, nat, z, r, grad, e_tot)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      real(dp), allocatable, intent(out) :: grad(:, :)
      real(dp), intent(out) :: e_tot
      type(trc_basis_t) :: bas, aux
      type(trc_scf_result_t) :: res
      type(trc_rimp2_grad_options_t) :: gopts
      type(trc_rimp2_grad_result_t) :: gres
      real(dp) :: e_corr
      integer :: nocc

      e_tot = 0.0_dp
      nocc = sum(z)/2
      allocate (grad(3, nat))
      call build_bases(error, c, nat, z, r, bas, aux)
      if (allocated(error)) return
      call converged_scf(error, c, bas, nocc, res)
      if (allocated(error)) return
      call rimp2_energy(error, c, bas, aux, nocc, res, e_corr)
      if (allocated(error)) return
      e_tot = res%energy + e_corr
      gopts%nocc = nocc
      gopts%nfrozen = c%nfrozen
      call trc_rimp2_gradient(bas, aux, res%cmo(:, :, 1), res%eps(:, 1), gopts, grad, gres)
      call check(error, gres%ok, trim(c%name)//": "//trim(gres%message))
      if (allocated(error)) return
      print "(a,a,es10.2,a,i0)", "  ", trim(c%name)//": E_corr (gradient) - E_corr (energy) ", &
         gres%e_corr - e_corr, "  Z-vector iterations ", gres%cphf_iterations
      call check(error, abs(gres%e_corr - e_corr) < ECORR_TOL, trim(c%name)// &
                 ": the gradient's amplitudes give a different correlation energy than trc_rimp2_run")
      call aux%release()
      call bas%release()
   end subroutine analytic_gradient

end module check_rimp2_grad_suite

program check_rimp2_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rimp2_grad_suite, only: collect_rimp2_grad
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_rimp2_grad, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rimp2_grad
