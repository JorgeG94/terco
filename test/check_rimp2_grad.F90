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
! C. The --df reference (`dfjk`/`jkaux`, the SCF's own J and K fitted):
!    metalquicha's "RI-MP2 gradient over a fitted reference NH3
!    def2-svp/def2-svp-rifit" (same validation file; its reference is a
!    Richardson-extrapolated finite difference of PySCF 2.14's own
!    density-fitted RHF + DF-MP2 energy, both fits in def2-SVP-RIFIT,
!    gen_cpu_validation.py's pyscf_df_ref_mp2_gradient), and central
!    finite differences of terco's own --df RHF + RI-MP2 energy, spherical
!    with def2-universal-jkfit for the SCF and frozen core, and Cartesian
!    with one set doing both fits. Plus an internal identity:
!    scf_grad_df_bilinear at X = D/2 is the --df SCF's own two-electron
!    gradient (scf_grad_df_two_electron).
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
   use trc_rimp2_grad, only: trc_rimp2_gradient, trc_rimp2_grad_t
   use trc_scf_df, only: df_jk_t
   use trc_scf_grad_df, only: scf_grad_df_two_electron, scf_grad_df_bilinear, scf_grad_df_bilinear_t
   use trc_linalg, only: trc_linalg_t
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
   real(dp), parameter :: MQC_DFREF_TOL = 2.0e-8_dp
   !! Against a finite-difference reference (metalquicha's own gate is
   !! 1e-7; its measured agreement is a few times 1e-9).
   real(dp), parameter :: BILINEAR_TOL = 1.0e-11_dp

   type :: case_t
      character(len=64) :: name = ""
      character(len=64) :: xyz = ""
      character(len=64) :: basis = ""
      character(len=64) :: aux = ""
      character(len=64) :: jkfit = ""
      !! The --df reference's fitting set; empty for the exact reference.
      integer :: nfrozen = 0
      logical :: cart = .false.
   end type case_t

   type :: system_t
      !! One converged reference and what the gradient needs of it.
      type(trc_basis_t) :: bas, aux, jk
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: lp
      type(df_jk_t) :: dfjk
      type(trc_scf_result_t) :: res
      integer :: nocc = 0
      logical :: df = .false.
   end type system_t

contains

   subroutine collect_rimp2_grad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water cc-pVDZ/cc-pVDZ-RIFIT (sph) vs metalquicha", test_mqc_water_ccpvdz), &
                  new_unittest("water 6-31G/cc-pVDZ-RIFIT, frozen 1, vs metalquicha", test_mqc_water_631g_fc1), &
                  new_unittest("NH3 def2-SVP/def2-SVP-RIFIT (sph) vs metalquicha", test_mqc_nh3_def2svp), &
                  new_unittest("formamide cc-pVDZ (Cartesian) vs finite differences", test_fd_formamide_cart), &
                  new_unittest("formamide cc-pVDZ (sph), frozen 3, vs finite differences", test_fd_formamide_sph_fc3), &
                  new_unittest("water cc-pVDZ (Cartesian), frozen 1, vs finite differences", test_fd_water_cart_fc1), &
                  new_unittest("--df: NH3 def2-SVP/def2-SVP-RIFIT (sph) vs metalquicha", test_mqc_dfref_nh3), &
                  new_unittest("--df: formamide cc-pVDZ (sph), universal-jkfit, frozen 3, vs finite differences", &
                               test_fd_df_formamide_sph), &
                  new_unittest("--df: water cc-pVDZ (Cartesian) vs finite differences", test_fd_df_water_cart), &
                  new_unittest("--df: the bilinear derivative at X = D/2 is the SCF's own", test_bilinear_identity), &
                  new_unittest("--df: the AO fallback of the fitted orbital Hessian agrees with the MO one", &
                               test_df_response_ao) &
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

   ! ---- C. the --df reference ----------------------------------------------

   subroutine test_mqc_dfref_nh3(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      real(dp), allocatable :: grad(:, :), r(:, :)
      integer, allocatable :: z(:)
      integer :: nat
      real(dp) :: e_tot, worst
      real(dp), parameter :: E_MQC = -56.3395970169_dp
      real(dp), parameter :: G_MQC(3, 4) = reshape([ &
                                           -4e-11_dp, -2.9e-11_dp, 0.008288657367_dp, &
                                           -0.002408829768_dp, -2.38e-10_dp, -0.002762885685_dp, &
                                           0.001204414513_dp, -0.00208610766_dp, -0.002762886289_dp, &
                                           0.001204414784_dp, 0.00208610774_dp, -0.002762885861_dp], [3, 4])
      c%name = "NH3 def2-SVP --df fc0"
      c%xyz = "nh3_mqc.xyz"
      c%basis = "basis_sets/def2-svp.json"
      c%aux = "basis_sets/def2-svp-rifit.json"
      c%jkfit = "basis_sets/def2-svp-rifit.json"
      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, grad, e_tot)
      if (allocated(error)) return
      worst = maxval(abs(grad - G_MQC))
      print "(a,a,es10.2,a,es10.2,a,es10.2)", "  ", trim(c%name)//": |E - E_mqc| ", abs(e_tot - E_MQC), &
         "  max |g - g_mqc| ", worst, "  |sum g| ", maxval(abs(sum(grad, dim=2)))
      call check(error, abs(e_tot - E_MQC) < MQC_ENERGY_TOL, trim(c%name)// &
                 ": the --df RHF + RI-MP2 energy differs from metalquicha's -- not the same energy")
      if (allocated(error)) return
      call check(error, worst < MQC_DFREF_TOL, trim(c%name)//": gradient differs from metalquicha's")
      if (allocated(error)) return
      call check(error, maxval(abs(sum(grad, dim=2))) < TRANS_TOL, trim(c%name)//": not translationally invariant")
   end subroutine test_mqc_dfref_nh3

   subroutine test_fd_df_formamide_sph(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "formamide sph --df fc3"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/def2-universal-jkfit.json"
      c%nfrozen = 3
      call against_fd(error, c, reshape([1, 2, 2, 3, 2, 6], [2, 3]))
   end subroutine test_fd_df_formamide_sph

   subroutine test_fd_df_water_cart(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      c%name = "water cart --df fc0"
      c%xyz = "water.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/cc-pvdz-rifit.json"
      c%cart = .true.
      call against_fd(error, c, reshape([3, 1, 2, 2, 3, 3], [2, 3]))
   end subroutine test_fd_df_water_cart

   subroutine test_bilinear_identity(error)
      !! d[Tr (D/2) G_fit(D)] through scf_grad_df_bilinear against the
      !! --df SCF gradient's own two-electron term, at a converged --df
      !! SCF (Cartesian, where scf_grad_df_two_electron runs).
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      type(system_t) :: sys
      real(dp), allocatable :: r(:, :), cmo(:, :), dens(:, :), g_ref(:, :), g_bil(:, :)
      type(scf_grad_df_bilinear_t) :: bil
      integer, allocatable :: z(:)
      integer :: nat, n
      real(dp) :: e_tot
      c%name = "formamide cart --df identity"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/cc-pvdz-rifit.json"
      c%cart = .true.
      call read_xyz(trim(c%xyz), nat, z, r)
      call setup_system(error, c, nat, z, r, sys, e_tot)
      if (allocated(error)) return
      n = sys%bas%nao
      cmo = sys%res%cmo(:, :, 1)
      allocate (dens(n, n), g_ref(3, nat), g_bil(3, nat))
      dens = 2.0_dp*matmul(cmo(:, 1:sys%nocc), transpose(cmo(:, 1:sys%nocc)))
      g_ref = 0.0_dp
      g_bil = 0.0_dp
      !$acc enter data copyin(cmo, dens, g_ref)
      call sys%jk%to_device()
      call scf_grad_df_two_electron(sys%bas, sys%jk, sys%pl, sys%dfjk, sys%lp, sys%nocc, cmo, dens, &
                                    1.0_dp, 1.0_dp, g_ref)
      !$acc update self(g_ref)
      !$acc exit data delete(cmo, dens, g_ref)
      call scf_grad_df_bilinear(sys%bas, sys%jk, sys%pl, sys%dfjk, sys%lp, cmo(:, 1:sys%nocc), 0.5_dp*dens, bil)
      g_bil = bil%grad
      print "(a,a,es10.2)", "  ", trim(c%name)//": max |bilinear(D/2) - two-electron| ", maxval(abs(g_bil - g_ref))
      call check(error, maxval(abs(g_bil - g_ref)) < BILINEAR_TOL, trim(c%name)// &
                 ": scf_grad_df_bilinear(X = D/2) is not the --df SCF's two-electron gradient")
      call release_system(sys)
   end subroutine test_bilinear_identity

   subroutine test_df_response_ao(error)
      !! The Z-vector's memory fallback (df_rotation_t's AO path) against
      !! its default MO path, formamide spherical --df.
      type(error_type), allocatable, intent(out) :: error
      type(case_t) :: c
      real(dp), allocatable :: r(:, :), g_mo(:, :), g_ao(:, :)
      integer, allocatable :: z(:)
      integer :: nat
      real(dp) :: e_tot
      c%name = "formamide sph --df AO fallback"
      c%xyz = "formamide.xyz"
      c%basis = "basis_sets/cc-pvdz.json"
      c%aux = "basis_sets/cc-pvdz-rifit.json"
      c%jkfit = "basis_sets/def2-universal-jkfit.json"
      c%nfrozen = 3
      call read_xyz(trim(c%xyz), nat, z, r)
      call analytic_gradient(error, c, nat, z, r, g_mo, e_tot)
      if (allocated(error)) return
      call analytic_gradient(error, c, nat, z, r, g_ao, e_tot, force_ao=.true.)
      if (allocated(error)) return
      print "(a,a,es10.2)", "  ", trim(c%name)//": max |g(AO) - g(MO)| ", maxval(abs(g_ao - g_mo))
      call check(error, maxval(abs(g_ao - g_mo)) < BILINEAR_TOL, trim(c%name)// &
                 ": the two forms of the fitted orbital Hessian give different gradients")
   end subroutine test_df_response_ao

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

   subroutine setup_system(error, c, nat, z, r, sys, e_tot)
      !! The bases, the converged reference (four-centre, or --df over
      !! `c%jkfit`) and E_RHF + E_RI-MP2 at geometry r, the way the energy
      !! path computes them.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      type(system_t), intent(inout) :: sys
      real(dp), intent(out) :: e_tot
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      real(dp) :: e_corr

      e_tot = 0.0_dp
      sys%nocc = sum(z)/2
      sys%df = len_trim(c%jkfit) > 0
      ! Unpruned: the AO space metalquicha's and PySCF's readers build from
      ! the same file (pruning is span-preserving, so the energy would not
      ! move, but the comparison should not depend on that).
      call trc_basis_from_json(trim(c%basis), nat, z, r, sys%bas, err, force_cart=c%cart, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      if (.not. err%has_error()) call trc_basis_from_json(trim(c%aux), nat, z, r, sys%aux, err, force_cart=c%cart, &
                                                           overrides=trc_basis_json_overrides(prune=.false.))
      if (sys%df .and. .not. err%has_error()) then
         call trc_basis_from_json(trim(c%jkfit), nat, z, r, sys%jk, err, force_cart=c%cart, &
                                  overrides=trc_basis_json_overrides(prune=.false.))
      end if
      call check(error, .not. err%has_error(), trim(c%name)//": basis: "//err%get_message())
      if (allocated(error)) return
      call sys%bas%to_device()
      call sys%aux%to_device()
      call sys%pl%build(sys%bas, 1.0e-12_dp)
      call sys%pl%to_device()
      call sys%lp%init(sys%bas%nao)
      opts%conv_energy = 1.0e-13_dp
      opts%conv_density = 1.0e-11_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      if (sys%df) then
         call sys%jk%to_device()
         call sys%dfjk%setup(sys%bas, sys%jk, sys%pl, sys%lp, err, nocc_hint=sys%nocc)
         call check(error, .not. err%has_error(), trim(c%name)//": fit: "//err%get_message())
         if (allocated(error)) return
         opts%df = .true.
         call trc_scf_run(sys%bas, sys%nocc, sys%nocc, opts, sys%res, la=sys%lp, scf_dfjk=sys%dfjk)
      else
         call trc_scf_run(sys%bas, sys%nocc, sys%nocc, opts, sys%res, la=sys%lp)
      end if
      call check(error, sys%res%converged, trim(c%name)//": the SCF did not converge")
      if (allocated(error)) return
      call check(error, sys%res%n_lindep == 0, trim(c%name)//": functions dropped; the gradient needs them all")
      if (allocated(error)) return
      call rimp2_energy(error, c, sys, e_corr)
      if (allocated(error)) return
      e_tot = sys%res%energy + e_corr
   end subroutine setup_system

   subroutine release_system(sys)
      type(system_t), intent(inout) :: sys
      if (sys%df) then
         call sys%dfjk%release()
         call sys%jk%release()
      end if
      call sys%lp%release()
      call sys%pl%release()
      call sys%aux%release()
      call sys%bas%release()
   end subroutine release_system

   subroutine total_energy(error, c, nat, z, r, e_tot)
      !! E_RHF + E_RI-MP2 at geometry r.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      real(dp), intent(out) :: e_tot
      type(system_t) :: sys
      call setup_system(error, c, nat, z, r, sys, e_tot)
      if (allocated(error)) return
      call release_system(sys)
   end subroutine total_energy

   subroutine rimp2_energy(error, c, sys, e_corr)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      type(system_t), intent(inout) :: sys
      real(dp), intent(out) :: e_corr
      type(trc_rimp2_result_t) :: mp
      type(trc_rimp2_options_t) :: ropts
      real(dp), allocatable :: cmo_c(:, :), eps_c(:)
      allocate (cmo_c(sys%bas%nao, sys%bas%nao), eps_c(sys%bas%nao))
      if (sys%bas%spherical) then
         call rimp2_driver_mo_to_cart(sys%bas, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), cmo_c, eps_c)
      else
         cmo_c = sys%res%cmo(:, :, 1)
         eps_c = sys%res%eps(:, 1)
      end if
      ropts%nfrozen = c%nfrozen
      call trc_rimp2_run(sys%bas, sys%aux, sys%pl, sys%nocc, cmo_c, eps_c, mp, opts=ropts, la=sys%lp)
      e_corr = mp%e_corr
      call check(error, len_trim(mp%message) == 0, trim(c%name)//": RI-MP2: "//trim(mp%message))
   end subroutine rimp2_energy

   subroutine analytic_gradient(error, c, nat, z, r, grad, e_tot, force_ao)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), intent(in) :: c
      integer, intent(in) :: nat, z(nat)
      real(dp), intent(in) :: r(3, nat)
      real(dp), allocatable, intent(out) :: grad(:, :)
      real(dp), intent(out) :: e_tot
      logical, intent(in), optional :: force_ao
      type(system_t) :: sys
      type(trc_rimp2_grad_t) :: job
      real(dp) :: e_corr

      call setup_system(error, c, nat, z, r, sys, e_tot)
      if (allocated(error)) return
      e_corr = e_tot - sys%res%energy
      job%opts%nocc = sys%nocc
      job%opts%nfrozen = c%nfrozen
      if (present(force_ao)) job%opts%df_response_ao = force_ao
      if (sys%df) then
         call trc_rimp2_gradient(sys%bas, sys%aux, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), job, sys%lp, &
                                 dfjk=sys%dfjk, jkaux=sys%jk)
      else
         call trc_rimp2_gradient(sys%bas, sys%aux, sys%res%cmo(:, :, 1), sys%res%eps(:, 1), job, sys%lp)
      end if
      call check(error, job%res%ok, trim(c%name)//": "//trim(job%res%message))
      if (allocated(error)) return
      grad = job%grad
      print "(a,a,es10.2,a,i0)", "  ", trim(c%name)//": E_corr (gradient) - E_corr (energy) ", &
         job%res%e_corr - e_corr, "  Z-vector iterations ", job%res%cphf_iterations
      call check(error, abs(job%res%e_corr - e_corr) < ECORR_TOL, trim(c%name)// &
                 ": the gradient's amplitudes give a different correlation energy than trc_rimp2_run")
      call release_system(sys)
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
