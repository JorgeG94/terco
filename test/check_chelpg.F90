module check_chelpg_suite
   !!
   !! CHELPG atomic charges (trc_chelpg.F90) against two references:
   !!
   !!  * metalquicha's recipe run in pyscf (test/pyscf_refs/chelpg.ref, written
   !!    offline by `generate.py chelpg`): the same Breneman-Wiberg grid rule
   !!    and constrained fit, with pyscf's own SCF and its `int1e_grids`
   !!    potential. Water RHF and B3LYP/6-31G, formamide and the acetate anion
   !!    RHF and B3LYP/def2-SVP, and the NH2 radical UHF/cc-pVDZ. Each case also
   !!    checks that terco's grid has exactly the reference's point count.
   !!  * metalquicha's OWN Fortran charges (mqc_czt_charges' chelpg_charges at
   !!    mqc commit 66efbeed96, the geometry of its validation/check_charges.f90,
   !!    the same 6-31G and cc-pVDZ files byte for byte), hardcoded below.
   !!
   !! And the properties any CHELPG fit owes: the charges sum to the molecular
   !! charge (it is the constraint), they do not move when the molecule is
   !! translated and they permute with the atoms, the grid lies outside every
   !! van der Waals sphere and inside the head space, the fit reproduces the
   !! potential (RRMS), and a ghost atom or an element without a radius is
   !! refused.
   !!
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_sph, only: sph_sph_to_cart
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_chelpg, only: fortran_chelpg_charges => trc_chelpg_charges, trc_chelpg_grid, trc_chelpg_vdw_radius, &
                         trc_chelpg_fit_t, CHELPG_HEAD_SPACE
   use trc_constants, only: ANGSTROM_TO_BOHR
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_elements, &
                               trc_set_basis_json, trc_set_method, trc_set_convergence, trc_run_scf, &
                               capi_chelpg_charges => trc_chelpg_charges, TRC_OK, TRC_ERR_BADARG, &
                               TRC_ERR_STATE
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_str, ref_read_real, &
                           ref_read_vector
   implicit none
   private
   public :: collect_chelpg_tests

   real(dp), parameter :: CHARGE_TOL = 3.0e-8_dp
      !! terco vs pyscf: two SCF codes converged to 1e-12 Eh on the same
      !! grid rule, with two independent potentials (the PCM kernel at the
      !! point-charge limit vs libcint's int1e_grids) and two fits.
      !! Measured worst atom (gfortran 15, host): 3.4e-9 over all seven cases,
      !! B3LYP included (the two codes' level-3 XC grids are the same
      !! construction: their energies agree to 9e-11 too); about ten times that.
   real(dp), parameter :: ENERGY_TOL = 1.0e-9_dp
      !! The SCF energy against pyscf's, so a charge mismatch is never a
      !! different state. Measured worst 9.4e-11 (acetate B3LYP).
   real(dp), parameter :: MQC_CHARGE_TOL = 1.0e-7_dp
      !! terco vs mqc's own Fortran: mqc's SCF (1e-12 Eh, 1e-10 DIIS) against
      !! terco's, same grid, same basis files. Measured worst 1.1e-8 (water/6-31G).
   real(dp), parameter :: SUM_TOL = 1.0e-12_dp
      !! The total-charge constraint is solved exactly; this is round-off in
      !! a (natm+1)-square LU solve. Measured worst 3.3e-16.
   real(dp), parameter :: INVARIANCE_TOL = 3.0e-12_dp
      !! Same density, moved or relabelled molecule: only round-off in the
      !! sums over grid points and the LU pivot order may differ. Measured
      !! 3.0e-13 (translation) and 2.8e-13 (reorder).

   ! mqc's OWN numbers (validation/check_charges.f90's geometry, Bohr; RHF,
   ! mqc's chelpg_charges at its defaults), printed to 16 digits from a
   ! build of mqc commit 66efbeed96: water/6-31G, water/cc-pVDZ, methane/6-31G.
   real(dp), parameter :: MQC_WATER_R(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   real(dp), parameter :: MQC_METHANE_R(3, 5) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           1.19_dp, 1.19_dp, 1.19_dp, &
                           -1.19_dp, -1.19_dp, 1.19_dp, &
                           -1.19_dp, 1.19_dp, -1.19_dp, &
                           1.19_dp, -1.19_dp, -1.19_dp], [3, 5])
   real(dp), parameter :: MQC_Q_WATER_631G(3) = [-9.4383200635358688e-01_dp, 4.7200801378211721e-01_dp, &
                                                 4.7182399257146967e-01_dp]
   real(dp), parameter :: MQC_Q_WATER_CCPVDZ(3) = [-7.5073713272426579e-01_dp, 3.7519725305819784e-01_dp, &
                                                   3.7553987966606805e-01_dp]
   real(dp), parameter :: MQC_Q_METHANE_631G(5) = [-4.4679143871066740e-01_dp, 1.1271072823146555e-01_dp, &
                                                   1.1136023682640390e-01_dp, 1.1136023682639418e-01_dp, &
                                                   1.1136023682640379e-01_dp]
   integer, parameter :: MQC_NPTS_WATER = 4365
   integer, parameter :: MQC_NPTS_METHANE = 5868

   type :: scf_case_t
      character(len=64) :: xyz = "", basis = "", functional = ""
      integer :: charge = 0, n_unpaired = 0
   end type scf_case_t

contains

   subroutine collect_chelpg_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water RHF/6-31G matches mqc's recipe in pyscf", test_ref_1), &
                  new_unittest("water B3LYP/6-31G matches mqc's recipe in pyscf", test_ref_2), &
                  new_unittest("formamide RHF/def2-SVP matches mqc's recipe in pyscf", test_ref_3), &
                  new_unittest("formamide B3LYP/def2-SVP matches mqc's recipe in pyscf", test_ref_4), &
                  new_unittest("acetate RHF/def2-SVP matches mqc's recipe in pyscf", test_ref_5), &
                  new_unittest("acetate B3LYP/def2-SVP matches mqc's recipe in pyscf", test_ref_6), &
                  new_unittest("NH2 UHF/cc-pVDZ matches mqc's recipe in pyscf", test_ref_7), &
                  new_unittest("water and methane match mqc's OWN charges", test_mqc), &
                  new_unittest("charges are invariant to translation", test_translation), &
                  new_unittest("charges permute with the atoms", test_reorder), &
                  new_unittest("grid is outside the vdW spheres, inside the head space", test_grid), &
                  new_unittest("ghost atom and radius-less element are refused", test_refused), &
                  new_unittest("the C ABI matches pyscf for the anion and the radical", test_capi), &
                  new_unittest("the C ABI refuses a ghost and a run with no SCF", test_capi_refused) &
                  ]
   end subroutine collect_chelpg_tests

   subroutine test_ref_1(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(1, scf_case_t("water.xyz", "basis_sets/6-31g.json", "", 0, 0), error)
   end subroutine test_ref_1

   subroutine test_ref_2(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(2, scf_case_t("water.xyz", "basis_sets/6-31g.json", "b3lyp", 0, 0), error)
   end subroutine test_ref_2

   subroutine test_ref_3(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(3, scf_case_t("formamide.xyz", "basis_sets/def2-svp.json", "", 0, 0), error)
   end subroutine test_ref_3

   subroutine test_ref_4(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(4, scf_case_t("formamide.xyz", "basis_sets/def2-svp.json", "b3lyp", 0, 0), error)
   end subroutine test_ref_4

   subroutine test_ref_5(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(5, scf_case_t("acetate.xyz", "basis_sets/def2-svp.json", "", -1, 0), error)
   end subroutine test_ref_5

   subroutine test_ref_6(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(6, scf_case_t("acetate.xyz", "basis_sets/def2-svp.json", "b3lyp", -1, 0), error)
   end subroutine test_ref_6

   subroutine test_ref_7(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_case(7, scf_case_t("nh2.xyz", "basis_sets/cc-pvdz.json", "", 0, 1), error)
   end subroutine test_ref_7

   subroutine compare_case(ic, case, error)
      !! Case `ic` of chelpg.ref: terco's SCF, terco's CHELPG, against pyscf's
      integer, intent(in) :: ic
      type(scf_case_t), intent(in) :: case
      type(error_type), allocatable, intent(out) :: error

      type(trc_basis_t) :: b
      type(trc_chelpg_fit_t) :: fit
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), d_cart(:, :), q(:), q_ref(:)
      integer, allocatable :: z(:)
      character(len=64) :: name
      real(dp) :: e_ref, e_terco, worst
      integer :: nat, npts_ref, grid_level

      call read_ref_case(ic, name, npts_ref, e_ref, q_ref, grid_level)
      call read_xyz(trim(case%xyz), nat, z, r)
      call check(error, nat == size(q_ref), trim(name)//": atom count differs from the reference")
      if (allocated(error)) return
      call converged_density(z, r, case, grid_level, b, d_cart, e_terco, error)
      if (allocated(error)) return

      call fortran_chelpg_charges(b, d_cart, real(case%charge, dp), q, err, fit=fit)
      call check(error,.not. err%has_error(), trim(name)//": "//err%get_message())
      if (allocated(error)) return

      worst = maxval(abs(q - q_ref))
      call check(error, fit%n_points == npts_ref, trim(name)//": terco's grid has "// &
                 to_string(fit%n_points)//" points, the reference's "//to_string(npts_ref))
      if (allocated(error)) return
      call check(error, abs(e_terco - e_ref) < ENERGY_TOL, trim(name)//": |E_terco - E_pyscf| = "// &
                 to_text(abs(e_terco - e_ref))//", a different state from the reference's")
      if (allocated(error)) return
      call check(error, worst < CHARGE_TOL, trim(name)//": worst |q_terco - q_pyscf| = "//to_text(worst)// &
                 " exceeds "//to_text(CHARGE_TOL))
      if (allocated(error)) return
      call check(error, abs(sum(q) - real(case%charge, dp)) < SUM_TOL, trim(name)// &
                 ": the charges do not sum to the molecular charge")
      if (allocated(error)) return
      ! Breneman-Wiberg's own quality bar for a polar molecule is RRMS < ~0.1;
      ! the radical's is larger (its potential is weaker, 0.21 in pyscf too).
      call check(error, fit%rrms < 0.25_dp, trim(name)//": RRMS "//to_text(fit%rrms)// &
                 ", the charges do not reproduce the potential")
   end subroutine compare_case

   subroutine test_mqc(error)
      !! mqc's own Fortran result, RHF at its own geometry
      type(error_type), allocatable, intent(out) :: error
      call compare_mqc("water/6-31G", [8, 1, 1], MQC_WATER_R, "basis_sets/6-31g.json", &
                       MQC_Q_WATER_631G, MQC_NPTS_WATER, error)
      if (allocated(error)) return
      call compare_mqc("water/cc-pVDZ", [8, 1, 1], MQC_WATER_R, "basis_sets/cc-pvdz.json", &
                       MQC_Q_WATER_CCPVDZ, MQC_NPTS_WATER, error)
      if (allocated(error)) return
      call compare_mqc("methane/6-31G", [6, 1, 1, 1, 1], MQC_METHANE_R, "basis_sets/6-31g.json", &
                       MQC_Q_METHANE_631G, MQC_NPTS_METHANE, error)
   end subroutine test_mqc

   subroutine compare_mqc(label, z, r, basis, q_mqc, npts_mqc, error)
      character(len=*), intent(in) :: label, basis
      integer, intent(in) :: z(:), npts_mqc
      real(dp), intent(in) :: r(:, :), q_mqc(:)
      type(error_type), allocatable, intent(out) :: error

      type(trc_basis_t) :: b
      type(trc_chelpg_fit_t) :: fit
      type(error_t) :: err
      real(dp), allocatable :: d_cart(:, :), q(:)
      real(dp) :: e_terco

      call converged_density(z, r, scf_case_t("", basis, "", 0, 0), 3, b, d_cart, e_terco, error)
      if (allocated(error)) return
      call fortran_chelpg_charges(b, d_cart, 0.0_dp, q, err, fit=fit)
      call check(error,.not. err%has_error(), label//": "//err%get_message())
      if (allocated(error)) return
      call check(error, fit%n_points == npts_mqc, label//": "//to_string(fit%n_points)// &
                 " grid points, mqc's grid has "//to_string(npts_mqc))
      if (allocated(error)) return
      call check(error, maxval(abs(q - q_mqc)) < MQC_CHARGE_TOL, label//": worst |q_terco - q_mqc| = "// &
                 to_text(maxval(abs(q - q_mqc)))//" exceeds "//to_text(MQC_CHARGE_TOL))
   end subroutine compare_mqc

   subroutine test_translation(error)
      !! The same density on a rigidly moved molecule: the lattice is anchored
      !! to the bounding box, so it moves with the atoms and so must nothing else
      type(error_type), allocatable, intent(out) :: error

      type(trc_basis_t) :: b, b_moved
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), d_cart(:, :), q(:), q_moved(:)
      integer, allocatable :: z(:)
      real(dp) :: e_terco
      integer :: nat, ia
      real(dp), parameter :: SHIFT(3) = [0.731_dp, -1.927_dp, 3.113_dp]

      call read_xyz("formamide.xyz", nat, z, r)
      call converged_density(z, r, scf_case_t("", "basis_sets/6-31g.json", "", 0, 0), 3, b, d_cart, &
                             e_terco, error)
      if (allocated(error)) return
      call fortran_chelpg_charges(b, d_cart, 0.0_dp, q, err)
      call check(error,.not. err%has_error(), "formamide: "//err%get_message())
      if (allocated(error)) return

      do ia = 1, nat
         r(:, ia) = r(:, ia) + SHIFT
      end do
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r, b_moved, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call fortran_chelpg_charges(b_moved, d_cart, 0.0_dp, q_moved, err)
      call check(error,.not. err%has_error(), "moved formamide: "//err%get_message())
      if (allocated(error)) return
      call check(error, maxval(abs(q - q_moved)) < INVARIANCE_TOL, "translation moved a charge by "// &
                 to_text(maxval(abs(q - q_moved))))
   end subroutine test_translation

   subroutine test_reorder(error)
      !! Formamide with its atoms listed in another order, and the SAME
      !! density relabelled to match: each atom keeps its charge
      type(error_type), allocatable, intent(out) :: error

      type(trc_basis_t) :: b, b_perm
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), d_cart(:, :), d_perm(:, :), q(:), q_perm(:)
      integer, allocatable :: z(:), ao_map(:)
      real(dp) :: e_terco
      integer :: nat, ia, i, j
      integer, parameter :: ORDER(6) = [5, 3, 1, 6, 2, 4]   !! new atom k is old atom ORDER(k)

      call read_xyz("formamide.xyz", nat, z, r)
      call converged_density(z, r, scf_case_t("", "basis_sets/6-31g.json", "", 0, 0), 3, b, d_cart, &
                             e_terco, error)
      if (allocated(error)) return
      call fortran_chelpg_charges(b, d_cart, 0.0_dp, q, err)
      call check(error,.not. err%has_error(), "formamide: "//err%get_message())
      if (allocated(error)) return

      call trc_basis_from_json("basis_sets/6-31g.json", nat, z(ORDER), r(:, ORDER), b_perm, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call atom_ao_map(b, b_perm, ORDER, ao_map)
      allocate (d_perm(b%nao, b%nao))
      do j = 1, b%nao
         do i = 1, b%nao
            d_perm(i, j) = d_cart(ao_map(i), ao_map(j))
         end do
      end do
      call fortran_chelpg_charges(b_perm, d_perm, 0.0_dp, q_perm, err)
      call check(error,.not. err%has_error(), "reordered formamide: "//err%get_message())
      if (allocated(error)) return
      do ia = 1, nat
         call check(error, abs(q_perm(ia) - q(ORDER(ia))) < INVARIANCE_TOL, "atom "// &
                    to_string(ORDER(ia))//" changed its charge by "//to_text(abs(q_perm(ia) - q(ORDER(ia))))// &
                    " when listed as atom "//to_string(ia))
         if (allocated(error)) return
      end do
   end subroutine test_reorder

   subroutine atom_ao_map(b, b_perm, order, ao_map)
      ! ao_map(i) = the AO of `b` that AO i of `b_perm` is: atom k of b_perm
      ! is atom order(k) of b, with the same shells in the same order.
      type(trc_basis_t), intent(in) :: b, b_perm
      integer, intent(in) :: order(:)
      integer, allocatable, intent(out) :: ao_map(:)
      integer :: k, s_new, s_old, m, n_cart

      allocate (ao_map(b_perm%nao))
      do k = 1, size(order)
         s_old = 0
         do s_new = 1, b_perm%nshell
            if (b_perm%sh_at(s_new) /= k) cycle
            s_old = next_shell_of(b, order(k), s_old)
            n_cart = (b_perm%sh_l(s_new) + 1)*(b_perm%sh_l(s_new) + 2)/2
            do m = 0, n_cart - 1
               ao_map(b_perm%sh_ao(s_new) + m) = b%sh_ao(s_old) + m
            end do
         end do
      end do
   end subroutine atom_ao_map

   pure function next_shell_of(b, atom, after) result(s)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: atom, after
      integer :: s
      do s = after + 1, b%nshell
         if (b%sh_at(s) == atom) return
      end do
      s = 0
   end function next_shell_of

   subroutine test_grid(error)
      !! mqc's test_grid_excludes: every point outside every sphere, every
      !! point within the head space of some nucleus, and enough of them
      type(error_type), allocatable, intent(out) :: error

      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), pts(:, :)
      integer, allocatable :: z(:)
      real(dp) :: dist, nearest, margin
      integer :: nat, g, ia

      call read_xyz("acetate.xyz", nat, z, r)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r, b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call trc_chelpg_grid(b, pts, err)
      call check(error,.not. err%has_error(), "grid: "//err%get_message())
      if (allocated(error)) return
      call check(error, size(pts, 2) > 10*nat, "the grid has only "//to_string(size(pts, 2))//" points")
      if (allocated(error)) return
      margin = CHELPG_HEAD_SPACE*ANGSTROM_TO_BOHR
      do g = 1, size(pts, 2)
         nearest = huge(1.0_dp)
         do ia = 1, nat
            dist = norm2(pts(:, g) - r(:, ia))
            call check(error, dist >= trc_chelpg_vdw_radius(z(ia)), "point "//to_string(g)// &
                       " is inside atom "//to_string(ia)//"'s van der Waals sphere")
            if (allocated(error)) return
            nearest = min(nearest, dist)
         end do
         call check(error, nearest <= margin, "point "//to_string(g)//" is beyond the head space")
         if (allocated(error)) return
      end do
   end subroutine test_grid

   subroutine test_refused(error)
      !! A ghost (charge zero) and Zr (no radius in the table) are refused
      !! rather than given a charge
      type(error_type), allocatable, intent(out) :: error

      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: r(:, :), q(:), d(:, :)
      integer, allocatable :: z(:)
      integer :: nat

      call read_xyz("water.xyz", nat, z, r)
      call trc_basis_from_json("basis_sets/6-31g.json", nat, z, r, b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      b%at_z(3) = 0.0_dp
      allocate (d(b%nao, b%nao))
      d = 0.0_dp
      call fortran_chelpg_charges(b, d, 0.0_dp, q, err)
      call check(error, err%has_error(), "a ghost atom was not refused")
      if (allocated(error)) return
      call err%clear()

      call trc_basis_from_json("basis_sets/sto-3g-rb-zr.json", 1, [40], reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]), &
                               b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      deallocate (d)
      allocate (d(b%nao, b%nao))
      d = 0.0_dp
      call fortran_chelpg_charges(b, d, 0.0_dp, q, err)
      call check(error, err%has_error(), "Zr, which has no van der Waals radius, was not refused")
   end subroutine test_refused

   subroutine test_capi(error)
      !! trc_chelpg_charges through the C ABI, on the cases that exercise its
      !! own conversion: a spherical basis (def2-SVP, cc-pVDZ, so the density
      !! goes through sph_sph_to_cart), the anion's charge as the constraint,
      !! and the radical's two spin densities summed
      type(error_type), allocatable, intent(out) :: error
      call capi_case(5, "acetate.xyz", "basis_sets/def2-svp.json", -1, 1, error)
      if (allocated(error)) return
      call capi_case(7, "nh2.xyz", "basis_sets/cc-pvdz.json", 0, 2, error)
   end subroutine test_capi

   subroutine capi_case(ic, xyz, basis, charge, multiplicity, error)
      integer, intent(in) :: ic, charge, multiplicity
      character(len=*), intent(in) :: xyz, basis
      type(error_type), allocatable, intent(out) :: error

      type(c_ptr) :: h
      real(dp), allocatable :: r(:, :), q_ref(:)
      real(c_double), allocatable :: q(:)
      integer, allocatable :: z(:)
      character(len=64) :: name
      real(dp) :: e_ref
      integer :: nat, npts_ref, grid_level
      integer(c_int) :: rc

      call read_ref_case(ic, name, npts_ref, e_ref, q_ref, grid_level)
      call read_xyz(xyz, nat, z, r)
      call capi_scf(h, z, r, basis, charge, multiplicity, error)
      if (allocated(error)) then
         rc = trc_destroy(h)
         return
      end if
      allocate (q(nat))
      rc = capi_chelpg_charges(h, q)
      call check(error, rc == TRC_OK, trim(name)//": trc_chelpg_charges failed")
      if (.not. allocated(error)) &
         call check(error, maxval(abs(real(q, dp) - q_ref)) < CHARGE_TOL, trim(name)// &
                    ": C ABI worst |q_terco - q_pyscf| = "//to_text(maxval(abs(real(q, dp) - q_ref))))
      if (.not. allocated(error)) &
         call check(error, abs(sum(real(q, dp)) - real(charge, dp)) < SUM_TOL, trim(name)// &
                    ": C ABI charges do not sum to the molecular charge")
      rc = trc_destroy(h)
   end subroutine capi_case

   subroutine test_capi_refused(error)
      !! Before any SCF: TRC_ERR_STATE. With a ghost hydrogen: TRC_ERR_BADARG.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      real(dp), allocatable :: r(:, :)
      real(c_double) :: q(3)
      integer, allocatable :: z(:)
      integer :: nat
      integer(c_int) :: rc

      call read_xyz("water.xyz", nat, z, r)
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(nat, c_int), real(z, c_double), reshape(r, [3*nat]), 0_c_int, 1_c_int)
      rc = capi_chelpg_charges(h, q)
      call check(error, rc == TRC_ERR_STATE, "trc_chelpg_charges before any SCF returned "//to_string(int(rc)))
      rc = trc_destroy(h)
      if (allocated(error)) return

      ! Water with its second hydrogen a ghost (basis functions, no nucleus):
      ! the OH radical, nine electrons, a doublet.
      rc = trc_create(h)
      rc = trc_set_molecule(h, int(nat, c_int), real([8, 1, 0], c_double), reshape(r, [3*nat]), 0_c_int, 2_c_int)
      rc = trc_set_elements(h, int(nat, c_int), [8_c_int, 1_c_int, 1_c_int])
      call check(error, rc == TRC_OK, "trc_set_elements failed")
      if (.not. allocated(error)) then
         rc = trc_set_basis_json(h, "basis_sets/6-31g.json"//c_null_char)
         rc = trc_run_scf(h)
         call check(error, rc == TRC_OK, "the ghosted SCF did not converge")
      end if
      if (.not. allocated(error)) then
         rc = capi_chelpg_charges(h, q)
         call check(error, rc == TRC_ERR_BADARG, "a ghost atom was not refused: "//to_string(int(rc)))
      end if
      rc = trc_destroy(h)
   end subroutine test_capi_refused

   subroutine capi_scf(h, z, r, basis, charge, multiplicity, error)
      type(c_ptr), intent(out) :: h
      integer, intent(in) :: z(:), charge, multiplicity
      real(dp), intent(in) :: r(:, :)
      character(len=*), intent(in) :: basis
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc

      rc = trc_create(h)
      rc = trc_set_molecule(h, int(size(z), c_int), real(z, c_double), reshape(r, [size(r)]), &
                            int(charge, c_int), int(multiplicity, c_int))
      call check(error, rc == TRC_OK, "trc_set_molecule failed")
      if (allocated(error)) return
      rc = trc_set_basis_json(h, basis//c_null_char)
      call check(error, rc == TRC_OK, "trc_set_basis_json failed")
      if (allocated(error)) return
      rc = trc_set_method(h, ""//c_null_char, 3_c_int)
      call check(error, rc == TRC_OK, "trc_set_method failed")
      if (allocated(error)) return
      rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 10_c_int)
      call check(error, rc == TRC_OK, "trc_set_convergence failed")
      if (allocated(error)) return
      rc = trc_run_scf(h)
      call check(error, rc == TRC_OK, "trc_run_scf did not converge")
   end subroutine capi_scf

   subroutine converged_density(z, r, case, grid_level, b, d_cart, energy, error)
      ! The SCF (RHF/UHF, or RKS/UKS with case%functional), and its total
      ! density in terco's Cartesian AO space -- the sph_to_cart step
      ! trc_capi's trc_chelpg_charges makes too.
      integer, intent(in) :: z(:), grid_level
      real(dp), intent(in) :: r(:, :)
      type(scf_case_t), intent(in) :: case
      type(trc_basis_t), intent(out) :: b
      real(dp), allocatable, intent(out) :: d_cart(:, :)
      real(dp), intent(out) :: energy
      type(error_type), allocatable, intent(out) :: error

      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      real(dp), allocatable :: d_scf(:, :)
      integer :: nelec, nalpha, nbeta

      call trc_basis_from_json(trim(case%basis), size(z), z, r, b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      nelec = sum(z) - case%charge
      nbeta = (nelec - case%n_unpaired)/2
      nalpha = nelec - nbeta
      opts%unrestricted = case%n_unpaired /= 0
      opts%functional = trim(case%functional)
      opts%grid_level = grid_level
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 200
      call trc_scf_run(b, nalpha, nbeta, opts, res)
      call check(error, res%converged, "the reference SCF did not converge")
      if (allocated(error)) return
      energy = res%energy

      allocate (d_scf(size(res%dmat, 1), size(res%dmat, 2)))
      if (res%nspin == 1) then
         d_scf = res%dmat(:, :, 1)
      else
         d_scf = res%dmat(:, :, 1) + res%dmat(:, :, 2)
      end if
      if (b%spherical) then
         allocate (d_cart(b%nao, b%nao))
         call sph_sph_to_cart(b, d_scf, d_cart)
      else
         d_cart = d_scf
      end if
   end subroutine converged_density

   subroutine read_ref_case(ic, name, npts, energy, q, grid_level)
      ! Case `ic` of chelpg.ref. Read in order: every case before `ic` is
      ! read and dropped.
      integer, intent(in) :: ic
      character(len=*), intent(out) :: name
      integer, intent(out) :: npts, grid_level
      real(dp), intent(out) :: energy
      real(dp), allocatable, intent(out) :: q(:)
      real(dp), allocatable :: v(:)
      integer :: u, ncase, k, nat

      call ref_require("chelpg.ref", u)
      call ref_read_int(u, "GRID_LEVEL", grid_level)
      call ref_read_int(u, "NCASE", ncase)
      do k = 1, min(ic, ncase)
         call ref_read_str(u, "NAME_"//to_string(k), name)
         call ref_read_int(u, "NATM_"//to_string(k), nat)
         call ref_read_int(u, "NPTS_"//to_string(k), npts)
         call ref_read_real(u, "E_"//to_string(k), energy)
         allocate (v(nat))
         call ref_read_vector(u, "Q_"//to_string(k), nat, v)
         if (k == ic) q = v
         deallocate (v)
      end do
      call ref_close(u)
   end subroutine read_ref_case

   function to_text(x) result(s)
      real(dp), intent(in) :: x
      character(len=32) :: s
      write (s, "(es12.4)") x
   end function to_text

end module check_chelpg_suite

program check_chelpg
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_chelpg_suite, only: collect_chelpg_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat
   stat = 0
   testsuites = [new_testsuite("chelpg", collect_chelpg_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_chelpg
