!
! Hirshfeld atomic charges (trc_hirshfeld.F90) against two references:
!
!  * metalquicha's OWN stored test numbers (its test/test_mqc_hirshfeld.f90:
!    water/cc-pVDZ, the neutral and the cation with spins), hardcoded below;
!  * metalquicha's recipe run in pyscf, for water RHF/6-31G and O2 triplet
!    UHF/cc-pVDZ: test/pyscf_refs/hirshfeld.ref, written offline by
!    `generate.py hirshfeld` (free atoms from mqc's
!    tools/spherical_atom/spherical_atom_reference.py, its commit recorded in
!    the .ref) and covered by the staleness hash like every other anchor.
!
! The pyscf integral runs on a (150, 974) atom grid and evaluates every free
! atom from its own orbitals; terco's grid (HIRSHFELD_GRID_LEVEL=5) and its
! cubic-Hermite proatom table are coarser, so agreement is bounded by the
! quadrature, not by the partition formula (CHARGE_TOL below).
!
! Also: the free-atom table integrates to Z, an element past Kr is refused,
! and the table loader recovers from a bad path and never silently keeps a
! different file than the one asked for.
!
module check_hirshfeld_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_error, only: error_t, ERROR_IO
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_sph, only: sph_sph_to_cart
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_hirshfeld, only: trc_hirshfeld_populations, trc_hirshfeld_free_atom_density
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_str, ref_read_vector
   implicit none
   private
   public :: collect_hirshfeld_tests


   ! mqc's OWN stored numbers (test/test_mqc_hirshfeld.f90, WATER/WATER_Z/
   ! Q_NEUTRAL/Q_CATION/S_CATION), not a pyscf cross-check: this is mqc's
   ! Fortran result itself (its free atoms solved at run time through its
   ! CPU HF/DFT backend, mqc_czt_spherical_atom), at mqc's own geometry --
   ! Bohr, hardcoded in that file, not read from an xyz (hence hardcoded
   ! here too rather than round-tripped through an Angstrom file).
   integer, parameter :: MQC_WATER_Z(3) = [8, 1, 1]
   real(dp), parameter :: MQC_WATER_R(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp], [3, 3])
   real(dp), parameter :: Q_NEUTRAL(3) = [-0.333925335674_dp, 0.166962667733_dp, &
                                          0.166962667733_dp]
   real(dp), parameter :: Q_CATION(3) = [0.335263944926_dp, 0.332368027475_dp, &
                                         0.332368027475_dp]
   real(dp), parameter :: S_CATION(3) = [0.956254823341_dp, 0.021872588392_dp, &
                                         0.021872588392_dp]
   real(dp), parameter :: MQC_CHARGE_TOL = 1.0e-7_dp
      !! mqc's own tolerance (test_mqc_hirshfeld.f90) for this EXACT
      !! comparison: mqc's free atoms solved through its CPU backend
      !! (libxc, mqc_czt_spherical_atom) vs terco's own proatom table
      !! (generated offline from the SAME recipe run in pyscf) and terco's
      !! own grid quadrature (HIRSHFELD_GRID_LEVEL=5) vs mqc's
      !! (HIRSHFELD_GRID_LEVEL=5 there too, but mqc's own radial/angular
      !! tables, not terco's). Three independent sources of small
      !! disagreement rather than the two in test_water/test_o2 above, so
      !! this uses mqc's own tolerance directly rather than tightening it.
      !! Measured (forcing MQC_CHARGE_TOL to 1e-12 to print the gap, as
      !! CHARGE_TOL's own comment does): water's worst atom differs from
      !! mqc's Q_NEUTRAL by 2.23e-9; the cation's worst atom differs from
      !! Q_CATION by 1.71e-9 and from S_CATION by 1.49e-9 -- all comfortably
      !! inside mqc's own 1e-7, and the same order of magnitude as
      !! test_water/test_o2's pyscf-side agreement, so the three extra
      !! sources of disagreement this comparison carries (mqc's own CPU
      !! SCF/grid rather than pyscf's) are not actually costing much here.

   real(dp), parameter :: CHARGE_TOL = 5.0e-8_dp
      !! Two different quadratures (terco's grid level 5 vs pyscf's
      !! 150x974) AND two different free-atom representations (terco's
      !! cubic-Hermite radial table vs pyscf's exact AOs) integrating the
      !! SAME partition formula, so this is not a self-consistency check.
      !! Measured (this file, with CHARGE_TOL set to 1e-12 to force the
      !! failure message to print the actual gap): water's worst atom
      !! differs by 1.30e-9, O2's (which should be exactly zero by
      !! symmetry) by 1.69e-9 -- both close to the 1.7e-9 pyscf itself
      !! could not do better than for O2 (the residual its own,
      !! unrenormalised, 150x974 grid leaves; hirshfeld.ref's Q_2). So the
      !! two quadratures and the two free-atom representations already
      !! agree far inside mqc's own CHARGE_TOL=1e-7 (test_mqc_hirshfeld.f90,
      !! a looser comparison of the SAME two things on a different
      !! molecule); 5e-8 leaves comfortable headroom above the measured
      !! 1.7e-9 without being so loose that a real regression (a wrong
      !! partition formula, a misread proatom table) would slip under it.
   real(dp), parameter :: SUM_TOL = 1.0e-9_dp
      !! The charges are renormalised to Tr(D S) (trc_hirshfeld's own
      !! `hirshfeld_renormalise`), so this is round-off, not quadrature.

contains

   subroutine collect_hirshfeld_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("free atom holds its electrons", test_free_atom_electrons), &
                  new_unittest("water matches mqc's recipe in pyscf", test_water), &
                  new_unittest("o2 triplet matches mqc's recipe in pyscf", test_o2), &
                  new_unittest("water matches mqc's OWN stored Q_NEUTRAL", test_mqc_water_neutral), &
                  new_unittest("cation matches mqc's OWN stored Q_CATION/S_CATION", test_mqc_water_cation), &
                  new_unittest("heavy element is refused", test_refused), &
                  new_unittest("table loader recovers and never keeps a stale file", test_loader) &
                  ]
   end subroutine collect_hirshfeld_tests

   subroutine test_free_atom_electrons(error)
      !! 4 pi int r^2 rho(r) dr = Z, on a radial grid uniform in ln r --
      !! the table itself, independent of any molecule or grid, for H and O
      !! (the elements the two molecular cases below need).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 4000
      type(error_t) :: err
      real(dp) :: r(N), rho(N), h, total
      integer :: k, iz
      integer, parameter :: ELEMENTS(2) = [1, 8]

      h = (log(80.0_dp) - log(1.0e-7_dp))/real(N - 1, dp)
      do k = 1, N
         r(k) = 1.0e-7_dp*exp(real(k - 1, dp)*h)
      end do
      do iz = 1, size(ELEMENTS)
         block
            real(dp), allocatable :: rho_a(:)
            call trc_hirshfeld_free_atom_density(ELEMENTS(iz), r, rho_a, err)
            call check(error,.not. err%has_error(), "free atom: "//err%get_message())
            if (allocated(error)) return
            rho = rho_a
         end block
         total = 0.0_dp
         do k = 1, N - 1
            total = total + 0.5_dp*(r(k + 1) - r(k))*(r(k)**2*rho(k) + r(k + 1)**2*rho(k + 1))
         end do
         total = 4.0_dp*3.14159265358979323846_dp*total
         call check(error, total, real(ELEMENTS(iz), dp), thr=5.0e-3_dp, &
                   message="free atom Z="//"integrates off its nuclear charge")
         if (allocated(error)) return
      end do
   end subroutine test_free_atom_electrons

   subroutine test_water(error)
      !! RHF/6-31G at mqc's own w1.xyz geometry (test/hirshfeld_water.xyz),
      !! against hirshfeld_reference.py's recipe run for this exact case.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: q(:)
      real(dp) :: worst
      integer :: ia, nat
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: q_ref(3), s_ref(3)

      call read_pyscf_case(1, 3, q_ref, s_ref)
      call read_xyz("hirshfeld_water.xyz", nat, z, r)
      call run_and_charge(z, r, "basis_sets/6-31g.json", 0, 0, q, error)
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         worst = max(worst, abs(q(ia) - q_ref(ia)))
      end do
      call check(error, worst < CHARGE_TOL, "water: worst |q_terco - q_pyscf| = "//to_text(worst)// &
                " exceeds CHARGE_TOL = "//to_text(CHARGE_TOL))
      if (allocated(error)) return
      call check(error, abs(sum(q)), 0.0_dp, thr=SUM_TOL, &
                message="water: charges do not sum to the molecular charge (0)")
   end subroutine test_water

   subroutine test_o2(error)
      !! UHF/cc-pVDZ triplet O2, against the same recipe.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: q(:)
      real(dp) :: worst
      integer :: ia, nat
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: sp(:)
      real(dp) :: q_ref(2), s_ref(2)

      call read_pyscf_case(2, 2, q_ref, s_ref)
      call read_xyz("o2.xyz", nat, z, r)
      call run_and_charge(z, r, "basis_sets/cc-pvdz.json", 0, 2, q, error, spin_q=sp)
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 2
         worst = max(worst, abs(q(ia) - q_ref(ia)), abs(sp(ia) - s_ref(ia)))
      end do
      call check(error, worst < CHARGE_TOL, "o2: worst |q or spin, terco - pyscf| = "//to_text(worst)// &
                " exceeds CHARGE_TOL = "//to_text(CHARGE_TOL))
      if (allocated(error)) return
      call check(error, abs(sum(q)), 0.0_dp, thr=SUM_TOL, &
                message="o2: charges do not sum to the molecular charge (0)")
   end subroutine test_o2

   subroutine test_mqc_water_neutral(error)
      !! RHF/cc-pVDZ at mqc's OWN water geometry (test_mqc_hirshfeld.f90's
      !! WATER, Bohr), against mqc's OWN stored Q_NEUTRAL -- mqc's Fortran
      !! result itself (its free atoms solved through its CPU HF/DFT
      !! backend), not a pyscf cross-check.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: q(:)
      real(dp) :: worst
      integer :: ia

      call run_and_charge(MQC_WATER_Z, MQC_WATER_R, "basis_sets/cc-pvdz.json", 0, 0, q, error)
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         worst = max(worst, abs(q(ia) - Q_NEUTRAL(ia)))
      end do
      call check(error, worst < MQC_CHARGE_TOL, "water: worst |q_terco - q_mqc| = "// &
                to_text(worst)//" exceeds mqc's own CHARGE_TOL = "//to_text(MQC_CHARGE_TOL))
      if (allocated(error)) return
      call check(error, abs(sum(q)), 0.0_dp, thr=SUM_TOL, &
                message="water: charges do not sum to the molecular charge (0)")
   end subroutine test_mqc_water_neutral

   subroutine test_mqc_water_cation(error)
      !! UHF/cc-pVDZ doublet, same geometry, against mqc's OWN stored
      !! Q_CATION and S_CATION (charges AND spin populations from the
      !! same weights, test_mqc_hirshfeld.f90's test_cation).
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: q(:), s(:)
      real(dp) :: worst_q, worst_s
      integer :: ia

      call run_and_charge(MQC_WATER_Z, MQC_WATER_R, "basis_sets/cc-pvdz.json", 1, 1, q, error, &
                          spin_q=s)
      if (allocated(error)) return

      worst_q = 0.0_dp
      worst_s = 0.0_dp
      do ia = 1, 3
         worst_q = max(worst_q, abs(q(ia) - Q_CATION(ia)))
         worst_s = max(worst_s, abs(s(ia) - S_CATION(ia)))
      end do
      call check(error, worst_q < MQC_CHARGE_TOL, "cation: worst |q_terco - q_mqc| = "// &
                to_text(worst_q)//" exceeds mqc's own CHARGE_TOL = "//to_text(MQC_CHARGE_TOL))
      if (allocated(error)) return
      call check(error, worst_s < MQC_CHARGE_TOL, "cation: worst |s_terco - s_mqc| = "// &
                to_text(worst_s)//" exceeds mqc's own CHARGE_TOL = "//to_text(MQC_CHARGE_TOL))
      if (allocated(error)) return
      call check(error, abs(sum(q) - 1.0_dp), 0.0_dp, thr=SUM_TOL, &
                message="cation: charges do not sum to the molecular charge (1)")
      if (allocated(error)) return
      call check(error, abs(sum(s) - 1.0_dp), 0.0_dp, thr=SUM_TOL, &
                message="cation: spin populations do not sum to one unpaired electron")
   end subroutine test_mqc_water_cation

   subroutine test_refused(error)
      !! Rubidium (Z=37) has no free atom (the recipe covers H-Kr): a
      !! one-atom "molecule" of it is refused, not silently given a zero
      !! charge.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: q(:), dzero(:, :)
      integer :: z(1)
      real(dp) :: r(3, 1)

      z = [37]
      r = 0.0_dp
      call trc_basis_from_json("basis_sets/def2-svp.json", 1, z, r, b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      allocate (dzero(b%nao, b%nao))
      dzero = 0.0_dp
      call trc_hirshfeld_populations(b, dzero, q, err)
      call check(error, err%has_error(), "Rb was not refused")
   end subroutine test_refused

   subroutine run_and_charge(z, r, basis_path, charge, n_unpaired, q, error, spin_q)
      ! Converge the SCF (RHF if n_unpaired==0, else UHF) and hand the
      ! total density, in terco's Cartesian AO space, to
      ! trc_hirshfeld_populations -- the same sph_to_cart step
      ! trc_capi's own trc_hirshfeld_charges makes (see its header). With
      ! `spin_q`, the alpha-minus-beta density goes through the same
      ! transform and the same partition weights, as mqc's own
      ! `hirshfeld_populations(..., spin_density=..., spin_populations=...)`.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: r(:, :)
      character(len=*), intent(in) :: basis_path
      integer, intent(in) :: charge, n_unpaired
      real(dp), allocatable, intent(out) :: q(:)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable, intent(out), optional :: spin_q(:)

      type(trc_basis_t) :: b
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer :: nat, nalpha, nbeta, nelec
      real(dp), allocatable :: d_scf(:, :), d_cart(:, :)
      real(dp), allocatable :: s_scf(:, :), s_cart(:, :)
      logical :: want_spin

      nat = size(z)
      call trc_basis_from_json(basis_path, nat, z, r, b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return

      nelec = sum(z) - charge
      nbeta = (nelec - n_unpaired)/2
      nalpha = nelec - nbeta
      opts%unrestricted = n_unpaired /= 0
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 200
      call trc_scf_run(b, nalpha, nbeta, opts, res)
      call check(error, res%converged, "reference SCF did not converge")
      if (allocated(error)) return

      want_spin = present(spin_q) .and. res%nspin == 2
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

      if (want_spin) then
         allocate (s_scf(size(res%dmat, 1), size(res%dmat, 2)))
         s_scf = res%dmat(:, :, 1) - res%dmat(:, :, 2)
         if (b%spherical) then
            allocate (s_cart(b%nao, b%nao))
            call sph_sph_to_cart(b, s_scf, s_cart)
         else
            s_cart = s_scf
         end if
         call trc_hirshfeld_populations(b, d_cart, q, err, spin_density=s_cart, &
                                        spin_populations=spin_q)
      else
         call trc_hirshfeld_populations(b, d_cart, q, err)
      end if
      call check(error,.not. err%has_error(), "hirshfeld: "//err%get_message())
   end subroutine run_and_charge

   function to_text(x) result(s)
      real(dp), intent(in) :: x
      character(len=32) :: s
      write (s, "(es12.4)") x
   end function to_text

   subroutine test_loader(error)
      !! The proatom table loader, through the real module cache: a path
      !! that does not exist is ERROR_IO and leaves nothing half-loaded, so
      !! the default file still loads afterwards; and once a file IS
      !! loaded, asking for a different one reloads (here: fails) instead
      !! of silently answering from the old one. One H atom in STO-3G with
      !! D = 1/S_11 = 1 is enough -- its charge is 0 by construction.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: b
      type(error_t) :: err
      real(dp), allocatable :: q(:), d(:, :)
      integer :: z(1)
      real(dp) :: r(3, 1)

      z = [1]
      r = 0.0_dp
      call trc_basis_from_json("basis_sets/sto-3g.json", 1, z, r, b, err)
      call check(error,.not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      allocate (d(b%nao, b%nao))
      d = 1.0_dp

      call trc_hirshfeld_populations(b, d, q, err, data_path="no/such/proatoms.json")
      call check(error, err%has_error() .and. err%get_code() == ERROR_IO, &
                 "a missing table must be ERROR_IO")
      if (allocated(error)) return
      call err%clear()

      call trc_hirshfeld_populations(b, d, q, err)
      call check(error,.not. err%has_error(), "the default table must load after a failed one: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, abs(q(1)), 0.0_dp, thr=SUM_TOL, message="a lone H atom must carry no charge")
      if (allocated(error)) return

      call trc_hirshfeld_populations(b, d, q, err, data_path="no/such/other.json")
      call check(error, err%has_error(), "a different, missing table was answered from the cached one")
      if (allocated(error)) return
      call err%clear()

      call trc_hirshfeld_populations(b, d, q, err)
      call check(error,.not. err%has_error(), "the default table must reload: "//err%get_message())
   end subroutine test_loader

   subroutine read_pyscf_case(ic, natm, q_ref, s_ref)
      ! Case `ic` of hirshfeld.ref (generate.py hirshfeld): its charges and
      ! spin populations (zero for a closed shell). The file is read in
      ! order, so every case before `ic` is read and dropped.
      integer, intent(in) :: ic, natm
      real(dp), intent(out) :: q_ref(natm), s_ref(natm)
      character(len=64) :: name
      real(dp), allocatable :: v(:)
      integer :: u, ncase, k, nat

      s_ref = 0.0_dp
      call ref_require("hirshfeld.ref", u)
      call ref_read_int(u, "NCASE", ncase)
      do k = 1, min(ic, ncase)
         call ref_read_str(u, "NAME_"//to_string(k), name)
         call ref_read_int(u, "NATM_"//to_string(k), nat)
         allocate (v(nat))
         call ref_read_vector(u, "Q_"//to_string(k), nat, v)
         if (k == ic) q_ref = v
         call ref_read_vector(u, "S_"//to_string(k), nat, v)
         if (k == ic) s_ref = v
         deallocate (v)
      end do
      call ref_close(u)
   end subroutine read_pyscf_case

end module check_hirshfeld_suite

program check_hirshfeld
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_hirshfeld_suite, only: collect_hirshfeld_tests
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat
   stat = 0
   testsuites = [new_testsuite("hirshfeld", collect_hirshfeld_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0,a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_hirshfeld
