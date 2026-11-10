!
! The SCF with the SKELETON Fock build (S14.3, docs/SYMMETRY_FOCK.md): the
! four-centre engine built with the point group, only the orbit
! representatives launched, G symmetrized after the fold.
!
! THE CENTRAL ORACLE, ONE MOLECULE AT A TIME, THREE RUNS:
!   skeleton   trc_sym_scf_t%skeleton = SYM_SKEL_ON
!   S4         skeleton = SYM_SKEL_OFF: the full Fock build, then F symmetrized
!   C1         no group
! The total energy (<= 1e-10 Eh HF, <= 1e-9 DFT), the orbital energies and
! the density of the skeleton run must match both of the others. A silent
! fallback to the full build fails two ways: the engine of the skeleton run
! must hold the group's orbits, and its last build (`eri%nwork`, the
! quartets the plan launched) must be smaller than the S4 run's, which is
! measured on the same engine type with the group left out.
!
! Cases cover every group the SCF supports (Cs, C2, Ci, C2v, C2h, D2, D2h and
! the rotated tops), RHF, UHF, a global hybrid (B3LYP), a range-separated
! hybrid (wB97X: the long-range exchange pass) and the incremental Fock build
! on and off. The molecules are those of the group's own check_sym_<group>.
!
module check_sym_scf_skeleton_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run, trc_pair_eri_cache_t
   use trc_error, only: error_t
   use trc_sym_ao, only: sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_run, only: SYM_SKEL_ON, SYM_SKEL_OFF
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_test_basis, only: read_xyz
   implicit none
   private
   public :: collect_sym_scf_skeleton

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the oracle's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the oracle's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density
   integer, parameter :: GRID_LEVEL = 3           !! DFT grid: the three runs share it, nothing is compared to a code

contains

   subroutine collect_sym_scf_skeleton(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Cs: RHF HOCl cc-pVDZ", test_cs_rhf), &
                  new_unittest("Cs: UHF HCO cc-pVDZ", test_cs_uhf), &
                  new_unittest("Cs: B3LYP HOCl cc-pVDZ", test_cs_b3lyp), &
                  new_unittest("Cs: wB97X (range-separated) HOCl 6-31G*", test_cs_wb97x), &
                  new_unittest("Cs: UHF wB97X HCO 6-31G* (several densities, long-range pass)", test_cs_uwb97x), &
                  new_unittest("Cs: RHF CHFCl2 6-31G*, incremental off", test_cs_full_builds), &
                  new_unittest("C2: RHF H2O2 cc-pVDZ", test_c2_rhf), &
                  new_unittest("Ci: RHF N2H2F2 6-31G*", test_ci_rhf), &
                  new_unittest("C2v: RHF water cc-pVDZ", test_c2v_rhf), &
                  new_unittest("C2v: UHF NH2 cc-pVDZ", test_c2v_uhf), &
                  new_unittest("C2v: B3LYP water cc-pVDZ", test_c2v_b3lyp), &
                  new_unittest("C2v: RHF CH2Cl2 6-31G*, incremental off", test_c2v_full_builds), &
                  new_unittest("C2h: RHF trans-N2H2 cc-pVDZ", test_c2h_rhf), &
                  new_unittest("D2: RHF twisted C2H4 cc-pVDZ", test_d2_rhf), &
                  new_unittest("D2h: RHF C2H4 cc-pVDZ", test_d2h_rhf), &
                  new_unittest("D2h: UHF C2H4+ cc-pVDZ", test_d2h_uhf), &
                  new_unittest("D2h: RHF B2H6 (rotated input) 6-31G*", test_d2h_b2h6), &
                  new_unittest("D2h: wB97X C2H4 6-31G*", test_d2h_wb97x), &
                  new_unittest("tops: RHF CH4 (rotated, D2 subgroup of Td) 6-31G*", test_tops_ch4), &
                  new_unittest("tops: UHF NH3+ (Cs) 6-31G*", test_tops_nh3_uhf), &
                  new_unittest("tops: RHF N2 (rotated, D2h subgroup of D-inf-h) 6-31G*", test_tops_n2) &
                  ]
   end subroutine collect_sym_scf_skeleton

   subroutine oracle(error, label, xyz, basis_json, spherical, functional, charge, mult, group, requested, incremental, min_gain)
      !! The three runs of one molecule, compared; `group` is the group the
      !! run must be in, `requested` (may be empty) the subgroup to ask for.
      !! `min_gain`: the least ratio of the full build's quartets to the
      !! skeleton's. A molecule whose atoms all lie on the mirror plane (every
      !! atom its own image) has trivial shell-pair orbits under Cs, so its
      !! skeleton plan is no smaller than the full one (docs/SYMMETRY_FOCK.md,
      !! planar molecules gain |G|/2 = 1, here at a few percent of overhead
      !! for the plan's rep x bin shape); those pass 0.9. Every other case
      !! passes 1.2 (measured 1.35 to 3.35).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label, xyz, basis_json, functional, group, requested
      logical, intent(in) :: spherical, incremental
      integer, intent(in) :: charge, mult
      real(dp), intent(in) :: min_gain

      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_scf_t) :: sk, s4
      type(trc_scf_result_t) :: rk, r4, r1
      type(trc_pair_eri_cache_t) :: ck, c4, c1
      type(trc_scf_options_t) :: opts
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      integer :: nat, nelec, nalpha, nbeta, dev
      integer(int64) :: nw_skel, nw_full
      real(dp) :: tol_e

      dev = trc_bind_device(0)
      tol_e = merge(TOL_E_DFT, TOL_E_HF, len_trim(functional) > 0)
      call read_xyz(xyz, nat, elem, r)
      if (len_trim(requested) > 0) then
         geom = sym_geom_detect(r, elem, err, basis_id=elem, requested=requested)
      else
         geom = sym_geom_detect(r, elem, err, basis_id=elem)
      end if
      call check(error, .not. err%has_error(), label//": sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, group, label//": group of "//xyz)
      if (allocated(error)) return
      call trc_basis_from_json(basis_json, nat, elem, geom%r_std, bas, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), label//": basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = spherical
      call bas%to_device()

      nelec = sum(elem) - charge
      nalpha = (nelec + mult - 1)/2
      nbeta = nelec - nalpha
      opts%functional = functional
      opts%grid_level = GRID_LEVEL
      opts%incremental = incremental
      ! Above the accumulation-order noise floor (check_sym_scf_cs); the
      ! energy is variational, so the dE gates are unaffected.
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 3.0e-9_dp
      opts%max_iteration = 100

      sk%ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), label//": sym_ao_build: "//err%get_message())
      if (allocated(error)) then
         call bas%release()
         return
      end if
      s4%ao = sk%ao
      ! Both runs are FORCED. The default (SYM_SKEL_AUTO) would let
      ! sym_run_skeleton_size_rule decline the skeleton build on most of these
      ! molecules -- they are small and compactly orbited on purpose, one per
      ! group -- and this oracle is about the two builds agreeing, not about
      ! which one a user gets.
      sk%skeleton = SYM_SKEL_ON
      s4%skeleton = SYM_SKEL_OFF

      call trc_scf_run(bas, nalpha, nbeta, opts, rk, sym=sk, pair_eri_cache=ck)
      call check(error, rk%converged, label//": the skeleton SCF did not converge: "//trim(rk%message))
      if (.not. allocated(error)) then
         nw_skel = ck%eri%nwork
         call check(error, ck%eri%orb%built .and. ck%orbits_built, &
                    label//": the skeleton run's engine was not built with the point group (silent fallback)")
      end if
      if (.not. allocated(error)) then
         call trc_scf_run(bas, nalpha, nbeta, opts, r4, sym=s4, pair_eri_cache=c4)
         call check(error, r4%converged, label//": the S4 SCF did not converge: "//trim(r4%message))
      end if
      if (.not. allocated(error)) then
         nw_full = c4%eri%nwork
         call check(error, .not. c4%eri%orb%built .and. .not. c4%orbits_built, &
                    label//": skeleton = .false. must leave the engine without the group")
      end if
      if (.not. allocated(error)) then
         call trc_scf_run(bas, nalpha, nbeta, opts, r1, pair_eri_cache=c1)
         call check(error, r1%converged, label//": the C1 SCF did not converge: "//trim(r1%message))
      end if
      ! The flag off is the S4 path: the engine and its launch plan are the
      ! C1 run's own (the plan is the plain one, only F is symmetrized).
      if (.not. allocated(error)) &
         call check(error, nw_full == c1%eri%nwork, label//": skeleton = .false. launched "//to_string(nw_full)// &
                    " quartets, the C1 run "//to_string(c1%eri%nwork))
      if (.not. allocated(error)) then
         print '(a,a,a,i0,a,i0,a,f5.2,a)', '  ', label, ': last build quartets skeleton ', nw_skel, ' full ', nw_full, &
            '  (x', real(nw_full, dp)/real(max(nw_skel, 1_int64), dp), ')'
         call check(error, nw_skel /= nw_full .and. real(nw_full, dp) >= min_gain*real(nw_skel, dp), &
                    label//": the skeleton build launched "//to_string(nw_skel)//" quartets against the full build's "// &
                    to_string(nw_full)//": expected a gain of at least "//to_string(min_gain, "f4.2"))
      end if
      if (.not. allocated(error)) call compare(error, label//" skeleton vs C1", rk, r1, tol_e)
      if (.not. allocated(error)) call compare(error, label//" skeleton vs S4", rk, r4, tol_e)
      if (.not. allocated(error) .and. incremental) &
         call check(error, rk%n_increments > 0, label//": the run must have used incremental Fock builds")
      if (.not. allocated(error) .and. .not. incremental) &
         call check(error, rk%n_increments == 0, label//": incremental builds were off")
      call ck%release()
      call c4%release()
      call c1%release()
      call bas%release()
   end subroutine oracle

   subroutine compare(error, label, a, b, tol_e)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      type(trc_scf_result_t), intent(in) :: a, b
      real(dp), intent(in) :: tol_e
      real(dp) :: de, deps, dd
      integer :: nmo

      de = abs(a%energy - b%energy)
      nmo = size(a%eps, 1) - a%n_lindep
      deps = maxval(abs(a%eps(1:nmo, :) - b%eps(1:nmo, :)))
      dd = maxval(abs(a%dmat - b%dmat))
      print '(a,a,a,es10.2,a,es10.2,a,es10.2,a,i0,a,i0,a)', '  ', label, ': |dE| ', de, '  max|d eps| ', deps, &
         '  max|dD| ', dd, '  iterations ', a%iterations, ' (', b%iterations, ')'
      call check(error, de <= tol_e, label//": |dE| = "//to_string(de, "es10.2")//" exceeds "//to_string(tol_e, "es10.2"))
      if (allocated(error)) return
      call check(error, deps <= TOL_EPS, label//": orbital energies differ by "//to_string(deps, "es10.2"))
      if (allocated(error)) return
      call check(error, dd <= TOL_D, label//": densities differ by "//to_string(dd, "es10.2"))
   end subroutine compare

   subroutine test_cs_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Cs RHF', 'hocl.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, 'Cs', '', .true., 0.9_dp)
   end subroutine test_cs_rhf

   subroutine test_cs_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Cs UHF', 'hco.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 2, 'Cs', '', .true., 0.9_dp)
   end subroutine test_cs_uhf

   subroutine test_cs_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Cs B3LYP', 'hocl.xyz', 'basis_sets/cc-pvdz.json', .true., 'b3lyp', 0, 1, 'Cs', '', .true., 0.9_dp)
   end subroutine test_cs_b3lyp

   subroutine test_cs_wb97x(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Cs wB97X', 'hocl.xyz', 'basis_sets/6-31gs.json', .false., 'wb97x', 0, 1, 'Cs', '', .true., 0.9_dp)
   end subroutine test_cs_wb97x

   subroutine test_cs_uwb97x(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Cs UwB97X', 'hco.xyz', 'basis_sets/6-31gs.json', .false., 'wb97x', 0, 2, 'Cs', '', .true., 0.9_dp)
   end subroutine test_cs_uwb97x

   subroutine test_cs_full_builds(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Cs RHF full builds', 'chfcl2.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, 'Cs', '', .false., 1.2_dp)
   end subroutine test_cs_full_builds

   subroutine test_c2_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'C2 RHF', 'h2o2.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, 'C2', '', .true., 1.2_dp)
   end subroutine test_c2_rhf

   subroutine test_ci_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'Ci RHF', 'n2h2f2_ci.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, 'Ci', '', .true., 1.2_dp)
   end subroutine test_ci_rhf

   subroutine test_c2v_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'C2v RHF', 'water.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, 'C2v', '', .true., 1.2_dp)
   end subroutine test_c2v_rhf

   subroutine test_c2v_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'C2v UHF', 'nh2.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 2, 'C2v', '', .true., 1.2_dp)
   end subroutine test_c2v_uhf

   subroutine test_c2v_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'C2v B3LYP', 'water.xyz', 'basis_sets/cc-pvdz.json', .true., 'b3lyp', 0, 1, 'C2v', '', .true., 1.2_dp)
   end subroutine test_c2v_b3lyp

   subroutine test_c2v_full_builds(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'C2v RHF full builds', 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, 'C2v', '', &
         .false., 1.2_dp)
   end subroutine test_c2v_full_builds

   subroutine test_c2h_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'C2h RHF', 'n2h2_trans.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, 'C2h', '', .true., 1.2_dp)
   end subroutine test_c2h_rhf

   subroutine test_d2_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'D2 RHF', 'c2h4_twisted.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, 'D2', '', .true., 1.2_dp)
   end subroutine test_d2_rhf

   subroutine test_d2h_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'D2h RHF', 'c2h4_d2h.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, 'D2h', '', .true., 1.2_dp)
   end subroutine test_d2h_rhf

   subroutine test_d2h_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'D2h UHF', 'c2h4_d2h.xyz', 'basis_sets/cc-pvdz.json', .true., '', 1, 2, 'D2h', '', .true., 1.2_dp)
   end subroutine test_d2h_uhf

   subroutine test_d2h_b2h6(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'D2h B2H6', 'b2h6_d2h_rot.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, 'D2h', '', .true., 1.2_dp)
   end subroutine test_d2h_b2h6

   subroutine test_d2h_wb97x(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'D2h wB97X', 'c2h4_d2h.xyz', 'basis_sets/6-31gs.json', .false., 'wb97x', 0, 1, 'D2h', '', .true., 1.2_dp)
   end subroutine test_d2h_wb97x

   subroutine test_tops_ch4(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'tops CH4', 'tops_ch4_rot.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, 'D2', 'D2', .true., 1.2_dp)
   end subroutine test_tops_ch4

   subroutine test_tops_nh3_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'tops NH3+', 'tops_nh3_rot.xyz', 'basis_sets/6-31gs.json', .false., '', 1, 2, 'Cs', '', .true., 1.2_dp)
   end subroutine test_tops_nh3_uhf

   subroutine test_tops_n2(error)
      type(error_type), allocatable, intent(out) :: error
      call oracle(error, 'tops N2', 'tops_n2_rot.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, 'D2h', '', .true., 1.2_dp)
   end subroutine test_tops_n2

end module check_sym_scf_skeleton_suite

program check_sym_scf_skeleton
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_sym_scf_skeleton_suite, only: collect_sym_scf_skeleton
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, is
   stat = 0
   testsuites = [new_testsuite("sym_scf_skeleton", collect_sym_scf_skeleton)]
   do is = 1, size(testsuites)
      write (error_unit, '(a)') "Testing: "//testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
   end do
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_scf_skeleton
