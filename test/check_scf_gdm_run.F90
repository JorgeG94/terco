!
! Geometric direct minimisation as a step mode of the SCF loop
! (src/scf/trc_scf_gdm_step.F90), end to end through `trc_scf_run`, as a
! test-drive suite. check_scf_gdm.F90 covers the numerical core (L-BFGS,
! preconditioner, step control) on its own; this is the loop around it.
!
! Most cases run GDM at the DEFAULT `gdm_start` (DIIS first, GDM once the
! commutator is below 0.1, which is how it is meant to be used); the water RHF,
! OH UHF and the reviewer's started-too-early case run it from the first plain
! diagonalisation (`gdm_start = huge`), where the guess is good (or, for the
! last, where the failure is the point). Each case asserts `gdm_steps > 0`, so a
! run that never left DIIS cannot pass, and that the accepted energies rose
! by no more than the accept rule allows: the plain slack (1e-10) between two
! full Fock builds, 1e-8 where an incremental build is involved
! (`check_gdm_ran`). Cases that say "exact" switch the incremental build off,
! which makes the energy deterministic and that bound sharp.
!
! References, and where each comes from:
!
!  - water RHF and OH UHF, cc-pVDZ: the pyscf energies check_spherical.F90
!    and check_scf_stabilise.F90 hold (cart=False, conv_tol 1e-11), 1e-8;
!  - RKS (B3LYP, water) and UKS (PBE, OH), 6-31G uncontracted: rks.ref
!    (pyscf on terco's own grid), at scf_rks's tolerances: 1e-8 closed shell;
!    2e-6 for the open-shell DFT converged energy, which scf_rks justifies by
!    the energy AT pyscf's density agreeing to 1e-8 (OH's pi pair is
!    degenerate and the orientation costs grid error): this check replicates
!    that at-density run and its 1e-8, and starts from the reference's own
!    symmetry-broken density as scf_rks does. The closed-shell case starts
!    from the converged Hartree-Fock density of the same basis;
!  - PCM (C-PCM, 6-31G*, water RHF and OH UHF): pcm_scf_main.ref, 1e-8,
!    symmetric, as scf_pcm.F90 holds it. OH is started from pyscf's own
!    converged density (scf_pcm's reason: its pi pair is split only by the
!    cavity and the energy along that mode is flat to 3e-8), with scf_pcm's
!    at-density check; a second case runs OH from the default guess with
!    exact energies as a regression of the slow flat mode (iteration ceiling);
!  - density fitting (water cc-pVDZ HF): scf_df.ref, 1e-8;
!  - a linearly dependent basis (O2- d-aug-cc-pVDZ, 5 functions dropped):
!    check_scf_lindep_gap's pyscf energy, 1e-8;
!  - water B3LYP/6-31G*: cpks_water_6-31gs_b3lyp.ref's E_SCF, 1e-8 (the
!    default-start run and the DIIS run; the started-too-early run is the
!    case that must NOT reach it).
!
! The mechanism cases (rejection, breakdown, the iteration ceiling, a rotation
! cap, MOM) have no external reference for the mechanism, only for the energy
! the run must still reach, or for the DIIS/MOM run of the same input.
!
module check_scf_gdm_run_suite
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_scf_gdm, only: scf_gdm_scale_min
   use trc_scf_gdm_step, only: scf_gdm_step_accept_slack, GDM_NOISE_BOUND
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz, uncontract_basis
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_real, ref_read_int, ref_read_matrix3
   implicit none
   private
   public :: collect_scf_gdm_run

   ! pyscf, cart=False, conv_tol 1e-11 (check_spherical.F90 / check_scf_stabilise.F90).
   real(dp), parameter :: E_WATER_RHF = -76.026772053394_dp
   real(dp), parameter :: E_OH_UHF = -75.393846033475_dp
   ! pyscf, check_scf_lindep_gap.F90.
   real(dp), parameter :: E_O2M_LINDEP = -74.147460338947_dp
   real(dp), parameter :: TOL_REF = 1.0e-8_dp
   real(dp), parameter :: EPS_WATER = 78.3553_dp
   character(len=*), parameter :: GDM_MESSAGE_PART = "no downhill step"
   integer, parameter :: MAX_ITER = 300
   real(dp), parameter :: SLACK = 1.0e-10_dp
   !! The default `gdm_slack`: the most an accepted point may sit above the
   !! one before when both Fock builds were full.
   real(dp), parameter :: NOISE_BOUND = 1.0e-8_dp
   !! The wider bound of src/scf/trc_scf_gdm_step.F90 (GDM_NOISE_BOUND), for a
   !! comparison in which an incremental build is involved.
   integer, parameter :: FLAT_CEILING = 250
   !! Iterations the flat-mode regression may take (set after measuring).
   real(dp), parameter :: ROUNDOFF = 1.0e-13_dp
   !! E_try <= E_acc + slack is exact in floating point; E_try - E_acc can
   !! still exceed the slack by an ulp of 76 Eh (1.4e-14).

contains

   subroutine collect_scf_gdm_run(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water RHF cc-pVDZ: GDM from iteration 2 reaches pyscf", test_water_rhf), &
                  new_unittest("OH UHF cc-pVDZ: GDM reaches pyscf, both spins", test_oh_uhf), &
                  new_unittest("water RKS B3LYP: GDM reaches pyscf (rks.ref)", test_water_rks), &
                  new_unittest("OH UKS PBE: GDM reaches pyscf (rks.ref), at-density control", test_oh_uks), &
                  new_unittest("water C-PCM RHF: GDM reaches pyscf (pcm_scf_main.ref)", test_water_pcm), &
                  new_unittest("OH C-PCM UHF from pyscf's density: GDM reaches pyscf at 1e-8", test_oh_pcm), &
                  new_unittest("OH C-PCM UHF flat mode from the guess: iteration ceiling, exact energies", &
                               test_oh_pcm_flat_mode), &
                  new_unittest("water density-fitted RHF: GDM reaches pyscf (scf_df.ref)", test_water_df), &
                  new_unittest("O2- with 5 functions dropped: GDM reaches pyscf, null states stay lifted", &
                               test_lindep), &
                  new_unittest("default handover: DIIS first, then GDM, same energy", test_default_handover), &
                  new_unittest("accept slack: plain between full builds, wide only when an incremental build is involved", &
                               test_accept_slack), &
                  new_unittest("adiis before the handover: GDM then reaches the same energy", test_adiis_handover), &
                  new_unittest("a too-large maximum rotation: steps are rejected, accepted energies never rise", &
                               test_rejected_steps), &
                  new_unittest("H atom UHF: an empty spin block is diagonalised whole, not rotated", test_empty_spin), &
                  new_unittest("H atom in a one-function basis: nothing to rotate in either spin", test_nothing_to_rotate), &
                  new_unittest("the maximum rotation caps every step", test_rotation_cap), &
                  new_unittest("returned orbitals are semicanonical and consistent with the density", &
                               test_orbitals_consistent), &
                  new_unittest("gdm with frac_occ or soscf is refused", test_refusals), &
                  new_unittest("a rejected trial never converges the run, however loose the thresholds", &
                               test_rejected_never_converges), &
                  new_unittest("no downhill step: the run stops after one round, not converged, with the GDM message", &
                               test_no_downhill), &
                  new_unittest("max_iteration reached on a rejected trial: the accepted point is returned", &
                               test_ceiling_on_rejected), &
                  new_unittest("an unusable maximum rotation stops the run on the accepted point", test_invalid_rotation), &
                  new_unittest("mom + gdm holds the excited state mom holds", test_mom_gdm), &
                  new_unittest("gdm started on the first diagonalisation can end non-aufbau and says so", &
                               test_started_too_early) &
                  ]
   end subroutine collect_scf_gdm_run

   ! --- helpers -----------------------------------------------------------

   function base_options() result(opts)
      !! The options every case starts from: tight enough to compare with
      !! pyscf at 1e-8, a generous iteration budget. GDM off.
      type(trc_scf_options_t) :: opts
      character(len=8) :: v
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-7_dp
      opts%max_iteration = MAX_ITER
      call get_environment_variable("TRC_GDM_TEST_VERBOSE", v)
      opts%verbose = len_trim(v) > 0
   end function base_options

   function exact_options() result(opts)
      !! `base_options` with the incremental Fock build off: every energy is a
      !! full build's, deterministic, so the plain accept slack is the bound.
      type(trc_scf_options_t) :: opts
      opts = base_options()
      opts%incremental = .false.
   end function exact_options

   function gdm_options(base, early) result(opts)
      !! `base` with GDM on, at the default `gdm_start`; with `early` from the
      !! first iteration that may (`gdm_start` huge: after the damped start,
      !! one plain diagonalisation).
      type(trc_scf_options_t), intent(in) :: base
      logical, intent(in), optional :: early
      type(trc_scf_options_t) :: opts
      opts = base
      opts%gdm = .true.
      if (present(early)) then
         if (early) opts%gdm_start = huge(1.0_dp)
      end if
   end function gdm_options

   subroutine make_basis(xyz, basis, na, nb, bas, error)
      !! Basis from `basis_sets/<basis>.json` on `<xyz>`; an odd electron count
      !! puts the extra one in alpha.
      character(len=*), intent(in) :: xyz, basis
      integer, intent(out) :: na, nb
      type(trc_basis_t), intent(out) :: bas
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm
      call read_xyz(xyz, natm, z, at_r)
      na = (sum(z) + 1)/2
      nb = sum(z) - na
      call trc_basis_from_json('basis_sets/'//basis//'.json', natm, z, at_r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
   end subroutine make_basis

   subroutine make_uncontracted(xyz, bas, na, nb, error)
      !! scf_rks's basis: 6-31G with every primitive its own shell (see that
      !! check's header for why pyscf can be compared with at 1e-8 only so).
      character(len=*), intent(in) :: xyz
      type(trc_basis_t), intent(out) :: bas
      integer, intent(out) :: na, nb
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: jbas
      type(error_t) :: err
      integer, allocatable :: sh_l(:), sh_np(:), zint(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), at_z(:), at_r(:, :)
      integer :: natm, nsh, maxnp
      call read_xyz(xyz, natm, zint, at_r)
      na = (sum(zint) + 1)/2
      nb = sum(zint) - na
      allocate (at_z(natm))
      at_z = real(zint, dp)
      call trc_basis_from_json("basis_sets/6-31g.json", natm, zint, at_r, jbas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call uncontract_basis(jbas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
      call bas%to_device()
   end subroutine make_uncontracted

   subroutine ref_seek(u, tag, found)
      !! Position unit `u` (a pyscf reference opened by `ref_open`) at the
      !! line that starts with `tag`, so `ref_read_*` of that tag reads it.
      !! The readers are sequential and a case needs one record out of many.
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      logical, intent(out) :: found
      character(len=1024) :: line
      integer :: ios
      found = .false.
      rewind (u)
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) return
         if (line(1:len(tag) + 1) == tag//' ') then
            backspace (u)
            found = .true.
            return
         end if
      end do
   end subroutine ref_seek

   subroutine ref_energy(file, tag, e, error)
      character(len=*), intent(in) :: file, tag
      real(dp), intent(out) :: e
      type(error_type), allocatable, intent(out) :: error
      integer :: u
      logical :: ok, found
      e = 0.0_dp
      call ref_open(file, u, ok)
      call check(error, ok, "missing pyscf reference "//file)
      if (allocated(error)) return
      call ref_seek(u, tag, found)
      call check(error, found, file//" has no record "//tag)
      if (.not. allocated(error)) call ref_read_real(u, tag, e)
      call ref_close(u)
   end subroutine ref_energy

   subroutine ref_density(file, tag, nao, nspin, d, error)
      character(len=*), intent(in) :: file, tag
      integer, intent(in) :: nao, nspin
      real(dp), intent(out) :: d(nao, nao, nspin)
      type(error_type), allocatable, intent(out) :: error
      integer :: u
      logical :: ok, found
      d = 0.0_dp
      call ref_open(file, u, ok)
      call check(error, ok, "missing pyscf reference "//file)
      if (allocated(error)) return
      call ref_seek(u, tag, found)
      call check(error, found, file//" has no record "//tag)
      if (.not. allocated(error)) call ref_read_matrix3(u, tag, nao, nao, nspin, d)
      call ref_close(u)
   end subroutine ref_density

   subroutine check_gdm_ran(error, res, what, exact, off_aufbau_ok)
      !! The run converged, took GDM steps, ended on an aufbau state (unless
      !! `off_aufbau_ok`), and the accepted energies rose by no more than the
      !! accept rule allows: the plain slack between full builds always, and
      !! overall that slack too when the run was `exact` (no incremental
      !! build), 1e-8 otherwise.
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_result_t), intent(in) :: res
      character(len=*), intent(in) :: what
      logical, intent(in), optional :: exact, off_aufbau_ok
      logical :: sharp, off_ok
      sharp = .false.
      off_ok = .false.
      if (present(exact)) sharp = exact
      if (present(off_aufbau_ok)) off_ok = off_aufbau_ok
      print "(a,a,a,i0,a,i0,a,i0,a,es9.2,a,es9.2,a,f20.12)", "    ", what, ": iterations ", res%iterations, "  gdm_steps ", &
         res%gdm_steps, "  rejected ", res%gdm_rejected, "  max rise ", res%gdm_max_rise, "  (full) ", res%gdm_max_rise_full, &
         "  E = ", res%energy
      call check(error, res%converged, what//": not converged ("//trim(res%message)//")")
      if (allocated(error)) return
      call check(error, res%gdm_steps > 0, what//": GDM never took a step ("//to_string(res%gdm_steps)//")")
      if (allocated(error)) return
      call check(error, res%gdm_max_rise_full <= SLACK + ROUNDOFF, what// &
                 ": two full-build points rose by "//to_string(res%gdm_max_rise_full))
      if (allocated(error)) return
      call check(error, res%gdm_max_rise <= merge(SLACK, NOISE_BOUND, sharp) + ROUNDOFF, what// &
                 ": an accepted point was higher than the last by "//to_string(res%gdm_max_rise))
      if (allocated(error)) return
      if (.not. off_ok) call check(error, res%gdm_off_aufbau == 0, what// &
                                   ": ended on a non-aufbau state ("//to_string(res%gdm_off_aufbau)//")")
      if (allocated(error)) return
   end subroutine check_gdm_ran

   ! --- Hartree-Fock ------------------------------------------------------

   subroutine test_water_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_result_t) :: res
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      call trc_scf_run(bas, na, nb, gdm_options(exact_options(), early=.true.), res)
      call bas%release()
      call check_gdm_ran(error, res, "water RHF", exact=.true.)
      if (allocated(error)) return
      call check(error, abs(res%energy - E_WATER_RHF) < TOL_REF, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - E_WATER_RHF)))
   end subroutine test_water_rhf

   subroutine test_oh_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_result_t) :: res
      integer :: na, nb
      call make_basis('oh.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      call trc_scf_run(bas, na, nb, gdm_options(exact_options(), early=.true.), res)
      call bas%release()
      call check_gdm_ran(error, res, "OH UHF", exact=.true.)
      if (allocated(error)) return
      call check(error, res%nspin == 2, "expected an unrestricted run")
      if (allocated(error)) return
      call check(error, abs(res%energy - E_OH_UHF) < TOL_REF, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - E_OH_UHF)))
   end subroutine test_oh_uhf

   ! --- Kohn-Sham ---------------------------------------------------------

   subroutine test_water_rks(error)
      !! rks.ref case 4: water, B3LYP, restricted, from the converged
      !! Hartree-Fock density of the same basis (a sane guess: the default GWH
      !! guess of this uncontracted basis starts 10000 Eh high, a poor test of
      !! anything but the guess).
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: base, hf_opts
      type(trc_scf_result_t) :: res, plain, hf
      real(dp) :: e_ref
      integer :: na, nb
      call ref_energy("rks.ref", "E_4", e_ref, error)
      if (allocated(error)) return
      call make_uncontracted('water.xyz', bas, na, nb, error)
      if (allocated(error)) return
      hf_opts = base_options()
      call trc_scf_run(bas, na, nb, hf_opts, hf)
      call check(error, hf%converged, "the Hartree-Fock guess did not converge: "//trim(hf%message))
      if (allocated(error)) then
         call bas%release()
         return
      end if
      base = base_options()
      base%functional = "b3lyp"
      base%grid_level = 3
      base%conv_diis = 1.0e-8_dp
      call trc_scf_run(bas, na, nb, gdm_options(base), res, dguess=hf%dmat)
      call trc_scf_run(bas, na, nb, base, plain, dguess=hf%dmat)
      call bas%release()
      call check_gdm_ran(error, res, "water B3LYP")
      if (allocated(error)) return
      call check(error, res%iterations < 40, "water B3LYP took "//to_string(res%iterations)//" iterations from a sane guess")
      if (allocated(error)) return
      call check(error, abs(res%energy - e_ref) < TOL_REF, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - e_ref)))
      if (allocated(error)) return
      call check(error, abs(res%energy - plain%energy) < 1.0e-8_dp, "GDM and DIIS energies differ by "// &
                 to_string(abs(res%energy - plain%energy)))
   end subroutine test_water_rks

   subroutine test_oh_uks(error)
      !! rks.ref case 10: OH, PBE, unrestricted, from the reference's own
      !! symmetry-broken density (scf_rks's reason: the degenerate pi pair).
      !! The converged energy gets scf_rks's 2e-6 for open-shell DFT, justified
      !! there by the energy AT pyscf's density agreeing to 1e-8: that
      !! at-density run (one iteration from DM_10, no GDM) is repeated here.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: base, opts2
      type(trc_scf_result_t) :: res, plain, at
      real(dp) :: e_ref
      real(dp), allocatable :: dm0(:, :, :), dm_ref(:, :, :)
      integer :: na, nb
      call ref_energy("rks.ref", "E_10", e_ref, error)
      if (allocated(error)) return
      call make_uncontracted('oh.xyz', bas, na, nb, error)
      if (allocated(error)) return
      allocate (dm0(bas%nao, bas%nao, 2), dm_ref(bas%nao, bas%nao, 2))
      call ref_density("rks.ref", "DM0_10", bas%nao, 2, dm0, error)
      if (.not. allocated(error)) call ref_density("rks.ref", "DM_10", bas%nao, 2, dm_ref, error)
      if (allocated(error)) then
         call bas%release()
         return
      end if
      base = base_options()
      base%functional = "pbe"
      base%grid_level = 3
      base%conv_diis = 1.0e-8_dp
      opts2 = base
      opts2%max_iteration = 1
      call trc_scf_run(bas, na, nb, opts2, at, dguess=dm_ref)
      call trc_scf_run(bas, na, nb, gdm_options(base), res, dguess=dm0)
      call trc_scf_run(bas, na, nb, base, plain, dguess=dm0)
      call bas%release()
      call check(error, abs(at%energy - e_ref) < TOL_REF, "terco's energy at pyscf's density differs from pyscf by "// &
                 to_string(abs(at%energy - e_ref)))
      if (allocated(error)) return
      call check_gdm_ran(error, res, "OH PBE")
      if (allocated(error)) return
      call check(error, res%nspin == 2, "expected an unrestricted run")
      if (allocated(error)) return
      call check(error, abs(res%energy - e_ref) < 2.0e-6_dp, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - e_ref)))
      if (allocated(error)) return
      call check(error, abs(res%energy - plain%energy) < 2.0e-6_dp, "GDM and DIIS energies differ by "// &
                 to_string(abs(res%energy - plain%energy)))
   end subroutine test_oh_uks

   ! --- continuum, density fitting ----------------------------------------

   subroutine pcm_case(error, xyz, tag, seeded, what)
      !! One pcm_scf_main.ref case (6-31G*, C-PCM, Hartree-Fock), against
      !! pyscf at scf_pcm.F90's 1e-8, symmetric. `seeded`: start from pyscf's
      !! converged density (OH, as scf_pcm does: its pi pair is split only by
      !! the cavity, the energy along that mode is flat to 3e-8 and has a
      !! stationary point at each end, so a run from the guess may end in
      !! either), with scf_pcm's at-density check of terco's energy at that
      !! density.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, tag, what
      logical, intent(in) :: seeded
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: base, opts2
      type(trc_scf_result_t) :: res, at
      real(dp) :: e_ref, e_dm
      real(dp), allocatable :: dm_ref(:, :, :)
      integer :: na, nb, nspin
      call ref_energy("pcm_scf_main.ref", "E_"//tag, e_ref, error)
      if (allocated(error)) return
      call make_basis(xyz, '6-31gs', na, nb, bas, error)
      if (allocated(error)) return
      base = base_options()
      base%grid_level = 3
      base%conv_diis = 1.0e-8_dp
      base%pcm_model = "cpcm"
      base%pcm_eps = EPS_WATER
      if (seeded) then
         nspin = merge(2, 1, na /= nb)
         allocate (dm_ref(bas%nao, bas%nao, nspin))
         call ref_energy("pcm_scf_main.ref", "EDM_"//tag, e_dm, error)
         if (.not. allocated(error)) call ref_density("pcm_scf_main.ref", "DM_"//tag, bas%nao, nspin, dm_ref, error)
         if (allocated(error)) then
            call bas%release()
            return
         end if
         opts2 = base
         opts2%max_iteration = 1
         call trc_scf_run(bas, na, nb, opts2, at, dguess=dm_ref)
         ! pyscf's density is converged, so the damped start (4 iterations for an
         ! open shell) would finish the run before GDM could take over:
         ! no damping, GDM from the first diagonalisation.
         base%ndamp = 0
         call trc_scf_run(bas, na, nb, gdm_options(base, early=.true.), res, dguess=dm_ref)
      else
         call trc_scf_run(bas, na, nb, gdm_options(base), res)
      end if
      call bas%release()
      if (seeded) then
         call check(error, abs(at%energy - e_dm) < TOL_REF, what//": terco's energy at pyscf's density differs by "// &
                    to_string(abs(at%energy - e_dm)))
         if (allocated(error)) return
      end if
      call check_gdm_ran(error, res, what)
      if (allocated(error)) return
      call check(error, res%e_pcm /= 0.0_dp, what//": no dielectric energy, the continuum was off")
      if (allocated(error)) return
      call check(error, abs(res%energy - e_ref) < TOL_REF, what//": energy differs from pyscf by "// &
                 to_string(res%energy - e_ref))
   end subroutine pcm_case

   subroutine test_water_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      call pcm_case(error, 'water.xyz', '0', .false., "water C-PCM")
   end subroutine test_water_pcm

   subroutine test_oh_pcm(error)
      type(error_type), allocatable, intent(out) :: error
      call pcm_case(error, 'oh.xyz', '12', .true., "OH C-PCM")
   end subroutine test_oh_pcm

   subroutine test_oh_pcm_flat_mode(error)
      !! OH C-PCM UHF from the default guess, every energy a full build's
      !! (`incremental` off), conv_diis 1e-8: the case whose pi pair makes the
      !! energy flat to 3e-8 along one mode (scf_pcm.F90's header). A rule that
      !! accepted real rises there let the energy random-walk along it (291
      !! iterations, 49 rejections, against DIIS's 118); with the plain slack
      !! the walk cannot happen. Asserts the exact bound on every rise, an
      !! iteration ceiling, and that the converged energy is within the flat
      !! mode's 3e-8 of pyscf's stationary point (either end of the mode is a
      !! stationary point, so this is two-sided by construction).
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: base
      type(trc_scf_result_t) :: res
      real(dp) :: e_ref
      integer :: na, nb
      call ref_energy("pcm_scf_main.ref", "E_12", e_ref, error)
      if (allocated(error)) return
      call make_basis('oh.xyz', '6-31gs', na, nb, bas, error)
      if (allocated(error)) return
      base = exact_options()
      base%grid_level = 3
      base%conv_diis = 1.0e-8_dp
      base%pcm_model = "cpcm"
      base%pcm_eps = EPS_WATER
      base%max_iteration = FLAT_CEILING + 100
      call trc_scf_run(bas, na, nb, gdm_options(base), res)
      call bas%release()
      call check_gdm_ran(error, res, "OH C-PCM flat mode", exact=.true.)
      if (allocated(error)) return
      call check(error, res%iterations <= FLAT_CEILING, "the flat mode took "//to_string(res%iterations)// &
                 " iterations, ceiling "//to_string(FLAT_CEILING))
      if (allocated(error)) return
      call check(error, abs(res%energy - e_ref) < 3.0e-8_dp, "energy is "//to_string(res%energy - e_ref)// &
                 " from pyscf's; the flat mode allows 3e-8")
   end subroutine test_oh_pcm_flat_mode

   subroutine test_water_df(error)
      !! scf_df.ref case 0: water, cc-pVDZ, Hartree-Fock, def2-universal-jkfit.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t) :: dfjk
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp) :: e_ref
      integer :: natm, na, nb
      call ref_energy("scf_df.ref", "E_0", e_ref, error)
      if (allocated(error)) return
      call read_xyz("water.xyz", natm, z, at_r)
      na = sum(z)/2
      nb = na
      call trc_basis_from_json("basis_sets/cc-pvdz.json", natm, z, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call trc_basis_from_json("basis_sets/def2-universal-jkfit.json", natm, z, at_r, aux, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "auxiliary basis: "//err%get_message())
      if (allocated(error)) return
      if (maxval(aux%sh_l) > DF_LMAX_AUX) then
         call skip_test(error, "def2-universal-jkfit needs l = "//to_string(maxval(aux%sh_l))// &
                        " and this build has DF_LMAX_AUX = "//to_string(DF_LMAX_AUX))
         return
      end if
      call bas%to_device()
      call aux%to_device()
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      call la%init(bas%nao)
      call dfjk%setup(bas, aux, pl, la, err, nocc_hint=max(na, 1))
      call check(error, .not. err%has_error(), "DF setup: "//err%get_message())
      if (.not. allocated(error)) then
         opts = gdm_options(base_options())
         opts%df = .true.
         opts%conv_diis = 1.0e-8_dp
         call trc_scf_run(bas, na, nb, opts, res, la=la, scf_dfjk=dfjk)
      end if
      call pl%release()
      call dfjk%release()
      call la%release()
      call aux%release()
      call bas%release()
      if (allocated(error)) return
      call check_gdm_ran(error, res, "water DF-RHF")
      if (allocated(error)) return
      call check(error, abs(res%energy - e_ref) < TOL_REF, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - e_ref)))
   end subroutine test_water_df

   ! --- linear dependence -------------------------------------------------

   subroutine test_lindep(error)
      !! check_scf_lindep_gap's O2- (d-aug-cc-pVDZ, threshold 0.2, 5 functions
      !! dropped). The retained levels are the DIIS run's; the dropped states
      !! stay above every real level, where `scf_driver_null_shift` put them.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: base
      type(trc_scf_result_t) :: res, plain
      integer :: z(1), n_mo, nao
      real(dp) :: r(3, 1)
      z = 8
      r = 0.0_dp
      call trc_basis_from_json("basis_sets/d-aug-cc-pvdz.json", 1, z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      base = base_options()
      base%lindep_thresh = 0.2_dp
      base%conv_diis = 1.0e-9_dp
      call trc_scf_run(bas, 5, 5, gdm_options(base), res)
      call trc_scf_run(bas, 5, 5, base, plain)
      call bas%release()
      call check(error, res%n_lindep == 5, "expected 5 functions dropped, got "//to_string(res%n_lindep))
      if (allocated(error)) return
      call check_gdm_ran(error, res, "O2- lindep")
      if (allocated(error)) return
      call check(error, abs(res%energy - E_O2M_LINDEP) < TOL_REF, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - E_O2M_LINDEP)))
      if (allocated(error)) return
      call check(error, plain%converged .and. abs(res%energy - plain%energy) < 1.0e-8_dp, &
                 "GDM and DIIS energies differ")
      if (allocated(error)) return
      nao = size(res%eps, 1)
      n_mo = nao - res%n_lindep
      call check(error, maxval(abs(res%eps(1:n_mo, 1) - plain%eps(1:n_mo, 1))) < 1.0e-6_dp, &
                 "retained orbital energies differ from the DIIS run's by "// &
                 to_string(maxval(abs(res%eps(1:n_mo, 1) - plain%eps(1:n_mo, 1)))))
      if (allocated(error)) return
      call check(error, minval(res%eps(n_mo + 1:nao, 1)) > maxval(res%eps(1:n_mo, 1)), &
                 "a dropped state fell below a real level")
   end subroutine test_lindep

   ! --- handover, rejection, orbitals -------------------------------------

   subroutine test_default_handover(error)
      !! `gdm_start` left at its default: DIIS (the first iterations, then the
      !! handover iteration's plain diagonalisation), GDM after.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res, plain
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = exact_options()
      opts%gdm = .true.
      call trc_scf_run(bas, na, nb, opts, res)
      opts%gdm = .false.
      call trc_scf_run(bas, na, nb, opts, plain)
      call bas%release()
      call check_gdm_ran(error, res, "default handover", exact=.true.)
      if (allocated(error)) return
      call check(error, res%iterations - res%gdm_steps >= 2, "no DIIS iteration before the handover: "// &
                 to_string(res%iterations)//" iterations, "//to_string(res%gdm_steps)//" GDM steps")
      if (allocated(error)) return
      call check(error, abs(res%energy - plain%energy) < 1.0e-9_dp .and. abs(res%energy - E_WATER_RHF) < TOL_REF, &
                 "energy differs from the DIIS run or from pyscf")
   end subroutine test_default_handover

   subroutine test_accept_slack(error)
      !! The accept rule's slack, all four kinds of comparison. A rule that
      !! widened it for a full-build comparison (or for none, or only for the
      !! trial) let the energy random-walk along OH's flat mode, accepting rises
      !! of 8e-9 with no noise source at all.
      type(error_type), allocatable, intent(out) :: error
      call check(error, scf_gdm_step_accept_slack(SLACK, .false., .false.) == SLACK, "full vs full is not the plain slack")
      if (allocated(error)) return
      call check(error, scf_gdm_step_accept_slack(SLACK, .true., .false.) == NOISE_BOUND, "an incremental trial is not widened")
      if (allocated(error)) return
      call check(error, scf_gdm_step_accept_slack(SLACK, .false., .true.) == NOISE_BOUND, &
                 "an incremental accepted point is not widened")
      if (allocated(error)) return
      call check(error, scf_gdm_step_accept_slack(SLACK, .true., .true.) == NOISE_BOUND, "incremental vs incremental is not widened")
      if (allocated(error)) return
      call check(error, scf_gdm_step_accept_slack(NOISE_BOUND*10.0_dp, .true., .true.) == NOISE_BOUND*10.0_dp, &
                 "a larger option was narrowed")
      if (allocated(error)) return
      call check(error, scf_gdm_step_accept_slack(-1.0_dp, .true., .true.) == -1.0_dp, "a negative slack (a demanded decrease) was widened")
   end subroutine test_accept_slack

   subroutine test_adiis_handover(error)
      !! ADIIS runs while the commutator is above `diis_start` (set to 1e-3 so
      !! it runs through the whole window above the handover at 0.1), then GDM
      !! takes over: the energy is the ADIIS/DIIS run's of the same input and
      !! pyscf's. Without the handover switching ADIIS off, its extrapolated
      !! Fock matrix would be the one GDM takes its gradient from.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res, plain
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = exact_options()
      opts%adiis = .true.
      opts%diis_start = 1.0e-3_dp
      opts%gdm = .true.
      call trc_scf_run(bas, na, nb, opts, res)
      opts%gdm = .false.
      call trc_scf_run(bas, na, nb, opts, plain)
      call bas%release()
      call check_gdm_ran(error, res, "adiis + gdm", exact=.true.)
      if (allocated(error)) return
      call check(error, res%iterations - res%gdm_steps >= 2, "no ADIIS iteration before the handover: "// &
                 to_string(res%iterations)//" iterations, "//to_string(res%gdm_steps)//" GDM steps")
      if (allocated(error)) return
      call check(error, plain%converged .and. abs(res%energy - plain%energy) < 1.0e-9_dp .and. &
                 abs(res%energy - E_WATER_RHF) < TOL_REF, "energy differs from the ADIIS run or from pyscf")
   end subroutine test_adiis_handover

   subroutine test_rejected_steps(error)
      !! GDM from the guess with a maximum rotation of 4 radians, far more than
      !! the energy surface is quadratic over: the first quasi-Newton steps
      !! overshoot and must be rejected and halved, and the run still lands on
      !! pyscf's energy with every accepted point below the one before.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = gdm_options(exact_options(), early=.true.)
      opts%gdm_max_rotation = 4.0_dp
      opts%guess = "core"
      call trc_scf_run(bas, na, nb, opts, res)
      call bas%release()
      call check_gdm_ran(error, res, "large rotation", exact=.true.)
      if (allocated(error)) return
      call check(error, res%gdm_rejected > 0, "no trial was rejected, so the rejection path was not exercised")
      if (allocated(error)) return
      call check(error, abs(res%energy - E_WATER_RHF) < TOL_REF, "energy differs from pyscf by "// &
                 to_string(abs(res%energy - E_WATER_RHF)))
   end subroutine test_rejected_steps

   subroutine test_rotation_cap(error)
      !! From the same poor start, the first GDM step lowers water's energy by
      !! 0.55 Eh at the default cap (0.2 rad). With a cap of 1e-3 rad it can
      !! lower it by at most the gradient's 1-norm times that, ~1e-2 Eh: the
      !! energy after the second step (iteration 3) sits below the first
      !! accepted point's, by less than 0.05. And every rotation element of
      !! EVERY trial of a longer run stays within the cap (`gdm_max_step`, the
      !! largest applied), while some step does reach it: a cap that held only
      !! the first step (or none, with nothing near it) fails one of those.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: one, two, many
      integer :: na, nb
      real(dp), parameter :: CAP = 1.0e-3_dp
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = gdm_options(base_options(), early=.true.)
      opts%gdm_max_rotation = CAP
      opts%max_iteration = 2
      call trc_scf_run(bas, na, nb, opts, one)
      opts%max_iteration = 3
      call trc_scf_run(bas, na, nb, opts, two)
      opts%max_iteration = 12
      call trc_scf_run(bas, na, nb, opts, many)
      call bas%release()
      call check(error, two%energy < one%energy, "the capped step did not lower the energy")
      if (allocated(error)) return
      call check(error, one%energy - two%energy < 0.05_dp, "a capped step lowered the energy by "// &
                 to_string(one%energy - two%energy)//" Eh: the cap was not applied")
      if (allocated(error)) return
      call check(error, many%gdm_steps >= 8, "too few GDM steps to test every step: "//to_string(many%gdm_steps))
      if (allocated(error)) return
      call check(error, many%gdm_max_step <= CAP*(1.0_dp + 1.0e-12_dp), "a step rotated by "// &
                 to_string(many%gdm_max_step)//" rad, above the cap "//to_string(CAP))
      if (allocated(error)) return
      call check(error, many%gdm_max_step >= 0.5_dp*CAP, "no step came near the cap ("//to_string(many%gdm_max_step)// &
                 " rad): the cap was never exercised")
   end subroutine test_rotation_cap

   subroutine test_orbitals_consistent(error)
      !! The orbitals and density a GDM run hands back belong together and are
      !! the DIIS run's: the density rebuilt from `cmo` and the occupations is
      !! `dmat`, and every orbital energy equals the canonical one (at the
      !! minimum the semicanonical orbitals ARE the canonical ones, up to the
      !! commutator). Without the final semicanonicalisation `eps` would be
      !! the handover iteration's.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: base
      type(trc_scf_result_t) :: res, plain
      real(dp), allocatable :: d(:, :)
      real(dp) :: occ
      integer :: na, nb, nao, s, nocc_s
      character(len=*), parameter :: SYSTEM(2) = [character(len=9) :: 'water.xyz', 'oh.xyz   ']
      integer :: isys
      do isys = 1, 2
         call make_basis(trim(SYSTEM(isys)), 'cc-pvdz', na, nb, bas, error)
         if (allocated(error)) return
         base = base_options()
         base%conv_diis = 1.0e-9_dp
         call trc_scf_run(bas, na, nb, gdm_options(base), res)
         call trc_scf_run(bas, na, nb, base, plain)
         ! The SCF's own dimension (spherical here), not bas%nao, which counts
         ! Cartesian functions. With `d` allocated at the wrong shape the
         ! matmul assignment below reallocates it, and nvfortran 24.11-26.9
         ! corrupts the heap doing that.
         nao = size(res%dmat, 1)
         call bas%release()
         call check_gdm_ran(error, res, trim(SYSTEM(isys)))
         if (allocated(error)) return
         allocate (d(nao, nao))
         do s = 1, res%nspin
            nocc_s = merge(na, nb, s == 1)
            occ = merge(2.0_dp, 1.0_dp, res%nspin == 1)
            d = occ*matmul(res%cmo(:, 1:nocc_s, s), transpose(res%cmo(:, 1:nocc_s, s)))
            call check(error, maxval(abs(d - res%dmat(:, :, s))) < 1.0e-8_dp, trim(SYSTEM(isys))// &
                       ": the density is not the one the returned orbitals make, spin "//to_string(s))
            if (allocated(error)) return
            call check(error, maxval(abs(res%eps(1:nao, s) - plain%eps(1:nao, s))) < 1.0e-6_dp, &
                       trim(SYSTEM(isys))//": orbital energies differ from the DIIS run's, spin "//to_string(s)// &
                       " by "//to_string(maxval(abs(res%eps(1:nao, s) - plain%eps(1:nao, s)))))
            if (allocated(error)) return
         end do
         deallocate (d)
      end do
   end subroutine test_orbitals_consistent

   subroutine h_atom(basis, bas, error)
      type(trc_basis_t), intent(out) :: bas
      character(len=*), intent(in) :: basis
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: z(1)
      real(dp) :: r(3, 1)
      z = 1
      r = 0.0_dp
      call trc_basis_from_json('basis_sets/'//basis//'.json', 1, z, r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
   end subroutine h_atom

   subroutine compare_with_diis(error, bas, what)
      !! One electron (alpha only), GDM from the first iteration it may take
      !! over (`ndamp` 0), against the DIIS run of the same input: energy to
      !! 1e-8 and every orbital energy of both spins to 1e-6. No external
      !! reference is involved; the one-electron answer needs none, the
      !! point is the spin block GDM has no occupied (or no virtual) orbital in.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t), intent(in) :: bas
      character(len=*), intent(in) :: what
      type(trc_scf_options_t) :: base
      type(trc_scf_result_t) :: res, plain
      base = base_options()
      base%ndamp = 0
      base%conv_diis = 1.0e-9_dp
      call trc_scf_run(bas, 1, 0, gdm_options(base, early=.true.), res)
      call trc_scf_run(bas, 1, 0, base, plain)
      call check_gdm_ran(error, res, what)
      if (allocated(error)) return
      call check(error, plain%converged .and. abs(res%energy - plain%energy) < 1.0e-8_dp, &
                 what//": GDM and DIIS energies differ by "//to_string(abs(res%energy - plain%energy)))
      if (allocated(error)) return
      call check(error, maxval(abs(res%eps - plain%eps)) < 1.0e-6_dp, what//": orbital energies differ from the DIIS run's by "// &
                 to_string(maxval(abs(res%eps - plain%eps))))
   end subroutine compare_with_diis

   subroutine test_empty_spin(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      call h_atom('cc-pvdz', bas, error)
      if (allocated(error)) return
      call compare_with_diis(error, bas, "H atom cc-pVDZ")
      call bas%release()
   end subroutine test_empty_spin

   subroutine test_nothing_to_rotate(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      call h_atom('sto-3g', bas, error)
      if (allocated(error)) return
      call compare_with_diis(error, bas, "H atom STO-3G")
      call bas%release()
   end subroutine test_nothing_to_rotate

   ! --- refusals, breakdown -----------------------------------------------

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = base_options()
      opts%gdm = .true.
      opts%frac_occ = .true.
      opts%nelec_frac = 10.0_dp
      call trc_scf_run(bas, na, nb, opts, res)
      call check(error, .not. res%converged .and. res%iterations == 0 .and. index(res%message, "gdm") > 0 &
                 .and. index(res%message, "fractional") > 0, "gdm + frac_occ was not refused: "//trim(res%message))
      if (allocated(error)) then
         call bas%release()
         return
      end if
      opts = base_options()
      opts%gdm = .true.
      opts%soscf = .true.
      call trc_scf_run(bas, na, nb, opts, res)
      call bas%release()
      call check(error, .not. res%converged .and. res%iterations == 0 .and. index(res%message, "gdm") > 0 &
                 .and. index(res%message, "soscf") > 0, "gdm + soscf was not refused: "//trim(res%message))
   end subroutine test_refusals

   subroutine test_rejected_never_converges(error)
      !! `gdm_slack` = -1 rejects every trial but the first, and the thresholds
      !! are set from this input so that the first rejected trial (iteration 3)
      !! meets them and the accepted one before it (iteration 2) does not:
      !! conv_diis is huge (always met), conv_energy is 1.1 times the energy the
      !! first trial differs from the handover point by, which is far below
      !! what the handover point differs from the guess by (asserted, so the
      !! setup cannot silently stop discriminating). If a rejected trial
      !! could converge the run it would, at iteration 3. It must not, so the
      !! run goes on to the GDM stop.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts, probe
      type(trc_scf_result_t) :: res, at_guess, at_handover, at_trial
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      probe = gdm_options(exact_options(), early=.true.)
      probe%max_iteration = 1
      call trc_scf_run(bas, na, nb, probe, at_guess)
      probe%max_iteration = 2
      call trc_scf_run(bas, na, nb, probe, at_handover)
      probe%max_iteration = 3
      probe%gdm_slack = 1.0_dp
      call trc_scf_run(bas, na, nb, probe, at_trial)
      opts = gdm_options(exact_options(), early=.true.)
      opts%gdm_slack = -1.0_dp
      opts%conv_energy = 1.1_dp*abs(at_trial%energy - at_handover%energy)
      opts%conv_diis = huge(1.0_dp)
      call check(error, abs(at_handover%energy - at_guess%energy) > opts%conv_energy .and. opts%conv_energy > 0.0_dp, &
                 "the thresholds do not separate iteration 2 from iteration 3: guess-to-handover "// &
                 to_string(abs(at_handover%energy - at_guess%energy))//", conv_energy "//to_string(opts%conv_energy))
      if (.not. allocated(error)) call trc_scf_run(bas, na, nb, opts, res)
      call bas%release()
      if (allocated(error)) return
      call check(error, .not. res%converged, "a rejected trial converged the run at iteration "// &
                 to_string(res%iterations))
      if (allocated(error)) return
      call check(error, res%iterations > 3 .and. res%gdm_rejected > 0 .and. index(res%message, GDM_MESSAGE_PART) > 0, &
                 "expected the GDM stop after the rejections, got: "//trim(res%message))
   end subroutine test_rejected_never_converges

   subroutine test_no_downhill(error)
      !! A negative `gdm_slack` demands a decrease of a hartree: no step gives
      !! it, so every trial is rejected. The first direction is the
      !! steepest descent (empty history), so after the one round of eleven
      !! halvings there is nothing else to try and the run stops, saying so
      !! -- not report the last (rejected) trial as converged, not repeat the
      !! round with the same direction -- and leaves a finite energy.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res, first
      integer :: na, nb, one_round
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = gdm_options(exact_options(), early=.true.)
      opts%gdm_slack = -1.0_dp
      call trc_scf_run(bas, na, nb, opts, res)
      ! The one accepted point of that run is its first GDM iteration's: the
      ! same input stopped after that iteration (default slack) reports it.
      opts%gdm_slack = 0.0_dp
      opts%max_iteration = 2
      call trc_scf_run(bas, na, nb, opts, first)
      call bas%release()
      print "(a,i0,a,i0,a,f20.12,a,a)", "    no downhill: iterations ", res%iterations, "  rejected ", res%gdm_rejected, &
         "  E = ", res%energy, "  ", trim(res%message)
      call check(error, .not. res%converged, "the run claims convergence with no downhill step")
      if (allocated(error)) return
      call check(error, index(res%message, GDM_MESSAGE_PART) > 0, "the message does not say GDM found no downhill "// &
                 "step: "//trim(res%message))
      if (allocated(error)) return
      ! A scale of 1 halved until it is below the minimum: 11 rejections.
      one_round = nint(log(1.0_dp/scf_gdm_scale_min)/log(2.0_dp)) + 1
      call check(error, res%gdm_rejected == one_round, "expected one round of "//to_string(one_round)// &
                 " rejections, got "//to_string(res%gdm_rejected))
      if (allocated(error)) return
      call check(error, abs(res%energy - first%energy) < 1.0e-10_dp, "the energy reported is not the last accepted "// &
                 "point's (the first GDM iteration's): it differs by "//to_string(abs(res%energy - first%energy)))
      if (allocated(error)) return
      call check(error, maxval(abs(res%dmat - first%dmat)) < 1.0e-10_dp, "the density reported is a rejected trial's, "// &
                 "not the last accepted point's")
      if (allocated(error)) return
      call check(error, maxval(abs(res%eps - first%eps)) < 1.0e-8_dp, "the orbital energies are not those of the "// &
                 "last accepted point")
   end subroutine test_no_downhill

   subroutine test_ceiling_on_rejected(error)
      !! `max_iteration` reached while the trial is a rejected one (water
      !! 6-31G* C-PCM, every trial after the first rejected by a negative
      !! slack, five iterations: three rejections, no breakdown yet): the
      !! energy, the dielectric energy, the apparent surface charge, the
      !! density and the orbital energies are those of the last ACCEPTED
      !! point, which a run stopped there (two iterations) reports as its own.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res, accepted
      integer :: na, nb
      call make_basis('water.xyz', '6-31gs', na, nb, bas, error)
      if (allocated(error)) return
      opts = gdm_options(exact_options(), early=.true.)
      opts%pcm_model = "cpcm"
      opts%pcm_eps = EPS_WATER
      opts%gdm_slack = -1.0_dp
      opts%max_iteration = 5
      call trc_scf_run(bas, na, nb, opts, res)
      opts%gdm_slack = 0.0_dp
      opts%max_iteration = 2
      call trc_scf_run(bas, na, nb, opts, accepted)
      call bas%release()
      call check(error, .not. res%converged .and. res%iterations == 5 .and. res%gdm_rejected == 3, &
                 "expected 5 iterations, 3 rejected, not converged: "//to_string(res%iterations)//", "// &
                 to_string(res%gdm_rejected)//", "//trim(res%message))
      if (allocated(error)) return
      call check(error, abs(res%energy - accepted%energy) < 1.0e-10_dp .and. abs(res%e_pcm - accepted%e_pcm) < 1.0e-10_dp, &
                 "the energies are not the accepted point's: "//to_string(res%energy - accepted%energy))
      if (allocated(error)) return
      call check(error, abs(res%pcm_qtotal - accepted%pcm_qtotal) < 1.0e-10_dp .and. accepted%pcm_qtotal /= 0.0_dp, &
                 "the surface charge is a rejected trial's: "//to_string(res%pcm_qtotal - accepted%pcm_qtotal))
      if (allocated(error)) return
      call check(error, maxval(abs(res%dmat - accepted%dmat)) < 1.0e-10_dp, "the density is a rejected trial's")
      if (allocated(error)) return
      call check(error, maxval(abs(res%eps - accepted%eps)) < 1.0e-8_dp, "the orbital energies are not the accepted point's")
   end subroutine test_ceiling_on_rejected

   subroutine test_invalid_rotation(error)
      !! An infinite maximum rotation is not a usable cap: the first direction
      !! is refused (`scf_gdm_step_limit` reports it invalid), which stops the
      !! run with the GDM message, not converged, on the accepted point (that of
      !! the same input stopped after two iterations). The stop on a
      !! non-finite GRADIENT shares this exit but cannot be provoked from the
      !! options, and is not covered.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res, accepted
      integer :: na, nb
      call make_basis('water.xyz', 'cc-pvdz', na, nb, bas, error)
      if (allocated(error)) return
      opts = gdm_options(exact_options(), early=.true.)
      opts%gdm_max_rotation = ieee_value(1.0_dp, ieee_positive_inf)
      call trc_scf_run(bas, na, nb, opts, res)
      opts%gdm_max_rotation = 0.0_dp
      opts%max_iteration = 2
      call trc_scf_run(bas, na, nb, opts, accepted)
      call bas%release()
      call check(error, .not. res%converged .and. index(res%message, GDM_MESSAGE_PART) > 0 .and. res%iterations == 2, &
                 "expected the GDM stop at iteration 2, got "//to_string(res%iterations)//": "//trim(res%message))
      if (allocated(error)) return
      call check(error, abs(res%energy - accepted%energy) < 1.0e-10_dp .and. &
                 maxval(abs(res%dmat - accepted%dmat)) < 1.0e-10_dp, "the result is not the accepted point")
   end subroutine test_invalid_rotation

   subroutine test_mom_gdm(error)
      !! check_scf_stabilise's test_mom_holds_excited_rhf set up (water
      !! cc-pVDZ, the ground state's HOMO and LUMO swapped in the given
      !! density): `mom` alone holds the excited state; `mom` + `gdm` must keep
      !! holding it, at the default `gdm_start` and with the handover at the
      !! first iteration. A minimiser cannot do that by itself (the held state
      !! is a saddle point and it runs down to the ground state: measured, the
      !! control below), so GDM must decline at the handover (`gdm_declined`,
      !! no GDM step) when MOM's occupation is not aufbau. When MOM merely
      !! latched on the aufbau occupation (given the ground density) GDM must
      !! not decline.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o
      type(trc_scf_result_t) :: ground, held, g_default, g_early, g_nomom, g_aufbau
      real(dp), allocatable :: dex(:, :, :)
      integer, allocatable :: occ_cols(:)
      integer :: nocc, k, nb
      call make_basis('water.xyz', 'cc-pvdz', nocc, nb, bas, error)
      if (allocated(error)) return
      o = exact_options()
      o%conv_energy = 1.0e-10_dp
      o%conv_diis = 1.0e-8_dp
      o%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, o, ground)
      call check(error, ground%converged, "ground state: "//trim(ground%message))
      if (allocated(error)) then
         call bas%release()
         return
      end if
      occ_cols = [(k, k=1, nocc - 1), nocc + 1]
      allocate (dex(size(ground%dmat, 1), size(ground%dmat, 2), 1))
      dex(:, :, 1) = 2.0_dp*matmul(ground%cmo(:, occ_cols, 1), transpose(ground%cmo(:, occ_cols, 1)))
      o%mom = .true.
      call trc_scf_run(bas, nocc, nocc, o, held, dguess=dex)
      call trc_scf_run(bas, nocc, nocc, gdm_options(o), g_default, dguess=dex)
      call trc_scf_run(bas, nocc, nocc, gdm_options(o, early=.true.), g_early, dguess=dex)
      call trc_scf_run(bas, nocc, nocc, gdm_options(o, early=.true.), g_aufbau, dguess=ground%dmat)
      o%mom = .false.
      call trc_scf_run(bas, nocc, nocc, gdm_options(o, early=.true.), g_nomom, dguess=dex)
      call bas%release()
      call check(error, held%converged .and. held%energy > ground%energy + 0.3_dp, "mom alone did not hold the excited state")
      if (allocated(error)) return
      call check(error, g_default%converged .and. g_early%converged, "mom + gdm did not converge: "// &
                 trim(g_default%message)//" / "//trim(g_early%message))
      if (allocated(error)) return
      call check(error, g_default%gdm_declined .and. g_early%gdm_declined .and. g_default%gdm_steps == 0 .and. &
                 g_early%gdm_steps == 0, "GDM did not decline against a mom-held non-aufbau occupation")
      if (allocated(error)) return
      call check(error, abs(g_default%energy - held%energy) < 1.0e-7_dp .and. &
                 abs(g_early%energy - held%energy) < 1.0e-7_dp, "mom + gdm left the mom state: E - E(mom) = "// &
                 to_string(g_default%energy - held%energy)//", "//to_string(g_early%energy - held%energy))
      if (allocated(error)) return
      call check(error, g_default%energy > ground%energy + 0.3_dp .and. g_early%energy > ground%energy + 0.3_dp, &
                 "mom + gdm fell to the ground state")
      if (allocated(error)) return
      call check(error, g_nomom%converged .and. g_nomom%gdm_steps > 0 .and. &
                 abs(g_nomom%energy - ground%energy) < 1.0e-7_dp, &
                 "without mom GDM should run down to the ground state (the control of this test)")
      if (allocated(error)) return
      call check_gdm_ran(error, g_aufbau, "mom + gdm on the aufbau occupation", exact=.true.)
      if (allocated(error)) return
      call check(error, .not. g_aufbau%gdm_declined .and. abs(g_aufbau%energy - ground%energy) < 1.0e-8_dp, &
                 "GDM declined, or missed the ground state, with mom on the aufbau occupation")
   end subroutine test_mom_gdm

   subroutine test_started_too_early(error)
      !! A minimiser keeps the occupation of the diagonalisation it starts
      !! from. Water B3LYP/6-31G* from a core guess with `gdm_start` huge
      !! (GDM from the first diagonalisation of the guess Fock matrix) at the
      !! default thresholds (conv_diis 1e-5) stops at a stationary point 0.82 Eh
      !! above the ground state, reported converged: `gdm_off_aufbau` > 0 (it
      !! logs a warning) is the only sign. (At conv_diis 1e-7 the same run
      !! slides off the saddle point and finds the ground state: the test
      !! needs the loose threshold, which is the default.) The same input at the
      !! default `gdm_start` (DIIS settles the occupation first) reaches
      !! pyscf's energy (cpks_water_6-31gs_b3lyp.ref) and is aufbau, at both
      !! thresholds, as is the DIIS run of the same input.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: base, loose
      type(trc_scf_result_t) :: early, normal, normal_loose, plain
      real(dp) :: e_ref
      integer :: na, nb
      call ref_energy("cpks_water_6-31gs_b3lyp.ref", "E_SCF", e_ref, error)
      if (allocated(error)) return
      call make_basis('water.xyz', '6-31gs', na, nb, bas, error)
      if (allocated(error)) return
      base = exact_options()
      base%functional = "b3lyp"
      base%grid_level = 3
      base%conv_diis = 1.0e-8_dp
      base%guess = "core"
      loose = base
      loose%conv_diis = 1.0e-5_dp
      loose%conv_energy = 1.0e-8_dp
      base%gdm = .false.
      base%guess = "gwh"
      call trc_scf_run(bas, na, nb, base, plain)
      base%guess = "core"
      call trc_scf_run(bas, na, nb, gdm_options(loose, early=.true.), early)
      call trc_scf_run(bas, na, nb, gdm_options(loose), normal_loose)
      call trc_scf_run(bas, na, nb, gdm_options(base), normal)
      call bas%release()
      print "(a,f18.9,a,f18.9,a,f18.9,a,i0,a,i0)", "    started early E = ", early%energy, "  default E = ", normal%energy, &
         "  DIIS E = ", plain%energy, "  off-aufbau early ", early%gdm_off_aufbau, " default ", normal%gdm_off_aufbau
      call check(error, plain%converged .and. abs(plain%energy - e_ref) < TOL_REF, "the DIIS run does not reach pyscf's energy")
      if (allocated(error)) return
      call check(error, early%converged .and. early%gdm_steps > 0 .and. early%energy > e_ref + 0.3_dp, &
                 "the started-too-early run did not end converged on a higher state: E - E(pyscf) = "// &
                 to_string(early%energy - e_ref))
      if (allocated(error)) return
      call check(error, early%gdm_off_aufbau > 0, "the higher state is not flagged non-aufbau (gdm_off_aufbau = "// &
                 to_string(early%gdm_off_aufbau)//")")
      if (allocated(error)) return
      call check_gdm_ran(error, normal_loose, "default start, core guess, loose", exact=.true.)
      if (allocated(error)) return
      call check_gdm_ran(error, normal, "default start, core guess", exact=.true.)
      if (allocated(error)) return
      call check(error, abs(normal%energy - e_ref) < TOL_REF .and. abs(normal_loose%energy - e_ref) < TOL_REF, &
                 "the default-start run differs from pyscf by "//to_string(normal%energy - e_ref)//", loose "// &
                 to_string(normal_loose%energy - e_ref))
   end subroutine test_started_too_early

end module check_scf_gdm_run_suite

program check_scf_gdm_run
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_scf_gdm_run_suite, only: collect_scf_gdm_run
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_gdm_run, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_gdm_run
