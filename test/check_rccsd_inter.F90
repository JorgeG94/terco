!
! The one-particle intermediates and the energy functional of RI-CCSD(T)
! (src/post_hf/ccsd/trc_rccsd_inter.F90) against pyscf's DF-RCCSD, on a water and
! on a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT, spherical, frozen core 1
! and 2 (the molecules, SCF and probes of check_rccsd_ints, shared through
! rccsd_test_support).
!
! The amplitudes are the FIXED, unconverged ones test/pyscf_refs/generate.py
! (gen_rccsd_ints) probes: t2 = MP2, t1_ia = (0.03 x + 0.02 y + 0.01 z)_ia / (e_i - e_a)
! with x, y, z the MO dipole blocks (origin 0) over the active orbitals. The
! reference holds, per tensor, sign-invariant probes: for a 2-index A,
! |A|, sum A_pq X_pq and sum A_pq (XY)_pq, XY = X @ Y; for a 4-index one
! |A|, sum A X_pq Y_rs, sum A X_pr Z_qs, sum A Y_ps Z_qr.
!
! INDEX ORDER. mqc's arrays are PySCF's, index for index: t1(i,a), t2(i,j,a,b),
! Foo(k,i), Fvv(a,c), Fov(k,c), Loo(k,i), Lvv(a,c). One difference of content:
! pyscf's cc_Foo, cc_Fvv (and so Loo, Lvv) add the bare Fock block, the diagonal
! eps_o, eps_v for canonical orbitals, which the module leaves out (its header
! says why). The probes below therefore add diag(eps) to the module's array
! first, and the reference's `fov` terms (zero up to 1e-12) are not in the module.
!
! Cases, per molecule:
!   1 t1, t2 (mp2_guess) probes against the reference; tau against an
!     independent host t2 + t1 t1; t2 against a host quotient
!   2 Foo and Fvv probes
!   3 Fov, Loo and Lvv probes
!   4 E_MP2 and the energy functional at the fixed amplitudes against the
!     reference, singles/doubles split against a host sum, repeatability
!   5 t1 = 0: Fov = 0, Loo = Foo, Lvv = Fvv, Foo and Fvv the t2-only forms
!   6 all five at the fixed amplitudes against host loops written term by term
!     as pyscf's einsums (t2 and t1 t1 terms separate, no tau, no L)
!   7 refusals: an empty integral type, mismatched sizes of amplitudes or intermediates
!     (status set, message non-empty, nothing changed)
!   8 stale scratch: after `ints%build` is called again the routines that take `ints` refuse
!     until the scratch is re-initialised, then reproduce the earlier numbers
!
module check_rccsd_inter_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t, rccsd_inter_one_particle, &
                              rccsd_inter_energy, rccsd_inter_mp2_guess, rccsd_inter_build_tau, rccsd_inter_foo, rccsd_inter_fvv
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, TOL_PROBE, PROBE_FLOOR, load_molecule, &
                                 run_scf, mo_dipoles, probe2, probe4
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_rccsd_inter

   integer, parameter :: NCASE = 8
   character(len=4), parameter :: SKIP_TAG(5) = [character(len=4) :: "oooo", "ooov", "oovv", "ovov", "ovvv"]
   real(dp), parameter :: TOL_ENERGY = 1.0e-10_dp   !! E_MP2, E_fixed against pyscf
   real(dp), parameter :: TOL_HOST = 1.0e-12_dp     !! against independent host loops, relative to max(1, |A|max)

   logical :: computed = .false.
   logical :: ok_case(NCASE, NMOL)
   character(len=700) :: msg_case(NCASE, NMOL)

contains

   subroutine collect_check_rccsd_inter(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: t1, t2, tau", w1), &
                  new_unittest("water: Foo, Fvv vs pyscf", w2), &
                  new_unittest("water: Fov, Loo, Lvv vs pyscf", w3), &
                  new_unittest("water: E_MP2 and energy functional vs pyscf", w4), &
                  new_unittest("water: t1 = 0 reductions", w5), &
                  new_unittest("water: intermediates vs host einsum loops", w6), &
                  new_unittest("water: refusals", w7), &
                  new_unittest("water: stale scratch", w8), &
                  new_unittest("h2o2: t1, t2, tau", h1), &
                  new_unittest("h2o2: Foo, Fvv vs pyscf", h2), &
                  new_unittest("h2o2: Fov, Loo, Lvv vs pyscf", h3), &
                  new_unittest("h2o2: E_MP2 and energy functional vs pyscf", h4), &
                  new_unittest("h2o2: t1 = 0 reductions", h5), &
                  new_unittest("h2o2: intermediates vs host einsum loops", h6), &
                  new_unittest("h2o2: refusals", h7), &
                  new_unittest("h2o2: stale scratch", h8) &
                  ]
   end subroutine collect_check_rccsd_inter

   subroutine w1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 1)
   end subroutine w1

   subroutine w2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 1)
   end subroutine w2

   subroutine w3(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 1)
   end subroutine w3

   subroutine w4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 1)
   end subroutine w4

   subroutine w5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 1)
   end subroutine w5

   subroutine w6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6, 1)
   end subroutine w6

   subroutine w7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7, 1)
   end subroutine w7

   subroutine w8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 1)
   end subroutine w8

   subroutine h1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 2)
   end subroutine h1

   subroutine h2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 2)
   end subroutine h2

   subroutine h3(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 2)
   end subroutine h3

   subroutine h4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 2)
   end subroutine h4

   subroutine h5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 2)
   end subroutine h5

   subroutine h6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 6, 2)
   end subroutine h6

   subroutine h7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 7, 2)
   end subroutine h7

   subroutine h8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 2)
   end subroutine h8

   subroutine run_case(error, ic, im)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic, im
      call ensure_computed()
      call check(error, ok_case(ic, im), trim(msg_case(ic, im)))
   end subroutine run_case

   subroutine ensure_computed()
      integer :: im, dev
      if (computed) return
      dev = trc_bind_device(0)
      ok_case = .false.
      msg_case = "not computed"
      do im = 1, NMOL
         call compute_molecule(im)
      end do
      computed = .true.
   end subroutine ensure_computed

   subroutine record(im, ic, passed, text)
      !! One (ok, message) per case, also printed so a run shows the measured numbers
      integer, intent(in) :: im, ic
      logical, intent(in) :: passed
      character(len=*), intent(in) :: text
      ok_case(ic, im) = passed
      msg_case(ic, im) = text
      print "(a,l1,a,a)", "  check_rccsd_inter ", passed, " ", text
   end subroutine record

   function rel_probe(p, ref) result(w)
      !! Worst |p - ref| / max(|ref|, PROBE_FLOOR)
      real(dp), intent(in) :: p(:), ref(:)
      real(dp) :: w
      w = maxval(abs(p - ref)/max(abs(ref), PROBE_FLOOR))
   end function rel_probe

   function with_diagonal(a, d) result(b)
      !! a + diag(d)
      real(dp), intent(in) :: a(:, :), d(:)
      real(dp), allocatable :: b(:, :)
      integer :: i
      allocate (b(size(a, 1), size(a, 2)))
      b = a
      do i = 1, size(d)
         b(i, i) = b(i, i) + d(i)
      end do
   end function with_diagonal

   function rel_diff(a, b) result(w)
      !! max |a - b| / max(1, max |b|)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp) :: w
      w = maxval(abs(a - b))/max(1.0_dp, maxval(abs(b)))
   end function rel_diff

   subroutine compute_molecule(im)
      !! SCF, the integrals, and every number the cases compare
      integer, intent(in) :: im

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_rccsd_ints_t) :: ints
      type(trc_rccsd_ints_options_t) :: bopts
      type(trc_rccsd_amps_t) :: amps
      type(trc_rccsd_inter_t) :: inter
      type(trc_rccsd_scratch_t) :: scr
      type(trc_linalg_t) :: la
      integer, allocatable :: zint(:)
      real(dp), allocatable :: dip(:, :, :), tauh(:, :, :, :), t2h(:, :, :, :), t1fix(:, :)
      real(dp), allocatable :: foo_fix(:, :), fvv_fix(:, :), fov_fix(:, :), loo_fix(:, :), lvv_fix(:, :)
      real(dp), allocatable :: h_foo(:, :), h_fvv(:, :), h_fov(:, :), h_loo(:, :), h_lvv(:, :)
      real(dp), allocatable :: foo0(:, :), fvv0(:, :), fov0(:, :), loo0(:, :), lvv0(:, :)
      real(dp), allocatable :: zeros(:, :)
      integer :: nocc, nfz, nmo, no, nv, nact, ref_unit, i, a, nocc_ref
      real(dp) :: ref2(3), ref4(4), e_scf_skip, e_mp2_ref, e_fixed_ref, e_mp2, e_s, e_d, e_s2, e_d2, e_mp2_ints
      real(dp) :: d_t1, d_t2, d_tau, d_t2h, d_foo, d_fvv, d_fov, d_loo, d_lvv, d_e1, d_e2, d_e3, d_split, d_rep
      real(dp) :: d_fov0, d_loo0, d_lvv0, d_foo0, d_fvv0
      real(dp) :: es_host, ed_host, ll
      character(len=:), allocatable :: label
      logical :: supported, ok

      label = trim(MOL_LABEL(im))//": "
      nfz = MOL_NFROZEN(im)
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_inter: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      nmo = bas%nao_sph
      call run_scf(bas, nocc, res)
      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = nfz
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_inter: build: "//ints%status%get_message()
      no = ints%no
      nv = ints%nv
      nact = no + nv
      allocate (dip(nact, nact, 4))
      call mo_dipoles(bas, pl, res%cmo(:, nfz + 1:nmo, 1), dip)
      call ref_require(trim(MOL_REF(im)), ref_unit)
      ! The reader is sequential: the tags in the file's order, the integral blocks' probes skipped.
      call ref_read_int(ref_unit, "NO", i)
      call ref_read_int(ref_unit, "NV", a)
      call ref_read_int(ref_unit, "NAUX", nocc_ref)
      call ref_read_int(ref_unit, "NFROZEN", nocc_ref)
      call ref_read_real(ref_unit, "E_SCF", e_scf_skip)
      call ref_read_real(ref_unit, "E_MP2", e_mp2_ref)
      call ref_read_real(ref_unit, "E_FIXED", e_fixed_ref)
      do i = 1, size(SKIP_TAG)
         call ref_read_vector(ref_unit, "P_"//SKIP_TAG(i), 4, ref4)
      end do

      call scr%init(ints)
      if (scr%status%has_error()) error stop "check_rccsd_inter: scratch: "//scr%status%get_message()
      call amps%init(no, nv)
      call inter%init(no, nv)

      ! ---- the amplitudes: t2 = MP2 (the module's guess), t1 from the dipole blocks ----------
      call rccsd_inter_mp2_guess(ints, amps, scr, e_mp2)
      allocate (t1fix(no, nv))
      do a = 1, nv
         do i = 1, no
            t1fix(i, a) = (0.03_dp*dip(i, no + a, 1) + 0.02_dp*dip(i, no + a, 2) + 0.01_dp*dip(i, no + a, 3))/ &
                          (ints%eps_o(i) - ints%eps_v(a))
         end do
      end do
      amps%t1 = t1fix
      call amps%push()

      call ref_read_vector(ref_unit, "P_t1", 3, ref2)
      d_t1 = rel_probe(probe2(amps%t1, "ov", dip, no), ref2)
      call ref_read_vector(ref_unit, "P_t2", 4, ref4)
      d_t2 = rel_probe(probe4(amps%t2, "oovv", dip, no), ref4)
      call host_t2(ints, t2h)
      d_t2h = maxval(abs(amps%t2 - t2h))
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      call inter%pull()
      call host_tau(amps%t1, amps%t2, tauh)
      d_tau = maxval(abs(inter%tau - tauh))
      call record(im, 1, d_t1 < TOL_PROBE .and. d_t2 < TOL_PROBE .and. d_t2h < 1.0e-14_dp .and. d_tau < 1.0e-15_dp &
                  .and. .not. scr%status%has_error(), &
                  label//"t1 probes vs pyscf "//to_string(d_t1, "es9.2")//", t2 probes vs pyscf "// &
                  to_string(d_t2, "es9.2")//", t2 vs host quotient "//to_string(d_t2h, "es9.2")// &
                  ", tau vs host t2 + t1 t1 "//to_string(d_tau, "es9.2"))

      ! ---- the intermediates at the fixed amplitudes ---------------------------------------
      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      call inter%pull()
      if (scr%status%has_error()) error stop "check_rccsd_inter: "//scr%status%get_message()
      foo_fix = inter%foo
      fvv_fix = inter%fvv
      fov_fix = inter%fov
      loo_fix = inter%loo
      lvv_fix = inter%lvv

      call ref_read_vector(ref_unit, "P_Foo", 3, ref2)
      d_foo = rel_probe(probe2(with_diagonal(foo_fix, ints%eps_o), "oo", dip, no), ref2)
      call ref_read_vector(ref_unit, "P_Fvv", 3, ref2)
      d_fvv = rel_probe(probe2(with_diagonal(fvv_fix, ints%eps_v), "vv", dip, no), ref2)
      call record(im, 2, d_foo < TOL_PROBE .and. d_fvv < TOL_PROBE, &
                  label//"Foo (+ eps_o on the diagonal) probes vs pyscf, worst relative (floor 0.1) "// &
                  to_string(d_foo, "es9.2")//", Fvv (+ eps_v) "//to_string(d_fvv, "es9.2"))

      call ref_read_vector(ref_unit, "P_Fov", 3, ref2)
      d_fov = rel_probe(probe2(fov_fix, "ov", dip, no), ref2)
      call ref_read_vector(ref_unit, "P_Loo", 3, ref2)
      d_loo = rel_probe(probe2(with_diagonal(loo_fix, ints%eps_o), "oo", dip, no), ref2)
      call ref_read_vector(ref_unit, "P_Lvv", 3, ref2)
      d_lvv = rel_probe(probe2(with_diagonal(lvv_fix, ints%eps_v), "vv", dip, no), ref2)
      call record(im, 3, d_fov < TOL_PROBE .and. d_loo < TOL_PROBE .and. d_lvv < TOL_PROBE, &
                  label//"Fov probes vs pyscf "//to_string(d_fov, "es9.2")//", Loo (+ eps_o) "// &
                  to_string(d_loo, "es9.2")//", Lvv (+ eps_v) "//to_string(d_lvv, "es9.2"))

      ! ---- energies -----------------------------------------------------------------------
      call rccsd_inter_energy(amps, scr, e_s, e_d)
      call rccsd_inter_energy(amps, scr, e_s2, e_d2)
      e_mp2_ints = ints%mp2_energy()
      d_e1 = abs(e_mp2 - e_mp2_ref)
      d_e2 = abs((e_s + e_d) - e_fixed_ref)
      d_e3 = abs(e_mp2 - e_mp2_ints)
      call host_energy(ints, amps%t1, amps%t2, es_host, ed_host)
      d_split = max(abs(e_s - es_host), abs(e_d - ed_host))
      d_rep = max(abs(e_s - e_s2), abs(e_d - e_d2))
      ok = d_e1 < TOL_ENERGY .and. d_e2 < TOL_ENERGY .and. d_e3 < 1.0e-13_dp .and. d_split < 1.0e-13_dp .and. d_rep == 0.0_dp
      call record(im, 4, ok, label//"E_MP2 vs pyscf "//to_string(d_e1, "es9.2")//" (vs ovov sum "// &
                  to_string(d_e3, "es9.2")//"), E(t1,t2) vs pyscf "//to_string(d_e2, "es9.2")// &
                  ", singles/doubles vs host sum "//to_string(d_split, "es9.2")//", repeat "//to_string(d_rep, "es9.2")// &
                  ", singles "//to_string(e_s, "es12.4"))

      ! ---- t1 = 0 --------------------------------------------------------------------------
      allocate (zeros(no, nv))
      zeros = 0.0_dp
      amps%t1 = zeros
      call amps%push()
      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      call inter%pull()
      foo0 = inter%foo
      fvv0 = inter%fvv
      fov0 = inter%fov
      loo0 = inter%loo
      lvv0 = inter%lvv
      call host_intermediates(ints, amps%t1, amps%t2, h_foo, h_fvv, h_fov, h_loo, h_lvv)
      d_fov0 = maxval(abs(fov0))
      d_foo0 = rel_diff(foo0, h_foo)
      d_fvv0 = rel_diff(fvv0, h_fvv)
      d_loo0 = maxval(abs(loo0 - foo0))
      d_lvv0 = maxval(abs(lvv0 - fvv0))
      ok = d_fov0 == 0.0_dp .and. d_foo0 < TOL_HOST .and. d_fvv0 < TOL_HOST .and. d_loo0 == 0.0_dp .and. d_lvv0 == 0.0_dp
      call record(im, 5, ok, label//"t1 = 0: max|Fov| "//to_string(d_fov0, "es9.2")//", Foo, Fvv vs t2-only host loops "// &
                  to_string(d_foo0, "es9.2")//" "//to_string(d_fvv0, "es9.2")//", max|Loo - Foo| "// &
                  to_string(d_loo0, "es9.2")//", max|Lvv - Fvv| "//to_string(d_lvv0, "es9.2"))

      ! ---- every intermediate at the fixed amplitudes against pyscf-shaped host loops ----------
      call host_intermediates(ints, t1fix, amps%t2, h_foo, h_fvv, h_fov, h_loo, h_lvv)
      d_foo = rel_diff(foo_fix, h_foo)
      d_fvv = rel_diff(fvv_fix, h_fvv)
      d_fov = rel_diff(fov_fix, h_fov)
      d_loo = rel_diff(loo_fix, h_loo)
      d_lvv = rel_diff(lvv_fix, h_lvv)
      ll = max(d_foo, d_fvv, d_fov, d_loo, d_lvv)
      call record(im, 6, ll < TOL_HOST, label//"Foo/Fvv/Fov/Loo/Lvv vs host einsum loops, worst relative "// &
                  to_string(d_foo, "es9.2")//" "//to_string(d_fvv, "es9.2")//" "//to_string(d_fov, "es9.2")//" "// &
                  to_string(d_loo, "es9.2")//" "//to_string(d_lvv, "es9.2"))

      call refusal_case(im, label, ints, amps, inter, scr, la)
      call stale_case(im, label, bas, aux, pl, res, la, bopts, ints, amps, inter, scr, t1fix)

      call ref_close(ref_unit)
      call inter%release()
      call amps%release()
      call scr%release()
      call ints%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine compute_molecule

   subroutine stale_case(im, label, bas, aux, pl, res, la, bopts, ints, amps, inter, scr, t1fix)
      !! After a second `ints%build` the packs of `scr` are stale: the routines that take `ints` refuse
      !! and change nothing; re-initialising the scratch restores the earlier numbers
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_ints_options_t), intent(in) :: bopts
      type(trc_rccsd_ints_t), intent(inout) :: ints
      type(trc_rccsd_amps_t), intent(inout) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      real(dp), intent(in) :: t1fix(:, :)

      real(dp), allocatable :: foo_b(:, :), fvv_b(:, :), fov_b(:, :), loo_b(:, :), lvv_b(:, :)
      real(dp), allocatable :: tau_b(:, :, :, :), t2_b(:, :, :, :), t1_b(:, :)
      integer :: stamp_before
      real(dp) :: d, e_mp2
      character(len=500) :: why

      why = ""
      amps%t1 = t1fix
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      if (scr%status%has_error()) error stop "check_rccsd_inter: stale: "//scr%status%get_message()
      call inter%pull()
      call amps%pull()
      foo_b = inter%foo
      fvv_b = inter%fvv
      fov_b = inter%fov
      loo_b = inter%loo
      lvv_b = inter%lvv
      tau_b = inter%tau
      t1_b = amps%t1
      t2_b = amps%t2
      stamp_before = ints%stamp
      call la%release()
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_inter: rebuild: "//ints%status%get_message()
      if (ints%stamp == stamp_before .or. ints%stamp == 0) why = trim(why)//" [the stamp did not change on rebuild]"

      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      call expect_refusal(scr, "one_particle with stale packs", why)
      call rccsd_inter_foo(ints, inter, scr, la)
      call expect_refusal(scr, "Foo with stale packs", why)
      call rccsd_inter_fvv(ints, inter, scr, la)
      call expect_refusal(scr, "Fvv with stale packs", why)
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      call expect_refusal(scr, "build_tau with stale packs", why)
      call rccsd_inter_mp2_guess(ints, amps, scr, e_mp2)
      call expect_refusal(scr, "mp2_guess with stale packs", why)
      call inter%pull()
      call amps%pull()
      if (any(inter%foo /= foo_b) .or. any(inter%fvv /= fvv_b) .or. any(inter%fov /= fov_b) .or. &
          any(inter%loo /= loo_b) .or. any(inter%lvv /= lvv_b) .or. any(inter%tau /= tau_b)) &
         why = trim(why)//" [a refused call changed the intermediates]"
      if (any(amps%t1 /= t1_b) .or. any(amps%t2 /= t2_b)) why = trim(why)//" [a refused call changed the amplitudes]"

      call scr%init(ints)
      if (scr%status%has_error()) error stop "check_rccsd_inter: stale: scratch: "//scr%status%get_message()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      if (scr%status%has_error()) why = trim(why)//" [one_particle after re-init fails: "//scr%status%get_message()//"]"
      call inter%pull()
      d = max(rel_diff(inter%foo, foo_b), rel_diff(inter%fvv, fvv_b), rel_diff(inter%fov, fov_b), &
              rel_diff(inter%loo, loo_b), rel_diff(inter%lvv, lvv_b))
      if (d > 1.0e-14_dp) why = trim(why)//" [numbers after re-init differ]"
      call record(im, 8, len_trim(why) == 0, label//"after a second ints%build the stale scratch is refused with a "// &
                  "message and nothing changes; after re-init the intermediates reproduce the earlier numbers, worst "// &
                  "relative "//to_string(d, "es9.2")//trim(why))
   end subroutine stale_case

   subroutine expect_refusal(scr, what, why)
      !! Appends to `why` when `scr` does not carry a refusal with a message
      type(trc_rccsd_scratch_t), intent(in) :: scr
      character(len=*), intent(in) :: what
      character(len=*), intent(inout) :: why
      if (.not. scr%status%has_error()) then
         why = trim(why)//" ["//what//" not refused]"
      else if (len_trim(scr%status%get_message()) == 0) then
         why = trim(why)//" ["//what//": refusal has no message]"
      end if
   end subroutine expect_refusal

   subroutine refusal_case(im, label, ints, amps, inter, scr, la)
      !! Refusals set `scr%status` and leave the arrays as they were
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la

      type(trc_rccsd_ints_t) :: empty
      type(trc_rccsd_scratch_t) :: s2
      type(trc_rccsd_amps_t) :: small
      type(trc_rccsd_inter_t) :: small_inter
      real(dp), allocatable :: before(:, :), before_fvv(:, :)
      character(len=500) :: why

      why = ""
      call s2%init(empty)
      if (.not. s2%status%has_error()) why = trim(why)//" [scratch of empty integrals not refused]"
      if (allocated(s2%lov)) why = trim(why)//" [scratch of empty integrals allocated]"
      call small%init(ints%no - 1, ints%nv)
      call small_inter%init(ints%no - 1, ints%nv)
      call inter%pull()
      before = inter%foo
      before_fvv = inter%fvv
      call rccsd_inter_one_particle(ints, small, inter, scr, la)
      if (.not. scr%status%has_error()) why = trim(why)//" [mismatched amplitudes not refused]"
      call inter%pull()
      if (any(inter%foo /= before)) why = trim(why)//" [a refused call changed Foo]"
      ! an inter built for other sizes: Foo and Fvv refuse, and one_particle returns that refusal
      call rccsd_inter_foo(ints, small_inter, scr, la)
      if (.not. scr%status%has_error()) why = trim(why)//" [Foo with a smaller inter not refused]"
      if (len_trim(scr%status%get_message()) == 0) why = trim(why)//" [Foo refusal has no message]"
      call rccsd_inter_fvv(ints, small_inter, scr, la)
      if (.not. scr%status%has_error()) why = trim(why)//" [Fvv with a smaller inter not refused]"
      if (len_trim(scr%status%get_message()) == 0) why = trim(why)//" [Fvv refusal has no message]"
      call rccsd_inter_one_particle(ints, amps, small_inter, scr, la)
      if (.not. scr%status%has_error()) why = trim(why)//" [one_particle with a smaller inter not refused]"
      call inter%pull()
      call small_inter%pull()
      if (any(inter%foo /= before) .or. any(inter%fvv /= before_fvv)) why = trim(why)//" [a refused call changed Foo or Fvv]"
      if (any(small_inter%foo /= 0.0_dp) .or. any(small_inter%fvv /= 0.0_dp) .or. any(small_inter%tau /= 0.0_dp)) &
         why = trim(why)//" [a refused call changed the smaller inter]"
      call rccsd_inter_foo(empty, small_inter, s2, la)
      if (.not. s2%status%has_error()) why = trim(why)//" [Foo without a scratch not refused]"
      call rccsd_inter_one_particle(ints, amps, inter, scr, la)
      if (scr%status%has_error()) why = trim(why)//" [a good call after a refusal still fails]"
      call small_inter%release()
      call small%release()
      call s2%release()
      call record(im, 7, len_trim(why) == 0, label//"empty integrals, mismatched sizes and a missing scratch are refused "// &
                  "with a status and change nothing"//trim(why))
   end subroutine refusal_case

   ! ---- independent host references -----------------------------------------------------------

   subroutine host_t2(ints, t2)
      !! t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), allocatable, intent(out) :: t2(:, :, :, :)
      integer :: i, j, a, b

      allocate (t2(ints%no, ints%no, ints%nv, ints%nv))
      do b = 1, ints%nv
         do a = 1, ints%nv
            do j = 1, ints%no
               do i = 1, ints%no
                  t2(i, j, a, b) = ints%ovov(i, a, j, b)/(ints%eps_o(i) + ints%eps_o(j) - ints%eps_v(a) - ints%eps_v(b))
               end do
            end do
         end do
      end do
   end subroutine host_t2

   subroutine host_tau(t1, t2, tau)
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), allocatable, intent(out) :: tau(:, :, :, :)
      integer :: i, j, a, b

      allocate (tau, mold=t2)
      do b = 1, size(t1, 2)
         do a = 1, size(t1, 2)
            do j = 1, size(t1, 1)
               do i = 1, size(t1, 1)
                  tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
   end subroutine host_tau

   subroutine host_energy(ints, t1, t2, es, ed)
      !! E = sum_ijab [2 (ia|jb) - (ib|ja)] (t2 + t1 t1), split into the t2 part and the t1 t1 part
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: es, ed
      integer :: i, j, a, b
      real(dp) :: l

      es = 0.0_dp
      ed = 0.0_dp
      do b = 1, ints%nv
         do a = 1, ints%nv
            do j = 1, ints%no
               do i = 1, ints%no
                  l = 2.0_dp*ints%ovov(i, a, j, b) - ints%ovov(i, b, j, a)
                  ed = ed + l*t2(i, j, a, b)
                  es = es + l*t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
   end subroutine host_energy

   subroutine host_intermediates(ints, t1, t2, foo, fvv, fov, loo, lvv)
      !! pyscf/cc/rintermediates.py's einsums, term by term (t2 and t1 t1 apart), without the
      !! bare Fock blocks and without the zero fov terms. Foo[k,i], Fvv[a,c], Fov[k,c], Loo[k,i], Lvv[a,c].
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), allocatable, intent(out) :: foo(:, :), fvv(:, :), fov(:, :), loo(:, :), lvv(:, :)

      integer :: no, nv, i, k, l, a, c, d
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (foo(no, no), fvv(nv, nv), fov(no, nv), loo(no, no), lvv(nv, nv))
      ! Fki = 2 kcld,ilcd - kdlc,ilcd + 2 kcld,ic,ld - kdlc,ic,ld
      do i = 1, no
         do k = 1, no
            s = 0.0_dp
            do d = 1, nv
               do c = 1, nv
                  do l = 1, no
                     s = s + 2.0_dp*ints%ovov(k, c, l, d)*t2(i, l, c, d) - ints%ovov(k, d, l, c)*t2(i, l, c, d)
                     s = s + 2.0_dp*ints%ovov(k, c, l, d)*t1(i, c)*t1(l, d) - ints%ovov(k, d, l, c)*t1(i, c)*t1(l, d)
                  end do
               end do
            end do
            foo(k, i) = s
         end do
      end do
      ! Fac = -2 kcld,klad + kdlc,klad - 2 kcld,ka,ld + kdlc,ka,ld
      do c = 1, nv
         do a = 1, nv
            s = 0.0_dp
            do d = 1, nv
               do l = 1, no
                  do k = 1, no
                     s = s - 2.0_dp*ints%ovov(k, c, l, d)*t2(k, l, a, d) + ints%ovov(k, d, l, c)*t2(k, l, a, d)
                     s = s - 2.0_dp*ints%ovov(k, c, l, d)*t1(k, a)*t1(l, d) + ints%ovov(k, d, l, c)*t1(k, a)*t1(l, d)
                  end do
               end do
            end do
            fvv(a, c) = s
         end do
      end do
      ! Fkc = 2 kcld,ld - kdlc,ld
      do c = 1, nv
         do k = 1, no
            s = 0.0_dp
            do d = 1, nv
               do l = 1, no
                  s = s + 2.0_dp*ints%ovov(k, c, l, d)*t1(l, d) - ints%ovov(k, d, l, c)*t1(l, d)
               end do
            end do
            fov(k, c) = s
         end do
      end do
      ! Lki = Fki + 2 lcki,lc - kcli,lc, with ovoo(l,c,k,i) = (lc|ki) = ooov(k,i,l,c)
      do i = 1, no
         do k = 1, no
            s = foo(k, i)
            do c = 1, nv
               do l = 1, no
                  s = s + 2.0_dp*ints%ooov(k, i, l, c)*t1(l, c) - ints%ooov(l, i, k, c)*t1(l, c)
               end do
            end do
            loo(k, i) = s
         end do
      end do
      ! Lac = Fac + 2 kdac,kd - kcad,kd, with ovvv(k,d,a,c) = (kd|ac)
      do c = 1, nv
         do a = 1, nv
            s = fvv(a, c)
            do d = 1, nv
               do k = 1, no
                  s = s + 2.0_dp*ints%ovvv(k, d, a, c)*t1(k, d) - ints%ovvv(k, c, a, d)*t1(k, d)
               end do
            end do
            lvv(a, c) = s
         end do
      end do
   end subroutine host_intermediates

end module check_rccsd_inter_suite

program check_rccsd_inter
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_inter_suite, only: collect_check_rccsd_inter
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_inter, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_inter
