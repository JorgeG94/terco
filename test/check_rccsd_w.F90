module check_rccsd_w_suite
   !!
   !! The two-particle intermediates Woooo, Wvoov and Wvovo of RI-CCSD(T)
   !! (src/post_hf/ccsd/trc_rccsd_w.F90) against the DF-RCCSD of pyscf, on a water and
   !! on a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT, spherical, frozen core 1
   !! and 2 (the molecules, SCF and probes of check_rccsd_ints, shared through
   !! rccsd_test_support).
   !!
   !! The amplitudes are the FIXED, unconverged ones test/pyscf_refs/generate.py
   !! (gen_rccsd_ints) probes: t2 = MP2, t1_ia = (0.03 x + 0.02 y + 0.01 z)_ia / (e_i - e_a).
   !! The reference holds, per tensor, the sign-invariant 4-index probes |A|,
   !! sum A X_pq Y_rs, sum A X_pr Z_qs, sum A Y_ps Z_qr (probe4 of rccsd_test_support),
   !! with the dipole blocks on the spaces of the indices the tensor carries.
   !!
   !! INDEX ORDER. The W arrays of mqc are those of PySCF, index for index: Woooo(k,l,i,j),
   !! Wvoov(a,k,i,c), Wvovo(a,k,c,i) (spaces "oooo", "voov", "vovo"), so the probes
   !! take the arrays of the module as they are.
   !!
   !! Cases, per molecule:
   !!   1 Woooo, Wvoov and Wvovo probes against the reference
   !!   2 each W against host loops written term by term as the einsums of pyscf (t2 and
   !!     t1 t1 terms apart, no tau, no helper of the module), reading only the integral
   !!     blocks and t1, t2
   !!   3 t1 = 0: the t1 terms vanish, each W against the t2-only host form
   !!   4 Woooo(k,l,i,j) = Woooo(l,k,j,i)
   !!   5 refusals: an empty integral type, mismatched sizes (also a scratch initialised for other
   !!     (no, nv)), a missing scratch or tau (status set, message non-empty, nothing changed)
   !!   6 stale scratch: after `ints%build` is called again every routine refuses until the
   !!     scratch is re-initialised, then reproduces the earlier numbers
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t, rccsd_inter_build_tau
   use trc_rccsd_w, only: trc_rccsd_w_t, trc_rccsd_wscratch_t, rccsd_w_woooo, rccsd_w_wvoov, rccsd_w_wvovo, rccsd_w_build
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, TOL_PROBE, PROBE_FLOOR, load_molecule, &
                                 run_scf, mo_dipoles, probe4
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_rccsd_w

   integer, parameter :: NCASE = 6
   character(len=4), parameter :: SKIP_TAG(5) = [character(len=4) :: "oooo", "ooov", "oovv", "ovov", "ovvv"]
   character(len=4), parameter :: SKIP_ONE(7) = [character(len=4) :: "t1", "t2", "Foo", "Fvv", "Fov", "Loo", "Lvv"]
   real(dp), parameter :: TOL_HOST = 1.0e-12_dp  !! against independent host loops, relative to max(1, |W|max)
   real(dp), parameter :: TOL_SYM = 1.0e-13_dp   !! Woooo(k,l,i,j) against Woooo(l,k,j,i), same scale

   logical :: computed = .false.
   logical :: ok_case(NCASE, NMOL)
   character(len=700) :: msg_case(NCASE, NMOL)

contains

   subroutine collect_check_rccsd_w(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: Woooo, Wvoov, Wvovo vs pyscf", w1), &
                  new_unittest("water: W vs host einsum loops", w2), &
                  new_unittest("water: t1 = 0 reductions", w3), &
                  new_unittest("water: Woooo symmetry", w4), &
                  new_unittest("water: refusals", w5), &
                  new_unittest("water: stale scratch", w6), &
                  new_unittest("h2o2: Woooo, Wvoov, Wvovo vs pyscf", h1), &
                  new_unittest("h2o2: W vs host einsum loops", h2), &
                  new_unittest("h2o2: t1 = 0 reductions", h3), &
                  new_unittest("h2o2: Woooo symmetry", h4), &
                  new_unittest("h2o2: refusals", h5), &
                  new_unittest("h2o2: stale scratch", h6) &
                  ]
   end subroutine collect_check_rccsd_w

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
      print "(a,l1,a,a)", "  check_rccsd_w ", passed, " ", text
   end subroutine record

   function rel_probe(p, ref) result(w)
      !! Worst |p - ref| / max(|ref|, PROBE_FLOOR)
      real(dp), intent(in) :: p(:), ref(:)
      real(dp) :: w
      w = maxval(abs(p - ref)/max(abs(ref), PROBE_FLOOR))
   end function rel_probe

   function rel_diff4(a, b) result(w)
      !! max |a - b| / max(1, max |b|)
      real(dp), intent(in) :: a(:, :, :, :), b(:, :, :, :)
      real(dp) :: w
      w = maxval(abs(a - b))/max(1.0_dp, maxval(abs(b)))
   end function rel_diff4

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
      type(trc_rccsd_w_t) :: w
      type(trc_rccsd_wscratch_t) :: ws, ws_other
      type(trc_linalg_t) :: la
      integer, allocatable :: zint(:)
      real(dp), allocatable :: dip(:, :, :), t1fix(:, :), t2h(:, :, :, :)
      real(dp), allocatable :: woooo_fix(:, :, :, :), wvoov_fix(:, :, :, :), wvovo_fix(:, :, :, :)
      real(dp), allocatable :: zeros(:, :)
      integer :: nocc, nfz, nmo, no, nv, nact, ref_unit, i, a, nskip
      real(dp) :: ref4(4), ref2(3), skip_real
      logical :: supported
      character(len=:), allocatable :: label

      label = trim(MOL_LABEL(im))//": "
      nfz = MOL_NFROZEN(im)
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_w: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      nmo = bas%nao_sph
      call run_scf(bas, nocc, res)
      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = nfz
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_w: build: "//ints%status%get_message()
      no = ints%no
      nv = ints%nv
      nact = no + nv
      allocate (dip(nact, nact, 4))
      call mo_dipoles(bas, pl, res%cmo(:, nfz + 1:nmo, 1), dip)

      call ref_require(trim(MOL_REF(im)), ref_unit)
      ! The reader is sequential: the tags in the order of the file, everything before Woooo skipped.
      call ref_read_int(ref_unit, "NO", nskip)
      call ref_read_int(ref_unit, "NV", nskip)
      call ref_read_int(ref_unit, "NAUX", nskip)
      call ref_read_int(ref_unit, "NFROZEN", nskip)
      call ref_read_real(ref_unit, "E_SCF", skip_real)
      call ref_read_real(ref_unit, "E_MP2", skip_real)
      call ref_read_real(ref_unit, "E_FIXED", skip_real)
      do i = 1, size(SKIP_TAG)
         call ref_read_vector(ref_unit, "P_"//SKIP_TAG(i), 4, ref4)
      end do
      do i = 1, size(SKIP_ONE)
         if (SKIP_ONE(i) == "t2") then
            call ref_read_vector(ref_unit, "P_"//trim(SKIP_ONE(i)), 4, ref4)
         else
            call ref_read_vector(ref_unit, "P_"//trim(SKIP_ONE(i)), 3, ref2)
         end if
      end do

      call scr%init(ints)
      if (scr%status%has_error()) error stop "check_rccsd_w: scratch: "//scr%status%get_message()
      call ws%init(ints)
      if (ws%status%has_error()) error stop "check_rccsd_w: W scratch: "//ws%status%get_message()
      call amps%init(no, nv)
      call inter%init(no, nv)
      call w%init(no, nv)

      ! ---- the fixed amplitudes: t2 = MP2 from a host quotient, t1 from the dipole blocks -----
      allocate (t1fix(no, nv))
      do a = 1, nv
         do i = 1, no
            t1fix(i, a) = (0.03_dp*dip(i, no + a, 1) + 0.02_dp*dip(i, no + a, 2) + 0.01_dp*dip(i, no + a, 3))/ &
                          (ints%eps_o(i) - ints%eps_v(a))
         end do
      end do
      call host_t2(ints, t2h)
      amps%t1 = t1fix
      amps%t2 = t2h
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_w: tau: "//scr%status%get_message()
      call rccsd_w_build(ints, amps, inter, w, ws, la)
      if (ws%status%has_error()) error stop "check_rccsd_w: "//ws%status%get_message()
      call w%pull()
      woooo_fix = w%woooo
      wvoov_fix = w%wvoov
      wvovo_fix = w%wvovo

      call probe_case(im, label, ref_unit, dip, no, woooo_fix, wvoov_fix, wvovo_fix)
      call host_case(im, label, ints, t1fix, t2h, woooo_fix, wvoov_fix, wvovo_fix)
      call symmetry_case(im, label, woooo_fix)

      ! ---- t1 = 0 ---------------------------------------------------------------------------
      allocate (zeros(no, nv))
      zeros = 0.0_dp
      amps%t1 = zeros
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_w: t1 = 0: tau: "//scr%status%get_message()
      call rccsd_w_build(ints, amps, inter, w, ws, la)
      if (ws%status%has_error()) error stop "check_rccsd_w: t1 = 0: "//ws%status%get_message()
      call w%pull()
      call zero_case(im, label, ints, t2h, w)

      ! ---- refusals, with the fixed amplitudes back ---------------------------------------
      amps%t1 = t1fix
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_w: refusals: tau: "//scr%status%get_message()
      call rccsd_w_build(ints, amps, inter, w, ws, la)
      if (ws%status%has_error()) error stop "check_rccsd_w: refusals: "//ws%status%get_message()
      call other_scratch(bas, aux, pl, res, la, bopts, ws_other)
      call refusal_case(im, label, ints, amps, inter, w, ws, ws_other, la, woooo_fix, wvoov_fix, wvovo_fix)
      call stale_case(im, label, bas, aux, pl, res, la, bopts, ints, amps, inter, scr, w, ws, woooo_fix, wvoov_fix, wvovo_fix)

      call ref_close(ref_unit)
      call w%release()
      call inter%release()
      call amps%release()
      call ws%release()
      call ws_other%release()
      call scr%release()
      call ints%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()
   end subroutine compute_molecule

   subroutine probe_case(im, label, ref_unit, dip, no, woooo, wvoov, wvovo)
      !! Case 1: the three probe sets against the reference
      integer, intent(in) :: im, ref_unit, no
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: dip(:, :, :), woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :)

      real(dp) :: ref4(4), d_oooo, d_voov, d_vovo

      call ref_read_vector(ref_unit, "P_Woooo", 4, ref4)
      d_oooo = rel_probe(probe4(woooo, "oooo", dip, no), ref4)
      call ref_read_vector(ref_unit, "P_Wvoov", 4, ref4)
      d_voov = rel_probe(probe4(wvoov, "voov", dip, no), ref4)
      call ref_read_vector(ref_unit, "P_Wvovo", 4, ref4)
      d_vovo = rel_probe(probe4(wvovo, "vovo", dip, no), ref4)
      call record(im, 1, d_oooo < TOL_PROBE .and. d_voov < TOL_PROBE .and. d_vovo < TOL_PROBE, &
                  label//"Woooo probes vs pyscf, worst relative (floor 0.1) "//to_string(d_oooo, "es9.2")// &
                  ", Wvoov "//to_string(d_voov, "es9.2")//", Wvovo "//to_string(d_vovo, "es9.2"))
   end subroutine probe_case

   subroutine host_case(im, label, ints, t1, t2, woooo, wvoov, wvovo)
      !! Case 2: every W against pyscf-shaped host loops
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :), woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :)

      real(dp), allocatable :: h_oooo(:, :, :, :), h_voov(:, :, :, :), h_vovo(:, :, :, :)
      real(dp) :: d_oooo, d_voov, d_vovo

      call host_w(ints, t1, t2, h_oooo, h_voov, h_vovo)
      d_oooo = rel_diff4(woooo, h_oooo)
      d_voov = rel_diff4(wvoov, h_voov)
      d_vovo = rel_diff4(wvovo, h_vovo)
      call record(im, 2, max(d_oooo, d_voov, d_vovo) < TOL_HOST, &
                  label//"Woooo/Wvoov/Wvovo vs host einsum loops, worst relative "//to_string(d_oooo, "es9.2")// &
                  " "//to_string(d_voov, "es9.2")//" "//to_string(d_vovo, "es9.2"))
   end subroutine host_case

   subroutine symmetry_case(im, label, woooo)
      !! Case 4: Woooo(k,l,i,j) = Woooo(l,k,j,i), which t2(i,j,a,b) = t2(j,i,b,a) and (kc|ld) = (ld|kc) give
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: woooo(:, :, :, :)

      integer :: i, j, k, l
      real(dp) :: dmax, asym

      dmax = 0.0_dp
      asym = 0.0_dp
      do j = 1, size(woooo, 4)
         do i = 1, size(woooo, 3)
            do l = 1, size(woooo, 2)
               do k = 1, size(woooo, 1)
                  dmax = max(dmax, abs(woooo(k, l, i, j) - woooo(l, k, j, i)))
                  asym = max(asym, abs(woooo(k, l, i, j) - woooo(k, l, j, i)))
               end do
            end do
         end do
      end do
      dmax = dmax/max(1.0_dp, maxval(abs(woooo)))
      ! asym is the unsymmetrised (i,j) swap alone: it must NOT vanish, or the test could not see a lost transposition
      call record(im, 4, dmax < TOL_SYM .and. asym > 1.0e-6_dp, label//"max|Woooo(k,l,i,j) - Woooo(l,k,j,i)| "// &
                  to_string(dmax, "es9.2")//" (the i<->j swap alone differs by "//to_string(asym, "es9.2")//")")
   end subroutine symmetry_case

   subroutine zero_case(im, label, ints, t2, w)
      !! Case 3: with t1 = 0 the W of the module against the t2-only forms
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t2(:, :, :, :)
      type(trc_rccsd_w_t), intent(in) :: w

      real(dp), allocatable :: h_oooo(:, :, :, :), h_voov(:, :, :, :), h_vovo(:, :, :, :)
      real(dp) :: d_oooo, d_voov, d_vovo

      call host_w_t2only(ints, t2, h_oooo, h_voov, h_vovo)
      d_oooo = rel_diff4(w%woooo, h_oooo)
      d_voov = rel_diff4(w%wvoov, h_voov)
      d_vovo = rel_diff4(w%wvovo, h_vovo)
      call record(im, 3, max(d_oooo, d_voov, d_vovo) < TOL_HOST, &
                  label//"t1 = 0: Woooo, Wvoov, Wvovo vs t2-only host loops "//to_string(d_oooo, "es9.2")// &
                  " "//to_string(d_voov, "es9.2")//" "//to_string(d_vovo, "es9.2"))
   end subroutine zero_case

   subroutine refusal_case(im, label, ints, amps, inter, w, ws, ws_other, la, ref_oooo, ref_voov, ref_vovo)
      !! Refusals set `ws%status` with a message and leave W as it was
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_w_t), intent(inout) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      type(trc_rccsd_wscratch_t), intent(inout) :: ws_other  !! initialised for integrals with one active occupied orbital fewer
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: ref_oooo(:, :, :, :), ref_voov(:, :, :, :), ref_vovo(:, :, :, :)

      type(trc_rccsd_ints_t) :: empty
      type(trc_rccsd_wscratch_t) :: s2
      type(trc_rccsd_amps_t) :: small
      type(trc_rccsd_inter_t) :: small_inter, no_tau
      type(trc_rccsd_w_t) :: small_w
      real(dp), allocatable :: b_oooo(:, :, :, :), b_voov(:, :, :, :), b_vovo(:, :, :, :)
      character(len=500) :: why

      why = ""
      call s2%init(empty)
      if (.not. s2%status%has_error()) why = trim(why)//" [scratch of empty integrals not refused]"
      if (len_trim(s2%status%get_message()) == 0) why = trim(why)//" [empty-integral refusal has no message]"
      if (allocated(s2%akl)) why = trim(why)//" [scratch of empty integrals allocated]"
      call small%init(ints%no - 1, ints%nv)
      call small_inter%init(ints%no - 1, ints%nv)
      call small_w%init(ints%no - 1, ints%nv)
      call w%pull()
      b_oooo = w%woooo
      b_voov = w%wvoov
      b_vovo = w%wvovo

      ! mismatched amplitudes, tau and W each refuse every routine
      call rccsd_w_woooo(ints, small, inter, w, ws, la)
      call expect_refusal(ws, "Woooo with smaller amplitudes", why)
      call rccsd_w_wvoov(ints, small, w, ws, la)
      call expect_refusal(ws, "Wvoov with smaller amplitudes", why)
      call rccsd_w_wvovo(ints, small, w, ws, la)
      call expect_refusal(ws, "Wvovo with smaller amplitudes", why)
      call rccsd_w_woooo(ints, amps, small_inter, w, ws, la)
      call expect_refusal(ws, "Woooo with a smaller tau", why)
      call rccsd_w_woooo(ints, amps, no_tau, w, ws, la)
      call expect_refusal(ws, "Woooo without tau", why)
      call rccsd_w_wvoov(ints, amps, small_w, ws, la)
      call expect_refusal(ws, "Wvoov into a smaller W", why)
      call rccsd_w_build(ints, amps, inter, small_w, ws, la)
      call expect_refusal(ws, "build into a smaller W", why)
      ! a scratch initialised for other (no, nv)
      if (ws_other%no == ints%no) why = trim(why)//" [the other scratch has the same no]"
      call rccsd_w_build(ints, amps, inter, w, ws_other, la)
      call expect_refusal(ws_other, "build with a scratch of other (no, nv)", why)
      call rccsd_w_woooo(ints, amps, inter, w, ws_other, la)
      call expect_refusal(ws_other, "Woooo with a scratch of other (no, nv)", why)
      call rccsd_w_wvoov(ints, amps, w, ws_other, la)
      call expect_refusal(ws_other, "Wvoov with a scratch of other (no, nv)", why)
      call rccsd_w_wvovo(ints, amps, w, ws_other, la)
      call expect_refusal(ws_other, "Wvovo with a scratch of other (no, nv)", why)
      ! a scratch that was never initialised, and empty integrals
      call rccsd_w_wvovo(ints, amps, w, s2, la)
      call expect_refusal(s2, "Wvovo without a scratch", why)
      call rccsd_w_build(empty, amps, inter, w, s2, la)
      call expect_refusal(s2, "build with empty integrals", why)
      call w%pull()
      if (any(w%woooo /= b_oooo) .or. any(w%wvoov /= b_voov) .or. any(w%wvovo /= b_vovo)) &
         why = trim(why)//" [a refused call changed W]"
      call small_w%pull()
      if (any(small_w%woooo /= 0.0_dp) .or. any(small_w%wvoov /= 0.0_dp) .or. any(small_w%wvovo /= 0.0_dp)) &
         why = trim(why)//" [a refused call changed the smaller W]"
      call rccsd_w_build(ints, amps, inter, w, ws, la)
      if (ws%status%has_error()) why = trim(why)//" [a good call after a refusal still fails]"
      call w%pull()
      if (rel_diff4(w%woooo, ref_oooo) > 1.0e-14_dp .or. rel_diff4(w%wvoov, ref_voov) > 1.0e-14_dp .or. &
          rel_diff4(w%wvovo, ref_vovo) > 1.0e-14_dp) why = trim(why)//" [the good call after a refusal gave other numbers]"
      call small_w%release()
      call small_inter%release()
      call small%release()
      call s2%release()
      call record(im, 5, len_trim(why) == 0, label//"empty integrals, mismatched sizes, a missing scratch and a missing "// &
                  "tau are refused with a status and change nothing"//trim(why))
   end subroutine refusal_case

   subroutine other_scratch(bas, aux, pl, res, la, bopts, ws_other)
      !! A W scratch for the same molecule with one more frozen orbital, so one active occupied orbital fewer
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_ints_options_t), intent(in) :: bopts
      type(trc_rccsd_wscratch_t), intent(inout) :: ws_other

      type(trc_rccsd_ints_t) :: ints2
      type(trc_rccsd_ints_options_t) :: o2

      o2 = bopts
      o2%nfrozen = bopts%nfrozen + 1
      call la%release()
      call ints2%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, o2)
      if (ints2%status%has_error()) error stop "check_rccsd_w: other build: "//ints2%status%get_message()
      call ws_other%init(ints2)
      if (ws_other%status%has_error()) error stop "check_rccsd_w: other scratch: "//ws_other%status%get_message()
      call ints2%release()
   end subroutine other_scratch

   subroutine stale_case(im, label, bas, aux, pl, res, la, bopts, ints, amps, inter, scr, w, ws, ref_oooo, ref_voov, ref_vovo)
      !! After a second `ints%build` the packs of `ws` are stale: every routine refuses and writes
      !! nothing; re-initialising the scratch restores the earlier numbers
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_ints_options_t), intent(in) :: bopts
      type(trc_rccsd_ints_t), intent(inout) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_rccsd_w_t), intent(inout) :: w
      type(trc_rccsd_wscratch_t), intent(inout) :: ws
      real(dp), intent(in) :: ref_oooo(:, :, :, :), ref_voov(:, :, :, :), ref_vovo(:, :, :, :)

      real(dp), allocatable :: b_oooo(:, :, :, :), b_voov(:, :, :, :), b_vovo(:, :, :, :)
      integer :: stamp_before
      real(dp) :: d
      character(len=500) :: why

      why = ""
      call w%pull()
      b_oooo = w%woooo
      b_voov = w%wvoov
      b_vovo = w%wvovo
      stamp_before = ints%stamp
      call la%release()
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_w: rebuild: "//ints%status%get_message()
      if (ints%stamp == stamp_before .or. ints%stamp == 0) why = trim(why)//" [the stamp did not change on rebuild]"
      call rccsd_w_build(ints, amps, inter, w, ws, la)
      call expect_refusal(ws, "build with stale packs", why)
      call rccsd_w_woooo(ints, amps, inter, w, ws, la)
      call expect_refusal(ws, "Woooo with stale packs", why)
      call rccsd_w_wvoov(ints, amps, w, ws, la)
      call expect_refusal(ws, "Wvoov with stale packs", why)
      call rccsd_w_wvovo(ints, amps, w, ws, la)
      call expect_refusal(ws, "Wvovo with stale packs", why)
      call w%pull()
      if (any(w%woooo /= b_oooo) .or. any(w%wvoov /= b_voov) .or. any(w%wvovo /= b_vovo)) &
         why = trim(why)//" [a refused call changed W]"
      ! re-init of both scratches, then the same numbers (the integrals are rebuilt from the same input)
      call scr%init(ints)
      if (scr%status%has_error()) error stop "check_rccsd_w: stale: scratch: "//scr%status%get_message()
      call ws%init(ints)
      if (ws%status%has_error()) error stop "check_rccsd_w: stale: W scratch: "//ws%status%get_message()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_w: stale: tau: "//scr%status%get_message()
      call rccsd_w_build(ints, amps, inter, w, ws, la)
      if (ws%status%has_error()) why = trim(why)//" [build after re-init fails: "//ws%status%get_message()//"]"
      call w%pull()
      d = max(rel_diff4(w%woooo, ref_oooo), rel_diff4(w%wvoov, ref_voov), rel_diff4(w%wvovo, ref_vovo))
      if (d > 1.0e-14_dp) why = trim(why)//" [numbers after re-init differ]"
      call record(im, 6, len_trim(why) == 0, label//"after a second ints%build the stale scratch is refused with a "// &
                  "message and W is unchanged; after re-init W reproduces the earlier numbers, worst relative "// &
                  to_string(d, "es9.2")//trim(why))
   end subroutine stale_case

   subroutine expect_refusal(ws, what, why)
      !! Appends to `why` when `ws` does not carry a refusal with a message
      type(trc_rccsd_wscratch_t), intent(in) :: ws
      character(len=*), intent(in) :: what
      character(len=*), intent(inout) :: why
      if (.not. ws%status%has_error()) then
         why = trim(why)//" ["//what//" not refused]"
      else if (len_trim(ws%status%get_message()) == 0) then
         why = trim(why)//" ["//what//": refusal has no message]"
      end if
   end subroutine expect_refusal

   subroutine host_t2(ints, t2)
      !! ---- independent host references -----------------------------------------------------------
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

   subroutine host_w(ints, t1, t2, woooo, wvoov, wvovo)
      !! The einsums of pyscf/cc/rintermediates.py cc_Woooo, cc_Wvoov, cc_Wvovo, term by term (t2 and
      !! t1 t1 apart), in the index order and with the integral orderings of pyscf
      !! (ovoo(l,c,k,i) = ooov(k,i,l,c), ovvo(k,c,a,i) = ovov(k,c,i,a)).
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), allocatable, intent(out) :: woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :)
      real(dp), allocatable :: ovvv(:, :, :, :)

      integer :: no, nv, i, j, k, l, a, c, d
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (ovvv(no, nv, nv, nv))
      ! (ia|bc) from the factors: trc_rccsd_ints no longer holds the dense block
      call ints%fill_ovvv(ovvv)
      allocate (woooo(no, no, no, no), wvoov(nv, no, no, nv), wvovo(nv, no, nv, no))
      do j = 1, no
         do i = 1, no
            do l = 1, no
               do k = 1, no
                  s = ints%oooo(k, i, l, j)
                  do c = 1, nv
                     s = s + ints%ooov(k, i, l, c)*t1(j, c)      ! lcki,jc->klij
                     s = s + ints%ooov(l, j, k, c)*t1(i, c)      ! kclj,ic->klij
                  end do
                  do d = 1, nv
                     do c = 1, nv
                        s = s + ints%ovov(k, c, l, d)*t2(i, j, c, d)                ! kcld,ijcd->klij
                        s = s + ints%ovov(k, c, l, d)*t1(i, c)*t1(j, d)             ! kcld,ic,jd->klij
                     end do
                  end do
                  woooo(k, l, i, j) = s
               end do
            end do
         end do
      end do
      do c = 1, nv
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  s = ints%ovov(k, c, i, a)
                  do d = 1, nv
                     s = s + ovvv(k, c, a, d)*t1(i, d)                      ! kcad,id->akic
                  end do
                  do l = 1, no
                     s = s - ints%ooov(l, i, k, c)*t1(l, a)                      ! kcli,la->akic
                  end do
                  do d = 1, nv
                     do l = 1, no
                        s = s - 0.5_dp*ints%ovov(l, d, k, c)*t2(i, l, d, a)      ! ldkc,ilda->akic
                        s = s - 0.5_dp*ints%ovov(l, c, k, d)*t2(i, l, a, d)      ! lckd,ilad->akic
                        s = s - ints%ovov(l, d, k, c)*t1(i, d)*t1(l, a)          ! ldkc,id,la->akic
                        s = s + ints%ovov(l, d, k, c)*t2(i, l, a, d)             ! ldkc,ilad->akic
                     end do
                  end do
                  wvoov(a, k, i, c) = s
               end do
            end do
         end do
      end do
      do i = 1, no
         do c = 1, nv
            do k = 1, no
               do a = 1, nv
                  s = ints%oovv(k, i, a, c)
                  do d = 1, nv
                     s = s + ovvv(k, d, a, c)*t1(i, d)                      ! kdac,id->akci
                  end do
                  do l = 1, no
                     s = s - ints%ooov(k, i, l, c)*t1(l, a)                      ! lcki,la->akci
                  end do
                  do d = 1, nv
                     do l = 1, no
                        s = s - 0.5_dp*ints%ovov(l, c, k, d)*t2(i, l, d, a)      ! lckd,ilda->akci
                        s = s - ints%ovov(l, c, k, d)*t1(i, d)*t1(l, a)          ! lckd,id,la->akci
                     end do
                  end do
                  wvovo(a, k, c, i) = s
               end do
            end do
         end do
      end do
   end subroutine host_w

   subroutine host_w_t2only(ints, t2, woooo, wvoov, wvovo)
      !! The same three with every t1 term dropped (not evaluated at t1 = 0): what is left is
      !! the t2 terms and the bare integrals.
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t2(:, :, :, :)
      real(dp), allocatable, intent(out) :: woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :)

      integer :: no, nv, i, j, k, l, a, c, d
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (woooo(no, no, no, no), wvoov(nv, no, no, nv), wvovo(nv, no, nv, no))
      do j = 1, no
         do i = 1, no
            do l = 1, no
               do k = 1, no
                  s = ints%oooo(k, i, l, j)
                  do d = 1, nv
                     do c = 1, nv
                        s = s + ints%ovov(k, c, l, d)*t2(i, j, c, d)
                     end do
                  end do
                  woooo(k, l, i, j) = s
               end do
            end do
         end do
      end do
      do c = 1, nv
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  s = ints%ovov(k, c, i, a)
                  do d = 1, nv
                     do l = 1, no
                        s = s - 0.5_dp*ints%ovov(l, d, k, c)*t2(i, l, d, a) - 0.5_dp*ints%ovov(l, c, k, d)*t2(i, l, a, d) &
                            + ints%ovov(l, d, k, c)*t2(i, l, a, d)
                     end do
                  end do
                  wvoov(a, k, i, c) = s
               end do
            end do
         end do
      end do
      do i = 1, no
         do c = 1, nv
            do k = 1, no
               do a = 1, nv
                  s = ints%oovv(k, i, a, c)
                  do d = 1, nv
                     do l = 1, no
                        s = s - 0.5_dp*ints%ovov(l, c, k, d)*t2(i, l, d, a)
                     end do
                  end do
                  wvovo(a, k, c, i) = s
               end do
            end do
         end do
      end do
   end subroutine host_w_t2only

end module check_rccsd_w_suite

program check_rccsd_w
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_w_suite, only: collect_check_rccsd_w
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_w, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_w
