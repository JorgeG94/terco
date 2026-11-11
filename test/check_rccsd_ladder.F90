!
! The particle-particle ladder of RI-CCSD(T) and its t1 dressing
! (src/post_hf/ccsd/trc_rccsd_ladder.F90) against the DF-RCCSD of pyscf, on a water and
! on a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT, spherical, frozen core 1 and 2
! (the molecules, SCF and probes of check_rccsd_ints, shared through rccsd_test_support).
!
! The amplitudes are the FIXED, unconverged ones test/pyscf_refs/generate.py
! (gen_rccsd_ints) probes: t2 = MP2, t1_ia = (0.03 x + 0.02 y + 0.01 z)_ia / (e_i - e_a).
! The reference holds three oovv tensors, each as the sign-invariant 4-index probes |A|,
! sum A X_pq Y_rs, sum A X_pr Z_qs, sum A Y_ps Z_qr:
!   P_ladder     pyscf's mycc._add_vvvv(t1, t2, eris) = sum_cd (ac|bd) tau(i,j,c,d), UNDRESSED
!   P_ladder_w   the einsum('abcd,ijcd->ijab', cc_Wvvvv, tau) of update_amps = ladder + dressing
!   P_dressing   the two ovvv terms of cc_Wvvvv contracted with tau
!
! Cases, per molecule (each form = the assembled one, the direct one and the half one, the
! last on both of its paths: `tau_symmetric` true, one contraction, and false, two):
!   1 ladder, ladder + dressing and dressing probes against the reference, every form
!   2 assembled == direct == half (both paths), elementwise
!   3 independence of the batch sizes: assembled cd_batch 1, 7, all; direct p_batch 1, 5, all;
!     half d_batch 1, 2, 3, all on both paths; dressing virtual_batch 1, 3, all
!   4 every form and the dressing against host loops written from the einsums (explicit
!     (ac|bd) from b_vv, explicit Wvvvv dressing terms, no helper of the module)
!   5 t1 = 0: the dressing is exactly zero and the ladder is the t2-only host form
!   6 R(i,j,a,b) = R(j,i,b,a) (t2 = MP2 has the symmetry), and the i<->j swap alone does not
!   7 refusals: empty integrals, mismatched sizes, a missing scratch, result or tau, a scratch of other
!     (no, nv), bad options, then a stale scratch after a second `ints%build` (status set, message
!     non-empty, nothing written; after re-init the earlier numbers come back)
!   8 a deliberately NON-symmetric tau (t2 scaled by 1 + (i - j)/(4 no)): assembled, direct and the
!     general half path give the host ladder of that tau at every d_batch; the symmetric half path
!     gives the host ladder of the symmetrised tau and visibly differs from the true one (the
!     reason `tau_symmetric` is an explicit argument)
!
module check_rccsd_ladder_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_rccsd_inter, only: trc_rccsd_amps_t, trc_rccsd_inter_t, trc_rccsd_scratch_t, rccsd_inter_build_tau
   use trc_rccsd_ladder, only: trc_rccsd_ladder_t, trc_rccsd_lscratch_t, trc_rccsd_ladder_options_t, &
                               rccsd_ladder_pp, rccsd_ladder_dressing, rccsd_ladder_build, &
                               rccsd_ladder_prefers_direct, rccsd_ladder_macs, &
                               RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_DIRECT, RCCSD_LADDER_AUTO, RCCSD_LADDER_HALF
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, TOL_PROBE, PROBE_FLOOR, load_molecule, &
                                 run_scf, mo_dipoles, probe4
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_rccsd_ladder

   integer, parameter :: NCASE = 8
   character(len=4), parameter :: SKIP_TAG(5) = [character(len=4) :: "oooo", "ooov", "oovv", "ovov", "ovvv"]
   character(len=4), parameter :: SKIP_ONE(7) = [character(len=4) :: "t1", "t2", "Foo", "Fvv", "Fov", "Loo", "Lvv"]
   character(len=6), parameter :: SKIP_W(3) = [character(len=6) :: "Woooo", "Wvoov", "Wvovo"]
   real(dp), parameter :: DRESSING_FLOOR = 1.0e-4_dp  !! the dressing is ~1e-4 in size, so its probes are held to this floor, not 0.1
   real(dp), parameter :: TOL_FORMS = 1.0e-11_dp  !! assembled against direct, max abs
   real(dp), parameter :: TOL_BATCH = 1.0e-12_dp  !! any batch size against the default one, max abs
   real(dp), parameter :: TOL_HOST = 1.0e-12_dp   !! against independent host loops, relative to max(1, |R|max)
   real(dp), parameter :: TOL_SYM = 1.0e-13_dp    !! R(i,j,a,b) against R(j,i,b,a), same scale

   logical :: computed = .false.
   logical :: ok_case(NCASE, NMOL)
   character(len=2000) :: msg_case(NCASE, NMOL)

contains

   subroutine collect_check_rccsd_ladder(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: ladder, dressing vs pyscf", w1), &
                  new_unittest("water: assembled == direct", w2), &
                  new_unittest("water: batch-size independence", w3), &
                  new_unittest("water: vs host einsum loops", w4), &
                  new_unittest("water: t1 = 0 reductions", w5), &
                  new_unittest("water: symmetry", w6), &
                  new_unittest("water: refusals and stale scratch", w7), &
                  new_unittest("water: non-symmetric tau", w8), &
                  new_unittest("h2o2: ladder, dressing vs pyscf", h1), &
                  new_unittest("h2o2: assembled == direct", h2), &
                  new_unittest("h2o2: batch-size independence", h3), &
                  new_unittest("h2o2: vs host einsum loops", h4), &
                  new_unittest("h2o2: t1 = 0 reductions", h5), &
                  new_unittest("h2o2: symmetry", h6), &
                  new_unittest("h2o2: refusals and stale scratch", h7), &
                  new_unittest("h2o2: non-symmetric tau", h8) &
                  ]
   end subroutine collect_check_rccsd_ladder

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

   subroutine h8(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 8, 2)
   end subroutine h8

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
      !! One (ok, message) per case, also printed so a run shows the measured numbers. A case recorded again
      !! (the half form's part of it) passes only if both parts do, and its message holds both.
      integer, intent(in) :: im, ic
      logical, intent(in) :: passed
      character(len=*), intent(in) :: text
      if (msg_case(ic, im) == "not computed") then
         ok_case(ic, im) = passed
         msg_case(ic, im) = text
      else
         ok_case(ic, im) = ok_case(ic, im) .and. passed
         msg_case(ic, im) = trim(msg_case(ic, im))//"; "//text
      end if
      print "(a,l1,a,a)", "  check_rccsd_ladder ", passed, " ", text
   end subroutine record

   function rel_probe(p, ref, floor) result(w)
      !! Worst |p - ref| / max(|ref|, floor)
      real(dp), intent(in) :: p(:), ref(:), floor
      real(dp) :: w
      w = maxval(abs(p - ref)/max(abs(ref), floor))
   end function rel_probe

   function rel_diff4(a, b) result(w)
      !! max |a - b| / max(1, max |b|)
      real(dp), intent(in) :: a(:, :, :, :), b(:, :, :, :)
      real(dp) :: w
      w = maxval(abs(a - b))/max(1.0_dp, maxval(abs(b)))
   end function rel_diff4

   function abs_diff4(a, b) result(w)
      !! max |a - b|
      real(dp), intent(in) :: a(:, :, :, :), b(:, :, :, :)
      real(dp) :: w
      w = maxval(abs(a - b))
   end function abs_diff4

   subroutine run_config(ints, amps, inter, la, opts, hpp, hdr, htot, bytes, tau_symmetric)
      !! One scratch of `opts`: the undressed ladder, the dressing alone (on a zeroed result) and the total.
      !! `tau_symmetric` is passed to the ladder as given (absent: the routines' default, general tau).
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_linalg_t), intent(in) :: la
      type(trc_rccsd_ladder_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: hpp(:, :, :, :), hdr(:, :, :, :), htot(:, :, :, :)
      integer(int64), intent(out) :: bytes
      logical, intent(in), optional :: tau_symmetric

      type(trc_rccsd_lscratch_t) :: ls
      type(trc_rccsd_ladder_t) :: lad

      call ls%init(ints, opts)
      if (ls%status%has_error()) error stop "check_rccsd_ladder: scratch: "//ls%status%get_message()
      bytes = ls%bytes
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_pp(ints, inter, lad, ls, la, tau_symmetric)
      if (ls%status%has_error()) error stop "check_rccsd_ladder: pp: "//ls%status%get_message()
      call lad%pull()
      hpp = lad%h
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_dressing(ints, amps, inter, lad, ls, la)
      if (ls%status%has_error()) error stop "check_rccsd_ladder: dressing: "//ls%status%get_message()
      call lad%pull()
      hdr = lad%h
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_build(ints, amps, inter, lad, ls, la, tau_symmetric=tau_symmetric)
      if (ls%status%has_error()) error stop "check_rccsd_ladder: build: "//ls%status%get_message()
      call lad%pull()
      htot = lad%h
      call lad%release()
      call ls%release()
   end subroutine run_config

   function make_options(form, cd_batch, p_batch, virtual_batch) result(opts)
      integer, intent(in) :: form, cd_batch, p_batch, virtual_batch
      type(trc_rccsd_ladder_options_t) :: opts
      opts%form = form
      opts%cd_batch = cd_batch
      opts%p_batch = p_batch
      opts%virtual_batch = virtual_batch
   end function make_options

   function make_half(d_batch, virtual_batch) result(opts)
      !! The half form with `d_batch` virtuals d per pass (0: from the budget)
      integer, intent(in) :: d_batch, virtual_batch
      type(trc_rccsd_ladder_options_t) :: opts
      opts%form = RCCSD_LADDER_HALF
      opts%d_batch = d_batch
      opts%virtual_batch = virtual_batch
   end function make_half

   subroutine skip_to_ladder(ref_unit)
      !! The reader is sequential: the tags in the order of the file, everything before P_ladder skipped
      integer, intent(in) :: ref_unit

      integer :: i, nskip
      real(dp) :: ref4(4), ref2(3), skip_real

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
      do i = 1, size(SKIP_W)
         call ref_read_vector(ref_unit, "P_"//trim(SKIP_W(i)), 4, ref4)
      end do
   end subroutine skip_to_ladder

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
      real(dp), allocatable :: dip(:, :, :), t1fix(:, :), t2h(:, :, :, :), zeros(:, :)
      real(dp), allocatable :: pp_a(:, :, :, :), dr_a(:, :, :, :), tot_a(:, :, :, :)
      real(dp), allocatable :: pp_d(:, :, :, :), dr_d(:, :, :, :), tot_d(:, :, :, :)
      real(dp), allocatable :: pp_s(:, :, :, :), dr_s(:, :, :, :), tot_s(:, :, :, :)
      real(dp), allocatable :: pp_g(:, :, :, :), dr_g(:, :, :, :), tot_g(:, :, :, :)
      integer :: nocc, nfz, nmo, no, nv, nact, ref_unit, i, a
      integer(int64) :: bytes_a, bytes_d, bytes_s, bytes_g
      logical :: supported
      character(len=:), allocatable :: label

      label = trim(MOL_LABEL(im))//": "
      nfz = MOL_NFROZEN(im)
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         print "(a)", "SKIP check_rccsd_ladder: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      nmo = bas%nao_sph
      call run_scf(bas, nocc, res)
      bopts%nalpha = nocc
      bopts%nbeta = nocc
      bopts%nfrozen = nfz
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_ladder: build: "//ints%status%get_message()
      no = ints%no
      nv = ints%nv
      nact = no + nv
      allocate (dip(nact, nact, 4))
      call mo_dipoles(bas, pl, res%cmo(:, nfz + 1:nmo, 1), dip)

      call ref_require(trim(MOL_REF(im)), ref_unit)
      call skip_to_ladder(ref_unit)
      call scr%init(ints)
      if (scr%status%has_error()) error stop "check_rccsd_ladder: scratch: "//scr%status%get_message()
      call amps%init(no, nv)
      call inter%init(no, nv)

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
      if (scr%status%has_error()) error stop "check_rccsd_ladder: tau: "//scr%status%get_message()

      call run_config(ints, amps, inter, la, make_options(RCCSD_LADDER_ASSEMBLED, 0, 0, 0), pp_a, dr_a, tot_a, bytes_a)
      call run_config(ints, amps, inter, la, make_options(RCCSD_LADDER_DIRECT, 0, 0, 0), pp_d, dr_d, tot_d, bytes_d)
      call run_config(ints, amps, inter, la, make_half(0, 0), pp_s, dr_s, tot_s, bytes_s, .true.)
      call run_config(ints, amps, inter, la, make_half(0, 0), pp_g, dr_g, tot_g, bytes_g, .false.)

      call probe_case(im, label, ref_unit, dip, no, pp_a, pp_d, dr_a, dr_d, tot_a, tot_d, pp_s, pp_g, tot_s, tot_g)
      call forms_case(im, label, ints, pp_a, pp_d, dr_a, dr_d, tot_a, tot_d, bytes_a, bytes_d)
      call forms_half(im, label, ints, pp_a, dr_a, tot_a, pp_s, dr_s, tot_s, pp_g, dr_g, tot_g, bytes_a, bytes_s)
      call batch_case(im, label, ints, amps, inter, la, pp_a, dr_a, tot_a, pp_d, tot_d)
      call batch_half(im, ints, amps, inter, la, pp_s, tot_s, pp_g, tot_g)
      call host_case(im, label, ints, t1fix, t2h, pp_a, pp_d, dr_a, dr_d, tot_a, tot_d)
      call host_half(im, ints, t1fix, t2h, pp_s, pp_g, tot_s, tot_g)
      call symmetry_case(im, label, tot_a, tot_d)
      call symmetry_half(im, tot_s, tot_g)

      ! ---- a non-symmetric tau, then t2 = MP2 back --------------------------------------------
      call nonsym_case(im, label, ints, amps, inter, scr, la, t1fix, t2h)
      amps%t2 = t2h

      ! ---- t1 = 0 ---------------------------------------------------------------------------
      allocate (zeros(no, nv))
      zeros = 0.0_dp
      amps%t1 = zeros
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_ladder: t1 = 0: tau: "//scr%status%get_message()
      call zero_case(im, label, ints, amps, inter, la, t2h)

      ! ---- refusals, with the fixed amplitudes back ---------------------------------------
      amps%t1 = t1fix
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_ladder: refusals: tau: "//scr%status%get_message()
      call refusal_case(im, label, bas, aux, pl, res, la, bopts, ints, amps, inter, scr, tot_a, tot_d, tot_s)

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

   subroutine probe_case(im, label, ref_unit, dip, no, pp_a, pp_d, dr_a, dr_d, tot_a, tot_d, pp_s, pp_g, tot_s, tot_g)
      !! Case 1: the three probe sets against the reference, every form (the half form's dressing is the
      !! shared routine, so its ladder and total are probed on both paths)
      integer, intent(in) :: im, ref_unit, no
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: dip(:, :, :)
      real(dp), intent(in) :: pp_a(:, :, :, :), pp_d(:, :, :, :), dr_a(:, :, :, :), dr_d(:, :, :, :)
      real(dp), intent(in) :: tot_a(:, :, :, :), tot_d(:, :, :, :)
      real(dp), intent(in) :: pp_s(:, :, :, :), pp_g(:, :, :, :), tot_s(:, :, :, :), tot_g(:, :, :, :)

      real(dp) :: ref_l(4), ref_w(4), ref_d(4), d_l_a, d_l_d, d_w_a, d_w_d, d_d_a, d_d_d
      real(dp) :: d_l_s, d_l_g, d_w_s, d_w_g

      call ref_read_vector(ref_unit, "P_ladder", 4, ref_l)
      call ref_read_vector(ref_unit, "P_ladder_w", 4, ref_w)
      call ref_read_vector(ref_unit, "P_dressing", 4, ref_d)
      d_l_a = rel_probe(probe4(pp_a, "oovv", dip, no), ref_l, PROBE_FLOOR)
      d_l_d = rel_probe(probe4(pp_d, "oovv", dip, no), ref_l, PROBE_FLOOR)
      d_w_a = rel_probe(probe4(tot_a, "oovv", dip, no), ref_w, PROBE_FLOOR)
      d_w_d = rel_probe(probe4(tot_d, "oovv", dip, no), ref_w, PROBE_FLOOR)
      d_d_a = rel_probe(probe4(dr_a, "oovv", dip, no), ref_d, DRESSING_FLOOR)
      d_d_d = rel_probe(probe4(dr_d, "oovv", dip, no), ref_d, DRESSING_FLOOR)
      call record(im, 1, max(d_l_a, d_l_d, d_w_a, d_w_d, d_d_a, d_d_d) < TOL_PROBE, &
                  label//"vs pyscf, worst relative (floor 0.1, dressing 1e-4): ladder assembled/direct "//to_string(d_l_a, "es9.2")// &
                  " "//to_string(d_l_d, "es9.2")//", ladder+dressing "//to_string(d_w_a, "es9.2")//" "// &
                  to_string(d_w_d, "es9.2")//", dressing "//to_string(d_d_a, "es9.2")//" "//to_string(d_d_d, "es9.2"))
      d_l_s = rel_probe(probe4(pp_s, "oovv", dip, no), ref_l, PROBE_FLOOR)
      d_l_g = rel_probe(probe4(pp_g, "oovv", dip, no), ref_l, PROBE_FLOOR)
      d_w_s = rel_probe(probe4(tot_s, "oovv", dip, no), ref_w, PROBE_FLOOR)
      d_w_g = rel_probe(probe4(tot_g, "oovv", dip, no), ref_w, PROBE_FLOOR)
      call record(im, 1, max(d_l_s, d_l_g, d_w_s, d_w_g) < TOL_PROBE, &
                  "half symmetric/general: ladder "//to_string(d_l_s, "es9.2")//" "//to_string(d_l_g, "es9.2")// &
                  ", ladder+dressing "//to_string(d_w_s, "es9.2")//" "//to_string(d_w_g, "es9.2"))
   end subroutine probe_case

   subroutine forms_half(im, label, ints, pp_a, dr_a, tot_a, pp_s, dr_s, tot_s, pp_g, dr_g, tot_g, bytes_a, bytes_s)
      !! Case 2, the half form: both paths against the assembled form, elementwise
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: pp_a(:, :, :, :), dr_a(:, :, :, :), tot_a(:, :, :, :)
      real(dp), intent(in) :: pp_s(:, :, :, :), dr_s(:, :, :, :), tot_s(:, :, :, :)
      real(dp), intent(in) :: pp_g(:, :, :, :), dr_g(:, :, :, :), tot_g(:, :, :, :)
      integer(int64), intent(in) :: bytes_a, bytes_s

      real(dp) :: d_s, d_g, d_split

      d_s = max(abs_diff4(pp_s, pp_a), abs_diff4(dr_s, dr_a), abs_diff4(tot_s, tot_a))
      d_g = max(abs_diff4(pp_g, pp_a), abs_diff4(dr_g, dr_a), abs_diff4(tot_g, tot_a))
      ! the half form's build is the dressed one (bt and the t1 t1 ovov correction); pp and the dressing are the
      ! separate routines, so build against pp + dressing compares the two routes
      d_split = max(abs_diff4(tot_s, pp_s + dr_s), abs_diff4(tot_g, pp_g + dr_g))
      call record(im, 2, max(d_s, d_g, d_split) < TOL_FORMS .and. bytes_s /= bytes_a, &
                  label//"half vs assembled, max abs over ladder, dressing, total: symmetric path "// &
                  to_string(d_s, "es9.2")//", general path "//to_string(d_g, "es9.2")//"; half dressed build vs "// &
                  "pp + dressing "//to_string(d_split, "es9.2")//" (scratch "// &
                  to_string(bytes_s)//" bytes; half MACs "// &
                  to_string(rccsd_ladder_macs(RCCSD_LADDER_HALF, ints%no, ints%nv, ints%naux), "es9.2")//")")
   end subroutine forms_half

   subroutine batch_half(im, ints, amps, inter, la, pp_s, tot_s, pp_g, tot_g)
      !! Case 3, the half form: d_batch 1 (one d per pass), 2, 3 and nv (all d in one pass) on both paths give the
      !! default run's numbers to round-off
      integer, intent(in) :: im
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: pp_s(:, :, :, :), tot_s(:, :, :, :), pp_g(:, :, :, :), tot_g(:, :, :, :)

      real(dp), allocatable :: hpp(:, :, :, :), hdr(:, :, :, :), htot(:, :, :, :)
      integer(int64) :: bytes
      integer :: k, db(4), vb(4)
      real(dp) :: worst_s, worst_g

      db = [1, 2, 3, ints%nv]
      vb = [1, 3, ints%nv, 2]
      worst_s = 0.0_dp
      worst_g = 0.0_dp
      do k = 1, 4
         call run_config(ints, amps, inter, la, make_half(db(k), vb(k)), hpp, hdr, htot, bytes, .true.)
         worst_s = max(worst_s, abs_diff4(hpp, pp_s), abs_diff4(htot, tot_s))
         call run_config(ints, amps, inter, la, make_half(db(k), vb(k)), hpp, hdr, htot, bytes, .false.)
         worst_g = max(worst_g, abs_diff4(hpp, pp_g), abs_diff4(htot, tot_g))
      end do
      call record(im, 3, max(worst_s, worst_g) < TOL_BATCH, &
                  "half d_batch 1/2/3/"//to_string(ints%nv)//" worst symmetric "//to_string(worst_s, "es9.2")// &
                  ", general "//to_string(worst_g, "es9.2"))
   end subroutine batch_half

   subroutine host_half(im, ints, t1, t2, pp_s, pp_g, tot_s, tot_g)
      !! Case 4, the half form: both paths against the host einsum loops
      integer, intent(in) :: im
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(in) :: pp_s(:, :, :, :), pp_g(:, :, :, :), tot_s(:, :, :, :), tot_g(:, :, :, :)

      real(dp), allocatable :: h_l(:, :, :, :), h_d(:, :, :, :)
      real(dp) :: d(4)

      call host_ladder_and_dressing(ints, t1, t2, .true., h_l, h_d)
      d = [rel_diff4(pp_s, h_l), rel_diff4(pp_g, h_l), rel_diff4(tot_s, h_l + h_d), rel_diff4(tot_g, h_l + h_d)]
      call record(im, 4, maxval(d) < TOL_HOST, &
                  "half symmetric/general: ladder "//to_string(d(1), "es9.2")//" "//to_string(d(2), "es9.2")// &
                  ", total "//to_string(d(3), "es9.2")//" "//to_string(d(4), "es9.2"))
   end subroutine host_half

   subroutine symmetry_half(im, tot_s, tot_g)
      !! Case 6, the half form. The symmetric path is symmetric by construction in the ladder (h = N + N^T), so
      !! the general path is the informative one here.
      integer, intent(in) :: im
      real(dp), intent(in) :: tot_s(:, :, :, :), tot_g(:, :, :, :)

      real(dp) :: sym(2), asym(2)

      call swap_measure(tot_s, sym(1), asym(1))
      call swap_measure(tot_g, sym(2), asym(2))
      call record(im, 6, maxval(sym) < TOL_SYM .and. minval(asym) > 1.0e-6_dp, &
                  "half symmetric/general "//to_string(sym(1), "es9.2")//" "//to_string(sym(2), "es9.2"))
   end subroutine symmetry_half

   subroutine nonsym_case(im, label, ints, amps, inter, scr, la, t1, t2)
      !! Case 8: a tau without the pair symmetry. t2ns(i,j,a,b) = t2(i,j,a,b) (1 + (i - j)/(4 no)) breaks
      !! tau(j,i,d,c) = tau(i,j,c,d). The forms that are exact for any tau (assembled, direct, the general half
      !! path at d_batch 1, 2 and all) must give the host ladder of tau; the symmetric half path must give the host
      !! ladder of the symmetrised tau, whose t2 is (t2ns(i,j,a,b) + t2ns(j,i,b,a))/2 (t1 t1 already has the
      !! symmetry), and must differ visibly from the true ladder. Leaves t2ns in `amps` and its tau in `inter`.
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(inout) :: amps
      type(trc_rccsd_inter_t), intent(inout) :: inter
      type(trc_rccsd_scratch_t), intent(inout) :: scr
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)

      real(dp), allocatable :: t2ns(:, :, :, :), t2s(:, :, :, :), h_l(:, :, :, :), h_d(:, :, :, :)
      real(dp), allocatable :: h_ls(:, :, :, :), h_ds(:, :, :, :)
      real(dp), allocatable :: hpp(:, :, :, :), hdr(:, :, :, :), htot(:, :, :, :)
      integer(int64) :: bytes
      integer :: no, nv, i, j, a, b, k, db(3)
      real(dp) :: d_exact, d_sym, gap, asym_t2, sym, asym_l

      no = ints%no
      nv = ints%nv
      allocate (t2ns(no, no, nv, nv), t2s(no, no, nv, nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2ns(i, j, a, b) = t2(i, j, a, b)*(1.0_dp + real(i - j, dp)/real(4*no, dp))
               end do
            end do
         end do
      end do
      asym_t2 = 0.0_dp
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2s(i, j, a, b) = 0.5_dp*(t2ns(i, j, a, b) + t2ns(j, i, b, a))
                  asym_t2 = max(asym_t2, abs(t2ns(i, j, a, b) - t2ns(j, i, b, a)))
               end do
            end do
         end do
      end do
      amps%t2 = t2ns
      call amps%push()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      if (scr%status%has_error()) error stop "check_rccsd_ladder: non-symmetric tau: "//scr%status%get_message()
      call host_ladder_and_dressing(ints, t1, t2ns, .true., h_l, h_d)
      call host_ladder_and_dressing(ints, t1, t2s, .true., h_ls, h_ds)

      d_exact = 0.0_dp
      call run_config(ints, amps, inter, la, make_options(RCCSD_LADDER_ASSEMBLED, 0, 0, 0), hpp, hdr, htot, bytes)
      d_exact = max(d_exact, rel_diff4(hpp, h_l), rel_diff4(htot, h_l + h_d))
      call run_config(ints, amps, inter, la, make_options(RCCSD_LADDER_DIRECT, 0, 0, 0), hpp, hdr, htot, bytes)
      d_exact = max(d_exact, rel_diff4(hpp, h_l), rel_diff4(htot, h_l + h_d))
      db = [1, 2, nv]
      d_sym = 0.0_dp
      gap = huge(1.0_dp)
      do k = 1, 3
         ! the default (general) path, by omission, then explicitly
         call run_config(ints, amps, inter, la, make_half(db(k), 0), hpp, hdr, htot, bytes)
         d_exact = max(d_exact, rel_diff4(hpp, h_l), rel_diff4(htot, h_l + h_d))
         call run_config(ints, amps, inter, la, make_half(db(k), 0), hpp, hdr, htot, bytes, .false.)
         d_exact = max(d_exact, rel_diff4(hpp, h_l), rel_diff4(htot, h_l + h_d))
         call run_config(ints, amps, inter, la, make_half(db(k), 0), hpp, hdr, htot, bytes, .true.)
         ! the dressed build of the symmetric path is L + D of the symmetrised tau as well
         d_sym = max(d_sym, rel_diff4(hpp, h_ls), rel_diff4(htot, h_ls + h_ds))
         gap = min(gap, maxval(abs(hpp - h_l)))
      end do
      call swap_measure(h_l, sym, asym_l)
      call record(im, 8, d_exact < TOL_HOST .and. d_sym < TOL_HOST .and. gap > 1.0e-6_dp .and. sym > 1.0e-6_dp, &
                  label//"non-symmetric tau (max|t2(i,j,a,b) - t2(j,i,b,a)| "//to_string(asym_t2, "es9.2")// &
                  ", of the host ladder "//to_string(sym, "es9.2")//"): assembled, direct and half general "// &
                  "(d_batch 1/2/"//to_string(nv)//") vs host, worst relative "//to_string(d_exact, "es9.2")// &
                  "; half symmetric vs host ladder (and ladder + dressing) of the symmetrised tau "// &
                  to_string(d_sym, "es9.2")// &
                  ", vs the true ladder at least "//to_string(gap, "es9.2")//" (must be visible)")
   end subroutine nonsym_case

   subroutine forms_case(im, label, ints, pp_a, pp_d, dr_a, dr_d, tot_a, tot_d, bytes_a, bytes_d)
      !! Case 2: assembled against direct, elementwise, for the ladder, the dressing and the total
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: pp_a(:, :, :, :), pp_d(:, :, :, :), dr_a(:, :, :, :), dr_d(:, :, :, :)
      real(dp), intent(in) :: tot_a(:, :, :, :), tot_d(:, :, :, :)
      integer(int64), intent(in) :: bytes_a, bytes_d

      real(dp) :: d_pp, d_dr, d_tot

      d_pp = abs_diff4(pp_a, pp_d)
      d_dr = abs_diff4(dr_a, dr_d)
      d_tot = abs_diff4(tot_a, tot_d)
      ! the two forms must differ in what they do, or this would compare a form with itself
      call record(im, 2, max(d_pp, d_dr, d_tot) < TOL_FORMS .and. maxval(abs(pp_a)) > 1.0e-3_dp .and. &
                  maxval(abs(dr_a)) > 1.0e-6_dp .and. bytes_a /= bytes_d, &
                  label//"assembled vs direct, max abs: ladder "//to_string(d_pp, "es9.2")//", dressing "// &
                  to_string(d_dr, "es9.2")//", total "//to_string(d_tot, "es9.2")//" (|ladder|max "// &
                  to_string(maxval(abs(pp_a)), "es9.2")//", |dressing|max "//to_string(maxval(abs(dr_a)), "es9.2")// &
                  "; scratch "//to_string(bytes_a)//" / "//to_string(bytes_d)//" bytes; auto picks "// &
                  merge("direct   ", "assembled", rccsd_ladder_prefers_direct(ints%no, ints%nv, ints%naux))// &
                  ", MACs "//to_string(rccsd_ladder_macs(RCCSD_LADDER_ASSEMBLED, ints%no, ints%nv, ints%naux), "es9.2")// &
                  " / "//to_string(rccsd_ladder_macs(RCCSD_LADDER_DIRECT, ints%no, ints%nv, ints%naux), "es9.2")//")")
   end subroutine forms_case

   subroutine batch_case(im, label, ints, amps, inter, la, pp_a, dr_a, tot_a, pp_d, tot_d)
      !! Case 3: other batch sizes give the default run's numbers to round-off
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: pp_a(:, :, :, :), dr_a(:, :, :, :), tot_a(:, :, :, :), pp_d(:, :, :, :), tot_d(:, :, :, :)

      real(dp), allocatable :: hpp(:, :, :, :), hdr(:, :, :, :), htot(:, :, :, :)
      integer(int64) :: bytes
      integer :: k, nv2, cdb(3), pb(3), vb(3)
      real(dp) :: worst_a, worst_d, worst_v, d
      character(len=300) :: seen

      nv2 = ints%nv*ints%nv
      cdb = [1, 7, nv2]
      pb = [1, 5, ints%naux]
      vb = [1, 3, ints%nv]
      worst_a = 0.0_dp
      worst_d = 0.0_dp
      worst_v = 0.0_dp
      seen = ""
      do k = 1, 3
         call run_config(ints, amps, inter, la, make_options(RCCSD_LADDER_ASSEMBLED, cdb(k), 0, vb(k)), hpp, hdr, htot, bytes)
         d = max(abs_diff4(hpp, pp_a), abs_diff4(htot, tot_a))
         worst_a = max(worst_a, d)
         worst_v = max(worst_v, abs_diff4(hdr, dr_a))
         seen = trim(seen)//" "//to_string(bytes)
         call run_config(ints, amps, inter, la, make_options(RCCSD_LADDER_DIRECT, 0, pb(k), vb(k)), hpp, hdr, htot, bytes)
         d = max(abs_diff4(hpp, pp_d), abs_diff4(htot, tot_d))
         worst_d = max(worst_d, d)
         worst_v = max(worst_v, abs_diff4(hdr, dr_a))
      end do
      call record(im, 3, max(worst_a, worst_d, worst_v) < TOL_BATCH, &
                  label//"batch sizes: assembled cd_batch 1/7/"//to_string(nv2)//" worst "//to_string(worst_a, "es9.2")// &
                  ", direct p_batch 1/5/"//to_string(ints%naux)//" worst "//to_string(worst_d, "es9.2")// &
                  ", dressing virtual_batch 1/3/"//to_string(ints%nv)//" worst "//to_string(worst_v, "es9.2"))
   end subroutine batch_case

   subroutine host_case(im, label, ints, t1, t2, pp_a, pp_d, dr_a, dr_d, tot_a, tot_d)
      !! Case 4: the undressed ladder and the dressing against host loops written from the einsums
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(in) :: pp_a(:, :, :, :), pp_d(:, :, :, :), dr_a(:, :, :, :), dr_d(:, :, :, :)
      real(dp), intent(in) :: tot_a(:, :, :, :), tot_d(:, :, :, :)

      real(dp), allocatable :: h_l(:, :, :, :), h_d(:, :, :, :)
      real(dp) :: d_la, d_ld, d_da, d_dd, d_ta, d_td

      call host_ladder_and_dressing(ints, t1, t2, .true., h_l, h_d)
      d_la = rel_diff4(pp_a, h_l)
      d_ld = rel_diff4(pp_d, h_l)
      d_da = rel_diff4(dr_a, h_d)
      d_dd = rel_diff4(dr_d, h_d)
      d_ta = rel_diff4(tot_a, h_l + h_d)
      d_td = rel_diff4(tot_d, h_l + h_d)
      call record(im, 4, max(d_la, d_ld, d_da, d_dd, d_ta, d_td) < TOL_HOST, &
                  label//"vs host einsum loops, worst relative: ladder assembled/direct "//to_string(d_la, "es9.2")// &
                  " "//to_string(d_ld, "es9.2")//", dressing "//to_string(d_da, "es9.2")//" "//to_string(d_dd, "es9.2")// &
                  ", total "//to_string(d_ta, "es9.2")//" "//to_string(d_td, "es9.2"))
   end subroutine host_case

   subroutine zero_case(im, label, ints, amps, inter, la, t2)
      !! Case 5: with t1 = 0 (amps and tau already so) the dressing vanishes and the ladder is the t2-only form
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_linalg_t), intent(in) :: la
      real(dp), intent(in) :: t2(:, :, :, :)

      real(dp), allocatable :: pp(:, :, :, :), dr(:, :, :, :), tot(:, :, :, :), h_l(:, :, :, :), h_d(:, :, :, :)
      real(dp), allocatable :: zero_t1(:, :)
      integer(int64) :: bytes
      real(dp) :: dmax_dr, d_l(4), d_t(4), lmax
      integer :: k, forms(4)
      logical :: sym(4)

      allocate (zero_t1(ints%no, ints%nv))
      zero_t1 = 0.0_dp
      call host_ladder_and_dressing(ints, zero_t1, t2, .false., h_l, h_d)
      forms = [RCCSD_LADDER_ASSEMBLED, RCCSD_LADDER_DIRECT, RCCSD_LADDER_HALF, RCCSD_LADDER_HALF]
      sym = [.false., .false., .true., .false.]
      dmax_dr = 0.0_dp
      do k = 1, 4
         call run_config(ints, amps, inter, la, make_options(forms(k), 0, 0, 0), pp, dr, tot, bytes, sym(k))
         dmax_dr = max(dmax_dr, maxval(abs(dr)))
         d_l(k) = rel_diff4(pp, h_l)
         d_t(k) = abs_diff4(tot, pp)
      end do
      lmax = maxval(abs(h_l))
      call record(im, 5, dmax_dr == 0.0_dp .and. maxval(d_l) < TOL_HOST .and. maxval(d_t) == 0.0_dp .and. lmax > 1.0e-3_dp, &
                  label//"t1 = 0: max|dressing| "//to_string(dmax_dr, "es9.2")//" (must be 0), ladder vs t2-only host loop "// &
                  "assembled/direct/half symmetric/half general "//to_string(d_l(1), "es9.2")//" "// &
                  to_string(d_l(2), "es9.2")//" "//to_string(d_l(3), "es9.2")//" "//to_string(d_l(4), "es9.2")// &
                  ", total - ladder "// &
                  to_string(maxval(d_t), "es9.2")//" (must be 0), |ladder|max "//to_string(lmax, "es9.2"))
   end subroutine zero_case

   subroutine symmetry_case(im, label, tot_a, tot_d)
      !! Case 6: R(i,j,a,b) = R(j,i,b,a) for the total of each form (t2 = MP2 and t1 t1 are symmetric under
      !! the pair swap); the i<->j swap alone must not leave R unchanged, or the test could not see a lost swap
      integer, intent(in) :: im
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: tot_a(:, :, :, :), tot_d(:, :, :, :)

      real(dp) :: sym(2), asym(2)

      call swap_measure(tot_a, sym(1), asym(1))
      call swap_measure(tot_d, sym(2), asym(2))
      call record(im, 6, maxval(sym) < TOL_SYM .and. minval(asym) > 1.0e-6_dp, &
                  label//"max|R(i,j,a,b) - R(j,i,b,a)| assembled/direct "//to_string(sym(1), "es9.2")//" "// &
                  to_string(sym(2), "es9.2")//" (the i<->j swap alone differs by "//to_string(asym(1), "es9.2")//")")
   end subroutine symmetry_case

   subroutine swap_measure(r, sym, asym)
      real(dp), intent(in) :: r(:, :, :, :)
      real(dp), intent(out) :: sym, asym

      integer :: i, j, a, b

      sym = 0.0_dp
      asym = 0.0_dp
      do b = 1, size(r, 4)
         do a = 1, size(r, 3)
            do j = 1, size(r, 2)
               do i = 1, size(r, 1)
                  sym = max(sym, abs(r(i, j, a, b) - r(j, i, b, a)))
                  asym = max(asym, abs(r(i, j, a, b) - r(j, i, a, b)))
               end do
            end do
         end do
      end do
      sym = sym/max(1.0_dp, maxval(abs(r)))
   end subroutine swap_measure

   subroutine expect_refusal(ls, what, why)
      !! Appends to `why` when `ls` does not carry a refusal with a message
      type(trc_rccsd_lscratch_t), intent(in) :: ls
      character(len=*), intent(in) :: what
      character(len=*), intent(inout) :: why
      if (.not. ls%status%has_error()) then
         why = trim(why)//" ["//what//" not refused]"
      else if (len_trim(ls%status%get_message()) == 0) then
         why = trim(why)//" ["//what//": refusal has no message]"
      end if
   end subroutine expect_refusal

   subroutine refuse_all(ints, amps, inter, lad, ls, la, what, why)
      !! The three routines on arguments that must be refused; every one must set `ls%status`
      type(trc_rccsd_ints_t), intent(in) :: ints
      type(trc_rccsd_amps_t), intent(in) :: amps
      type(trc_rccsd_inter_t), intent(in) :: inter
      type(trc_rccsd_ladder_t), intent(inout) :: lad
      type(trc_rccsd_lscratch_t), intent(inout) :: ls
      type(trc_linalg_t), intent(in) :: la
      character(len=*), intent(in) :: what
      character(len=*), intent(inout) :: why

      call rccsd_ladder_pp(ints, inter, lad, ls, la)
      call expect_refusal(ls, "pp "//what, why)
      call rccsd_ladder_dressing(ints, amps, inter, lad, ls, la)
      call expect_refusal(ls, "dressing "//what, why)
      call rccsd_ladder_build(ints, amps, inter, lad, ls, la)
      call expect_refusal(ls, "build "//what, why)
   end subroutine refuse_all

   subroutine refusal_case(im, label, bas, aux, pl, res, la, bopts, ints, amps, inter, scr, tot_a, tot_d, tot_s)
      !! Case 7: refusals set `ls%status` with a message and write nothing; a stale scratch is refused after a
      !! second `ints%build` and, re-initialised, reproduces the earlier numbers
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
      real(dp), intent(in) :: tot_a(:, :, :, :), tot_d(:, :, :, :), tot_s(:, :, :, :)

      type(trc_rccsd_ints_t) :: empty
      type(trc_rccsd_lscratch_t) :: ls_a, ls_d, ls_h, ls_none, ls_other, ls_bad
      type(trc_rccsd_ladder_t) :: lad, lad_small, lad_none
      type(trc_rccsd_amps_t) :: small, no_amps
      type(trc_rccsd_inter_t) :: small_inter, no_tau
      real(dp), allocatable :: before_a(:, :, :, :), before_stale(:, :, :, :)
      real(dp) :: d_a, d_d, d_h
      integer :: stamp_before
      character(len=700) :: why
      type(trc_rccsd_ladder_options_t) :: bad

      why = ""
      call ls_a%init(ints, make_options(RCCSD_LADDER_ASSEMBLED, 11, 0, 5))
      call ls_d%init(ints, make_options(RCCSD_LADDER_DIRECT, 0, 4, 5))
      call ls_h%init(ints, make_half(2, 5))
      if (ls_a%status%has_error() .or. ls_d%status%has_error() .or. ls_h%status%has_error()) &
         error stop "check_rccsd_ladder: refusals: scratch"
      call lad%init(ints%no, ints%nv)
      call rccsd_ladder_build(ints, amps, inter, lad, ls_a, la)
      if (ls_a%status%has_error()) error stop "check_rccsd_ladder: refusals: build: "//ls_a%status%get_message()
      call lad%pull()
      before_a = lad%h
      call small%init(ints%no - 1, ints%nv)
      call small_inter%init(ints%no - 1, ints%nv)
      call lad_small%init(ints%no - 1, ints%nv)
      call other_scratch(bas, aux, pl, res, la, bopts, ls_other)

      ! scratch, result, tau and amplitudes that are missing or of other sizes, for both forms
      call refuse_all(empty, amps, inter, lad, ls_a, la, "with empty integrals", why)
      call refuse_all(ints, amps, inter, lad_small, ls_a, la, "into a smaller result", why)
      call refuse_all(ints, amps, inter, lad_none, ls_a, la, "into an empty result", why)
      call refuse_all(ints, amps, small_inter, lad, ls_d, la, "with a smaller tau", why)
      call refuse_all(ints, amps, no_tau, lad, ls_d, la, "without tau", why)
      call refuse_all(ints, amps, small_inter, lad, ls_h, la, "with a smaller tau (half)", why)
      call refuse_all(ints, amps, inter, lad_small, ls_h, la, "into a smaller result (half)", why)
      call refuse_all(ints, amps, inter, lad, ls_none, la, "without a scratch", why)
      if (ls_other%no == ints%no) why = trim(why)//" [the other scratch has the same no]"
      call refuse_all(ints, amps, inter, lad, ls_other, la, "with a scratch of other (no, nv)", why)
      ! a smaller or missing t1 refuses the dressing and the build, not the undressed ladder
      call rccsd_ladder_dressing(ints, small, inter, lad, ls_a, la)
      call expect_refusal(ls_a, "dressing with smaller amplitudes", why)
      call rccsd_ladder_build(ints, small, inter, lad, ls_d, la)
      call expect_refusal(ls_d, "build with smaller amplitudes", why)
      call rccsd_ladder_dressing(ints, no_amps, inter, lad, ls_a, la)
      call expect_refusal(ls_a, "dressing without amplitudes", why)
      ! options and integrals the scratch itself refuses
      call ls_bad%init(empty)
      call expect_refusal(ls_bad, "scratch for empty integrals", why)
      if (allocated(ls_bad%taut)) why = trim(why)//" [scratch of empty integrals allocated]"
      bad = make_options(7, 0, 0, 0)
      call ls_bad%init(ints, bad)
      call expect_refusal(ls_bad, "scratch with an unknown form", why)
      bad = make_options(RCCSD_LADDER_AUTO, -1, 0, 0)
      call ls_bad%init(ints, bad)
      call expect_refusal(ls_bad, "scratch with a negative batch", why)
      bad = make_half(-1, 0)
      call ls_bad%init(ints, bad)
      call expect_refusal(ls_bad, "scratch with a negative d_batch", why)
      if (allocated(ls_bad%taut) .or. allocated(ls_bad%gtau) .or. allocated(ls_bad%wbatch) .or. &
          allocated(ls_bad%cbuf) .or. allocated(ls_bad%w2)) &
         why = trim(why)//" [a refused scratch allocated]"
      ! nothing written by any of it
      call lad%pull()
      if (any(lad%h /= before_a)) why = trim(why)//" [a refused call changed the result]"
      call lad_small%pull()
      if (any(lad_small%h /= 0.0_dp)) why = trim(why)//" [a refused call changed the smaller result]"
      call rccsd_ladder_build(ints, amps, inter, lad, ls_a, la)
      if (ls_a%status%has_error()) why = trim(why)//" [a good call after a refusal still fails]"
      call lad%pull()
      if (abs_diff4(lad%h, before_a) > 0.0_dp) why = trim(why)//" [the good call after a refusal gave other numbers]"

      ! ---- stale scratch -------------------------------------------------------------------
      before_stale = before_a
      stamp_before = ints%stamp
      call la%release()
      call ints%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
      if (ints%status%has_error()) error stop "check_rccsd_ladder: rebuild: "//ints%status%get_message()
      if (ints%stamp == stamp_before .or. ints%stamp == 0) why = trim(why)//" [the stamp did not change on rebuild]"
      call refuse_all(ints, amps, inter, lad, ls_a, la, "with stale assembled scratch", why)
      call refuse_all(ints, amps, inter, lad, ls_d, la, "with stale direct scratch", why)
      call refuse_all(ints, amps, inter, lad, ls_h, la, "with stale half scratch", why)
      call lad%pull()
      if (any(lad%h /= before_stale)) why = trim(why)//" [a stale call changed the result]"
      call scr%init(ints)
      if (scr%status%has_error()) error stop "check_rccsd_ladder: stale: scratch: "//scr%status%get_message()
      call rccsd_inter_build_tau(ints, amps, inter, scr)
      call ls_a%init(ints, make_options(RCCSD_LADDER_ASSEMBLED, 0, 0, 0))
      call ls_d%init(ints, make_options(RCCSD_LADDER_DIRECT, 0, 0, 0))
      call rccsd_ladder_build(ints, amps, inter, lad, ls_a, la)
      if (ls_a%status%has_error()) why = trim(why)//" [build after re-init fails: "//ls_a%status%get_message()//"]"
      call lad%pull()
      d_a = abs_diff4(lad%h, tot_a)
      call rccsd_ladder_build(ints, amps, inter, lad, ls_d, la)
      if (ls_d%status%has_error()) why = trim(why)//" [direct build after re-init fails: "//ls_d%status%get_message()//"]"
      call lad%pull()
      d_d = abs_diff4(lad%h, tot_d)
      call ls_h%init(ints, make_half(0, 0))
      call rccsd_ladder_build(ints, amps, inter, lad, ls_h, la, tau_symmetric=.true.)
      if (ls_h%status%has_error()) why = trim(why)//" [half build after re-init fails: "//ls_h%status%get_message()//"]"
      call lad%pull()
      d_h = abs_diff4(lad%h, tot_s)
      if (max(d_a, d_d, d_h) > 1.0e-14_dp) why = trim(why)//" [numbers after re-init differ]"
      call lad%release()
      call lad_small%release()
      call small_inter%release()
      call small%release()
      call ls_a%release()
      call ls_d%release()
      call ls_h%release()
      call ls_other%release()
      call record(im, 7, len_trim(why) == 0, label//"empty integrals, mismatched sizes, a missing scratch, result, "// &
                  "tau and t1, bad options and a stale scratch are refused with a status and change nothing; after "// &
                  "re-init the earlier numbers return, worst abs "//to_string(max(d_a, d_d, d_h), "es9.2")//trim(why))
   end subroutine refusal_case

   subroutine other_scratch(bas, aux, pl, res, la, bopts, ls_other)
      !! A scratch for the same molecule with one more frozen orbital, so one active occupied orbital fewer
      type(trc_basis_t), intent(in) :: bas, aux
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_scf_result_t), intent(in) :: res
      type(trc_linalg_t), intent(inout) :: la
      type(trc_rccsd_ints_options_t), intent(in) :: bopts
      type(trc_rccsd_lscratch_t), intent(inout) :: ls_other

      type(trc_rccsd_ints_t) :: ints2
      type(trc_rccsd_ints_options_t) :: o2

      o2 = bopts
      o2%nfrozen = bopts%nfrozen + 1
      call la%release()
      call ints2%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, o2)
      if (ints2%status%has_error()) error stop "check_rccsd_ladder: other build: "//ints2%status%get_message()
      call ls_other%init(ints2)
      if (ls_other%status%has_error()) error stop "check_rccsd_ladder: other scratch: "//ls_other%status%get_message()
      call ints2%release()
   end subroutine other_scratch

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

   subroutine host_ladder_and_dressing(ints, t1, t2, with_t1, h_l, h_d)
      !! The einsums of pyscf, term by term, from the integral blocks and t1, t2 only (no tau of the module,
      !! no helper): (ac|bd) = sum_P b_vv(a,c;P) b_vv(b,d;P) formed explicitly,
      !!   h_l(i,j,a,b) = sum_cd (ac|bd) tau(i,j,c,d)                         'ijcd,acbd->ijab'
      !!   h_d(i,j,a,b) = sum_cd Wd(a,b,c,d) tau(i,j,c,d)
      !!   Wd(a,b,c,d)  = - sum_k ovvv(k,d,a,c) t1(k,b) - sum_k ovvv(k,c,b,d) t1(k,a)
      !! with tau = t2 + t1 t1 when `with_t1` and tau = t2 else (then h_d is not evaluated and is zero)
      type(trc_rccsd_ints_t), intent(in) :: ints
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      logical, intent(in) :: with_t1
      real(dp), allocatable, intent(out) :: h_l(:, :, :, :), h_d(:, :, :, :)
      real(dp), allocatable :: ovvv(:, :, :, :)

      real(dp), allocatable :: vvvv(:, :, :, :), wd(:, :, :, :), tau(:, :, :, :)
      integer :: no, nv, i, j, a, b, c, d, k, p
      real(dp) :: s

      no = ints%no
      nv = ints%nv
      allocate (ovvv(no, nv, nv, nv))
      ! (ia|bc) from the factors: trc_rccsd_ints no longer holds the dense block
      call ints%fill_ovvv(ovvv)
      allocate (vvvv(nv, nv, nv, nv), wd(nv, nv, nv, nv), tau(no, no, nv, nv))
      allocate (h_l(no, no, nv, nv), h_d(no, no, nv, nv))
      do d = 1, nv
         do b = 1, nv
            do c = 1, nv
               do a = 1, nv
                  s = 0.0_dp
                  do p = 1, ints%naux
                     s = s + ints%b_vv(a + nv*(c - 1), p)*ints%b_vv(b + nv*(d - 1), p)
                  end do
                  vvvv(a, c, b, d) = s
               end do
            end do
         end do
      end do
      do d = 1, nv
         do c = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, c, d) = t2(i, j, c, d)
                  if (with_t1) tau(i, j, c, d) = tau(i, j, c, d) + t1(i, c)*t1(j, d)
               end do
            end do
         end do
      end do
      do d = 1, nv
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  s = 0.0_dp
                  if (with_t1) then
                     do k = 1, no
                        s = s - ovvv(k, d, a, c)*t1(k, b) - ovvv(k, c, b, d)*t1(k, a)
                     end do
                  end if
                  wd(a, b, c, d) = s
               end do
            end do
         end do
      end do
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  h_l(i, j, a, b) = 0.0_dp
                  h_d(i, j, a, b) = 0.0_dp
                  do d = 1, nv
                     do c = 1, nv
                        h_l(i, j, a, b) = h_l(i, j, a, b) + vvvv(a, c, b, d)*tau(i, j, c, d)
                        h_d(i, j, a, b) = h_d(i, j, a, b) + wd(a, b, c, d)*tau(i, j, c, d)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine host_ladder_and_dressing

end module check_rccsd_ladder_suite

program check_rccsd_ladder
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_ladder_suite, only: collect_check_rccsd_ladder
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_ladder, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_ladder
