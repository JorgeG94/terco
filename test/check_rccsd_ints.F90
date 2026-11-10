!
! The fitted integral blocks of RI-CCSD(T) (src/post_hf/ccsd/trc_rccsd_ints.F90)
! against the DF-RCCSD integrals of pyscf, on a water and on a C1 hydrogen peroxide
! (no probe vanishes by symmetry), cc-pVDZ with cc-pVDZ-RIFIT, spherical,
! frozen core 1 and 2.
!
! pyscf never runs here: test/pyscf_refs/generate.py (case gen_rccsd_ints) wrote
! rccsd_ints_<molecule>.ref. A 4-index block cannot be compared element by
! element across the orbital phases of two codes, so the reference holds sign-
! invariant probes of each block A(p,q,r,s): its Frobenius norm and three
! contractions with the MO dipole integrals (origin 0) over the active
! orbitals,
!     sum A X_pq Y_rs,   sum A X_pr Z_qs,   sum A Y_ps Z_qr,
! each dipole block restricted to the orbital spaces (o or v) of the indices
! it carries. A transposed or mis-packed block changes them.
!
! Each molecule is computed once (SCF to energy 1e-13 and commutator 1e-9,
! since RI-MP2 is 3e-9 off at the default thresholds) and every case below
! only reads what was recorded.
!
module check_rccsd_ints_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_result_t
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t, trc_rimp2_options_t, rimp2_driver_mo_to_cart
   use trc_rccsd_ints, only: trc_rccsd_ints_t, trc_rccsd_ints_options_t
   use trc_linalg, only: trc_linalg_t
   use rccsd_test_support, only: NMOL, MOL_LABEL, MOL_REF, MOL_NFROZEN, TOL_PROBE, PROBE_FLOOR, load_molecule, run_scf, &
                                 mo_dipoles, probe4
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector
   implicit none
   private
   public :: collect_check_rccsd_ints

   integer, parameter :: NCASE = 5
   integer, parameter :: NBLOCK = 5
   character(len=4), parameter :: BLOCK_NAME(NBLOCK) = [character(len=4) :: "oooo", "ooov", "oovv", "ovov", "ovvv"]

   real(dp), parameter :: TOL_SCF = 1.0e-10_dp     !! E_SCF against pyscf
   real(dp), parameter :: TOL_MP2 = 1.0e-10_dp     !! E_MP2 against pyscf
   real(dp), parameter :: TOL_RIMP2 = 1.0e-11_dp   !! E_MP2 against trc_rimp2_run, same orbitals
   real(dp), parameter :: TOL_ALGEBRA = 1.0e-12_dp

   logical :: computed = .false.
   logical :: ok_case(NCASE, NMOL)
   character(len=500) :: msg_case(NCASE, NMOL)
   logical :: refusals_ok = .false.
   character(len=500) :: refusals_msg = "not run"

contains

   subroutine collect_check_rccsd_ints(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water: E_SCF vs pyscf", w_scf), &
                  new_unittest("water: E_MP2 from ovov vs pyscf and RI-MP2", w_mp2), &
                  new_unittest("water: block probes vs pyscf", w_probes), &
                  new_unittest("water: B consistency and aux blocking", w_algebra), &
                  new_unittest("water: block symmetries", w_symmetry), &
                  new_unittest("h2o2: E_SCF vs pyscf", h_scf), &
                  new_unittest("h2o2: E_MP2 from ovov vs pyscf and RI-MP2", h_mp2), &
                  new_unittest("h2o2: block probes vs pyscf", h_probes), &
                  new_unittest("h2o2: B consistency and aux blocking", h_algebra), &
                  new_unittest("h2o2: block symmetries", h_symmetry), &
                  new_unittest("refusals: ECP, ranks, open shell, no occupied, no virtual, frozen, shapes", refusals) &
                  ]
   end subroutine collect_check_rccsd_ints

   ! One wrapper per case: a unittest must not depend on which ran before it.
   subroutine w_scf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 1)
   end subroutine w_scf

   subroutine w_mp2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 1)
   end subroutine w_mp2

   subroutine w_probes(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 1)
   end subroutine w_probes

   subroutine w_algebra(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 1)
   end subroutine w_algebra

   subroutine w_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 1)
   end subroutine w_symmetry

   subroutine h_scf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1, 2)
   end subroutine h_scf

   subroutine h_mp2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 2, 2)
   end subroutine h_mp2

   subroutine h_probes(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 3, 2)
   end subroutine h_probes

   subroutine h_algebra(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 4, 2)
   end subroutine h_algebra

   subroutine h_symmetry(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5, 2)
   end subroutine h_symmetry

   subroutine refusals(error)
      type(error_type), allocatable, intent(out) :: error
      call ensure_computed()
      call check(error, refusals_ok, trim(refusals_msg))
   end subroutine refusals

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
      print "(a,l1,a,a)", "  check_rccsd_ints ", passed, " ", text
   end subroutine record

   subroutine compute_molecule(im)
      !! SCF, the blocks, and every number the cases compare
      integer, intent(in) :: im

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_result_t) :: res
      type(trc_rimp2_result_t) :: mp
      type(trc_rimp2_options_t) :: mp_opts
      type(trc_rccsd_ints_t) :: ints, blocked
      type(trc_linalg_t) :: la
      integer, allocatable :: zint(:)
      real(dp), allocatable :: cmo_c(:, :), eps_c(:), dip(:, :, :), rebuilt(:, :)
      integer :: nocc, nfz, nmo, ref_unit, no, nv, nact, ib, k
      integer :: ref_no, ref_nv, ref_naux, ref_nfz
      real(dp) :: e_scf_ref, e_mp2_ref, e_fixed_ref, e_mp2, probe(4), ref(4), d, dscf, dmp2, drimp2
      real(dp) :: worst_block(NBLOCK), b_ovov, b_oovv, blocking
      character(len=:), allocatable :: label
      logical :: ok, supported

      label = trim(MOL_LABEL(im))//": "
      nfz = MOL_NFROZEN(im)
      call load_molecule(im, bas, aux, pl, zint, nocc, supported)
      if (.not. supported) then
         ! Same self-skip convention as check_df_vs_pyscf (SKIP_REGULAR_EXPRESSION in CMake).
         print "(a)", "SKIP check_rccsd_ints: this build's density-fitting kernels stop below the cc-pVDZ-RIFIT shells"
         stop 0
      end if
      nmo = bas%nao_sph
      call run_scf(bas, nocc, res)

      ! ---- the blocks, whole and in auxiliary blocks of 20 -------------------
      call out_build(ints, 0)
      call out_build(blocked, 20)
      no = ints%no
      nv = ints%nv
      nact = no + nv

      ! ---- reference file -----------------------------------------------------
      call ref_require(trim(MOL_REF(im)), ref_unit)
      call ref_read_int(ref_unit, "NO", ref_no)
      call ref_read_int(ref_unit, "NV", ref_nv)
      call ref_read_int(ref_unit, "NAUX", ref_naux)
      call ref_read_int(ref_unit, "NFROZEN", ref_nfz)
      call ref_read_real(ref_unit, "E_SCF", e_scf_ref)
      call ref_read_real(ref_unit, "E_MP2", e_mp2_ref)
      call ref_read_real(ref_unit, "E_FIXED", e_fixed_ref)

      ! ---- E_SCF and E_MP2 ----------------------------------------------------
      dscf = abs(res%energy - e_scf_ref)
      ok = res%converged .and. dscf < TOL_SCF
      call record(im, 1, ok, label//"E_SCF vs pyscf, |dE| = "//to_string(dscf, "es10.2")// &
                  trim(merge("                     ", "  (SCF NOT CONVERGED)", res%converged)))

      e_mp2 = ints%mp2_energy()
      mp_opts%nfrozen = nfz
      allocate (cmo_c(bas%nao, bas%nao), eps_c(bas%nao))
      call rimp2_driver_mo_to_cart(bas, res%cmo(:, :, 1), res%eps(:, 1), cmo_c, eps_c)
      call trc_rimp2_run(bas, aux, pl, nocc, cmo_c, eps_c, mp, opts=mp_opts)
      dmp2 = abs(e_mp2 - e_mp2_ref)
      drimp2 = abs(e_mp2 - mp%e_corr)
      ok = len_trim(mp%message) == 0 .and. dmp2 < TOL_MP2 .and. drimp2 < TOL_RIMP2 .and. &
           ref_no == no .and. ref_nv == nv .and. ref_naux == ints%naux .and. ref_nfz == ints%nfrozen
      call record(im, 2, ok, label//"E_MP2 from ovov: vs pyscf "//to_string(dmp2, "es10.2")// &
                  ", vs trc_rimp2_run "//to_string(drimp2, "es10.2")//"; no/nv/naux/nfrozen "// &
                  to_string(ints%no)//"/"//to_string(ints%nv)//"/"//to_string(ints%naux)//"/"// &
                  to_string(ints%nfrozen)//" vs ref "//to_string(ref_no)//"/"//to_string(ref_nv)//"/"// &
                  to_string(ref_naux)//"/"//to_string(ref_nfz)//" "//trim(mp%message))

      ! ---- the probes ---------------------------------------------------------
      allocate (dip(nact, nact, 4))
      call mo_dipoles(bas, pl, res%cmo(:, nfz + 1:nmo, 1), dip)
      worst_block = 0.0_dp
      do ib = 1, NBLOCK
         call ref_read_vector(ref_unit, "P_"//BLOCK_NAME(ib), 4, ref)
         select case (ib)
         case (1); probe = probe4(ints%oooo, BLOCK_NAME(ib), dip, no)
         case (2); probe = probe4(ints%ooov, BLOCK_NAME(ib), dip, no)
         case (3); probe = probe4(ints%oovv, BLOCK_NAME(ib), dip, no)
         case (4); probe = probe4(ints%ovov, BLOCK_NAME(ib), dip, no)
         case default; probe = probe4(ints%ovvv, BLOCK_NAME(ib), dip, no)
         end select
         do k = 1, 4
            d = abs(probe(k) - ref(k))/max(abs(ref(k)), PROBE_FLOOR)
            worst_block(ib) = max(worst_block(ib), d)
         end do
         print "(a,a,a,a,a,4es11.3)", "  check_rccsd_ints ", label, BLOCK_NAME(ib), " absolute probe diffs", &
            " vs pyscf: ", abs(probe - ref)
      end do
      call ref_close(ref_unit)
      ok = all(worst_block < TOL_PROBE)
      call record(im, 3, ok, label//"block probes vs pyscf, worst relative (floor 0.1) oooo/ooov/oovv/ovov/ovvv = "// &
                  to_string(worst_block(1), "es9.2")//" "//to_string(worst_block(2), "es9.2")//" "// &
                  to_string(worst_block(3), "es9.2")//" "//to_string(worst_block(4), "es9.2")//" "// &
                  to_string(worst_block(5), "es9.2"))

      ! ---- B reproduces the blocks, and blocking changes nothing --------------
      allocate (rebuilt(no*nv, no*nv))
      rebuilt(:, :) = matmul(ints%b_ov, transpose(ints%b_ov))
      b_ovov = maxval(abs(rebuilt - reshape(ints%ovov, [no*nv, no*nv])))
      deallocate (rebuilt)
      allocate (rebuilt(no*no, nv*nv))
      rebuilt(:, :) = matmul(ints%b_oo, transpose(ints%b_vv))
      b_oovv = maxval(abs(rebuilt - reshape(ints%oovv, [no*no, nv*nv])))
      deallocate (rebuilt)
      blocking = max(maxval(abs(ints%ovvv - blocked%ovvv)), maxval(abs(ints%oooo - blocked%oooo)), &
                     maxval(abs(ints%b_vv - blocked%b_vv)))
      ok = b_ovov < TOL_ALGEBRA .and. b_oovv < TOL_ALGEBRA .and. blocking < TOL_ALGEBRA
      call record(im, 4, ok, label//"ovov from b_ov: "//to_string(b_ovov, "es9.2")// &
                  ", oovv from b_oo b_vv: "//to_string(b_oovv, "es9.2")// &
                  ", aux blocks of 20 vs whole: "//to_string(blocking, "es9.2"))

      call symmetry_case(im, ints)
      if (im == 1) call refusal_cases()

      call ints%release()
      call blocked%release()
      call la%release()
      call pl%release()
      call aux%release()
      call bas%release()

   contains

      subroutine out_build(target, aux_block)
         type(trc_rccsd_ints_t), intent(inout) :: target
         integer, intent(in) :: aux_block
         type(trc_rccsd_ints_options_t) :: bopts
         bopts%nalpha = nocc
         bopts%nbeta = nocc
         bopts%nfrozen = nfz
         bopts%aux_block = aux_block
         call target%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, bopts)
         if (target%status%has_error()) error stop "check_rccsd_ints: build: "//target%status%get_message()
      end subroutine out_build

      subroutine refusal_cases()
         !! Each refusal returns its status and message and leaves nothing allocated.
         type(trc_rccsd_ints_t) :: r
         type(trc_rccsd_ints_options_t) :: ro
         character(len=500) :: why
         real(dp), allocatable :: cshort(:, :)

         why = ""
         ro = good_options()
         ro%has_ecp = .true.
         call r%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, ro)
         call note_refusal("an ECP", "ECP", r, why)
         ro = good_options()
         ro%nranks = 2
         call r%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, ro)
         call note_refusal("two ranks", "one rank only", r, why)
         ro = good_options()
         ro%nalpha = nocc + 1
         call r%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, ro)
         call note_refusal("an open shell", "open shell", r, why)
         ro = good_options()
         ro%nfrozen = nocc
         call r%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, ro)
         call note_refusal("every occupied orbital frozen", "no active occupied", r, why)
         ro = good_options()
         ro%nalpha = nmo
         ro%nbeta = nmo
         ro%nfrozen = 0
         call r%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, ro)
         call note_refusal("no virtual orbital", "no virtual", r, why)
         ro = good_options()
         ro%nfrozen = -1
         call r%build(bas, aux, pl, res%cmo(:, :, 1), res%eps(:, 1), la, ro)
         call note_refusal("a negative frozen count", "frozen-core count", r, why)
         allocate (cshort(nmo - 1, nmo - 1))
         cshort = res%cmo(1:nmo - 1, 1:nmo - 1, 1)
         call r%build(bas, aux, pl, cshort, res%eps(1:nmo - 1, 1), la, good_options())
         call note_refusal("orbitals of the wrong size", "sized by the basis", r, why)
         call r%release()
         refusals_ok = len_trim(why) == 0
         refusals_msg = "every refusal returns its message and allocates nothing"
         if (.not. refusals_ok) refusals_msg = "refusals: "//trim(why)
      end subroutine refusal_cases

      function good_options() result(o)
         type(trc_rccsd_ints_options_t) :: o
         o%nalpha = nocc
         o%nbeta = nocc
         o%nfrozen = nfz
      end function good_options

   end subroutine compute_molecule

   subroutine note_refusal(what, keyword, r, why)
      !! Appends to `why` unless `r%status` carries `keyword` and `r` holds nothing
      character(len=*), intent(in) :: what, keyword
      type(trc_rccsd_ints_t), intent(in) :: r
      character(len=*), intent(inout) :: why
      logical :: good
      good = r%status%has_error()
      if (good) good = index(r%status%get_message(), keyword) > 0
      if (good) good = .not. allocated(r%ovov) .and. .not. allocated(r%b_vv) .and. .not. allocated(r%eps_o)
      if (.not. good) why = trim(why)//" ["//what//" not refused with '"//keyword//"']"
   end subroutine note_refusal

   subroutine symmetry_case(im, ints)
      !! oooo(ijkl) = oooo(jikl) = oooo(klij); ovvv(iabc) = ovvv(iacb); ovov(iajb) = ovov(jbia)
      integer, intent(in) :: im
      type(trc_rccsd_ints_t), intent(in) :: ints
      integer :: i, j, k, l, a, b, c, no, nv
      real(dp) :: s_ij, s_pair, s_cb, s_ovov

      no = ints%no
      nv = ints%nv
      s_ij = 0.0_dp
      s_pair = 0.0_dp
      s_cb = 0.0_dp
      s_ovov = 0.0_dp
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  s_ij = max(s_ij, abs(ints%oooo(i, j, k, l) - ints%oooo(j, i, k, l)))
                  s_pair = max(s_pair, abs(ints%oooo(i, j, k, l) - ints%oooo(k, l, i, j)))
               end do
            end do
         end do
      end do
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do i = 1, no
                  s_cb = max(s_cb, abs(ints%ovvv(i, a, b, c) - ints%ovvv(i, a, c, b)))
               end do
            end do
         end do
      end do
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  s_ovov = max(s_ovov, abs(ints%ovov(i, a, j, b) - ints%ovov(j, b, i, a)))
               end do
            end do
         end do
      end do
      call record(im, 5, max(s_ij, s_pair, s_cb, s_ovov) < TOL_ALGEBRA, trim(MOL_LABEL(im))// &
                  ": oooo(ijkl)-oooo(jikl) "//to_string(s_ij, "es9.2")//", oooo(ijkl)-oooo(klij) "// &
                  to_string(s_pair, "es9.2")//", ovvv(iabc)-ovvv(iacb) "//to_string(s_cb, "es9.2")// &
                  ", ovov(iajb)-ovov(jbia) "//to_string(s_ovov, "es9.2"))
   end subroutine symmetry_case

end module check_rccsd_ints_suite

program check_rccsd_ints
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_rccsd_ints_suite, only: collect_check_rccsd_ints
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_check_rccsd_ints, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_rccsd_ints
