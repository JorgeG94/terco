!
! SPHERICAL-basis correlation: RI-MP2 (restricted and unrestricted) and
! density-fitted SCF (`--df`) on a spherical orbital basis, against pyscf's
! own `cart=False` on the identical sets -- the spherical counterpart of
! scf_rimp2.F90/scf_rimp2_uhf.F90/scf_df.F90/scf_dh.F90, which all force
! Cartesian.
!
! FOUR GROUPS:
!
!   A. Restricted RI-MP2: water and formamide, cc-pVDZ + cc-pVDZ-RIFIT
!      (both spherical in BSE), full core correlated -- scf_rimp2.F90's own
!      water/cc-pVDZ case's convention, just spherical.
!
!   B. Unrestricted RI-MP2: OH, cc-pVDZ + cc-pVDZ-RIFIT, frozen core (O),
!      from a FIXED initial density -- scf_rimp2_uhf.F90's own pattern:
!      OH's unpaired electron sits in a rotationally degenerate pair of pi
!      orbitals (a genuine zero-curvature direction for the SCF, and for
!      which particular rotation within that pair ends up occupied), so
!      both sides are seeded from the SAME converged density rather than
!      two that merely happen to share an energy. Unlike that file's own
!      Cartesian case, no AO-normalisation rescale is needed here: a
!      spherical AO's normalisation already matches pyscf's directly, with
!      nothing GAMESS-like to diverge on (scripts/gen_sph.py's own "NOTHING
!      EXTRA NEEDED" header) -- this suite's own S dump still lets
!      generate.py CHECK that empirically rather than assume it.
!
!   C. `--df` HF and B3LYP: water and formamide, cc-pVDZ, def2-universal-
!      jkfit -- scf_df.F90's own cases, spherical.
!
!   D. B2PLYP (a double hybrid): water, cc-pVDZ, four-centre AND `--df`
!      (cc-pVDZ-RIFIT doing double duty as both --df's jkfit and RI-MP2's
!      aux, chained -- scf_dh.F90's own pattern), frozen core (O).
!
! Writes "sph_corr_oh_probe.bin": NAO and S in terco's SPHERICAL AO space
! (trc_1e's own Cartesian S, folded once with sph_cart_to_sph -- the same
! fold trc_scf_driver itself applies to S every geometry) -- geometry- and
! basis-derived only, dumped unconditionally, group B's own "fixed input"
! (scf_rimp2_uhf.F90's term for the same category).
!
! Nine independent cases below (A: 2, B: 1, C: 4, D: 2), each its own
! test-drive unittest; the ONE pyscf reference file is opened once and
! read sequentially, in the SAME A->B->C->D, within-group order the
! original program used -- see check_xc_energy_suite's own `computed`
! comment for why, and why this suite stays `parallel=.false.`.
!
module scf_sph_corr_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, trc_1e
   use trc_sph, only: sph_cart_to_sph
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_rimp2_driver, only: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, rimp2_driver_mo_to_cart, &
                               trc_rimp2_options_t
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3
   implicit none
   private
   public :: collect_scf_sph_corr

   real(dp), parameter :: TOL_SCF = 1.0e-9_dp    !! scf_rimp2.ref's own closed-shell tolerance
   real(dp), parameter :: TOL_DF = 1.0e-8_dp     !! scf_df.ref's own fitted-SCF tolerance
   real(dp), parameter :: TOL_CORR = 1.0e-9_dp   !! rimp2.ref's own tolerance

   integer, parameter :: NCASE = 9
   logical :: computed = .false.
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)
   integer :: next_slot = 0
   !! Which ok_case/msg_case slot run_all's inner ref_report-style calls
   !! write to next, incremented by `record` -- run_all fills every slot
   !! in one pass; each case_NN wrapper only reads its own slot back.

contains

   subroutine collect_scf_sph_corr(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("A water/cc-pvdz(sph) RI-MP2 vs pyscf", case_01), &
                  new_unittest("A formamide/cc-pvdz(sph) RI-MP2 vs pyscf", case_02), &
                  new_unittest("B oh/cc-pvdz(sph) UHF RI-MP2 vs pyscf", case_03), &
                  new_unittest("C water/cc-pvdz(sph) --df HF vs pyscf", case_04), &
                  new_unittest("C water/cc-pvdz(sph) --df B3LYP vs pyscf", case_05), &
                  new_unittest("C formamide/cc-pvdz(sph) --df HF vs pyscf", case_06), &
                  new_unittest("C formamide/cc-pvdz(sph) --df B3LYP vs pyscf", case_07), &
                  new_unittest("D b2plyp/cc-pvdz(sph) 4c vs pyscf", case_08), &
                  new_unittest("D b2plyp/cc-pvdz(sph) df vs pyscf", case_09) &
                  ]
   end subroutine collect_scf_sph_corr

   ! One wrapper per case, each passing its own index: a unittest must not
   ! depend on which unittests ran before it (a name-selected run, or a
   ! reordered suite, would otherwise check the wrong case).
   subroutine case_01(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 1)
   end subroutine case_01

   subroutine case_02(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 2)
   end subroutine case_02

   subroutine case_03(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 3)
   end subroutine case_03

   subroutine case_04(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 4)
   end subroutine case_04

   subroutine case_05(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 5)
   end subroutine case_05

   subroutine case_06(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 6)
   end subroutine case_06

   subroutine case_07(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 7)
   end subroutine case_07

   subroutine case_08(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 8)
   end subroutine case_08

   subroutine case_09(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 9)
   end subroutine case_09

   subroutine test_case(error, ic)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic
      call ensure_computed()
      call check(error, ok_case(ic), trim(msg_case(ic)))
   end subroutine test_case

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   ! Appends one (ok, message) pair to the next ok_case/msg_case slot.
   subroutine record(ok, msg)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: msg
      next_slot = next_slot + 1
      ok_case(next_slot) = ok
      msg_case(next_slot) = msg
   end subroutine record

   subroutine run_all()
      integer :: dev, refu
      logical :: has_ref

      dev = trc_bind_device(0)
      ok_case = .false.
      msg_case = "pyscf reference missing"
      next_slot = 0
      call ref_open("sph_corr.ref", refu, has_ref)

      call part_a(refu, has_ref)
      call part_b(refu, has_ref)
      call part_c(refu, has_ref)
      call part_d(refu, has_ref)

      if (has_ref) call ref_close(refu)
      if (.not. has_ref) error stop "FAIL: missing reference pyscf_refs/sph_corr.ref -- run " &
         //"test/pyscf_refs/generate.py sph_corr (see test/pyscf_refs/README.md)"
   end subroutine run_all

   !
   ! A. Restricted RI-MP2, spherical: water and formamide, cc-pVDZ +
   ! cc-pVDZ-RIFIT, full core.
   !
   subroutine part_a(refu, has_ref)
      integer, intent(in) :: refu
      logical, intent(in) :: has_ref

      integer, parameter :: NC = 2
      character(len=16), parameter :: gmol(NC) = [character(len=16) :: "water", "formamide"]

      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_rimp2_result_t) :: mp
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), cmo_c(:, :), eps_c(:)
      integer :: ic, natm, nelec, nocc
      character(len=32) :: tag, name_ref, name_this
      real(dp) :: e_ref, e_os_ref, e_ss_ref, de, worst
      logical :: gates_ok
      character(len=:), allocatable :: gate_msg

      do ic = 1, NC
         name_this = trim(gmol(ic))//"_ccpvdz_sph"
         gates_ok = .true.
         gate_msg = ""
         call read_xyz(trim(gmol(ic))//".xyz", natm, zint, at_r)
         nelec = sum(zint)
         nocc = nelec/2

         call trc_basis_from_json("basis_sets/cc-pvdz.json", natm, zint, at_r, bas, err)
         if (err%has_error()) error stop "scf_sph_corr/A: "//err%get_message()
         if (.not. bas%spherical) then
            print "(a)", "scf_sph_corr/A: cc-pVDZ did not come back spherical"
            gates_ok = .false.
            gate_msg = gate_msg//"cc-pVDZ did not come back spherical; "
         end if
         call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err)
         if (err%has_error()) error stop "scf_sph_corr/A: "//err%get_message()
         call bas%to_device()
         call aux%to_device()
         call pl%build(bas, 1.0e-12_dp)
         call pl%to_device()

         opts = trc_scf_options_t()
         opts%conv_energy = 1.0e-12_dp
         opts%conv_density = 1.0e-9_dp
         call trc_scf_run(bas, nocc, nocc, opts, res)
         if (.not. res%converged) print "(a,a)", "  NOT CONVERGED: ", trim(res%message)

         allocate (cmo_c(bas%nao, bas%nao), eps_c(bas%nao))
         call rimp2_driver_mo_to_cart(bas, res%cmo(:, :, 1), res%eps(:, 1), cmo_c, eps_c)
         call trc_rimp2_run(bas, aux, pl, nocc, cmo_c, eps_c, mp)
         if (len_trim(mp%message) > 0) then
            print "(a)", "scf_sph_corr/A: "//trim(mp%message)
            gates_ok = .false.
            gate_msg = gate_msg//"RI-MP2: "//trim(mp%message)//"; "
         end if
         deallocate (cmo_c, eps_c)

         print "(a,a10,a,i0,a,f20.12,a,f16.12,a,f16.12)", "  A ", trim(gmol(ic)), "  nao_sph ", bas%nao_sph, &
            "  E_RHF ", res%energy, "  E_os ", mp%e_os, "  E_ss ", mp%e_ss

         if (has_ref) then
            tag = to_string(ic - 1)
            call ref_read_str(refu, "RMP2_NAME_"//trim(tag), name_ref)
            if (trim(name_ref) /= trim(name_this)) error stop "scf_sph_corr/A: case order mismatch (" &
               //trim(name_ref)//" vs "//trim(name_this)//")"
            call ref_read_real(refu, "RMP2_ESCF_"//trim(tag), e_ref)
            call ref_read_real(refu, "RMP2_EOS_"//trim(tag), e_os_ref)
            call ref_read_real(refu, "RMP2_ESS_"//trim(tag), e_ss_ref)
            de = abs(res%energy - e_ref)
            worst = max(de, max(abs(mp%e_os - e_os_ref), abs(mp%e_ss - e_ss_ref)))
            call record(gates_ok .and. res%converged .and. de < TOL_SCF .and. abs(mp%e_os - e_os_ref) < TOL_CORR &
                        .and. abs(mp%e_ss - e_ss_ref) < TOL_CORR, &
                        gate_msg//trim(name_this)//": E_RHF/E_os/E_ss vs pyscf, worst = "//to_string(worst, "es10.2"))
         else
            call record(.false., gate_msg//trim(name_this)//": pyscf reference missing")
         end if

         call pl%release(); call aux%release(); call bas%release()
         deallocate (zint, at_r)
      end do
   end subroutine part_a

   !
   ! B. Unrestricted RI-MP2, spherical: OH, cc-pVDZ + cc-pVDZ-RIFIT,
   ! frozen core, from a fixed initial density.
   !
   subroutine part_b(refu, has_ref)
      integer, intent(in) :: refu
      logical, intent(in) :: has_ref

      integer, parameter :: NA = 5, NB = 4, NFZ = 1
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_rimp2_result_t) :: mp
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp), allocatable :: s_c(:, :), t_c(:, :), v_c(:, :), s_sph(:, :)
      real(dp), allocatable :: cmo_ab(:, :, :), eps_ab(:, :)
      integer :: natm, probe_unit
      real(dp) :: e_scf_ref, e_os_ref, e_ss_ref, worst
      character(len=:), allocatable :: gate_msg

      call read_xyz("oh.xyz", natm, zint, at_r)
      call trc_basis_from_json("basis_sets/cc-pvdz.json", natm, zint, at_r, bas, err)
      if (err%has_error()) error stop "scf_sph_corr/B: "//err%get_message()
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err)
      if (err%has_error()) error stop "scf_sph_corr/B: "//err%get_message()
      call bas%to_device()
      call aux%to_device()
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()

      ! S alone, folded into the SPHERICAL AO space (trc_1e stays
      ! Cartesian; sph_cart_to_sph is the same fold trc_scf_driver itself
      ! applies to S every geometry) -- never from an SCF, dumped
      ! unconditionally (scf_rimp2_uhf.F90's own pattern).
      allocate (s_c(bas%nao, bas%nao), t_c(bas%nao, bas%nao), v_c(bas%nao, bas%nao))
      call trc_1e(bas, pl, s_c, t_c, v_c)
      allocate (s_sph(bas%nao_sph, bas%nao_sph))
      call sph_cart_to_sph(bas, s_c, s_sph)
      open (newunit=probe_unit, file="sph_corr_oh_probe.bin", access="stream", form="unformatted", &
            status="replace", action="write")
      write (probe_unit) int(bas%nao_sph, int64)
      write (probe_unit) s_sph
      close (probe_unit)
      deallocate (s_c, t_c, v_c, s_sph)

      opts = trc_scf_options_t()
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-9_dp
      opts%conv_diis = 1.0e-11_dp   ! the commutator norm gates, not the density; see scf_rimp2_uhf.F90

      if (has_ref) then
         block
            integer :: nspin_chk, nao_chk
            real(dp), allocatable :: dm0(:, :, :)
            call ref_read_int(refu, "OH_NSPIN", nspin_chk)
            call ref_read_int(refu, "OH_NAO", nao_chk)
            if (nspin_chk /= 2 .or. nao_chk /= bas%nao_sph) &
               error stop "scf_sph_corr/B: pyscf reference shape does not match this build"
            allocate (dm0(bas%nao_sph, bas%nao_sph, 2))
            call ref_read_matrix3(refu, "OH_DM0", bas%nao_sph, bas%nao_sph, 2, dm0)
            call trc_scf_run(bas, NA, NB, opts, res, dguess=dm0)
         end block
      else
         call trc_scf_run(bas, NA, NB, opts, res)
      end if
      if (.not. res%converged) print "(a,a)", "  NOT CONVERGED (OH): ", trim(res%message)

      allocate (cmo_ab(bas%nao, bas%nao, 2), eps_ab(bas%nao, 2))
      call rimp2_driver_mo_to_cart(bas, res%cmo(:, :, 1), res%eps(:, 1), cmo_ab(:, :, 1), eps_ab(:, 1))
      call rimp2_driver_mo_to_cart(bas, res%cmo(:, :, 2), res%eps(:, 2), cmo_ab(:, :, 2), eps_ab(:, 2))
      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%nfrozen = NFZ
         call trc_urimp2_run(bas, aux, pl, NA, NB, cmo_ab, eps_ab, mp, opts=rimp2_opts)
      end block
      gate_msg = ""
      if (len_trim(mp%message) > 0) then
         print "(a)", "scf_sph_corr/B: "//trim(mp%message)
         gate_msg = "RI-MP2: "//trim(mp%message)//"; "
      end if
      deallocate (cmo_ab, eps_ab)

      print "(a,i0,a,f20.12,a,f16.12,a,f16.12)", "  B oh  nao_sph ", bas%nao_sph, "  E_UHF ", res%energy, &
         "  E_os ", mp%e_os, "  E_ss ", mp%e_ss

      if (has_ref) then
         call ref_read_real(refu, "OH_ESCF", e_scf_ref)
         call ref_read_real(refu, "OH_EOS", e_os_ref)
         call ref_read_real(refu, "OH_ESS", e_ss_ref)
         worst = max(abs(res%energy - e_scf_ref), max(abs(mp%e_os - e_os_ref), abs(mp%e_ss - e_ss_ref)))
         call record(res%converged .and. abs(res%energy - e_scf_ref) < TOL_SCF .and. &
                     abs(mp%e_os - e_os_ref) < TOL_CORR .and. abs(mp%e_ss - e_ss_ref) < TOL_CORR .and. &
                     len_trim(mp%message) == 0, &
                     gate_msg//"oh/cc-pvdz(sph): E_UHF/E_os/E_ss vs pyscf, worst = "//to_string(worst, "es10.2"))
      else
         call record(.false., gate_msg//"oh/cc-pvdz(sph): pyscf reference missing")
      end if

      call pl%release(); call aux%release(); call bas%release()
   end subroutine part_b

   !
   ! C. `--df` HF and B3LYP, spherical: water and formamide, cc-pVDZ,
   ! def2-universal-jkfit.
   !
   subroutine part_c(refu, has_ref)
      integer, intent(in) :: refu
      logical, intent(in) :: has_ref

      integer, parameter :: NGROUP = 2, NMETHOD = 2
      character(len=16), parameter :: gmol(NGROUP) = [character(len=16) :: "water", "formamide"]
      character(len=8), parameter :: method(NMETHOD) = [character(len=8) :: "", "b3lyp"]

      type(trc_basis_t) :: bas, auxb
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t), target :: dfjk
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: ig, im, ic, natm, nalpha, nbeta, nelec
      character(len=32) :: tag, name_ref, name_this
      real(dp) :: e_ref, de

      do ig = 1, NGROUP
         call read_xyz(trim(gmol(ig))//".xyz", natm, zint, at_r)
         nelec = sum(zint)
         nalpha = (nelec + mod(nelec, 2))/2
         nbeta = nelec - nalpha

         call trc_basis_from_json("basis_sets/cc-pvdz.json", natm, zint, at_r, bas, err)
         if (err%has_error()) error stop "scf_sph_corr/C: "//err%get_message()
         call bas%to_device()
         call trc_basis_from_json("basis_sets/def2-universal-jkfit.json", natm, zint, at_r, auxb, err)
         if (err%has_error()) error stop "scf_sph_corr/C: "//err%get_message()
         call auxb%to_device()

         call pl%build(bas, 1.0e-12_dp)
         call pl%to_device()
         call la%init(bas%nao)
         call dfjk%setup(bas, auxb, pl, la, err, nocc_hint=max(nalpha, 1))
         if (err%has_error()) error stop "scf_sph_corr/C: "//err%get_message()
         call pl%release()

         do im = 1, NMETHOD
            opts = trc_scf_options_t()
            opts%df = .true.
            opts%functional = method(im)
            opts%grid_level = 3
            opts%conv_energy = 1.0e-11_dp
            opts%conv_diis = 1.0e-8_dp
            call trc_scf_run(bas, nalpha, nbeta, opts, res, la=la, scf_dfjk=dfjk)
            if (.not. res%converged) print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
            print "(a,a10,1x,a6,a,i3,a,f22.12)", "  C ", trim(gmol(ig)), method(im), " iter ", res%iterations, &
               "  E = ", res%energy

            ic = (ig - 1)*NMETHOD + im
            if (len_trim(method(im)) == 0) then
               name_this = trim(gmol(ig))//"_ccpvdz_sph_hf"
            else
               name_this = trim(gmol(ig))//"_ccpvdz_sph_"//trim(method(im))
            end if
            if (has_ref) then
               tag = to_string(ic - 1)
               call ref_read_str(refu, "DF_NAME_"//trim(tag), name_ref)
               if (trim(name_ref) /= trim(name_this)) error stop "scf_sph_corr/C: case order mismatch (" &
                  //trim(name_ref)//" vs "//trim(name_this)//")"
               call ref_read_real(refu, "DF_ESCF_"//trim(tag), e_ref)
               de = abs(res%energy - e_ref)
               call record(res%converged .and. de < TOL_DF, &
                           trim(name_this)//": E vs pyscf, |dE| = "//to_string(de, "es10.2"))
            else
               call record(.false., trim(name_this)//": pyscf reference missing")
            end if
         end do

         call dfjk%release(); call auxb%release(); call bas%release(); call la%release()
         deallocate (zint, at_r)
      end do
   end subroutine part_c

   !
   ! D. B2PLYP, spherical: water, cc-pVDZ, four-centre AND `--df` (chained
   ! RI-MP2 off --df's own fit, cc-pVDZ-RIFIT as both jkfit and aux),
   ! frozen core.
   !
   subroutine part_d(refu, has_ref)
      integer, intent(in) :: refu
      logical, intent(in) :: has_ref

      integer, parameter :: NFZ = 1
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t), target :: dfjk
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res_4c, res_df
      type(trc_rimp2_result_t) :: mp_4c, mp_df
      type(error_t) :: err
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), cmo_c(:, :), eps_c(:)
      integer :: natm, nalpha, nbeta, nelec
      real(dp) :: e_ref, e_os_ref, e_ss_ref, worst

      call read_xyz("water.xyz", natm, zint, at_r)
      nelec = sum(zint)
      nalpha = nelec/2
      nbeta = nelec - nalpha

      call trc_basis_from_json("basis_sets/cc-pvdz.json", natm, zint, at_r, bas, err)
      if (err%has_error()) error stop "scf_sph_corr/D: "//err%get_message()
      call bas%to_device()
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err)
      if (err%has_error()) error stop "scf_sph_corr/D: "//err%get_message()
      call aux%to_device()
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      call la%init(bas%nao)
      call dfjk%setup(bas, aux, pl, la, err, nocc_hint=max(nalpha, 1))
      if (err%has_error()) error stop "scf_sph_corr/D: "//err%get_message()

      ! --- four-centre ------------------------------------------------------
      opts = trc_scf_options_t()
      opts%functional = "b2plyp"
      opts%grid_level = 3
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-8_dp
      opts%df = .false.
      call trc_scf_run(bas, nalpha, nbeta, opts, res_4c, la=la)
      allocate (cmo_c(bas%nao, bas%nao), eps_c(bas%nao))
      call rimp2_driver_mo_to_cart(bas, res_4c%cmo(:, :, 1), res_4c%eps(:, 1), cmo_c, eps_c)
      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%nfrozen = NFZ
         call trc_rimp2_run(bas, aux, pl, nalpha, cmo_c, eps_c, mp_4c, opts=rimp2_opts, la=la)
      end block
      deallocate (cmo_c, eps_c)
      if (.not. res_4c%converged) print "(a,a)", "  NOT CONVERGED (4c): ", trim(res_4c%message)
      print "(a,i3,a,f20.12,a,f16.12,a,f16.12)", "  D b2plyp 4c  iter ", res_4c%iterations, "  E_SCF ", &
         res_4c%energy, "  E_os ", mp_4c%e_os, "  E_ss ", mp_4c%e_ss

      ! --- --df, chained RI-MP2 off its own fit -----------------------------
      opts%df = .true.
      call trc_scf_run(bas, nalpha, nbeta, opts, res_df, la=la, scf_dfjk=dfjk)
      allocate (cmo_c(bas%nao, bas%nao), eps_c(bas%nao))
      call rimp2_driver_mo_to_cart(bas, res_df%cmo(:, :, 1), res_df%eps(:, 1), cmo_c, eps_c)
      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%nfrozen = NFZ
         rimp2_opts%packed_b => dfjk%packed_b()
         call trc_rimp2_run(bas, aux, pl, nalpha, cmo_c, eps_c, mp_df, opts=rimp2_opts, la=la)
      end block
      deallocate (cmo_c, eps_c)
      if (.not. res_df%converged) print "(a,a)", "  NOT CONVERGED (df): ", trim(res_df%message)
      print "(a,i3,a,f20.12,a,f16.12,a,f16.12)", "  D b2plyp df  iter ", res_df%iterations, "  E_SCF ", &
         res_df%energy, "  E_os ", mp_df%e_os, "  E_ss ", mp_df%e_ss

      if (has_ref) then
         call ref_read_real(refu, "B2PLYP_ESCF_4C", e_ref)
         call ref_read_real(refu, "B2PLYP_EOS_4C", e_os_ref)
         call ref_read_real(refu, "B2PLYP_ESS_4C", e_ss_ref)
         worst = max(abs(res_4c%energy - e_ref), max(abs(mp_4c%e_os - e_os_ref), abs(mp_4c%e_ss - e_ss_ref)))
         call record(res_4c%converged .and. abs(res_4c%energy - e_ref) < TOL_SCF .and. &
                     abs(mp_4c%e_os - e_os_ref) < TOL_CORR .and. abs(mp_4c%e_ss - e_ss_ref) < TOL_CORR, &
                     "b2plyp/cc-pvdz(sph) 4c: E_SCF/E_os/E_ss vs pyscf, worst = "//to_string(worst, "es10.2"))

         call ref_read_real(refu, "B2PLYP_ESCF_DF", e_ref)
         call ref_read_real(refu, "B2PLYP_EOS_DF", e_os_ref)
         call ref_read_real(refu, "B2PLYP_ESS_DF", e_ss_ref)
         worst = max(abs(res_df%energy - e_ref), max(abs(mp_df%e_os - e_os_ref), abs(mp_df%e_ss - e_ss_ref)))
         call record(res_df%converged .and. abs(res_df%energy - e_ref) < TOL_DF .and. &
                     abs(mp_df%e_os - e_os_ref) < TOL_CORR .and. abs(mp_df%e_ss - e_ss_ref) < TOL_CORR, &
                     "b2plyp/cc-pvdz(sph) df: E_SCF/E_os/E_ss vs pyscf, worst = "//to_string(worst, "es10.2"))
      else
         call record(.false., "b2plyp/cc-pvdz(sph) 4c: pyscf reference missing")
         call record(.false., "b2plyp/cc-pvdz(sph) df: pyscf reference missing")
      end if

      call dfjk%release(); call pl%release(); call aux%release(); call bas%release(); call la%release()
   end subroutine part_d

end module scf_sph_corr_suite

program scf_sph_corr
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_sph_corr_suite, only: collect_scf_sph_corr
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_sph_corr, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_sph_corr
