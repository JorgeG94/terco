!
! Double hybrids (trc_xc_functional.F90's b2plyp/b2gp-plyp/mpw2plyp): the SCF
! through trc_scf_run (four-centre AND --df, trc_scf_df.F90), then RI-MP2
! on the converged KS orbitals through trc_rimp2_run -- chained off --df's
! own fit when the aux basis given to RI-MP2 is the SAME one --df was
! fitted on (cc-pVDZ-RIFIT doing double duty as both --df's jkfit and
! RI-MP2's aux, the one case trc_rimp2_run's own `packed_b` argument exists
! for) -- against pyscf's KS with the SAME hybrid composition
! (`0.53*HF + 0.47*B88, 0.73*LYP` for b2plyp, etc. -- test/pyscf_refs/
! generate.py's own gen_dh) plus DF-MP2/DF-UMP2 on the SAME aux basis and
! frozen core: the committed reference (dh.ref) carries E_SCF, the RAW
! (unscaled) E_os/E_ss, AND pyscf's own weighted total (mf.e_tot +
! c_os*e_os + c_ss*e_ss) as E_DH -- this suite's own c_os*E_os +
! c_ss*E_ss weighting (from LITERAL coefficients, PT2_COS/PT2_CSS below,
! not read back from trc_xc_functional's own table) is checked against
! that committed E_DH directly, never left as an unchecked header claim.
!
! 24 closed-shell cases: water and formamide, cc-pVDZ and 6-31G*,
! b2plyp, b2gp-plyp and mpw2plyp, four-centre and --df -- plus one open-shell case,
! OH/b2plyp UKS/6-31G*, from a FIXED initial density (rks.ref's and
! rimp2_uhf.ref's own pattern: pyscf's own converged density, rescaled
! into this build's AO normalisation via this suite's own dumped
! overlap matrix, so both sides settle into the SAME orientation rather
! than two that merely share an energy -- see scf_rimp2_uhf.F90's header
! for the full story and the measured consequence of skipping it).
!
! Frozen core: core_orbitals(Z) summed per molecule's atoms (trc_capi.F90's
! own RI-MP2 default rule) -- water 1 (O), formamide 3 (C, N, O), OH 1 (O).
!
! Writes "scf_dh_probe.bin": per closed-shell case, nao/natm/nalpha/nbeta
! and the molecule/basis/functional/kind names (a sanity check only --
! generate.py's gen_dh rebuilds the identical mol itself, straight off the
! same xyz/JSON files, scf_df.F90's own pattern); for OH, "scf_dh_oh_probe.bin"
! carries nao and the overlap matrix S, unconditionally (geometry- and
! basis-derived only, never from an SCF -- scf_rimp2_uhf.F90's own "fixed
! input" category).
!
! Twenty-five independent cases below (24 closed-shell + OH), each its
! own test-drive unittest; the shared per-group setup, the probe writes
! and the sequential pyscf-reference read run once -- see
! check_xc_energy_suite's own `computed` comment for why, and why this
! suite stays `parallel=.false.`.
!
module scf_dh_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, trc_1e
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_rimp2_driver, only: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, trc_rimp2_options_t
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3
   implicit none
   private
   public :: collect_scf_dh

   integer, parameter :: NGROUP = 4, NFUNC = 3, NCASE = NGROUP*NFUNC*2 + 1
   character(len=16), parameter :: GMOL(NGROUP) = [character(len=16) :: &
                                                   "water", "water", "formamide", "formamide"]
   character(len=16), parameter :: GBASIS(NGROUP) = [character(len=16) :: &
                                                      "cc-pvdz", "6-31gs", "cc-pvdz", "6-31gs"]
   integer, parameter :: GFROZEN(NGROUP) = [1, 1, 3, 3]   !! core_orbitals(Z), summed: water O, formamide C+N+O
   character(len=8), parameter :: FUNC(NFUNC) = [character(len=8) :: "b2plyp", "b2gpplyp", "mpw2plyp"]
   ! LITERAL coefficients, one pair per FUNC entry above, never read back
   ! from trc_xc_functional's own table: trc_xc_functional.F90's own
   ! comments at each case ("b2plyp")/("b2gp-plyp", "b2gpplyp")/
   ! ("mpw2plyp", "mpw2-plyp") give these same numbers.
   real(dp), parameter :: PT2_COS(NFUNC) = [0.27_dp, 0.36_dp, 0.25_dp]
   real(dp), parameter :: PT2_CSS(NFUNC) = [0.27_dp, 0.36_dp, 0.25_dp]
   character(len=2), parameter :: KIND(2) = [character(len=2) :: "4c", "df"]
   real(dp), parameter :: E_TOL_4C = 1.0e-9_dp    !! scf_rimp2.ref's own closed-shell tolerance
   real(dp), parameter :: E_TOL_DF = 1.0e-8_dp    !! scf_df.ref's own fitted-SCF tolerance
   real(dp), parameter :: CORR_TOL = 1.0e-9_dp    !! rimp2.ref's own tolerance, both SCF kinds

   logical :: computed = .false.
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)
   character(len=32) :: name_case(NCASE)

contains

   subroutine collect_scf_dh(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest(case_name(1), case_01), &
                  new_unittest(case_name(2), case_02), &
                  new_unittest(case_name(3), case_03), &
                  new_unittest(case_name(4), case_04), &
                  new_unittest(case_name(5), case_05), &
                  new_unittest(case_name(6), case_06), &
                  new_unittest(case_name(7), case_07), &
                  new_unittest(case_name(8), case_08), &
                  new_unittest(case_name(9), case_09), &
                  new_unittest(case_name(10), case_10), &
                  new_unittest(case_name(11), case_11), &
                  new_unittest(case_name(12), case_12), &
                  new_unittest(case_name(13), case_13), &
                  new_unittest(case_name(14), case_14), &
                  new_unittest(case_name(15), case_15), &
                  new_unittest(case_name(16), case_16), &
                  new_unittest(case_name(17), case_17), &
                  new_unittest(case_name(18), case_18), &
                  new_unittest(case_name(19), case_19), &
                  new_unittest(case_name(20), case_20), &
                  new_unittest(case_name(21), case_21), &
                  new_unittest(case_name(22), case_22), &
                  new_unittest(case_name(23), case_23), &
                  new_unittest(case_name(24), case_24), &
                  new_unittest(case_name(25), case_25) &
                ]
   end subroutine collect_scf_dh

   function case_name(ic) result(name)
      integer, intent(in) :: ic
      character(len=:), allocatable :: name
      integer :: ig, im, ik
      if (ic == NCASE) then
         name = "oh/b2plyp UKS/6-31G* (fixed density)"
      else
         ! ic = ((ig - 1)*NFUNC + (im - 1))*2 + ik, run_all's own order.
         ig = (ic - 1)/(2*NFUNC) + 1
         im = mod((ic - 1)/2, NFUNC) + 1
         ik = mod(ic - 1, 2) + 1
         name = trim(GMOL(ig))//"/"//trim(GBASIS(ig))//"/"//trim(FUNC(im))//"/"//trim(KIND(ik))
      end if
   end function case_name

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

   subroutine case_10(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 10)
   end subroutine case_10

   subroutine case_11(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 11)
   end subroutine case_11

   subroutine case_12(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 12)
   end subroutine case_12

   subroutine case_13(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 13)
   end subroutine case_13

   subroutine case_14(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 14)
   end subroutine case_14

   subroutine case_15(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 15)
   end subroutine case_15

   subroutine case_16(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 16)
   end subroutine case_16

   subroutine case_17(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 17)
   end subroutine case_17

   subroutine case_18(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 18)
   end subroutine case_18

   subroutine case_19(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 19)
   end subroutine case_19

   subroutine case_20(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 20)
   end subroutine case_20

   subroutine case_21(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 21)
   end subroutine case_21

   subroutine case_22(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 22)
   end subroutine case_22

   subroutine case_23(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 23)
   end subroutine case_23

   subroutine case_24(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 24)
   end subroutine case_24

   subroutine case_25(error)
      type(error_type), allocatable, intent(out) :: error
      call test_case(error, 25)
   end subroutine case_25

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

   subroutine run_all()
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t), target :: dfjk
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_rimp2_result_t) :: mp
      type(error_t) :: err
      logical :: has_ref
      integer :: ig, im, ik, unit, natm, nalpha, nbeta, nelec, ic, refu, ncase_ref
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      character(len=32) :: tag, name_ref, name_this
      real(dp) :: e_ref, de, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
      integer :: dev

      dev = trc_bind_device(0)
      ok_case = .false.
      msg_case = "pyscf reference missing"

      call ref_open("dh.ref", refu, has_ref)
      if (has_ref) then
         call ref_read_int(refu, "NCASE", ncase_ref)
         if (ncase_ref /= NGROUP*NFUNC*2) error stop "scf_dh: pyscf reference case count does not match this build"
      end if

      open (newunit=unit, file="scf_dh_probe.bin", access="stream", form="unformatted", status="replace", &
            action="write")
      write (unit) int(NGROUP*NFUNC*2, int64)

      do ig = 1, NGROUP
         call read_xyz(trim(GMOL(ig))//".xyz", natm, zint, at_r)
         nelec = sum(zint)
         nalpha = nelec/2
         nbeta = nelec - nalpha

         call trc_basis_from_json("basis_sets/"//trim(GBASIS(ig))//".json", natm, zint, at_r, bas, err, &
                                  force_cart=.true.)
         if (err%has_error()) error stop "scf_dh: "//err%get_message()
         call bas%to_device()
         call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err, force_cart=.true.)
         if (err%has_error()) error stop "scf_dh: "//err%get_message()
         call aux%to_device()

         call pl%build(bas, 1.0e-12_dp)
         call pl%to_device()
         call la%init(bas%nao)
         ! --df's own fit, cc-pVDZ-RIFIT doing double duty as the jkfit too:
         ! set up once per group, reused by both functionals' --df case and
         ! by RI-MP2's chained path (trc_rimp2_run's `packed_b` argument).
         call dfjk%setup(bas, aux, pl, la, err, nocc_hint=max(nalpha, 1))
         if (err%has_error()) error stop "scf_dh: "//err%get_message()

         print "(a)", ""
         print "(a,a,a,a,a,i0,a)", "scf_dh: ", trim(GMOL(ig)), "/", trim(GBASIS(ig)), ", ", bas%nao, " AOs"

         do im = 1, NFUNC
            do ik = 1, 2
               opts = trc_scf_options_t()
               opts%functional = trim(FUNC(im))
               opts%grid_level = 3
               opts%conv_energy = 1.0e-11_dp
               opts%conv_diis = 1.0e-8_dp
               if (ik == 1) then
                  opts%df = .false.
                  call trc_scf_run(bas, nalpha, nbeta, opts, res, la=la)
                  block
                     type(trc_rimp2_options_t) :: rimp2_opts
                     rimp2_opts%nfrozen = GFROZEN(ig)
                     call trc_rimp2_run(bas, aux, pl, nalpha, res%cmo(:, :, 1), res%eps(:, 1), mp, &
                                        opts=rimp2_opts, la=la)
                  end block
               else
                  opts%df = .true.
                  call trc_scf_run(bas, nalpha, nbeta, opts, res, la=la, scf_dfjk=dfjk)
                  block
                     type(trc_rimp2_options_t) :: rimp2_opts
                     rimp2_opts%nfrozen = GFROZEN(ig)
                     rimp2_opts%packed_b => dfjk%packed_b()
                     call trc_rimp2_run(bas, aux, pl, nalpha, res%cmo(:, :, 1), res%eps(:, 1), mp, &
                                        opts=rimp2_opts, la=la)
                  end block
               end if
               if (.not. res%converged) print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
               if (len_trim(mp%message) > 0) print "(a)", "scf_dh: "//trim(mp%message)
               print "(a,a8,1x,a2,a,i3,a,f20.12,a,f16.12,a,f16.12)", "  ", FUNC(im), KIND(ik), " iter ", &
                  res%iterations, "  E_SCF ", res%energy, "  E_os ", mp%e_os, "  E_ss ", mp%e_ss

               ic = ((ig - 1)*NFUNC + (im - 1))*2 + ik
               name_this = trim(GMOL(ig))//"_"//trim(GBASIS(ig))//"_"//trim(FUNC(im))//"_"//trim(KIND(ik))
               name_case(ic) = name_this
               if (has_ref) then
                  write (tag, "(i0)") ic - 1
                  call ref_read_str(refu, "NAME_"//trim(tag), name_ref)
                  call ref_read_real(refu, "E_SCF_"//trim(tag), e_ref)
                  if (trim(name_ref) /= trim(name_this)) error stop "scf_dh: pyscf reference case order does " &
                     //"not match this build ("//trim(name_ref)//" vs "//trim(name_this)//")"
                  de = abs(res%energy - e_ref)
                  block
                     real(dp) :: e_os_ref, e_ss_ref, e_dh_ref, e_dh_this, worst
                     logical :: ok1, ok2, ok3, ok4
                     call ref_read_real(refu, "E_OS_"//trim(tag), e_os_ref)
                     call ref_read_real(refu, "E_SS_"//trim(tag), e_ss_ref)
                     ok1 = res%converged .and. de < merge(E_TOL_4C, E_TOL_DF, ik == 1)
                     ok2 = abs(mp%e_os - e_os_ref) < CORR_TOL
                     ok3 = abs(mp%e_ss - e_ss_ref) < CORR_TOL
                     ! The weighted double-hybrid total, against pyscf's OWN
                     ! mf.e_tot + c_os*e_os + c_ss*e_ss (generate.py's gen_dh),
                     ! not merely this suite's own E_SCF/E_os/E_ss checked
                     ! separately above -- PT2_COS/PT2_CSS are the SAME
                     ! literals context_scf.F90's double-hybrid cases assert.
                     call ref_read_real(refu, "E_DH_"//trim(tag), e_dh_ref)
                     e_dh_this = res%energy + PT2_COS(im)*mp%e_os + PT2_CSS(im)*mp%e_ss
                     ok4 = abs(e_dh_this - e_dh_ref) < CORR_TOL
                     worst = max(de, max(abs(mp%e_os - e_os_ref), max(abs(mp%e_ss - e_ss_ref), &
                                 abs(e_dh_this - e_dh_ref))))
                     ! A refused RI-MP2 (mp%message set) fails the case on its
                     ! own, not only through the e_os = 0 it leaves behind.
                     ok_case(ic) = ok1 .and. ok2 .and. ok3 .and. ok4 .and. len_trim(mp%message) == 0
                     write (msg_case(ic), "(a,es10.2)") trim(name_this)// &
                        ": E_SCF/E_os/E_ss/E_dh vs pyscf, worst = ", worst
                     if (len_trim(mp%message) > 0) msg_case(ic) = "RI-MP2: "//trim(mp%message)//"; "// &
                                                                  trim(msg_case(ic))
                  end block
                  ! Provenance only -- scf_rimp2.F90's own comment on the same tags.
                  call ref_read_real(refu, "CONV_TOL_"//trim(tag), pyscf_conv_tol)
                  call ref_read_real(refu, "CONV_TOL_GRAD_"//trim(tag), pyscf_conv_tol_grad)
                  call ref_read_real(refu, "GRAD_NORM_"//trim(tag), pyscf_grad_norm)
               end if

               write (unit) int(bas%nao, int64), int(natm, int64), int(nalpha, int64), int(nbeta, int64)
               write (unit) GMOL(ig)
               write (unit) GBASIS(ig)
               write (unit) FUNC(im)
               write (unit) KIND(ik)
               write (unit) res%energy, mp%e_os, mp%e_ss
            end do
         end do

         call dfjk%release()
         call pl%release()
         call aux%release()
         call bas%release()
         call la%release()
         deallocate (zint, at_r)
      end do
      close (unit)

      ! --- the open-shell case: OH/b2plyp UKS/6-31G*, from a FIXED density ---
      block
         type(trc_basis_t) :: wbas, waux
         type(trc_pairlist_t) :: wpl
         type(trc_scf_options_t) :: wopts
         type(trc_scf_result_t) :: wres
         type(trc_rimp2_result_t) :: wmp
         type(error_t) :: werr
         integer, allocatable :: wz(:)
         real(dp), allocatable :: at_r_oh(:, :), smat(:, :), tmat(:, :), vmat(:, :)
         integer :: wnatm, probe_unit
         integer, parameter :: NA_OH = 5, NB_OH = 4, NFZ_OH = 1
         logical :: wh_converged

         call read_xyz("oh.xyz", wnatm, wz, at_r_oh)
         call trc_basis_from_json("basis_sets/6-31gs.json", wnatm, wz, at_r_oh, wbas, werr, force_cart=.true.)
         if (werr%has_error()) error stop "scf_dh: "//werr%get_message()
         call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", wnatm, wz, at_r_oh, waux, werr, force_cart=.true.)
         if (werr%has_error()) error stop "scf_dh: "//werr%get_message()
         call wbas%to_device()
         call waux%to_device()
         call wpl%build(wbas, 1.0e-12_dp)
         call wpl%to_device()

         ! S alone, this build's own AO order/normalisation -- never from an
         ! SCF, so dumped unconditionally (scf_rimp2_uhf.F90's own pattern).
         allocate (smat(wbas%nao, wbas%nao), tmat(wbas%nao, wbas%nao), vmat(wbas%nao, wbas%nao))
         call trc_1e(wbas, wpl, smat, tmat, vmat)
         open (newunit=probe_unit, file="scf_dh_oh_probe.bin", access="stream", form="unformatted", &
               status="replace", action="write")
         write (probe_unit) int(wbas%nao, int64)
         write (probe_unit) smat
         close (probe_unit)

         wopts = trc_scf_options_t()
         wopts%functional = "b2plyp"
         wopts%grid_level = 3
         wopts%conv_energy = 1.0e-12_dp
         wopts%conv_density = 1.0e-9_dp
         wopts%conv_diis = 1.0e-11_dp

         if (has_ref) then
            block
               integer :: nspin_chk, nao_chk
               real(dp), allocatable :: dm0(:, :, :)
               call ref_read_int(refu, "NSPIN_OH", nspin_chk)
               call ref_read_int(refu, "NAO_OH", nao_chk)
               if (nspin_chk /= 2 .or. nao_chk /= wbas%nao) &
                  error stop "scf_dh: pyscf reference shape (OH) does not match this build"
               allocate (dm0(wbas%nao, wbas%nao, 2))
               call ref_read_matrix3(refu, "DM0_OH", wbas%nao, wbas%nao, 2, dm0)
               call trc_scf_run(wbas, NA_OH, NB_OH, wopts, wres, dguess=dm0)
            end block
         else
            call trc_scf_run(wbas, NA_OH, NB_OH, wopts, wres)
         end if
         wh_converged = wres%converged
         if (.not. wh_converged) print "(a,a)", "  NOT CONVERGED (OH): ", trim(wres%message)
         block
            type(trc_rimp2_options_t) :: rimp2_opts
            rimp2_opts%nfrozen = NFZ_OH
            call trc_urimp2_run(wbas, waux, wpl, NA_OH, NB_OH, wres%cmo, wres%eps, wmp, opts=rimp2_opts)
         end block
         if (len_trim(wmp%message) > 0) print "(a)", "scf_dh: "//trim(wmp%message)
         print "(a)", ""
         print "(a,f20.12,a,f16.12,a,f16.12)", "  oh b2plyp UKS  E_SCF ", wres%energy, "  E_os ", wmp%e_os, &
            "  E_ss ", wmp%e_ss

         name_case(NCASE) = "oh_b2plyp_uks"
         if (has_ref) then
            block
               real(dp) :: e_scf_ref, e_os_ref, e_ss_ref, e_dh_ref, e_dh_this, worst
               call ref_read_real(refu, "E_SCF_OH", e_scf_ref)
               call ref_read_real(refu, "E_OS_OH", e_os_ref)
               call ref_read_real(refu, "E_SS_OH", e_ss_ref)
               call ref_read_real(refu, "E_DH_OH", e_dh_ref)
               call ref_read_real(refu, "CONV_TOL_OH", pyscf_conv_tol)
               call ref_read_real(refu, "CONV_TOL_GRAD_OH", pyscf_conv_tol_grad)
               call ref_read_real(refu, "GRAD_NORM_OH", pyscf_grad_norm)
               ! PT2_COS(1)/PT2_CSS(1) is b2plyp's literal 0.27/0.27 (FUNC(1)).
               e_dh_this = wres%energy + PT2_COS(1)*wmp%e_os + PT2_CSS(1)*wmp%e_ss
               worst = max(abs(wres%energy - e_scf_ref), max(abs(wmp%e_os - e_os_ref), &
                           max(abs(wmp%e_ss - e_ss_ref), abs(e_dh_this - e_dh_ref))))
               ok_case(NCASE) = wh_converged .and. abs(wres%energy - e_scf_ref) < CORR_TOL .and. &
                                abs(wmp%e_os - e_os_ref) < CORR_TOL .and. abs(wmp%e_ss - e_ss_ref) < CORR_TOL .and. &
                                abs(e_dh_this - e_dh_ref) < CORR_TOL .and. len_trim(wmp%message) == 0
               write (msg_case(NCASE), "(a,es10.2)") "oh b2plyp: E_SCF/E_os/E_ss/E_dh vs pyscf, worst = ", worst
               if (len_trim(wmp%message) > 0) msg_case(NCASE) = "RI-MP2: "//trim(wmp%message)//"; "// &
                                                                trim(msg_case(NCASE))
            end block
         end if

         call wpl%release()
         call waux%release()
         call wbas%release()
      end block

      if (has_ref) call ref_close(refu)
      if (.not. has_ref) error stop "FAIL: missing reference pyscf_refs/dh.ref -- run " &
         //"test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)"
   end subroutine run_all

end module scf_dh_suite

program scf_dh
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_dh_suite, only: collect_scf_dh
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_dh, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_dh
