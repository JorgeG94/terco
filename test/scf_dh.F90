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
! c_os*e_os + c_ss*e_ss) as E_DH -- this program's own c_os*E_os +
! c_ss*E_ss weighting (from LITERAL coefficients, PT2_COS/PT2_CSS below,
! not read back from trc_xc_functional's own table) is checked against
! that committed E_DH directly, never left as an unchecked header claim.
!
! 24 closed-shell cases: water and formamide, cc-pVDZ and 6-31G*,
! b2plyp, b2gp-plyp and mpw2plyp, four-centre and --df -- plus one open-shell case,
! OH/b2plyp UKS/6-31G*, from a FIXED initial density (rks.ref's and
! rimp2_uhf.ref's own pattern: pyscf's own converged density, rescaled
! into this build's AO normalisation via this program's own dumped
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
! same xyz/JSON files, scf_df.F90's own pattern); for OH, nao and the
! overlap matrix S, unconditionally (geometry- and basis-derived only,
! never from an SCF -- scf_rimp2_uhf.F90's own "fixed input" category).
!
program scf_dh
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, trc_1e
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_rimp2_driver, only: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, trc_rimp2_options_t
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3, ref_report
   implicit none

   integer, parameter :: NGROUP = 4, NFUNC = 3
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

   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_linalg_t) :: la
   type(df_jk_t), target :: dfjk
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(trc_rimp2_result_t) :: mp
   type(error_t) :: err
   logical :: all_ok, has_ref
   integer :: ig, im, ik, unit, natm, nalpha, nbeta, nelec, ic, nfail, refu, ncase_ref
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   character(len=32) :: tag, name_ref, name_this
   real(dp) :: e_ref, de, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
   integer :: dev

   dev = trc_bind_device(0)
   all_ok = .true.
   nfail = 0

   call ref_open("dh.ref", refu, has_ref)
   if (has_ref) then
      call ref_read_int(refu, "NCASE", ncase_ref)
      if (ncase_ref /= NGROUP*NFUNC*2) then
         print "(a)", "scf_dh: pyscf reference case count does not match this build"
         stop 1
      end if
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
      if (err%has_error()) then
         print "(a)", "scf_dh: "//err%get_message()
         stop 1
      end if
      call bas%to_device()
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err, force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_dh: "//err%get_message()
         stop 1
      end if
      call aux%to_device()

      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      call la%init(bas%nao)
      ! --df's own fit, cc-pVDZ-RIFIT doing double duty as the jkfit too:
      ! set up once per group, reused by both functionals' --df case and
      ! by RI-MP2's chained path (trc_rimp2_run's `packed_b` argument).
      call dfjk%setup(bas, aux, pl, la, err, nocc_hint=max(nalpha, 1))
      if (err%has_error()) then
         print "(a)", "scf_dh: "//err%get_message()
         stop 1
      end if

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
            if (.not. res%converged) then
               print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
               all_ok = .false.
            end if
            if (len_trim(mp%message) > 0) then
               print "(a)", "scf_dh: "//trim(mp%message)
               all_ok = .false.
            end if
            print "(a,a8,1x,a2,a,i3,a,f20.12,a,f16.12,a,f16.12)", "  ", FUNC(im), KIND(ik), " iter ", &
               res%iterations, "  E_SCF ", res%energy, "  E_os ", mp%e_os, "  E_ss ", mp%e_ss

            ic = ((ig - 1)*NFUNC + (im - 1))*2 + ik
            name_this = trim(GMOL(ig))//"_"//trim(GBASIS(ig))//"_"//trim(FUNC(im))//"_"//trim(KIND(ik))
            if (has_ref) then
               write (tag, "(i0)") ic - 1
               call ref_read_str(refu, "NAME_"//trim(tag), name_ref)
               call ref_read_real(refu, "E_SCF_"//trim(tag), e_ref)
               if (trim(name_ref) /= trim(name_this)) then
                  print "(a)", "scf_dh: pyscf reference case order does not match this build "// &
                     "("//trim(name_ref)//" vs "//trim(name_this)//")"
                  stop 1
               end if
               de = abs(res%energy - e_ref)
               call ref_report(nfail, res%converged .and. de < merge(E_TOL_4C, E_TOL_DF, ik == 1), &
                               trim(name_this)//" E_SCF vs pyscf", de)
               block
                  real(dp) :: e_os_ref, e_ss_ref, e_dh_ref, e_dh_this
                  call ref_read_real(refu, "E_OS_"//trim(tag), e_os_ref)
                  call ref_read_real(refu, "E_SS_"//trim(tag), e_ss_ref)
                  call ref_report(nfail, abs(mp%e_os - e_os_ref) < CORR_TOL, &
                                  trim(name_this)//" E_os vs pyscf", abs(mp%e_os - e_os_ref))
                  call ref_report(nfail, abs(mp%e_ss - e_ss_ref) < CORR_TOL, &
                                  trim(name_this)//" E_ss vs pyscf", abs(mp%e_ss - e_ss_ref))
                  ! The weighted double-hybrid total, against pyscf's OWN
                  ! mf.e_tot + c_os*e_os + c_ss*e_ss (generate.py's gen_dh),
                  ! not merely this program's own E_SCF/E_os/E_ss checked
                  ! separately above -- PT2_COS/PT2_CSS are the SAME
                  ! literals context_scf.F90's check_double_hybrid asserts.
                  call ref_read_real(refu, "E_DH_"//trim(tag), e_dh_ref)
                  e_dh_this = res%energy + PT2_COS(im)*mp%e_os + PT2_CSS(im)*mp%e_ss
                  call ref_report(nfail, abs(e_dh_this - e_dh_ref) < CORR_TOL, &
                                  trim(name_this)//" weighted E_dh vs pyscf", abs(e_dh_this - e_dh_ref))
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

      call read_xyz("oh.xyz", wnatm, wz, at_r_oh)
      call trc_basis_from_json("basis_sets/6-31gs.json", wnatm, wz, at_r_oh, wbas, werr, force_cart=.true.)
      if (werr%has_error()) then
         print "(a)", "scf_dh: "//werr%get_message()
         stop 1
      end if
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", wnatm, wz, at_r_oh, waux, werr, force_cart=.true.)
      if (werr%has_error()) then
         print "(a)", "scf_dh: "//werr%get_message()
         stop 1
      end if
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
            if (nspin_chk /= 2 .or. nao_chk /= wbas%nao) then
               print "(a)", "scf_dh: pyscf reference shape (OH) does not match this build"
               stop 1
            end if
            allocate (dm0(wbas%nao, wbas%nao, 2))
            call ref_read_matrix3(refu, "DM0_OH", wbas%nao, wbas%nao, 2, dm0)
            call trc_scf_run(wbas, NA_OH, NB_OH, wopts, wres, dguess=dm0)
         end block
      else
         call trc_scf_run(wbas, NA_OH, NB_OH, wopts, wres)
      end if
      if (.not. wres%converged) then
         print "(a,a)", "  NOT CONVERGED (OH): ", trim(wres%message)
         all_ok = .false.
      end if
      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%nfrozen = NFZ_OH
         call trc_urimp2_run(wbas, waux, wpl, NA_OH, NB_OH, wres%cmo, wres%eps, wmp, opts=rimp2_opts)
      end block
      if (len_trim(wmp%message) > 0) then
         print "(a)", "scf_dh: "//trim(wmp%message)
         all_ok = .false.
      end if
      print "(a)", ""
      print "(a,f20.12,a,f16.12,a,f16.12)", "  oh b2plyp UKS  E_SCF ", wres%energy, "  E_os ", wmp%e_os, &
         "  E_ss ", wmp%e_ss

      if (has_ref) then
         block
            real(dp) :: e_scf_ref, e_os_ref, e_ss_ref, e_dh_ref, e_dh_this
            call ref_read_real(refu, "E_SCF_OH", e_scf_ref)
            call ref_read_real(refu, "E_OS_OH", e_os_ref)
            call ref_read_real(refu, "E_SS_OH", e_ss_ref)
            call ref_read_real(refu, "E_DH_OH", e_dh_ref)
            call ref_read_real(refu, "CONV_TOL_OH", pyscf_conv_tol)
            call ref_read_real(refu, "CONV_TOL_GRAD_OH", pyscf_conv_tol_grad)
            call ref_read_real(refu, "GRAD_NORM_OH", pyscf_grad_norm)
            call ref_report(nfail, abs(wres%energy - e_scf_ref) < CORR_TOL, "oh b2plyp E_SCF vs pyscf", &
                            abs(wres%energy - e_scf_ref))
            call ref_report(nfail, abs(wmp%e_os - e_os_ref) < CORR_TOL, "oh b2plyp E_os vs pyscf", &
                            abs(wmp%e_os - e_os_ref))
            call ref_report(nfail, abs(wmp%e_ss - e_ss_ref) < CORR_TOL, "oh b2plyp E_ss vs pyscf", &
                            abs(wmp%e_ss - e_ss_ref))
            ! PT2_COS(1)/PT2_CSS(1) is b2plyp's literal 0.27/0.27 (FUNC(1)).
            e_dh_this = wres%energy + PT2_COS(1)*wmp%e_os + PT2_CSS(1)*wmp%e_ss
            call ref_report(nfail, abs(e_dh_this - e_dh_ref) < CORR_TOL, "oh b2plyp weighted E_dh vs pyscf", &
                            abs(e_dh_this - e_dh_ref))
         end block
      end if

      call wpl%release()
      call waux%release()
      call wbas%release()
   end block

   if (has_ref) call ref_close(refu)
   if (.not. has_ref) then
      print "(a)", "FAIL: missing reference pyscf_refs/dh.ref -- run test/pyscf_refs/generate.py "// &
         "(see test/pyscf_refs/README.md)"
      stop 1
   end if
   if (.not. all_ok .or. nfail > 0) stop 1
   print "(a)", ""
   print "(a)", "scf_dh: every case converged and agrees with the pyscf reference; RESULT: PASS"

end program scf_dh
