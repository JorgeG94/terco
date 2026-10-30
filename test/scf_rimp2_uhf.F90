!
! UNRESTRICTED RI-MP2: UHF through trc_scf_run, then the correlation
! energy through trc_urimp2_run (src/post_hf/rimp2/trc_rimp2_driver.F90),
! on three open-shell radicals -- OH, triplet O2, NH2 -- in 6-31G* and
! cc-pVDZ with cc-pVDZ-RIFIT, six cases against pyscf's own UHF + DF-UMP2
! on the identical sets, plus one closed-shell-as-UHF self-check at the
! end (water, forced unrestricted with trc_scf_options_t%unrestricted,
! compared to the plain restricted trc_rimp2_run on the same molecule).
!
! OH and O2 are linear: the unpaired electron(s) sit in a rotationally
! degenerate pair of pi orbitals, a genuine zero-curvature direction of
! the energy -- the same thing test/scf_rks.F90's own OH cases hit first
! (see that file's header). test/pyscf_refs/generate.py's gen_rimp2_uhf
! computes a symmetry-broken initial guess (core Hamiltonian plus a field
! along one fixed axis) straight from the geometry and basis, converges
! pyscf's own UHF from it, and commits the CONVERGED density as DM0_i;
! this suite reads DM0_i back (if the reference already exists) and
! starts ITS OWN SCF from the SAME density, so both sides settle into the
! same orientation rather than two that merely share an energy. A missing
! reference (the first, bootstrap run, before generate.py has written
! one) falls back to trc_scf_run's own default guess, same as
! scf_rks.F90's identical fallback.
!
! This suite ALSO writes "rimp2_uhf_probe.bin": NAO and the overlap
! matrix S (trc_1e), in THIS build's own AO order, for each case --
! geometry- and basis-derived only, never from an SCF, so it is written
! unconditionally, the same "fixed input" category as rys_fock's density
! or pcm_ints's grid (test/pyscf_refs/README.md's own term for it). Its
! ONLY use is test/pyscf_refs/generate.py's gen_rimp2_uhf, which needs
! this build's own d-shell Cartesian normalisation (ao_norm_from_overlap)
! to scale the UHF density it commits as DM0_i correctly -- d and higher
! Cartesian cross components (dxy, dxz, dyz, ...) are NOT unit-normalised
! in this build the way pyscf's cart=True ones are (see
! test/pyscf_refs/README.md's own "File format" section), and skipping
! this rescale measured landing terco's own SCF in a DIFFERENT UHF
! stationary point on O2/6-31G* (|dE| up to 0.26 Ha) even though both
! sides started from "the same" density in name only.
!
! Seven independent cases below (six radical/basis pairs plus the
! closed-shell-as-UHF self-check), each its own test-drive unittest; the
! probe write and the sequential pyscf-reference read run once -- see
! check_xc_energy_suite's own `next_ifn` comment for why, and why this
! suite stays `parallel=.false.`.
!
module scf_rimp2_uhf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, trc_1e
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, trc_rimp2_options_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_matrix3
   implicit none
   private
   public :: collect_scf_rimp2_uhf

   integer, parameter :: NCASE = 6, NCASE_TOTAL = NCASE + 1
   character(len=8), parameter :: mol(NCASE) = [character(len=8) :: &
                                                "oh", "oh", "o2", "o2", "nh2", "nh2"]
   character(len=16), parameter :: basname(NCASE) = [character(len=16) :: &
                                                     "6-31gs", "cc-pvdz", "6-31gs", "cc-pvdz", &
                                                     "6-31gs", "cc-pvdz"]
   integer, parameter :: nalpha_c(NCASE) = [5, 5, 9, 9, 5, 5]
   integer, parameter :: nbeta_c(NCASE) = [4, 4, 7, 7, 4, 4]
   ! core_orbitals(Z), summed per molecule's atoms -- the SAME per-element
   ! rule trc_capi.F90's own default (nfrozen = -1) counts RI-MP2's frozen
   ! core by: O (Z=8) and N (Z=7) both 1, H 0. Explicit here because this
   ! suite calls trc_urimp2_run directly, under the driver's own
   ! contract (nfrozen defaults to none), not through that context default.
   integer, parameter :: nfroz_c(NCASE) = [1, 1, 2, 2, 1, 1]

   logical :: computed = .false.
   logical :: ok_case(NCASE_TOTAL)
   character(len=256) :: msg_case(NCASE_TOTAL)

contains

   subroutine collect_scf_rimp2_uhf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      integer :: ic
      allocate (testsuite(NCASE_TOTAL))
      do ic = 1, NCASE
         testsuite(ic) = new_unittest(trim(mol(ic))//"/"//trim(basname(ic)), test_case)
      end do
      testsuite(NCASE_TOTAL) = new_unittest("closed-shell-as-UHF vs restricted", test_case)
   end subroutine collect_scf_rimp2_uhf

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      integer, save :: next_ic = 0
      call ensure_computed()
      next_ic = next_ic + 1
      call check(error, ok_case(next_ic), trim(msg_case(next_ic)))
   end subroutine test_case

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   subroutine run_all()
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_rimp2_result_t) :: mp
      type(error_t) :: err
      logical :: has_ref
      integer :: ic, dev, refu, probe_unit
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), smat(:, :), tmat(:, :), vmat(:, :)
      integer :: natm
      character(len=16) :: tag
      logical :: converged_ok

      dev = trc_bind_device(0)
      ok_case = .false.
      msg_case = "pyscf reference missing"

      open (newunit=probe_unit, file="rimp2_uhf_probe.bin", access="stream", form="unformatted", status="replace", &
            action="write")
      write (probe_unit) int(NCASE, int64)

      call ref_open("rimp2_uhf.ref", refu, has_ref)
      if (has_ref) then
         block
            integer :: ncase_ref
            call ref_read_int(refu, "NCASE", ncase_ref)
            if (ncase_ref /= NCASE) error stop "scf_rimp2_uhf: pyscf reference case count does not match this build"
         end block
      end if

      do ic = 1, NCASE
         call read_xyz(trim(mol(ic))//".xyz", natm, zint, at_r)
         call trc_basis_from_json("basis_sets/"//trim(basname(ic))//".json", natm, zint, at_r, bas, err, &
                                  force_cart=.true.)
         if (err%has_error()) error stop "scf_rimp2_uhf: "//err%get_message()
         call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err, force_cart=.true.)
         if (err%has_error()) error stop "scf_rimp2_uhf: "//err%get_message()
         call bas%to_device()
         call aux%to_device()
         call pl%build(bas, 1.0e-12_dp)
         call pl%to_device()

         ! S alone, this build's own AO order -- see this module's own
         ! header for why it is dumped unconditionally.
         if (allocated(smat)) deallocate (smat, tmat, vmat)
         allocate (smat(bas%nao, bas%nao), tmat(bas%nao, bas%nao), vmat(bas%nao, bas%nao))
         call trc_1e(bas, pl, smat, tmat, vmat)
         write (probe_unit) int(bas%nao, int64)
         write (probe_unit) smat

         opts = trc_scf_options_t()
         opts%conv_energy = 1.0e-12_dp
         opts%conv_density = 1.0e-9_dp
         opts%conv_diis = 1.0e-11_dp   ! the commutator norm gates, not the density; see below
         write (tag, "(i0)") ic - 1

         if (has_ref) then
            block
               integer :: nspin_chk, nao_chk
               real(dp), allocatable :: dm0(:, :, :)
               call ref_read_int(refu, "NSPIN_"//trim(tag), nspin_chk)
               call ref_read_int(refu, "NAO_"//trim(tag), nao_chk)
               if (nspin_chk /= 2 .or. nao_chk /= bas%nao) &
                  error stop "scf_rimp2_uhf: pyscf reference shape does not match this build"
               allocate (dm0(bas%nao, bas%nao, 2))
               call ref_read_matrix3(refu, "DM0_"//trim(tag), bas%nao, bas%nao, 2, dm0)
               call trc_scf_run(bas, nalpha_c(ic), nbeta_c(ic), opts, res, dguess=dm0)
            end block
         else
            call trc_scf_run(bas, nalpha_c(ic), nbeta_c(ic), opts, res)
         end if
         converged_ok = res%converged
         if (.not. converged_ok) print "(a,a)", "  NOT CONVERGED: ", trim(res%message)

         block
            type(trc_rimp2_options_t) :: rimp2_opts
            rimp2_opts%nfrozen = nfroz_c(ic)
            call trc_urimp2_run(bas, aux, pl, nalpha_c(ic), nbeta_c(ic), res%cmo, res%eps, mp, opts=rimp2_opts)
         end block
         if (len_trim(mp%message) > 0) print "(a)", "scf_rimp2_uhf: "//trim(mp%message)

         print "(a,a4,1x,a8,a,f20.12,a,f16.12,a,f16.12)", "  ", trim(mol(ic)), trim(basname(ic)), &
            "  E_UHF ", res%energy, "  E_os ", mp%e_os, "  E_ss ", mp%e_ss

         if (has_ref) then
            block
               real(dp) :: e_scf_ref, e_os_ref, e_ss_ref, worst
               real(dp) :: pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
               call ref_read_real(refu, "E_SCF_"//trim(tag), e_scf_ref)
               call ref_read_real(refu, "E_OS_"//trim(tag), e_os_ref)
               call ref_read_real(refu, "E_SS_"//trim(tag), e_ss_ref)
               ! Provenance only -- see scf_rimp2.F90's own comment on the same tags.
               call ref_read_real(refu, "CONV_TOL_"//trim(tag), pyscf_conv_tol)
               call ref_read_real(refu, "CONV_TOL_GRAD_"//trim(tag), pyscf_conv_tol_grad)
               call ref_read_real(refu, "GRAD_NORM_"//trim(tag), pyscf_grad_norm)
               worst = max(abs(res%energy - e_scf_ref), max(abs(mp%e_os - e_os_ref), abs(mp%e_ss - e_ss_ref)))
               ok_case(ic) = converged_ok .and. len_trim(mp%message) == 0 .and. &
                             abs(res%energy - e_scf_ref) < 1.0e-9_dp .and. abs(mp%e_os - e_os_ref) < 1.0e-9_dp .and. &
                             abs(mp%e_ss - e_ss_ref) < 1.0e-9_dp
               write (msg_case(ic), "(a,es10.2)") trim(mol(ic))//"/"//trim(basname(ic))// &
                  ": E_UHF/E_os/E_ss vs pyscf, worst = ", worst
            end block
         end if

         call pl%release()
         call aux%release()
         call bas%release()
      end do
      if (has_ref) call ref_close(refu)
      close (probe_unit)

      ! --- closed-shell-as-UHF: water, forced unrestricted, against the ----
      ! --- plain restricted trc_rimp2_run on the same molecule/orbitals ----
      block
         type(trc_basis_t) :: wbas, waux
         type(trc_pairlist_t) :: wpl
         type(trc_scf_options_t) :: wopts
         type(trc_scf_result_t) :: wres_r, wres_u
         type(trc_rimp2_result_t) :: mp_r, mp_u
         type(error_t) :: werr
         integer, allocatable :: wz(:)
         real(dp), allocatable :: wr(:, :)
         integer :: wnatm
         real(dp) :: de_scf, de_corr

         call read_xyz("water.xyz", wnatm, wz, wr)
         call trc_basis_from_json("basis_sets/6-31gs.json", wnatm, wz, wr, wbas, werr, force_cart=.true.)
         if (werr%has_error()) error stop "scf_rimp2_uhf: "//werr%get_message()
         call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", wnatm, wz, wr, waux, werr, force_cart=.true.)
         if (werr%has_error()) error stop "scf_rimp2_uhf: "//werr%get_message()
         call wbas%to_device()
         call waux%to_device()
         call wpl%build(wbas, 1.0e-12_dp)
         call wpl%to_device()

         wopts = trc_scf_options_t()
         wopts%conv_energy = 1.0e-12_dp
         wopts%conv_density = 1.0e-9_dp
         wopts%conv_diis = 1.0e-11_dp
         ! The DIIS commutator norm, not the density, is what gates
         ! convergence (trc_scf_settings.F90's own comment on conv_diis); its
         ! 1e-5 default left ~1e-9 between this check's two independently
         ! converged SCFs, past the 1e-10 gate below -- measured.
         call trc_scf_run(wbas, 5, 5, wopts, wres_r)
         block
            type(trc_rimp2_options_t) :: rimp2_opts
            rimp2_opts%nfrozen = 1
            call trc_rimp2_run(wbas, waux, wpl, 5, wres_r%cmo(:, :, 1), wres_r%eps(:, 1), mp_r, opts=rimp2_opts)

            wopts%unrestricted = .true.
            call trc_scf_run(wbas, 5, 5, wopts, wres_u)
            call trc_urimp2_run(wbas, waux, wpl, 5, 5, wres_u%cmo, wres_u%eps, mp_u, opts=rimp2_opts)
         end block

         de_scf = abs(wres_r%energy - wres_u%energy)
         de_corr = abs(mp_r%e_corr - mp_u%e_corr)
         print "(a,f20.12,a,f20.12,a,es9.2)", "  closed-shell-as-UHF: E_corr(R) ", mp_r%e_corr, &
            "  E_corr(U) ", mp_u%e_corr, "  |diff| ", de_corr

         ok_case(NCASE_TOTAL) = de_scf <= 1.0e-10_dp .and. de_corr <= 1.0e-10_dp
         write (msg_case(NCASE_TOTAL), "(a,es9.2,a,es9.2)") &
            "closed-shell-as-UHF vs restricted: |dE_SCF| = ", de_scf, ", |dE_corr| = ", de_corr

         call wpl%release()
         call waux%release()
         call wbas%release()
      end block
   end subroutine run_all

end module scf_rimp2_uhf_suite

program scf_rimp2_uhf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_rimp2_uhf_suite, only: collect_scf_rimp2_uhf
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_rimp2_uhf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_rimp2_uhf
