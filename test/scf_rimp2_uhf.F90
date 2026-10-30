!
! UNRESTRICTED RI-MP2: UHF through trc_scf_run, then the correlation
! energy through trc_urimp2_run (src/post_hf/rimp2/trc_rimp2_driver.F90),
! on three open-shell radicals -- OH, triplet O2, NH2 -- in 6-31G* and
! cc-pVDZ with cc-pVDZ-RIFIT, six cases against pyscf's own UHF + DF-UMP2
! on the identical sets, plus one closed-shell-as-UHF self-check at the
! bottom (water, forced unrestricted with trc_scf_options_t%unrestricted,
! compared to the plain restricted trc_rimp2_run on the same molecule).
!
! OH and O2 are linear: the unpaired electron(s) sit in a rotationally
! degenerate pair of pi orbitals, a genuine zero-curvature direction of
! the energy -- the same thing test/scf_rks.F90's own OH cases hit first
! (see that file's header). test/pyscf_refs/generate.py's gen_rimp2_uhf
! computes a symmetry-broken initial guess (core Hamiltonian plus a field
! along one fixed axis) straight from the geometry and basis, converges
! pyscf's own UHF from it, and commits the CONVERGED density as DM0_i;
! this program reads DM0_i back (if the reference already exists) and
! starts ITS OWN SCF from the SAME density, so both sides settle into the
! same orientation rather than two that merely share an energy. A missing
! reference (the first, bootstrap run, before generate.py has written
! one) falls back to trc_scf_run's own default guess, same as
! scf_rks.F90's identical fallback.
!
! This program ALSO writes "rimp2_uhf_probe.bin": NAO and the overlap
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
program scf_rimp2_uhf
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device, trc_1e
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_urimp2_run, trc_rimp2_result_t, trc_rimp2_options_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_matrix3, ref_report
   implicit none

   integer, parameter :: NCASE = 6
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
   ! program calls trc_urimp2_run directly, under the driver's own
   ! contract (nfrozen defaults to none), not through that context default.
   integer, parameter :: nfroz_c(NCASE) = [1, 1, 2, 2, 1, 1]

   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(trc_rimp2_result_t) :: mp
   type(error_t) :: err
   logical :: ok, has_ref
   integer :: ic, dev, refu, nfail, probe_unit
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :), smat(:, :), tmat(:, :), vmat(:, :)
   integer :: natm
   character(len=16) :: tag

   dev = trc_bind_device(0)
   ok = .true.
   nfail = 0

   open (newunit=probe_unit, file="rimp2_uhf_probe.bin", access="stream", form="unformatted", status="replace", &
         action="write")
   write (probe_unit) int(NCASE, int64)

   call ref_open("rimp2_uhf.ref", refu, has_ref)
   if (has_ref) then
      block
         integer :: ncase_ref
         call ref_read_int(refu, "NCASE", ncase_ref)
         if (ncase_ref /= NCASE) then
            print "(a)", "scf_rimp2_uhf: pyscf reference case count does not match this build"
            stop 1
         end if
      end block
   end if

   do ic = 1, NCASE
      call read_xyz(trim(mol(ic))//".xyz", natm, zint, at_r)
      call trc_basis_from_json("basis_sets/"//trim(basname(ic))//".json", natm, zint, at_r, bas, err, &
                               force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_rimp2_uhf: "//err%get_message()
         stop 1
      end if
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err, force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_rimp2_uhf: "//err%get_message()
         stop 1
      end if
      call bas%to_device()
      call aux%to_device()
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()

      ! S alone, this build's own AO order -- see this program's own
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
            if (nspin_chk /= 2 .or. nao_chk /= bas%nao) then
               print "(a)", "scf_rimp2_uhf: pyscf reference shape does not match this build"
               stop 1
            end if
            allocate (dm0(bas%nao, bas%nao, 2))
            call ref_read_matrix3(refu, "DM0_"//trim(tag), bas%nao, bas%nao, 2, dm0)
            call trc_scf_run(bas, nalpha_c(ic), nbeta_c(ic), opts, res, dguess=dm0)
         end block
      else
         call trc_scf_run(bas, nalpha_c(ic), nbeta_c(ic), opts, res)
      end if
      if (.not. res%converged) then
         print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
         ok = .false.
      end if

      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%nfrozen = nfroz_c(ic)
         call trc_urimp2_run(bas, aux, pl, nalpha_c(ic), nbeta_c(ic), res%cmo, res%eps, mp, opts=rimp2_opts)
      end block
      if (len_trim(mp%message) > 0) then
         print "(a)", "scf_rimp2_uhf: "//trim(mp%message)
         ok = .false.
      end if

      print "(a,a4,1x,a8,a,f20.12,a,f16.12,a,f16.12)", "  ", trim(mol(ic)), trim(basname(ic)), &
         "  E_UHF ", res%energy, "  E_os ", mp%e_os, "  E_ss ", mp%e_ss

      if (has_ref) then
         block
            real(dp) :: e_scf_ref, e_os_ref, e_ss_ref
            real(dp) :: pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
            call ref_read_real(refu, "E_SCF_"//trim(tag), e_scf_ref)
            call ref_read_real(refu, "E_OS_"//trim(tag), e_os_ref)
            call ref_read_real(refu, "E_SS_"//trim(tag), e_ss_ref)
            ! Provenance only -- see scf_rimp2.F90's own comment on the same tags.
            call ref_read_real(refu, "CONV_TOL_"//trim(tag), pyscf_conv_tol)
            call ref_read_real(refu, "CONV_TOL_GRAD_"//trim(tag), pyscf_conv_tol_grad)
            call ref_read_real(refu, "GRAD_NORM_"//trim(tag), pyscf_grad_norm)
            call ref_report(nfail, abs(res%energy - e_scf_ref) < 1.0e-9_dp, &
                           trim(mol(ic))//"/"//trim(basname(ic))//" E_UHF vs pyscf", abs(res%energy - e_scf_ref))
            call ref_report(nfail, abs(mp%e_os - e_os_ref) < 1.0e-9_dp, &
                           trim(mol(ic))//"/"//trim(basname(ic))//" E_os vs pyscf", abs(mp%e_os - e_os_ref))
            call ref_report(nfail, abs(mp%e_ss - e_ss_ref) < 1.0e-9_dp, &
                           trim(mol(ic))//"/"//trim(basname(ic))//" E_ss vs pyscf", abs(mp%e_ss - e_ss_ref))
         end block
      end if

      call pl%release()
      call aux%release()
      call bas%release()
   end do
   if (has_ref) call ref_close(refu)
   close (probe_unit)
   if (nfail > 0) ok = .false.

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
      if (werr%has_error()) then
         print "(a)", "scf_rimp2_uhf: "//werr%get_message()
         stop 1
      end if
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", wnatm, wz, wr, waux, werr, force_cart=.true.)
      if (werr%has_error()) then
         print "(a)", "scf_rimp2_uhf: "//werr%get_message()
         stop 1
      end if
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
      if (de_scf > 1.0e-10_dp) then
         print "(a,es9.2)", "  FAIL: closed-shell-as-UHF SCF energy vs restricted, |diff| ", de_scf
         ok = .false.
      end if
      if (de_corr > 1.0e-10_dp) then
         print "(a,es9.2)", "  FAIL: closed-shell-as-UHF E_corr vs restricted, |diff| ", de_corr
         ok = .false.
      end if

      call wpl%release()
      call waux%release()
      call wbas%release()
   end block

   print "(a)", merge("scf_rimp2_uhf: PASS", "scf_rimp2_uhf: FAIL", ok)
   if (.not. ok) stop 1
end program scf_rimp2_uhf
