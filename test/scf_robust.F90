!
! The hard-case set for SCF robustness (M2, Part B.4): four small, fast
! systems chosen because plain damped DIIS is exactly where each one is
! weak, each for a different reason --
!
!   n2_stretched  UHF/6-31G, charge +1   N2+ (13 electrons, genuinely
!                              open shell) with the bond stretched to 2x
!                              Re: the pi/pi* gap that is large at Re is
!                              small here, and the commutator stays large
!                              long after a normal molecule would have
!                              reached the DIIS window.
!   o3            RHF/6-31G    ozone: a small HOMO-LUMO gap and
!                              biradicaloid character, a textbook case
!                              simple damping oscillates on.
!   w5            RHF/d-aug-cc-pVDZ   five waters, doubly diffuse-augmented:
!                              S's smallest eigenvalues measured at 6.5e-7
!                              to 2.3e-5 (dumped with a throwaway syev
!                              driver, not shown here), so canonical
!                              orthogonalisation is asked for at 1e-5 --
!                              looser than the type's 1e-8 default, which
!                              a single aug-cc-pVDZ or def2-TZVPD calc at
!                              this size does not reach (measured: aug-cc-
!                              pVDZ on the same five waters bottoms out at
!                              4e-4, three orders above where this engages)
!                              but a doubly-augmented or larger diffuse
!                              cluster does.
!   o2_anion      UHF/6-31G, charge -1   superoxide: charged AND open
!                              shell, with the pi*/pi* near-degeneracy
!                              that makes it a standard hard DFT/HF case.
!
! This runs each with the SETTINGS THIS MILESTONE ADDS ON (A/EDIIS ahead
! of CDIIS where the case needs it, canonical orthogonalisation always
! live at its measured default threshold) and asserts convergence --
! that is the regression ctest protects. The BEFORE/AFTER comparison
! against plain damped DIIS and against an uncorrected S is numbers for
! the milestone report, not something this binary can run safely: the
! old, unconditional `sym_orthog` is an `error stop` on a singular S, and
! `error stop` takes the whole process with it, which a ctest case must
! not risk. That sweep is a separate, disposable driver, one process per
! configuration.
!
! Writes scf_robust_probe.bin for scf_robust_ref.py, same layout idea as
! scf_rks.F90's rks_probe.bin but on terco's own (non-uncontracted) basis,
! since what is being checked here is "does it converge to the right
! energy on a hard case", not a bit-level fold convention.
!
! Four independent cases below, each its own test-drive unittest; the
! probe write and the sequential pyscf-reference read run once -- see
! check_xc_energy_suite's own `next_ifn` comment for why, and why this
! suite stays `parallel=.false.`.
!
module scf_robust_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none
   private
   public :: collect_scf_robust

   integer, parameter :: NCASE = 4
   character(len=16), parameter :: xyzf(NCASE) = [character(len=16) :: &
      "n2_stretched", "o3", "w5", "o2_anion"]
   character(len=16), parameter :: basf(NCASE) = [character(len=16) :: &
      "6-31g", "6-31g", "d-aug-cc-pvdz", "6-31g"]
   logical, parameter :: uhf(NCASE) = [.true., .false., .false., .true.]
   integer, parameter :: charge(NCASE) = [1, 0, 0, -1]
   real(dp), parameter :: lindep(NCASE) = [1.0e-8_dp, 1.0e-8_dp, 1.0e-5_dp, 1.0e-8_dp]
   real(dp), parameter :: pyscf_tol(NCASE) = [2.0e-5_dp, 2.0e-5_dp, 3.0e-4_dp, 2.0e-5_dp]
   !! Agreement with pyscf's own converged energy (test/pyscf_refs/
   !! scf_robust.ref). Looser than the uncontracted comparisons elsewhere in
   !! this tree because the basis here is CONTRACTED and read by both codes
   !! from a name, not terco's exact exponents -- see this module's header.
   !! w5/d-aug-cc-pVDZ gets the loosest gate: canonical orthogonalisation
   !! drops functions below `lindep_thresh` there, and pyscf's default RHF
   !! applies no such cut, so the two are not solving the same, numerically
   !! well-posed problem to begin with.

   logical :: computed = .false.
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)

contains

   subroutine collect_scf_robust(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      integer :: ic
      allocate (testsuite(NCASE))
      do ic = 1, NCASE
         testsuite(ic) = new_unittest(trim(xyzf(ic))//"/"//trim(basf(ic)), test_case)
      end do
   end subroutine collect_scf_robust

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
      integer :: ic, unit, natm, nelec, nalpha, nbeta, dev
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer :: refu, ncase_ref
      character(len=16) :: tag
      logical :: has_ref

      dev = trc_bind_device(0)
      open (newunit=unit, file='scf_robust_probe.bin', access='stream', form='unformatted', status='replace')
      write (unit) int(NCASE, int64)

      ! Opened (not required) before any case runs -- a missing reference
      ! must not stop this suite before scf_robust_probe.bin exists.
      call ref_open('scf_robust.ref', refu, has_ref)
      call golden_begin('scf_robust')
      if (has_ref) then
         call ref_read_int(refu, 'NCASE', ncase_ref)
         if (ncase_ref /= NCASE) error stop 'scf_robust: pyscf reference case count does not match this build'
      end if
      ok_case = .false.
      msg_case = "pyscf reference missing"

      do ic = 1, NCASE
         call read_xyz(trim(xyzf(ic))//'.xyz', natm, zint, at_r)
         nelec = sum(zint) - charge(ic)
         if (uhf(ic)) then
            nalpha = (nelec + 1)/2
            nbeta = nelec - nalpha
         else
            nalpha = nelec/2
            nbeta = nalpha
         end if
         ! Cartesian, to match scf_robust_ref.py's cart=True and the lindep
         ! threshold measured on it: d-aug-cc-pVDZ tags shells gto_spherical.
         ! overrides%prune = .false.: w5/d-aug-cc-pVDZ's lindep_thresh (1e-5, this
         ! module's own header) drops near-linearly-dependent diffuse
         ! combinations off the overlap matrix's eigenvalues, and
         ! prune_block's basis rotation moves those eigenvalues enough to
         ! change WHICH combinations get dropped -- a real difference, not
         ! round-off (measured 9.66e-07 against the golden, four orders of
         ! magnitude over its 2e-10 gate), because this basis is deliberately
         ! near-redundant and lindep thresholding is not invariant to a
         ! rotation the way a well-conditioned SCF energy is.
         call trc_basis_from_json('basis_sets/'//trim(basf(ic))//'.json', natm, zint, at_r, bas, err, &
                                  force_cart=.true., overrides=trc_basis_json_overrides(prune=.false.))
         if (err%has_error()) error stop 'scf_robust: '//trim(xyzf(ic))//': '//err%get_message()
         call bas%to_device()

         opts = trc_scf_options_t()
         opts%max_iteration = 150
         opts%unrestricted = uhf(ic)
         opts%adiis = .true.
         opts%lindep_thresh = lindep(ic)

         call trc_scf_run(bas, nalpha, nbeta, opts, res)
         print '(a,a14,a,a12,a,i4,a,l1,a,f20.10,a,i0)', &
            '  ', trim(xyzf(ic)), '  ', trim(basf(ic)), '  iter ', res%iterations, &
            '  conv ', res%converged, '  E = ', res%energy, '  ndrop ', res%n_lindep
         if (.not. res%converged) print '(a,a)', '    NOT CONVERGED: ', trim(res%message)

         write (unit) int(bas%nao, int64), int(natm, int64), int(res%nspin, int64), &
            int(nalpha, int64), int(nbeta, int64), int(charge(ic), int64)
         write (unit) xyzf(ic), basf(ic)
         write (unit) res%energy
         call golden_put('E_'//trim(xyzf(ic))//'_'//trim(basf(ic)), res%energy)
         call bas%release()

         ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
         if (has_ref) then
            block
               character(len=16) :: xyz_ref, bas_ref
               real(dp) :: e_ref, de, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
               write (tag, '(i0)') ic - 1
               call ref_read_str(refu, 'XYZ_'//trim(tag), xyz_ref)
               call ref_read_str(refu, 'BAS_'//trim(tag), bas_ref)
               call ref_read_real(refu, 'E_'//trim(tag), e_ref)
               ! Provenance only -- see scf_rks.F90's own comment on the same tags.
               call ref_read_real(refu, 'CONV_TOL_'//trim(tag), pyscf_conv_tol)
               call ref_read_real(refu, 'CONV_TOL_GRAD_'//trim(tag), pyscf_conv_tol_grad)
               call ref_read_real(refu, 'GRAD_NORM_'//trim(tag), pyscf_grad_norm)
               if (trim(xyz_ref) /= trim(xyzf(ic)) .or. trim(bas_ref) /= trim(basf(ic))) &
                  error stop 'scf_robust: pyscf reference case order does not match this build'
               de = abs(res%energy - e_ref)
               ok_case(ic) = res%converged .and. de < pyscf_tol(ic)
               write (msg_case(ic), '(a,es10.2)') trim(xyzf(ic))//'/'//trim(basf(ic))//' vs pyscf, |dE| = ', de
            end block
         end if
      end do
      call golden_end()
      close (unit)
      if (has_ref) call ref_close(refu)
      if (.not. has_ref) error stop 'FAIL: missing reference pyscf_refs/scf_robust.ref -- run ' &
         //'test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)'
   end subroutine run_all

end module scf_robust_suite

program scf_robust
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_robust_suite, only: collect_scf_robust
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_robust, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_robust
