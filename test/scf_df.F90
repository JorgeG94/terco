!
! Density-fitted RHF/RKS (`--df`'s library entry, trc_scf_df.F90) against
! pyscf's own `density_fit`, same auxiliary basis (def2-universal-jkfit) on
! both sides -- water and formamide, cc-pVDZ and 6-31G*, Hartree-Fock, PBE
! and B3LYP. Converged energy only, at rks.ref's closed-shell tolerance
! (1e-8): no at-density control (unlike rks.ref itself), since that needs
! an AO-norm-corrected density transfer this case does not have.
!
! The dump carries the molecule/basis NAMES (and nao/natm/nalpha/nbeta as
! a sanity check), not shell data: test/pyscf_refs/generate.py's
! gen_scf_df rebuilds the identical basis itself, straight from the same
! xyz file and the same BSE JSON this program reads (bse_json_basis, the
! scf_robust.F90/check_rys_fock_vs_pyscf.F90 pattern) -- not from terco's
! OWN per-primitive sh_c/sh_e, which fold in a primitive normalisation
! pyscf's custom-basis-dict path applies a second time if handed those
! values directly (found debugging this very case: every contraction with
! more than one primitive came out wrong, confirmed against the exact
! same basis built via bse_json_basis instead).
!
! NEEDS DF_LMAX_AUX >= 4: def2-universal-jkfit has g auxiliary functions
! on O/C/N. Every case self-skips (test-drive's skip_test) rather than
! failing when `maxval(aux%sh_l) > DF_LMAX_AUX` -- a build below that
! ceiling.
!
! Twelve independent (molecule, basis, method) cases below, each its own
! test-drive unittest; the shared setup (per-group basis/DF-JK build, the
! probe write, the golden dump and the pyscf-reference read) runs once --
! see check_xc_energy_suite's own `next_ifn` comment for why, and why
! this suite stays `parallel=.false.`.
!
module scf_df_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none
   private
   public :: collect_scf_df

   integer, parameter :: NGROUP = 4, NMETHOD = 3, NCASE = NGROUP*NMETHOD
   character(len=16), parameter :: GMOL(NGROUP) = [character(len=16) :: &
                                                   "water", "water", "formamide", "formamide"]
   character(len=16), parameter :: GBASIS(NGROUP) = [character(len=16) :: &
                                                      "cc-pvdz", "6-31gs", "cc-pvdz", "6-31gs"]
   character(len=8), parameter :: METHOD(NMETHOD) = [character(len=8) :: "", "pbe", "b3lyp"]
   real(dp), parameter :: E_TOL = 1.0e-8_dp   !! rks.ref's closed-shell tolerance; every case here is nspin=1

   logical :: computed = .false.
   logical :: skip_all = .false.
   character(len=256) :: skip_message = ""
   logical :: ok_case(NCASE)
   character(len=256) :: msg_case(NCASE)
   character(len=32) :: name_case(NCASE)

contains

   subroutine collect_scf_df(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      integer :: ig, im, ic
      allocate (testsuite(NCASE))
      do ig = 1, NGROUP
         do im = 1, NMETHOD
            ic = (ig - 1)*NMETHOD + im
            if (len_trim(METHOD(im)) > 0) then
               testsuite(ic) = new_unittest(trim(GMOL(ig))//"/"//trim(GBASIS(ig))//"/"//trim(METHOD(im)), &
                                            test_case)
            else
               testsuite(ic) = new_unittest(trim(GMOL(ig))//"/"//trim(GBASIS(ig))//"/hf", test_case)
            end if
         end do
      end do
   end subroutine collect_scf_df

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      integer, save :: next_ic = 0
      call ensure_computed()
      next_ic = next_ic + 1
      if (skip_all) then
         call skip_test(error, trim(skip_message))
         return
      end if
      call check(error, ok_case(next_ic), trim(msg_case(next_ic)))
   end subroutine test_case

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   subroutine run_all()
      type(trc_basis_t) :: bas, auxb
      type(trc_pairlist_t) :: pl
      type(trc_linalg_t) :: la
      type(df_jk_t) :: dfjk
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer :: ig, im, unit, natm, nalpha, nbeta, nelec
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: refu, ncase_ref, ic
      logical :: has_ref
      character(len=32) :: tag, name_ref, name_this
      real(dp) :: e_ref, de, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm

      call ref_open("scf_df.ref", refu, has_ref)
      call golden_begin("scf_df")
      if (.not. has_ref) error stop "FAIL: missing reference pyscf_refs/scf_df.ref -- run " &
         //"test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)"
      call ref_read_int(refu, "NCASE", ncase_ref)
      if (ncase_ref /= NCASE) error stop "scf_df: pyscf reference case count does not match this build"

      open (newunit=unit, file="scf_df_probe.bin", access="stream", form="unformatted", status="replace", &
            action="write")
      write (unit) int(NCASE, int64)

      ok_case = .false.
      msg_case = ""

      do ig = 1, NGROUP
         call read_xyz(trim(GMOL(ig))//".xyz", natm, zint, at_r)
         nelec = sum(zint)
         nalpha = (nelec + mod(nelec, 2))/2
         nbeta = nelec - nalpha

         call trc_basis_from_json("basis_sets/"//trim(GBASIS(ig))//".json", natm, zint, at_r, bas, err, &
                                  force_cart=.true.)
         if (err%has_error()) error stop "scf_df: "//err%get_message()
         call bas%to_device()
         call trc_basis_from_json("basis_sets/def2-universal-jkfit.json", natm, zint, at_r, auxb, err, &
                                  force_cart=.true.)
         if (err%has_error()) error stop "scf_df: "//err%get_message()
         if (maxval(auxb%sh_l) > DF_LMAX_AUX) then
            write (skip_message, "(a,i0,a,i0,a)") "def2-universal-jkfit needs l = ", maxval(auxb%sh_l), &
               " and this build has DF_LMAX_AUX = ", DF_LMAX_AUX, "."
            skip_all = .true.
            call auxb%release()
            call bas%release()
            call ref_close(refu)
            close (unit)
            return
         end if
         call auxb%to_device()

         call pl%build(bas, 1.0e-12_dp)
         call pl%to_device()
         call la%init(bas%nao)
         call dfjk%setup(bas, auxb, pl, la, err, nocc_hint=max(nalpha, 1))
         if (err%has_error()) error stop "scf_df: "//err%get_message()
         call pl%release()

         print "(a)", ""
         print "(a,a,a,a,a,i0,a)", "scf_df: ", trim(GMOL(ig)), "/", trim(GBASIS(ig)), ", ", bas%nao, " AOs"

         do im = 1, NMETHOD
            opts = trc_scf_options_t()
            opts%df = .true.
            opts%functional = METHOD(im)
            opts%grid_level = 3
            opts%conv_energy = 1.0e-11_dp
            opts%conv_diis = 1.0e-8_dp
            call trc_scf_run(bas, nalpha, nbeta, opts, res, la=la, scf_dfjk=dfjk)
            if (.not. res%converged) print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
            if (len_trim(METHOD(im)) > 0) then
               print "(a,a8,a,i3,a,f22.12)", "  ", METHOD(im), " iter ", res%iterations, "  E = ", res%energy
            else
               print "(a,a8,a,i3,a,f22.12)", "  ", "hf", " iter ", res%iterations, "  E = ", res%energy
            end if

            ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
            ic = (ig - 1)*NMETHOD + im
            if (len_trim(METHOD(im)) > 0) then
               name_this = trim(GMOL(ig))//"_"//trim(GBASIS(ig))//"_"//trim(METHOD(im))
            else
               name_this = trim(GMOL(ig))//"_"//trim(GBASIS(ig))//"_hf"
            end if
            name_case(ic) = name_this
            call golden_put("E_"//trim(name_this), res%energy)

            write (tag, "(i0)") ic - 1
            call ref_read_str(refu, "NAME_"//trim(tag), name_ref)
            call ref_read_real(refu, "E_"//trim(tag), e_ref)
            ! Provenance only, as rks.ref's own CONV_TOL_i/CONV_TOL_GRAD_i/
            ! GRAD_NORM_i are: the pyscf convergence this energy was
            ! committed at.
            call ref_read_real(refu, "CONV_TOL_"//trim(tag), pyscf_conv_tol)
            call ref_read_real(refu, "CONV_TOL_GRAD_"//trim(tag), pyscf_conv_tol_grad)
            call ref_read_real(refu, "GRAD_NORM_"//trim(tag), pyscf_grad_norm)
            if (trim(name_ref) /= trim(name_this)) error stop "scf_df: pyscf reference case order does not " &
               //"match this build ("//trim(name_ref)//" vs "//trim(name_this)//")"
            de = abs(res%energy - e_ref)
            ok_case(ic) = res%converged .and. de < E_TOL
            write (msg_case(ic), "(a,es10.2)") trim(name_this)//": E vs pyscf, |dE| = ", de

            write (unit) int(bas%nao, int64), int(natm, int64), int(nalpha, int64), int(nbeta, int64)
            write (unit) GMOL(ig)
            write (unit) GBASIS(ig)
            write (unit) METHOD(im)
            write (unit) res%energy
         end do

         call dfjk%release()
         call auxb%release()
         call bas%release()
         call la%release()
         deallocate (zint, at_r)
      end do
      call golden_end()
      close (unit)
      call ref_close(refu)
   end subroutine run_all

end module scf_df_suite

program scf_df
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_df_suite, only: collect_scf_df
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_scf_df, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_df
