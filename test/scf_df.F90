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
! on O/C/N. Every case prints SKIP and this exits 0 rather than failing
! when `maxval(aux%sh_l) > DF_LMAX_AUX` -- a build below that ceiling
! self-skips cleanly rather than failing.
!
program scf_df
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_df, only: df_jk_t
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, ref_report
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none

   integer, parameter :: NGROUP = 4, NMETHOD = 3
   character(len=16), parameter :: GMOL(NGROUP) = [character(len=16) :: &
                                                   "water", "water", "formamide", "formamide"]
   character(len=16), parameter :: GBASIS(NGROUP) = [character(len=16) :: &
                                                      "cc-pvdz", "6-31gs", "cc-pvdz", "6-31gs"]
   character(len=8), parameter :: METHOD(NMETHOD) = [character(len=8) :: "", "pbe", "b3lyp"]
   real(dp), parameter :: E_TOL = 1.0e-8_dp   !! rks.ref's closed-shell tolerance; every case here is nspin=1

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
   logical :: all_ok
   integer :: refu, ncase_ref, nfail, ic
   logical :: has_ref
   character(len=32) :: tag, name_ref, name_this
   real(dp) :: e_ref, de, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm

   call ref_open("scf_df.ref", refu, has_ref)
   call golden_begin("scf_df")
   nfail = 0
   if (has_ref) then
      call ref_read_int(refu, "NCASE", ncase_ref)
      if (ncase_ref /= NGROUP*NMETHOD) then
         print "(a)", "scf_df: pyscf reference case count does not match this build"
         stop 1
      end if
   end if

   open (newunit=unit, file="scf_df_probe.bin", access="stream", form="unformatted", status="replace", &
         action="write")
   write (unit) int(NGROUP*NMETHOD, int64)
   all_ok = .true.

   do ig = 1, NGROUP
      call read_xyz(trim(GMOL(ig))//".xyz", natm, zint, at_r)
      nelec = sum(zint)
      nalpha = (nelec + mod(nelec, 2))/2
      nbeta = nelec - nalpha

      call trc_basis_from_json("basis_sets/"//trim(GBASIS(ig))//".json", natm, zint, at_r, bas, err, &
                               force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_df: "//err%get_message()
         stop 1
      end if
      call bas%to_device()
      call trc_basis_from_json("basis_sets/def2-universal-jkfit.json", natm, zint, at_r, auxb, err, &
                               force_cart=.true.)
      if (err%has_error()) then
         print "(a)", "scf_df: "//err%get_message()
         stop 1
      end if
      if (maxval(auxb%sh_l) > DF_LMAX_AUX) then
         print "(a,i0,a,i0,a)", "  SKIP: def2-universal-jkfit needs l = ", maxval(auxb%sh_l), &
            " and this build has DF_LMAX_AUX = ", DF_LMAX_AUX, "."
         call auxb%release()
         call bas%release()
         if (has_ref) call ref_close(refu)
         close (unit)
         stop 0
      end if
      call auxb%to_device()

      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
      call la%init(bas%nao)
      call dfjk%setup(bas, auxb, pl, la, err, nocc_hint=max(nalpha, 1))
      if (err%has_error()) then
         print "(a)", "scf_df: "//err%get_message()
         stop 1
      end if
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
         call trc_scf_run(bas, nalpha, nbeta, opts, res, la=la, dfjk=dfjk)
         if (.not. res%converged) then
            print "(a,a)", "  NOT CONVERGED: ", trim(res%message)
            all_ok = .false.
         end if
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
         call golden_put("E_"//trim(name_this), res%energy)
         if (has_ref) then
            write (tag, "(i0)") ic - 1
            call ref_read_str(refu, "NAME_"//trim(tag), name_ref)
            call ref_read_real(refu, "E_"//trim(tag), e_ref)
            ! Provenance only, as rks.ref's own CONV_TOL_i/CONV_TOL_GRAD_i/
            ! GRAD_NORM_i are: the pyscf convergence this energy was
            ! committed at.
            call ref_read_real(refu, "CONV_TOL_"//trim(tag), pyscf_conv_tol)
            call ref_read_real(refu, "CONV_TOL_GRAD_"//trim(tag), pyscf_conv_tol_grad)
            call ref_read_real(refu, "GRAD_NORM_"//trim(tag), pyscf_grad_norm)
            if (trim(name_ref) /= trim(name_this)) then
               print "(a)", "scf_df: pyscf reference case order does not match this build "// &
                  "("//trim(name_ref)//" vs "//trim(name_this)//")"
               stop 1
            end if
            de = abs(res%energy - e_ref)
            call ref_report(nfail, res%converged .and. de < E_TOL, trim(name_this)//": E vs pyscf", de)
         end if

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
   if (has_ref) call ref_close(refu)
   if (.not. has_ref) then
      print "(a)", "FAIL: missing reference pyscf_refs/scf_df.ref -- run test/pyscf_refs/generate.py "// &
         "(see test/pyscf_refs/README.md)"
      stop 1
   end if
   if (.not. all_ok .or. nfail > 0) stop 1
   print "(a)", ""
   print "(a)", "scf_df: every case converged and agrees with the pyscf reference; RESULT: PASS"

end program scf_df
