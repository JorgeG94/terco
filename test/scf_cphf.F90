!
! RHF static dipole polarizability through trc_cphf, against pyscf's
! pyscf.prop.polarizability.rhf (test/pyscf_refs/generate.py's gen_cphf).
!
! Converges the SCF independently and tightly (conv_energy 1e-12,
! conv_diis 1e-10) rather than reading back a committed pyscf density
! (scf_rks.F90's "at-density" pattern): water and formamide have no
! near-degenerate frontier orbitals or symmetry-breaking ambiguity, so two
! independent tight RHF solves land on the same stationary point to well
! past CPHF's own 1e-6 a.u. gate -- the same assumption scf_rimp2.F90/
! scf_df.F90 already make comparing SCF energies to pyscf.
!
! Builds its OWN pair list, ERI engine and linalg handles: trc_scf_run
! releases everything it built internally once it returns, and trc_cphf
! needs `eri%fock_resident` (device-resident, density_screen=.false. --
! see trc_cphf.F90's own header) and the dipole integrals, neither of
! which the converged trc_scf_result_t carries.
!
! Arguments: xyz [basis.json [molname [basname [probe]]]]. Defaults water/
! 6-31G*. `molname`/`basname` are what generate.py rebuilds the SAME
! pyscf molecule from (bse_json_basis on basis_sets/<basname>.json) --
! plain tags, not filenames, so they travel through the probe as fixed-
! width strings with no path to re-derive. ONE (molecule, basis) case per
! process (four ctest targets, same source, different ARGS), so this
! suite has a single test-drive case.
!
module scf_cphf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_scf_polarizability, trc_cphf_solve_rhf
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix
   implicit none
   private
   public :: collect_scf_cphf, xyz_arg, bas_json_arg, molname_arg, basname_arg, probe_arg

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = 'basis_sets/6-31gs.json'
   character(len=16) :: molname_arg = 'water'
   character(len=16) :: basname_arg = '6-31gs'
   character(len=256) :: probe_arg = 'cphf_probe.bin'

contains

   subroutine collect_scf_cphf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF polarizability vs pyscf", test_cphf) &
                  ]
   end subroutine collect_scf_cphf

   subroutine test_cphf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm, nocc, unit, dev, x, y
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_cphf_options_t) :: cphf_opts
      type(trc_cphf_result_t) :: cphf_res
      type(error_t) :: err
      real(dp) :: alpha(3, 3)
      real(dp) :: worst_alpha

      dev = trc_bind_device(0)

      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      call trc_basis_from_json(trim(bas_json_arg), natm, zint, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call pl%build(bas, ERI_THRESH)
      call pl%to_device()

      nocc = nint(sum(real(zint, dp)))/2
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-10_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      block
         character(len=8) :: envv
         integer :: st
         call get_environment_variable("TERCO_VERBOSE", envv, status=st)
         opts%verbose = st == 0
         call get_environment_variable("TERCO_CPHF_VERBOSE", envv, status=st)
         cphf_opts%verbose = st == 0
      end block
      call trc_scf_run(bas, nocc, nocc, opts, res)
      if (.not. res%converged) print '(a)', "scf_cphf: "//trim(res%message)

      call eri%build(bas, ERI_THRESH)
      call la%init(bas%nao)
      cphf_opts%conv_residual = 1.0e-8_dp
      cphf_opts%max_iteration = 100

      block
         type(timer_type) :: tm
         real(dp) :: t_cphf
         integer :: niter_tot
         call tm%start()
         call trc_scf_polarizability(bas, pl, eri, nocc, nocc, res, alpha, cphf_opts, cphf_res, la)
         t_cphf = tm%get_elapsed_time()
         niter_tot = sum(cphf_res%iterations)
         print '(a,f10.3,a,i0,a,f10.4,a)', "  CPHF wall time ", t_cphf, " s over ", niter_tot, &
            " total CG iterations: ", t_cphf/real(max(niter_tot, 1), dp), " s/iteration"
      end block
      if (.not. cphf_res%converged) print '(a)', "scf_cphf: "//trim(cphf_res%message)

      print '(a,i0,a,i0)', "scf_cphf: nao ", bas%nao, "  nocc ", nocc
      print '(a,f20.12,a,i0)', "  E_RHF = ", res%energy, "   SCF iter ", res%iterations
      print '(a,3(i0,1x))', "  CPHF iterations (x,y,z RHS): ", cphf_res%iterations
      print '(a)', "  polarizability tensor (a.u.):"
      do x = 1, 3
         print '(3f16.8)', (alpha(x, y), y=1, 3)
      end do

      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(bas%nao, int64), int(natm, int64), int(nocc, int64)
      write (unit) molname_arg, basname_arg
      write (unit) res%energy
      write (unit) alpha
      close (unit)

      call check(error, res%converged, "RHF SCF did not converge")
      if (allocated(error)) then
         call eri%release(); call la%release(); call pl%release(); call bas%release()
         return
      end if
      call check(error, cphf_res%converged, "CPHF solve did not converge")
      if (allocated(error)) then
         call eri%release(); call la%release(); call pl%release(); call bas%release()
         return
      end if

      block
         integer :: u
         real(dp) :: e_scf_ref, alpha_ref(3, 3)
         character(len=256) :: refname
         refname = 'cphf_'//trim(molname_arg)//'_'//trim(basname_arg)//'.ref'
         call ref_require(trim(refname), u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
         call ref_close(u)

         call check(error, abs(res%energy - e_scf_ref) < 1.0e-9_dp, &
                    "E_RHF vs pyscf, |dE| = "//to_string(abs(res%energy - e_scf_ref), "es10.2"))
         if (allocated(error)) then
            call eri%release(); call la%release(); call pl%release(); call bas%release()
            return
         end if
         worst_alpha = maxval(abs(alpha - alpha_ref))
         call check(error, all(abs(alpha - alpha_ref) < 1.0e-6_dp), &
                    "alpha vs pyscf, worst|diff| = "//to_string(worst_alpha, "es10.2"))
         if (allocated(error)) then
            call eri%release(); call la%release(); call pl%release(); call bas%release()
            return
         end if
      end block

      ! A NaN right-hand side must report NOT converged, never "converged"
      ! with zero iterations -- see trc_cphf_solve_rhf's own comment on why
      ! its CG loop condition reads a NaN resid as "not yet converged".
      block
         use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
         real(dp) :: cmo_nan(bas%nao, bas%nao)
         real(dp), allocatable :: bvec_nan(:, :, :), u_nan(:, :, :)
         type(trc_cphf_result_t) :: res_nan
         integer :: nvirt_t
         real(dp) :: nan_val

         nan_val = ieee_value(nan_val, ieee_quiet_nan)
         nvirt_t = bas%nao - nocc
         cmo_nan = res%cmo(:, :, 1)
         allocate (bvec_nan(nvirt_t, nocc, 1), u_nan(nvirt_t, nocc, 1))
         bvec_nan = nan_val
         !$acc enter data copyin(cmo_nan, bvec_nan) create(u_nan)
         call trc_cphf_solve_rhf(bas, eri, la, bas%nao, nocc, nvirt_t, cmo_nan, res%eps(:, 1), &
                                 1, bvec_nan, u_nan, cphf_opts, res_nan)
         !$acc exit data delete(cmo_nan, bvec_nan, u_nan)
         if (.not. res_nan%converged) print '(a)', "scf_cphf: NaN right-hand side correctly refused: "// &
            trim(res_nan%message)
         deallocate (bvec_nan, u_nan)

         call eri%release(); call la%release(); call pl%release(); call bas%release()
         call check(error, .not. res_nan%converged, "a NaN right-hand side must not report converged")
      end block
   end subroutine test_cphf

end module scf_cphf_suite

program scf_cphf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_cphf_suite, only: collect_scf_cphf, xyz_arg, bas_json_arg, molname_arg, basname_arg, probe_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, molname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, basname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, probe_arg)

   stat = 0
   call run_testsuite(collect_scf_cphf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_cphf
