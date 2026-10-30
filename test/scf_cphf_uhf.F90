!
! UHF static dipole polarizability through trc_cphf (the coupled spin
! channels, trc_cphf_solve_uhf), against pyscf's
! pyscf.prop.polarizability.uhf -- test/pyscf_refs/generate.py's
! gen_cphf_uhf. OH radical (9 electrons, nalpha=5, nbeta=4).
!
! COMPARES EIGENVALUES OF alpha, NOT ITS RAW COMPONENTS. OH's unpaired
! electron sits in the (in exact C-infinity-v symmetry, degenerate) pi*
! pair; UHF breaks that degeneracy by picking ONE real combination of the
! two perpendicular p-like orbitals, and which combination it picks is a
! numerical accident of the iteration path, not something either code
! controls. Measured: two independent tight UHF solves (terco's and
! pyscf's, same geometry/basis, energies agreeing to 1e-11) landed on
! DIFFERENT rotations within that plane -- alpha_xx/alpha_yy/alpha_xy
! disagreed by up to 0.17 a.u. while alpha_xx+alpha_yy (the trace) and
! alpha_zz (along the bond, not degenerate) agreed to 1e-6. The physically
! meaningful, rotation-invariant comparison is alpha's eigenvalues, which
! this module computes (LAPACK dsyev on the host, via pic_syev -- a
! plain 3x3, nothing device-resident) and writes instead of trusting the
! raw tensor to line up component-by-component. A molecule with no such
! orbital degeneracy (every RHF case in scf_cphf.F90) has no such freedom,
! which is why that suite's own cases compare the raw 3x3 tensor directly.
!
! Same pattern as scf_cphf.F90 otherwise: its own pair list, ERI engine
! and linalg handles, built and released here rather than reused from the
! (already-released) SCF context.
!
module scf_cphf_uhf_suite
   use pic_types, only: dp, int64, default_int
   use pic_strings, only: to_string
   use pic_lapack_interfaces, only: pic_syev
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_cphf, only: trc_cphf_options_t, trc_cphf_result_t, trc_scf_polarizability, trc_cphf_solve_uhf
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix, ref_read_vector
   implicit none
   private
   public :: collect_scf_cphf_uhf, xyz_arg, bas_json_arg, molname_arg, basname_arg, probe_arg

   real(dp), parameter :: ERI_THRESH = 1.0e-12_dp

   character(len=256) :: xyz_arg = 'oh.xyz'
   character(len=256) :: bas_json_arg = 'basis_sets/6-31gs.json'
   character(len=16) :: molname_arg = 'oh'
   character(len=16) :: basname_arg = '6-31gs'
   character(len=256) :: probe_arg = 'cphf_uhf_probe.bin'

contains

   subroutine collect_scf_cphf_uhf(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("UHF polarizability vs pyscf", test_cphf_uhf) &
                  ]
   end subroutine collect_scf_cphf_uhf

   subroutine test_cphf_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm, nalpha, nbeta, unit, dev, x, y
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
      real(dp) :: alpha(3, 3), alpha_eig(3)
      integer(default_int) :: eig_info

      dev = trc_bind_device(0)

      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      call trc_basis_from_json(trim(bas_json_arg), natm, zint, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call pl%build(bas, ERI_THRESH)
      call pl%to_device()

      ! OH, neutral, doublet: 9 electrons, nalpha - nbeta = 1.
      nalpha = (nint(sum(real(zint, dp))) + 1)/2
      nbeta = nint(sum(real(zint, dp))) - nalpha

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
      call trc_scf_run(bas, nalpha, nbeta, opts, res)
      if (.not. res%converged) print '(a)', "scf_cphf_uhf: "//trim(res%message)

      call eri%build(bas, ERI_THRESH)
      call la%init(bas%nao)
      cphf_opts%conv_residual = 1.0e-8_dp
      cphf_opts%max_iteration = 100

      call trc_scf_polarizability(bas, pl, eri, nalpha, nbeta, res, alpha, cphf_opts, cphf_res, la)
      if (.not. cphf_res%converged) print '(a)', "scf_cphf_uhf: "//trim(cphf_res%message)

      print '(a,i0,a,i0,a,i0)', "scf_cphf_uhf: nao ", bas%nao, "  nalpha ", nalpha, "  nbeta ", nbeta
      print '(a,f20.12,a,i0)', "  E_UHF = ", res%energy, "   SCF iter ", res%iterations
      print '(a,3(i0,1x))', "  CPHF iterations (x,y,z RHS): ", cphf_res%iterations
      print '(a)', "  polarizability tensor (a.u.):"
      do x = 1, 3
         print '(3f16.8)', (alpha(x, y), y=1, 3)
      end do

      ! Eigenvalues of the symmetric alpha, ascending (LAPACK dsyev's own
      ! order) -- the rotation-invariant quantities this module's own
      ! header explains are the fair comparison for OH's degenerate pi*
      ! open shell. `alpha` is consumed (overwritten with eigenvectors, of
      ! no further use here); work from a copy so the printed/written tensor
      ! above stays what trc_scf_polarizability actually returned.
      block
         real(dp) :: alpha_copy(3, 3)
         alpha_copy = alpha
         call pic_syev(alpha_copy, alpha_eig, "N", "U", eig_info)
      end block
      print '(a,3f16.8)', "  alpha eigenvalues: ", alpha_eig

      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(bas%nao, int64), int(natm, int64), int(nalpha, int64), int(nbeta, int64)
      write (unit) molname_arg, basname_arg
      write (unit) res%energy
      write (unit) alpha
      write (unit) alpha_eig
      close (unit)

      call check(error, res%converged, "UHF SCF did not converge")
      if (allocated(error)) then
         call eri%release(); call la%release(); call pl%release(); call bas%release()
         return
      end if
      call check(error, res%nspin == 2, "expected an unrestricted run (nalpha /= nbeta)")
      if (allocated(error)) then
         call eri%release(); call la%release(); call pl%release(); call bas%release()
         return
      end if
      call check(error, cphf_res%converged, "CPHF solve did not converge")
      if (allocated(error)) then
         call eri%release(); call la%release(); call pl%release(); call bas%release()
         return
      end if
      call check(error, eig_info == 0, "pic_syev failed on the polarizability tensor")
      if (allocated(error)) then
         call eri%release(); call la%release(); call pl%release(); call bas%release()
         return
      end if

      block
         real(dp) :: e_scf_ref, alpha_ref(3, 3), eig_ref(3), worst_eig, trace_diff, zz_diff
         integer :: u
         character(len=256) :: refname
         refname = 'cphf_uhf_'//trim(molname_arg)//'_'//trim(basname_arg)//'.ref'
         call ref_require(trim(refname), u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_matrix(u, 'ALPHA', 3, 3, alpha_ref)
         call ref_read_vector(u, 'ALPHA_EIG', 3, eig_ref)
         call ref_close(u)

         call check(error, abs(res%energy - e_scf_ref) < 1.0e-9_dp, &
                    "E_UHF vs pyscf, |dE| = "//to_string(abs(res%energy - e_scf_ref), "es10.2"))
         if (allocated(error)) then
            call eri%release(); call la%release(); call pl%release(); call bas%release()
            return
         end if

         worst_eig = maxval(abs(alpha_eig - eig_ref))
         call check(error, all(abs(alpha_eig - eig_ref) < 1.0e-6_dp), &
                    "alpha eigenvalue vs pyscf, worst|diff| = "//to_string(worst_eig, "es10.2"))
         if (allocated(error)) then
            call eri%release(); call la%release(); call pl%release(); call bas%release()
            return
         end if

         ! alpha_ref is read above but was never otherwise checked: trace(xx+yy)
         ! and zz alone, not every raw element, because those two combinations
         ! are the ones that stay well-defined under OH's degenerate pi*
         ! mixing (see this module's own header) the way the eigenvalues
         ! above already are. 1e-6, same as the eigenvalues: measured 4.3e-7
         ! (trace) and 9.3e-8 (zz), both already under it.
         trace_diff = abs((alpha(1, 1) + alpha(2, 2)) - (alpha_ref(1, 1) + alpha_ref(2, 2)))
         call check(error, trace_diff < 1.0e-6_dp, &
                    "alpha xx+yy (trace) vs pyscf, |diff| = "//to_string(trace_diff, "es10.2"))
         if (allocated(error)) then
            call eri%release(); call la%release(); call pl%release(); call bas%release()
            return
         end if
         zz_diff = abs(alpha(3, 3) - alpha_ref(3, 3))
         call check(error, zz_diff < 1.0e-6_dp, "alpha zz vs pyscf, |diff| = "//to_string(zz_diff, "es10.2"))
         if (allocated(error)) then
            call eri%release(); call la%release(); call pl%release(); call bas%release()
            return
         end if
      end block

      ! A NaN right-hand side must report NOT converged, never "converged"
      ! with zero iterations -- same reason as trc_cphf_solve_rhf's own
      ! comment, for trc_cphf_solve_uhf's CG loop condition.
      block
         use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
         real(dp) :: cmo_a_nan(bas%nao, bas%nao), cmo_b_nan(bas%nao, bas%nao)
         real(dp), allocatable :: bvec_a_nan(:, :, :), bvec_b_nan(:, :, :)
         real(dp), allocatable :: u_a_nan(:, :, :), u_b_nan(:, :, :)
         type(trc_cphf_result_t) :: res_nan
         integer :: nva_t, nvb_t
         real(dp) :: nan_val

         nan_val = ieee_value(nan_val, ieee_quiet_nan)
         nva_t = bas%nao - nalpha
         nvb_t = bas%nao - nbeta
         cmo_a_nan = res%cmo(:, :, 1)
         cmo_b_nan = res%cmo(:, :, 2)
         allocate (bvec_a_nan(nva_t, nalpha, 1), bvec_b_nan(nvb_t, nbeta, 1))
         allocate (u_a_nan(nva_t, nalpha, 1), u_b_nan(nvb_t, nbeta, 1))
         bvec_a_nan = nan_val
         bvec_b_nan = 0.0_dp
         !$acc enter data copyin(cmo_a_nan, cmo_b_nan, bvec_a_nan, bvec_b_nan) create(u_a_nan, u_b_nan)
         call trc_cphf_solve_uhf(bas, eri, la, bas%nao, nalpha, nva_t, nbeta, nvb_t, cmo_a_nan, cmo_b_nan, &
                                 res%eps(:, 1), res%eps(:, 2), 1, bvec_a_nan, bvec_b_nan, u_a_nan, u_b_nan, &
                                 cphf_opts, res_nan)
         !$acc exit data delete(cmo_a_nan, cmo_b_nan, bvec_a_nan, bvec_b_nan, u_a_nan, u_b_nan)
         if (.not. res_nan%converged) print '(a)', "scf_cphf_uhf: NaN right-hand side correctly refused: "// &
            trim(res_nan%message)
         deallocate (bvec_a_nan, bvec_b_nan, u_a_nan, u_b_nan)

         call eri%release(); call la%release(); call pl%release(); call bas%release()
         call check(error, .not. res_nan%converged, "a NaN right-hand side must not report converged")
      end block
   end subroutine test_cphf_uhf

end module scf_cphf_uhf_suite

program scf_cphf_uhf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_cphf_uhf_suite, only: collect_scf_cphf_uhf, xyz_arg, bas_json_arg, molname_arg, basname_arg, probe_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, molname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, basname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, probe_arg)

   stat = 0
   call run_testsuite(collect_scf_cphf_uhf, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_cphf_uhf
