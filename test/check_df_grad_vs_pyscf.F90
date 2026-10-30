!
! The FULL density-fitted (`--df`) RHF/RKS nuclear gradient
! (trc_scf_grad%trc_scf_gradient with `dfjk`/`aux`: fitted J and K, 1e, XC,
! nuclear repulsion) against pyscf's own density-fitted gradient
! (pyscf.df.grad.rhf / pyscf.df.grad.rks, auxiliary basis
! def2-universal-jkfit, mol.cart=True) AT TERCO'S OWN converged orbitals and,
! for a functional, TERCO'S OWN grid -- gen_dft_grad's recipe
! (test/pyscf_refs/generate.py's gen_df_grad), so the comparison is of the
! gradient formula at one fixed state, not of two SCFs.
!
! pyscf's DF gradient includes the auxiliary-basis response (both the
! (mu nu|P) aux-centre derivative and the metric's) by default; terco's
! always does. Grid response is OFF on both sides for this comparison
! (pyscf's default); check_df_grad_fd.F90 checks the grid-response-ON
! gradient against a finite difference instead.
!
! Writes the dump gen_df_grad reads: geometry, grid, orbitals, terco's
! gradient and the Cartesian overlap (for the per-AO normalisation).
!
! SKIPS CLEANLY (print SKIP, stop 0, before run_testsuite) when the
! auxiliary basis needs more angular momentum than this build's
! DF_LMAX_AUX, same convention as check_nosym.F90.
!
module check_df_grad_vs_pyscf_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_scf_df, only: df_jk_t
   use trc_df_kernels, only: DF_LMAX_AUX
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_linalg, only: trc_linalg_t
   use trc_test_basis, only: read_xyz
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix
   implicit none
   private
   public :: collect_df_grad_vs_pyscf_tests, check_df_grad_vs_pyscf_prepare

   real(dp), parameter :: THRESH = 1.0e-12_dp
   !! Per-case tolerance (`tol`, below) comes from test/CMakeLists.txt (the
   !! 8th argument): 1e-10 on water, 5e-9 on formamide. The formamide gap
   !! (9.4e-10 RHF, 6.2e-10 B3LYP) is the (P|Q) VALUES through an
   !! ill-conditioned fit, PROVEN, not assumed. With terco's converged
   !! orbitals fixed, the fitted gradient was re-evaluated in numpy from
   !! each code's raw tensors (terco's own (mn|P), (P|Q), d(mn|P)/dA,
   !! d(P|Q)/dA dumped from the same kernels this program runs),
   !! formamide/cc-pVDZ, Cartesian def2-universal-jkfit, cond(M) = 8.6e13:
   !!   - raw tensors, terco vs pyscf, max |diff| / max |element|: (mn|P)
   !!     5.2e-14, (P|Q) 2.8e-14, d(mn|P)/dA 2.5e-14, d(P|Q)/dA 5.6e-14;
   !!   - numpy with terco's tensors reproduces terco's gradient to 4.1e-12,
   !!     with pyscf's reproduces pyscf's df.grad to 4.3e-12 (so neither
   !!     contraction is off), and the two differ by 9.5e-10, the gap;
   !!   - terco's tensors with ONLY pyscf's (P|Q) values swapped in land
   !!     4.3e-11 from pyscf's gradient (9.1e-10 from terco's), and both
   !!     derivative tensors swapped together move it only 3.3e-11.
   !! So a 3e-12 (2.8e-14 relative) difference in the metric, amplified by
   !! cond(M), is the whole gap. Water (cond 4.6e11) agrees to 1.6e-13.
   !! Independently of pyscf, the same formamide gradient matches a
   !! Richardson finite difference of terco's reconverged --df energy to
   !! 7.1e-10 (check_df_grad_fd formamide).
   !! Revisit once spherical auxiliary functions land (PR #82): the
   !! Cartesian auxiliary set's non-spherical contaminants likely drive
   !! cond(M), and with them gone the formamide gate should tighten to
   !! water's.

   character(len=256) :: g_xyzfile, g_basis_json, g_aux_json, g_functional, g_dump, g_case_name
   integer :: g_grid_level
   real(dp) :: g_tol
   logical :: g_is_dft

contains

   subroutine check_df_grad_vs_pyscf_prepare(xyzfile, basis_json, aux_json, functional, grid_level, &
                                             dump, case_name, tol, skip)
      character(*), intent(in) :: xyzfile, basis_json, aux_json, functional, dump, case_name
      integer, intent(in) :: grid_level
      real(dp), intent(in) :: tol
      logical, intent(out) :: skip
      type(trc_basis_t) :: aux
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)

      g_xyzfile = xyzfile; g_basis_json = basis_json; g_aux_json = aux_json
      g_functional = functional; g_grid_level = grid_level; g_dump = dump
      g_case_name = case_name; g_tol = tol
      g_is_dft = trim(functional) /= 'hf'

      ! Only the aux basis's own ceiling decides the skip, and building it
      ! once here (then rebuilding inside the test) is cheap next to the
      ! SCF -- simpler than threading a cached basis through, and this
      ! program only ever builds one test.
      call read_xyz(trim(xyzfile), nat, at_z, r)
      call trc_basis_from_json(trim(aux_json), nat, at_z, r, aux, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
      end if
      skip = maxval(aux%sh_l) > DF_LMAX_AUX
      call aux%release()
   end subroutine check_df_grad_vs_pyscf_prepare

   subroutine collect_df_grad_vs_pyscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("--df gradient vs pyscf df.grad", test_df_grad_vs_pyscf)]
   end subroutine collect_df_grad_vs_pyscf_tests

   subroutine test_df_grad_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nat, nao, nalpha, nbeta, unit, dev, npts
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), grad(:, :)
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: lp
      type(df_jk_t) :: dfjk
      type(error_t) :: err
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      integer :: u, nat_ref, nao_ref
      real(dp), allocatable :: grad_ref(:, :)
      real(dp) :: worst
      logical :: shape_ok

      if (g_is_dft) then
         if (.not. xc_functional_by_name(trim(g_functional), func)) then
            print '(a)', 'check_df_grad_vs_pyscf: functional not available'; stop 1
         end if
      end if

      dev = trc_bind_device(0)
      call read_xyz(trim(g_xyzfile), nat, at_z, r)
      nalpha = sum(at_z)/2
      nbeta = nalpha
      call trc_basis_from_json(trim(g_basis_json), nat, at_z, r, bas, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
      end if
      call trc_basis_from_json(trim(g_aux_json), nat, at_z, r, aux, err, force_cart=.true.)
      if (err%has_error()) then
         print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
      end if
      nao = bas%nao
      call bas%to_device()
      call aux%to_device()
      call pl%build(bas, THRESH)
      call pl%to_device()
      call lp%init(nao)
      call dfjk%setup(bas, aux, pl, lp, err, nocc_hint=nalpha)
      if (err%has_error()) then
         print '(a)', 'check_df_grad_vs_pyscf: '//err%get_message(); stop 1
      end if

      opts%df = .true.
      opts%guess = "gwh"
      if (g_is_dft) opts%functional = trim(g_functional)
      opts%grid_level = g_grid_level
      opts%conv_energy = 1.0e-12_dp; opts%conv_diis = 1.0e-10_dp; opts%max_iteration = 200
      call trc_scf_run(bas, nalpha, nbeta, opts, res, la=lp, scf_dfjk=dfjk)
      call check(error, res%converged, "check_df_grad_vs_pyscf: reference SCF did not converge")
      if (allocated(error)) then
         call dfjk%release(); call pl%release(); call lp%release(); call aux%release(); call bas%release()
         return
      end if

      allocate (grad(3, nat))
      if (g_is_dft) then
         call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                               func=func, grid_level=g_grid_level, grid_max_pts=opts%grid_max_pts, &
                               rho_tol=opts%rho_tol, grid_response=.false., dfjk=dfjk, aux=aux)
      else
         call trc_scf_gradient(bas, lp, res%nspin, nalpha, nbeta, res%cmo, res%eps, THRESH, grad, &
                               dfjk=dfjk, aux=aux)
      end if

      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(smat)
      !$acc exit data delete(smat, tmat, vmat)
      call dfjk%release()
      call pl%release()
      call lp%release()
      call aux%release()
      call bas%release()

      npts = 0
      if (g_is_dft) npts = size(res%grid_r, 2)
      open (newunit=unit, file=trim(g_dump), access='stream', form='unformatted', status='replace')
      write (unit) int(nat, int64), int(nao, int64), int(nalpha, int64)
      write (unit) int(at_z, int64)
      write (unit) r
      write (unit) int(npts, int64)
      if (g_is_dft) then
         write (unit) res%grid_r
         write (unit) res%grid_w
      end if
      write (unit) res%cmo
      write (unit) res%eps
      write (unit) grad
      write (unit) smat
      close (unit)
      print '(a)', 'check_df_grad_vs_pyscf: wrote '//trim(g_dump)

      call ref_require('dfg_'//trim(g_case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      shape_ok = nat_ref == nat .and. nao_ref == nao
      call check(error, shape_ok, "check_df_grad_vs_pyscf: pyscf reference shape does not match this build")
      if (allocated(error)) then
         call ref_close(u)
         deallocate (at_z, r, grad, smat, tmat, vmat)
         return
      end if
      allocate (grad_ref(3, nat))
      call ref_read_matrix(u, 'GRAD', 3, nat, grad_ref)
      call ref_close(u)
      worst = maxval(abs(grad - grad_ref))
      print '(a,es10.2)', '  translational-invariance residual: ', maxval(abs(sum(grad, dim=2)))
      call check(error, worst < g_tol, "--df gradient vs pyscf df.grad, worst|diff| = "// &
                 to_string(worst, "es10.2")//" exceeds gate "//to_string(g_tol, "es10.2"))
      deallocate (at_z, r, grad, smat, tmat, vmat, grad_ref)
   end subroutine test_df_grad_vs_pyscf

end module check_df_grad_vs_pyscf_suite

program check_df_grad_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use pic_types, only: dp
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_df_grad_vs_pyscf_suite, only: collect_df_grad_vs_pyscf_tests, check_df_grad_vs_pyscf_prepare
   implicit none
   integer :: stat, grid_level
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=256) :: xyzfile, basis_json, aux_json, functional, dump, case_name, glevel_str, tol_str
   real(dp) :: tol
   logical :: skip

   if (command_argument_count() < 8) then
      print '(a)', 'usage: check_df_grad_vs_pyscf <xyz> <basis.json> <aux.json> <functional|hf> '// &
         '<grid_level> <dump.bin> <case_name> <tol>'; stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, aux_json)
   call get_command_argument(4, functional)
   call get_command_argument(5, glevel_str)
   call get_command_argument(6, dump)
   call get_command_argument(7, case_name)
   call get_command_argument(8, tol_str)
   read (glevel_str, *) grid_level
   read (tol_str, *) tol

   call check_df_grad_vs_pyscf_prepare(trim(xyzfile), trim(basis_json), trim(aux_json), trim(functional), &
                                       grid_level, trim(dump), trim(case_name), tol, skip)
   if (skip) then
      print '(a)', '  SKIP: the auxiliary basis needs more angular momentum than this build''s DF_LMAX_AUX.'
      stop 0
   end if

   stat = 0
   testsuites = [new_testsuite("check_df_grad_vs_pyscf", collect_df_grad_vs_pyscf_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_df_grad_vs_pyscf
