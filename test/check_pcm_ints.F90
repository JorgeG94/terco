!
! The PCM solute-surface integrals (trc_pcm_potential, trc_pcm_fock; the
! --pcm kernels of scripts/gen_1e.py) against pyscf, on fixed inputs.
!
! The points are water's own cavity at 302 per sphere plus eight more
! placed where the cavity never goes -- on and between the nuclei, with
! exponents from 0.3 to 300 -- so the comparison reaches the regimes a
! cavity does not: theta near one (a sharp charge on top of the density)
! and near zero, and T from zero up. The density is the fixed, dense,
! symmetric recipe check_rys_fock_vs_pyscf uses; the charges a fixed
! function of the point index. Nothing here comes from an SCF.
!
!   v(k)       = sum_{mu nu} D_mu,nu <mu| erf(zeta_k|r - C_k|)/|r - C_k| |nu>
!   F(mu, nu)  = sum_k q_k <mu| erf(zeta_k|r - C_k|)/|r - C_k| |nu>
!
! test/pyscf_refs/pcm_ints_<case>.ref (generate.py's gen_pcm_ints) holds
! pyscf's int3c2e against fakemol_for_charges -- the path pyscf's own PCM
! takes -- contracted the same two ways, after the generator has checked it
! against int1e_grids under set_range_coulomb(zeta), one zeta group at a
! time. Matrices already in this build's AO order and normalisation.
!
! SKIPS CLEANLY (prints SKIP, exits 0) when the basis needs more angular
! momentum than this build's PCM kernels reach, as check_rys_fock_vs_pyscf
! does, in run_case before run_testsuite is ever called -- past that guard
! a missing reference is a FAIL.
!
module check_pcm_ints_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_pcm_potential, trc_pcm_fock, PCM_LMAX
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_pcm, only: pcm_surface_t, pcm_build_surface
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector, ref_read_matrix, ref_report
   implicit none
   private
   public :: collect_pcm_ints_tests, run_case

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: TOL = 1.0e-10_dp
   integer, parameter :: NEXTRA = 8

   real(dp) :: g_s, g_v, g_v_off, g_f
      !! max|diff| against pyscf for the overlap (normalisation check), v_k
      !! at every point, v_k restricted to the eight off-cavity points, and
      !! F -- filled in by run_case.

contains

   subroutine collect_pcm_ints_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("overlap, v_k and F match pyscf", test_vs_pyscf) &
                  ]
   end subroutine collect_pcm_ints_tests

   !
   ! Everything the old program did: argument parsing, the PCM_LMAX self-
   ! skip, building the cavity plus the eight off-cavity points, the fixed
   ! density/charges, the potential/Fock/overlap computation, writing the
   ! probe dump, and reading the pyscf reference back into the four
   ! max|diff| values test_vs_pyscf checks.
   !
   subroutine run_case(skip, message)
      logical, intent(out) :: skip
      character(len=*), intent(out) :: message

      character(len=256) :: xyzfile, basis_json, dump, case_name
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(pcm_surface_t) :: surf
      type(error_t) :: err
      integer :: nat, nao, npts, i, k, unit
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :), pts(:, :), zeta(:), q(:), d(:, :), v(:), fmat(:, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      integer :: u, nao_ref, npts_ref, nfail
      real(dp), allocatable :: s_ref(:, :), f_ref(:, :), v_ref(:)

      skip = .false.
      message = ''

      if (command_argument_count() < 4) then
         print '(a)', 'usage: check_pcm_ints <xyz> <basis.json> <dump.bin> <case>'; stop 1
      end if
      call get_command_argument(1, xyzfile)
      call get_command_argument(2, basis_json)
      call get_command_argument(3, dump)
      call get_command_argument(4, case_name)

      call read_xyz(trim(xyzfile), nat, at_z, r)
      call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
      if (err%has_error()) then
         print '(a)', 'check_pcm_ints: '//err%get_message(); stop 1
      end if
      ! The kernels are Cartesian whatever the file says; so is the reference.
      bas%spherical = .false.
      if (maxval(bas%sh_l) > PCM_LMAX) then
         write (message, '(a,i0,a,i0,a)') '  SKIP: basis needs l = ', maxval(bas%sh_l), &
            ' and this build has PCM_LMAX = ', PCM_LMAX, '.'
         skip = .true.; return
      end if
      nao = bas%nao

      ! --- the points: the cavity, then the ones it never reaches ---
      call pcm_build_surface(at_z, r, 302, 1.2_dp, surf, err)
      if (err%has_error()) then
         print '(a)', 'check_pcm_ints: '//err%get_message(); stop 1
      end if
      npts = surf%n_points + NEXTRA
      allocate (pts(3, npts), zeta(npts), q(npts), v(npts))
      pts(:, 1:surf%n_points) = surf%coords
      zeta(1:surf%n_points) = surf%zeta
      do i = 1, NEXTRA
         k = surf%n_points + i
         ! on atom 1, on atom 2, between them, and off-centre, at widths
         ! spanning three decades
         select case (mod(i - 1, 4))
         case (0); pts(:, k) = r(:, 1)
         case (1); pts(:, k) = r(:, 2)
         case (2); pts(:, k) = 0.5_dp*(r(:, 1) + r(:, 2))
         case default; pts(:, k) = r(:, 1) + [0.31_dp, -0.17_dp, 0.23_dp]
         end select
         zeta(k) = 0.3_dp*10.0_dp**(real(i - 1, dp)*3.0_dp/real(NEXTRA - 1, dp))
      end do
      do k = 1, npts
         q(k) = 0.01_dp*cos(0.37_dp*real(k, dp)) + 0.002_dp
      end do

      allocate (d(nao, nao), fmat(nao, nao))
      call synth_density(nao, d)
      fmat = 0.0_dp

      call bas%to_device()
      call pl%build(bas, THRESH)
      call pl%to_device()
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      !$acc enter data copyin(pts, zeta, q, d, fmat) create(v, smat, tmat, vmat)
      call trc_pcm_potential(bas, pl, npts, pts, zeta, d, v)
      call trc_pcm_fock(bas, pl, npts, pts, zeta, q, fmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(v, fmat, smat)
      !$acc exit data delete(pts, zeta, q, d, fmat, v, smat, tmat, vmat)
      call pl%release()
      call bas%release()

      open (newunit=unit, file=trim(dump), access='stream', form='unformatted', status='replace')
      write (unit) int(nat, int64), int(nao, int64), int(npts, int64)
      write (unit) int(at_z, int64)
      write (unit) r
      write (unit) pts, zeta, q
      write (unit) d
      write (unit) smat
      close (unit)
      print '(a,i0,a,i0,a)', 'check_pcm_ints: wrote '//trim(dump)//' (', nao, ' AOs, ', npts, ' points)'

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      call ref_require('pcm_ints_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAO', nao_ref)
      call ref_read_int(u, 'NPTS', npts_ref)
      if (nao_ref /= nao .or. npts_ref /= npts) then
         print '(a)', 'check_pcm_ints: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (s_ref(nao, nao), f_ref(nao, nao), v_ref(npts))
      call ref_read_matrix(u, 'S', nao, nao, s_ref)
      call ref_read_vector(u, 'V', npts, v_ref)
      call ref_read_matrix(u, 'F', nao, nao, f_ref)
      call ref_close(u)
      nfail = 0
      g_s = maxval(abs(smat - s_ref))
      call ref_report(nfail, g_s < TOL, 'overlap check vs pyscf (normalisation)', g_s)
      g_v = maxval(abs(v - v_ref))
      call ref_report(nfail, g_v < TOL, 'v_k = D . <|erf(zeta_k r)/r|> at every point', g_v)
      g_v_off = maxval(abs(v(npts - NEXTRA + 1:) - v_ref(npts - NEXTRA + 1:)))
      call ref_report(nfail, g_v_off < TOL, '  of which the eight off-cavity points', g_v_off)
      g_f = maxval(abs(fmat - f_ref))
      call ref_report(nfail, g_f < TOL, 'F = sum_k q_k <|erf(zeta_k r)/r|>', g_f)

      deallocate (pts, zeta, q, v, d, fmat, smat, tmat, vmat, s_ref, f_ref, v_ref)
   end subroutine run_case

   subroutine synth_density(n, dd)
      !! check_rys_fock_vs_pyscf's: dense, symmetric, not an SCF density
      integer, intent(in) :: n
      real(dp), intent(out) :: dd(n, n)
      integer :: ii, jj
      do jj = 1, n
         do ii = 1, n
            dd(ii, jj) = 0.4_dp + 0.1_dp*cos(real(3*ii + 5*jj, dp))
         end do
      end do
      dd = 0.5_dp*(dd + transpose(dd))
   end subroutine synth_density

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, g_s < TOL, "overlap check vs pyscf (normalisation), max|diff| = "// &
                 to_string(g_s, "es12.4")//" exceeds "//to_string(TOL, "es12.4"))
      if (allocated(error)) return
      call check(error, g_v < TOL, "v_k = D . <|erf(zeta_k r)/r|> at every point, max|diff| = "// &
                 to_string(g_v, "es12.4")//" exceeds "//to_string(TOL, "es12.4"))
      if (allocated(error)) return
      call check(error, g_v_off < TOL, "v_k at the eight off-cavity points, max|diff| = "// &
                 to_string(g_v_off, "es12.4")//" exceeds "//to_string(TOL, "es12.4"))
      if (allocated(error)) return
      call check(error, g_f < TOL, "F = sum_k q_k <|erf(zeta_k r)/r|>, max|diff| = "// &
                 to_string(g_f, "es12.4")//" exceeds "//to_string(TOL, "es12.4"))
   end subroutine test_vs_pyscf

end module check_pcm_ints_suite

program check_pcm_ints
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_pcm_ints_suite, only: collect_pcm_ints_tests, run_case
   implicit none
   integer :: stat
   logical :: skip
   character(len=256) :: message
   type(testsuite_type), allocatable :: testsuites(:)

   call run_case(skip, message)
   if (skip) then
      print '(a)', trim(message)
      stop 0
   end if

   stat = 0
   testsuites = [new_testsuite("check_pcm_ints", collect_pcm_ints_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_pcm_ints
