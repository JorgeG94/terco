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
! does; past that guard a missing reference is a FAIL.
!
program check_pcm_ints
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_pcm_potential, trc_pcm_fock, PCM_LMAX
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_pcm, only: pcm_surface_t, build_pcm_surface
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector, ref_read_matrix, ref_report
   implicit none

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: TOL = 1.0e-10_dp
   integer, parameter :: NEXTRA = 8
   character(len=256) :: xyzfile, basis_json, dump, case_name
   type(trc_basis_t) :: bas
   type(trc_pairlist_t) :: pl
   type(pcm_surface_t) :: surf
   type(error_t) :: err
   integer :: nat, nao, npts, i, k, unit
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :), pts(:, :), zeta(:), q(:), d(:, :), v(:), fmat(:, :)
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)

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
      print '(a,i0,a,i0,a)', '  SKIP: basis needs l = ', maxval(bas%sh_l), &
         ' and this build has PCM_LMAX = ', PCM_LMAX, '.'
      stop 0
   end if
   nao = bas%nao

   ! --- the points: the cavity, then the ones it never reaches ---
   call build_pcm_surface(at_z, r, 302, 1.2_dp, surf, err)
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
   block
      integer :: u, nao_ref, npts_ref, nfail
      real(dp), allocatable :: s_ref(:, :), f_ref(:, :), v_ref(:)
      real(dp) :: worst
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
      worst = maxval(abs(smat - s_ref))
      call ref_report(nfail, worst < TOL, 'overlap check vs pyscf (normalisation)', worst)
      worst = maxval(abs(v - v_ref))
      call ref_report(nfail, worst < TOL, 'v_k = D . <|erf(zeta_k r)/r|> at every point', worst)
      worst = maxval(abs(v(npts - NEXTRA + 1:) - v_ref(npts - NEXTRA + 1:)))
      call ref_report(nfail, worst < TOL, '  of which the eight off-cavity points', worst)
      worst = maxval(abs(fmat - f_ref))
      call ref_report(nfail, worst < TOL, 'F = sum_k q_k <|erf(zeta_k r)/r|>', worst)
      if (nfail > 0) then
         print '(a)', 'check_pcm_ints: FAILED against the pyscf reference'; stop 1
      end if
      print '(a)', 'check_pcm_ints: RESULT: PASS'
   end block

contains

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

end program check_pcm_ints
