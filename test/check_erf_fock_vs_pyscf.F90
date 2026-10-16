!
! eri%fock / eri%fock_nosym with the range-separated (erf) operator against
! pyscf's own range-separated J/K, on a fixed, deterministic, non-SCF
! density -- the omega > 0 analogue of check_rys_fock_vs_pyscf.F90.
!
! eri%fock has no entry point that returns J and K separately: it always
! returns the combined g = j_scale*J - k_scale*K/2. Two calls at the same
! omega isolate them instead of adding a second interface just for this
! check: (j_scale=1, k_scale=0) leaves g = J, and (j_scale=0, k_scale=1)
! leaves g = -K/2, so K = -2*g. This is exact, not an approximation -- the
! kernels are linear in jfac/kfac by construction (see
! scripts/gen_perclass.py's digestion), and the same trick works whichever
! path g went through.
!
! The fourth argument picks eri%fock (0, the folded/symmetric path -- HGP by
! default, Rys when the caller sets TRC_RYS_CLASSES, see test/CMakeLists.txt)
! or eri%fock_nosym (1, the enumerated/general path trc_binkernel.F90's
! one_bin_item_nosym runs). Same AO-order and AO-normalisation story as
! check_rys_fock_vs_pyscf.F90 (and the same dump layout plus one more
! matrix); erf_fock_ref.py does the pyscf side.
!
program check_erf_fock_vs_pyscf
   use pic_types, only: dp, int64
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_tables, only: LMAX
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_report
   implicit none

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: OMEGA_TEST = 0.3_dp
   real(dp), parameter :: TOL = 1.0e-9_dp
   character(len=256) :: xyzfile, basis_json, dump, nosym_arg, case_name
   logical :: nosym
   type(trc_basis_t) :: bas
   type(trc_eri_t)   :: eri
   type(trc_pairlist_t) :: pl
   type(error_t) :: err
   integer :: nat, nao, unit
   integer, allocatable :: at_z(:)
   real(dp), allocatable :: r(:, :)
   real(dp), allocatable :: d(:, :), gj(:, :), gk(:, :)
   real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)

   if (command_argument_count() < 4) then
      print '(a)', 'usage: check_erf_fock_vs_pyscf <xyz> <basis.json> <dump.bin> <nosym 0|1>'
      stop 1
   end if
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, basis_json)
   call get_command_argument(3, dump)
   call get_command_argument(4, nosym_arg)
   nosym = (trim(nosym_arg) == '1')
   case_name = 'water_ccpvdz'
   if (command_argument_count() >= 5) call get_command_argument(5, case_name)

   call read_xyz(trim(xyzfile), nat, at_z, r)
   call trc_basis_from_json(trim(basis_json), nat, at_z, r, bas, err)
   if (err%has_error()) then
      print '(a)', 'check_erf_fock_vs_pyscf: '//err%get_message(); stop 1
   end if

   ! Same reasoning as check_rys_fock_vs_pyscf.F90: nothing to test yet on a
   ! build whose kernels stop short of this basis's highest shell.
   if (maxval(bas%sh_l) > LMAX) then
      print '(a,i0,a,i0,a)', '  SKIP: basis needs l = ', maxval(bas%sh_l), &
         ' and this build has LMAX = ', LMAX, '.'
      stop 0
   end if

   nao = bas%nao
   allocate (d(nao, nao), gj(nao, nao), gk(nao, nao))
   call synth_density(nao, d)

   call bas%to_device()
   call eri%build(bas, THRESH)
   if (nosym) then
      call eri%fock_nosym(bas, d, gj, density_screen=.false., &
                          j_scale=1.0_dp, k_scale=0.0_dp, omega=OMEGA_TEST)
      call eri%fock_nosym(bas, d, gk, density_screen=.false., &
                          j_scale=0.0_dp, k_scale=1.0_dp, omega=OMEGA_TEST)
   else
      call eri%fock(bas, d, gj, density_screen=.false., &
                    j_scale=1.0_dp, k_scale=0.0_dp, omega=OMEGA_TEST)
      call eri%fock(bas, d, gk, density_screen=.false., &
                    j_scale=0.0_dp, k_scale=1.0_dp, omega=OMEGA_TEST)
   end if
   call eri%release()

   ! K = -2*g from the second call; see this program's own header.
   gk = -2.0_dp*gk

   ! The overlap, for erf_fock_ref.py's per-AO normalisation check -- same
   ! reasoning as check_rys_fock_vs_pyscf.F90's own smat dump.
   allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
   call pl%build(bas, THRESH)
   call pl%to_device()
   !$acc enter data create(smat, tmat, vmat)
   call trc_1e(bas, pl, smat, tmat, vmat)
   !$acc update self(smat)
   !$acc exit data delete(smat, tmat, vmat)
   call pl%release()
   call bas%release()

   open (newunit=unit, file=trim(dump), access='stream', &
         form='unformatted', status='replace')
   write (unit) int(nat, int64), int(nao, int64)
   write (unit) int(at_z, int64)
   write (unit) r
   write (unit) d
   write (unit) gj
   write (unit) gk
   write (unit) smat
   close (unit)
   print '(a,i0,a,i0,a)', 'check_erf_fock_vs_pyscf: wrote '//trim(dump)//' (', &
      nat, ' atoms, ', nao, ' AOs) for the pyscf-reference comparison'

   ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
   ! One reference file per (geometry, basis), shared by every kernel path
   ! (hgp/nosym/rys) this program is built for -- see
   ! test/pyscf_refs/generate.py's gen_erf_fock for why.
   block
      integer :: u, nat_ref, nao_ref, nfail
      real(dp), allocatable :: s_ref(:, :), j_ref(:, :), k_ref(:, :)
      real(dp) :: ds, worst_j, worst_k
      call ref_require('erf_fock_'//trim(case_name)//'.ref', u)
      call ref_read_int(u, 'NAT', nat_ref)
      call ref_read_int(u, 'NAO', nao_ref)
      if (nat_ref /= nat .or. nao_ref /= nao) then
         print '(a)', 'check_erf_fock_vs_pyscf: pyscf reference shape does not match this build'; stop 1
      end if
      allocate (s_ref(nao, nao), j_ref(nao, nao), k_ref(nao, nao))
      call ref_read_matrix(u, 'S', nao, nao, s_ref)
      call ref_read_matrix(u, 'J', nao, nao, j_ref)
      call ref_read_matrix(u, 'K', nao, nao, k_ref)
      call ref_close(u)
      nfail = 0
      ds = maxval(abs(smat - s_ref))
      call ref_report(nfail, ds < TOL, 'overlap check vs pyscf, max|diff|', ds)
      worst_j = maxval(abs(gj - j_ref))
      worst_k = maxval(abs(gk - k_ref))
      call ref_report(nfail, worst_j < TOL, 'J(omega) vs pyscf, worst|diff|', worst_j)
      call ref_report(nfail, worst_k < TOL, 'K(omega) vs pyscf, worst|diff|', worst_k)
      if (nfail > 0) then
         print '(a)', 'check_erf_fock_vs_pyscf: FAILED against the pyscf reference'; stop 1
      end if
      print '(a)', 'check_erf_fock_vs_pyscf: RESULT: PASS'
   end block

   deallocate (at_z, r, d, gj, gk, smat, tmat, vmat)

contains

   !
   ! Same deterministic, dense, symmetric construction as
   ! check_rys_fock_vs_pyscf.F90.
   !
   subroutine synth_density(n, dd)
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

end program check_erf_fock_vs_pyscf
