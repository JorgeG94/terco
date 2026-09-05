!
! RI-MP2 on water: RHF through the driver, then the correlation energy
! through trc_rimp2 with a synthetic auxiliary basis, once as a single
! auxiliary block and once in blocks of 20 functions (bit-identical
! required), on every rank of the communicator when there is one. Writes
! `rimp2_probe.bin` for rimp2_ref.py, which runs pyscf's RHF and DF-MP2 on
! the same orbital and auxiliary sets and compares E_os, E_ss and E_corr.
!
program scf_rimp2
   use trc_boys, only: dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t
   use trc_test_basis, only: read_xyz, build_631g, build_aux
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none

   integer :: natm, nsh, maxnp, nsha, maxnpa, unit, ish, dev, nocc
   integer, allocatable :: zint(:), sh_l(:), sh_np(:), ax_l(:), ax_np(:)
   real(dp), allocatable :: at_r(:, :), at_z(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(dp), allocatable :: ax_e(:, :), ax_c(:, :), ax_r(:, :)
   type(trc_basis_t) :: bas, aux
   type(trc_pairlist_t) :: pl
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   type(trc_rimp2_result_t) :: mp, mp2
   type(comm_t) :: comm
   type(error_t) :: err
   logical :: ok
   character(len=256) :: xyz, bas_json, aux_json, probe

   call pic_mpi_init()
   comm = comm_world()
   dev = trc_bind_device(comm%rank())

   ! Arguments: xyz [basis.json aux.json [probe]]. Without files, primitive
   ! 6-31G and the synthetic auxiliary ladder; with them, both sets from
   ! MolSSI BSE JSON through trc_basis_json, which is how mqc's sets arrive.
   xyz = 'water.xyz'; bas_json = ''; aux_json = ''; probe = 'rimp2_probe.bin'
   if (command_argument_count() >= 1) call get_command_argument(1, xyz)
   if (command_argument_count() >= 3) then
      call get_command_argument(2, bas_json)
      call get_command_argument(3, aux_json)
   end if
   if (command_argument_count() >= 4) call get_command_argument(4, probe)

   call read_xyz(trim(xyz), natm, zint, at_r)
   allocate (at_z(natm)); at_z = real(zint, dp)
   if (len_trim(bas_json) > 0) then
      call trc_basis_from_json(trim(bas_json), natm, zint, at_r, bas, err)
      if (.not. err%has_error()) call trc_basis_from_json(trim(aux_json), natm, zint, at_r, aux, err)
      if (err%has_error()) then
         print '(a)', "scf_rimp2: "//err%get_message()
         stop 1
      end if
      call unpack_basis(bas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      call unpack_basis(aux, nsha, ax_l, ax_np, ax_e, ax_c, ax_r, maxnpa)
   else
      call build_631g(natm, zint, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp, uncontracted=.true.)
      call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
      call build_aux(natm, zint, at_r, 2, nsha, ax_l, ax_np, ax_e, ax_c, ax_r, maxnpa)
      call aux%build(nsha, ax_l, ax_np, ax_e, ax_c, ax_r, natm, at_z, at_r, maxnpa)
   end if
   call bas%to_device()
   call aux%to_device()
   call pl%build(bas, 1.0e-12_dp)
   call pl%to_device()

   nocc = nint(sum(at_z))/2
   opts%conv_energy = 1.0e-12_dp
   opts%conv_density = 1.0e-9_dp
   if (comm%size() > 1) then
      call trc_scf_run(bas, nocc, nocc, opts, res, comm=comm)
      call trc_rimp2_run(bas, aux, pl, nocc, res%cmo(:, :, 1), res%eps(:, 1), mp, comm=comm)
      call trc_rimp2_run(bas, aux, pl, nocc, res%cmo(:, :, 1), res%eps(:, 1), mp2, aux_block=20, comm=comm)
   else
      call trc_scf_run(bas, nocc, nocc, opts, res)
      call trc_rimp2_run(bas, aux, pl, nocc, res%cmo(:, :, 1), res%eps(:, 1), mp)
      call trc_rimp2_run(bas, aux, pl, nocc, res%cmo(:, :, 1), res%eps(:, 1), mp2, aux_block=20)
   end if
   ok = res%converged .and. len_trim(mp%message) == 0
   if (len_trim(mp%message) > 0 .and. comm%rank() == 0) print '(a)', "scf_rimp2: "//trim(mp%message)
   if (abs(mp%e_corr - mp2%e_corr) > 1.0e-13_dp .or. abs(mp%e_os - mp2%e_os) > 1.0e-13_dp) ok = .false.

   ! The same through the C entries mqc uses: trc_scf_mpi then trc_rimp2 on
   ! the orbitals it kept, collectively when there are ranks.
   block
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
      use trc_c_interfaces, only: trc_basis_create, trc_basis_destroy, trc_scf_mpi, trc_rimp2, TRC_OK
#ifdef TERCO_HAVE_MPI
      use mpi_f08, only: MPI_COMM_WORLD
#endif
      integer(c_int) :: fcomm, rc, niter
      integer(c_int), allocatable :: lc(:), npc(:), lca(:), npca(:)
      real(c_double), allocatable :: zat(:), d_lib(:, :), eps_lib(:)
      real(c_double) :: e_lib, e_xc, c_os, c_ss
      type(c_ptr) :: hbas, haux
      fcomm = -1
#ifdef TERCO_HAVE_MPI
      if (comm%size() > 1) fcomm = MPI_COMM_WORLD%mpi_val
#endif
      allocate (zat(natm), lc(nsh), npc(nsh), lca(nsha), npca(nsha))
      zat = real(at_z, c_double); lc = int(sh_l, c_int); npc = int(sh_np, c_int)
      lca = int(ax_l, c_int); npca = int(ax_np, c_int)
      rc = trc_basis_create(int(nsh, c_int), int(maxnp, c_int), lc, npc, sh_e, sh_c, sh_r, &
                            int(natm, c_int), zat, at_r, hbas)
      if (rc /= TRC_OK) ok = .false.
      rc = trc_basis_create(int(nsha, c_int), int(maxnpa, c_int), lca, npca, ax_e, ax_c, ax_r, &
                            int(natm, c_int), zat, at_r, haux)
      if (rc /= TRC_OK) ok = .false.
      allocate (d_lib(bas%nao, bas%nao), eps_lib(bas%nao))
      rc = trc_scf_mpi(fcomm, hbas, int(nocc, c_int), int(nocc, c_int), c_null_char, 3_c_int, &
                       1.0e-12_c_double, 1.0e-9_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                       e_lib, e_xc, d_lib, eps_lib, niter)
      if (rc /= TRC_OK) ok = .false.
      rc = trc_rimp2(fcomm, hbas, haux, 0_c_int, 20_c_int, c_os, c_ss)
      if (rc /= TRC_OK) ok = .false.
      if (comm%rank() == 0) print '(a,es9.2,a,es9.2)', "  C entries vs driver:  E_RHF ", abs(e_lib - res%energy), &
         "   E_corr ", abs(c_os + c_ss - mp%e_corr)
      if (abs(e_lib - res%energy) > 1.0e-10_dp .or. abs(c_os + c_ss - mp%e_corr) > 1.0e-10_dp) ok = .false.
      rc = trc_basis_destroy(haux); rc = trc_basis_destroy(hbas)
   end block

   if (comm%rank() == 0) then
      print '(a,i0,a,i0,a,i0)', "scf_rimp2: nao ", bas%nao, "  naux ", aux%nao, "  ranks ", comm%size()
      print '(a,f20.12,a,i0)', "  E_RHF   = ", res%energy, "   iter ", res%iterations
      print '(a,f20.12,a,f20.12)', "  E_os    = ", mp%e_os, "   E_ss = ", mp%e_ss
      print '(a,f20.12,a,es9.2)', "  E_corr  = ", mp%e_corr, "   blocked vs whole ", abs(mp%e_corr - mp2%e_corr)
      open (newunit=unit, file=trim(probe), access='stream', form='unformatted', status='replace')
      write (unit) int(bas%nao, kind=8), int(natm, kind=8), int(nsh, kind=8), int(nsha, kind=8), int(nocc, kind=8)
      write (unit) at_z, at_r
      ! Shells as the JSON gives them: raw coefficients when read from a
      ! file, so pyscf normalises them its own way; unit-weight primitives
      ! for the synthetic sets, where the coefficient is the normalisation.
      do ish = 1, nsh
         write (unit) int(centre_of(sh_r(:, ish)), kind=8), int(sh_l(ish), kind=8), int(sh_np(ish), kind=8)
         write (unit) sh_e(1:sh_np(ish), ish), raw_coef(sh_l(ish), sh_np(ish), sh_e(1:sh_np(ish), ish), &
            sh_c(1:sh_np(ish), ish), len_trim(bas_json) > 0)
      end do
      do ish = 1, nsha
         write (unit) int(centre_of(ax_r(:, ish)), kind=8), int(ax_l(ish), kind=8), int(ax_np(ish), kind=8)
         write (unit) ax_e(1:ax_np(ish), ish), raw_coef(ax_l(ish), ax_np(ish), ax_e(1:ax_np(ish), ish), &
            ax_c(1:ax_np(ish), ish), len_trim(aux_json) > 0)
      end do
      write (unit) res%energy, mp%e_os, mp%e_ss
      close (unit)
      print '(a)', merge("scf_rimp2: PASS", "scf_rimp2: FAIL", ok)
   end if
   call pl%release(); call aux%release(); call bas%release()
   call pic_mpi_finalize()
   if (.not. ok) stop 1

contains

   ! The flat arrays back out of a built basis, so the C-entry check and the
   ! probe can be written from one set of names. sh_c has common_fac_sp
   ! folded in; take it out so trc_basis_create folds it once.
   subroutine unpack_basis(b, n, l, np, e, c, r, mnp)
      type(trc_basis_t), intent(in) :: b
      integer, intent(out) :: n, mnp
      integer, allocatable, intent(out) :: l(:), np(:)
      real(dp), allocatable, intent(out) :: e(:, :), c(:, :), r(:, :)
      integer :: i
      n = b%nshell; mnp = b%maxnp
      allocate (l(n), np(n), e(mnp, n), c(mnp, n), r(3, n))
      l = b%sh_l; np = b%sh_np; e = b%sh_e; c = b%sh_c; r = b%sh_r
      do i = 1, n
         if (l(i) == 0) c(:, i) = c(:, i)/0.282094791773878143_dp
         if (l(i) == 1) c(:, i) = c(:, i)/0.488602511902919921_dp
      end do
   end subroutine unpack_basis

   ! What pyscf is given as the coefficient. pyscf normalises a contraction
   ! itself from raw coefficients, so a file-read shell hands over the file
   ! shape: the normalised coefficients with the primitive norm divided out
   ! (the overall scale is pyscf's to fix). A synthetic shell is a single
   ! primitive whose coefficient IS the norm: pyscf gets 1.
   function raw_coef(l, np, e, c, from_file) result(raw)
      integer, intent(in) :: l, np
      real(dp), intent(in) :: e(np), c(np)
      logical, intent(in) :: from_file
      real(dp) :: raw(np)
      real(dp), parameter :: PI = 3.14159265358979323846_dp
      real(dp) :: dfac
      integer :: k
      if (.not. from_file) then
         raw = 1.0_dp
         return
      end if
      dfac = 1.0_dp
      do k = 1, l
         dfac = dfac*real(2*k - 1, dp)
      end do
      do k = 1, np
         raw(k) = c(k)/sqrt(2.0_dp**l/(PI**1.5_dp*dfac)*(2.0_dp*e(k))**(real(l, dp) + 1.5_dp))
      end do
   end function raw_coef

   integer function centre_of(r)
      real(dp), intent(in) :: r(3)
      integer :: a
      centre_of = 0
      do a = 1, natm
         if (all(abs(r - at_r(:, a)) < 1.0e-12_dp)) centre_of = a
      end do
   end function centre_of

end program scf_rimp2
