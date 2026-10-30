!
! RI-MP2 on water: RHF through the driver, then the correlation energy
! through trc_rimp2 with a synthetic auxiliary basis, once as a single
! auxiliary block and once in blocks of 20 functions (bit-identical
! required). One process: RI-MP2 over ranks is checked through the
! command line, by app/cli_mpi_check.py. Writes
! `rimp2_probe.bin` for rimp2_ref.py, which runs pyscf's RHF and DF-MP2 on
! the same orbital and auxiliary sets and compares E_os, E_ss and E_corr.
!
! Arguments: xyz [basis.json aux.json [probe [sad.json]]]. Without files,
! primitive 6-31G and the synthetic auxiliary ladder; with them, both sets
! from MolSSI BSE JSON through trc_basis_json, which is how mqc's sets
! arrive, and optionally mqc's SAD guess to start from. ONE (basis, aux)
! case per process (the SAME executable is registered twice, with
! different ARGS -- terco/scf_rimp2 and terco/scf_rimp2_ccpvdz), so this
! suite has a single test-drive case.
!
module scf_rimp2_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_rimp2_driver, only: trc_rimp2_run, trc_rimp2_result_t, trc_rimp2_options_t
   use trc_test_basis, only: read_xyz, build_aux, uncontract_basis
   use trc_basis_json, only: trc_basis_from_json
   use trc_guess_json, only: trc_guess_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none
   private
   public :: collect_scf_rimp2, xyz_arg, bas_json_arg, aux_json_arg, probe_arg, guess_json_arg

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = ''
   character(len=256) :: aux_json_arg = ''
   character(len=256) :: probe_arg = 'rimp2_probe.bin'
   character(len=256) :: guess_json_arg = ''

contains

   subroutine collect_scf_rimp2(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RI-MP2 vs pyscf, aux blocking, and the C entries", test_rimp2) &
                  ]
   end subroutine collect_scf_rimp2

   subroutine test_rimp2(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: natm, nsh, maxnp, nsha, maxnpa, unit, ish, dev, nocc
      integer, allocatable :: zint(:), sh_l(:), sh_np(:), ax_l(:), ax_np(:)
      real(dp), allocatable :: at_r(:, :), at_z(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
      real(dp), allocatable :: ax_e(:, :), ax_c(:, :), ax_r(:, :)
      type(trc_basis_t) :: bas, aux
      type(trc_pairlist_t) :: pl
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(trc_rimp2_result_t) :: mp, mp2
      type(error_t) :: err
      real(dp), allocatable :: dguess(:, :), dg3(:, :, :)
      logical :: c_ok

      dev = trc_bind_device(0)

      call read_xyz(trim(xyz_arg), natm, zint, at_r)
      allocate (at_z(natm)); at_z = real(zint, dp)
      if (len_trim(bas_json_arg) > 0) then
         ! Forced Cartesian: RI-MP2 (below) is Cartesian-only (trc_rimp2_driver.F90's
         ! own refusal of a spherical basis), and mqc's SAD guess fixture, when
         ! one is given, is in terco's Cartesian AO convention regardless.
         call trc_basis_from_json(trim(bas_json_arg), natm, zint, at_r, bas, err, force_cart=.true.)
         if (.not. err%has_error()) call trc_basis_from_json(trim(aux_json_arg), natm, zint, at_r, aux, err, &
                                                              force_cart=.true.)
         call check(error, .not. err%has_error(), err%get_message())
         if (allocated(error)) return
         call unpack_basis(bas, natm, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
         if (len_trim(guess_json_arg) > 0) then
            call trc_guess_from_json(trim(guess_json_arg), bas, zint, dguess, err)
            call check(error, .not. err%has_error(), err%get_message())
            if (allocated(error)) return
            allocate (dg3(bas%nao, bas%nao, 1)); dg3(:, :, 1) = dguess
         end if
         call unpack_basis(aux, natm, nsha, ax_l, ax_np, ax_e, ax_c, ax_r, maxnpa)
      else
         ! 6-31G, through the one basis-construction path, then uncontracted
         ! to primitives exactly as the old `build_631g(..., uncontracted=
         ! .true.)` call did (see uncontract_basis's own comment). `bas` is
         ! reused as scratch for the JSON load -- `%build` below releases
         ! and rebuilds it, so nothing of the intermediate load survives.
         call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, bas, err)
         call check(error, .not. err%has_error(), err%get_message())
         if (allocated(error)) return
         call uncontract_basis(bas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
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
      ! TERCO_VERBOSE in the environment: the SCF iterations and the RI-MP2
      ! steps with their times.
      block
         character(len=8) :: env
         integer :: st
         call get_environment_variable("TERCO_VERBOSE", env, status=st)
         opts%verbose = st == 0
      end block
      if (allocated(dg3)) then
         call trc_scf_run(bas, nocc, nocc, opts, res, dguess=dg3)
      else
         call trc_scf_run(bas, nocc, nocc, opts, res)
      end if
      block
         type(trc_rimp2_options_t) :: rimp2_opts
         rimp2_opts%verbose = opts%verbose
         call trc_rimp2_run(bas, aux, pl, nocc, res%cmo(:, :, 1), res%eps(:, 1), mp, opts=rimp2_opts)
         rimp2_opts%verbose = .false.
         rimp2_opts%aux_block = 20
         call trc_rimp2_run(bas, aux, pl, nocc, res%cmo(:, :, 1), res%eps(:, 1), mp2, opts=rimp2_opts)
      end block
      if (len_trim(mp%message) > 0) print '(a)', "scf_rimp2: "//trim(mp%message)

      print '(a,i0,a,i0)', "scf_rimp2: nao ", bas%nao, "  naux ", aux%nao
      print '(a,f20.12,a,i0)', "  E_RHF   = ", res%energy, "   iter ", res%iterations
      print '(a,f20.12,a,f20.12)', "  E_os    = ", mp%e_os, "   E_ss = ", mp%e_ss
      print '(a,f20.12,a,es9.2)', "  E_corr  = ", mp%e_corr, "   blocked vs whole ", abs(mp%e_corr - mp2%e_corr)

      ! The same through the C entries mqc uses: trc_scf_mpi then trc_rimp2 on
      ! the orbitals it kept, with no communicator (fcomm = -1).
      block
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
         use trc_c_interfaces, only: trc_basis_create, trc_basis_destroy, trc_scf_mpi, trc_rimp2, TRC_OK
         integer(c_int) :: fcomm, rc, niter
         integer(c_int), allocatable :: lc(:), npc(:), lca(:), npca(:)
         real(c_double), allocatable :: zat(:), d_lib(:, :), eps_lib(:)
         real(c_double) :: e_lib, e_xc, c_os, c_ss
         type(c_ptr) :: hbas, haux
         c_ok = .true.
         fcomm = -1
         allocate (zat(natm), lc(nsh), npc(nsh), lca(nsha), npca(nsha))
         zat = real(at_z, c_double); lc = int(sh_l, c_int); npc = int(sh_np, c_int)
         lca = int(ax_l, c_int); npca = int(ax_np, c_int)
         rc = trc_basis_create(int(nsh, c_int), int(maxnp, c_int), lc, npc, sh_e, sh_c, sh_r, &
                               int(natm, c_int), zat, at_r, hbas)
         if (rc /= TRC_OK) c_ok = .false.
         rc = trc_basis_create(int(nsha, c_int), int(maxnpa, c_int), lca, npca, ax_e, ax_c, ax_r, &
                               int(natm, c_int), zat, at_r, haux)
         if (rc /= TRC_OK) c_ok = .false.
         allocate (d_lib(bas%nao, bas%nao), eps_lib(bas%nao))
         rc = trc_scf_mpi(fcomm, hbas, int(nocc, c_int), int(nocc, c_int), c_null_char, 3_c_int, &
                          real(opts%conv_energy, c_double), real(opts%conv_density, c_double), 100_c_int, &
                          c_null_ptr, 0_c_int, e_lib, e_xc, d_lib, eps_lib, niter)
         if (rc /= TRC_OK) c_ok = .false.
         rc = trc_rimp2(fcomm, hbas, haux, 0_c_int, 20_c_int, c_os, c_ss)
         if (rc /= TRC_OK) c_ok = .false.
         print '(a,es9.2,a,es9.2)', "  C entries vs driver:  E_RHF ", abs(e_lib - res%energy), &
            "   E_corr ", abs(c_os + c_ss - mp%e_corr)
         ! Two SCFs converged separately, then two RI-MP2s over their orbitals:
         ! the atomic accumulation order alone puts ~1e-10 between them.
         if (abs(e_lib - res%energy) > 1.0e-10_dp .or. abs(c_os + c_ss - mp%e_corr) > 1.0e-9_dp) c_ok = .false.
         rc = trc_basis_destroy(haux); rc = trc_basis_destroy(hbas)
      end block

      ! Golden case named after the probe, as the pyscf reference is below:
      ! rimp2_probe.bin -> scf_rimp2, rimp2_ccpvdz_probe.bin -> scf_rimp2_ccpvdz.
      if (index(trim(probe_arg), '_probe.bin') > 0) then
         call golden_begin('scf_'//probe_arg(1:index(trim(probe_arg), '_probe.bin') - 1))
         call golden_put('E_RHF', res%energy)
         call golden_put('E_OS', mp%e_os)
         call golden_put('E_SS', mp%e_ss)
         call golden_put('E_CORR', mp%e_corr)
         call golden_end()
      end if
      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(bas%nao, int64), int(natm, int64), int(nsh, int64), int(nsha, int64), int(nocc, int64)
      write (unit) at_z, at_r
      ! Shells as the JSON gives them: raw coefficients when read from a
      ! file, so pyscf normalises them its own way; unit-weight primitives
      ! for the synthetic sets, where the coefficient is the normalisation.
      do ish = 1, nsh
         write (unit) int(centre_of(sh_r(:, ish), at_r, natm), int64), int(sh_l(ish), int64), int(sh_np(ish), int64)
         write (unit) sh_e(1:sh_np(ish), ish), raw_coef(sh_l(ish), sh_np(ish), sh_e(1:sh_np(ish), ish), &
            sh_c(1:sh_np(ish), ish), len_trim(bas_json_arg) > 0)
      end do
      do ish = 1, nsha
         write (unit) int(centre_of(ax_r(:, ish), at_r, natm), int64), int(ax_l(ish), int64), int(ax_np(ish), int64)
         write (unit) ax_e(1:ax_np(ish), ish), raw_coef(ax_l(ish), ax_np(ish), ax_e(1:ax_np(ish), ish), &
            ax_c(1:ax_np(ish), ish), len_trim(aux_json_arg) > 0)
      end do
      write (unit) res%energy, mp%e_os, mp%e_ss
      close (unit)

      ! The golden and the probe are written above whatever happens below:
      ! test/pyscf_refs/generate.py's gen_rimp2 reads the probe even when a
      ! check here fails (e.g. a regeneration where the C entries disagree).
      call pl%release(); call aux%release(); call bas%release()
      call check(error, res%converged, "SCF did not converge")
      if (allocated(error)) return
      call check(error, len_trim(mp%message) == 0, "RI-MP2: "//trim(mp%message))
      if (allocated(error)) return
      call check(error, abs(mp%e_corr - mp2%e_corr) <= 1.0e-13_dp .and. abs(mp%e_os - mp2%e_os) <= 1.0e-13_dp, &
                 "aux blocking changed the correlation energy")
      if (allocated(error)) return
      call check(error, c_ok, "the C entries (trc_scf_mpi/trc_rimp2) disagree with the Fortran driver")
      if (allocated(error)) return

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      ! test/pyscf_refs/generate.py's gen_rimp2 names each case's reference
      ! after its own probe file, "<probe base>.ref" in place of
      ! "<probe base>_probe.bin" -- rimp2_probe.bin -> rimp2.ref,
      ! rimp2_ccpvdz_probe.bin -> rimp2_ccpvdz.ref.
      block
         integer :: u, iprobe
         real(dp) :: e_scf_ref, e_os_ref, e_ss_ref
         real(dp) :: pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
         character(len=256) :: refname
         iprobe = index(trim(probe_arg), '_probe.bin')
         if (iprobe == 0) error stop 'scf_rimp2: cannot derive a reference name from '//trim(probe_arg)
         refname = probe_arg(1:iprobe - 1)//'.ref'
         call ref_require(trim(refname), u)
         call ref_read_real(u, 'E_SCF', e_scf_ref)
         call ref_read_real(u, 'E_OS', e_os_ref)
         call ref_read_real(u, 'E_SS', e_ss_ref)
         ! Provenance only -- see scf_rks.F90's own comment on the same tags.
         call ref_read_real(u, 'CONV_TOL', pyscf_conv_tol)
         call ref_read_real(u, 'CONV_TOL_GRAD', pyscf_conv_tol_grad)
         call ref_read_real(u, 'GRAD_NORM', pyscf_grad_norm)
         call ref_close(u)

         call check(error, abs(res%energy - e_scf_ref) < 1.0e-9_dp, &
                    "E_RHF vs pyscf, |dE| = "//to_string(abs(res%energy - e_scf_ref), "es10.2"))
         if (allocated(error)) return
         call check(error, abs(mp%e_os - e_os_ref) < 1.0e-9_dp, &
                    "E_os vs pyscf, |dE| = "//to_string(abs(mp%e_os - e_os_ref), "es10.2"))
         if (allocated(error)) return
         call check(error, abs(mp%e_ss - e_ss_ref) < 1.0e-9_dp, &
                    "E_ss vs pyscf, |dE| = "//to_string(abs(mp%e_ss - e_ss_ref), "es10.2"))
      end block
   end subroutine test_rimp2

   ! The flat arrays back out of a built basis, so the C-entry check and the
   ! probe can be written from one set of names. sh_c has common_fac_sp
   ! folded in; take it out so trc_basis_create folds it once.
   subroutine unpack_basis(b, natm, n, l, np, e, c, r, mnp)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: natm
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

   integer function centre_of(r, at_r, natm)
      real(dp), intent(in) :: r(3), at_r(:, :)
      integer, intent(in) :: natm
      integer :: a
      centre_of = 0
      do a = 1, natm
         if (all(abs(r - at_r(:, a)) < 1.0e-12_dp)) centre_of = a
      end do
   end function centre_of

end module scf_rimp2_suite

program scf_rimp2
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_rimp2_suite, only: collect_scf_rimp2, xyz_arg, bas_json_arg, aux_json_arg, probe_arg, guess_json_arg
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 3) then
      call get_command_argument(2, bas_json_arg)
      call get_command_argument(3, aux_json_arg)
   end if
   if (command_argument_count() >= 4) call get_command_argument(4, probe_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, guess_json_arg)

   stat = 0
   call run_testsuite(collect_scf_rimp2, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_rimp2
