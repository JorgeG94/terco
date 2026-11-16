module link_gfortran_suite
   !!
   !! terco driven through its C ABI from Fortran, via the shipped interface.
   !!
   !! WHY THIS TEST EXISTS
   !! --------------------
   !! A host code built with another compiler (mqc is built with gfortran, terco
   !! with nvfortran) cannot read terco's `.mod` files, so the only usable
   !! interface is `trc_capi.F90`'s `bind(c)` skin and its Fortran declaration,
   !! `trc_c_interfaces`. This test compiles that declaration with the project
   !! compiler and gets the right number out of the library through it. It uses
   !! terco's modules only for set-up (`trc_test_basis`'s xyz reader and basis
   !! unpacking, `trc_basis_json`'s basis construction); every call it checks
   !! goes through `trc_c_interfaces`.
   !!
   !! It runs a complete RHF SCF: basis, pair list, Schwarz bounds, one-electron
   !! matrices and every Fock build come from terco across the C boundary, and
   !! only the eigenproblem is local. Water at 6-31G must give
   !!
   !!     -75.983974465717   (trc_basis_from_json's own value, 2026-10-08)
   !!
   !! which is what `test/standalone/scf_rhf.F90` gets natively and what pyscf gets (both
   !! to 1e-8, through the SAME basis-construction path, `trc_basis_from_json`,
   !! this file now also uses -- see trc_test_basis.F90's header for why
   !! `build_631g` is no longer that path; it used to be build_631g's own
   !! 7-figure-literal value, -75.983974469890, 4.2e-9 Eh away from this one).
   !! Anything
   !! else means the ABI is misdescribing an array somewhere -- and a transposed
   !! or mis-strided matrix usually still converges, to the wrong energy, so the
   !! energy is the check and "it ran" is not.
   !!
   !! It also exercises EVERY entry point in the ABI, not just the ones an SCF
   !! needs, so that a signature change anywhere is caught here.
   !!
   !! Built as the `link_gfortran` target (test/CMakeLists.txt) and run as
   !!   link_gfortran <xyz>
   !! A test-drive suite of two cases: the SCF and its energy, and every
   !! remaining entry point.
   !!
   use, intrinsic :: iso_c_binding
   !
   ! The interfaces come from terco, not from a copy written here. A copy is
   ! a contract kept in step by hand, and nothing catches it drifting: add an
   ! argument to `trc_capi.F90`, forget the copy, and it still compiles and
   ! still links, because a C symbol carries no signature. It fails at run
   ! time by reading whatever was in the register.
   !
   ! Compiling the shipped declaration here and then asserting
   ! an energy is what makes that arrangement self-checking -- the drift shows
   ! up as a wrong number in terco's own test suite rather than in a host
   ! code months later.
   !
   use pic_strings, only: to_string
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   use trc_c_interfaces
   use trc_test_basis, only: read_xyz, unpack_basis
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_link_gfortran, release_link_gfortran

   ! The SCF is converged once, by whichever case runs first, and shared:
   ! the cases read the converged density and the handles from here.
   logical :: ready = .false.
   logical :: failed = .false.
   character(len=256) :: xyzfile
   integer :: nat, nsh, maxnp, nocc, nao
   integer,  allocatable :: at_z(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), sh_e(:, :), sh_c(:, :), sh_r(:, :)
   real(c_double), allocatable :: zat(:)
   type(c_ptr) :: hbas = c_null_ptr, hpair = c_null_ptr, heri = c_null_ptr
   real(dp), allocatable :: s(:, :), d(:, :)
   real(dp) :: escf
   integer :: iters
   real(dp), parameter :: thresh = 1.0e-10_dp

contains

   subroutine collect_link_gfortran(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF through the C ABI matches the native energy", test_scf_energy), &
                  new_unittest("every other entry point", test_entry_points) &
                  ]
   end subroutine collect_link_gfortran

   subroutine release_link_gfortran()
      !! Destroy the handles the first case created.
      call release_handles()
      ready = .false.
   end subroutine release_link_gfortran

   subroutine release_handles()
      !! Destroy whichever handles were created (a null handle was not).
      integer(c_int) :: rc
      if (c_associated(heri)) rc = trc_eri_destroy(heri)
      if (c_associated(hpair)) rc = trc_pairs_destroy(hpair)
      if (c_associated(hbas)) rc = trc_basis_destroy(hbas)
      heri = c_null_ptr; hpair = c_null_ptr; hbas = c_null_ptr
   end subroutine release_handles

   subroutine check_rc(error, code, what)
      !! A nonzero ABI status fails the case.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: code
      character(len=*), intent(in) :: what
      if (code /= 0) call test_failed(error, what//" failed with status "//to_string(code))
   end subroutine check_rc

   subroutine setup(error)
      !! Converge the shared SCF once. A failure releases what was created and
      !! is remembered, so a later case fails with a message instead of
      !! re-running into already-allocated state.
      type(error_type), allocatable, intent(out) :: error
      if (ready) return
      if (failed) then
         call test_failed(error, "an earlier case failed during set-up"); return
      end if
      call setup_build(error)
      if (allocated(error)) then
         failed = .true.
         call release_handles()
      end if
   end subroutine setup

   subroutine setup_build(error)
      !! Water 6-31G through the C ABI, then a plain RHF loop with a local
      !! eigensolver: only the Fock builds and one-electron matrices are terco's.
      type(error_type), allocatable, intent(out) :: error
      integer :: i, j, it, rc
      integer(c_int) :: nao_c
      integer(c_int), allocatable :: lc(:), npc(:)
      type(trc_basis_t) :: jbas
      type(error_t) :: err
      real(dp), allocatable :: t(:, :), v(:, :), h(:, :)
      real(dp), allocatable :: g(:, :), f(:, :), x(:, :), c(:, :), e(:)
      real(dp) :: enuc, eold

      if (command_argument_count() < 1) then
         call test_failed(error, "usage: link_gfortran <xyz>"); return
      end if
      call get_command_argument(1, xyzfile)

      call read_xyz(xyzfile, nat, at_z, at_r)
      call trc_basis_from_json('basis_sets/6-31g.json', nat, at_z, at_r, jbas, err)
      if (err%has_error()) then
         call test_failed(error, "link_gfortran: "//err%get_message()); return
      end if
      call unpack_basis(jbas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
      allocate (zat(nat), lc(nsh), npc(nsh))
      zat = real(at_z, c_double)
      lc = int(sh_l, c_int); npc = int(sh_np, c_int)

      rc = trc_basis_create(int(nsh, c_int), int(maxnp, c_int), lc, npc, &
                            sh_e, sh_c, sh_r, int(nat, c_int), zat, at_r, hbas)
      call check_rc(error, rc, 'trc_basis_create'); if (allocated(error)) return
      rc = trc_basis_nao(hbas, nao_c)
      call check_rc(error, rc, 'trc_basis_nao'); if (allocated(error)) return
      nao = int(nao_c)
      rc = trc_pairs_create(hbas, thresh, hpair)
      call check_rc(error, rc, 'trc_pairs_create'); if (allocated(error)) return
      rc = trc_eri_create(hbas, thresh, heri)
      call check_rc(error, rc, 'trc_eri_create'); if (allocated(error)) return

      allocate (s(nao, nao), t(nao, nao), v(nao, nao), h(nao, nao))
      rc = trc_compute_1e(hbas, hpair, s, t, v)
      call check_rc(error, rc, 'trc_compute_1e'); if (allocated(error)) return
      h = t + v

      enuc = 0.0_dp
      do i = 2, nat
         do j = 1, i - 1
            enuc = enuc + zat(i)*zat(j)/norm2(at_r(:, i) - at_r(:, j))
         end do
      end do

      nocc = sum(at_z)/2
      allocate (d(nao, nao), g(nao, nao), f(nao, nao), x(nao, nao))
      allocate (c(nao, nao), e(nao))
      call sym_orthog(nao, s, x)
      call diag(nao, h, x, c, e)
      d = 2.0_dp*matmul(c(:, 1:nocc), transpose(c(:, 1:nocc)))
      eold = 0.0_dp
      do it = 1, 200
         rc = trc_fock(heri, hbas, d, g, 1.0_c_double, 1.0_c_double, 1_c_int, 0.0_c_double)
         call check_rc(error, rc, 'trc_fock'); if (allocated(error)) return
         f = h + g
         escf = enuc + 0.5_dp*sum(d*(h + f))
         call diag(nao, f, x, c, e)
         d = 2.0_dp*matmul(c(:, 1:nocc), transpose(c(:, 1:nocc)))
         if (abs(escf - eold) < 1.0e-11_dp) exit
         eold = escf
      end do
      iters = it
      ready = .true.
   end subroutine setup_build

   subroutine test_scf_energy(error)
      !! The converged energy, the multipoles, and the library's own driver
      !! across the same ABI.
      type(error_type), allocatable, intent(out) :: error
      integer :: rc, nmult
      real(dp), allocatable :: mm(:, :, :)

      call setup(error); if (allocated(error)) return

      nmult = int(trc_multipole_count())
      allocate (mm(nao, nao, nmult))
      rc = trc_compute_multipoles(hbas, hpair, [0.0_c_double, 0.0_c_double, &
                                  0.0_c_double], mm)
      call check_rc(error, rc, 'trc_compute_multipoles'); if (allocated(error)) return

      print '(a)', '  terco through the C ABI, via trc_c_interfaces'
      print '(a,i0,a,i0)', '    nao ', nao, '   iterations ', iters
      print '(a,f18.10)',  '    RHF energy   ', escf
      print '(a,f18.10)',  '    sum(S)       ', sum(s)
      print '(a,3f14.8)',  '    dipole <r>   ', -sum(d*mm(:, :, 1)), &
                                                -sum(d*mm(:, :, 2)), &
                                                -sum(d*mm(:, :, 3))
      ! The same SCF through the library's own driver, across the same ABI.
      ! It must land on the energy the loop above found, and PBE must converge.
      block
         real(c_double) :: e_lib, e_xc
         real(c_double), allocatable :: d_lib(:, :), eps_lib(:)
         integer(c_int) :: niter
         allocate (d_lib(nao, nao), eps_lib(nao))
         rc = trc_scf(hbas, int(nocc, c_int), int(nocc, c_int), c_null_char, 3_c_int, &
                      1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                      e_lib, e_xc, d_lib, eps_lib, niter)
         call check_rc(error, rc, 'trc_scf (HF)'); if (allocated(error)) return
         print '(a,f18.10,a,i0,a)', '    trc_scf RHF  ', e_lib, '   (', niter, ' iterations)'
         if (.not. (abs(e_lib - escf) < 1.0e-8_dp)) then
            call test_failed(error, "trc_scf disagrees with the loop above"); return
         end if
         rc = trc_scf(hbas, int(nocc, c_int), int(nocc, c_int), 'pbe'//c_null_char, 3_c_int, &
                      1.0e-10_c_double, 1.0e-7_c_double, 100_c_int, c_null_ptr, 0_c_int, &
                      e_lib, e_xc, d_lib, eps_lib, niter)
         call check_rc(error, rc, 'trc_scf (PBE)'); if (allocated(error)) return
         print '(a,f18.10,a,f16.10,a,i0,a)', '    trc_scf PBE  ', e_lib, '   E_xc ', e_xc, '   (', niter, ' iterations)'
      end block
      ! trc_basis_from_json's own 6-31G water RHF value, 2026-10-08.
      if (.not. (abs(escf + 75.983974465717_dp) < 1.0e-8_dp)) then
         call test_failed(error, "RHF energy does not match the native build and pyscf: "//to_string(escf))
      end if
   end subroutine test_scf_energy

   subroutine test_entry_points(error)
      type(error_type), allocatable, intent(out) :: error
      call setup(error); if (allocated(error)) return
      call cover_rest(error)
   end subroutine test_entry_points

   subroutine cover_rest(error)
      !!
      !! Every remaining entry point, each with an invariant that a wrong
      !! signature would break.
      !!
      !! Coverage is the point. An SCF alone touches eight of the eighteen calls,
      !! so a changed argument list on any of the other ten would sail through
      !! this test and surface in a host code instead.
      !!
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: maxl, npair, nao_c
      real(dp), allocatable :: g1(:, :), gn(:, :), gm(:, :, :)
      real(dp), allocatable :: j2c(:, :), t3c(:, :, :)
      real(dp) :: worst
      integer :: i, j, rc
      type(c_ptr) :: hbas2, hpair2

      ! --- basis_maxl: 6-31G for H and O is s and p, nothing higher --------
      rc = trc_basis_maxl(hbas, maxl)
      call check_rc(error, rc, 'trc_basis_maxl'); if (allocated(error)) return
      if (maxl /= 1) then
         call test_failed(error, "basis_maxl should be 1 for 6-31G, got "//to_string(int(maxl)))
         return
      end if

      ! --- pairs_count -----------------------------------------------------
      rc = trc_pairs_count(hpair, npair)
      call check_rc(error, rc, 'trc_pairs_count'); if (allocated(error)) return
      if (npair <= 0) then
         call test_failed(error, "pairs_count returned nothing"); return
      end if

      ! --- the libcint-layout constructor ----------------------------------
      !
      ! This is the path mqc uses, so it gets the strongest check available:
      ! build the SAME basis both ways and require the same overlap. A
      ! mistake in the packed layout -- a slot off by one, a pointer that is
      ! 1-based where libcint is 0-based -- changes the basis rather than
      ! failing, so comparing against the array constructor is what catches
      ! it.
      call build_libcint_arrays(error, hbas2, hpair2); if (allocated(error)) return
      rc = trc_basis_nao(hbas2, nao_c)
      call check_rc(error, rc, 'nao of libcint basis'); if (allocated(error)) return
      if (int(nao_c) /= nao) then
         call test_failed(error, "libcint-layout basis has nao "//to_string(int(nao_c))// &
                          " against "//to_string(nao))
         return
      end if
      block
         real(dp), allocatable :: s2(:, :), t2(:, :), v2(:, :)
         allocate (s2(nao, nao), t2(nao, nao), v2(nao, nao))
         rc = trc_compute_1e(hbas2, hpair2, s2, t2, v2)
         call check_rc(error, rc, '1e over the libcint-layout basis'); if (allocated(error)) return
         worst = maxval(abs(s2 - s))
         print '(a,es10.2)', '    libcint-layout basis: overlap agrees to ', worst
         if (.not. (worst < 1.0e-12_dp)) then
            call test_failed(error, "the two constructors describe different bases")
            return
         end if
         deallocate (s2, t2, v2)
      end block

      ! --- fock_many against fock ------------------------------------------
      !
      ! Two copies of the converged density must give two copies of the
      ! answer `fock` already gave. That checks the batched entry point AND
      ! the (ndens, nao, nao) layout at once: with the density index in the
      ! wrong position this reads across the wrong stride and disagrees.
      allocate (g1(nao, nao), gn(nao, nao), gm(2, nao, nao))
      rc = trc_fock(heri, hbas, d, g1, 1.0_c_double, 1.0_c_double, 1_c_int, 0.0_c_double)
      call check_rc(error, rc, 'trc_fock'); if (allocated(error)) return
      block
         real(dp), allocatable :: dm(:, :, :)
         allocate (dm(2, nao, nao))
         do j = 1, nao
            do i = 1, nao
               dm(1, i, j) = d(i, j)
               dm(2, i, j) = d(i, j)
            end do
         end do
         rc = trc_fock_many(heri, hbas, 2_c_int, dm, gm, 1.0_c_double, &
                            1.0_c_double, 1_c_int, 0.0_c_double)
         call check_rc(error, rc, 'trc_fock_many'); if (allocated(error)) return
         worst = 0.0_dp
         do j = 1, nao
            do i = 1, nao
               worst = max(worst, abs(gm(1, i, j) - g1(i, j)))
               worst = max(worst, abs(gm(2, i, j) - g1(i, j)))
            end do
         end do
         deallocate (dm)
      end block
      print '(a,es10.2)', '    fock_many vs fock                    ', worst
      if (.not. (worst < 1.0e-10_dp)) then
         call test_failed(error, "the batched build disagrees with the single one")
         return
      end if

      ! --- fock_nosym against fock -----------------------------------------
      !
      ! On a SYMMETRIC density the folded and enumerated digestions must
      ! agree. That is the control that makes a disagreement on an
      ! asymmetric density mean the asymmetry rather than a bug.
      rc = trc_fock_nosym(heri, hbas, d, gn, 1_c_int)
      call check_rc(error, rc, 'trc_fock_nosym'); if (allocated(error)) return
      worst = maxval(abs(gn - g1))
      print '(a,es10.2)', '    fock_nosym vs fock (symmetric D)     ', worst
      print '(a,es24.16,a,es24.16)', '    sum folded ', sum(g1), &
         '   sum enumerated ', sum(gn)
      print '(a,f18.10)', '    ratio ', sum(gn)/sum(g1)
      if (.not. (worst < 1.0e-10_dp)) then
         call test_failed(error, "the enumerated digestion disagrees on a symmetric density")
         return
      end if

      ! --- density fitting, with the orbital basis as its own auxiliary ----
      allocate (j2c(nao, nao), t3c(nao, nao, nao))
      rc = trc_compute_df2c(hbas, j2c)
      call check_rc(error, rc, 'trc_compute_df2c'); if (allocated(error)) return
      rc = trc_compute_df3c(hbas, hpair, hbas, t3c)
      call check_rc(error, rc, 'trc_compute_df3c'); if (allocated(error)) return
      ! (P|Q) is symmetric and positive definite, so its diagonal is positive
      if (.not. ieee_ok(j2c) .or. .not. (minval([(j2c(i, i), i=1, nao)]) > 0.0_dp)) then
         call test_failed(error, "the two-centre metric is not positive on its diagonal")
         return
      end if
      print '(a,f16.8,a,f16.8)', '    df 2c trace ', &
         sum([(j2c(i, i), i=1, nao)]), '   3c sum ', sum(t3c)

      rc = trc_pairs_destroy(hpair2); rc = trc_basis_destroy(hbas2)
      deallocate (g1, gn, gm, j2c, t3c)
      allocate (g1(nao, nao))

      ! --- omega refusal: negative and NaN must be TRC_ERR_BADARG, not a
      ! crash or a silently wrong G ------------------------------------
      block
         use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
         real(c_double) :: nan_omega
         nan_omega = ieee_value(nan_omega, ieee_quiet_nan)
         rc = trc_fock(heri, hbas, d, g1, 1.0_c_double, 1.0_c_double, &
                       1_c_int, -1.0_c_double)
         if (rc /= TRC_ERR_BADARG) then
            call test_failed(error, "trc_fock(omega=-1) should refuse with "// &
                             to_string(int(TRC_ERR_BADARG))//", got "//to_string(rc))
            return
         end if
         rc = trc_fock(heri, hbas, d, g1, 1.0_c_double, 1.0_c_double, &
                       1_c_int, nan_omega)
         if (rc /= TRC_ERR_BADARG) then
            call test_failed(error, "trc_fock(omega=NaN) should refuse with "// &
                             to_string(int(TRC_ERR_BADARG))//", got "//to_string(rc))
            return
         end if
         block
            real(dp), allocatable :: dm(:, :, :), gm2(:, :, :)
            allocate (dm(2, nao, nao), gm2(2, nao, nao))
            dm(1, :, :) = d; dm(2, :, :) = d
            rc = trc_fock_many(heri, hbas, 2_c_int, dm, gm2, 1.0_c_double, &
                               1.0_c_double, 1_c_int, -1.0_c_double)
            if (rc /= TRC_ERR_BADARG) then
               call test_failed(error, "trc_fock_many(omega=-1) should refuse with "// &
                                to_string(int(TRC_ERR_BADARG))//", got "//to_string(rc))
               return
            end if
            rc = trc_fock_many(heri, hbas, 2_c_int, dm, gm2, 1.0_c_double, &
                               1.0_c_double, 1_c_int, nan_omega)
            if (rc /= TRC_ERR_BADARG) then
               call test_failed(error, "trc_fock_many(omega=NaN) should refuse with "// &
                                to_string(int(TRC_ERR_BADARG))//", got "//to_string(rc))
               return
            end if
            deallocate (dm, gm2)
         end block
         print '(a)', '    trc_fock/trc_fock_many refuse omega<0 and omega=NaN'
      end block
   end subroutine cover_rest

   logical function ieee_ok(a)
      real(dp), intent(in) :: a(:, :)
      ieee_ok = all(a == a) .and. all(abs(a) < huge(1.0_dp))
   end function ieee_ok

   subroutine build_libcint_arrays(error, hbas2, hpair2)
      !!
      !! Pack the same basis into libcint's atm/bas/env.
      !!
      !! The slot numbers are libcint's, 0-based, so a Fortran array indexed from
      !! one carries them at index+1. `env` pointers are 0-based offsets into
      !! `env` itself, which is the part that is easy to get wrong from Fortran
      !! and the reason the overlap is compared afterwards rather than trusted.
      !!
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(out) :: hbas2, hpair2
      integer, parameter :: PTR_ENV_START = 20
      integer(c_int), allocatable :: atm(:, :), bas(:, :)
      real(c_double), allocatable :: env(:)
      integer :: off, ia, is, k, ncoord, rc

      ncoord = 3*nat
      allocate (atm(6, nat), bas(8, nsh))
      allocate (env(PTR_ENV_START + ncoord + 2*sum(sh_np)))
      env = 0.0_c_double
      atm = 0_c_int; bas = 0_c_int

      off = PTR_ENV_START
      do ia = 1, nat
         atm(1, ia) = int(at_z(ia), c_int)      ! CHARGE_OF
         atm(2, ia) = int(off, c_int)           ! PTR_COORD
         env(off + 1:off + 3) = real(at_r(:, ia), c_double)
         off = off + 3
      end do

      do is = 1, nsh
         bas(1, is) = int(atom_of(is) - 1, c_int)   ! ATOM_OF, 0-based
         bas(2, is) = int(sh_l(is), c_int)          ! ANG_OF
         bas(3, is) = int(sh_np(is), c_int)         ! NPRIM_OF
         bas(4, is) = 1_c_int                       ! NCTR_OF
         bas(5, is) = 0_c_int                       ! KAPPA_OF
         bas(6, is) = int(off, c_int)               ! PTR_EXP
         do k = 1, sh_np(is)
            env(off + k) = real(sh_e(k, is), c_double)
         end do
         off = off + sh_np(is)
         bas(7, is) = int(off, c_int)               ! PTR_COEFF
         do k = 1, sh_np(is)
            env(off + k) = real(sh_c(k, is), c_double)
         end do
         off = off + sh_np(is)
      end do

      rc = trc_basis_create_libcint(atm, int(nat, c_int), bas, &
                                    int(nsh, c_int), env, &
                                    int(size(env), c_int), 1_c_int, hbas2)
      call check_rc(error, rc, 'trc_basis_create_libcint'); if (allocated(error)) return
      rc = trc_pairs_create(hbas2, thresh, hpair2)
      call check_rc(error, rc, 'pairs over the libcint-layout basis')
      deallocate (atm, bas, env)
   end subroutine build_libcint_arrays

   integer function atom_of(is)
      !! Which atom a shell sits on, by matching its centre.
      integer, intent(in) :: is
      integer :: ia
      atom_of = 1
      do ia = 1, nat
         if (maxval(abs(sh_r(:, is) - at_r(:, ia))) < 1.0e-12_dp) atom_of = ia
      end do
   end function atom_of

   subroutine sym_orthog(n, sm, xo)
      integer, intent(in) :: n
      real(dp), intent(in) :: sm(n, n)
      real(dp), intent(out) :: xo(n, n)
      real(dp), allocatable :: u(:, :), w(:), work(:)
      integer :: info, lw, ii, jj
      allocate (u(n, n), w(n), work(64*n)); u = sm; lw = 64*n
      call dsyev('V', 'U', n, u, n, w, work, lw, info)
      do jj = 1, n
         do ii = 1, n
            xo(ii, jj) = u(ii, jj)/sqrt(w(jj))
         end do
      end do
      xo = matmul(xo, transpose(u))
   end subroutine sym_orthog

   subroutine diag(n, fm, xo, cm, ev)
      integer, intent(in) :: n
      real(dp), intent(in) :: fm(n, n), xo(n, n)
      real(dp), intent(out) :: cm(n, n), ev(n)
      real(dp), allocatable :: fp(:, :), work(:)
      integer :: info, lw
      allocate (fp(n, n), work(64*n)); lw = 64*n
      fp = matmul(transpose(xo), matmul(fm, xo))
      call dsyev('V', 'U', n, fp, n, ev, work, lw, info)
      cm = matmul(xo, fp)
   end subroutine diag

end module link_gfortran_suite

program link_gfortran
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use link_gfortran_suite, only: collect_link_gfortran, release_link_gfortran
   implicit none
   integer :: stat
   stat = 0
   if (command_argument_count() < 1) then
      write (error_unit, "(a)") "usage: link_gfortran <xyz>"
      error stop 1
   end if
   call run_testsuite(collect_link_gfortran, error_unit, stat, parallel=.false.)
   call release_link_gfortran()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program link_gfortran
