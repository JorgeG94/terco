module check_fock_nosym_many_suite
   !!
   !! eri%fock_nosym_many: J and K of several NON-SYMMETRIC densities in one
   !! pass over the integrals, the matvec of a TDDFT (TDA/RPA) Davidson solve,
   !! whose transition densities have an antisymmetric part and need the
   !! enumerated digestion (see eri%fock_nosym) for every root at once.
   !!
   !! Four checks, on water/6-31G* by default (argv: xyz, basis json, probe
   !! dump, reference case; test/CMakeLists.txt):
   !!
   !!  1. a batch equals the same densities through eri%fock_nosym one at a
   !!     time, across several j_scale/k_scale mixes;
   !!  2. a purely antisymmetric density has J = 0 (the integral is symmetric
   !!     in the pair it is contracted with) and an antisymmetric K, the part
   !!     a folded build cannot represent;
   !!  3. J, K and J - K/2 of three non-symmetric densities against pyscf's
   !!     get_jk(dm, hermi=0), committed as test/pyscf_refs/fock_nosym_many_*.ref
   !!     (generate.py's gen_fock_nosym_many; the densities are the closed form
   !!     in nonsym_densities, so the reference never depends on terco's own
   !!     output);
   !!  4. the same agreement with the batch at omega > 0.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_eri, only: trc_eri_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_matrix, ref_read_matrix3
   implicit none
   private
   public :: collect_fock_nosym_many_tests, g_xyzfile, g_basis_json, g_dump, g_case

   real(dp), parameter :: THRESH = 1.0e-14_dp
   real(dp), parameter :: OMEGA_TEST = 0.3_dp
   real(dp), parameter :: TOL_BATCH = 1.0e-13_dp
      !! Batch against single-density calls. The kernel is the same code and
      !! the same quartet order, so the two differ only by the order atomic
      !! updates land in; measured 0 on the host build.
   real(dp), parameter :: TOL_ZERO = 1.0e-13_dp
   real(dp), parameter :: TOL_PYSCF = 1.0e-9_dp
      !! The gate every other pyscf Fock comparison uses (erf_fock, rys_fock).
   integer, parameter :: NPROBE = 3
      !! Densities in the pyscf reference.
   integer, parameter :: NBATCH = 4
      !! Densities in the self-consistency checks.

   character(len=256) :: g_xyzfile = "water.xyz"
   character(len=256) :: g_basis_json = "basis_sets/6-31gs.json"
   character(len=256) :: g_dump = ""
      !! Where test_vs_pyscf writes the fixed inputs generate.py reads
      !! (empty: nothing written).
   character(len=256) :: g_case = "water_631gs"
      !! Names the reference file, fock_nosym_many_<case>.ref.

contains

   subroutine collect_fock_nosym_many_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("batch equals one fock_nosym call per density", test_batch_vs_single), &
                  new_unittest("antisymmetric density: J vanishes, K is antisymmetric", test_antisymmetric), &
                  new_unittest("J, K and J - K/2 match pyscf get_jk(hermi=0)", test_vs_pyscf), &
                  new_unittest("omega > 0 batch equals single-density fock_nosym", test_omega) &
                  ]
   end subroutine collect_fock_nosym_many_tests

   subroutine nonsym_densities(n, nd, dd)
      !! Dense, deterministic, NOT symmetric densities: the closed form the
      !! pyscf reference is generated from. In a module routine with explicit
      !! shapes, as check_nosym's fill_density (an inline one-liner crashed
      !! nvfortran).
      integer, intent(in) :: n, nd
      real(dp), intent(out) :: dd(n, n, nd)
      integer :: i, j, k
      do k = 1, nd
         do j = 1, n
            do i = 1, n
               dd(i, j, k) = 0.5_dp*cos(real(3*i + 5*j + 7*k, dp)) + 0.3_dp*sin(real(2*i - 3*j + k, dp))
            end do
         end do
      end do
   end subroutine nonsym_densities

   subroutine setup(bas, eri, error)
      type(trc_basis_t), intent(out) :: bas
      type(trc_eri_t), intent(out) :: eri
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: nat
      integer, allocatable :: at_z(:)
      real(dp), allocatable :: r(:, :)
      call read_xyz(trim(g_xyzfile), nat, at_z, r)
      ! prune=.false.: the pyscf reference is built from the file's literal
      ! columns, as in check_erf_fock_vs_pyscf.
      call trc_basis_from_json(trim(g_basis_json), nat, at_z, r, bas, err, &
                               overrides=trc_basis_json_overrides(prune=.false.))
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call eri%build(bas, THRESH)
   end subroutine setup

   subroutine test_batch_vs_single(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      real(dp), allocatable :: d(:, :, :), gm(:, :, :), g1(:, :)
      real(dp) :: worst, scl(2, 3)
      integer :: nao, id, imix
      call setup(bas, eri, error)
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao, NBATCH), gm(nao, nao, NBATCH), g1(nao, nao))
      call nonsym_densities(nao, NBATCH, d)
      call check(error, maxval(abs(d(:, :, 1) - transpose(d(:, :, 1)))) > 0.1_dp, "densities must be non-symmetric")
      if (allocated(error)) return
      scl = reshape([1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3])
      worst = 0.0_dp
      do imix = 1, 3
         call eri%fock_nosym_many(bas, NBATCH, d, gm, j_scale=scl(1, imix), k_scale=scl(2, imix))
         do id = 1, NBATCH
            call eri%fock_nosym(bas, d(:, :, id), g1, j_scale=scl(1, imix), k_scale=scl(2, imix), &
                                density_screen=.false.)
            worst = max(worst, maxval(abs(gm(:, :, id) - g1)))
         end do
      end do
      print "(a,es12.4)", "  batch vs single nosym, worst|diff| ", worst
      call check(error, maxval(abs(gm)) > 1.0e-3_dp, "batch result is trivially zero")
      if (allocated(error)) return
      call check(error, worst < TOL_BATCH, "batch vs single nosym, worst|diff| = "// &
                 to_string(worst, "es12.4")//" exceeds "//to_string(TOL_BATCH, "es12.4"))
      call eri%release(); call bas%release()
   end subroutine test_batch_vs_single

   subroutine test_antisymmetric(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      real(dp), allocatable :: d(:, :, :), a(:, :, :), gj(:, :, :), gk(:, :, :)
      real(dp) :: jmax, kasym, kmax
      integer :: nao, id
      call setup(bas, eri, error)
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao, NBATCH), a(nao, nao, NBATCH), gj(nao, nao, NBATCH), gk(nao, nao, NBATCH))
      call nonsym_densities(nao, NBATCH, d)
      do id = 1, NBATCH
         a(:, :, id) = d(:, :, id) - transpose(d(:, :, id))
      end do
      ! g = J with (j, k) = (1, 0), g = -K/2 with (0, 1).
      call eri%fock_nosym_many(bas, NBATCH, a, gj, j_scale=1.0_dp, k_scale=0.0_dp)
      call eri%fock_nosym_many(bas, NBATCH, a, gk, j_scale=0.0_dp, k_scale=1.0_dp)
      jmax = maxval(abs(gj))
      kmax = maxval(abs(gk))
      kasym = 0.0_dp
      do id = 1, NBATCH
         kasym = max(kasym, maxval(abs(gk(:, :, id) + transpose(gk(:, :, id)))))
      end do
      print "(a,es12.4,a,es12.4,a,es12.4)", "  antisymmetric density: max|J| ", jmax, "  max|K| ", kmax, "  max|K+K^T| ", kasym
      call check(error, kmax > 1.0e-2_dp, "K of an antisymmetric density is trivially zero")
      if (allocated(error)) return
      call check(error, jmax < TOL_ZERO, "J of an antisymmetric density, max|J| = "// &
                 to_string(jmax, "es12.4")//" exceeds "//to_string(TOL_ZERO, "es12.4"))
      if (allocated(error)) return
      call check(error, kasym < TOL_ZERO*max(1.0_dp, kmax), "K + K^T for an antisymmetric density = "// &
                 to_string(kasym, "es12.4")//" (max|K| "//to_string(kmax, "es12.4")//")")
      call eri%release(); call bas%release()
   end subroutine test_antisymmetric

   subroutine test_omega(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      real(dp), allocatable :: d(:, :, :), gm(:, :, :), g1(:, :), gc(:, :)
      real(dp) :: worst, shift
      integer :: nao, id
      call setup(bas, eri, error)
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao, NPROBE), gm(nao, nao, NPROBE), g1(nao, nao), gc(nao, nao))
      call nonsym_densities(nao, NPROBE, d)
      call eri%fock_nosym_many(bas, NPROBE, d, gm, omega=OMEGA_TEST)
      worst = 0.0_dp
      shift = 0.0_dp
      do id = 1, NPROBE
         call eri%fock_nosym(bas, d(:, :, id), g1, omega=OMEGA_TEST, density_screen=.false.)
         call eri%fock_nosym(bas, d(:, :, id), gc, density_screen=.false.)
         worst = max(worst, maxval(abs(gm(:, :, id) - g1)))
         shift = max(shift, maxval(abs(g1 - gc)))
      end do
      print "(a,es12.4,a,es12.4)", "  omega batch vs single nosym, worst|diff| ", worst, "  omega shift ", shift
      ! omega must actually have changed the answer, or this compares nothing.
      call check(error, shift > 1.0e-2_dp, "omega did not change the Fock matrix")
      if (allocated(error)) return
      call check(error, worst < TOL_BATCH, "omega batch vs single nosym, worst|diff| = "// &
                 to_string(worst, "es12.4")//" exceeds "//to_string(TOL_BATCH, "es12.4"))
      call eri%release(); call bas%release()
   end subroutine test_omega

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_eri_t) :: eri
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: d(:, :, :), gj(:, :, :), gk(:, :, :), gc(:, :, :)
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      real(dp), allocatable :: s_ref(:, :), j_ref(:, :, :), k_ref(:, :, :), r(:, :)
      integer, allocatable :: at_z(:)
      real(dp) :: ds, wj, wk, wg
      integer :: nao, nat, nat_ref, nao_ref, u, unit
      call setup(bas, eri, error)
      if (allocated(error)) return
      nao = bas%nao
      allocate (d(nao, nao, NPROBE), gj(nao, nao, NPROBE), gk(nao, nao, NPROBE), gc(nao, nao, NPROBE))
      call nonsym_densities(nao, NPROBE, d)
      ! J, then K = -2 g from (0, 1) -- exact, the kernels are linear in the
      ! scales (as check_erf_fock_vs_pyscf) -- then the combined build.
      call eri%fock_nosym_many(bas, NPROBE, d, gj, j_scale=1.0_dp, k_scale=0.0_dp)
      call eri%fock_nosym_many(bas, NPROBE, d, gk, j_scale=0.0_dp, k_scale=1.0_dp)
      gk = -2.0_dp*gk
      call eri%fock_nosym_many(bas, NPROBE, d, gc, j_scale=1.0_dp, k_scale=1.0_dp)

      ! The overlap, for generate.py's per-AO normalisation check.
      allocate (smat(nao, nao), tmat(nao, nao), vmat(nao, nao))
      call pl%build(bas, THRESH)
      call pl%to_device()
      !$acc enter data create(smat, tmat, vmat)
      call trc_1e(bas, pl, smat, tmat, vmat)
      !$acc update self(smat)
      !$acc exit data delete(smat, tmat, vmat)
      call pl%release()

      if (len_trim(g_dump) > 0) then
         call read_xyz(trim(g_xyzfile), nat, at_z, r)
         open (newunit=unit, file=trim(g_dump), access="stream", form="unformatted", status="replace")
         write (unit) int(nat, int64), int(nao, int64)
         write (unit) int(at_z, int64)
         write (unit) r
         write (unit) d
         write (unit) smat
         close (unit)
      end if

      call ref_require("fock_nosym_many_"//trim(g_case)//".ref", u)
      call ref_read_int(u, "NAT", nat_ref)
      call ref_read_int(u, "NAO", nao_ref)
      call check(error, nao_ref == nao, "pyscf reference has "//to_string(nao_ref)//" AOs, this build "//to_string(nao))
      if (allocated(error)) return
      allocate (s_ref(nao, nao), j_ref(nao, nao, NPROBE), k_ref(nao, nao, NPROBE))
      call ref_read_matrix(u, "S", nao, nao, s_ref)
      call ref_read_matrix3(u, "J", nao, nao, NPROBE, j_ref)
      call ref_read_matrix3(u, "K", nao, nao, NPROBE, k_ref)
      call ref_close(u)

      ds = maxval(abs(smat - s_ref))
      wj = maxval(abs(gj - j_ref))
      wk = maxval(abs(gk - k_ref))
      wg = maxval(abs(gc - (j_ref - 0.5_dp*k_ref)))
      print "(a,es12.4,a,es12.4,a,es12.4,a,es12.4)", "  vs pyscf: S ", ds, "  J ", wj, "  K ", wk, "  J-K/2 ", wg
      call check(error, ds < TOL_PYSCF, "overlap vs pyscf, max|diff| = "//to_string(ds, "es12.4"))
      if (allocated(error)) return
      call check(error, wj < TOL_PYSCF, "J vs pyscf get_jk, worst|diff| = "// &
                 to_string(wj, "es12.4")//" exceeds "//to_string(TOL_PYSCF, "es12.4"))
      if (allocated(error)) return
      call check(error, wk < TOL_PYSCF, "K vs pyscf get_jk, worst|diff| = "// &
                 to_string(wk, "es12.4")//" exceeds "//to_string(TOL_PYSCF, "es12.4"))
      if (allocated(error)) return
      call check(error, wg < TOL_PYSCF, "J - K/2 vs pyscf get_jk, worst|diff| = "// &
                 to_string(wg, "es12.4")//" exceeds "//to_string(TOL_PYSCF, "es12.4"))
      call eri%release(); call bas%release()
   end subroutine test_vs_pyscf

end module check_fock_nosym_many_suite

program check_fock_nosym_many
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_fock_nosym_many_suite, only: collect_fock_nosym_many_tests, g_xyzfile, g_basis_json, g_dump, g_case
   implicit none
   integer :: stat
   type(testsuite_type), allocatable :: testsuites(:)

   if (command_argument_count() >= 1) call get_command_argument(1, g_xyzfile)
   if (command_argument_count() >= 2) call get_command_argument(2, g_basis_json)
   if (command_argument_count() >= 3) call get_command_argument(3, g_dump)
   if (command_argument_count() >= 4) call get_command_argument(4, g_case)

   stat = 0
   testsuites = [new_testsuite("check_fock_nosym_many", collect_fock_nosym_many_tests)]
   call run_testsuite(testsuites(1)%collect, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop 1
   end if
end program check_fock_nosym_many
