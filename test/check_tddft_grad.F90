module check_tddft_grad_suite
   !!
   !! Nuclear gradients of TDA (CIS for Hartree-Fock) singlet excited states
   !! through trc_tddft_grad, for several roots at once.
   !!
   !! One (molecule, basis, functional) per process; one SCF, one solve and one
   !! gradient shared by the cases (built by the first, released by the last).
   !! Arguments: xyz basis.json molname basname funcname nroots probe.bin [fd].
   !!  1. every root's FULL gradient (nuclear repulsion and ground state
   !!     included) against pyscf's grad.tdrks/tdrhf at terco's own converged
   !!     orbitals, amplitudes and grid (test/pyscf_refs/generate.py,
   !!     gen_tddft_grad: pyscf's amplitudes come from an exact diagonalisation
   !!     of its own TDA matrix, so only the gradient code is compared). Pulay-only
   !!     XC terms, as pyscf's tdrks (`grid_response` off), gate GATE_PYSCF;
   !!  2. chunking: one to five roots per sweep (two is the default, four the cap)
   !!     agree, and the number of sweeps is the number of chunks;
   !!  3. the gradient of a root does not depend on the other roots solved with
   !!     it (nroots = 1 against the same root of the full solve);
   !!  4. refusals: a non-converged result, a meta-GGA, a range-separated
   !!     hybrid, VV10 and a missing grid leave an error and no gradient.
   !! With the argument `fd` the one test run is the finite-difference check:
   !! the analytic gradient of every root against the Richardson-extrapolated
   !! central difference of E_SCF + omega over every nuclear coordinate, the
   !! basis, pair list, grid and SCF rebuilt at each displacement. For a
   !! functional the analytic gradient runs with `grid_response` on (a grid
   !! from `dft_grid_build`, which is the SCF's own, carrying its partition),
   !! since only then is it the exact derivative of terco's energy; the
   !! Pulay-only gradient is reported against the same differences and must
   !! miss them by more than the gate, so the test can tell the two apart.
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_bind_device
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_FAMILY_MGGA
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_solve_rks, trc_tddft_reference_t
   use trc_tddft_grad, only: trc_tddft_grad_options_t, trc_tddft_gradient_rks
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector, ref_read_matrix3
   implicit none
   private
   public :: collect_check_tddft_grad, xyz_arg, bas_json_arg, molname_arg, basname_arg, funcname_arg, nroots_arg, &
             probe_arg, mode_arg

   real(dp), parameter :: ERI_THRESH = 1.0e-14_dp
   real(dp), parameter :: CONV_TOL = 1.0e-10_dp
   real(dp), parameter :: CPKS_TOL = 1.0e-11_dp
   real(dp), parameter :: GATE_PYSCF = 5.0e-9_dp
   !! Gradient against pyscf, max over roots, atoms and components, for every
   !! functional and Hartree-Fock. pyscf's own Z-vector (cphf.solve, a Krylov
   !! solve) stops at a residual near 1e-6 whatever tolerance it is given, which
   !! is 1e-7 off the exact solution, so the references use a dense exact solve
   !! (generate.py); measured agreement is at most 1.5e-9.
   real(dp), parameter :: GATE_CHUNK = 1.0e-11_dp
   !! Chunked against unchunked: the same sums in another order.
   real(dp), parameter :: FD_STEP = 2.5e-3_dp
   real(dp), parameter :: FD_ABSTOL = 1.0e-8_dp
   !! Hartree-Fock: Richardson at h = 2.5e-3 measures 2.6e-9 on H2O2/6-31G (at
   !! h = 0.01 the h^4 term alone is 2.2e-7 there; water, stiffer, gave 1.2e-9).
   real(dp), parameter :: FD_ABSTOL_DFT = 2.0e-7_dp
   !! Kohn-Sham (H2O2 PBE, 2 roots, grid response on, the grid rebuilt at every
   !! displacement): measured 3.6e-8 at h = 2.5e-3; the Pulay-only gradient
   !! misses the same differences by 1.5e-5.

   character(len=256) :: xyz_arg = 'water.xyz'
   character(len=256) :: bas_json_arg = 'basis_sets/6-31g.json'
   character(len=16) :: molname_arg = 'water'
   character(len=16) :: basname_arg = '6-31g'
   character(len=16) :: funcname_arg = 'hf'
   integer :: nroots_arg = 3
   character(len=256) :: probe_arg = 'tddft_grad_probe.bin'
   character(len=16) :: mode_arg = ''

   type :: fixture_t
      logical :: ready = .false.
      logical :: ok = .false.
      !! ready: the objects are built (and must be released); ok: the SCF, the
      !! solve and the gradient all succeeded, so a later test may use them.
      logical :: dft = .false.
      integer :: nocc = 0, natm = 0
      integer, allocatable :: zint(:)
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_xc_grid_t) :: xg
      type(trc_xc_functional_t) :: func
      type(trc_scf_result_t) :: scf
      type(trc_tddft_result_t) :: td
      type(trc_tddft_reference_t) :: ref
      type(trc_tddft_grad_options_t) :: gopts
      real(dp) :: t_grad = 0.0_dp
   end type fixture_t

   type(fixture_t), save :: fx

contains

   subroutine collect_check_tddft_grad(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (trim(mode_arg) == 'fd') then
         testsuite = [new_unittest("analytic gradient of every root against finite differences", test_fd), &
                      new_unittest("release", test_release)]
      else
         testsuite = [ &
                     new_unittest("every root's gradient vs pyscf", test_vs_pyscf), &
                     new_unittest("roots per sweep 1..5 (capped at 4): same gradients, one sweep per chunk", test_chunking), &
                     new_unittest("a root's gradient does not depend on the other roots", test_single_root), &
                     new_unittest("refusals leave an error and no gradient", test_refusals), &
                     new_unittest("release", test_release) &
                     ]
      end if
   end subroutine collect_check_tddft_grad

   subroutine scf_options(sopts)
      type(trc_scf_options_t), intent(out) :: sopts
      ! The gradient is first order in the orbital error, so converge well past
      ! the gate (1e-8), but not to the arithmetic's floor. Incremental Fock builds
      ! are off: their accumulated screening error holds |FDS-SDF| at 1e-11..4e-11
      ! and dE at 1e-11 for H2O2 (|E| = 150), so a threshold of 1e-11 was met
      ! only by luck on one toolchain and never on the CI builds, which ran out of
      ! iterations. Full builds fall smoothly (dE ~ 2e-13). With them, 1e-11 and
      ! 1e-10 leave the finite difference at 8e-10 against 1e-8 (1e-12/1e-11 gave
      ! 9e-10: the residual is the h^4 and screening error, not the SCF).
      sopts%incremental = .false.
      sopts%conv_energy = 1.0e-11_dp
      sopts%conv_density = 1.0e-10_dp
      sopts%conv_diis = 1.0e-10_dp
      sopts%max_iteration = 200
      sopts%grid_level = 3
      if (trim(funcname_arg) /= 'hf') sopts%functional = trim(funcname_arg)
   end subroutine scf_options

   subroutine fixture_build(error)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(trc_scf_options_t) :: sopts
      type(trc_tddft_options_t) :: topts
      integer :: dev, natm
      real(dp), allocatable :: at_r(:, :)
      logical :: ok

      if (fx%ready) then
         call check(error, fx%ok, "the fixture (SCF, TDA solve or gradient) failed in an earlier test")
         return
      end if
      dev = trc_bind_device(0)
      fx%dft = trim(funcname_arg) /= 'hf'
      if (fx%dft) then
         ok = xc_functional_by_name(trim(funcname_arg), fx%func)
         call check(error, ok, "unrecognised functional '"//trim(funcname_arg)//"'")
         if (allocated(error)) return
      end if
      call read_xyz(trim(xyz_arg), natm, fx%zint, at_r)
      fx%natm = natm
      call trc_basis_from_json(trim(bas_json_arg), natm, fx%zint, at_r, fx%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%pl%build(fx%bas, ERI_THRESH)
      call fx%pl%to_device()
      call fx%eri%build(fx%bas, ERI_THRESH)
      call fx%la%init(fx%bas%nao)
      fx%ready = .true.
      fx%nocc = nint(sum(real(fx%zint, dp)))/2

      call scf_options(sopts)
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, sopts, fx%scf)
      call check(error, fx%scf%converged, "SCF did not converge")
      if (allocated(error)) return
      if (fx%dft) then
         call fx%xg%build(size(fx%scf%grid_w), fx%scf%grid_r, fx%scf%grid_w, fx%bas, max_pts=sopts%grid_max_pts)
         call fx%xg%to_device()
         call fx%ref%set(fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1), fx%func)
      else
         call fx%ref%set(fx%nocc, fx%scf%cmo(:, :, 1), fx%scf%eps(:, 1), fx%scf%dmat(:, :, 1))
      end if
      topts%nroots = nroots_arg
      topts%conv_tol = CONV_TOL
      call trc_tddft_solve_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, topts, fx%td, fx%xg)
      call check(error, fx%td%converged, "TDA solve failed: "//fx%td%err%get_message())
      if (allocated(error)) return
      call gradient(fx%gopts, 0, err)
      call check(error, .not. err%has_error(), "gradient failed: "//err%get_message())
      if (allocated(error)) return
      fx%ok = .true.
   end subroutine fixture_build

   subroutine gradient(gopts, roots_per_sweep, err)
      !! The gradient of every root of the fixture's solve.
      type(trc_tddft_grad_options_t), intent(inout) :: gopts
      integer, intent(in) :: roots_per_sweep
      type(error_t), intent(out) :: err
      gopts%roots_per_sweep = roots_per_sweep
      gopts%cpks_tol = CPKS_TOL
      gopts%grid_response = .false.
      call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, fx%td, gopts, fx%xg)
      err = gopts%err
   end subroutine gradient

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: gref(:, :, :)
      real(dp) :: e_ref, worst, d, scale, gate
      integer :: u, nr, unit, k
      character(len=256) :: refname

      call fixture_build(error)
      if (allocated(error)) return
      nr = nroots_arg

      ! The probe pyscf's generator reads: terco's grid, energies and gradients.
      open (newunit=unit, file=trim(probe_arg), access='stream', form='unformatted', status='replace')
      write (unit) int(fx%bas%nao, int64), int(fx%natm, int64), int(fx%nocc, int64), int(nr, int64)
      write (unit) molname_arg, basname_arg, funcname_arg
      write (unit) fx%scf%energy
      write (unit) fx%td%omega
      write (unit) fx%gopts%grad
      ! Hartree-Fock leaves the grid unallocated: write a zero count and no arrays.
      if (allocated(fx%scf%grid_w)) then
         write (unit) int(size(fx%scf%grid_w), int64)
         write (unit) fx%scf%grid_r, fx%scf%grid_w
      else
         write (unit) 0_int64
      end if
      close (unit)

      refname = 'tddft_grad_'//trim(molname_arg)//'_'//trim(basname_arg)//'_'//trim(funcname_arg)//'.ref'
      call ref_require(trim(refname), u)
      allocate (gref(3, fx%natm, nr))
      call ref_read_real(u, 'E_SCF', e_ref)
      call ref_read_matrix3(u, 'GRAD', 3, fx%natm, nr, gref)
      call ref_close(u)

      call check(error, abs(fx%scf%energy - e_ref) < 1.0e-9_dp, &
                 "E vs pyscf, |dE| = "//to_string(abs(fx%scf%energy - e_ref), "es10.2"))
      if (allocated(error)) return
      worst = 0.0_dp
      do k = 1, nr
         d = maxval(abs(fx%gopts%grad(:, :, k) - gref(:, :, k)))
         scale = maxval(abs(gref(:, :, k)))
         call global_logger%info("tddft_grad "//trim(molname_arg)//"/"//trim(funcname_arg)//" root "//to_string(k)// &
                                 ": omega "//to_string(fx%td%omega(k), "f12.8")//", max|d grad| vs pyscf "// &
                                 to_string(d, "es10.2")//" (max|grad| "//to_string(scale, "es9.2")//")")
         worst = max(worst, d)
      end do
      gate = GATE_PYSCF
      call global_logger%info("tddft_grad "//trim(molname_arg)//"/"//trim(funcname_arg)//": worst "// &
                              to_string(worst, "es10.2")//" (gate "//to_string(gate, "es8.1")//"), "// &
                              to_string(fx%gopts%nsweep)//" sweep(s), CPKS iterations "// &
                              to_string(maxval(fx%gopts%cpks_iterations))//" max")
      call check(error, worst < gate, "gradient vs pyscf, max|diff| = "//to_string(worst, "es10.2"))
      if (allocated(error)) return
      ! The excitation must matter: a ground-state gradient alone would pass a
      ! check that compared nothing but nuclear repulsion and the ground state.
      call check(error, maxval(abs(fx%gopts%grad(:, :, 1) - fx%gopts%grad(:, :, nr))) > 1.0e-4_dp, &
                 "the roots' gradients should differ")
   end subroutine test_vs_pyscf

   subroutine test_chunking(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_grad_options_t) :: o
      type(error_t) :: err
      real(dp) :: d, worst
      integer :: nr, per, nchunk

      call fixture_build(error)
      if (allocated(error)) return
      nr = nroots_arg
      call check(error, fx%gopts%nsweep == (nr + 1)/2, "the default is two roots per sweep, got "// &
                 to_string(fx%gopts%nsweep)//" sweeps")
      if (allocated(error)) return
      worst = 0.0_dp
      do per = 1, 5
         call gradient(o, per, err)
         call check(error, .not. err%has_error(), "roots_per_sweep = "//to_string(per)//": "//err%get_message())
         if (allocated(error)) return
         nchunk = (nr + min(per, 4) - 1)/min(per, 4)
         call check(error, o%nsweep == nchunk, "roots_per_sweep = "//to_string(per)//" made "//to_string(o%nsweep)// &
                    " sweeps, expected "//to_string(nchunk))
         if (allocated(error)) return
         d = maxval(abs(o%grad - fx%gopts%grad))
         worst = max(worst, d)
         call check(error, d < GATE_CHUNK, "roots_per_sweep = "//to_string(per)//" differs from one sweep by "// &
                    to_string(d, "es10.2"))
         if (allocated(error)) return
      end do
      call global_logger%info("tddft_grad "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                              ": chunked against one sweep, worst "//to_string(worst, "es10.2")//" (gate "// &
                              to_string(GATE_CHUNK, "es8.1")//")")
   end subroutine test_chunking

   subroutine test_single_root(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_result_t) :: one
      type(trc_tddft_grad_options_t) :: o
      real(dp) :: d
      integer :: k
      call fixture_build(error)
      if (allocated(error)) return
      k = nroots_arg
      ! The last root alone: same amplitudes (the solver's own, so the same up
      ! to its convergence), but a Z-vector solved for it by itself.
      one%converged = .true.
      one%nroots = 1
      allocate (one%xamp(size(fx%td%xamp, 1), size(fx%td%xamp, 2), 1), one%omega(1))
      one%xamp(:, :, 1) = fx%td%xamp(:, :, k)
      one%omega(1) = fx%td%omega(k)
      o%grid_response = .false.
      o%cpks_tol = CPKS_TOL
      call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, one, o, fx%xg)
      call check(error, .not. o%err%has_error(), "single root: "//o%err%get_message())
      if (allocated(error)) return
      d = maxval(abs(o%grad(:, :, 1) - fx%gopts%grad(:, :, k)))
      call global_logger%info("tddft_grad "//trim(molname_arg)//"/"//trim(funcname_arg)//": root "//to_string(k)// &
                              " alone against the batch, "//to_string(d, "es10.2")//" (gate "// &
                              to_string(GATE_CHUNK, "es8.1")//")")
      call check(error, d < GATE_CHUNK, "root "//to_string(k)//" alone differs from the batch by "//to_string(d, "es10.2"))
   end subroutine test_single_root

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_result_t) :: bad
      type(trc_tddft_reference_t) :: r
      type(trc_tddft_grad_options_t) :: o
      type(trc_xc_functional_t) :: f_mgga, f_rsh, f_vv10
      logical :: ok

      call fixture_build(error)
      if (allocated(error)) return

      ! A solve that did not converge.
      bad = fx%td
      bad%converged = .false.
      call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, bad, o, fx%xg)
      call check(error, o%err%has_error() .and. .not. allocated(o%grad), "an unconverged result must be refused")
      if (allocated(error)) return
      call check(error, index(o%err%get_message(), "converged") > 0, "message: "//o%err%get_message())
      if (allocated(error)) return

      ok = xc_functional_by_name('pbe', f_mgga)
      f_mgga%family = XC_FAMILY_MGGA
      ok = xc_functional_by_name('pbe', f_rsh)
      f_rsh%omega = 0.3_dp
      f_rsh%exx_lr = 0.8_dp
      ok = xc_functional_by_name('pbe', f_vv10)
      f_vv10%nlc_b = 6.0_dp
      f_vv10%nlc_c = 0.01_dp
      call check(error, ok, "pbe is registered")
      if (allocated(error)) return

      call refuse_with(f_mgga, "meta-GGA", error)
      if (allocated(error)) return
      call refuse_with(f_rsh, "range-separated", error)
      if (allocated(error)) return
      call refuse_with(f_vv10, "VV10", error)
      if (allocated(error)) return

      ! The amplitudes do not match the reference's orbital counts.
      bad = fx%td
      deallocate (bad%xamp)
      allocate (bad%xamp(1, 1, 1))
      bad%xamp = 0.0_dp
      o%cpks_iterations = [7]
      call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, bad, o, fx%xg)
      call check(error, o%err%has_error() .and. index(o%err%get_message(), "do not match") > 0, &
                 "an amplitude shape mismatch must be refused: "//o%err%get_message())
      if (allocated(error)) return
      call check(error, .not. allocated(o%grad) .and. .not. allocated(o%cpks_iterations), &
                 "a refusal must not leave out fields from an earlier call")
      if (allocated(error)) return

      call refuse_spherical(error)
      if (allocated(error)) return

      ! Kohn-Sham without a grid.
      if (fx%dft) then
         r = fx%ref
         call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, r, fx%td, o)
         call check(error, o%err%has_error() .and. index(o%err%get_message(), "grid") > 0, &
                    "a functional without a grid must be refused")
      end if
   end subroutine test_refusals

   subroutine refuse_spherical(error)
      !! Orbitals left in the spherical space (cc-pVDZ loads spherical, nao_sph < nao
      !! rows) are refused: the gradient takes Cartesian AO coefficients, which
      !! trc_tddft_run folds from a spherical SCF.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: sph
      type(trc_pairlist_t) :: spl
      type(trc_eri_t) :: seri
      type(trc_tddft_reference_t) :: r
      type(trc_tddft_result_t) :: td
      type(trc_tddft_grad_options_t) :: o
      type(error_t) :: err
      real(dp), allocatable :: cmo(:, :), eps(:), dmat(:, :), at_r(:, :)
      integer, allocatable :: zint(:)
      integer :: n, natm, i

      call read_xyz('water.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/cc-pvdz.json', natm, zint, at_r, sph, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, sph%spherical, "cc-pvdz must load spherical")
      if (allocated(error)) return
      ! Orbitals still in the spherical space (nao_sph rows): the gradient wants them folded to Cartesian.
      n = sph%nao_sph
      allocate (cmo(n, n), eps(n), dmat(n, n))
      cmo = 0.0_dp
      do i = 1, n
         cmo(i, i) = 1.0_dp
         eps(i) = real(i, dp)
      end do
      dmat = 0.0_dp
      call r%set(5, cmo, eps, dmat)
      td%converged = .true.
      allocate (td%xamp(5, n - 5, 1), td%omega(1))
      td%xamp = 0.0_dp
      td%omega = 1.0_dp
      call trc_tddft_gradient_rks(sph, spl, seri, fx%la, r, td, o)
      call check(error, o%err%has_error() .and. index(o%err%get_message(), "Cartesian basis") > 0 .and. &
                 .not. allocated(o%grad), "spherical-space orbitals must be refused: "//o%err%get_message())
      call sph%release()
   end subroutine refuse_spherical

   subroutine refuse_with(func, what, error)
      !! The gradient with `func` in place of the reference's functional is refused.
      type(trc_xc_functional_t), intent(in) :: func
      character(len=*), intent(in) :: what
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_reference_t) :: r
      type(trc_tddft_grad_options_t) :: o
      r = fx%ref
      r%dft = .true.
      r%func = func
      call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, r, fx%td, o, fx%xg)
      call check(error, o%err%has_error() .and. .not. allocated(o%grad), what//" must be refused")
      if (allocated(error)) return
      call check(error, index(o%err%get_message(), what) > 0, what//" message: "//o%err%get_message())
   end subroutine refuse_with

   subroutine excited_energies(at_r, nr, energy, error)
      !! E_SCF + omega_k of the lowest nr Hartree-Fock roots with the nuclei at at_r.
      real(dp), intent(in) :: at_r(:, :)
      integer, intent(in) :: nr
      real(dp), intent(out) :: energy(nr)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_linalg_t) :: la
      type(trc_scf_result_t) :: scf
      type(trc_scf_options_t) :: sopts
      type(trc_tddft_options_t) :: topts
      type(trc_tddft_result_t) :: td
      type(trc_tddft_reference_t) :: ref
      type(trc_xc_grid_t) :: xgrid
      type(error_t) :: err

      energy = 0.0_dp
      call trc_basis_from_json(trim(bas_json_arg), fx%natm, fx%zint, at_r, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call pl%build(bas, ERI_THRESH)
      call pl%to_device()
      call eri%build(bas, ERI_THRESH)
      call la%init(bas%nao)
      call scf_options(sopts)
      call trc_scf_run(bas, fx%nocc, fx%nocc, sopts, scf)
      if (scf%converged) then
         if (fx%dft) then
            call xgrid%build(size(scf%grid_w), scf%grid_r, scf%grid_w, bas, max_pts=sopts%grid_max_pts)
            call xgrid%to_device()
            call ref%set(fx%nocc, scf%cmo(:, :, 1), scf%eps(:, 1), scf%dmat(:, :, 1), fx%func)
         else
            call ref%set(fx%nocc, scf%cmo(:, :, 1), scf%eps(:, 1), scf%dmat(:, :, 1))
         end if
         topts%nroots = nr
         topts%conv_tol = 1.0e-10_dp
         call trc_tddft_solve_rks(bas, pl, eri, la, ref, topts, td, xgrid)
         if (td%converged) energy = scf%energy + td%omega
         if (fx%dft) call xgrid%release()
      end if
      call check(error, scf%converged, "displaced SCF did not converge: "//to_string(scf%iterations)// &
                 " iterations, "//trim(scf%message))
      if (.not. allocated(error)) call check(error, td%converged, "displaced TDA solve: "//td%err%get_message())
      call eri%release()
      call la%release()
      call pl%release()
      call bas%release()
   end subroutine excited_energies

   subroutine test_fd(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: rp(:, :), ep(:), em(:), fd(:, :, :), coarse(:)
      type(trc_tddft_grad_options_t) :: exact, pulay
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: gxg
      type(error_t) :: err
      real(dp) :: worst_abs, worst_pulay, gate
      integer, allocatable :: fd_atoms(:)
      integer :: nr, atom, dir, natom_fd, k

      call fixture_build(error)
      if (allocated(error)) return
      nr = nroots_arg
      ! The analytic gradient that is the exact derivative of the energy: for a
      ! functional, the grid response on a grid that carries its partition.
      exact%cpks_tol = CPKS_TOL
      pulay = exact
      if (fx%dft) then
         call dft_grid_build(fx%bas%at_r, fx%bas%at_elem, grid, err, level=3)
         call check(error, .not. err%has_error(), "grid: "//err%get_message())
         if (allocated(error)) return
         call check(error, grid%n_points == size(fx%scf%grid_w), &
                    "the FD grid must be the SCF's own, "//to_string(grid%n_points)//" against "// &
                    to_string(size(fx%scf%grid_w))//" points")
         if (allocated(error)) return
         call gxg%build_from_grid(grid, fx%bas, max_pts=512)
         call gxg%to_device()
         exact%grid_response = .true.
         exact%partition = grid%partition
         call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, fx%td, exact, gxg)
         call trc_tddft_gradient_rks(fx%bas, fx%pl, fx%eri, fx%la, fx%ref, fx%td, pulay, gxg)
         call gxg%release()
         call check(error, .not. exact%err%has_error() .and. .not. pulay%err%has_error(), "gradient failed")
         if (allocated(error)) return
      else
         exact%grad = fx%gopts%grad
      end if
      gate = merge(FD_ABSTOL_DFT, FD_ABSTOL, fx%dft)

      allocate (ep(nr), em(nr), coarse(nr), fd(3, fx%natm, nr))
      rp = fx%bas%at_r
      ! Kohn-Sham differences are expensive (a grid and an SCF per point): the
      ! first and the last atom only, which for H2O2 (an O and an H) still
      ! gives every Cartesian component a non-zero gradient.
      natom_fd = fx%natm
      if (fx%dft) natom_fd = 2
      allocate (fd_atoms(natom_fd))
      fd_atoms = [(k, k=1, natom_fd)]
      if (fx%dft) fd_atoms = [1, fx%natm]
      ! Richardson: (4 D(h/2) - D(h))/3 removes the h^2 truncation term of the
      ! central difference D(h), leaving h^4.
      fd = exact%grad
      do k = 1, natom_fd
         atom = fd_atoms(k)
         do dir = 1, 3
            rp(dir, atom) = fx%bas%at_r(dir, atom) + FD_STEP
            call excited_energies(rp, nr, ep, error)
            if (allocated(error)) return
            rp(dir, atom) = fx%bas%at_r(dir, atom) - FD_STEP
            call excited_energies(rp, nr, em, error)
            if (allocated(error)) return
            coarse = (ep - em)/(2.0_dp*FD_STEP)
            rp(dir, atom) = fx%bas%at_r(dir, atom) + 0.5_dp*FD_STEP
            call excited_energies(rp, nr, ep, error)
            if (allocated(error)) return
            rp(dir, atom) = fx%bas%at_r(dir, atom) - 0.5_dp*FD_STEP
            call excited_energies(rp, nr, em, error)
            if (allocated(error)) return
            rp(dir, atom) = fx%bas%at_r(dir, atom)
            fd(dir, atom, :) = (4.0_dp*(ep - em)/FD_STEP - coarse)/3.0_dp
         end do
      end do
      worst_abs = maxval(abs(exact%grad - fd))
      call global_logger%info("tddft_grad "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                              " finite difference (Richardson, h = "//to_string(FD_STEP, "es8.1")//"): worst abs "// &
                              to_string(worst_abs, "es10.2")//" (gate "//to_string(gate, "es8.1")//"), max|grad| "// &
                              to_string(maxval(abs(fd)), "es9.2")//", smallest component "// &
                              to_string(minval(maxval(maxval(abs(fd(:, fd_atoms, :)), dim=3), dim=2)), "es9.2"))
      call check(error, worst_abs < gate, "analytic gradient vs finite difference, worst abs "// &
                 to_string(worst_abs, "es10.2"))
      if (allocated(error)) return
      ! Every Cartesian component carries a gradient: no zero one by symmetry.
      call check(error, minval(maxval(maxval(abs(fd(:, fd_atoms, :)), dim=3), dim=2)) > 1.0e-3_dp, &
                 "a Cartesian component of the gradient vanishes: the test molecule is planar")
      if (allocated(error)) return
      if (fx%dft) then
         worst_pulay = maxval(abs(pulay%grad - fd))
         call global_logger%info("tddft_grad "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                                 ": Pulay-only (grid_response off) misses the differences by "// &
                                 to_string(worst_pulay, "es10.2"))
         call check(error, worst_pulay > 10.0_dp*gate, "the Pulay-only gradient should miss the FD by more than "// &
                    "ten gates, so the grid response is being tested, got "//to_string(worst_pulay, "es10.2"))
      end if
   end subroutine test_fd

   subroutine test_release(error)
      type(error_type), allocatable, intent(out) :: error
      if (fx%ready) then
         call fx%eri%release()
         call fx%la%release()
         if (fx%dft) call fx%xg%release()
         call fx%pl%release()
         call fx%bas%release()
         fx%ready = .false.
         fx%ok = .false.
      end if
      if (.false.) allocate (error)
   end subroutine test_release

end module check_tddft_grad_suite

program check_tddft_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_tddft_grad_suite, only: collect_check_tddft_grad, xyz_arg, bas_json_arg, molname_arg, basname_arg, &
                                      funcname_arg, nroots_arg, probe_arg, mode_arg
   implicit none
   integer :: stat
   character(len=16) :: argument

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, bas_json_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, molname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, basname_arg)
   if (command_argument_count() >= 5) call get_command_argument(5, funcname_arg)
   if (command_argument_count() >= 6) then
      call get_command_argument(6, argument)
      read (argument, *) nroots_arg
   end if
   if (command_argument_count() >= 7) call get_command_argument(7, probe_arg)
   if (command_argument_count() >= 8) call get_command_argument(8, mode_arg)

   stat = 0
   call run_testsuite(collect_check_tddft_grad, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft_grad
