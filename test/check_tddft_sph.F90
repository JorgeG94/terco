module check_tddft_sph_suite
   !!
   !! TDA and full RPA excited states of a SPHERICAL orbital basis (cc-pVDZ loads
   !! spherical) through the Fortran entry trc_tddft_run, against pyscf run with
   !! mol.cart = False (test/pyscf_refs/generate.py, `tddft_sph_*` cases):
   !!  1. TDA excitation energies and oscillator strengths of the lowest three
   !!     singlets;
   !!  2. RPA energies, oscillator strengths and ||Y||^2 (pyscf's (X, Y) are
   !!     normalised 2 (X.X - Y.Y) = 1, so ||Y||^2(terco) = 2 ||Y||^2(pyscf));
   !!  3. (argument `grad`) the Pulay-only TDA and RPA gradient of every one of
   !!     the three roots, against grad.tdrhf / grad.tdrks (pyscf's amplitudes
   !!     are the exact eigenvectors of its own TDA / RPA matrix and its Z-vector
   !!     is a dense solve, so only the gradient formulas are compared).
   !! The kernels run over Cartesian AOs and the run folds the spherical orbitals
   !! to Cartesian coefficients once (trc_tddft_run's SPHERICAL BASES), so what
   !! this pins is that fold and the nmo < nao bookkeeping of the solver, the
   !! Z-vector and the gradient. Tolerances are those of the Cartesian checks
   !! (check_tddft_basis, check_tddft_grad). One (molecule, functional) per
   !! process; the SCF uses pyscf's convergence and writes terco's numbers and its
   !! grid to the probe file for generate.py.
   !!
   !! Arguments: xyz molname funcname probe.bin [grad]. Run from test/. With `grad` and a Kohn-Sham functional it
   !! also checks the grid response (test_grid_response).
   !!
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_tddft_run, only: trc_tddft_run_options_t, trc_tddft_run_result_t, trc_tddft_run_solve, &
                            trc_tddft_run_gradient
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_vector, ref_read_matrix3
   implicit none
   private
   public :: collect_check_tddft_sph
   public :: xyz_arg, molname_arg, funcname_arg, probe_arg, grad_arg

   ! TOL_YY and GATE_GRAD are the Cartesian checks' gates against pyscf (check_tddft_basis, check_tddft_grad).
   real(dp), parameter :: TOL_OMEGA = 1.0e-8_dp
   real(dp), parameter :: TOL_OSC = 1.0e-7_dp
   real(dp), parameter :: TOL_YY = 1.0e-6_dp
   real(dp), parameter :: GATE_GRAD = 5.0e-9_dp
   ! Host-measured over the six cases: max|d omega| 2.8e-10 and max|d f| 2.1e-9 (the gates keep 35x and 50x for
   ! the device's 1e-12 atomics noise); the grid response moves the water gradient by 4.5e-6 to 5.5e-6 and the atom
   ! sum is 2e-13 with it and 4e-6 to 5e-6 without.
   real(dp), parameter :: GRID_MOVE_MIN = 1.0e-6_dp
   real(dp), parameter :: GRID_SUM_MAX = 1.0e-9_dp
   integer, parameter :: NROOTS = 3
   character(len=*), parameter :: BASIS_JSON = "basis_sets/cc-pvdz.json"

   character(len=256) :: xyz_arg = "water.xyz"
   character(len=16) :: molname_arg = "water"
   character(len=16) :: funcname_arg = "hf"
   character(len=256) :: probe_arg = "tddft_sph_probe.bin"
   logical :: grad_arg = .false.

   type :: fixture_t
      logical :: ready = .false.   ! the basis exists: test_release frees it
      logical :: solved = .false.  ! the SCF, both solves and the gradients all succeeded
      integer :: nocc = 0, natm = 0
      type(trc_basis_t) :: bas
      type(trc_linalg_t) :: la
      type(trc_scf_options_t) :: sopts
      type(trc_scf_result_t) :: scf
      type(trc_tddft_run_result_t) :: tda, rpa
      real(dp), allocatable :: g_tda(:, :, :), g_rpa(:, :, :)
   end type fixture_t

   type(fixture_t), save :: fx

contains

   subroutine collect_check_tddft_sph(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (grad_arg) then
         testsuite = [ &
                     new_unittest("spherical TDA energies vs pyscf", test_tda), &
                     new_unittest("spherical RPA energies and ||Y||^2 vs pyscf", test_rpa), &
                     new_unittest("spherical TDA gradients vs pyscf (Pulay only)", test_grad_tda), &
                     new_unittest("spherical RPA gradients vs pyscf (Pulay only)", test_grad_rpa), &
                     new_unittest("spherical grid response (Kohn-Sham)", test_grid_response), &
                     new_unittest("release", test_release) &
                     ]
      else
         testsuite = [ &
                     new_unittest("spherical TDA energies vs pyscf", test_tda), &
                     new_unittest("spherical RPA energies and ||Y||^2 vs pyscf", test_rpa), &
                     new_unittest("release", test_release) &
                     ]
      end if
   end subroutine collect_check_tddft_sph

   function ref_name(prefix) result(name)
      character(len=*), intent(in) :: prefix
      character(len=256) :: name
      name = prefix//trim(molname_arg)//"_cc-pvdz_"//trim(funcname_arg)//".ref"
   end function ref_name

   subroutine fixture_build(error)
      !! The spherical basis and one tightly converged SCF, then the TDA and RPA solves (and the gradients).
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(trc_tddft_run_options_t) :: opts
      real(dp), allocatable :: at_r(:, :)
      integer, allocatable :: zint(:)
      integer :: dev, k

      if (fx%solved) return
      call check(error, .not. fx%ready, "fixture setup failed in an earlier test; nothing to compare")
      if (allocated(error)) return
      dev = trc_bind_device(0)
      call read_xyz(trim(xyz_arg), fx%natm, zint, at_r)
      ! Not force_cart: cc-pVDZ is a spherical basis, which is the point.
      call trc_basis_from_json(BASIS_JSON, fx%natm, zint, at_r, fx%bas, err)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, fx%bas%spherical .and. fx%bas%nao_sph < fx%bas%nao, "cc-pvdz must load spherical")
      if (allocated(error)) return
      call fx%bas%to_device()
      call fx%la%init(fx%bas%nao_sph)
      fx%ready = .true.
      fx%nocc = nint(sum(real(zint, dp)))/2
      fx%sopts%conv_energy = 1.0e-12_dp
      fx%sopts%conv_density = 1.0e-10_dp
      fx%sopts%conv_diis = 1.0e-11_dp
      fx%sopts%max_iteration = 200
      fx%sopts%grid_level = 3
      if (trim(funcname_arg) /= "hf") fx%sopts%functional = trim(funcname_arg)
      call trc_scf_run(fx%bas, fx%nocc, fx%nocc, fx%sopts, fx%scf)
      call check(error, fx%scf%converged, "SCF did not converge")
      if (allocated(error)) return

      opts%nroots = NROOTS
      opts%conv_tol = 1.0e-10_dp
      opts%cpks_tol = 1.0e-11_dp
      opts%grid_response = .false.
      opts%tda = .true.
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, fx%tda)
      call check(error, fx%tda%converged, "TDA solve failed: "//fx%tda%err%get_message())
      if (allocated(error)) return
      opts%tda = .false.
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, fx%rpa)
      call check(error, fx%rpa%converged, "RPA solve failed: "//fx%rpa%err%get_message())
      if (allocated(error)) return

      if (grad_arg) then
         allocate (fx%g_tda(3, fx%natm, NROOTS), fx%g_rpa(3, fx%natm, NROOTS))
         do k = 1, NROOTS
            opts%target_root = k
            opts%tda = .true.
            call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, fx%tda)
            call check(error, .not. fx%tda%err%has_error(), "TDA gradient of root "//to_string(k)//": "// &
                       fx%tda%err%get_message())
            if (allocated(error)) return
            fx%g_tda(:, :, k) = fx%tda%grad
            opts%tda = .false.
            call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, fx%rpa)
            call check(error, .not. fx%rpa%err%has_error(), "RPA gradient of root "//to_string(k)//": "// &
                       fx%rpa%err%get_message())
            if (allocated(error)) return
            fx%g_rpa(:, :, k) = fx%rpa%grad
         end do
      end if
      fx%solved = .true.
      call write_probe()
   end subroutine fixture_build

   subroutine write_probe()
      !! Terco's numbers and its grid, for generate.py's tddft_sph cases.
      integer :: unit, ngrad
      character(len=16) :: bas_tag
      bas_tag = "cc-pvdz"
      ngrad = 0
      if (grad_arg) ngrad = NROOTS
      open (newunit=unit, file=trim(probe_arg), access="stream", form="unformatted", status="replace")
      write (unit) int(fx%bas%nao_sph, int64), int(fx%natm, int64), int(fx%nocc, int64), int(NROOTS, int64)
      write (unit) molname_arg, bas_tag, funcname_arg
      write (unit) fx%scf%energy
      write (unit) fx%tda%omega
      write (unit) fx%tda%osc
      write (unit) fx%rpa%omega
      write (unit) fx%rpa%osc
      write (unit) fx%rpa%ynorm
      write (unit) int(ngrad, int64)
      if (grad_arg) write (unit) fx%g_tda, fx%g_rpa
      ! Hartree-Fock has no grid; the generator reads zero points.
      if (allocated(fx%scf%grid_w)) then
         write (unit) int(size(fx%scf%grid_w), int64)
         write (unit) fx%scf%grid_r, fx%scf%grid_w
      else
         write (unit) 0_int64
      end if
      close (unit)
   end subroutine write_probe

   subroutine compare_roots(error, res, prefix, with_y, tag)
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_run_result_t), intent(in) :: res
      character(len=*), intent(in) :: prefix, tag
      logical, intent(in) :: with_y
      real(dp), allocatable :: om_ref(:), osc_ref(:), yy_ref(:)
      real(dp) :: dom, dosc, dyy, e_ref
      integer :: u

      call fixture_build(error)
      if (allocated(error)) return
      call ref_require(ref_name(prefix), u)
      allocate (om_ref(NROOTS), osc_ref(NROOTS), yy_ref(NROOTS))
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_vector(u, "OMEGA", NROOTS, om_ref)
      call ref_read_vector(u, "OSC", NROOTS, osc_ref)
      if (with_y) call ref_read_vector(u, "YY", NROOTS, yy_ref)
      call ref_close(u)
      dom = maxval(abs(res%omega - om_ref))
      dosc = maxval(abs(res%osc - osc_ref))
      call global_logger%info("tddft sph "//tag//" "//trim(molname_arg)//"/"//trim(funcname_arg)//": max|d omega| "// &
                              to_string(dom, "es10.2")//", max|d f| "//to_string(dosc, "es10.2"))
      call check(error, abs(fx%scf%energy - e_ref) < 1.0e-9_dp, "E_SCF vs pyscf")
      if (allocated(error)) return
      call check(error, dom < TOL_OMEGA, "omega vs pyscf, "//to_string(dom, "es10.2"))
      if (allocated(error)) return
      call check(error, dosc < TOL_OSC, "oscillator strength vs pyscf, "//to_string(dosc, "es10.2"))
      if (allocated(error)) return
      if (with_y) then
         dyy = maxval(abs(res%ynorm - yy_ref))
         call global_logger%info("tddft sph "//tag//": max|d ||Y||^2| "//to_string(dyy, "es10.2"))
         call check(error, dyy < TOL_YY, "||Y||^2 vs pyscf, "//to_string(dyy, "es10.2"))
      else
         call check(error, maxval(res%ynorm) == 0.0_dp, "the TDA has Y = 0")
      end if
   end subroutine compare_roots

   subroutine test_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call fixture_build(error)
      if (allocated(error)) return
      call compare_roots(error, fx%tda, "tddft_sph_", .false., "TDA")
   end subroutine test_tda

   subroutine test_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call fixture_build(error)
      if (allocated(error)) return
      call compare_roots(error, fx%rpa, "tddft_sph_rpa_", .true., "RPA")
   end subroutine test_rpa

   subroutine compare_gradients(error, prefix, g, tag)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: prefix, tag
      real(dp), intent(in) :: g(:, :, :)
      real(dp), allocatable :: gref(:, :, :)
      real(dp) :: d, e_ref
      integer :: u, k

      call fixture_build(error)
      if (allocated(error)) return
      allocate (gref(3, fx%natm, NROOTS))
      call ref_require(ref_name(prefix), u)
      call ref_read_real(u, "E_SCF", e_ref)
      call ref_read_matrix3(u, "GRAD", 3, fx%natm, NROOTS, gref)
      call ref_close(u)
      do k = 1, NROOTS
         d = maxval(abs(g(:, :, k) - gref(:, :, k)))
         call global_logger%info("tddft sph "//tag//" gradient "//trim(molname_arg)//"/"//trim(funcname_arg)//", root "// &
                                 to_string(k)//": max|d grad| vs pyscf "//to_string(d, "es10.2")//" (gate "// &
                                 to_string(GATE_GRAD, "es8.1")//")")
         call check(error, d < GATE_GRAD, tag//" gradient of root "//to_string(k)//" vs pyscf, "//to_string(d, "es10.2"))
         if (allocated(error)) return
      end do
   end subroutine compare_gradients

   subroutine test_grad_tda(error)
      type(error_type), allocatable, intent(out) :: error
      call fixture_build(error)
      if (allocated(error)) return
      call compare_gradients(error, "tddft_sph_grad_", fx%g_tda, "TDA")
   end subroutine test_grad_tda

   subroutine test_grad_rpa(error)
      type(error_type), allocatable, intent(out) :: error
      call fixture_build(error)
      if (allocated(error)) return
      call compare_gradients(error, "tddft_sph_rpa_grad_", fx%g_rpa, "RPA")
   end subroutine test_grad_rpa

   subroutine test_grid_response(error)
      !! The grid response of a spherical basis, which pyscf has no reference for (its tdrks has none): the TDA
      !! gradient of root 1 with it on, against the Pulay-only one of the fixture. The grid moves with the nuclei,
      !! so (a) the response changes the gradient by the grid's own derivative, well above the solver noise, and
      !! (b) the full gradient is translation invariant (its atom sum vanishes) while the Pulay-only one is not.
      !! Hartree-Fock has no grid: the test passes trivially.
      type(error_type), allocatable, intent(out) :: error
      type(trc_tddft_run_options_t) :: opts
      type(trc_tddft_run_result_t) :: res
      real(dp) :: diff, sum_full, sum_pulay

      if (trim(funcname_arg) == "hf") return
      call fixture_build(error)
      if (allocated(error)) return
      opts%nroots = 1
      opts%target_root = 1
      opts%tda = .true.
      opts%conv_tol = 1.0e-10_dp
      opts%cpks_tol = 1.0e-11_dp
      call check(error, opts%grid_response, "grid_response is on by default")
      if (allocated(error)) return
      call trc_tddft_run_solve(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call check(error, res%converged, "solve failed: "//res%err%get_message())
      if (allocated(error)) return
      call trc_tddft_run_gradient(fx%bas, fx%la, fx%scf, fx%nocc, fx%sopts, opts, res)
      call check(error, .not. res%err%has_error(), "gradient: "//res%err%get_message())
      if (allocated(error)) return
      diff = maxval(abs(res%grad - fx%g_tda(:, :, 1)))
      sum_full = maxval(abs(sum(res%grad, dim=2)))
      sum_pulay = maxval(abs(sum(fx%g_tda(:, :, 1), dim=2)))
      call global_logger%info("tddft sph grid response "//trim(molname_arg)//"/"//trim(funcname_arg)// &
                              ": moves the gradient by "//to_string(diff, "es10.2")//"; atom sum "// &
                              to_string(sum_full, "es10.2")//" with it, "//to_string(sum_pulay, "es10.2")//" without")
      call check(error, diff > GRID_MOVE_MIN, "the grid response must move the gradient, got "//to_string(diff, "es10.2"))
      if (allocated(error)) return
      call check(error, sum_full < GRID_SUM_MAX, "the full gradient must sum to zero over the atoms, got "// &
                 to_string(sum_full, "es10.2"))
      if (allocated(error)) return
      call check(error, sum_pulay > 10.0_dp*sum_full, "the Pulay-only gradient must fail the atom sum by far more")
   end subroutine test_grid_response

   subroutine test_release(error)
      type(error_type), allocatable, intent(out) :: error
      if (fx%ready) then
         call fx%la%release()
         call fx%bas%release()
         fx%ready = .false.
         fx%solved = .false.
      end if
      call check(error, .not. fx%ready, "released")
   end subroutine test_release

end module check_tddft_sph_suite

program check_tddft_sph
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_tddft_sph_suite, only: collect_check_tddft_sph, xyz_arg, molname_arg, funcname_arg, probe_arg, grad_arg
   implicit none
   integer :: stat
   character(len=16) :: grad_text

   if (command_argument_count() >= 1) call get_command_argument(1, xyz_arg)
   if (command_argument_count() >= 2) call get_command_argument(2, molname_arg)
   if (command_argument_count() >= 3) call get_command_argument(3, funcname_arg)
   if (command_argument_count() >= 4) call get_command_argument(4, probe_arg)
   if (command_argument_count() >= 5) then
      call get_command_argument(5, grad_text)
      if (trim(grad_text) /= "grad") then
         write (error_unit, "(a)") "check_tddft_sph: the fifth argument must be 'grad', got '"//trim(grad_text)//"'"
         error stop 2
      end if
      grad_arg = .true.
   end if
   stat = 0
   call run_testsuite(collect_check_tddft_sph, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_tddft_sph
