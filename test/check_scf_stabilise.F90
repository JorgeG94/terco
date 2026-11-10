!
! The two SCF convergence aids for small-gap systems, `mom` (maximum overlap
! occupation) and `level_shift_until` (the level shift kept on through DIIS),
! as a test-drive suite. Every case fails if the aid it names is broken.
!
! 1. scf_driver_mom_weights + scf_driver_mom_permute on a synthetic problem
!    (S = 1, C a Householder reflector, D_old built from chosen columns of C,
!    occ = 2 and occ = 1): the previously occupied columns come first, eps
!    follows, each block keeps its order, nswap counts the columns aufbau
!    would have left empty, an aufbau D_old and nocc = 0 or nocc >= nao move
!    nothing.
!
! 2. Ground-state SCFs with both aids on (water RHF, OH UHF) land on the plain
!    run's energy and on the pyscf energies check_spherical.F90 already holds
!    for these systems (cart=False, conv_tol 1e-11), and water's HOMO and LUMO
!    equal the plain run's: the converged iteration was not level shifted
!    (level_shift_until below conv_diis is clamped to it).
!
! 3. The level shift really stays on through DIIS. Started from a converged
!    density the Fock matrix is block diagonal, so one iteration with the
!    shift on moves exactly the virtual eigenvalues up by the shift mu and
!    leaves the occupied ones; with level_shift_until = 0 nothing moves
!    (DIIS is on from iteration 1 here). Water RHF, OH UHF (both spins), and
!    through the C ABI. The non-converged result says its virtuals are shifted.
!
! 4. MOM keeps a non-aufbau occupation. A converged water density is rebuilt
!    with the HOMO (1b1) and LUMO (4a1, another symmetry) swapped and given as
!    the guess: with mom the SCF converges to that doubly excited stationary
!    point, well above the ground state; without it aufbau falls back to the
!    ground state. Fortran and C ABI. The excited energy has no external
!    reference (a pyscf mom_occ anchor still has to be generated), so only the
!    relations are asserted.
!
! 5. trc_set_scf_stabilise validates its arguments.
!
module check_scf_stabilise_suite
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char, c_loc, c_associated, c_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run, scf_driver_mom_weights, &
                             scf_driver_mom_permute
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, &
                               trc_set_method, trc_set_convergence, trc_set_level_shift, &
                               trc_set_scf_stabilise, trc_set_guess, trc_run_scf, trc_energy, trc_nao, &
                               trc_density, trc_mo_coefficients, trc_mo_energies, trc_message, &
                               trc_converged, TRC_OK, TRC_ERR_BADARG, TRC_ERR_NOCONV, TRC_GUESS_GIVEN
   implicit none
   private
   public :: collect_scf_stabilise

   ! pyscf, cart=False, conv_tol 1e-11; the literals check_spherical.F90 holds
   ! for water.xyz/cc-pvdz RHF and oh.xyz/cc-pvdz UHF.
   real(dp), parameter :: E_WATER_RHF = -76.026772053394_dp
   real(dp), parameter :: E_OH_UHF = -75.393846033475_dp
   real(dp), parameter :: TOL_REF = 1.0e-8_dp
   real(dp), parameter :: MU = 0.5_dp
   integer, parameter :: UHF_NDAMP = 4
   !! The damped start of an unrestricted run (`ndamp` < 0 resolves to 4).
   !! The level shift of the probes. A shifted virtual sits MU above the
   !! unshifted one; the tolerance on that is 1e-5, the unshifted 1e-6.
   real(dp), parameter :: MIN_EXCITATION = 0.3_dp
   !! The doubly excited water state must sit at least this far above the
   !! ground state (measured: see the cases), or it is not a different state.

contains

   subroutine collect_scf_stabilise(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("mom weights+permute: previously occupied columns first", test_mom_order_swaps), &
                  new_unittest("mom weights+permute: aufbau, nocc 0 and nocc >= nao move nothing", &
                               test_mom_order_noops), &
                  new_unittest("mom weights+permute: occ = 1 (unrestricted) density", test_mom_order_occ1), &
                  new_unittest("water RHF ground state: mom and level_shift_until keep energy and HOMO/LUMO", &
                               test_water_ground_state), &
                  new_unittest("OH UHF ground state: mom and level_shift_until keep the energy", test_oh_uhf), &
                  new_unittest("water RHF: level_shift_until keeps the shift on through DIIS", test_shift_rhf), &
                  new_unittest("OH UHF: level_shift_until shifts the virtuals of both spins", test_shift_uhf), &
                  new_unittest("OH UHF: a restart that converges inside the damped start is unshifted", &
                               test_converging_iteration_unshifted), &
                  new_unittest("water RHF: SOSCF's eps are not called level shifted", test_soscf_clears_shift), &
                  new_unittest("water RHF 6-31G* core guess: mom does not lock a poor guess", test_mom_poor_guess), &
                  new_unittest("water RHF: mom keeps a swapped HOMO/LUMO occupation, aufbau does not", &
                               test_mom_holds_excited_rhf), &
                  new_unittest("C ABI trc_set_scf_stabilise: bad arguments are refused", test_capi_validation), &
                  new_unittest("C ABI trc_set_scf_stabilise: level_shift_until reaches the SCF", test_capi_shift), &
                  new_unittest("C ABI trc_set_scf_stabilise: mom reaches the SCF", test_capi_mom) &
                  ]
   end subroutine collect_scf_stabilise

   ! --- 1. the MOM ordering on a synthetic problem ------------------------

   subroutine reflector(n, c)
      !! An orthogonal matrix with every entry non-zero: 1 - 2 v v^T/(v^T v).
      integer, intent(in) :: n
      real(dp), intent(out) :: c(n, n)
      real(dp) :: v(n)
      integer :: i, j
      v = [(real(i, dp), i=1, n)]
      do j = 1, n
         do i = 1, n
            c(i, j) = merge(1.0_dp, 0.0_dp, i == j) - 2.0_dp*v(i)*v(j)/sum(v*v)
         end do
      end do
   end subroutine reflector

   subroutine run_mom_order(cols, nocc, occ, cmo, eps, nswap)
      !! D_old = occ sum_{k in cols} c_k c_k^T from the reflector, S = 1;
      !! returns the reordered cmo/eps (starting eps -1.0, -0.5, 0.2, 0.9).
      integer, intent(in) :: cols(:), nocc
      real(dp), intent(in) :: occ
      real(dp), intent(out) :: cmo(4, 4), eps(4)
      integer, intent(out) :: nswap
      integer, parameter :: n = 4
      type(trc_linalg_t) :: la
      real(dp) :: smat(n, n), d_old(n, n), t(n, n), sc(n, n), p(n)
      integer :: i, k
      call la%init(n)
      call reflector(n, cmo)
      smat = 0.0_dp
      do i = 1, n
         smat(i, i) = 1.0_dp
      end do
      d_old = 0.0_dp
      do k = 1, size(cols)
         d_old = d_old + occ*spread(cmo(:, cols(k)), 2, n)*spread(cmo(:, cols(k)), 1, n)
      end do
      eps = [-1.0_dp, -0.5_dp, 0.2_dp, 0.9_dp]
      t = 0.0_dp; sc = 0.0_dp; p = 0.0_dp
      !$acc enter data copyin(smat, d_old, cmo, eps) create(t, sc, p)
      call scf_driver_mom_weights(la, n, smat, d_old, cmo, t, sc, p)
      call scf_driver_mom_permute(n, nocc, p, cmo, eps, t, nswap)
      !$acc update self(cmo, eps)
      !$acc exit data delete(smat, d_old, cmo, eps, t, sc, p)
      call la%release()
   end subroutine run_mom_order

   subroutine test_mom_order_swaps(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: cmo(4, 4), c0(4, 4), eps(4)
      integer :: nswap
      call reflector(4, c0)
      ! Previously occupied: columns 1 and 3 -> order 1, 3, 2, 4.
      call run_mom_order([1, 3], 2, 2.0_dp, cmo, eps, nswap)
      call check(error, nswap == 1, "nswap for occupied {1,3} should be 1, got "//to_string(nswap))
      if (allocated(error)) return
      call check(error, maxval(abs(cmo - c0(:, [1, 3, 2, 4]))) < 1.0e-12_dp, "columns not ordered 1,3,2,4")
      if (allocated(error)) return
      call check(error, maxval(abs(eps - [-1.0_dp, 0.2_dp, -0.5_dp, 0.9_dp])) < 1.0e-14_dp, &
                 "eps not permuted with the columns")
      if (allocated(error)) return
      ! Occupied {1,4}: order 1, 4, 2, 3 (energy order kept inside the empty block).
      call run_mom_order([1, 4], 2, 2.0_dp, cmo, eps, nswap)
      call check(error, nswap == 1, "nswap for occupied {1,4} should be 1, got "//to_string(nswap))
      if (allocated(error)) return
      call check(error, maxval(abs(cmo - c0(:, [1, 4, 2, 3]))) < 1.0e-12_dp, "columns not ordered 1,4,2,3")
      if (allocated(error)) return
      call check(error, maxval(abs(eps - [-1.0_dp, 0.9_dp, -0.5_dp, 0.2_dp])) < 1.0e-14_dp, &
                 "eps not permuted with the columns (1,4)")
      if (allocated(error)) return
      ! Occupied {3,4}: both swapped.
      call run_mom_order([3, 4], 2, 2.0_dp, cmo, eps, nswap)
      call check(error, nswap == 2, "nswap for occupied {3,4} should be 2, got "//to_string(nswap))
      if (allocated(error)) return
      call check(error, maxval(abs(cmo - c0(:, [3, 4, 1, 2]))) < 1.0e-12_dp, "columns not ordered 3,4,1,2")
   end subroutine test_mom_order_swaps

   subroutine test_mom_order_noops(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: cmo(4, 4), c0(4, 4), eps(4)
      real(dp), parameter :: eps0(4) = [-1.0_dp, -0.5_dp, 0.2_dp, 0.9_dp]
      integer :: nswap
      call reflector(4, c0)
      call run_mom_order([1, 2], 2, 2.0_dp, cmo, eps, nswap)
      call check(error, nswap == 0, "nswap for the aufbau set should be 0, got "//to_string(nswap))
      if (allocated(error)) return
      call check(error, maxval(abs(cmo - c0)) == 0.0_dp, "an aufbau occupation must not move a column")
      if (allocated(error)) return
      call check(error, all(eps == eps0), "an aufbau occupation must not move eps")
      if (allocated(error)) return
      ! nocc = 0 and nocc >= nao: nothing to choose, nothing moves even with a non-aufbau D_old.
      call run_mom_order([1, 3], 0, 2.0_dp, cmo, eps, nswap)
      call check(error, nswap == 0 .and. maxval(abs(cmo - c0)) == 0.0_dp .and. all(eps == eps0), &
                 "nocc = 0 must return untouched")
      if (allocated(error)) return
      call run_mom_order([1, 3], 4, 2.0_dp, cmo, eps, nswap)
      call check(error, nswap == 0 .and. maxval(abs(cmo - c0)) == 0.0_dp .and. all(eps == eps0), &
                 "nocc = nao must return untouched")
   end subroutine test_mom_order_noops

   subroutine test_mom_order_occ1(error)
      !! The unrestricted path: one electron per orbital, D_old = sum c c^T.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: cmo(4, 4), c0(4, 4), eps(4)
      integer :: nswap
      call reflector(4, c0)
      call run_mom_order([2, 4], 2, 1.0_dp, cmo, eps, nswap)
      call check(error, nswap == 1, "occ = 1: nswap for occupied {2,4} should be 1, got "//to_string(nswap))
      if (allocated(error)) return
      call check(error, maxval(abs(cmo - c0(:, [2, 4, 1, 3]))) < 1.0e-12_dp, "occ = 1: columns not ordered 2,4,1,3")
      if (allocated(error)) return
      call check(error, maxval(abs(eps - [-0.5_dp, 0.9_dp, -1.0_dp, 0.2_dp])) < 1.0e-14_dp, "occ = 1: eps not permuted")
   end subroutine test_mom_order_occ1

   ! --- Fortran-level SCF helpers ------------------------------------------

   subroutine make_basis(xyz, basis, nelec, bas, error)
      character(len=*), intent(in) :: xyz, basis
      integer, intent(out) :: nelec
      type(trc_basis_t), intent(out) :: bas
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm
      call read_xyz(xyz, natm, z, at_r)
      nelec = sum(z)
      call trc_basis_from_json('basis_sets/'//basis//'.json', natm, z, at_r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call bas%to_device()
   end subroutine make_basis

   subroutine test_water_ground_state(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o0, o1
      type(trc_scf_result_t) :: plain, aided
      integer :: nelec, nocc
      call make_basis('water.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      nocc = nelec/2
      o0%level_shift = MU
      o0%conv_energy = 1.0e-11_dp
      o0%conv_diis = 1.0e-7_dp
      o0%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, o0, plain)
      ! level_shift_until far below conv_diis: the shift is on through DIIS
      ! up to the last iteration that can converge, which is never shifted.
      o1 = o0
      o1%level_shift_until = 1.0e-12_dp
      o1%mom = .true.
      call trc_scf_run(bas, nocc, nocc, o1, aided)
      call bas%release()
      call check(error, plain%converged, "plain run: "//trim(plain%message))
      if (allocated(error)) return
      call check(error, abs(plain%energy - E_WATER_RHF) < TOL_REF, "plain energy differs from pyscf")
      if (allocated(error)) return
      call check(error, aided%converged, "mom + level_shift_until: "//trim(aided%message))
      if (allocated(error)) return
      call check(error, abs(aided%energy - plain%energy) < 1.0e-9_dp, "energy moved by the aids")
      if (allocated(error)) return
      call check(error, abs(aided%energy - E_WATER_RHF) < TOL_REF, "aided energy differs from pyscf")
      if (allocated(error)) return
      call check(error, abs(aided%eps(nocc, 1) - plain%eps(nocc, 1)) < 1.0e-7_dp, "HOMO eps carries the shift")
      if (allocated(error)) return
      call check(error, abs(aided%eps(nocc + 1, 1) - plain%eps(nocc + 1, 1)) < 1.0e-7_dp, &
                 "LUMO eps carries the shift (converged iteration was shifted)")
   end subroutine test_water_ground_state

   subroutine test_oh_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o0, o1
      type(trc_scf_result_t) :: plain, aided
      integer :: nelec, na, nb
      call make_basis('oh.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      na = (nelec + 1)/2
      nb = nelec - na
      o0%level_shift = MU
      o0%conv_energy = 1.0e-11_dp
      o0%conv_diis = 1.0e-7_dp
      o0%max_iteration = 100
      call trc_scf_run(bas, na, nb, o0, plain)
      o1 = o0
      o1%level_shift_until = 1.0e-2_dp
      o1%mom = .true.
      call trc_scf_run(bas, na, nb, o1, aided)
      call bas%release()
      call check(error, plain%converged, "plain UHF: "//trim(plain%message))
      if (allocated(error)) return
      call check(error, abs(plain%energy - E_OH_UHF) < TOL_REF, "plain UHF energy differs from pyscf")
      if (allocated(error)) return
      call check(error, aided%converged, "UHF with level_shift_until + mom: "//trim(aided%message))
      if (allocated(error)) return
      call check(error, abs(aided%energy - plain%energy) < 1.0e-9_dp, "UHF energy moved by the aids")
      if (allocated(error)) return
      call check(error, abs(aided%energy - E_OH_UHF) < TOL_REF, "aided UHF energy differs from pyscf")
   end subroutine test_oh_uhf

   subroutine one_iteration(bas, na, nb, base, d0, until, res)
      !! One SCF iteration from the converged density `d0`, DIIS on from the
      !! start (no damped phase), with `level_shift` = MU and the given
      !! `level_shift_until`. `conv_diis` is 1e-30: the commutator at the
      !! converged density (8.5e-11 measured) is above it by construction, so
      !! the iteration is one that may shift (an iteration under `conv_diis`
      !! never does), whatever the reference converged to.
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: na, nb
      type(trc_scf_options_t), intent(in) :: base
      real(dp), intent(in) :: d0(:, :, :), until
      type(trc_scf_result_t), intent(out) :: res
      type(trc_scf_options_t) :: o
      o = base
      o%max_iteration = 1
      o%conv_diis = 1.0e-30_dp
      o%level_shift_until = until
      call trc_scf_run(bas, na, nb, o, res, dguess=d0)
   end subroutine one_iteration

   subroutine test_shift_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o
      type(trc_scf_result_t) :: ref, off, on
      integer :: nelec, nocc
      call make_basis('water.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      nocc = nelec/2
      o%level_shift = MU
      o%conv_energy = 1.0e-12_dp
      o%conv_diis = 1.0e-10_dp
      o%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, o, ref)
      call one_iteration(bas, nocc, nocc, o, ref%dmat, 0.0_dp, off)
      call one_iteration(bas, nocc, nocc, o, ref%dmat, 1.0e-30_dp, on)
      call bas%release()
      call check(error, ref%converged, "reference run: "//trim(ref%message))
      if (allocated(error)) return
      call check(error, abs(off%eps(nocc + 1, 1) - ref%eps(nocc + 1, 1)) < 1.0e-6_dp, &
                 "without level_shift_until the LUMO must not move")
      if (allocated(error)) return
      call check(error, abs(on%eps(nocc, 1) - ref%eps(nocc, 1)) < 1.0e-6_dp, "the shift moved the HOMO")
      if (allocated(error)) return
      call check(error, abs(on%eps(nocc + 1, 1) - ref%eps(nocc + 1, 1) - MU) < 1.0e-5_dp, &
                 "LUMO with level_shift_until is not mu above the unshifted one: "// &
                 to_string(on%eps(nocc + 1, 1) - ref%eps(nocc + 1, 1), "f12.8"))
      if (allocated(error)) return
      call check(error, index(on%message, "level shifted") > 0, &
                 "the non-converged message should say the virtuals are shifted: "//trim(on%message))
      if (allocated(error)) return
      call check(error, index(off%message, "level shifted") == 0, "an unshifted run must not claim a shift")
   end subroutine test_shift_rhf

   subroutine test_shift_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o
      type(trc_scf_result_t) :: ref, off, on
      integer :: nelec, na, nb
      call make_basis('oh.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      na = (nelec + 1)/2
      nb = nelec - na
      o%level_shift = MU
      o%ndamp = 0
      o%diis_start = 1.0e10_dp
      o%conv_energy = 1.0e-12_dp
      o%conv_diis = 1.0e-10_dp
      o%max_iteration = 100
      call trc_scf_run(bas, na, nb, o, ref)
      call one_iteration(bas, na, nb, o, ref%dmat, 0.0_dp, off)
      call one_iteration(bas, na, nb, o, ref%dmat, 1.0e-30_dp, on)
      call bas%release()
      call check(error, ref%converged, "reference UHF: "//trim(ref%message))
      if (allocated(error)) return
      call check(error, abs(off%eps(na + 1, 1) - ref%eps(na + 1, 1)) < 1.0e-6_dp .and. &
                 abs(off%eps(nb + 1, 2) - ref%eps(nb + 1, 2)) < 1.0e-6_dp, &
                 "without level_shift_until neither LUMO may move")
      if (allocated(error)) return
      call check(error, abs(on%eps(na, 1) - ref%eps(na, 1)) < 1.0e-6_dp .and. &
                 abs(on%eps(nb, 2) - ref%eps(nb, 2)) < 1.0e-6_dp, "the shift moved an occupied level")
      if (allocated(error)) return
      call check(error, abs(on%eps(na + 1, 1) - ref%eps(na + 1, 1) - MU) < 1.0e-5_dp, &
                 "alpha LUMO is not mu above the unshifted one: "// &
                 to_string(on%eps(na + 1, 1) - ref%eps(na + 1, 1), "f12.8"))
      if (allocated(error)) return
      call check(error, abs(on%eps(nb + 1, 2) - ref%eps(nb + 1, 2) - MU) < 1.0e-5_dp, &
                 "beta LUMO is not mu above the unshifted one: "// &
                 to_string(on%eps(nb + 1, 2) - ref%eps(nb + 1, 2), "f12.8"))
   end subroutine test_shift_uhf

   subroutine test_converging_iteration_unshifted(error)
      !! OH UHF with a level shift, restarted from its own converged density,
      !! converges at iteration 2 still inside the damped start (ndamp = 4)
      !! with F shifted: before the fix the returned alpha LUMO sat MU
      !! (+0.49999999802) too high. Now an iteration that can converge is
      !! never shifted.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o, o_shifted
      type(trc_scf_result_t) :: ref, restart
      integer :: nelec, na, nb
      call make_basis('oh.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      na = (nelec + 1)/2
      nb = nelec - na
      o%conv_energy = 1.0e-10_dp
      o%conv_diis = 1.0e-7_dp
      o%max_iteration = 100
      call trc_scf_run(bas, na, nb, o, ref)
      o_shifted = o
      o_shifted%level_shift = MU
      call trc_scf_run(bas, na, nb, o_shifted, restart, dguess=ref%dmat)
      call bas%release()
      call check(error, ref%converged, "unshifted reference: "//trim(ref%message))
      if (allocated(error)) return
      call check(error, restart%converged .and. restart%iterations <= UHF_NDAMP, &
                 "restart converges inside the damped start: "//to_string(restart%iterations)//" iterations")
      if (allocated(error)) return
      call check(error, abs(restart%energy - ref%energy) < 1.0e-8_dp, "restart energy moved")
      if (allocated(error)) return
      call check(error, abs(restart%eps(na + 1, 1) - ref%eps(na + 1, 1)) < 1.0e-6_dp, &
                 "alpha LUMO carries the shift: "//to_string(restart%eps(na + 1, 1) - ref%eps(na + 1, 1), "f12.8"))
      if (allocated(error)) return
      call check(error, abs(restart%eps(nb + 1, 2) - ref%eps(nb + 1, 2)) < 1.0e-6_dp, "beta LUMO carries the shift")
   end subroutine test_converging_iteration_unshifted

   subroutine test_soscf_clears_shift(error)
      !! SOSCF takes over at the end of iteration 1, which was level shifted
      !! (ndamp = 2). An unreachable conv_diis keeps the run unconverged, and
      !! the message must not blame a shift SOSCF's own semicanonical eps do
      !! not carry.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o
      type(trc_scf_result_t) :: res
      integer :: nelec, nocc
      call make_basis('water.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      nocc = nelec/2
      o%level_shift = MU
      o%ndamp = 2
      o%soscf = .true.
      o%soscf_start = 1.0e6_dp
      o%conv_energy = 1.0e-30_dp
      o%conv_diis = 1.0e-30_dp
      o%max_iteration = 5
      call trc_scf_run(bas, nocc, nocc, o, res)
      call bas%release()
      call check(error, .not. res%converged .and. res%iterations >= 1, "the run must stay unconverged")
      if (allocated(error)) return
      call check(error, index(res%message, "level shifted") == 0, &
                 "SOSCF's eps blamed on the shift: "//trim(res%message))
   end subroutine test_soscf_clears_shift

   subroutine test_mom_poor_guess(error)
      !! Restricted runs have DIIS on from iteration 1, so MOM used to lock the
      !! occupation of the core guess: water/6-31G* converged 7.93 Eh above
      !! the aufbau run (HCN 3.31 Eh). MOM now waits for a settled commutator.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o
      type(trc_scf_result_t) :: aufbau, followed
      integer :: nelec, nocc
      call make_basis('water.xyz', '6-31gs', nelec, bas, error)
      if (allocated(error)) return
      nocc = nelec/2
      o%guess = "core"
      o%conv_energy = 1.0e-10_dp
      o%conv_diis = 1.0e-8_dp
      o%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, o, aufbau)
      o%mom = .true.
      call trc_scf_run(bas, nocc, nocc, o, followed)
      call bas%release()
      call check(error, aufbau%converged, "aufbau run: "//trim(aufbau%message))
      if (allocated(error)) return
      call check(error, followed%converged, "mom run: "//trim(followed%message))
      if (allocated(error)) return
      call check(error, abs(followed%energy - aufbau%energy) < 1.0e-8_dp, &
                 "mom from a core guess locked a wrong occupation: E - E0 = "// &
                 to_string(followed%energy - aufbau%energy, "es10.2"))
   end subroutine test_mom_poor_guess

   subroutine test_mom_holds_excited_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: o
      type(trc_scf_result_t) :: ground, held, fallen
      real(dp), allocatable :: dex(:, :, :)
      integer, allocatable :: occ_cols(:)
      integer :: nelec, nocc, k
      call make_basis('water.xyz', 'cc-pvdz', nelec, bas, error)
      if (allocated(error)) return
      nocc = nelec/2
      o%conv_energy = 1.0e-10_dp
      o%conv_diis = 1.0e-8_dp
      o%max_iteration = 100
      call trc_scf_run(bas, nocc, nocc, o, ground)
      call check(error, ground%converged, "ground state: "//trim(ground%message))
      if (allocated(error)) then
         call bas%release()
         return
      end if
      ! HOMO and LUMO swapped: columns 1..nocc-1 and nocc+1.
      occ_cols = [(k, k=1, nocc - 1), nocc + 1]
      allocate (dex(size(ground%dmat, 1), size(ground%dmat, 2), 1))
      dex(:, :, 1) = 2.0_dp*matmul(ground%cmo(:, occ_cols, 1), transpose(ground%cmo(:, occ_cols, 1)))
      o%mom = .true.
      call trc_scf_run(bas, nocc, nocc, o, held, dguess=dex)
      o%mom = .false.
      call trc_scf_run(bas, nocc, nocc, o, fallen, dguess=dex)
      call bas%release()
      call check(error, held%converged, "mom run from the swapped occupation: "//trim(held%message))
      if (allocated(error)) return
      call check(error, held%energy > ground%energy + MIN_EXCITATION, &
                 "mom did not hold the excited occupation: E - E0 = "//to_string(held%energy - ground%energy, "f12.8"))
      if (allocated(error)) return
      call check(error, fallen%converged, "aufbau run from the swapped occupation: "//trim(fallen%message))
      if (allocated(error)) return
      call check(error, abs(fallen%energy - ground%energy) < 1.0e-8_dp, &
                 "without mom the SCF should fall back to the ground state: E - E0 = "// &
                 to_string(fallen%energy - ground%energy, "es10.2"))
   end subroutine test_mom_holds_excited_rhf

   ! --- C ABI --------------------------------------------------------------

   subroutine capi_free(hs)
      type(c_ptr), intent(inout) :: hs(:)
      integer :: i
      integer(c_int) :: rc
      do i = 1, size(hs)
         if (c_associated(hs(i))) then
            rc = trc_destroy(hs(i))
            hs(i) = c_null_ptr
         end if
      end do
   end subroutine capi_free

   subroutine capi_water(h, maxit, conv_diis, error)
      !! A water/6-31G* handle, every setup call checked. On failure `h` is
      !! destroyed here and the error set.
      type(c_ptr), intent(out) :: h
      integer, intent(in) :: maxit
      real(dp), intent(in) :: conv_diis
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: z(:)
      real(dp), allocatable :: at_r(:, :)
      integer :: natm
      integer(c_int) :: rc
      h = c_null_ptr
      call read_xyz('water.xyz', natm, z, at_r)
      rc = trc_create(h)
      call check(error, rc == TRC_OK, "trc_create")
      if (allocated(error)) return
      rc = trc_set_molecule(h, int(natm, c_int), real(z, c_double), reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
      call check(error, rc == TRC_OK, "trc_set_molecule")
      if (allocated(error)) go to 90
      rc = trc_set_basis_json(h, 'basis_sets/6-31gs.json'//c_null_char)
      call check(error, rc == TRC_OK, "trc_set_basis_json")
      if (allocated(error)) go to 90
      rc = trc_set_method(h, ''//c_null_char, 3_c_int)
      call check(error, rc == TRC_OK, "trc_set_method")
      if (allocated(error)) go to 90
      rc = trc_set_convergence(h, 1.0e-11_c_double, real(conv_diis, c_double), int(maxit, c_int), 10_c_int)
      call check(error, rc == TRC_OK, "trc_set_convergence")
      if (allocated(error)) go to 90
      rc = trc_set_level_shift(h, real(MU, c_double))
      call check(error, rc == TRC_OK, "trc_set_level_shift")
      if (allocated(error)) go to 90
      return
90    rc = trc_destroy(h)
      h = c_null_ptr
   end subroutine capi_water

   subroutine capi_ground(hs, nao, dmat, cmo, eps, e0, error)
      !! A converged ground-state water handle in hs(1), and its density, orbitals and energy.
      type(c_ptr), intent(inout) :: hs(:)
      integer, intent(out) :: nao
      real(dp), allocatable, intent(out) :: dmat(:, :), cmo(:, :), eps(:)
      real(dp), intent(out) :: e0
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc, n
      real(c_double) :: e
      call capi_water(hs(1), 100, 1.0e-10_dp, error)
      if (allocated(error)) return
      rc = trc_run_scf(hs(1))
      call check(error, rc == TRC_OK, "ground state trc_run_scf")
      if (allocated(error)) return
      rc = trc_nao(hs(1), n)
      call check(error, rc == TRC_OK, "trc_nao")
      if (allocated(error)) return
      nao = n
      allocate (dmat(nao, nao), cmo(nao, nao), eps(nao))
      rc = trc_density(hs(1), dmat)
      call check(error, rc == TRC_OK, "trc_density")
      if (allocated(error)) return
      rc = trc_mo_coefficients(hs(1), cmo)
      call check(error, rc == TRC_OK, "trc_mo_coefficients")
      if (allocated(error)) return
      rc = trc_mo_energies(hs(1), eps)
      call check(error, rc == TRC_OK, "trc_mo_energies")
      if (allocated(error)) return
      rc = trc_energy(hs(1), e)
      call check(error, rc == TRC_OK, "trc_energy")
      e0 = e
   end subroutine capi_ground

   subroutine test_capi_validation(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hs(1)
      integer(c_int) :: rc
      real(c_double) :: nan
      hs = c_null_ptr
      nan = 0.0_c_double
      nan = nan/nan
      call capi_water(hs(1), 100, 1.0e-8_dp, error)
      if (allocated(error)) return
      rc = trc_set_scf_stabilise(hs(1), 2_c_int, 0.0_c_double)
      call check(error, rc == TRC_ERR_BADARG, "mom = 2 must be refused")
      if (.not. allocated(error)) then
         rc = trc_set_scf_stabilise(hs(1), 0_c_int, -1.0_c_double)
         call check(error, rc == TRC_ERR_BADARG, "a negative level_shift_until must be refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_stabilise(hs(1), 0_c_int, nan)
         call check(error, rc == TRC_ERR_BADARG, "a NaN level_shift_until must be refused")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_stabilise(hs(1), 0_c_int, 0.0_c_double)
         call check(error, rc == TRC_OK, "(0, 0) must be accepted")
      end if
      if (.not. allocated(error)) then
         rc = trc_set_scf_stabilise(hs(1), 1_c_int, 1.0e-2_c_double)
         call check(error, rc == TRC_OK, "(1, 1e-2) must be accepted")
      end if
      call capi_free(hs)
   end subroutine test_capi_validation

   subroutine capi_probe(hs, i, dmat, mom, until, e_lumo_homo, error)
      !! A fresh handle hs(i): one iteration from the converged density with the
      !! given stabilise settings; returns eps(nocc), eps(nocc+1).
      type(c_ptr), intent(inout) :: hs(:)
      integer, intent(in) :: i
      real(dp), intent(in), target :: dmat(:, :)
      integer, intent(in) :: mom
      real(dp), intent(in) :: until
      real(dp), intent(out) :: e_lumo_homo(2)
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc
      real(dp), allocatable :: eps(:)
      call capi_water(hs(i), 1, 1.0e-30_dp, error)
      if (allocated(error)) return
      rc = trc_set_scf_stabilise(hs(i), int(mom, c_int), real(until, c_double))
      call check(error, rc == TRC_OK, "trc_set_scf_stabilise")
      if (allocated(error)) return
      rc = trc_set_guess(hs(i), TRC_GUESS_GIVEN, c_loc(dmat(1, 1)), 1_c_int)
      call check(error, rc == TRC_OK, "trc_set_guess")
      if (allocated(error)) return
      rc = trc_run_scf(hs(i))
      call check(error, rc == TRC_ERR_NOCONV, "one iteration cannot be converged (status "//to_string(int(rc))//")")
      if (allocated(error)) return
      allocate (eps(size(dmat, 1)))
      rc = trc_mo_energies(hs(i), eps)
      call check(error, rc == TRC_OK, "trc_mo_energies")
      if (allocated(error)) return
      e_lumo_homo = eps(5:6)
   end subroutine capi_probe

   subroutine test_capi_shift(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hs(3)
      real(dp), allocatable :: dmat(:, :), cmo(:, :), eps(:)
      real(dp) :: e0, off(2), on(2)
      integer :: nao
      hs = c_null_ptr
      call capi_ground(hs, nao, dmat, cmo, eps, e0, error)
      if (.not. allocated(error)) call capi_probe(hs, 2, dmat, 0, 0.0_dp, off, error)
      if (.not. allocated(error)) call capi_probe(hs, 3, dmat, 0, 1.0e-30_dp, on, error)
      if (.not. allocated(error)) call check(error, abs(off(2) - eps(6)) < 1.0e-6_dp, &
                                             "level_shift_until = 0: the LUMO must not move")
      if (.not. allocated(error)) call check(error, abs(on(1) - eps(5)) < 1.0e-6_dp, "the shift moved the HOMO")
      if (.not. allocated(error)) call check(error, abs(on(2) - eps(6) - MU) < 1.0e-5_dp, &
                                             "level_shift_until through the C ABI: LUMO not mu above: "// &
                                             to_string(on(2) - eps(6), "f12.8"))
      call capi_free(hs)
   end subroutine test_capi_shift

   subroutine capi_swapped_run(hs, i, dex, mom, e, conv, error)
      type(c_ptr), intent(inout) :: hs(:)
      integer, intent(in) :: i, mom
      real(dp), intent(in), target :: dex(:, :)
      real(dp), intent(out) :: e
      integer, intent(out) :: conv
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: rc, flag
      real(c_double) :: ce
      call capi_water(hs(i), 100, 1.0e-8_dp, error)
      if (allocated(error)) return
      rc = trc_set_scf_stabilise(hs(i), int(mom, c_int), 0.0_c_double)
      call check(error, rc == TRC_OK, "trc_set_scf_stabilise")
      if (allocated(error)) return
      rc = trc_set_guess(hs(i), TRC_GUESS_GIVEN, c_loc(dex(1, 1)), 1_c_int)
      call check(error, rc == TRC_OK, "trc_set_guess")
      if (allocated(error)) return
      rc = trc_run_scf(hs(i))
      call check(error, rc == TRC_OK, "trc_run_scf from the swapped occupation (status "//to_string(int(rc))//")")
      if (allocated(error)) return
      rc = trc_converged(hs(i), flag)
      conv = flag
      rc = trc_energy(hs(i), ce)
      call check(error, rc == TRC_OK, "trc_energy")
      e = ce
   end subroutine capi_swapped_run

   subroutine test_capi_mom(error)
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: hs(3)
      real(dp), allocatable :: dmat(:, :), cmo(:, :), eps(:), dex(:, :)
      real(dp) :: e0, e_mom, e_aufbau
      integer :: nao, conv, k
      integer, allocatable :: cols(:)
      hs = c_null_ptr
      call capi_ground(hs, nao, dmat, cmo, eps, e0, error)
      if (allocated(error)) then
         call capi_free(hs)
         return
      end if
      cols = [(k, k=1, 4), 6]
      allocate (dex(nao, nao))
      dex = 2.0_dp*matmul(cmo(:, cols), transpose(cmo(:, cols)))
      call capi_swapped_run(hs, 2, dex, 1, e_mom, conv, error)
      if (.not. allocated(error)) call check(error, conv == 1, "mom run did not converge")
      if (.not. allocated(error)) call check(error, e_mom > e0 + MIN_EXCITATION, &
                                             "mom through the C ABI did not hold the excited occupation: E - E0 = "// &
                                             to_string(e_mom - e0, "f12.8"))
      if (.not. allocated(error)) call capi_swapped_run(hs, 3, dex, 0, e_aufbau, conv, error)
      if (.not. allocated(error)) call check(error, abs(e_aufbau - e0) < 1.0e-8_dp, &
                                             "without mom the swapped start should fall back to the ground state: E - E0 = "// &
                                             to_string(e_aufbau - e0, "es10.2"))
      call capi_free(hs)
   end subroutine test_capi_mom

end module check_scf_stabilise_suite

program check_scf_stabilise
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_api, only: trc_bind_device
   use check_scf_stabilise_suite, only: collect_scf_stabilise
   implicit none
   integer :: stat, dev
   dev = trc_bind_device(0)
   stat = 0
   ! Serial: the C-ABI case builds a context, which is not thread-safe, and
   ! an -fopenmp build runs test-drive cases at once.
   call run_testsuite(collect_scf_stabilise, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_scf_stabilise
