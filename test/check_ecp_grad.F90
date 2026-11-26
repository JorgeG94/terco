module check_ecp_grad_suite
   !!
   !! The analytic effective-core-potential term of the nuclear gradient.
   !!
   !! Three layers, each against something independent of the code under test:
   !!
   !!   1. The derivative integrals. For HI, a bent TeH2 and an HI dimer (two
   !!      ECP atoms) in def2-SVP (Cartesian), dU/dR_X, assembled from
   !!      `trc_ecp_deriv_cart`'s one-ECP-atom bra derivatives (basis-centre
   !!      terms on the bra and, by symmetry, the ket, plus the ECP centre by
   !!      translational invariance), against a five-point difference of
   !!      `trc_ecp_matrix_cart` over every atom X and component: the basis
   !!      centres AND the ECP centre move. `trc_ecp_gradient` on a fixed
   !!      synthetic density against the same differences contracted with it:
   !!      with two ECP atoms that fails if a centre term lands on the wrong one.
   !!      On a pull request the dimer's contraction is one directional
   !!      derivative, every atom moving at once (four matrices); element by
   !!      element it is LONG. HI with two more iodines, 9 and 60 bohr off
   !!      (three ECP atoms, far enough for the screen to drop tasks): screened
   !!      against unscreened, one chunk of ECP atoms against three, both
   !!      against the density contracted with `trc_ecp_deriv_cart` (the
   !!      dense, one-atom-a-time reference, a different derivative basis);
   !!      the chunks the free device memory forces, and the refusal when it
   !!      cannot hold one atom with the kernels' local memory (both through
   !!      the memory query's test cap, so on the host too).
   !!   2. Whole gradients against pyscf 2.14 (the same BSE-JSON basis and ECP,
   !!      tight SCF), read from pyscf_refs/ecp_grad.ref: RHF of HI and TeH2
   !!      (spherical and Cartesian) and of the HI dimer, UHF of the HI cation,
   !!      RI-MP2 of HI (pyscf's side a difference of its DF-MP2 energy).
   !!   3. The whole gradient through the C API, at the same two molecules in
   !!      def2-SVP's own spherical functions:
   !!      `trc_run_gradient` for RHF, UHF, RKS/PBE and RI-MP2 against a central
   !!      difference of the reconverged total energy, a handful of
   !!      coordinates each. Without the ECP term every one of these
   !!      gradients is off by 0.02 to 0.06 Eh/bohr, against tolerances of
   !!      1e-6 or so.
   !!
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_basis_types, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_json, only: trc_ecp_from_json
   use trc_ecp, only: trc_ecp_matrix_cart, trc_ecp_deriv_cart, trc_ecp_gradient, trc_ecp_screen_stats_t
   use trc_device_memory, only: device_memory_set_cap
   use trc_error, only: error_t
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_convergence, trc_set_method, trc_set_rimp2, trc_set_force_cartesian, &
                               trc_run_scf, trc_run_rimp2, trc_run_gradient, trc_gradient, trc_energy, trc_message, &
                               TRC_OK, TRC_ERR_UNSUPPORTED
   use trc_capi, only: trc_check_gradient_upfront
   implicit none
   private
   public :: collect_ecp_grad_tests

   interface
      function c_setenv(name, value, overwrite) bind(c, name="setenv") result(r)
         import :: c_int, c_char
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
      end function c_setenv
      function c_unsetenv(name) bind(c, name="unsetenv") result(r)
         import :: c_int, c_char
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int) :: r
      end function c_unsetenv
   end interface

   character(len=*), parameter :: DEF2_SVP = 'basis_sets/def2-svp.json'
   real(dp), parameter :: INT_TOL = 1.0e-10_dp
   !! The derivative matrix against its five-point difference, absolute: measured 8e-13.
   real(dp), parameter :: CONTRACT_TOL = 1.0e-10_dp
   !! trc_ecp_gradient against the contracted difference of U, Eh/bohr: measured 3e-12.
   real(dp), parameter :: DIRECTION_TOL = 1.0e-8_dp
   !! The dimer's directional derivative, a sum over all twelve components: measured 7e-11 (gfortran), 3e-11 (nvfortran), on a value of 1.2.
   real(dp), parameter :: SCREEN_TOL = 1.0e-12_dp
   !! Screened against unscreened, Eh/bohr: the dropped blocks are below 1e-20.
   real(dp), parameter :: CHUNK_TOL = 1.0e-13_dp
   !! Three chunks against one: the same tasks in the same order, each with its own slot.
   real(dp), parameter :: PYSCF_TOL = 1.0e-7_dp
   !! Whole gradients against pyscf_refs/ecp_grad.ref: measured at most 1e-9 (RHF, UHF and RI-MP2).
   real(dp), parameter :: GRAD_TOL_HF = 2.0e-6_dp
   real(dp), parameter :: GRAD_TOL_KS = 2.0e-6_dp
   !! The DFT difference carries the quadrature's own noise (the grid moves
   !! with the atoms and is finite).
   real(dp), parameter :: GRAD_TOL_MP2 = 5.0e-6_dp

   type :: pyscf_refs_t
      !! pyscf_refs/ecp_grad.ref.
      real(dp) :: hi_sph(3, 2), hi_cart(3, 2), teh2_sph(3, 3), teh2_cart(3, 3), hi_bond(3, 2)
      real(dp) :: dimer(3, 4), uhf(3, 2), mp2(3, 2), e_uhf, e_mp2, e_skip
   end type pyscf_refs_t

contains

   subroutine collect_ecp_grad_tests(testsuite)
      !! No argument: what a pull request runs (about half a minute). `long`: the
      !! heavier integral and pyscf cases. `fd`: the gradients against energy
      !! differences, which reconverge about fifty SCFs (minutes). The last two
      !! are registered LONG.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      character(len=8) :: arg

      arg = ""
      if (command_argument_count() >= 1) call get_command_argument(1, arg)
      if (trim(arg) == "fd") then
         testsuite = [ &
                     new_unittest("RHF gradient of HI vs the energy difference", test_rhf_hi), &
                     new_unittest("RHF gradient of bent TeH2 vs the energy difference", test_rhf_teh2), &
                     new_unittest("RKS/PBE gradient of HI vs the energy difference", test_pbe_hi), &
                     new_unittest("UHF gradient of the HI cation (doublet) vs the energy difference", test_uhf_hi), &
                     new_unittest("RI-MP2 gradient of HI vs the energy difference", test_mp2_hi), &
                     new_unittest("RI-MP2 gradient of bent TeH2 vs the energy difference", test_mp2_teh2) &
                     ]
      else if (trim(arg) == "long") then
         testsuite = [ &
                     new_unittest("dU/dR of bent TeH2 vs differences of the ECP matrix", test_ints_teh2), &
                     new_unittest("dU/dR of an HI dimer (two ECP atoms) vs differences of the ECP matrix", test_ints_dimer), &
                     new_unittest("RHF gradients of HI, Cartesian TeH2 and an HI dimer vs pyscf", test_rhf_pyscf_long) &
                     ]
      else
         testsuite = [ &
                     new_unittest("dU/dR of HI vs differences of the ECP matrix", test_ints_hi), &
                     new_unittest("ECP gradient of an HI dimer (two ECP atoms) along a direction", test_dimer_direction), &
                     new_unittest("ECP gradient with three ECP atoms: screened, chunked and the dense reference", &
                                  test_screen_chunk), &
                     new_unittest("RHF gradients of bent TeH2 and Cartesian HI vs pyscf", test_rhf_pyscf), &
                     new_unittest("UHF gradient of the HI cation (doublet) vs pyscf", test_uhf_pyscf), &
                     new_unittest("RI-MP2 gradient of HI vs pyscf", test_mp2_pyscf), &
                     new_unittest("HI/def2-SVP with an ECP: the gradient pre-check accepts it", test_upfront_ok), &
                     new_unittest("an f basis function with an ECP: the gradient is refused up front", test_f_refused) &
                     ]
      end if
   end subroutine collect_ecp_grad_tests

   subroutine geometry(name, nat, elements, r)
      character(len=*), intent(in) :: name
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: elements(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      select case (name)
      case ("hi")
         ! Off-axis, so every component of both atoms is live.
         nat = 2
         allocate (elements(nat), r(3, nat))
         elements = [1, 53]
         r(:, 1) = [0.1_dp, -0.2_dp, 0.0_dp]
         r(:, 2) = [0.3_dp, 0.2_dp, 3.04_dp]
      case ("hi_dimer")
         ! Two ECP atoms, neither on an axis of the other.
         nat = 4
         allocate (elements(nat), r(3, nat))
         elements = [1, 53, 1, 53]
         r(:, 1) = [0.1_dp, -0.2_dp, 0.0_dp]
         r(:, 2) = [0.3_dp, 0.2_dp, 3.04_dp]
         r(:, 3) = [4.6_dp, 1.1_dp, 0.9_dp]
         r(:, 4) = [4.9_dp, 4.0_dp, 1.9_dp]
      case default
         ! Te at the apex, 2.2 / 1.8 bohr to each H: a bent, non-symmetric frame.
         nat = 3
         allocate (elements(nat), r(3, nat))
         elements = [52, 1, 1]
         r(:, 1) = [0.05_dp, 0.0_dp, 0.3_dp]
         r(:, 2) = [2.2_dp, 0.1_dp, -1.8_dp]
         r(:, 3) = [-2.1_dp, -0.1_dp, -1.7_dp]
      end select
   end subroutine geometry

   subroutine ecp_u(nat, elements, r, ecp, u, err)
      !! The Cartesian ECP matrix at `r`; `ecp` is read once by the caller.
      integer, intent(in) :: nat, elements(:)
      real(dp), intent(in) :: r(3, nat)
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), allocatable, intent(out) :: u(:, :)
      type(error_t), intent(inout) :: err
      type(trc_basis_t) :: bas

      call trc_basis_from_json(DEF2_SVP, nat, elements, r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      if (err%has_error()) return
      call trc_ecp_matrix_cart(bas, ecp, u, err)
      call bas%release()
   end subroutine ecp_u

   subroutine synthetic_density(dens)
      !! Symmetric and nowhere special: a density the ECP gradient can be
      !! contracted with, no SCF needed.
      real(dp), intent(out) :: dens(:, :)
      integer :: j, k

      do j = 1, size(dens, 2)
         do k = 1, size(dens, 1)
            dens(k, j) = cos(0.37_dp*real(k + j, dp))/(1.0_dp + 0.05_dp*real(abs(k - j), dp))
         end do
      end do
   end subroutine synthetic_density

   subroutine test_dimer_direction(error)
      !! trc_ecp_gradient of the HI dimer against a five-point difference of
      !! sum(density * U) with every atom moving at once: four matrices, and
      !! all twelve components enter. A centre term put on the wrong ECP atom
      !! shifts the sum by that term times the two atoms' relative motion.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      integer, parameter :: STEP(4) = [1, -1, 2, -2]
      real(dp), parameter :: H = 1.0e-3_dp
      real(dp), parameter :: WEIGHT(4) = [8.0_dp, -8.0_dp, -1.0_dp, 1.0_dp]
      integer :: nat, ia, c, k
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :), v(:, :), g(:, :), dens(:, :), u(:, :)
      real(dp) :: fd, ana, diff

      call geometry("hi_dimer", nat, elements, r)
      call trc_basis_from_json(DEF2_SVP, nat, elements, r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      call check(error,.not. err%has_error(), "dimer basis must build: "//err%get_message())
      if (allocated(error)) return
      call trc_ecp_from_json(DEF2_SVP, nat, elements, ecp, err)
      call check(error,.not. err%has_error(), "dimer ECP must read: "//err%get_message())
      if (allocated(error)) return

      allocate (v(3, nat), g(3, nat), dens(bas%nao, bas%nao))
      do ia = 1, nat
         do c = 1, 3
            v(c, ia) = cos(1.3_dp*real(c, dp) + 2.1_dp*real(ia, dp))
         end do
      end do
      call synthetic_density(dens)
      g = 0.0_dp
      call trc_ecp_gradient(bas, ecp, dens, g, err)
      call check(error,.not. err%has_error(), "dimer gradient: "//err%get_message())
      if (allocated(error)) return
      ana = sum(g*v)

      fd = 0.0_dp
      do k = 1, 4
         call ecp_u(nat, elements, r + real(STEP(k), dp)*H*v, ecp, u, err)
         call check(error,.not. err%has_error(), "dimer displaced matrix: "//err%get_message())
         if (allocated(error)) return
         fd = fd + WEIGHT(k)*sum(dens*u)
      end do
      fd = fd/(12.0_dp*H)
      diff = abs(ana - fd)
      call global_logger%info("  hi_dimer: directional derivative "//to_string(fd, "es14.6")// &
                              ", |trc_ecp_gradient - difference| = "//to_string(diff, "es10.3"))
      call check(error, abs(fd) > 1.0e-2_dp, "dimer: the directional difference is ~0, the test is blind")
      if (allocated(error)) return
      call check(error, diff < DIRECTION_TOL, "dimer: trc_ecp_gradient along a direction vs the difference, |diff| = "// &
                 to_string(diff, "es12.4"))
      call bas%release()
      call ecp%release()
   end subroutine test_dimer_direction

   subroutine dense_reference(bas, ecp, dens, g, err)
      !! The same contraction from `trc_ecp_deriv_cart`, one ECP atom at a
      !! time: its bra derivative times the density, twice for the ket, on the
      !! bra's atom, and minus that on the ECP atom.
      type(trc_basis_t), intent(in) :: bas
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), intent(in) :: dens(:, :)
      real(dp), intent(out) :: g(:, :)
      type(error_t), intent(inout) :: err
      real(dp), allocatable :: d(:, :, :)
      real(dp) :: w
      integer :: ia, ish, mu, c, x

      g = 0.0_dp
      do ia = 1, bas%natm
         if (.not. any(ecp%ecp_at(1:ecp%nrow) == ia)) cycle
         call trc_ecp_deriv_cart(bas, ecp, ia, d, err)
         if (err%has_error()) return
         do ish = 1, bas%nshell
            x = bas%sh_at(ish)
            do mu = bas%sh_ao(ish), bas%sh_ao(ish) + (bas%sh_l(ish) + 1)*(bas%sh_l(ish) + 2)/2 - 1
               do c = 1, 3
                  w = 2.0_dp*dot_product(dens(mu, :), d(mu, :, c))
                  g(c, x) = g(c, x) + w
                  g(c, ia) = g(c, ia) - w
               end do
            end do
         end do
         deallocate (d)
      end do
   end subroutine dense_reference

   subroutine test_screen_chunk(error)
      !! HI, an iodine 9 bohr off it and another 60 bohr off: three ECP
      !! atoms (and a hydrogen without one), tasks the screen drops,
      !! (TERCO_ECP_GRAD_ATOM_CHUNK=1) one ECP atom a chunk, and a capped
      !! device memory that forces chunks by itself, then refuses.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err, refused
      type(trc_ecp_screen_stats_t) :: st, st_all, st_chunk, st_fit, st_short
      integer, parameter :: NAT = 4
      integer, parameter :: ELEMENTS(NAT) = [1, 53, 53, 53]
      real(dp), parameter :: R(3, NAT) = reshape([0.1_dp, -0.2_dp, 0.0_dp, &
                                                  0.3_dp, 0.2_dp, 3.04_dp, &
                                                  9.0_dp, 1.0_dp, 2.0_dp, &
                                                  60.0_dp, 0.5_dp, 1.0_dp], [3, NAT])
      integer :: rc_set, rc_unset
      integer(int64) :: whole
      real(dp), allocatable :: dens(:, :)
      real(dp) :: g(3, NAT), g_all(3, NAT), g_chunk(3, NAT), g_ref(3, NAT)

      call trc_basis_from_json(DEF2_SVP, NAT, ELEMENTS, R, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      call check(error,.not. err%has_error(), "three I: basis must build: "//err%get_message())
      if (allocated(error)) return
      call trc_ecp_from_json(DEF2_SVP, NAT, ELEMENTS, ecp, err)
      call check(error,.not. err%has_error(), "three I: ECP must read: "//err%get_message())
      if (allocated(error)) return
      allocate (dens(bas%nao, bas%nao))
      call synthetic_density(dens)

      g = 0.0_dp
      call trc_ecp_gradient(bas, ecp, dens, g, err, stats=st)
      g_all = 0.0_dp
      call trc_ecp_gradient(bas, ecp, dens, g_all, err, screen=.false., stats=st_all)
      rc_set = int(c_setenv("TERCO_ECP_GRAD_ATOM_CHUNK"//c_null_char, "1"//c_null_char, 1_c_int))
      g_chunk = 0.0_dp
      call trc_ecp_gradient(bas, ecp, dens, g_chunk, err, stats=st_chunk)
      rc_unset = int(c_unsetenv("TERCO_ECP_GRAD_ATOM_CHUNK"//c_null_char))
      call check(error, rc_set == 0 .and. rc_unset == 0, "three I: setenv and unsetenv must succeed")
      if (allocated(error)) return
      call dense_reference(bas, ecp, dens, g_ref, err)
      call check(error,.not. err%has_error(), "three I: "//err%get_message())
      if (allocated(error)) return

      call global_logger%info("  three I: kept "//to_string(st%nkept)//" of "//to_string(st%ntask)// &
                              " tasks ("//to_string(st_all%nkept)//" unscreened), chunks "//to_string(st%nchunk)// &
                              " and "//to_string(st_chunk%nchunk))
      call global_logger%info("  three I: max|screened - unscreened| = "//to_string(maxval(abs(g - g_all)), "es10.3")// &
                              ", |chunked - whole| = "//to_string(maxval(abs(g - g_chunk)), "es10.3")// &
                              ", |gradient - dense reference| = "//to_string(maxval(abs(g - g_ref)), "es10.3"))
      call check(error, maxval(abs(g_ref)) > 1.0e-2_dp, "three I: the reference is ~0, the test is blind")
      if (allocated(error)) return
      call check(error, st%nkept > 0_int64 .and. st%nkept < st_all%nkept, "three I: the screen must drop tasks")
      if (allocated(error)) return
      call check(error, st%nchunk == 1 .and. st_chunk%nchunk == 3, "three I: one chunk, and one per ECP atom when capped")
      if (allocated(error)) return
      call check(error, maxval(abs(g - g_all)) < SCREEN_TOL, "three I: screened vs unscreened, max|diff| = "// &
                 to_string(maxval(abs(g - g_all)), "es12.4"))
      if (allocated(error)) return
      call check(error, maxval(abs(g - g_chunk)) < CHUNK_TOL, "three I: chunked vs whole, max|diff| = "// &
                 to_string(maxval(abs(g - g_chunk)), "es12.4"))
      if (allocated(error)) return
      call check(error, maxval(abs(g - g_ref)) < CONTRACT_TOL, "three I: vs the dense reference, max|diff| = "// &
                 to_string(maxval(abs(g - g_ref)), "es12.4"))
      if (allocated(error)) return

      ! The free memory decides the chunks. `st_chunk` ran last, so its
      ! reserve is the one the next call computes: a first launch has taken
      ! the local memory a device had still to take. Exactly enough for every
      ! atom at once is one chunk; one byte less is more than one.
      whole = st_chunk%chunk_reserve + st_chunk%cache_bytes
      call check(error, st_chunk%atom_bytes_max > 0_int64 .and. st_chunk%cache_bytes > st_chunk%atom_bytes_max, &
                 "three I: more than one ECP atom must hold memory, or no cap can chunk them")
      if (allocated(error)) return
      g_chunk = 0.0_dp
      call device_memory_set_cap(whole)
      call trc_ecp_gradient(bas, ecp, dens, g_chunk, err, stats=st_fit)
      g_chunk = 0.0_dp
      call device_memory_set_cap(whole - 1_int64)
      call trc_ecp_gradient(bas, ecp, dens, g_chunk, err, stats=st_short)
      call device_memory_set_cap(0_int64)
      call check(error,.not. err%has_error(), "three I: capped memory: "//err%get_message())
      if (allocated(error)) return
      call global_logger%info("  three I: capped at every atom's bytes and the reserve "//to_string(st_fit%nchunk)// &
                              " chunk(s), one byte less "//to_string(st_short%nchunk))
      call check(error, st_fit%nchunk == 1, "three I: memory for every ECP atom must give one chunk")
      if (allocated(error)) return
      call check(error, st_short%nchunk > 1, "three I: memory one byte short of every ECP atom must chunk")
      if (allocated(error)) return
      call check(error, maxval(abs(g - g_chunk)) < CHUNK_TOL, "three I: chunked by memory vs whole, max|diff| = "// &
                 to_string(maxval(abs(g - g_chunk)), "es12.4"))
      if (allocated(error)) return

      ! One byte short of the largest atom with the reserve: a status naming
      ! the sizes, not a launch that aborts.
      call device_memory_set_cap(st_chunk%chunk_reserve + st_chunk%atom_bytes_max - 1_int64)
      call trc_ecp_gradient(bas, ecp, dens, g_chunk, refused)
      call device_memory_set_cap(0_int64)
      call check(error, refused%has_error(), "three I: no room for one ECP atom with the reserve must be refused")
      if (allocated(error)) return
      call check(error, index(refused%get_message(), "device memory") > 0, &
                 "three I: the refusal names the memory: "//refused%get_message())
      call bas%release()
      call ecp%release()
   end subroutine test_screen_chunk

   subroutine ints_case(error, name)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      integer :: nat, nao, ix, c, ia, ish, a0, k
      integer, allocatable :: elements(:), ao_at(:)
      real(dp), allocatable :: r(:, :), rr(:, :), up(:, :), d(:, :, :), um(:, :), fd(:, :), ana(:, :)
      real(dp), allocatable :: u2p(:, :), u2m(:, :), dens(:, :), g(:, :)
      real(dp) :: worst, refmax, gworst
      real(dp), parameter :: H = 1.0e-3_dp

      call geometry(name, nat, elements, r)
      call trc_basis_from_json(DEF2_SVP, nat, elements, r, bas, err, force_cart=.true., &
                               overrides=trc_basis_json_overrides(allow_ecp=.true.))
      call check(error,.not. err%has_error(), name//": basis must build: "//err%get_message())
      if (allocated(error)) return
      call trc_ecp_from_json(DEF2_SVP, nat, elements, ecp, err)
      call check(error,.not. err%has_error(), name//": ECP must read: "//err%get_message())
      if (allocated(error)) return
      nao = bas%nao
      allocate (ao_at(nao))
      do ish = 1, bas%nshell
         a0 = bas%sh_ao(ish)
         ao_at(a0:a0 + (bas%sh_l(ish) + 1)*(bas%sh_l(ish) + 2)/2 - 1) = bas%sh_at(ish)
      end do

      allocate (dens(nao, nao), g(3, nat))
      call synthetic_density(dens)
      g = 0.0_dp
      call trc_ecp_gradient(bas, ecp, dens, g, err)
      call check(error,.not. err%has_error(), name//": gradient: "//err%get_message())
      if (allocated(error)) return

      worst = 0.0_dp
      refmax = 0.0_dp
      gworst = 0.0_dp
      do ix = 1, nat
         ! dU/dR_X over every ECP atom's own potential, from the bra
         ! derivatives: rows on X, the transpose for the ket, minus both for
         ! the potential sitting on X.
         allocate (ana(nao, nao))
         do c = 1, 3
            ana = 0.0_dp
            do ia = 1, nat
               if (.not. any(ecp%ecp_at(1:ecp%nrow) == ia)) cycle
               call trc_ecp_deriv_cart(bas, ecp, ia, d, err)
               call check(error,.not. err%has_error(), name//": deriv: "//err%get_message())
               if (allocated(error)) return
               do k = 1, nao
                  if (ao_at(k) == ix) then
                     ana(k, :) = ana(k, :) + d(k, :, c)
                     ana(:, k) = ana(:, k) + d(k, :, c)
                  end if
                  if (ia == ix) then
                     ana(k, :) = ana(k, :) - d(k, :, c)
                     ana(:, k) = ana(:, k) - d(k, :, c)
                  end if
               end do
               deallocate (d)
            end do
            ! Five-point central difference of the whole matrix.
            rr = r
            rr(c, ix) = r(c, ix) + H
            call ecp_u(nat, elements, rr, ecp, up, err)
            rr(c, ix) = r(c, ix) - H
            call ecp_u(nat, elements, rr, ecp, um, err)
            rr(c, ix) = r(c, ix) + 2.0_dp*H
            call ecp_u(nat, elements, rr, ecp, u2p, err)
            rr(c, ix) = r(c, ix) - 2.0_dp*H
            call ecp_u(nat, elements, rr, ecp, u2m, err)
            call check(error,.not. err%has_error(), name//": displaced matrix: "//err%get_message())
            if (allocated(error)) return
            fd = (8.0_dp*(up - um) - (u2p - u2m))/(12.0_dp*H)
            gworst = max(gworst, abs(g(c, ix) - sum(dens*fd)))
            worst = max(worst, maxval(abs(fd - ana)))
            refmax = max(refmax, maxval(abs(fd)))
         end do
         deallocate (ana)
      end do
      call global_logger%info("  "//name//": max|dU/dR analytic - difference| = "//to_string(worst, "es10.3")// &
                              " (largest element "//to_string(refmax, "es10.3")//")")
      call check(error, refmax > 1.0e-2_dp, name//": the differences are all ~0, the test is blind")
      if (allocated(error)) return
      call check(error, worst < INT_TOL, name//": dU/dR vs the difference of U, max|diff| = "//to_string(worst, "es12.4"))
      if (allocated(error)) return
      call global_logger%info("  "//name//": max|trc_ecp_gradient - contracted difference| = "// &
                              to_string(gworst, "es10.3")//" (largest component "//to_string(maxval(abs(g)), "es10.3")//")")
      call check(error, gworst < CONTRACT_TOL, name//": trc_ecp_gradient vs the contracted difference, max|diff| = "// &
                 to_string(gworst, "es12.4"))
      call bas%release()
      call ecp%release()
   end subroutine ints_case

   subroutine test_ints_hi(error)
      type(error_type), allocatable, intent(out) :: error
      call ints_case(error, "hi")
   end subroutine test_ints_hi

   subroutine test_ints_teh2(error)
      type(error_type), allocatable, intent(out) :: error
      call ints_case(error, "teh2")
   end subroutine test_ints_teh2

   subroutine test_ints_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      call ints_case(error, "hi_dimer")
   end subroutine test_ints_dimer

   subroutine run_method(nat, elements, r, method, want_grad, cart, e, grad, status)
      !! SCF (+ RI-MP2 for "mp2") at `r` through the C API; the total energy,
      !! and the gradient when asked.
      integer, intent(in) :: nat, elements(:)
      real(dp), intent(in) :: r(3, nat)
      character(len=*), intent(in) :: method
      logical, intent(in) :: want_grad, cart
      !! Cartesian AOs (def2-SVP's own are spherical).
      real(dp), intent(out) :: e, grad(3, nat)
      integer, intent(out) :: status
      type(c_ptr) :: h
      integer(c_int) :: rc, charge, mult
      real(c_double) :: en

      ! "uhf": the cation, a doublet, so the unrestricted (spin-summed) density.
      charge = merge(1_c_int, 0_c_int, method == "uhf")
      mult = merge(2_c_int, 1_c_int, method == "uhf")
      e = 0.0_dp
      grad = 0.0_dp
      status = -1
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, int(nat, c_int), real(elements, c_double), &
                                              reshape(real(r, c_double), [3*nat]), charge, mult)
      ! Before the basis load: it does not reach a basis already loaded.
      if (rc == TRC_OK .and. cart) rc = trc_set_force_cartesian(h, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, DEF2_SVP//c_null_char)
      if (rc == TRC_OK) rc = trc_set_convergence(h, 1.0e-12_c_double, 1.0e-9_c_double, 200_c_int, 10_c_int)
      if (rc == TRC_OK .and. method == "pbe") rc = trc_set_method(h, "pbe"//c_null_char, 4_c_int)
      if (rc == TRC_OK .and. method == "mp2") rc = trc_set_aux_json(h, 'basis_sets/def2-svp-rifit.json'//c_null_char)
      if (rc == TRC_OK .and. method == "mp2") rc = trc_set_rimp2(h, 0_c_int, 0_c_int)
      if (rc == TRC_OK) rc = trc_run_scf(h)
      if (rc == TRC_OK .and. method == "mp2") rc = trc_run_rimp2(h)
      if (rc == TRC_OK .and. want_grad) then
         rc = trc_run_gradient(h)
         if (rc == TRC_OK) rc = trc_gradient(h, grad)
      end if
      if (rc == TRC_OK) rc = trc_energy(h, en)
      if (rc == TRC_OK) e = en
      status = int(rc)
      rc = trc_destroy(h)
   end subroutine run_method

   subroutine grad_case(error, name, method, tol, ncoord)
      !! The analytic gradient against a central difference of the energy,
      !! `ncoord` coordinates (atom, component) spread over the molecule.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name, method
      real(dp), intent(in) :: tol
      integer, intent(in) :: ncoord
      integer :: nat, status, k, ia, c, nsel
      character(len=:), allocatable :: tag
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :), rr(:, :), g0(:, :), gd(:, :)
      real(dp) :: e0, ep, em, fd, worst
      real(dp), parameter :: H = 2.0e-3_dp

      tag = name//"/"//method
      call geometry(name, nat, elements, r)
      allocate (g0(3, nat), gd(3, nat))
      call run_method(nat, elements, r, method, .true., .false., e0, g0, status)
      call check(error, status == 0, name//"/"//method//": the SCF and gradient must run, status "//to_string(status))
      if (allocated(error)) return
      worst = 0.0_dp
      nsel = min(ncoord, 3*nat)
      do k = 1, nsel
         ! Coordinate k of the 3*nat list, visited with a stride (5) coprime
         ! to 6 and 9 so the first ones cover every atom and component.
         ia = mod((k - 1)*5, 3*nat)
         c = mod(ia, 3) + 1
         ia = ia/3 + 1
         rr = r
         rr(c, ia) = r(c, ia) + H
         call run_method(nat, elements, rr, method, .false., .false., ep, gd, status)
         call check(error, status == 0, name//"/"//method//": displaced SCF must run")
         if (allocated(error)) return
         rr(c, ia) = r(c, ia) - H
         call run_method(nat, elements, rr, method, .false., .false., em, gd, status)
         call check(error, status == 0, name//"/"//method//": displaced SCF must run")
         if (allocated(error)) return
         fd = (ep - em)/(2.0_dp*H)
         worst = max(worst, abs(fd - g0(c, ia)))
         call global_logger%info("  "//tag//": atom "//to_string(ia)//" comp "//to_string(c)//" analytic "// &
                                 to_string(g0(c, ia), "f14.8")//" difference "//to_string(fd, "f14.8"))
      end do
      call global_logger%info("  "//tag//": max|analytic - difference| = "//to_string(worst, "es10.3"))
      call check(error, worst < tol, name//"/"//method//": gradient vs the energy difference, max|diff| = "// &
                 to_string(worst, "es12.4"))
   end subroutine grad_case

   subroutine read_refs(refs)
      !! Every case of pyscf_refs/ecp_grad.ref (generate.py's gen_ecp_grad), in
      !! the file's order.
      type(pyscf_refs_t), intent(out) :: refs
      integer :: u

      call ref_require('ecp_grad.ref', u)
      call ref_read_real(u, 'E_hi_sph', refs%e_skip)
      call ref_read_matrix(u, 'G_hi_sph', 3, 2, refs%hi_sph)
      call ref_read_real(u, 'E_hi_cart', refs%e_skip)
      call ref_read_matrix(u, 'G_hi_cart', 3, 2, refs%hi_cart)
      call ref_read_real(u, 'E_teh2_sph', refs%e_skip)
      call ref_read_matrix(u, 'G_teh2_sph', 3, 3, refs%teh2_sph)
      call ref_read_real(u, 'E_teh2_cart', refs%e_skip)
      call ref_read_matrix(u, 'G_teh2_cart', 3, 3, refs%teh2_cart)
      ! check_ecp_scf's case.
      call ref_read_real(u, 'E_hi_bond_sph', refs%e_skip)
      call ref_read_matrix(u, 'G_hi_bond_sph', 3, 2, refs%hi_bond)
      call ref_read_real(u, 'E_hi_dimer_sph', refs%e_skip)
      call ref_read_matrix(u, 'G_hi_dimer_sph', 3, 4, refs%dimer)
      call ref_read_real(u, 'E_hi_uhf', refs%e_uhf)
      call ref_read_matrix(u, 'G_hi_uhf', 3, 2, refs%uhf)
      call ref_read_real(u, 'E_hi_mp2', refs%e_mp2)
      call ref_read_matrix(u, 'G_hi_mp2', 3, 2, refs%mp2)
      call ref_close(u)
   end subroutine read_refs

   subroutine test_rhf_pyscf(error)
      !! pyscf's RHF gradients (conv_tol 1e-13-ish via converge_tight, the basis
      !! and ECP read from def2-svp.json by bse_json_basis/bse_json_ecp),
      !! spherical and Cartesian, Eh/bohr, (x, y, z) per atom. Spherical HI is
      !! also check_ecp_scf's and the UHF and RI-MP2 cases' molecule.
      type(error_type), allocatable, intent(out) :: error
      type(pyscf_refs_t) :: refs

      call read_refs(refs)
      call vs_pyscf(error, "teh2", "hf", .false., refs%teh2_sph)
      if (allocated(error)) return
      call vs_pyscf(error, "hi", "hf", .true., refs%hi_cart)
   end subroutine test_rhf_pyscf

   subroutine test_rhf_pyscf_long(error)
      type(error_type), allocatable, intent(out) :: error
      type(pyscf_refs_t) :: refs

      call read_refs(refs)
      call vs_pyscf(error, "hi", "hf", .false., refs%hi_sph)
      if (allocated(error)) return
      call vs_pyscf(error, "teh2", "hf", .true., refs%teh2_cart)
      if (allocated(error)) return
      call vs_pyscf(error, "hi_dimer", "hf", .false., refs%dimer)
   end subroutine test_rhf_pyscf_long

   subroutine test_uhf_pyscf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pyscf_refs_t) :: refs

      call read_refs(refs)
      call vs_pyscf(error, "hi", "uhf", .false., refs%uhf, refs%e_uhf)
   end subroutine test_uhf_pyscf

   subroutine test_mp2_pyscf(error)
      !! The ECP term takes the relaxed density D + P: with the reference's D
      !! alone this misses by far more than PYSCF_TOL.
      type(error_type), allocatable, intent(out) :: error
      type(pyscf_refs_t) :: refs

      call read_refs(refs)
      call vs_pyscf(error, "hi", "mp2", .false., refs%mp2, refs%e_mp2)
   end subroutine test_mp2_pyscf

   subroutine vs_pyscf(error, name, method, cart, ref, e_ref)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name, method
      logical, intent(in) :: cart
      real(dp), intent(in) :: ref(:, :)
      real(dp), intent(in), optional :: e_ref
      character(len=:), allocatable :: tag
      integer :: nat, status
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :), g(:, :)
      real(dp) :: e, worst

      tag = name//"/"//method//merge("/cartesian", "/spherical", cart)
      call geometry(name, nat, elements, r)
      allocate (g(3, nat))
      call run_method(nat, elements, r, method, .true., cart, e, g, status)
      call check(error, status == 0, tag//": the SCF and gradient must run, status "//to_string(status))
      if (allocated(error)) return
      if (present(e_ref)) then
         ! The same electronic state as pyscf's, before its gradient is asked about.
         call check(error, abs(e - e_ref) < 1.0e-8_dp, tag//": energy vs pyscf, |diff| = "// &
                    to_string(abs(e - e_ref), "es12.4"))
         if (allocated(error)) return
      end if
      worst = maxval(abs(g - ref))
      call global_logger%info("  "//tag//" vs pyscf: max|diff| = "//to_string(worst, "es10.3"))
      call check(error, worst < PYSCF_TOL, tag//": gradient vs pyscf, max|diff| = "//to_string(worst, "es12.4"))
   end subroutine vs_pyscf

   subroutine test_upfront_ok(error)
      !! What --opt and --hessian ask before their first step: nothing about
      !! this molecule, basis and ECP refuses a gradient.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      integer :: nat
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :)

      call geometry("hi", nat, elements, r)
      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, int(nat, c_int), real(elements, c_double), &
                                              reshape(real(r, c_double), [3*nat]), 0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, DEF2_SVP//c_null_char)
      call check(error, rc == TRC_OK, "HI/def2-SVP setup must succeed")
      if (allocated(error)) return
      rc = trc_check_gradient_upfront(h)
      call check(error, rc == TRC_OK, "trc_check_gradient_upfront refuses HI/def2-SVP, status "//to_string(int(rc)))
      rc = trc_destroy(h)
   end subroutine test_upfront_ok

   subroutine test_rhf_hi(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, "hi", "hf", GRAD_TOL_HF, 6)
   end subroutine test_rhf_hi

   subroutine test_rhf_teh2(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, "teh2", "hf", GRAD_TOL_HF, 6)
   end subroutine test_rhf_teh2

   subroutine test_pbe_hi(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, "hi", "pbe", GRAD_TOL_KS, 4)
   end subroutine test_pbe_hi

   subroutine test_uhf_hi(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, "hi", "uhf", GRAD_TOL_HF, 6)
   end subroutine test_uhf_hi

   subroutine test_f_refused(error)
      !! Xe in def2-SVP has an f shell and an ECP: the ECP derivative raises
      !! each shell by one l, past the ECP kernels' f ceiling. The refusal is
      !! the context's own, before any SCF, and says so; an LMAX=2 build could
      !! refuse the same basis for its integral kernels first, which this
      !! message check would still tell apart.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h
      integer(c_int) :: rc
      character(len=1024) :: msg

      rc = trc_create(h)
      if (rc == TRC_OK) rc = trc_set_molecule(h, 1_c_int, [54.0_c_double], [0.0_c_double, 0.0_c_double, 0.0_c_double], &
                                              0_c_int, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_basis_json(h, DEF2_SVP//c_null_char)
      call check(error, rc == TRC_OK, "Xe/def2-SVP setup must succeed")
      if (allocated(error)) return
      rc = trc_check_gradient_upfront(h)
      call check(error, rc == TRC_ERR_UNSUPPORTED, "an f basis function over an ECP: the gradient pre-check must refuse")
      if (.not. allocated(error)) then
         msg = ""
         call trc_message_text(h, msg)
         call check(error, index(msg, "f basis functions") > 0, "the refusal names the f functions: "//trim(msg))
      end if
      rc = trc_destroy(h)
   end subroutine test_f_refused

   subroutine trc_message_text(h, text)
      !! trc_message's NUL-terminated buffer as a Fortran string.
      type(c_ptr), intent(in) :: h
      character(len=*), intent(out) :: text
      character(kind=c_char) :: buf(1024)
      integer(c_int) :: rc
      integer :: i

      buf = c_null_char
      rc = trc_message(h, buf, int(size(buf), c_int))
      text = ""
      do i = 1, min(size(buf), len(text))
         if (buf(i) == c_null_char) exit
         text(i:i) = buf(i)
      end do
   end subroutine trc_message_text

   subroutine test_mp2_hi(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, "hi", "mp2", GRAD_TOL_MP2, 4)
   end subroutine test_mp2_hi

   subroutine test_mp2_teh2(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, "teh2", "mp2", GRAD_TOL_MP2, 4)
   end subroutine test_mp2_teh2

end module check_ecp_grad_suite

program check_ecp_grad
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_ecp_grad_suite, only: collect_ecp_grad_tests
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_ecp_grad_tests, error_unit, stat, parallel=.false.)
   if (stat > 0) error stop 1
end program check_ecp_grad
