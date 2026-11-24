module check_ecp_grad_suite
   !!
   !! The analytic effective-core-potential term of the nuclear gradient.
   !!
   !! Three layers, each against something independent of the code under test:
   !!
   !!   1. The derivative integrals. For HI and a bent TeH2 in def2-SVP
   !!      (Cartesian), dU/dR_X, assembled from `trc_ecp_deriv_cart`'s
   !!      one-ECP-atom bra derivatives (basis-centre terms on the bra and, by
   !!      symmetry, the ket, plus the ECP centre by translational invariance),
   !!      against a five-point difference of `trc_ecp_matrix_cart` over every
   !!      atom X and component: the basis centres AND the ECP centre move.
   !!   2. The RHF gradient of both molecules against pyscf 2.14 (the same
   !!      BSE-JSON basis and ECP, spherical and Cartesian, tight SCF), read from
   !!      pyscf_refs/ecp_grad.ref. A doublet (the HI cation, UHF) is checked
   !!      against energy differences only.
   !!   3. The whole gradient through the C API, at the same two molecules in
   !!      def2-SVP's own spherical functions:
   !!      `trc_run_gradient` for RHF, UHF, RKS/PBE and RI-MP2 against a central
   !!      difference of the reconverged total energy, a handful of
   !!      coordinates each. Without the ECP term every one of these
   !!      gradients is off by 0.02 to 0.06 Eh/bohr, against tolerances of
   !!      1e-6 or so.
   !!
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_real, ref_read_matrix
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_basis_types, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json, trc_basis_json_overrides
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_json, only: trc_ecp_from_json
   use trc_ecp, only: trc_ecp_matrix_cart, trc_ecp_deriv_cart
   use trc_error, only: error_t
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_json, trc_set_aux_json, &
                               trc_set_convergence, trc_set_method, trc_set_rimp2, trc_set_force_cartesian, &
                               trc_run_scf, trc_run_rimp2, trc_run_gradient, trc_gradient, trc_energy, trc_message, &
                               TRC_OK, TRC_ERR_UNSUPPORTED
   use trc_capi, only: trc_check_gradient_upfront
   implicit none
   private
   public :: collect_ecp_grad_tests

   character(len=*), parameter :: DEF2_SVP = 'basis_sets/def2-svp.json'
   real(dp), parameter :: INT_TOL = 1.0e-10_dp
   !! The derivative matrix against its five-point difference, absolute: measured 8e-13.
   real(dp), parameter :: GRAD_TOL_HF = 2.0e-6_dp
   real(dp), parameter :: GRAD_TOL_KS = 2.0e-6_dp
   !! The DFT difference carries the quadrature's own noise (the grid moves
   !! with the atoms and is finite).
   real(dp), parameter :: GRAD_TOL_MP2 = 5.0e-6_dp

contains

   subroutine collect_ecp_grad_tests(testsuite)
      !! No argument: the integrals and the pyscf values (seconds). `fd`: the
      !! gradients against energy differences, which reconverge about fifty
      !! SCFs (minutes; registered LONG).
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
      else
         testsuite = [ &
                     new_unittest("dU/dR of HI vs differences of the ECP matrix", test_ints_hi), &
                     new_unittest("dU/dR of bent TeH2 vs differences of the ECP matrix", test_ints_teh2), &
                     new_unittest("RHF gradient of HI and bent TeH2 vs pyscf", test_rhf_pyscf), &
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

   subroutine ints_case(error, name)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name
      type(trc_basis_t) :: bas
      type(trc_ecp_t) :: ecp
      type(error_t) :: err
      integer :: nat, nao, ix, c, ia, ish, a0, k
      integer, allocatable :: elements(:), ao_at(:)
      real(dp), allocatable :: r(:, :), rr(:, :), up(:, :), d(:, :, :), um(:, :), fd(:, :), ana(:, :), dcen(:, :, :)
      real(dp), allocatable :: u2p(:, :), u2m(:, :)
      real(dp) :: worst, refmax
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

      worst = 0.0_dp
      refmax = 0.0_dp
      allocate (dcen(nao, nao, 3))
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

   subroutine test_rhf_pyscf(error)
      !! pyscf's RHF gradients from test/pyscf_refs/ecp_grad.ref (generate.py's
      !! gen_ecp_grad: conv_tol 1e-13-ish via converge_tight, the basis and ECP
      !! read from def2-svp.json by bse_json_basis/bse_json_ecp), spherical and
      !! Cartesian, Eh/bohr, (x, y, z) per atom. The file's fifth case is
      !! check_ecp_scf's.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: hi_sph(3, 2), hi_cart(3, 2), teh2_sph(3, 3), teh2_cart(3, 3), hi_bond(3, 2), e
      integer :: u

      call ref_require('ecp_grad.ref', u)
      call ref_read_real(u, 'E_hi_sph', e)
      call ref_read_matrix(u, 'G_hi_sph', 3, 2, hi_sph)
      call ref_read_real(u, 'E_hi_cart', e)
      call ref_read_matrix(u, 'G_hi_cart', 3, 2, hi_cart)
      call ref_read_real(u, 'E_teh2_sph', e)
      call ref_read_matrix(u, 'G_teh2_sph', 3, 3, teh2_sph)
      call ref_read_real(u, 'E_teh2_cart', e)
      call ref_read_matrix(u, 'G_teh2_cart', 3, 3, teh2_cart)
      call ref_read_real(u, 'E_hi_bond_sph', e)
      call ref_read_matrix(u, 'G_hi_bond_sph', 3, 2, hi_bond)
      call ref_close(u)
      call vs_pyscf(error, "hi", .false., hi_sph)
      if (allocated(error)) return
      call vs_pyscf(error, "teh2", .false., teh2_sph)
      if (allocated(error)) return
      call vs_pyscf(error, "hi", .true., hi_cart)
      if (allocated(error)) return
      call vs_pyscf(error, "teh2", .true., teh2_cart)
   end subroutine test_rhf_pyscf

   subroutine vs_pyscf(error, name, cart, ref)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name
      logical, intent(in) :: cart
      real(dp), intent(in) :: ref(:, :)
      character(len=:), allocatable :: tag
      integer :: nat, status
      integer, allocatable :: elements(:)
      real(dp), allocatable :: r(:, :), g(:, :)
      real(dp) :: e, worst

      tag = name//merge("/cartesian", "/spherical", cart)
      call geometry(name, nat, elements, r)
      allocate (g(3, nat))
      call run_method(nat, elements, r, "hf", .true., cart, e, g, status)
      call check(error, status == 0, tag//": the SCF and gradient must run, status "//to_string(status))
      if (allocated(error)) return
      worst = maxval(abs(g - ref))
      call global_logger%info("  "//tag//" vs pyscf: max|diff| = "//to_string(worst, "es10.3"))
      call check(error, worst < 1.0e-7_dp, tag//": RHF gradient vs pyscf, max|diff| = "//to_string(worst, "es12.4"))
   end subroutine vs_pyscf

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
