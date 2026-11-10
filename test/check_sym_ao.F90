!
! trc_sym_ao: AO signed permutations and the SO basis, on terco's own 1e
! integrals.
!
! For HOCl (planar, so every atom is its own mirror image: only the signs are
! exercised) and nonplanar CHFCl2 (the mirror swaps the two Cl), both Cs (C2v
! has its own check, check_sym_c2v), in 6-31G* and cc-pVDZ, in the Cartesian
! and the spherical AO space, with the basis built on the standard frame of trc_sym_geom:
!   - U^T U = I, and the irrep dimensions add up to the AO count;
!   - the dimensions equal those of pyscf 2.14.0's mol.symm_orb
!     (pyscf_refs/sym_so_dims.ref);
!   - U^T S U and U^T H U (H = T + V) are block-diagonal in the irreps,
!     off-block elements <= 1e-12 (scaled by the largest element for the
!     uncontracted set, whose exponents run to 1e3);
!   - S is invariant under every signed permutation.
! Also: a decontracted (primitive) basis keeps the shell map; equivalent
! atoms with different shells, and a basis not on the standard frame, are
! refused.
!
! Registered twice, as check_nosym is: `check_sym_ao` (mode "main") and
! `check_sym_ao_f` (mode "f", the uncontracted l = 0..3 set, which prints
! SKIP and exits when this build's LMAX is below 3). The Cartesian component
! order is whatever TERCO_CART_ORDER this build has; the tests run under
! both orders (see the stage report).
!
module check_sym_ao_suite
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_basis_json, only: trc_basis_from_json
   use trc_decontract, only: decon_t, decontract_basis, decontract_release
   use trc_constants, only: ANGSTROM_TO_BOHR
   use trc_error, only: error_t
   use trc_sph, only: sph_cart_to_sph
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_test_basis, only: read_xyz, build_aux
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   implicit none
   private
   public :: collect_sym_ao, g_mode

   character(len=8) :: g_mode = 'main'
      !! "main" or "f", set from argv before run_testsuite.

   integer, parameter :: NMOL = 2, NBAS = 2
   character(len=8), parameter :: MOL(NMOL) = [character(len=8) :: "hocl", "chfcl2"]
   character(len=8), parameter :: BAS_NAME(NBAS) = [character(len=8) :: "631gs", "ccpvdz"]
   character(len=32), parameter :: BAS_FILE(NBAS) = [character(len=32) :: &
                                                      "basis_sets/6-31gs.json", "basis_sets/cc-pvdz.json"]
   character(len=4), parameter :: KIND_NAME(2) = [character(len=4) :: "cart", "sph"]

contains

   subroutine collect_sym_ao(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == 'f') then
         testsuite = [new_unittest("uncontracted l = 0..3: U^T U, block-diagonal S and H", test_lmax3)]
      else
         testsuite = [ &
                     new_unittest("SO basis vs pyscf dimensions, block-diagonal S and H", test_main), &
                     new_unittest("decontracted basis keeps the shell map", test_decontracted), &
                     new_unittest("different shells on equivalent atoms are refused", test_refuse_shells), &
                     new_unittest("basis off the standard frame is refused", test_refuse_frame) &
                     ]
      end if
   end subroutine collect_sym_ao

   subroutine setup_geom(name, elem, geom, error)
      !! Read molecule `name`, detect its group, and return the geometry result.
      character(len=*), intent(in) :: name
      integer, allocatable, intent(out) :: elem(:)
      type(trc_sym_geom_t), intent(out) :: geom
      type(error_type), allocatable, intent(out) :: error
      integer :: nat
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      call read_xyz(trim(name)//'.xyz', nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, "Cs", "group of "//trim(name))
   end subroutine setup_geom

   subroutine make_basis(path, elem, geom, spherical, bas, error)
      character(len=*), intent(in) :: path
      integer, intent(in) :: elem(:)
      type(trc_sym_geom_t), intent(in) :: geom
      logical, intent(in) :: spherical
      type(trc_basis_t), intent(out) :: bas
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      call trc_basis_from_json(path, size(elem), elem, geom%r_std, bas, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      ! A Cartesian-tagged file (6-31G*) is also used spherically: the flag
      ! is all that selects the AO space, nao_sph is always computed.
      bas%spherical = spherical
   end subroutine make_basis

   subroutine one_electron(bas, s, h)
      !! S and H = T + V in the AO space of `bas` (spherical when it says so).
      type(trc_basis_t), intent(inout) :: bas
      real(dp), allocatable, intent(out) :: s(:, :), h(:, :)
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: sc(:, :), tc(:, :), vc(:, :), hc(:, :)
      integer :: n, ns
      n = bas%nao
      ns = merge(bas%nao_sph, bas%nao, bas%spherical)
      call pl%build(bas, 1.0e-14_dp)
      call bas%to_device(); call pl%to_device()
      allocate (sc(n, n), tc(n, n), vc(n, n), hc(n, n), s(ns, ns), h(ns, ns))
      call trc_1e(bas, pl, sc, tc, vc)
      hc = tc + vc
      if (bas%spherical) then
         call sph_cart_to_sph(bas, sc, s)
         call sph_cart_to_sph(bas, hc, h)
      else
         s = sc
         h = hc
      end if
      call pl%release()
   end subroutine one_electron

   subroutine analyse(bas, geom, label, scaled, error, so_dims)
      !! The identities of the header on one basis; returns the SO dimensions.
      type(trc_basis_t), intent(inout) :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      character(len=*), intent(in) :: label
      logical, intent(in) :: scaled
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable, intent(out), optional :: so_dims(:)
      type(trc_sym_ao_t) :: ao
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), utu(:, :)
      integer, allocatable :: col_irrep(:)
      integer :: i, j, k, ir, n
      real(dp) :: worst_s, worst_h, worst_u, worst_r, tol, scale

      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), label//": sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      n = ao%nao
      call check(error, sum(ao%irrep_dim), n, label//": irrep dimensions add up to nao")
      if (allocated(error)) return

      utu = matmul(transpose(ao%u), ao%u)
      do i = 1, n
         utu(i, i) = utu(i, i) - 1.0_dp
      end do
      worst_u = maxval(abs(utu))
      call check(error, worst_u < 1.0e-13_dp, label//": U^T U = I")
      if (allocated(error)) return

      allocate (col_irrep(n))
      do ir = 1, ao%nirrep
         col_irrep(ao%irrep_off(ir) + 1:ao%irrep_off(ir) + ao%irrep_dim(ir)) = ir
      end do

      call one_electron(bas, s, h)
      ! S invariant under the signed permutations: S(map a, map b) sa sb = S(a, b)
      worst_r = 0.0_dp
      do k = 1, ao%nop
         do j = 1, n
            do i = 1, n
               worst_r = max(worst_r, abs(real(ao%ao_sign(i, k)*ao%ao_sign(j, k), dp)* &
                                          s(ao%ao_map(i, k), ao%ao_map(j, k)) - s(i, j)))
            end do
         end do
      end do

      worst_s = offblock(ao%u, s, col_irrep)
      worst_h = offblock(ao%u, h, col_irrep)
      scale = merge(max(1.0_dp, maxval(abs(h))), 1.0_dp, scaled)
      tol = 1.0e-12_dp*scale
      print '(2a,i5,a,es9.2,a,es9.2,a,es9.2,a,es9.2)', '  ', label, n, '  |UtU-1| ', worst_u, &
         '  off(S) ', worst_s, '  off(H) ', worst_h, '  R-inv(S) ', worst_r
      call check(error, worst_s <= tol, label//": U^T S U off-block elements")
      if (allocated(error)) return
      call check(error, worst_h <= tol, label//": U^T H U off-block elements")
      if (allocated(error)) return
      call check(error, worst_r <= tol, label//": S invariant under the signed permutations")
      if (allocated(error)) return
      if (present(so_dims)) so_dims = ao%irrep_dim
   end subroutine analyse

   function offblock(u, a, col_irrep) result(worst)
      !! Largest element of U^T A U between columns of different irreps.
      real(dp), intent(in) :: u(:, :), a(:, :)
      integer, intent(in) :: col_irrep(:)
      real(dp) :: worst
      real(dp), allocatable :: m(:, :)
      integer :: i, j
      m = matmul(transpose(u), matmul(a, u))
      worst = 0.0_dp
      do j = 1, size(m, 2)
         do i = 1, size(m, 1)
            if (col_irrep(i) /= col_irrep(j)) worst = max(worst, abs(m(i, j)))
         end do
      end do
   end function offblock

   subroutine test_main(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), dims(:)
      real(dp) :: ref(2)
      integer :: im, ib, ik, u, nao_ref, nao
      character(len=64) :: label

      call ref_require('sym_so_dims.ref', u)
      do im = 1, NMOL
         call setup_geom(MOL(im), elem, geom, error)
         if (allocated(error)) return
         do ib = 1, NBAS
            do ik = 1, 2
               label = trim(MOL(im))//' '//trim(BAS_NAME(ib))//' '//KIND_NAME(ik)
               call make_basis(BAS_FILE(ib), elem, geom, ik == 2, bas, error)
               if (allocated(error)) return
               call analyse(bas, geom, trim(label), .false., error, dims)
               nao = merge(bas%nao_sph, bas%nao, bas%spherical)
               call bas%release()
               if (allocated(error)) return
               call ref_read_int(u, 'NAO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib))//'_'//KIND_NAME(ik), nao_ref)
               call ref_read_vector(u, 'SO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib))//'_'//KIND_NAME(ik), 2, ref)
               call check(error, nao, nao_ref, trim(label)//": AO count against pyscf")
               if (allocated(error)) return
               call check(error, size(dims), 2, trim(label)//": number of irreps")
               if (allocated(error)) return
               call check(error, all(dims == nint(ref)), trim(label)//": SO dimensions against pyscf")
               if (allocated(error)) return
            end do
         end do
      end do
      call ref_close(u)
   end subroutine test_main

   subroutine test_decontracted(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas, pb
      type(trc_sym_geom_t) :: geom
      type(decon_t) :: m
      integer, allocatable :: elem(:)
      call setup_geom('hocl', elem, geom, error)
      if (allocated(error)) return
      call make_basis('basis_sets/cc-pvdz.json', elem, geom, .false., bas, error)
      if (allocated(error)) return
      call decontract_basis(bas, pb, m)
      call check(error, m%active, "cc-pVDZ HOCl has something to decontract")
      if (allocated(error)) return
      call analyse(pb, geom, 'hocl cc-pVDZ primitive', .false., error)
      call decontract_release(m)
      call pb%release(); call bas%release()
   end subroutine test_decontracted

   subroutine test_refuse_shells(error)
      !! CH2Cl2 (nonplanar, C2v): the operations swap pairs of equivalent
      !! atoms, so a change to one of them must be refused. (Planar HOCl lies
      !! in the mirror plane, so every atom is its own image.)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(error_t) :: err
      integer, parameter :: elem(5) = [6, 17, 17, 1, 1]
      real(dp) :: r(3, 5)
      integer :: s, ia, ja
      r = reshape([0.3000000000_dp, -0.2000000000_dp, 0.5000000000_dp, &
                   -0.9303249111_dp, -0.1525786544_dp, 1.7202671139_dp, &
                   0.7904300607_dp, -1.8070483858_dp, 0.0735333056_dp, &
                   1.1684143839_dp, 0.3327600172_dp, 0.8721870790_dp, &
                   -0.0855356395_dp, 0.2851018405_dp, -0.3902463002_dp], [3, 5])*ANGSTROM_TO_BOHR
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect CH2Cl2")
      if (allocated(error)) return
      call make_basis('basis_sets/6-31gs.json', elem, geom, .false., bas, error)
      if (allocated(error)) return
      ! the whole thing works with identical shells ...
      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), "CH2Cl2 SO basis: "//err%get_message())
      if (allocated(error)) return
      call check(error, any(ao%ao_map(:, 2) /= [(s, s=1, ao%nao)]), "the mirror moves some AOs")
      if (allocated(error)) return
      ! ... and nudging one exponent on an atom that has a distinct image is refused
      ia = 0
      do s = 1, 5
         if (geom%perm(s, 2) /= s) then
            ia = s
            exit
         end if
      end do
      call check(error, ia > 0, "the mirror moves an atom")
      if (allocated(error)) return
      ja = ia
      do s = 1, bas%nshell
         if (bas%sh_at(s) == ja) exit
      end do
      bas%sh_e(1, s) = bas%sh_e(1, s)*(1.0_dp + 1.0e-6_dp)
      ao = sym_ao_build(bas, geom, err)
      call check(error, err%has_error(), "perturbed exponent on an equivalent atom is refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "exponents") > 0, &
                 "refusal names the exponents: "//err%get_message())
      call bas%release()
   end subroutine test_refuse_shells

   subroutine test_refuse_frame(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(error_t) :: err
      integer, allocatable :: elem(:)
      integer :: nat
      real(dp), allocatable :: r(:, :)
      call setup_geom('hocl', elem, geom, error)
      if (allocated(error)) return
      call read_xyz('hocl.xyz', nat, elem, r)  ! the input frame, not the standard one
      call trc_basis_from_json('basis_sets/6-31gs.json', nat, elem, r, bas, err)
      call check(error, .not. err%has_error(), "basis in input frame")
      if (allocated(error)) return
      ao = sym_ao_build(bas, geom, err)
      call check(error, err%has_error(), "basis off the standard frame is refused")
      call bas%release()
   end subroutine test_refuse_frame

   subroutine test_lmax3(error)
      !! Uncontracted l = 0..3 shells on every atom, Cartesian and spherical.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer :: im, ik, nsh, maxnp
      do im = 1, NMOL
         call setup_geom(MOL(im), elem, geom, error)
         if (allocated(error)) return
         do ik = 1, 2
            call build_aux(size(elem), elem, geom%r_std, 3, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
            call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, size(elem), real(elem, dp), geom%r_std, maxnp)
            bas%spherical = (ik == 2)
            call analyse(bas, geom, trim(MOL(im))//' l<=3 '//KIND_NAME(ik), .true., error)
            call bas%release()
            if (allocated(error)) return
         end do
      end do
   end subroutine test_lmax3

end module check_sym_ao_suite

program check_sym_ao
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use check_sym_ao_suite, only: collect_sym_ao, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if (trim(g_mode) == 'f' .and. LMAX < 3) then
      print '(a,i0,a)', '  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ', LMAX, '.'
      stop 0
   end if
   stat = 0
   call run_testsuite(collect_sym_ao, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_ao
