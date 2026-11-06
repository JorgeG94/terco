!
! The point group Ci (an inversion centre, nothing else): trc_sym_geom,
! trc_sym_ao and the SO-blocked SCF of trc_sym_scf on two molecules whose only
! symmetry element is the inversion. Both are staggered, anti conformers that
! pyscf 2.14.0 detects as Ci (see the comments of the .xyz files):
!   chfcl2_ci   meso-CHFCl-CHFCl, 8 atoms (C, H, F, Cl and their images)
!   n2h2f2_ci   N2H2F2, 6 atoms
! No atom is fixed by the inversion and every AO moves to another atom, with
! the sign (-1)^l, so d shells (cc-pVDZ) and the l = 0..3 set test l-parity.
!
! Detection: group Ci, irreps Ag, Au; the permutation has no fixed point, the
! symmetrised geometry is inversion symmetric to round-off, a broken geometry
! falls back to C1. AO: U^T U = I, the sign of every AO under the inversion is
! (-1)^l, U^T S U and U^T H U are block-diagonal (<= 1e-12), the dimensions
! per irrep equal pyscf's mol.symm_orb (pyscf_refs/sym_ci.ref), Cartesian and
! spherical, 6-31G* and cc-pVDZ. Registered twice, like check_sym_ao:
! `check_sym_ci` (mode "main") and `check_sym_ci_f` (mode "f", the uncontracted
! l = 0..3 set, which prints SKIP when this build's LMAX is below 3).
! SCF: the central oracle (symmetry vs C1 on the symmetrised geometry: |dE|
! <= 1e-10 Eh for HF, 1e-9 for B3LYP, orbital energies, density) for RHF,
! UHF (radical cations) and B3LYP, each against the pyscf anchors (electrons
! per irrep, orbital energies per irrep); the gradient symmetrisation.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_ci_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_sph, only: sph_cart_to_sph
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_sym_scf, only: trc_sym_scf_t, sym_grad_symmetrize
   use trc_test_basis, only: read_xyz, build_aux
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   use golden_io, only: golden_put
   implicit none
   private
   public :: collect_sym_ci, g_mode

   character(len=8) :: g_mode = 'main'
      !! "main" or "f", set from argv before run_testsuite.

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the brief's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the brief's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies, symmetry vs C1
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density, symmetry vs C1
   real(dp), parameter :: TOL_ANCHOR_HF = 1.0e-6_dp   !! terco vs pyscf, HF orbital energies
   real(dp), parameter :: TOL_ANCHOR_DFT = 1.0e-4_dp
      !! B3LYP: the two codes use different grids; measured 2.0e-5 here (|dE| 5.7e-6); a mislabelled
      !! orbital would be off by an orbital level spacing, orders of magnitude more

   integer, parameter :: NMOL = 2, NBAS = 2
   character(len=12), parameter :: MOL(NMOL) = [character(len=12) :: "chfcl2_ci", "n2h2f2_ci"]
   character(len=8), parameter :: BAS_NAME(NBAS) = [character(len=8) :: "631gs", "ccpvdz"]
   character(len=32), parameter :: BAS_FILE(NBAS) = [character(len=32) :: &
                                                      "basis_sets/6-31gs.json", "basis_sets/cc-pvdz.json"]
   character(len=4), parameter :: KIND_NAME(2) = [character(len=4) :: "cart", "sph"]

   type :: pair_t
      !! One molecule run with and without symmetry.
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_scf_result_t) :: rs, r1
      integer :: nalpha = 0, nbeta = 0
   end type pair_t

contains

   subroutine collect_sym_ci(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == 'f') then
         testsuite = [new_unittest("uncontracted l = 0..3: parity, U^T U, block-diagonal S and H", test_lmax3)]
      else
         testsuite = [ &
                     new_unittest("detection: Ci, Ag/Au, no fixed atom, symmetrisation", test_detect), &
                     new_unittest("SO basis vs pyscf dimensions, block-diagonal S and H", test_main), &
                     new_unittest("RHF meso-CHFCl-CHFCl cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_chfcl2_rhf), &
                     new_unittest("B3LYP meso-CHFCl-CHFCl 6-31G* (Cartesian): symmetry vs C1 and pyscf", test_chfcl2_b3lyp), &
                     new_unittest("UHF meso-CHFCl-CHFCl cation 6-31G* (Cartesian): symmetry vs C1 and pyscf", test_chfcl2_uhf), &
                     new_unittest("RHF N2H2F2 cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_n2h2f2_rhf), &
                     new_unittest("UHF N2H2F2 cation cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_n2h2f2_uhf), &
                     new_unittest("gradient: odd under the inversion, equals C1", test_gradient) &
                     ]
      end if
   end subroutine collect_sym_ci

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
      geom = sym_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, "Ci", "group of "//trim(name))
   end subroutine setup_geom

   subroutine check_inversion(bas, ao, label, error)
      !! The inversion sends every AO to an AO on another atom, with the sign
      !! (-1)^l of its shell (op 2 of Ci is the inversion).
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_ao_t), intent(in) :: ao
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      integer :: s, l, n, a0, a, nbad, nfix
      nbad = 0
      nfix = 0
      do s = 1, bas%nshell
         l = bas%sh_l(s)
         n = merge(2*l + 1, (l + 1)*(l + 2)/2, bas%spherical)
         a0 = merge(bas%sh_ao_sph(s), bas%sh_ao(s), bas%spherical)
         do a = a0, a0 + n - 1
            if (ao%ao_sign(a, 2) /= (-1)**l) nbad = nbad + 1
            if (ao%ao_map(a, 2) == a) nfix = nfix + 1
         end do
      end do
      call check(error, nbad == 0, label//": "//to_string(nbad)//" AOs have a sign other than (-1)^l under the inversion")
      if (allocated(error)) return
      call check(error, nfix == 0, label//": the inversion fixes an AO")
   end subroutine check_inversion

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
      call check_inversion(bas, ao, label, error)
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

   subroutine test_detect(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :), rb(:, :)
      integer :: im, i, nat, n
      real(dp) :: worst, det, orth, pos

      do im = 1, NMOL
         call setup_geom(MOL(im), elem, geom, error)
         if (allocated(error)) return
         nat = size(elem)
         call check(error, geom%detected, "Ci", trim(MOL(im))//": detected group")
         if (allocated(error)) return
         call check(error, geom%group%order, 2, trim(MOL(im))//": order")
         if (allocated(error)) return
         call check(error, trim(geom%group%irrep_names(1))//" "//trim(geom%group%irrep_names(2)), "Ag Au", &
                    trim(MOL(im))//": irrep names")
         if (allocated(error)) return
         ! the inversion permutation: an involution without a fixed atom
         n = 0
         worst = 0.0_dp
         do i = 1, nat
            if (geom%perm(i, 2) == i) n = n + 1
            if (geom%perm(geom%perm(i, 2), 2) /= i) n = n + 1
            worst = max(worst, maxval(abs(geom%r_std(:, geom%perm(i, 2)) + geom%r_std(:, i))))
            if (elem(geom%perm(i, 2)) /= elem(i)) n = n + 1
         end do
         call check(error, n == 0, trim(MOL(im))//": the inversion permutation is not a fixed-point-free involution")
         if (allocated(error)) return
         call check(error, worst <= 1.0e-14_dp, trim(MOL(im))//": symmetrised r_std is not inversion symmetric: "// &
                    to_string(worst, "es10.2"))
         if (allocated(error)) return
         ! the rotation is orthogonal with determinant +1
         det = geom%rot(1, 1)*(geom%rot(2, 2)*geom%rot(3, 3) - geom%rot(2, 3)*geom%rot(3, 2)) &
               - geom%rot(1, 2)*(geom%rot(2, 1)*geom%rot(3, 3) - geom%rot(2, 3)*geom%rot(3, 1)) &
               + geom%rot(1, 3)*(geom%rot(2, 1)*geom%rot(3, 2) - geom%rot(2, 2)*geom%rot(3, 1))
         orth = maxval(abs(matmul(geom%rot, transpose(geom%rot)) - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])))
         call check(error, abs(det - 1.0_dp) <= 1.0e-14_dp .and. orth <= 1.0e-14_dp, &
                    trim(MOL(im))//": rot is not a proper rotation")
         if (allocated(error)) return
         ! the standard frame's origin is the inversion centre: the centroid of the charge is 0
         pos = maxval(abs(sum(geom%r_std, dim=2)))
         call check(error, pos <= 1.0e-12_dp, trim(MOL(im))//": the inversion centre is not the origin")
         if (allocated(error)) return

         ! a geometry broken by 1e-3 Angstrom on one atom (well above the 1e-5 Bohr
         ! tolerance) has no symmetry; one broken by 1e-8 is symmetrised away
         call read_xyz(trim(MOL(im))//'.xyz', nat, elem, r)
         rb = r
         rb(1, 1) = rb(1, 1) + 1.0e-3_dp*1.8897261246257702_dp
         geom = sym_detect(rb, elem, err, basis_id=elem)
         call check(error, .not. err%has_error() .and. geom%group%name == "C1", &
                    trim(MOL(im))//": a broken geometry is not C1")
         if (allocated(error)) return
         rb = r
         rb(2, 3) = rb(2, 3) + 1.0e-8_dp
         geom = sym_detect(rb, elem, err, basis_id=elem)
         call check(error, .not. err%has_error() .and. geom%group%name == "Ci", &
                    trim(MOL(im))//": a 1e-8 Bohr distortion is not symmetrised to Ci")
         if (allocated(error)) return
         call check(error, geom%residual > 0.0_dp .and. geom%residual < 1.0e-7_dp, &
                    trim(MOL(im))//": residual of the symmetrisation")
         if (allocated(error)) return
      end do
   end subroutine test_detect

   subroutine test_main(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), dims(:)
      real(dp) :: ref(2)
      integer :: im, ib, ik, u, nao_ref, nao
      character(len=64) :: label

      call ref_require('sym_ci.ref', u)
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

   subroutine make_pair(error, xyz, basis_json, spherical, functional, charge, mult, p, incremental, with_c1)
      !! Detect Ci, build the basis on the standard frame, and run the SCF
      !! with symmetry and (unless `with_c1` is false) as C1.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, functional
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t), intent(inout) :: p
      logical, intent(in), optional :: incremental, with_c1
      integer :: nat, nelec, dev
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_sym_scf_t) :: sc

      dev = trc_bind_device(0)
      call read_xyz(xyz, nat, elem, r)
      p%geom = sym_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, p%geom%group%name, "Ci", "group of "//xyz)
      if (allocated(error)) return
      call trc_basis_from_json(basis_json, nat, elem, p%geom%r_std, p%bas, err, force_cart=.not. spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      p%bas%spherical = spherical
      call p%bas%to_device()
      p%ao = sym_ao_build(p%bas, p%geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (allocated(error)) return

      nelec = sum(elem) - charge
      p%nalpha = (nelec + mult - 1)/2
      p%nbeta = nelec - p%nalpha
      opts%functional = functional
      ! Convergence just above the accumulation-order noise floor, not below
      ! it. Under OpenMP (and on the device) the C1 UHF CH2Cl2+ run settles by
      ! iteration ~40 and then sits at |FDS-SDF| 1.7e-10..4e-10 and |dE|
      ! 0..9e-12 (|E| ~ 958 Eh) until max_iteration, so 1e-10/1e-12 converged
      ! only by luck (2 of 6 OpenMP runs failed). The energy is variational,
      ! so a 1e-9 commutator moves E by ~1e-18 and the 1e-10 oracle gates on
      ! dE, d eps and dD below are unaffected.
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      if (present(incremental)) opts%incremental = incremental

      sc%ao = p%ao
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%rs, sym=sc)
      call check(error, p%rs%converged, xyz//": the symmetric SCF did not converge: "//trim(p%rs%message))
      if (allocated(error)) return
      call sc%release()
      if (present(with_c1)) then
         if (.not. with_c1) return
      end if
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%r1)
      call check(error, p%r1%converged, xyz//": the C1 SCF did not converge: "//trim(p%r1%message))
   end subroutine make_pair

   subroutine compare_pair(error, p, label, tol_e)
      !! The central oracle: energy, orbital energies and density.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: tol_e
      real(dp) :: de, deps, dd
      integer :: nmo

      de = abs(p%rs%energy - p%r1%energy)
      nmo = size(p%rs%eps, 1)
      deps = maxval(abs(p%rs%eps(1:nmo - p%rs%n_lindep, :) - p%r1%eps(1:nmo - p%rs%n_lindep, :)))
      dd = maxval(abs(p%rs%dmat - p%r1%dmat))
      print '(a,a,a,es10.2,a,es10.2,a,es10.2,a,i0,a,i0,a)', '  ', label, ': |dE| ', de, '  max|d eps| ', deps, &
         '  max|dD| ', dd, '  iterations ', p%rs%iterations, ' (C1 ', p%r1%iterations, ')'
      call check(error, de <= tol_e, label//": |E(sym) - E(C1)| = "//to_string(de, "es10.2")// &
                 " exceeds "//to_string(tol_e, "es10.2"))
      if (allocated(error)) return
      call check(error, deps <= TOL_EPS, label//": orbital energies differ by "//to_string(deps, "es10.2"))
      if (allocated(error)) return
      call check(error, dd <= TOL_D, label//": densities differ by "//to_string(dd, "es10.2"))
   end subroutine compare_pair

   subroutine compare_anchor(error, p, tag, tol)
      !! Electrons per irrep and orbital energies per irrep against pyscf.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: tag
      real(dp), intent(in) :: tol
      integer :: u, nspin, ir, s, m, nocc(2), nref, k, nirr
      real(dp) :: eref_e, nelec(4), worst
      real(dp), allocatable :: ev(:), eir(:)
      integer :: cnt(2, 2)

      nirr = p%ao%nirrep
      call ref_require('sym_ci.ref', u)
      ! the file is sequential: skip to this case
      call find_case(u, tag, nspin, eref_e)
      call ref_read_vector(u, 'NELEC_'//tag, 4, nelec)
      nocc = [p%nalpha, p%nbeta]
      cnt = 0
      do s = 1, p%rs%nspin
         do m = 1, nocc(s)
            ir = p%rs%mo_irrep(m, s)
            cnt(ir, s) = cnt(ir, s) + 1
         end do
      end do
      if (p%rs%nspin == 1) cnt(:, 2) = cnt(:, 1)
      call check(error, all(cnt(:, 1) == nint(nelec(1:2))) .and. all(cnt(:, 2) == nint(nelec(3:4))), &
                 tag//": electrons per irrep "//to_string(cnt(1, 1))//","//to_string(cnt(2, 1))//" / "// &
                 to_string(cnt(1, 2))//","//to_string(cnt(2, 2))//" differ from pyscf")
      if (allocated(error)) then
         call ref_close(u)
         return
      end if
      worst = 0.0_dp
      do s = 1, nspin
         do ir = 1, nirr
            nref = p%ao%irrep_dim(ir)
            allocate (eir(nref))
            call ref_read_vector(u, 'EPS'//to_string(ir)//'_'//tag//'_'//to_string(s), nref, eir)
            ev = pack(p%rs%eps(1:p%ao%nao, s), p%rs%mo_irrep(1:p%ao%nao, s) == ir)
            call check(error, size(ev) == nref, tag//": orbitals in irrep "//to_string(ir)//" differ from pyscf")
            if (allocated(error)) then
               call ref_close(u)
               return
            end if
            do k = 1, nref
               worst = max(worst, abs(ev(k) - eir(k)))
            end do
            deallocate (eir)
         end do
      end do
      call ref_close(u)
      print '(a,a,a,es10.2,a,es10.2)', '  ', tag, ': max |eps - eps(pyscf)| per irrep ', worst, '  |E - E(pyscf)| ', &
         abs(p%rs%energy - eref_e)
      call check(error, worst <= tol, tag//": orbital energies per irrep differ from pyscf by "// &
                 to_string(worst, "es10.2"))
   end subroutine compare_anchor

   subroutine find_case(u, tag, nspin, e)
      !! Position the reader on the NSPIN record of `tag`. The generator writes
      !! cases in a fixed order, so records of the other cases are consumed.
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(out) :: nspin
      real(dp), intent(out) :: e
      character(len=512) :: line
      integer :: ios
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            print '(a)', 'check_sym_ci: no case '//tag//' in sym_ci.ref'
            stop 1
         end if
         if (index(line, 'NSPIN_'//tag//' ') == 1) exit
      end do
      read (line(len('NSPIN_'//tag) + 1:), *) nspin
      read (u, '(a)') line
      read (line(len('E_'//tag) + 1:), *) e
      ! positioned at the NELEC record, which the caller reads
   end subroutine find_case

   subroutine test_chfcl2_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'chfcl2_ci.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF CHFCl-CHFCl cc-pVDZ', TOL_E_HF)
      if (.not. allocated(error)) call check(error, p%rs%n_increments > 0 .and. p%r1%n_increments > 0, &
                                             "the default runs must have used incremental Fock builds")
      call golden_put('E_chfcl2_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'chfcl2_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_chfcl2_rhf

   subroutine test_chfcl2_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'chfcl2_ci.xyz', 'basis_sets/6-31gs.json', .false., 'b3lyp', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'B3LYP CHFCl-CHFCl 6-31G*', TOL_E_DFT)
      call golden_put('E_chfcl2_b3lyp', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'chfcl2_b3lyp', TOL_ANCHOR_DFT)
      call p%bas%release()
   end subroutine test_chfcl2_b3lyp

   subroutine test_chfcl2_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'chfcl2_ci.xyz', 'basis_sets/6-31gs.json', .false., '', 1, 2, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF CHFCl-CHFCl+ 6-31G*', TOL_E_HF)
      if (.not. allocated(error)) call compare_anchor(error, p, 'chfcl2_cation_uhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_chfcl2_uhf

   subroutine test_n2h2f2_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'n2h2f2_ci.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF N2H2F2 cc-pVDZ', TOL_E_HF)
      call golden_put('E_n2h2f2_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'n2h2f2_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_n2h2f2_rhf

   subroutine test_n2h2f2_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'n2h2f2_ci.xyz', 'basis_sets/cc-pvdz.json', .true., '', 1, 2, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF N2H2F2+ cc-pVDZ', TOL_E_HF)
      call golden_put('E_n2h2f2_cation_uhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'n2h2f2_cation_uhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_n2h2f2_uhf

   subroutine test_gradient(error)
      !! The inversion makes the gradient odd: g(image of a) = -g(a), and no atom
      !! is its own image, so no component is forced to vanish by itself.
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, 'n2h2f2_ci.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'chfcl2_ci.xyz', 'basis_sets/6-31gs.json', .true., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'n2h2f2_ci.xyz', 'basis_sets/6-31gs.json', .false., 1, 2)
      if (allocated(error)) return
      call noise_case(error)
   end subroutine test_gradient

   subroutine noise_case(error)
      !! The projection itself, on a gradient that is NOT symmetric: the
      !! result is odd under the inversion, unchanged by a second application,
      !! and keeps the part of an already odd gradient.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp), allocatable :: g(:, :), gs(:, :), g2(:, :), g3(:, :)
      real(dp) :: asym, idem, keep
      integer :: a, ic
      call setup_geom('chfcl2_ci', elem, geom, error)
      if (allocated(error)) return
      allocate (g(3, size(elem)))
      do a = 1, size(elem)
         do ic = 1, 3
            g(ic, a) = sin(1.7_dp*real(3*a + ic, dp)) + 0.3_dp
         end do
      end do
      gs = g
      call sym_grad_symmetrize(geom%group, geom%perm, gs)
      g2 = gs
      call sym_grad_symmetrize(geom%group, geom%perm, g2)
      idem = maxval(abs(g2 - gs))
      asym = 0.0_dp
      do a = 1, size(elem)
         asym = max(asym, maxval(abs(gs(:, geom%perm(a, 2)) + gs(:, a))))
      end do
      ! an odd gradient is kept as it is
      g3 = g
      do a = 1, size(elem)
         g3(:, a) = g(:, a) - g(:, geom%perm(a, 2))
      end do
      g2 = g3
      call sym_grad_symmetrize(geom%group, geom%perm, g2)
      keep = 0.0_dp
      do a = 1, size(elem)
         keep = max(keep, maxval(abs(g2(:, a) - g3(:, a))))
      end do
      call check(error, asym <= 1.0e-14_dp .and. idem <= 1.0e-14_dp .and. keep <= 1.0e-14_dp, &
                 "gradient projection: oddness "//to_string(asym, "es10.2")//", idempotence "// &
                 to_string(idem, "es10.2")//", odd part "//to_string(keep, "es10.2"))
   end subroutine noise_case

   subroutine grad_case(error, xyz, basis_json, spherical, charge, mult)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t) :: p
      type(trc_linalg_t) :: lp
      real(dp), allocatable :: graw(:, :), gsym(:, :), g1(:, :)
      real(dp) :: dg, asym, raw
      integer :: a, nat

      call make_pair(error, xyz, basis_json, spherical, '', charge, mult, p)
      if (allocated(error)) return
      nat = p%bas%natm
      allocate (graw(3, nat), gsym(3, nat), g1(3, nat))
      call lp%init(p%bas%nao)
      call trc_scf_gradient(p%bas, lp, p%rs%nspin, p%nalpha, p%nbeta, p%rs%cmo, p%rs%eps, 1.0e-12_dp, graw)
      call trc_scf_gradient(p%bas, lp, p%rs%nspin, p%nalpha, p%nbeta, p%rs%cmo, p%rs%eps, 1.0e-12_dp, gsym, &
                            sym_geom=p%geom)
      call trc_scf_gradient(p%bas, lp, p%r1%nspin, p%nalpha, p%nbeta, p%r1%cmo, p%r1%eps, 1.0e-12_dp, g1)
      call lp%release()
      call p%bas%release()

      asym = 0.0_dp
      raw = 0.0_dp
      do a = 1, nat
         asym = max(asym, maxval(abs(gsym(:, p%geom%perm(a, 2)) + gsym(:, a))))
         raw = max(raw, maxval(abs(graw(:, p%geom%perm(a, 2)) + graw(:, a))))
      end do
      dg = maxval(abs(gsym - g1))
      print '(a,a,a,es10.2,a,es10.2,a,es10.2)', '  ', xyz, ': oddness raw ', raw, '  after symmetrisation ', asym, &
         '  |sym - C1| ', dg
      call check(error, raw <= 1.0e-10_dp, xyz//": the raw gradient is not odd under the inversion: "// &
                 to_string(raw, "es10.2"))
      if (allocated(error)) return
      call check(error, asym <= 1.0e-14_dp, xyz//": symmetrised gradient is not odd: "//to_string(asym, "es10.2"))
      if (allocated(error)) return
      call check(error, dg <= 1.0e-9_dp, xyz//": symmetrised gradient differs from C1 by "//to_string(dg, "es10.2"))
   end subroutine grad_case

end module check_sym_ci_suite

program check_sym_ci
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use golden_io, only: golden_begin, golden_end
   use check_sym_ci_suite, only: collect_sym_ci, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if (trim(g_mode) == 'f' .and. LMAX < 3) then
      print '(a,i0,a)', '  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ', LMAX, '.'
      stop 0
   end if
   stat = 0
   if (trim(g_mode) /= 'f') call golden_begin('sym_ci')
   call run_testsuite(collect_sym_ci, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
   if (trim(g_mode) /= 'f') call golden_end()
end program check_sym_ci
