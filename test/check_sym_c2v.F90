module check_sym_c2v_suite
   !!
   !! The point group C2v: detection and orientation, the SO basis, the
   !! SO-blocked SCF, and the gradient, against pyscf 2.14.0 (symmetry=True) and
   !! against terco's own C1 runs.
   !!
   !! Molecules: water (planar, in the yz plane on input), formaldehyde (planar,
   !! C2 along x, in the xy plane), the NH2 radical (planar, 2B1) and nonplanar
   !! CH2Cl2 (both mirrors move atoms, in a general orientation).
   !!
   !! Orientation. pyscf puts the C2 axis on z and, for a planar molecule, the
   !! plane normal on x (detect_symm / _adjust_planar_c2v), so the out-of-plane
   !! orbitals are B1 and the in-plane ones B2: the HOMO of water is 1b1 and its
   !! 1b2 is the in-plane O-H bonding orbital. The checks that pin this:
   !!   - the frame (rows x, y, z, up to sign) and the origin that terco finds are
   !!     those pyscf prints for the four geometries;
   !!   - the irrep of the highest occupied orbital of each spin, the electrons
   !!     per irrep and the orbital energies of each irrep equal pyscf's
   !!     (pyscf_refs/sym_c2v.ref); a swap of B1 and B2 changes every one of them;
   !!   - the HOMO (b1) and the lowest b2 of water are asserted by name, also for
   !!     a copy of water that is rotated arbitrarily on input.
   !!
   !! THE CENTRAL ORACLE. A molecule on the symmetrised, standard-frame geometry
   !! of trc_sym_geom is run with the irrep blocks and as plain C1; they must
   !! agree in the total energy (<= 1e-10 Eh for HF, <= 1e-9 for B3LYP), the
   !! orbital energies and the density. RHF water (cc-pVDZ, spherical), RHF
   !! formaldehyde and RHF CH2Cl2 (6-31G*, Cartesian), UHF NH2 (cc-pVDZ), B3LYP
   !! water.
   !!
   !! The SO basis: U^T U = I, U^T S U and U^T H U block-diagonal (<= 1e-12), S
   !! invariant under the signed permutations, and the irrep dimensions equal
   !! pyscf's mol.symm_orb (pyscf_refs/sym_c2v_so_dims.ref), Cartesian and
   !! spherical, 6-31G* and cc-pVDZ.
   !!
   !! The gradient: the components that symmetry forces to zero are <= 1e-10 in
   !! the raw gradient of the symmetric run, vanish after symmetrisation, the
   !! symmetrised gradient equals the C1 one, and the projection is equivariant
   !! and idempotent on a gradient that is not symmetric.
   !!
   !! WHICH FOCK BUILD EACH CASE RUNS, DELIBERATELY. Every molecule here is below
   !! the skeleton build's size rule (`sym_run_skeleton_size_rule`): water and NH2 are
   !! 3 atoms in 2 orbits, formaldehyde 4 in 3 and CH2Cl2 5 in 3, because the C2
   !! axis and the mirrors fix most of their atoms. A default symmetric run
   !! therefore takes the FULL build with F symmetrized after it, and the skeleton
   !! build would not be exercised here at all unless it is asked for. So
   !! `make_pair` takes `skeleton`, and every case states what it wants and
   !! asserts (`trc_scf_result_t%sym_skeleton`) that it got it:
   !!   skeleton = .true.   CH2Cl2 RHF (nonplanar: both mirrors move atoms, so the
   !!                       orbit fold is nontrivial) and the NH2 radical (UHF).
   !!                       This is the group's skeleton-build coverage against
   !!                       both C1 and the pyscf anchor.
   !!   absent              the rule in charge, which for these molecules means
   !!                       the full build -- the path a user actually gets. The
   !!                       assertion is that the rule declined, so if the rule
   !!                       ever starts accepting them this test says so rather
   !!                       than silently changing what it covers.
   !!
   !! NO INTERNAL PROCEDURES: the helpers are module procedures.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient, trc_scf_grad_sym_t
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_sph, only: sph_cart_to_sph
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_run, only: SYM_SKEL_ON
   use trc_sym_scf, only: trc_sym_scf_t, sym_scf_grad_symmetrize
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none
   private
   public :: collect_sym_c2v

   integer, parameter :: NIRREP = 4
   integer, parameter :: A1 = 1, A2 = 2, B1 = 3, B2 = 4  !! the C2v irrep order of trc_sym_group (pyscf's)
   character(len=2), parameter :: IRREP_NAME(NIRREP) = [character(len=2) :: "A1", "A2", "B1", "B2"]

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the brief's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the brief's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies, symmetry vs C1
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density, symmetry vs C1
   real(dp), parameter :: TOL_ANCHOR_HF = 1.0e-6_dp   !! terco vs pyscf, HF orbital energies
   real(dp), parameter :: TOL_ANCHOR_DFT = 1.0e-5_dp  !! B3LYP: the two codes use different grids
   real(dp), parameter :: TOL_GRAD = 1.0e-10_dp

   integer, parameter :: NMOL = 4, NBAS = 2
   character(len=12), parameter :: MOL(NMOL) = [character(len=12) :: "water", "formaldehyde", "nh2", "ch2cl2"]
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

   subroutine collect_sym_c2v(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("C2v detected, frame and origin as pyscf's", test_frame), &
                  new_unittest("SO basis vs pyscf dimensions, block-diagonal S and H", test_so_basis), &
                  new_unittest("RHF water cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_water_rhf), &
                  new_unittest("RHF formaldehyde 6-31G* (Cartesian): vs C1 and pyscf", test_formaldehyde_rhf), &
                  new_unittest("RHF CH2Cl2 6-31G* (Cartesian, mirrors move atoms): vs C1 and pyscf", test_ch2cl2_rhf), &
                  new_unittest("UHF NH2 cc-pVDZ (2B1): symmetry vs C1 and pyscf", test_nh2_uhf), &
                  new_unittest("B3LYP water cc-pVDZ: symmetry vs C1 and pyscf", test_water_b3lyp), &
                  new_unittest("B1 and B2 of water, input rotated arbitrarily", test_rotated_water), &
                  new_unittest("gradient: forbidden components vanish, equals C1", test_gradient), &
                  new_unittest("fixed per-irrep occupations", test_fixed_occupations) &
                  ]
   end subroutine collect_sym_c2v

   subroutine test_frame(error)
      !! ---------------------------------------------------------------- frame
      !! pyscf 2.14.0 _symm_axes (rows x, y, z, compared up to sign) and
      !! _symm_orig (Bohr) for the four geometries as read from the .xyz files.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: ax(3, 3, NMOL), orig(3, NMOL), g(3, 3)
      integer :: im, k, nat, i, j, kk, m
      logical, allocatable :: seen(:)

      ! water and NH2: C2 along z, plane normal x
      ax(:, :, 1) = reshape([-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3], order=[2, 1])
      ! formaldehyde: C2 along x, in the xy plane, so the standard x is the input z
      ax(:, :, 2) = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], [3, 3], order=[2, 1])
      ax(:, :, 3) = ax(:, :, 1)
      ! CH2Cl2: nonplanar, x and y by the principal moments
      ax(:, :, 4) = reshape([0.59336378_dp, -0.5705068_dp, -0.56783924_dp, &
                            0.70446631_dp, 0.02677426_dp, 0.70923224_dp, &
                            -0.38941834_dp, -0.82085634_dp, 0.41778969_dp], [3, 3], order=[2, 1])
      orig(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      orig(:, 2) = [1.00006454_dp, 0.0_dp, 0.0_dp]
      orig(:, 3) = [0.0_dp, 0.0_dp, -0.26651594_dp]
      orig(:, 4) = [0.02270669_dp, -1.52508987_dp, 1.5287231_dp]
      do im = 1, NMOL
         call read_xyz(trim(MOL(im))//'.xyz', nat, elem, r)
         geom = sym_geom_detect(r, elem, err, basis_id=elem)
         call check(error, .not. err%has_error(), "sym_geom_detect "//trim(MOL(im))//": "//err%get_message())
         if (allocated(error)) return
         call check(error, geom%detected, "C2v", "detected group of "//trim(MOL(im)))
         if (allocated(error)) return
         call check(error, geom%group%name, "C2v", "reported group of "//trim(MOL(im)))
         if (allocated(error)) return
         call check(error, geom%group%order, 4, "four operations")
         if (allocated(error)) return
         call check(error, maxval(abs(geom%shift - orig(:, im))) < 1.0e-6_dp, "origin of "//trim(MOL(im)))
         if (allocated(error)) return
         do k = 1, 3
            call check(error, abs(abs(dot_product(geom%rot(k, :), ax(k, :, im))) - 1.0_dp) < 1.0e-7_dp, &
                       "axis "//to_string(k)//" of "//trim(MOL(im))//" against pyscf")
            if (allocated(error)) return
         end do
         ! proper rotation
         g = matmul(geom%rot, transpose(geom%rot))
         do i = 1, 3
            g(i, i) = g(i, i) - 1.0_dp
         end do
         call check(error, maxval(abs(g)) < 1.0e-14_dp, "orthogonal frame of "//trim(MOL(im)))
         if (allocated(error)) return
         call check(error, abs(det3(geom%rot) - 1.0_dp) < 1.0e-14_dp, "determinant +1 of "//trim(MOL(im)))
         if (allocated(error)) return
         ! every operation is a bijection of the atoms and E is the identity
         call check(error, all(geom%perm(:, 1) == [(i, i=1, nat)]), "E is the identity")
         if (allocated(error)) return
         allocate (seen(nat))
         do kk = 1, 4
            seen = .false.
            do i = 1, nat
               seen(geom%perm(i, kk)) = .true.
            end do
            call check(error, all(seen), "perm is a bijection for "//trim(MOL(im)))
            if (allocated(error)) return
            ! g x = x[perm], exactly, in the standard frame
            do i = 1, nat
               do j = 1, 3
                  m = merge(-1, 1, btest(geom%group%ops(kk), j - 1))
                  call check(error, abs(real(m, dp)*geom%r_std(j, i) - geom%r_std(j, geom%perm(i, kk))) < 1.0e-14_dp, &
                             "symmetrised coordinates are invariant")
                  if (allocated(error)) return
               end do
            end do
         end do
         deallocate (seen)
      end do
      ! the planar molecules lie in the yz plane (x = 0) of the standard frame; CH2Cl2 does not
      call read_xyz('water.xyz', nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, maxval(abs(geom%r_std(1, :))) < 1.0e-12_dp, "water lies in the yz plane")
      if (allocated(error)) return
      call read_xyz('ch2cl2.xyz', nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, maxval(abs(geom%r_std(1, :))) > 0.5_dp .and. maxval(abs(geom%r_std(2, :))) > 0.5_dp, &
                 "CH2Cl2 is not planar")
   end subroutine test_frame

   pure function det3(m) result(d)
      real(dp), intent(in) :: m(3, 3)
      real(dp) :: d
      d = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) - m(1, 2)*(m(2, 1)*m(3, 3) - m(2, 3)*m(3, 1)) &
          + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))
   end function det3

   subroutine setup_geom(name, elem, geom, error)
      !! ------------------------------------------------------------- SO basis
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
      call check(error, geom%group%name, "C2v", "group of "//trim(name))
   end subroutine setup_geom

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

   subroutine analyse(bas, geom, label, error, so_dims)
      type(trc_basis_t), intent(inout) :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable, intent(out) :: so_dims(:)
      type(trc_sym_ao_t) :: ao
      type(error_t) :: err
      real(dp), allocatable :: s(:, :), h(:, :), utu(:, :)
      integer, allocatable :: col_irrep(:)
      integer :: i, j, k, ir, n
      real(dp) :: worst_s, worst_h, worst_u, worst_r

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
      print '(2a,i5,a,es9.2,a,es9.2,a,es9.2,a,es9.2)', '  ', label, n, '  |UtU-1| ', worst_u, &
         '  off(S) ', worst_s, '  off(H) ', worst_h, '  R-inv(S) ', worst_r
      call check(error, worst_s <= 1.0e-12_dp, label//": U^T S U off-block elements")
      if (allocated(error)) return
      call check(error, worst_h <= 1.0e-12_dp, label//": U^T H U off-block elements")
      if (allocated(error)) return
      call check(error, worst_r <= 1.0e-12_dp, label//": S invariant under the signed permutations")
      if (allocated(error)) return
      so_dims = ao%irrep_dim
   end subroutine analyse

   subroutine test_so_basis(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer, allocatable :: elem(:), dims(:)
      real(dp) :: ref(NIRREP)
      integer :: im, ib, ik, u, nao_ref, nao
      character(len=64) :: label

      call ref_require('sym_c2v_so_dims.ref', u)
      do im = 1, NMOL
         call setup_geom(MOL(im), elem, geom, error)
         if (allocated(error)) return
         do ib = 1, NBAS
            do ik = 1, 2
               label = trim(MOL(im))//' '//trim(BAS_NAME(ib))//' '//KIND_NAME(ik)
               call trc_basis_from_json(BAS_FILE(ib), size(elem), elem, geom%r_std, bas, err, force_cart=ik == 1)
               call check(error, .not. err%has_error(), "basis: "//err%get_message())
               if (allocated(error)) return
               bas%spherical = ik == 2
               call analyse(bas, geom, trim(label), error, dims)
               nao = merge(bas%nao_sph, bas%nao, bas%spherical)
               call bas%release()
               if (allocated(error)) return
               call ref_read_int(u, 'NAO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib))//'_'//KIND_NAME(ik), nao_ref)
               call ref_read_vector(u, 'SO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib))//'_'//KIND_NAME(ik), NIRREP, ref)
               call check(error, nao, nao_ref, trim(label)//": AO count against pyscf")
               if (allocated(error)) return
               call check(error, size(dims), NIRREP, trim(label)//": number of irreps")
               if (allocated(error)) return
               call check(error, all(dims == nint(ref)), trim(label)//": SO dimensions (A1 A2 B1 B2) against pyscf")
               if (allocated(error)) return
            end do
         end do
      end do
      call ref_close(u)
   end subroutine test_so_basis

   subroutine run_pair(error, p, elem, bas_file, spherical, functional, charge, mult, with_c1, skeleton)
      !! ------------------------------------------------------------------ SCF
      !! Build the basis on the standard frame of `p%geom`, run the SCF with
      !! symmetry and (unless `with_c1` is false) as C1. `skeleton` picks the
      !! Fock build: present and true forces the skeleton build, absent leaves
      !! the size rule in charge (which declines every molecule of this file --
      !! see the header). Either way the path the run took is asserted, so no
      !! case can change its coverage without changing this file.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(inout) :: p
      integer, intent(in) :: elem(:)
      character(len=*), intent(in) :: bas_file, functional
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      logical, intent(in), optional :: with_c1, skeleton
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_sym_scf_t) :: sc
      logical :: want_skeleton
      integer :: nelec, dev

      dev = trc_bind_device(0)
      call trc_basis_from_json(bas_file, size(elem), elem, p%geom%r_std, p%bas, err, force_cart=.not. spherical)
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
      ! Just above the accumulation-order noise floor (see check_sym_scf_cs).
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      sc%ao = p%ao
      want_skeleton = .false.
      if (present(skeleton)) want_skeleton = skeleton
      if (want_skeleton) then
         sc%geom = p%geom
         sc%skeleton = SYM_SKEL_ON
      end if
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%rs, sym=sc)
      call check(error, p%rs%converged, "the symmetric SCF did not converge: "//trim(p%rs%message))
      if (allocated(error)) return
      if (want_skeleton) then
         call check(error, p%rs%sym_skeleton, "the skeleton Fock build was forced on and was not used")
      else
         call check(error, .not. p%rs%sym_skeleton, "the size rule should decline the skeleton build here ("// &
                    to_string(p%rs%sym_atom_orbits)//" atom orbits); it used it")
      end if
      if (allocated(error)) return
      call sc%release()
      if (present(with_c1)) then
         if (.not. with_c1) return
      end if
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%r1)
      call check(error, p%r1%converged, "the C1 SCF did not converge: "//trim(p%r1%message))
   end subroutine run_pair

   subroutine make_pair(error, xyz, bas_file, spherical, functional, charge, mult, p, with_c1, skeleton)
      !! Detect C2v, then `run_pair`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, bas_file, functional
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t), intent(inout) :: p
      logical, intent(in), optional :: with_c1, skeleton
      integer :: nat
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err

      call read_xyz(xyz, nat, elem, r)
      p%geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, p%geom%group%name, "C2v", "group of "//xyz)
      if (allocated(error)) return
      call run_pair(error, p, elem, bas_file, spherical, functional, charge, mult, with_c1, skeleton)
   end subroutine make_pair

   subroutine compare_pair(error, p, label, tol_e, key)
      !! The central oracle: energy, orbital energies and density.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: label, key
      real(dp), intent(in) :: tol_e
      real(dp) :: de, deps, dd
      integer :: nmo

      call golden_put('E_'//key, p%rs%energy)  ! the symmetric run's energy, for the golden tier
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

   subroutine find_case(u, tag, nspin, e)
      !! Position the reader on the NSPIN record of `tag`; the NELEC record
      !! is the next one.
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(out) :: nspin
      real(dp), intent(out) :: e
      character(len=512) :: line
      integer :: ios
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            print '(a)', 'check_sym_c2v: no case '//tag//' in sym_c2v.ref'
            stop 1
         end if
         if (index(line, 'NSPIN_'//tag//' ') == 1) exit
      end do
      read (line(len('NSPIN_'//tag) + 1:), *) nspin
      read (u, '(a)') line
      read (line(len('E_'//tag) + 1:), *) e
   end subroutine find_case

   subroutine compare_anchor(error, p, tag, tol)
      !! HOMO irrep, electrons per irrep and orbital energies per irrep
      !! against pyscf.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: tag
      real(dp), intent(in) :: tol
      integer :: u, nspin, ir, s, m, nocc(2), nref, k
      real(dp) :: eref_e, nelec(2*NIRREP), homo(2), worst
      real(dp), allocatable :: ev(:), eir(:)
      integer :: cnt(NIRREP, 2), my_homo(2)

      call ref_require('sym_c2v.ref', u)
      call find_case(u, tag, nspin, eref_e)
      call ref_read_vector(u, 'NELEC_'//tag, 2*NIRREP, nelec)
      call ref_read_vector(u, 'HOMO_'//tag, 2, homo)
      nocc = [p%nalpha, p%nbeta]
      cnt = 0
      my_homo = 0
      do s = 1, p%rs%nspin
         do m = 1, nocc(s)
            ir = p%rs%mo_irrep(m, s)
            cnt(ir, s) = cnt(ir, s) + 1
         end do
         if (nocc(s) > 0) my_homo(s) = p%rs%mo_irrep(nocc(s), s)
      end do
      if (p%rs%nspin == 1) then
         cnt(:, 2) = cnt(:, 1)
         my_homo(2) = my_homo(1)
      end if
      call check(error, all(cnt(:, 1) == nint(nelec(1:NIRREP))) .and. all(cnt(:, 2) == nint(nelec(NIRREP + 1:))), &
                 tag//": electrons per irrep (A1 A2 B1 B2) "//to_string(cnt(1, 1))//to_string(cnt(2, 1))// &
                 to_string(cnt(3, 1))//to_string(cnt(4, 1))//" / "//to_string(cnt(1, 2))//to_string(cnt(2, 2))// &
                 to_string(cnt(3, 2))//to_string(cnt(4, 2))//" differ from pyscf")
      if (allocated(error)) then
         call ref_close(u)
         return
      end if
      call check(error, all(my_homo == nint(homo)), tag//": the HOMO is "//IRREP_NAME(max(1, my_homo(1)))//"/"// &
                 IRREP_NAME(max(1, my_homo(2)))//", pyscf has "//IRREP_NAME(nint(homo(1)))//"/"//IRREP_NAME(nint(homo(2))))
      if (allocated(error)) then
         call ref_close(u)
         return
      end if
      worst = 0.0_dp
      do s = 1, nspin
         do ir = 1, NIRREP
            nref = p%ao%irrep_dim(ir)
            allocate (eir(nref))
            call ref_read_vector(u, 'EPS'//to_string(ir)//'_'//tag//'_'//to_string(s), nref, eir)
            ev = pack(p%rs%eps(1:p%ao%nao, s), p%rs%mo_irrep(1:p%ao%nao, s) == ir)
            call check(error, size(ev) == nref, tag//": orbitals in irrep "//IRREP_NAME(ir)//" differ from pyscf")
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

   subroutine test_water_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'water.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF water cc-pVDZ', TOL_E_HF, 'water_rhf')
      if (.not. allocated(error)) call compare_anchor(error, p, 'water_rhf', TOL_ANCHOR_HF)
      ! By name (pyscf 2.14.0): 1a1 2a1 1b2 3a1 1b1, the HOMO being the out-of-plane lone pair b1
      if (.not. allocated(error)) call check(error, all(p%rs%mo_irrep(1:5, 1) == [A1, A1, B2, A1, B1]), &
                                             "water occupied orbitals are a1 a1 b2 a1 b1")
      call p%bas%release()
   end subroutine test_water_rhf

   subroutine test_formaldehyde_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'formaldehyde.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF formaldehyde 6-31G*', TOL_E_HF, 'formaldehyde_rhf')
      if (.not. allocated(error)) call compare_anchor(error, p, 'formaldehyde_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_formaldehyde_rhf

   subroutine test_ch2cl2_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      ! The skeleton Fock build, forced: nonplanar, so both mirrors move atoms
      ! and the orbit fold is nontrivial. See the header.
      call make_pair(error, 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, skeleton=.true.)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF CH2Cl2 6-31G*, skeleton build', TOL_E_HF, 'ch2cl2_rhf')
      if (.not. allocated(error)) call check(error, p%rs%n_increments > 0 .and. p%r1%n_increments > 0, &
                                             "the default runs must have used incremental Fock builds")
      if (.not. allocated(error)) call compare_anchor(error, p, 'ch2cl2_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_ch2cl2_rhf

   subroutine test_nh2_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      ! The skeleton Fock build, forced, for an open shell. See the header.
      call make_pair(error, 'nh2.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 2, p, skeleton=.true.)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF NH2 cc-pVDZ, skeleton build', TOL_E_HF, 'nh2_uhf')
      if (.not. allocated(error)) call compare_anchor(error, p, 'nh2_uhf', TOL_ANCHOR_HF)
      ! 2B1: the unpaired electron is in the out-of-plane 1b1; the beta HOMO is 3a1
      if (.not. allocated(error)) call check(error, p%rs%mo_irrep(4, 1) == B1 .and. p%rs%mo_irrep(5, 1) == A1 .and. &
                                             p%rs%mo_irrep(4, 2) == A1, "NH2 is 2B1: alpha 1b1 is occupied, beta 1b1 is not")
      call p%bas%release()
   end subroutine test_nh2_uhf

   subroutine test_water_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'water.xyz', 'basis_sets/cc-pvdz.json', .true., 'b3lyp', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'B3LYP water cc-pVDZ', TOL_E_DFT, 'water_b3lyp')
      if (.not. allocated(error)) call compare_anchor(error, p, 'water_b3lyp', TOL_ANCHOR_DFT)
      call p%bas%release()
   end subroutine test_water_b3lyp

   subroutine test_rotated_water(error)
      !! The labels do not depend on the orientation of the input: water
      !! rotated by an arbitrary rotation and shifted gives the same orbitals
      !! (a1 a1 b2 a1 b1) and energies as the standard one.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p, q
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: rot(3, 3)
      integer :: nat, i

      call make_pair(error, 'water.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, with_c1=.false.)
      if (allocated(error)) return
      call p%bas%release()
      call read_xyz('water.xyz', nat, elem, r)
      ! a proper rotation about the axis (1, 2, -2)/3 by 1.1 rad
      rot = rotation([1.0_dp, 2.0_dp, -2.0_dp]/3.0_dp, 1.1_dp)
      r = matmul(rot, r)
      do i = 1, nat
         r(:, i) = r(:, i) + [0.4_dp, -0.7_dp, 1.3_dp]
      end do
      q%geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error() .and. q%geom%group%name == "C2v", "rotated water is C2v")
      if (allocated(error)) return
      call run_pair(error, q, elem, 'basis_sets/6-31gs.json', .false., '', 0, 1, with_c1=.false.)
      if (allocated(error)) return
      call check(error, all(q%rs%mo_irrep(1:5, 1) == [A1, A1, B2, A1, B1]), "rotated water: occupied a1 a1 b2 a1 b1")
      if (.not. allocated(error)) call check(error, abs(q%rs%energy - p%rs%energy) <= 1.0e-8_dp .and. &
                                             maxval(abs(q%rs%eps(1:q%ao%nao, 1) - p%rs%eps(1:p%ao%nao, 1))) <= 1.0e-7_dp, &
                                             "rotated water: same energy and orbital energies")
      call q%bas%release()
   end subroutine test_rotated_water

   pure function rotation(u, theta) result(m)
      !! Rotation by theta about the unit vector u (Rodrigues).
      real(dp), intent(in) :: u(3), theta
      real(dp) :: m(3, 3)
      real(dp) :: c, s
      integer :: i, j
      c = cos(theta)
      s = sin(theta)
      do j = 1, 3
         do i = 1, 3
            m(i, j) = (1.0_dp - c)*u(i)*u(j)
         end do
         m(j, j) = m(j, j) + c
      end do
      m(1, 2) = m(1, 2) - s*u(3)
      m(2, 1) = m(2, 1) + s*u(3)
      m(1, 3) = m(1, 3) + s*u(2)
      m(3, 1) = m(3, 1) - s*u(2)
      m(2, 3) = m(2, 3) - s*u(1)
      m(3, 2) = m(3, 2) + s*u(1)
   end function rotation

   subroutine test_gradient(error)
      !! -------------------------------------------------------------- gradient
      !! Water (planar), NH2 (UHF, planar) and nonplanar CH2Cl2, then the
      !! projection itself on a gradient that is not symmetric.
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, 'water.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'nh2.xyz', 'basis_sets/6-31gs.json', .false., 0, 2)
      if (allocated(error)) return
      call grad_case(error, 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call noise_case(error)
   end subroutine test_gradient

   subroutine forbidden_and_equivariance(geom, g, forbidden, asym)
      !! The largest component of `g` that symmetry forces to zero (an atom
      !! fixed by an operation has no force along an axis the operation
      !! flips), and the largest violation of g(image) = S g(atom).
      type(trc_sym_geom_t), intent(in) :: geom
      real(dp), intent(in) :: g(:, :)
      real(dp), intent(out) :: forbidden, asym
      integer :: a, k, ic
      real(dp) :: sg
      forbidden = 0.0_dp
      asym = 0.0_dp
      do a = 1, size(g, 2)
         do k = 2, geom%group%order
            do ic = 1, 3
               sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
               if (geom%perm(a, k) == a .and. sg < 0.0_dp) forbidden = max(forbidden, abs(g(ic, a)))
               asym = max(asym, abs(g(ic, geom%perm(a, k)) - sg*g(ic, a)))
            end do
         end do
      end do
   end subroutine forbidden_and_equivariance

   subroutine grad_case(error, xyz, basis_json, spherical, charge, mult)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t) :: p
      type(trc_linalg_t) :: lp
      real(dp), allocatable :: graw(:, :), gsym(:, :), g1(:, :)
      real(dp) :: forbidden, dg, asym, after, dummy
      integer :: nat

      call make_pair(error, xyz, basis_json, spherical, '', charge, mult, p)
      if (allocated(error)) return
      nat = p%bas%natm
      allocate (graw(3, nat), gsym(3, nat), g1(3, nat))
      call lp%init(p%bas%nao)
      call trc_scf_gradient(p%bas, lp, p%rs%nspin, p%nalpha, p%nbeta, p%rs%cmo, p%rs%eps, 1.0e-12_dp, graw)
      call trc_scf_gradient(p%bas, lp, p%rs%nspin, p%nalpha, p%nbeta, p%rs%cmo, p%rs%eps, 1.0e-12_dp, gsym, &
                            sym=trc_scf_grad_sym_t(geometry=p%geom))
      call trc_scf_gradient(p%bas, lp, p%r1%nspin, p%nalpha, p%nbeta, p%r1%cmo, p%r1%eps, 1.0e-12_dp, g1)
      call lp%release()
      call p%bas%release()

      call forbidden_and_equivariance(p%geom, graw, forbidden, dummy)
      call forbidden_and_equivariance(p%geom, gsym, after, asym)
      dg = maxval(abs(gsym - g1))
      print '(a,a,a,es10.2,a,es10.2,a,es10.2,a,es10.2)', '  ', xyz, ': forbidden raw ', forbidden, &
         '  after symmetrisation ', after, '  equivariance ', asym, '  |sym - C1| ', dg
      call check(error, forbidden <= TOL_GRAD, xyz//": forbidden gradient component "//to_string(forbidden, "es10.2"))
      if (allocated(error)) return
      call check(error, after <= 1.0e-14_dp, xyz//": forbidden component survives symmetrisation: "// &
                 to_string(after, "es10.2"))
      if (allocated(error)) return
      call check(error, asym <= 1.0e-14_dp, xyz//": symmetrised gradient is not equivariant: "//to_string(asym, "es10.2"))
      if (allocated(error)) return
      call check(error, dg <= 1.0e-9_dp, xyz//": symmetrised gradient differs from C1 by "//to_string(dg, "es10.2"))
   end subroutine grad_case

   subroutine noise_case(error)
      !! The projection on a gradient that is NOT symmetric: equivariant,
      !! free of forbidden components, and idempotent.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :), g(:, :), gs(:, :), g2(:, :)
      real(dp) :: asym, forbidden, idem
      integer :: a, ic, nat
      call read_xyz('ch2cl2.xyz', nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      allocate (g(3, nat))
      do a = 1, nat
         do ic = 1, 3
            g(ic, a) = sin(1.7_dp*real(3*a + ic, dp)) + 0.3_dp
         end do
      end do
      gs = g
      call sym_scf_grad_symmetrize(geom%group, geom%perm, gs)
      g2 = gs
      call sym_scf_grad_symmetrize(geom%group, geom%perm, g2)
      idem = maxval(abs(g2 - gs))
      call forbidden_and_equivariance(geom, gs, forbidden, asym)
      call check(error, asym <= 1.0e-14_dp .and. forbidden <= 1.0e-14_dp .and. idem <= 1.0e-14_dp, &
                 "gradient projection: equivariance "//to_string(asym, "es10.2")//", forbidden "// &
                 to_string(forbidden, "es10.2")//", idempotence "//to_string(idem, "es10.2"))
   end subroutine noise_case

   subroutine test_fixed_occupations(error)
      !! ------------------------------------------------------ fixed occupations
      !! pyscf's irrep_nelec for four irreps: the aufbau counts reproduce the
      !! aufbau energy; moving the b1 pair into a1 gives a higher converged
      !! state; a wrong sum is refused.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      type(trc_sym_scf_t) :: sc
      type(trc_scf_result_t) :: rf
      type(trc_scf_options_t) :: opts
      integer :: m, cnt(NIRREP)

      call make_pair(error, 'water.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, with_c1=.false.)
      if (allocated(error)) return
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      cnt = 0
      do m = 1, p%nalpha
         cnt(p%rs%mo_irrep(m, 1)) = cnt(p%rs%mo_irrep(m, 1)) + 1
      end do
      call check(error, all(cnt == [3, 0, 1, 1]), "water aufbau counts per irrep are 3 0 1 1")
      if (allocated(error)) then
         call p%bas%release()
         return
      end if

      sc%ao = p%ao
      allocate (sc%irrep_nocc(NIRREP, 1))
      sc%irrep_nocc(:, 1) = cnt
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rf, sym=sc)
      call check(error, rf%converged, "fixed aufbau occupations: not converged: "//trim(rf%message))
      if (.not. allocated(error)) call check(error, abs(rf%energy - p%rs%energy) <= TOL_E_HF, &
         "fixed aufbau occupations: energy differs from aufbau by "//to_string(abs(rf%energy - p%rs%energy), "es10.2"))
      if (allocated(error)) then
         call p%bas%release()
         return
      end if

      ! b1 -> a1: a higher state (the 1b1 -> 4a1 double excitation)
      sc%irrep_nocc(:, 1) = [cnt(A1) + 1, cnt(A2), cnt(B1) - 1, cnt(B2)]
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rf, sym=sc)
      call check(error, rf%converged, "excited occupations: not converged: "//trim(rf%message))
      if (.not. allocated(error)) call check(error, rf%energy > p%rs%energy + 1.0e-2_dp, &
         "excited occupations: energy "//to_string(rf%energy)//" is not above the ground state "// &
         to_string(p%rs%energy))
      if (allocated(error)) then
         call p%bas%release()
         return
      end if

      ! a wrong sum is refused, not run
      sc%irrep_nocc(:, 1) = [cnt(A1), cnt(A2), cnt(B1) + 1, cnt(B2)]
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rf, sym=sc)
      call check(error, .not. rf%converged .and. index(rf%message, "irrep_nocc") > 0, &
                 "a wrong irrep_nocc sum must be refused, got: "//trim(rf%message))
      call p%bas%release()
   end subroutine test_fixed_occupations

end module check_sym_c2v_suite

program check_sym_c2v
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_sym_c2v_suite, only: collect_sym_c2v
   use golden_io, only: golden_begin, golden_end
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, is
   stat = 0
   call golden_begin("sym_c2v")
   testsuites = [new_testsuite("sym_c2v", collect_sym_c2v)]
   do is = 1, size(testsuites)
      write (error_unit, '(a)') "Testing: "//testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
   end do
   call golden_end()
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_c2v
