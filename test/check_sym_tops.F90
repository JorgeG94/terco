module check_sym_tops_suite
   !!
   !! Symmetric and spherical tops, linear molecules and atoms (trc_sym_geom +
   !! trc_sym_tops). The group and frame are those pyscf 2.14.0 gives with
   !! symmetry=True (its `detect_symm` full group, reduced by `get_subgroup` to the
   !! Abelian `mol.groupname`). Molecules (sources in the .xyz comment lines; each
   !! non-linear one and each linear one also rotated and shifted: *_rot):
   !!   tops_benzene   D6h -> D2h      tops_nh3  C3v -> Cs      tops_ch4  Td -> D2
   !!   tops_sf6       Oh  -> D2h      tops_n2   Dooh -> D2h    tops_co   Coov -> C2v
   !!   tops_hcn       Coov -> C2v     o2 (UHF triplet) Dooh -> D2h
   !!   tops_ne, o_atom   SO3 -> D2h (pyscf: D2h for a Cartesian basis, SO3 for a spherical one)
   !!
   !!   - detection: group, origin and frame against pyscf. The frame is compared up
   !!     to the sign of each axis. Three kinds: (1) all three axes must equal
   !!     pyscf's (everything deterministic in pyscf's search); (2) only z (linear
   !!     molecules: pyscf takes x and y from LAPACK, arbitrary in the plane); (3) the
   !!     two frames are related by an operation of the molecule (NH3 rotated:
   !!     `search_mirrorx` picks one of three equivalent mirrors from the sign of
   !!     LAPACK's unique axis), checked by comparing the atoms in the two frames.
   !!     The labels and SO dimensions of two such frames are the same;
   !!   - the symmetrised geometry is invariant under every operation to 1e-14;
   !!   - SO dimensions per irrep (cart 6-31G*, sph cc-pVDZ) against pyscf; for the
   !!     linear molecules and atoms pyscf's anchor is the Abelian group requested on
   !!     the aligned file (the spherical basis would use Dooh/Coov/SO3 irreps);
   !!     U^T U = I, U^T S U and U^T H U block-diagonal to 1e-12;
   !!   - THE CENTRAL ORACLE (symmetry vs C1 on the symmetrised geometry) for RHF,
   !!     UHF (NH3+, O2 triplet, O atom) and B3LYP (NH3, grid level 5), electrons per
   !!     irrep and orbital energies per irrep against pyscf (the O atom: the sorted
   !!     orbital energies only, since which degenerate p orbital the beta electron
   !!     takes is arbitrary);
   !!   - requested subgroups C2v and Cs of benzene, in pyscf's frame;
   !!   - the symmetrised gradient; `auto` and explicit group choice.
   !!
   !! NO INTERNAL PROCEDURES: the helpers are module procedures.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use pic_ascii, only: to_lower
   use pic_sorting, only: sort
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
   use trc_sym_run, only: trc_sym_request_t, sym_run_choose_group, SYM_AUTO, SYM_EXPLICIT
   use trc_sym_scf, only: trc_sym_scf_t, sym_scf_grad_symmetrize
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   use golden_io, only: golden_put
   implicit none
   private
   public :: collect_sym_tops

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the brief's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the brief's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies, symmetry vs C1
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density, symmetry vs C1
   real(dp), parameter :: TOL_ANCHOR_HF = 1.0e-6_dp   !! terco vs pyscf, HF orbital energies
   real(dp), parameter :: TOL_ANCHOR_DFT = 1.0e-5_dp
      !! B3LYP orbital energies against pyscf: different grids; terco at level 5 (see the plan, "Known issues")
   integer, parameter :: DFT_GRID_LEVEL = 5

   integer, parameter :: NMOL = 18
   character(len=16), parameter :: MOL(NMOL) = [character(len=16) :: &
      "tops_benzene", "tops_benzene_rot", "tops_nh3", "tops_nh3_rot", "tops_ch4", "tops_ch4_rot", &
      "tops_sf6", "tops_sf6_rot", "tops_n2", "tops_n2_rot", "tops_co", "tops_co_rot", "tops_hcn", &
      "tops_hcn_rot", "o2", "tops_o2_rot", "tops_ne", "o_atom"]
   character(len=3), parameter :: GROUP(NMOL) = [character(len=3) :: &
      "D2h", "D2h", "Cs", "Cs", "D2", "D2", "D2h", "D2h", "D2h", "D2h", "C2v", "C2v", "C2v", "C2v", "D2h", "D2h", &
      "D2h", "D2h"]
   integer, parameter :: FKIND(NMOL) = [1, 1, 1, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1]
      !! 1: all axes pyscf's; 2: only z (linear); 3: related by a symmetry operation
   character(len=16), parameter :: DIMREF(NMOL) = [character(len=16) :: &
      "tops_benzene", "tops_benzene_rot", "tops_nh3", "tops_nh3_rot", "tops_ch4", "tops_ch4_rot", &
      "tops_sf6", "tops_sf6_rot", "tops_n2", "tops_n2", "tops_co", "tops_co", "tops_hcn", &
      "tops_hcn", "o2", "o2", "tops_ne", "o_atom"]
      !! the molecule whose pyscf SO dimensions apply (the aligned file for a rotated linear molecule)
   integer, parameter :: NBAS = 2
   character(len=8), parameter :: BAS_NAME(NBAS) = [character(len=8) :: "631gs", "ccpvdz"]
   character(len=32), parameter :: BAS_FILE(NBAS) = [character(len=32) :: &
                                                      "basis_sets/6-31gs.json", "basis_sets/cc-pvdz.json"]
   character(len=4), parameter :: KIND_NAME(2) = [character(len=4) :: "cart", "sph"]
   character(len=*), parameter :: REF = 'sym_tops.ref'

   type :: pair_t
      !! One molecule run with and without symmetry.
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_scf_result_t) :: rs, r1
      integer :: nalpha = 0, nbeta = 0
   end type pair_t

contains

   subroutine collect_sym_tops(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("detection: group, origin and pyscf frame (aligned and rotated)", test_detection), &
                  new_unittest("SO basis vs pyscf dimensions, block-diagonal S and H", test_salc), &
                  new_unittest("RHF benzene 6-31G* (D2h, rotated): symmetry vs C1 and pyscf", test_benzene), &
                  new_unittest("RHF NH3 cc-pVDZ (Cs, rotated): symmetry vs C1 and pyscf", test_nh3_rhf), &
                  new_unittest("UHF NH3+ cc-pVDZ (Cs): symmetry vs C1 and pyscf", test_nh3_uhf), &
                  new_unittest("B3LYP NH3 cc-pVDZ (Cs, rotated): symmetry vs C1 and pyscf", test_nh3_b3lyp), &
                  new_unittest("RHF CH4 cc-pVDZ (D2, rotated): symmetry vs C1 and pyscf", test_ch4), &
                  new_unittest("RHF SF6 6-31G* (D2h from Oh, rotated): symmetry vs C1 and pyscf", test_sf6), &
                  new_unittest("RHF N2 cc-pVDZ (D2h, rotated): symmetry vs C1 and pyscf", test_n2), &
                  new_unittest("RHF CO cc-pVDZ (C2v, rotated): symmetry vs C1 and pyscf", test_co), &
                  new_unittest("RHF HCN cc-pVDZ (C2v, rotated): symmetry vs C1 and pyscf", test_hcn), &
                  new_unittest("UHF O2 triplet cc-pVDZ (D2h, rotated): symmetry vs C1 and pyscf", test_o2), &
                  new_unittest("RHF Ne cc-pVDZ (D2h atom): symmetry vs C1 and pyscf", test_ne), &
                  new_unittest("UHF O atom cc-pVDZ (D2h atom): symmetry vs C1 and pyscf", test_o_atom), &
                  new_unittest("requested C2v and Cs of benzene: pyscf frame, dimensions, energy", test_benzene_sub), &
                  new_unittest("gradient: equivariant, equals C1, projection", test_gradient), &
                  new_unittest("auto and explicit group choice", test_choose) &
                  ]
   end subroutine collect_sym_tops

   subroutine goto_tag(u, tag)
      !! Position the reader on the record `tag` (anywhere in the file).
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      character(len=512) :: line
      integer :: ios
      rewind (u)
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            print '(a)', 'check_sym_tops: no record '//tag//' in '//REF
            stop 1
         end if
         if (index(line, tag//' ') == 1) exit
      end do
      backspace (u)
   end subroutine goto_tag

   subroutine pyscf_frame(tag, ax, orig)
      !! pyscf's frame (rows x, y, z) and origin from the reference.
      character(len=*), intent(in) :: tag
      real(dp), intent(out) :: ax(3, 3), orig(3)
      real(dp) :: v(9)
      integer :: u
      call ref_require(REF, u)
      call goto_tag(u, 'AXES_'//trim(tag))
      call ref_read_vector(u, 'AXES_'//trim(tag), 9, v)
      ax = transpose(reshape(v, [3, 3]))
      call goto_tag(u, 'ORIG_'//trim(tag))
      call ref_read_vector(u, 'ORIG_'//trim(tag), 3, orig)
      call ref_close(u)
   end subroutine pyscf_frame

   function frames_equivalent(r, elem, rot1, sh1, rot2, sh2) result(ok)
      !! Whether the atoms, in the frame (rot1, sh1) with some axis signs flipped, are
      !! the atoms of the frame (rot2, sh2): the two frames differ by a symmetry
      !! operation of the molecule (and axis signs), so every irrep holds the same functions.
      real(dp), intent(in) :: r(:, :), rot1(3, 3), sh1(3), rot2(3, 3), sh2(3)
      integer, intent(in) :: elem(:)
      logical :: ok
      real(dp) :: q1(3), q2(3), sg(3)
      integer :: mask, i, j, d
      logical :: found, all_found
      ok = .false.
      do mask = 0, 7
         do d = 1, 3
            sg(d) = merge(-1.0_dp, 1.0_dp, btest(mask, d - 1))
         end do
         all_found = .true.
         do i = 1, size(elem)
            q1 = sg*matmul(rot1, r(:, i) - sh1)
            found = .false.
            do j = 1, size(elem)
               if (elem(j) /= elem(i)) cycle
               q2 = matmul(rot2, r(:, j) - sh2)
               if (maxval(abs(q1 - q2)) < 1.0e-5_dp) found = .true.
            end do
            if (.not. found) all_found = .false.
         end do
         if (all_found) then
            ok = .true.
            return
         end if
      end do
   end function frames_equivalent

   subroutine detect_mol(name, nat, elem, r, geom, error, requested)
      character(len=*), intent(in) :: name
      integer, intent(out) :: nat
      integer, allocatable, intent(out) :: elem(:)
      real(dp), allocatable, intent(out) :: r(:, :)
      type(trc_sym_geom_t), intent(out) :: geom
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in), optional :: requested
      type(error_t) :: err
      call read_xyz(trim(name)//'.xyz', nat, elem, r)
      if (present(requested)) then
         geom = sym_geom_detect(r, elem, err, basis_id=elem, requested=requested)
      else
         geom = sym_geom_detect(r, elem, err, basis_id=elem)
      end if
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
   end subroutine detect_mol

   subroutine test_detection(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: ax(3, 3), orig(3), dev, sg
      integer :: im, i, j, k, nat, ic
      do im = 1, NMOL
         call detect_mol(MOL(im), nat, elem, r, geom, error)
         if (allocated(error)) return
         call check(error, geom%detected, GROUP(im), "detected group of "//trim(MOL(im)))
         if (allocated(error)) return
         call check(error, geom%group%name, GROUP(im), "group of "//trim(MOL(im)))
         if (allocated(error)) return
         call check(error, geom%degenerate, trim(MOL(im))//": flagged as a top, linear molecule or atom")
         if (allocated(error)) return
         call pyscf_frame(MOL(im), ax, orig)
         call check(error, maxval(abs(geom%shift - orig)) < 1.0e-6_dp, "origin of "//trim(MOL(im)))
         if (allocated(error)) return
         call check(error, maxval(abs(matmul(geom%rot, transpose(geom%rot)) - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))) &
                    < 1.0e-12_dp, trim(MOL(im))//": rotation is orthogonal")
         if (allocated(error)) return
         call check(error, abs(geom%rot(1, 1)*(geom%rot(2, 2)*geom%rot(3, 3) - geom%rot(2, 3)*geom%rot(3, 2)) &
                               - geom%rot(1, 2)*(geom%rot(2, 1)*geom%rot(3, 3) - geom%rot(2, 3)*geom%rot(3, 1)) &
                               + geom%rot(1, 3)*(geom%rot(2, 1)*geom%rot(3, 2) - geom%rot(2, 2)*geom%rot(3, 1)) &
                               - 1.0_dp) < 1.0e-12_dp, trim(MOL(im))//": det(rot) = +1")
         if (allocated(error)) return
         select case (FKIND(im))
         case (1)
            do k = 1, 3
               call check(error, abs(abs(dot_product(geom%rot(k, :), ax(k, :))) - 1.0_dp) < 1.0e-6_dp, &
                          trim(MOL(im))//": axis "//to_string(k)//" is pyscf's")
               if (allocated(error)) return
            end do
         case (2)
            call check(error, abs(abs(dot_product(geom%rot(3, :), ax(3, :))) - 1.0_dp) < 1.0e-6_dp, &
                       trim(MOL(im))//": z is the molecular axis, as pyscf's")
            if (allocated(error)) return
            call check(error, maxval(abs(geom%r_std(1:2, :))) < 1.0e-8_dp, trim(MOL(im))//": atoms on the z axis")
            if (allocated(error)) return
         case default
            call check(error, frames_equivalent(r, elem, geom%rot, geom%shift, ax, orig), &
                       trim(MOL(im))//": frame is pyscf's up to a symmetry operation")
            if (allocated(error)) return
         end select
         ! the image of an atom under operation k is the mask-signed coordinate of its partner
         dev = 0.0_dp
         do k = 2, geom%group%order
            do i = 1, nat
               j = geom%perm(i, k)
               do ic = 1, 3
                  sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
                  dev = max(dev, abs(geom%r_std(ic, j) - sg*geom%r_std(ic, i)))
               end do
            end do
         end do
         call check(error, dev < 1.0e-14_dp, trim(MOL(im))//": images are the signed coordinates, deviation "// &
                    to_string(dev, "es10.2"))
         if (allocated(error)) return
         call check(error, geom%residual < 1.0e-8_dp, trim(MOL(im))//": the input geometry is symmetric")
         if (allocated(error)) return
      end do
   end subroutine test_detection


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

   subroutine compare_pair(error, p, label, tol_e, density)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: tol_e
      logical, intent(in) :: density  !! false: which degenerate orbital is occupied is arbitrary
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
      if (density) call check(error, dd <= TOL_D, label//": densities differ by "//to_string(dd, "es10.2"))
   end subroutine compare_pair

   subroutine check_equivariance(error, geom, g, label, asym, forbidden)
      !! g(image) = s g(atom) for every operation, s the sign of the component;
      !! on an atom that is its own image a component with s = -1 is zero.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t), intent(in) :: geom
      real(dp), intent(in) :: g(:, :)
      character(len=*), intent(in) :: label
      real(dp), intent(out) :: asym, forbidden
      integer :: a, k, ic
      real(dp) :: sg
      asym = 0.0_dp
      forbidden = 0.0_dp
      do a = 1, size(g, 2)
         do k = 2, geom%group%order
            do ic = 1, 3
               sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
               asym = max(asym, abs(g(ic, geom%perm(a, k)) - sg*g(ic, a)))
               if (geom%perm(a, k) == a .and. sg < 0.0_dp) forbidden = max(forbidden, abs(g(ic, a)))
            end do
         end do
      end do
      call check(error, asym <= 1.0e-14_dp .and. forbidden <= 1.0e-14_dp, label//": equivariance "// &
                 to_string(asym, "es10.2")//", forbidden "//to_string(forbidden, "es10.2"))
   end subroutine check_equivariance

   subroutine noise_case(error, xyz)
      !! The projection on a gradient that is NOT symmetric: the result is
      !! equivariant, idempotent, and an already symmetric gradient is kept.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :), g(:, :), gs(:, :), g2(:, :)
      real(dp) :: asym, forbidden, idem, moved
      integer :: a, ic, nat
      call read_xyz(xyz, nat, elem, r)
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
      moved = maxval(abs(gs - g))
      call check(error, moved > 0.1_dp, xyz//": the projection changed a non-symmetric gradient")
      if (allocated(error)) return
      g2 = gs
      call sym_scf_grad_symmetrize(geom%group, geom%perm, g2)
      idem = maxval(abs(g2 - gs))
      call check_equivariance(error, geom, gs, xyz//" projection", asym, forbidden)
      if (allocated(error)) return
      call check(error, idem <= 1.0e-14_dp, xyz//": projection is idempotent, "//to_string(idem, "es10.2"))
   end subroutine noise_case

   subroutine analyse(bas, geom, label, error, so_dims)
      !! U^T U = I, block-diagonal S and H, S invariant under the signed permutations;
      !! returns the SO dimensions.
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
      real(dp) :: worst_s, worst_h, worst_u, worst_r, tol

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
      tol = 1.0e-12_dp
      call check(error, worst_s <= tol, label//": U^T S U off-block elements "//to_string(worst_s, "es10.2"))
      if (allocated(error)) return
      call check(error, worst_h <= tol, label//": U^T H U off-block elements "//to_string(worst_h, "es10.2"))
      if (allocated(error)) return
      call check(error, worst_r <= tol, label//": S invariant under the signed permutations")
      if (allocated(error)) return
      so_dims = ao%irrep_dim
   end subroutine analyse

   subroutine check_dims(u, im, tagname, geom, elem, label, error)
      !! SALC identities and the SO dimensions against the anchor `tagname`, both bases.
      integer, intent(in) :: u, im
      character(len=*), intent(in) :: tagname, label
      type(trc_sym_geom_t), intent(in) :: geom
      integer, intent(in) :: elem(:)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      integer, allocatable :: dims(:)
      real(dp), allocatable :: ref(:)
      integer :: ib, nao, nao_ref, nirr
      character(len=64) :: lab
      do ib = 1, NBAS
         lab = trim(label)//' '//trim(BAS_NAME(ib))//' '//KIND_NAME(ib)
         call make_basis(BAS_FILE(ib), elem, geom, ib == 2, bas, error)
         if (allocated(error)) return
         call analyse(bas, geom, trim(lab), error, dims)
         nao = merge(bas%nao_sph, bas%nao, bas%spherical)
         call bas%release()
         if (allocated(error)) return
         call goto_tag(u, 'NAO_'//trim(tagname)//'_'//trim(BAS_NAME(ib)))
         call ref_read_int(u, 'NAO_'//trim(tagname)//'_'//trim(BAS_NAME(ib)), nao_ref)
         nirr = size(dims)
         allocate (ref(nirr))
         call ref_read_vector(u, 'SO_'//trim(tagname)//'_'//trim(BAS_NAME(ib)), nirr, ref)
         call check(error, nao, nao_ref, trim(lab)//": AO count against pyscf")
         if (.not. allocated(error)) call check(error, all(dims == nint(ref)), &
                                                trim(lab)//": SO dimensions against pyscf")
         deallocate (ref)
         if (allocated(error)) return
      end do
      if (im < 0) return
   end subroutine check_dims

   subroutine test_salc(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      integer :: im, u, nat
      call ref_require(REF, u)
      do im = 1, NMOL
         call detect_mol(MOL(im), nat, elem, r, geom, error)
         if (allocated(error)) exit
         call check_dims(u, im, trim(DIMREF(im)), geom, elem, trim(MOL(im)), error)
         if (allocated(error)) exit
      end do
      call ref_close(u)
   end subroutine test_salc

   subroutine make_pair(error, xyz, basis_json, spherical, functional, charge, mult, group, requested, run_c1, p)
      !! Detect the group, build the basis on the standard frame, and run the SCF
      !! with symmetry (and as C1 when `run_c1`).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, functional, group
      logical, intent(in) :: spherical, run_c1
      integer, intent(in) :: charge, mult
      character(len=*), intent(in) :: requested  !! "" for the full group
      type(pair_t), intent(inout) :: p
      integer :: nat, nelec, dev
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_sym_scf_t) :: sc

      dev = trc_bind_device(0)
      call read_xyz(xyz//'.xyz', nat, elem, r)
      if (len_trim(requested) > 0) then
         p%geom = sym_geom_detect(r, elem, err, basis_id=elem, requested=requested)
      else
         p%geom = sym_geom_detect(r, elem, err, basis_id=elem)
      end if
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, p%geom%group%name, group, "group of "//xyz)
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
      opts%grid_level = DFT_GRID_LEVEL
      ! Just above the OpenMP/device accumulation-order noise floor (see check_sym_scf_cs).
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100

      sc%ao = p%ao
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%rs, sym=sc)
      call check(error, p%rs%converged, xyz//": the symmetric SCF did not converge: "//trim(p%rs%message))
      if (allocated(error)) return
      call sc%release()
      call check(error, allocated(p%rs%mo_irrep), xyz//": the symmetric run carries irrep labels")
      if (allocated(error)) return
      if (.not. run_c1) return
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%r1)
      call check(error, p%r1%converged, xyz//": the C1 SCF did not converge: "//trim(p%r1%message))
   end subroutine make_pair


   subroutine find_case(u, tag, nspin, e)
      !! Position the reader on the NSPIN record of `tag`.
      integer, intent(in) :: u
      character(len=*), intent(in) :: tag
      integer, intent(out) :: nspin
      real(dp), intent(out) :: e
      character(len=512) :: line
      integer :: ios
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            print '(a)', 'check_sym_tops: no case '//tag//' in '//REF
            stop 1
         end if
         if (index(line, 'NSPIN_'//tag//' ') == 1) exit
      end do
      read (line(len('NSPIN_'//tag) + 1:), *) nspin
      read (u, '(a)') line
      read (line(len('E_'//tag) + 1:), *) e
   end subroutine find_case

   subroutine compare_anchor(error, p, tag, tol, tol_e, by_irrep)
      !! Electrons per irrep and orbital energies per irrep against pyscf (with
      !! `by_irrep` false: only the sorted orbital energies of each spin; the energy is
      !! compared in both cases).
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: tag
      real(dp), intent(in) :: tol, tol_e
      logical, intent(in) :: by_irrep
      integer :: u, nspin, ir, s, m, nocc(2), nref, k, nirr, ioff
      integer, allocatable :: cnt(:, :), want(:, :)
      real(dp) :: eref_e, worst
      real(dp), allocatable :: ev(:), eir(:), nelec(:), eall(:), tall(:)

      nirr = p%ao%nirrep
      allocate (cnt(nirr, 2), want(nirr, 2), nelec(2*nirr))
      call ref_require(REF, u)
      call goto_tag(u, 'NSPIN_'//tag)
      call find_case(u, tag, nspin, eref_e)
      call ref_read_vector(u, 'NELEC_'//tag, 2*nirr, nelec)
      nocc = [p%nalpha, p%nbeta]
      cnt = 0
      do s = 1, p%rs%nspin
         do m = 1, nocc(s)
            ir = p%rs%mo_irrep(m, s)
            cnt(ir, s) = cnt(ir, s) + 1
         end do
      end do
      if (p%rs%nspin == 1) cnt(:, 2) = cnt(:, 1)
      do s = 1, 2
         ioff = (s - 1)*nirr
         do ir = 1, nirr
            want(ir, s) = nint(nelec(ioff + ir))
         end do
      end do
      print '(a,a,a,es10.2)', '  ', tag, ': |E - E(pyscf)| ', abs(p%rs%energy - eref_e)
      call check(error, abs(p%rs%energy - eref_e) <= tol_e, tag//": energy differs from pyscf by "// &
                 to_string(abs(p%rs%energy - eref_e), "es10.2"))
      if (allocated(error)) then
         call ref_close(u)
         return
      end if
      if (by_irrep) then
         call check(error, all(cnt == want), tag//": electrons per irrep differ from pyscf")
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
      end if
      worst = 0.0_dp
      do s = 1, nspin
         allocate (eall(0))
         do ir = 1, nirr
            nref = p%ao%irrep_dim(ir)
            allocate (eir(nref))
            call goto_tag(u, 'EPS'//to_string(ir)//'_'//tag//'_'//to_string(s))
            call ref_read_vector(u, 'EPS'//to_string(ir)//'_'//tag//'_'//to_string(s), nref, eir)
            if (by_irrep) then
               ev = pack(p%rs%eps(1:p%ao%nao, s), p%rs%mo_irrep(1:p%ao%nao, s) == ir)
               call check(error, size(ev) == nref, tag//": orbitals in irrep "//to_string(ir)//" differ from pyscf")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               do k = 1, nref
                  worst = max(worst, abs(ev(k) - eir(k)))
               end do
            else
               eall = [eall, eir]
            end if
            deallocate (eir)
         end do
         if (.not. by_irrep) then
            tall = p%rs%eps(1:p%ao%nao, s)
            call sort(eall)
            call sort(tall)
            worst = max(worst, maxval(abs(eall - tall)))
         end if
         deallocate (eall)
      end do
      call ref_close(u)
      print '(a,a,a,es10.2)', '  ', tag, ': max |eps - eps(pyscf)| ', worst
      call check(error, worst <= tol, tag//": orbital energies differ from pyscf by "//to_string(worst, "es10.2"))
   end subroutine compare_anchor

   subroutine run_case(error, tag, xyz, ib, functional, charge, mult, group, by_irrep)
      !! The central oracle and the anchors for one molecule.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: tag, xyz, functional, group
      integer, intent(in) :: ib, charge, mult
      logical, intent(in) :: by_irrep
      type(pair_t) :: p
      real(dp) :: tol_e, tol_a
      logical :: dft
      dft = len_trim(functional) > 0
      tol_e = merge(TOL_E_DFT, TOL_E_HF, dft)
      tol_a = merge(TOL_ANCHOR_DFT, TOL_ANCHOR_HF, dft)
      call make_pair(error, xyz, BAS_FILE(ib), ib == 2, functional, charge, mult, group, '', .true., p)
      if (allocated(error)) return
      call compare_pair(error, p, tag, tol_e, by_irrep)
      call golden_put('E_'//tag, p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, tag, tol_a, merge(1.0e-6_dp, 1.0e-8_dp, dft), by_irrep)
      call p%bas%release()
   end subroutine run_case

   subroutine test_benzene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'benzene_rhf', 'tops_benzene_rot', 1, '', 0, 1, 'D2h', .true.)
   end subroutine test_benzene

   subroutine test_nh3_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'nh3_rhf', 'tops_nh3_rot', 2, '', 0, 1, 'Cs', .true.)
   end subroutine test_nh3_rhf

   subroutine test_nh3_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'nh3_uhf', 'tops_nh3', 2, '', 1, 2, 'Cs', .true.)
   end subroutine test_nh3_uhf

   subroutine test_nh3_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'nh3_b3lyp', 'tops_nh3_rot', 2, 'b3lyp', 0, 1, 'Cs', .true.)
   end subroutine test_nh3_b3lyp

   subroutine test_ch4(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'ch4_rhf', 'tops_ch4_rot', 2, '', 0, 1, 'D2', .true.)
   end subroutine test_ch4

   subroutine test_sf6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'sf6_rhf', 'tops_sf6_rot', 1, '', 0, 1, 'D2h', .true.)
   end subroutine test_sf6

   subroutine test_n2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'n2_rhf', 'tops_n2_rot', 2, '', 0, 1, 'D2h', .true.)
   end subroutine test_n2

   subroutine test_co(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'co_rhf', 'tops_co_rot', 2, '', 0, 1, 'C2v', .true.)
   end subroutine test_co

   subroutine test_hcn(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'hcn_rhf', 'tops_hcn_rot', 2, '', 0, 1, 'C2v', .true.)
   end subroutine test_hcn

   subroutine test_o2(error)
      !! The triplet: both codes must land in the same state (energy, and one alpha
      !! electron in each of the degenerate pi* pair, B2g and B3g).
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'o2_uhf', 'tops_o2_rot', 2, '', 0, 3, 'D2h', .true.)
   end subroutine test_o2

   subroutine test_ne(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'ne_rhf', 'tops_ne', 2, '', 0, 1, 'D2h', .true.)
   end subroutine test_ne

   subroutine test_o_atom(error)
      !! 3P oxygen: which degenerate p orbital takes the beta electron is arbitrary in
      !! both codes, so only the energy and the sorted orbital energies are anchored.
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 'o_uhf', 'o_atom', 2, '', 0, 3, 'D2h', .false.)
   end subroutine test_o_atom

   subroutine test_benzene_sub(error)
      !! symmetry='C2v' and 'Cs' on benzene: the frame pyscf uses (input axes when the
      !! group holds there, else the D2h frame, which `as_subgroup` leaves unchanged), the
      !! SO dimensions and the energy.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      real(dp) :: ax(3, 3), orig(3)
      integer, parameter :: NS = 2
      character(len=3), parameter :: SUBS(NS) = [character(len=3) :: "C2v", "Cs"]
      character(len=16), parameter :: MB(NS) = [character(len=16) :: "tops_benzene", "tops_benzene_rot"]
      integer :: ib, is, k, nat, u
      type(pair_t) :: p
      do ib = 1, NS
         do is = 1, NS
            call detect_mol(MB(ib), nat, elem, r, geom, error, requested=SUBS(is))
            if (allocated(error)) return
            call check(error, geom%detected, "D2h", "benzene detected group")
            if (allocated(error)) return
            call check(error, geom%group%name, SUBS(is), "benzene as "//SUBS(is))
            if (allocated(error)) return
            call pyscf_frame(trim(MB(ib))//'_'//trim(SUBS(is)), ax, orig)
            call check(error, maxval(abs(geom%shift - orig)) < 1.0e-6_dp, "origin of "//trim(MB(ib))//' '//SUBS(is))
            if (allocated(error)) return
            do k = 1, 3
               call check(error, abs(abs(dot_product(geom%rot(k, :), ax(k, :))) - 1.0_dp) < 1.0e-6_dp, &
                          trim(MB(ib))//' '//trim(SUBS(is))//": axis "//to_string(k)//" is pyscf's")
               if (allocated(error)) return
            end do
            call ref_require(REF, u)
            call check_dims(u, 0, trim(MB(ib))//'_'//trim(SUBS(is)), geom, elem, trim(MB(ib))//' '//SUBS(is), error)
            call ref_close(u)
            if (allocated(error)) return
         end do
      end do
      do is = 1, NS
         call make_pair(error, 'tops_benzene_rot', BAS_FILE(1), .false., '', 0, 1, SUBS(is), SUBS(is), .false., p)
         if (allocated(error)) return
         call golden_put('E_benzene_'//to_lower(trim(SUBS(is))), p%rs%energy)
         call compare_anchor(error, p, 'benzene_'//to_lower(trim(SUBS(is)))//'_rhf', TOL_ANCHOR_HF, 1.0e-8_dp, .true.)
         call p%bas%release()
         if (allocated(error)) return
      end do
   end subroutine test_benzene_sub

   subroutine grad_case(error, xyz, basis_json, spherical, group)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, group
      logical, intent(in) :: spherical
      type(pair_t) :: p
      type(trc_linalg_t) :: lp
      real(dp), allocatable :: gsym(:, :), g1(:, :)
      real(dp) :: dg, asym, forbidden

      call make_pair(error, xyz, basis_json, spherical, '', 0, 1, group, '', .true., p)
      if (allocated(error)) return
      allocate (gsym(3, p%bas%natm), g1(3, p%bas%natm))
      call lp%init(p%bas%nao)
      call trc_scf_gradient(p%bas, lp, p%rs%nspin, p%nalpha, p%nbeta, p%rs%cmo, p%rs%eps, 1.0e-12_dp, gsym, &
                            sym=trc_scf_grad_sym_t(geometry=p%geom))
      call trc_scf_gradient(p%bas, lp, p%r1%nspin, p%nalpha, p%nbeta, p%r1%cmo, p%r1%eps, 1.0e-12_dp, g1)
      call lp%release()
      call p%bas%release()
      dg = maxval(abs(gsym - g1))
      print '(a,a,a,es10.2)', '  ', xyz, ': |sym gradient - C1 gradient| ', dg
      call check_equivariance(error, p%geom, gsym, xyz//" symmetrised gradient", asym, forbidden)
      if (allocated(error)) return
      call check(error, dg <= 1.0e-9_dp, xyz//": symmetrised gradient differs from C1 by "//to_string(dg, "es10.2"))
   end subroutine grad_case

   subroutine test_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: im
      call grad_case(error, 'tops_nh3_rot', 'basis_sets/6-31gs.json', .false., 'Cs')
      if (allocated(error)) return
      call grad_case(error, 'tops_n2_rot', 'basis_sets/6-31gs.json', .false., 'D2h')
      if (allocated(error)) return
      do im = 1, NMOL
         call noise_case(error, trim(MOL(im))//'.xyz')
         if (allocated(error)) return
      end do
   end subroutine test_gradient

   subroutine test_choose(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_sym_request_t) :: req
      type(error_t) :: e1, e2, e3
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      character(len=3) :: g
      integer :: nat, im
      do im = 1, NMOL
         call detect_mol(MOL(im), nat, elem, r, geom, error)
         if (allocated(error)) return
         req%mode = SYM_AUTO
         g = sym_run_choose_group(req, geom, e1)
         call check(error, g == GROUP(im) .and. .not. e1%has_error(), 'auto picks '//GROUP(im)//' for '//trim(MOL(im)))
         if (allocated(error)) return
         req%mode = SYM_EXPLICIT
         req%group = GROUP(im)
         g = sym_run_choose_group(req, geom, e2)
         call check(error, g == GROUP(im) .and. .not. e2%has_error(), 'explicit '//GROUP(im)//' is accepted for '// &
                    trim(MOL(im)))
         if (allocated(error)) return
      end do
      ! a subgroup of a linear molecule or atom has no frame rule: refused rather than guessed
      call detect_mol('tops_n2', nat, elem, r, geom, error, requested='C2v')
      if (allocated(error)) return
      req%group = 'C2v'
      g = sym_run_choose_group(req, geom, e3)
      call check(error, e3%has_error(), 'C2v requested on N2 is refused: '//e3%get_message())
   end subroutine test_choose

end module check_sym_tops_suite

program check_sym_tops
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use golden_io, only: golden_begin, golden_end
   use check_sym_tops_suite, only: collect_sym_tops
   implicit none
   integer :: stat

   stat = 0
   call golden_begin('sym_tops')
   call run_testsuite(collect_sym_tops, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
   call golden_end()
end program check_sym_tops
