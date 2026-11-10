!
! Point group D2h (trc_sym_geom, trc_sym_ao, trc_sym_scf), the last Abelian
! group. Molecules, with the source of each geometry in its .xyz comment line:
!   c2h4_d2h          planar ethylene, general orientation, shifted
!   c2h4_d2h_aligned  planar ethylene, axis-aligned (xy plane, C=C along x)
!   naphthalene_d2h   planar naphthalene, axis-aligned (xy plane)
!   b2h6_d2h          diborane, NONPLANAR (the bridging H atoms leave the BH2
!                     plane), axis-aligned
!   b2h6_d2h_rot      diborane, general orientation, shifted
! Before the group was enabled a D2h molecule was capped to Cs with a
! provisional principal-axis order that had never been checked against pyscf.
!
!   - detection: group D2h, the frame and origin as pyscf 2.14.0 gives them
!     (`_symm_axes`, `_symm_orig`, read from sym_d2h.ref). pyscf names the D2h
!     axes from the principal axes in ascending order of the second moment (no
!     `alias_axes`, unlike D2), so the frame does not depend on LAPACK's
!     eigenvector signs and no x/y exchange of the D2 kind can happen: the rows
!     must equal pyscf's up to sign, for every molecule. For a planar molecule
!     pyscf's `_adjust_planar_d2h` makes x the plane normal; checked here;
!   - nonplanar B2H6: no mirror leaves every atom fixed, so the operations
!     move atoms off every plane; planar molecules have exactly one mirror
!     that fixes all atoms;
!   - SALC: U^T U = I, U^T S U and U^T H U block-diagonal (<= 1e-12), S
!     invariant under the signed permutations, the dimensions per irrep (Ag,
!     B1g, B2g, B3g, Au, B1u, B2u, B3u) equal to pyscf `mol.symm_orb`; 6-31G*
!     Cartesian and cc-pVDZ spherical; and (mode "f", `check_sym_d2h_f`, which
!     prints SKIP when this build LMAX is below 3) the uncontracted l = 0..3 set;
!   - THE CENTRAL ORACLE: symmetry vs C1 on the symmetrised geometry, RHF, UHF
!     (the ethylene cation) and B3LYP at grid level 5, |dE| <= 1e-10 (HF) /
!     1e-9 (DFT), orbital energies, density;
!   - anchors (sym_d2h.ref): electrons per irrep and orbital energies per irrep;
!   - gradient: symmetrised gradient equivariant and equal to C1, and the
!     projection itself on an arbitrary (non-symmetric) gradient;
!   - `symmetry = auto` and `symmetry = D2h` choose D2h; D2 is refused
!     (pyscf has no D2 subgroup of D2h either).
! The requested subgroups are in check_sym_d2h_subgroup.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_d2h_suite
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
   use trc_sym_run, only: trc_sym_request_t, sym_run_choose_group, SYM_AUTO, SYM_EXPLICIT
   use trc_sym_scf, only: trc_sym_scf_t, sym_scf_grad_symmetrize
   use trc_test_basis, only: read_xyz, build_aux
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   use golden_io, only: golden_put
   implicit none
   private
   public :: collect_sym_d2h, g_mode

   character(len=8) :: g_mode = 'main'
      !! "main" or "f", set from argv before run_testsuite.

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the brief's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the brief's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies, symmetry vs C1
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density, symmetry vs C1
   real(dp), parameter :: TOL_ANCHOR_HF = 1.0e-6_dp   !! terco vs pyscf, HF orbital energies
   real(dp), parameter :: TOL_ANCHOR_DFT = 1.0e-5_dp
      !! B3LYP orbital energies against pyscf: the two codes use different grids; the gate
      !! rests on the grid study in check_sym_ci (terco grid level 5 against pyscf level 6).
   integer, parameter :: DFT_GRID_LEVEL = 5  !! terco grid level of the B3LYP case (see the plan, "Known issues")

   integer, parameter :: NMOL = 5, NBAS = 2, NIR = 8
   character(len=20), parameter :: MOL(NMOL) = [character(len=20) :: "c2h4_d2h", "c2h4_d2h_aligned", "naphthalene_d2h", &
                                                                     "b2h6_d2h", "b2h6_d2h_rot"]
   logical, parameter :: PLANAR(NMOL) = [.true., .true., .true., .false., .false.]
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

   subroutine collect_sym_d2h(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == 'f') then
         testsuite = [new_unittest("uncontracted l = 0..3: parity, U^T U, block-diagonal S and H", test_lmax3)]
      else
         testsuite = [ &
                     new_unittest("detection: D2h, pyscf frame (rotated, shifted and aligned input)", test_detection), &
                     new_unittest("SO basis vs pyscf dimensions, block-diagonal S and H", test_salc), &
                     new_unittest("RHF C2H4 cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_c2h4_rhf), &
                     new_unittest("RHF C2H4 6-31G* (Cartesian, aligned): symmetry vs C1 and pyscf", test_c2h4a_rhf), &
                     new_unittest("RHF naphthalene 6-31G* (Cartesian): symmetry vs C1 and pyscf", test_nap_rhf), &
                     new_unittest("RHF B2H6 cc-pVDZ (spherical, rotated): symmetry vs C1 and pyscf", test_b2h6_rhf), &
                     new_unittest("RHF B2H6 6-31G* (Cartesian, aligned): symmetry vs C1 and pyscf", test_b2h6a_rhf), &
                     new_unittest("UHF C2H4 cation cc-pVDZ: symmetry vs C1 and pyscf", test_c2h4_uhf), &
                     new_unittest("B3LYP C2H4 cc-pVDZ: symmetry vs C1 and pyscf", test_c2h4_b3lyp), &
                     new_unittest("gradient: equivariant, equals C1, projection", test_gradient), &
                     new_unittest("auto and explicit D2h choose D2h, D2 is refused", test_choose) &
                     ]
      end if
   end subroutine collect_sym_d2h

   subroutine detect_mol(name, elem, geom, error)
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
      call check(error, geom%detected, "D2h", "detected group of "//trim(name))
      if (allocated(error)) return
      call check(error, geom%group%name, "D2h", "group of "//trim(name))
   end subroutine detect_mol

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
            print '(a)', 'check_sym_d2h: no record '//tag//' in sym_d2h.ref'
            stop 1
         end if
         if (index(line, tag//' ') == 1) exit
      end do
      backspace (u)
   end subroutine goto_tag

   subroutine pyscf_frame(mol, ax, orig)
      !! pyscf's frame of `mol` (rows x, y, z) and origin from the reference.
      character(len=*), intent(in) :: mol
      real(dp), intent(out) :: ax(3, 3)
      real(dp), intent(out), optional :: orig(3)
      real(dp) :: v(9)
      integer :: u
      call ref_require('sym_d2h.ref', u)
      call goto_tag(u, 'AXES_'//trim(mol))
      call ref_read_vector(u, 'AXES_'//trim(mol), 9, v)
      ax = transpose(reshape(v, [3, 3]))
      if (present(orig)) then
         call goto_tag(u, 'ORIG_'//trim(mol))
         call ref_read_vector(u, 'ORIG_'//trim(mol), 3, orig)
      end if
      call ref_close(u)
   end subroutine pyscf_frame

   subroutine test_detection(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp) :: ax(3, 3), orig(3), dev, sg
      integer :: im, i, j, k, nat, ic, ntot_fixed_mirror, nfix
      logical :: moved
      do im = 1, NMOL
         call detect_mol(MOL(im), elem, geom, error)
         if (allocated(error)) return
         call check(error, geom%group%order, 8, "eight operations")
         if (allocated(error)) return
         call check(error, all(geom%group%ops == [0, 6, 5, 3, 7, 1, 2, 4]), "operations E, C2x, C2y, C2z, i, sx, sy, sz")
         if (allocated(error)) return
         call pyscf_frame(MOL(im), ax, orig)
         do k = 1, 3
            call check(error, abs(abs(dot_product(geom%rot(k, :), ax(k, :))) - 1.0_dp) < 1.0e-6_dp, &
                       trim(MOL(im))//": axis "//to_string(k)//" is pyscf's")
            if (allocated(error)) return
         end do
         call check(error, maxval(abs(geom%shift - orig)) < 1.0e-6_dp, "origin of "//trim(MOL(im)))
         if (allocated(error)) return
         nat = size(elem)
         if (MOL(im) == 'c2h4_d2h_aligned' .or. MOL(im) == 'naphthalene_d2h' .or. MOL(im) == 'b2h6_d2h') then
            ! axes aligned with the input: the frame is a signed permutation of the input axes
            call check(error, all(abs(abs(geom%rot) - 1.0_dp) < 1.0e-6_dp .or. abs(geom%rot) < 1.0e-6_dp) .and. &
                       all(count(abs(geom%rot) > 0.5_dp, dim=1) == 1), trim(MOL(im))//": aligned input, permuted axes")
            if (allocated(error)) return
         end if
         ! the image of an atom under operation k is the mask-signed coordinate of its partner
         dev = 0.0_dp
         ntot_fixed_mirror = 0
         do k = 2, 8
            moved = .false.
            nfix = 0
            do i = 1, nat
               j = geom%perm(i, k)
               if (j /= i) moved = .true.
               if (j == i) nfix = nfix + 1
               do ic = 1, 3
                  sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
                  dev = max(dev, abs(geom%r_std(ic, j) - sg*geom%r_std(ic, i)))
               end do
            end do
            if (k >= 6 .and. nfix == nat) ntot_fixed_mirror = ntot_fixed_mirror + 1
            if (k == 5) call check(error, nfix == 0, trim(MOL(im))//": the inversion fixes an atom")
            if (allocated(error)) return
            call check(error, moved .or. (k >= 6 .and. PLANAR(im)), trim(MOL(im))//": operation "//to_string(k)//" moves atoms")
            if (allocated(error)) return
         end do
         ! planar: exactly one mirror fixes every atom; nonplanar: none does
         call check(error, ntot_fixed_mirror == merge(1, 0, PLANAR(im)), trim(MOL(im))//": mirrors fixing every atom: "// &
                    to_string(ntot_fixed_mirror))
         if (allocated(error)) return
         if (PLANAR(im)) then
            ! pyscf _adjust_planar_d2h: x is the normal of the plane
            call check(error, maxval(abs(geom%r_std(1, :))) < 1.0e-8_dp, trim(MOL(im))//": x is the plane normal")
            if (allocated(error)) return
         end if
         call check(error, dev < 1.0e-14_dp, trim(MOL(im))//": images are the signed coordinates, deviation "// &
                    to_string(dev, "es10.2"))
         if (allocated(error)) return
         call check(error, geom%residual < 1.0e-5_dp, trim(MOL(im))//": the input geometry is symmetric to 1e-5 Bohr")
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

   subroutine check_inversion(bas, ao, label, error)
      !! The inversion sends every AO to an AO on another atom, with the sign
      !! (-1)^l of its shell (op 5 of D2h is the inversion).
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
            if (ao%ao_sign(a, 5) /= (-1)**l) nbad = nbad + 1
            if (ao%ao_map(a, 5) == a) nfix = nfix + 1
         end do
      end do
      call check(error, nbad == 0, label//": "//to_string(nbad)//" AOs have a sign other than (-1)^l under the inversion")
      if (allocated(error)) return
      call check(error, nfix == 0, label//": the inversion fixes an AO")
   end subroutine check_inversion

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

   subroutine test_salc(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), dims(:)
      real(dp) :: ref(NIR)
      integer :: im, ib, ik, u, nao_ref, nao
      character(len=64) :: label

      call ref_require('sym_d2h.ref', u)
      do im = 1, NMOL
         call detect_mol(MOL(im), elem, geom, error)
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
         do ib = 1, NBAS
            do ik = 1, 2
               if (ik /= ib) cycle  ! 6-31G* Cartesian, cc-pVDZ spherical: the anchored pairs
               label = trim(MOL(im))//' '//trim(BAS_NAME(ib))//' '//KIND_NAME(ik)
               call make_basis(BAS_FILE(ib), elem, geom, ik == 2, bas, error)
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call analyse(bas, geom, trim(label), .false., error, dims)
               nao = merge(bas%nao_sph, bas%nao, bas%spherical)
               call bas%release()
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call goto_tag(u, 'NAO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib)))
               call ref_read_int(u, 'NAO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib)), nao_ref)
               call ref_read_vector(u, 'SO_'//trim(MOL(im))//'_'//trim(BAS_NAME(ib)), NIR, ref)
               call check(error, nao, nao_ref, trim(label)//": AO count against pyscf")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call check(error, size(dims), NIR, trim(label)//": number of irreps")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call check(error, all(dims == nint(ref)), trim(label)//": SO dimensions (Ag, B1g, ... B3u) against pyscf")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
            end do
         end do
      end do
      call ref_close(u)
   end subroutine test_salc

   subroutine test_lmax3(error)
      !! Uncontracted l = 0..3 shells on every atom, Cartesian and spherical.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: LMAX3_MOL(3) = [1, 4, 5]  !! ethylene, B2H6 aligned and rotated
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
      integer :: im, ii, ik, nsh, maxnp
      do ii = 1, size(LMAX3_MOL)
         im = LMAX3_MOL(ii)
         call detect_mol(MOL(im), elem, geom, error)
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

   subroutine make_pair(error, xyz, basis_json, spherical, functional, charge, mult, p, conv_diis)
      !! Detect D2h, build the basis on the standard frame, and run the SCF
      !! with symmetry and as C1. `conv_diis` tightens the commutator
      !! convergence of both runs when given (default 1e-9).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, functional
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t), intent(inout) :: p
      real(dp), intent(in), optional :: conv_diis
      integer :: nat, nelec, dev
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_sym_scf_t) :: sc

      dev = trc_bind_device(0)
      call read_xyz(xyz, nat, elem, r)
      p%geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, p%geom%group%name, "D2h", "group of "//xyz)
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
      ! Just above the OpenMP/device accumulation-order noise floor (see
      ! check_sym_scf_cs); the energy is variational, so the 1e-10 gates hold.
      opts%conv_energy = 1.0e-10_dp
      opts%conv_diis = 1.0e-9_dp
      if (present(conv_diis)) opts%conv_diis = conv_diis
      opts%max_iteration = 100

      sc%ao = p%ao
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%rs, sym=sc)
      call check(error, p%rs%converged, xyz//": the symmetric SCF did not converge: "//trim(p%rs%message))
      if (allocated(error)) return
      call sc%release()
      call check(error, allocated(p%rs%mo_irrep), xyz//": the symmetric run carries irrep labels")
      if (allocated(error)) return
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, p%r1)
      call check(error, p%r1%converged, xyz//": the C1 SCF did not converge: "//trim(p%r1%message))
   end subroutine make_pair

   subroutine compare_pair(error, p, label, tol_e)
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
            print '(a)', 'check_sym_d2h: no case '//tag//' in sym_d2h.ref'
            stop 1
         end if
         if (index(line, 'NSPIN_'//tag//' ') == 1) exit
      end do
      read (line(len('NSPIN_'//tag) + 1:), *) nspin
      read (u, '(a)') line
      read (line(len('E_'//tag) + 1:), *) e
   end subroutine find_case

   subroutine compare_anchor(error, p, tag, tol)
      !! Electrons per irrep and orbital energies per irrep against pyscf
      !! .
      type(error_type), allocatable, intent(out) :: error
      type(pair_t), intent(in) :: p
      character(len=*), intent(in) :: tag
      real(dp), intent(in) :: tol
      integer :: u, nspin, ir, s, m, nocc(2), nref, k, nirr, ioff
      integer, allocatable :: cnt(:, :), want(:, :)
      real(dp) :: eref_e, worst
      real(dp), allocatable :: ev(:), eir(:), nelec(:)

      nirr = p%ao%nirrep
      allocate (cnt(nirr, 2), want(nirr, 2), nelec(2*nirr))
      call ref_require('sym_d2h.ref', u)
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
      call check(error, all(cnt == want), tag//": electrons per irrep differ from pyscf: "// &
                 to_string(cnt(1, 1))//","//to_string(cnt(2, 1))//" / "//to_string(want(1, 1))//","//to_string(want(2, 1)))
      if (allocated(error)) then
         call ref_close(u)
         return
      end if
      worst = 0.0_dp
      do s = 1, nspin
         do ir = 1, nirr
            nref = p%ao%irrep_dim(ir)
            allocate (eir(nref))
            call goto_tag(u, 'EPS'//to_string(ir)//'_'//tag//'_'//to_string(s))
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

   subroutine test_c2h4_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2h4_d2h.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF C2H4 cc-pVDZ', TOL_E_HF)
      call golden_put('E_c2h4_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_c2h4_rhf

   subroutine test_c2h4a_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2h4_d2h_aligned.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF C2H4 6-31G* aligned', TOL_E_HF)
      call golden_put('E_c2h4a_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4a_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_c2h4a_rhf

   subroutine test_nap_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      ! The full Fock build's noise on 18 atoms leaves the orbital energies
      ! 1.2e-8 apart at the 1e-9 commutator gate (51 iterations against C1's
      ! 38), over the 1e-8 TOL_EPS; the gates stay, the SCF converges tighter.
      call make_pair(error, 'naphthalene_d2h.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, conv_diis=1.0e-10_dp)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF naphthalene 6-31G*', TOL_E_HF)
      call golden_put('E_nap_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'nap_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_nap_rhf

   subroutine test_b2h6_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'b2h6_d2h_rot.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF B2H6 cc-pVDZ', TOL_E_HF)
      call golden_put('E_b2h6_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'b2h6_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_b2h6_rhf

   subroutine test_b2h6a_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'b2h6_d2h.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF B2H6 6-31G* aligned', TOL_E_HF)
      call golden_put('E_b2h6a_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'b2h6a_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_b2h6a_rhf

   subroutine test_c2h4_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2h4_d2h.xyz', 'basis_sets/cc-pvdz.json', .true., '', 1, 2, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF C2H4+ cc-pVDZ', TOL_E_HF)
      call golden_put('E_c2h4_uhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4_uhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_c2h4_uhf

   subroutine test_c2h4_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2h4_d2h.xyz', 'basis_sets/cc-pvdz.json', .true., 'b3lyp', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'B3LYP C2H4 cc-pVDZ', TOL_E_DFT)
      call golden_put('E_c2h4_b3lyp', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4_b3lyp', TOL_ANCHOR_DFT)
      call p%bas%release()
   end subroutine test_c2h4_b3lyp

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

   subroutine grad_case(error, xyz, basis_json, spherical, charge, mult)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t) :: p
      type(trc_linalg_t) :: lp
      real(dp), allocatable :: gsym(:, :), g1(:, :)
      real(dp) :: dg, asym, forbidden

      call make_pair(error, xyz, basis_json, spherical, '', charge, mult, p)
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
      call grad_case(error, 'c2h4_d2h.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'b2h6_d2h_rot.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'c2h4_d2h.xyz', 'basis_sets/6-31gs.json', .false., 1, 2)
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
      type(error_t) :: e1, e2, e3, e4
      integer, allocatable :: elem(:)
      character(len=3) :: g
      call detect_mol('c2h4_d2h', elem, geom, error)
      if (allocated(error)) return
      req%mode = SYM_AUTO
      g = sym_run_choose_group(req, geom, e1)
      call check(error, g == 'D2h' .and. .not. e1%has_error(), 'auto picks D2h')
      if (allocated(error)) return
      req%mode = SYM_EXPLICIT
      req%group = 'D2h'
      g = sym_run_choose_group(req, geom, e2)
      call check(error, g == 'D2h' .and. .not. e2%has_error(), 'explicit D2h is accepted')
      if (allocated(error)) return
      ! D2 is a subgroup of D2h as a group, but pyscf refuses symmetry='D2' on D2h and no frame is defined
      req%group = 'D2'
      g = sym_run_choose_group(req, geom, e3)
      call check(error, e3%has_error() .and. index(e3%get_message(), 'no frame is defined') > 0, &
                 'D2 on a D2h molecule is refused: '//e3%get_message())
      if (allocated(error)) return
      ! C2v needs sym_geom_detect(requested='C2v'): not passed here, so the full group stands and the request is refused
      req%group = 'C2v'
      g = sym_run_choose_group(req, geom, e4)
      call check(error, e4%has_error(), 'C2v is refused unless it is passed to sym_geom_detect')
   end subroutine test_choose

end module check_sym_d2h_suite

program check_sym_d2h
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use golden_io, only: golden_begin, golden_end
   use check_sym_d2h_suite, only: collect_sym_d2h, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if (trim(g_mode) == 'f' .and. LMAX < 3) then
      print '(a,i0,a)', '  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ', LMAX, '.'
      stop 0
   end if
   stat = 0
   if (trim(g_mode) /= 'f') call golden_begin('sym_d2h')
   call run_testsuite(collect_sym_d2h, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
   if (trim(g_mode) /= 'f') call golden_end()
end program check_sym_d2h
