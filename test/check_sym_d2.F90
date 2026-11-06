!
! Point group D2 (trc_sym_geom, trc_sym_ao, trc_sym_scf): twisted ethylene
! (general orientation, shifted: the standard frame is not the input frame) and
! twisted tetrafluoroethylene (axes aligned with the input axes). Before the
! group was enabled both were capped to C1 (D2 has no Cs below it).
!
!   - detection: group D2, the frame and origin as pyscf 2.14.0 gives them
!     (`_symm_axes`, `_symm_orig`, read from sym_d2.ref). pyscf swaps x and y
!     when the LAPACK signs of its principal axes are left-handed, a sign that
!     is not defined: the test accepts terco's frame or the one with x and y
!     exchanged, and then maps terco's B2 and B3 to pyscf's B3 and B2 in every
!     anchor below. For the aligned molecule pyscf's axes are the input axes
!     and no swap is possible: the test requires that;
!   - SALC: U^T U = I, U^T S U and U^T H U block-diagonal (<= 1e-12), the
!     dimensions per irrep (A, B1, B2, B3) equal to pyscf's `mol.symm_orb`;
!     6-31G* and cc-pVDZ, Cartesian and spherical;
!   - THE CENTRAL ORACLE: symmetry vs C1 on the symmetrised geometry, RHF,
!     UHF (the ethylene cation) and B3LYP at grid level 5, |dE| <= 1e-10 (HF)
!     / 1e-9 (DFT), orbital energies, density;
!   - anchors (sym_d2.ref): electrons per irrep and orbital energies per irrep;
!   - the requested subgroup D2 -> C2 (pyscf `symmetry='C2'`): reported group
!     C2 with `detected` D2, the pyscf frame (the D2 z axis; the input axes
!     when they hold), SO dimensions, electrons per irrep, orbital energies,
!     and the energy equal to the D2 and C1 runs;
!   - gradient: symmetrised gradient equivariant and equal to C1, and the
!     projection itself on an arbitrary (non-symmetric) gradient;
!   - `symmetry = auto` and `symmetry = D2` choose D2, `C2` is accepted.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_d2_suite
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
   use trc_sym_run, only: trc_sym_request_t, sym_choose_group, SYM_AUTO, SYM_EXPLICIT
   use trc_sym_scf, only: trc_sym_scf_t, sym_grad_symmetrize
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   use golden_io, only: golden_put
   implicit none
   private
   public :: collect_sym_d2

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the brief's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the brief's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies, symmetry vs C1
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density, symmetry vs C1
   real(dp), parameter :: TOL_ANCHOR_HF = 1.0e-6_dp   !! terco vs pyscf, HF orbital energies
   real(dp), parameter :: TOL_ANCHOR_DFT = 1.0e-5_dp  !! B3LYP: the two codes use different grids

   integer, parameter :: NMOL = 2, NBAS = 2
   integer, parameter :: DFT_GRID_LEVEL = 5  !! terco grid level of the B3LYP case (see the plan, "Known issues")
   character(len=16), parameter :: MOL(NMOL) = [character(len=16) :: "c2h4_twisted", "c2f4_twisted"]
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
      integer :: pm(4) = [1, 2, 3, 4]  !! pyscf irrep index of terco's irrep (B2, B3 exchanged when pyscf swapped x, y)
   end type pair_t

contains

   subroutine collect_sym_d2(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("detection: D2, pyscf frame (rotated, shifted and aligned input)", test_detection), &
                  new_unittest("SO basis vs pyscf dimensions, block-diagonal S and H", test_salc), &
                  new_unittest("RHF C2H4 cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_c2h4_rhf), &
                  new_unittest("RHF C2F4 6-31G* (Cartesian): symmetry vs C1 and pyscf", test_c2f4_rhf), &
                  new_unittest("UHF C2H4 cation cc-pVDZ: symmetry vs C1 and pyscf", test_c2h4_uhf), &
                  new_unittest("B3LYP C2H4 cc-pVDZ: symmetry vs C1 and pyscf", test_c2h4_b3lyp), &
                  new_unittest("D2 -> C2: group, frame, SO dimensions vs pyscf", test_sub_frame), &
                  new_unittest("D2 -> C2: energy equals D2 and C1, irreps vs pyscf", test_sub_scf), &
                  new_unittest("gradient: equivariant, equals C1, projection", test_gradient), &
                  new_unittest("auto and explicit D2 choose D2, C2 is accepted", test_choose) &
                  ]
   end subroutine collect_sym_d2

   subroutine detect_mol(name, elem, geom, error, requested)
      character(len=*), intent(in) :: name
      integer, allocatable, intent(out) :: elem(:)
      type(trc_sym_geom_t), intent(out) :: geom
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in), optional :: requested
      integer :: nat
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      call read_xyz(trim(name)//'.xyz', nat, elem, r)
      geom = sym_detect(r, elem, err, basis_id=elem, requested=requested)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%detected, "D2", "detected group of "//trim(name))
      if (allocated(error)) return
      if (present(requested)) then
         call check(error, geom%group%name, requested, "group of "//trim(name))
      else
         call check(error, geom%group%name, "D2", "group of "//trim(name))
      end if
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
            print '(a)', 'check_sym_d2: no record '//tag//' in sym_d2.ref'
            stop 1
         end if
         if (index(line, tag//' ') == 1) exit
      end do
      backspace (u)
   end subroutine goto_tag

   subroutine pyscf_frame(mol, suffix, ax, orig)
      !! pyscf's frame of `mol` (rows x, y, z) and origin from the reference.
      character(len=*), intent(in) :: mol, suffix
      real(dp), intent(out) :: ax(3, 3)
      real(dp), intent(out), optional :: orig(3)
      real(dp) :: v(9)
      integer :: u
      call ref_require('sym_d2.ref', u)
      call goto_tag(u, 'AXES_'//trim(mol)//suffix)
      call ref_read_vector(u, 'AXES_'//trim(mol)//suffix, 9, v)
      ax = transpose(reshape(v, [3, 3]))
      if (present(orig)) then
         call goto_tag(u, 'ORIG_'//trim(mol))
         call ref_read_vector(u, 'ORIG_'//trim(mol), 3, orig)
      end if
      call ref_close(u)
   end subroutine pyscf_frame

   function frame_map(rot, ax, swapped) result(ok)
      !! Whether terco's frame `rot` is pyscf's `ax` (rows, up to sign), or
      !! pyscf's with x and y exchanged (`swapped`).
      real(dp), intent(in) :: rot(3, 3), ax(3, 3)
      logical, intent(out) :: swapped
      logical :: ok
      real(dp) :: m(3, 3)
      m = abs(matmul(rot, transpose(ax)))
      swapped = .false.
      ok = all(abs(m(1, :) - [1, 0, 0]) < 1.0e-6_dp) .and. all(abs(m(2, :) - [0, 1, 0]) < 1.0e-6_dp) .and. &
           all(abs(m(3, :) - [0, 0, 1]) < 1.0e-6_dp)
      if (ok) return
      swapped = .true.
      ok = all(abs(m(1, :) - [0, 1, 0]) < 1.0e-6_dp) .and. all(abs(m(2, :) - [1, 0, 0]) < 1.0e-6_dp) .and. &
           all(abs(m(3, :) - [0, 0, 1]) < 1.0e-6_dp)
   end function frame_map

   function irrep_map(geom, mol, error) result(pm)
      !! pyscf irrep index of each D2 irrep of terco (B2 and B3 exchanged
      !! when pyscf's x and y are terco's y and x).
      type(trc_sym_geom_t), intent(in) :: geom
      character(len=*), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error
      integer :: pm(4)
      real(dp) :: ax(3, 3)
      logical :: swapped, ok
      call pyscf_frame(mol, '', ax)
      ok = frame_map(geom%rot, ax, swapped)
      call check(error, ok, trim(mol)//": the frame is neither pyscf's nor pyscf's with x, y exchanged")
      pm = [1, 2, 3, 4]
      if (swapped) pm = [1, 2, 4, 3]
   end function irrep_map

   subroutine test_detection(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp) :: ax(3, 3), orig(3), dev, resid
      integer :: im, i, j, k, nat, ic
      logical :: swapped, ok, moved
      real(dp) :: sg
      do im = 1, NMOL
         call detect_mol(MOL(im), elem, geom, error)
         if (allocated(error)) return
         call check(error, geom%group%order, 4, "four operations")
         if (allocated(error)) return
         call check(error, all(geom%group%ops(2:4) == [6, 5, 3]), "operations C2x, C2y, C2z (masks 6, 5, 3)")
         if (allocated(error)) return
         call pyscf_frame(MOL(im), '', ax, orig)
         ok = frame_map(geom%rot, ax, swapped)
         call check(error, ok, trim(MOL(im))//": frame is pyscf's (or pyscf's with x, y exchanged)")
         if (allocated(error)) return
         print '(3a,l2)', '  ', trim(MOL(im)), ': pyscf x and y exchanged relative to terco: ', swapped
         if (MOL(im) == 'c2f4_twisted') then
            ! axes aligned with the input: pyscf keeps the input axes, and so must terco
            call check(error, .not. swapped .and. maxval(abs(abs(geom%rot) - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))) &
                       < 1.0e-6_dp, "aligned input: the frame is the input axes")
            if (allocated(error)) return
         end if
         call check(error, maxval(abs(geom%shift - orig)) < 1.0e-6_dp, "origin of "//trim(MOL(im)))
         if (allocated(error)) return
         ! the image of an atom under operation k is the mask-signed coordinate of its partner
         nat = size(elem)
         dev = 0.0_dp
         do k = 2, 4
            moved = .false.
            do i = 1, nat
               j = geom%perm(i, k)
               if (j /= i) moved = .true.
               do ic = 1, 3
                  sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
                  dev = max(dev, abs(geom%r_std(ic, j) - sg*geom%r_std(ic, i)))
               end do
            end do
            call check(error, moved, trim(MOL(im))//": operation "//to_string(k)//" moves atoms")
            if (allocated(error)) return
         end do
         call check(error, dev < 1.0e-14_dp, trim(MOL(im))//": images are the signed coordinates, deviation "// &
                    to_string(dev, "es10.2"))
         if (allocated(error)) return
         resid = geom%residual
         call check(error, resid < 1.0e-5_dp, trim(MOL(im))//": the input geometry is symmetric to 1e-5 Bohr")
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
      integer :: i, ir, n
      real(dp) :: worst_s, worst_h, worst_u

      ao = sym_ao_build(bas, geom, err)
      call check(error, .not. err%has_error(), label//": sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      n = ao%nao
      call check(error, sum(ao%irrep_dim), n, label//": irrep dimensions add up to nao")
      if (allocated(error)) return
      call check(error, any(ao%ao_map(:, 2) /= [(i, i=1, n)]), label//": the C2z rotation moves AOs")
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
      worst_s = offblock(ao%u, s, col_irrep)
      worst_h = offblock(ao%u, h, col_irrep)
      print '(2a,i5,a,es9.2,a,es9.2,a,es9.2)', '  ', label, n, '  |UtU-1| ', worst_u, &
         '  off(S) ', worst_s, '  off(H) ', worst_h
      call check(error, worst_s <= 1.0e-12_dp, label//": U^T S U off-block elements")
      if (allocated(error)) return
      call check(error, worst_h <= 1.0e-12_dp, label//": U^T H U off-block elements")
      if (allocated(error)) return
      so_dims = ao%irrep_dim
   end subroutine analyse

   subroutine test_salc(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), dims(:)
      real(dp) :: ref(4)
      integer :: im, ib, ik, u, nao_ref, nao, pm(4)
      character(len=64) :: label

      call ref_require('sym_d2.ref', u)
      do im = 1, NMOL
         call detect_mol(MOL(im), elem, geom, error)
         if (allocated(error)) then
            call ref_close(u)
            return
         end if
         pm = irrep_map(geom, MOL(im), error)
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
               call analyse(bas, geom, trim(label), error, dims)
               nao = merge(bas%nao_sph, bas%nao, bas%spherical)
               call bas%release()
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call goto_tag(u, 'NAO_'//trim(MOL(im))//'_D2_'//trim(BAS_NAME(ib)))
               call ref_read_int(u, 'NAO_'//trim(MOL(im))//'_D2_'//trim(BAS_NAME(ib)), nao_ref)
               call ref_read_vector(u, 'SO_'//trim(MOL(im))//'_D2_'//trim(BAS_NAME(ib)), 4, ref)
               call check(error, nao, nao_ref, trim(label)//": AO count against pyscf")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call check(error, size(dims), 4, trim(label)//": number of irreps")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
               call check(error, all(dims == nint(ref(pm))), trim(label)//": SO dimensions (A, B1, B2, B3) against pyscf")
               if (allocated(error)) then
                  call ref_close(u)
                  return
               end if
            end do
         end do
      end do
      call ref_close(u)
   end subroutine test_salc

   subroutine make_pair(error, xyz, basis_json, spherical, functional, charge, mult, p, sub, sym_first)
      !! Detect D2 (or the requested subgroup `sub`), build the basis on the
      !! standard frame, and run the SCF with symmetry and as C1. With
      !! `sym_first` the molecule is first symmetrised to D2, so that D2, C2
      !! and C1 runs see one geometry.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json, functional
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t), intent(inout) :: p
      character(len=*), intent(in), optional :: sub
      logical, intent(in), optional :: sym_first
      integer :: nat, nelec, dev
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_sym_scf_t) :: sc
      type(trc_sym_geom_t) :: g0
      character(len=3) :: want

      dev = trc_bind_device(0)
      call read_xyz(xyz, nat, elem, r)
      if (present(sym_first)) then
         if (sym_first) then
            g0 = sym_detect(r, elem, err, basis_id=elem)
            call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
            if (allocated(error)) return
            r = g0%r_sym
         end if
      end if
      p%geom = sym_detect(r, elem, err, basis_id=elem, requested=sub)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
      if (allocated(error)) return
      want = "D2"
      if (present(sub)) want = sub
      call check(error, p%geom%group%name, want, "group of "//xyz)
      if (allocated(error)) return
      p%pm = [1, 2, 3, 4]
      if (.not. present(sub)) then
         p%pm = irrep_map(p%geom, xyz(1:index(xyz, '.') - 1), error)
         if (allocated(error)) return
      end if
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
            print '(a)', 'check_sym_d2: no case '//tag//' in sym_d2.ref'
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
      !! (terco's irrep ir is pyscf's irrep p%pm(ir)).
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
      call ref_require('sym_d2.ref', u)
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
            want(ir, s) = nint(nelec(ioff + p%pm(ir)))
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
            call goto_tag(u, 'EPS'//to_string(p%pm(ir))//'_'//tag//'_'//to_string(s))
            call ref_read_vector(u, 'EPS'//to_string(p%pm(ir))//'_'//tag//'_'//to_string(s), nref, eir)
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
      call make_pair(error, 'c2h4_twisted.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF C2H4 cc-pVDZ', TOL_E_HF)
      call golden_put('E_c2h4_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_c2h4_rhf

   subroutine test_c2f4_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2f4_twisted.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF C2F4 6-31G*', TOL_E_HF)
      call golden_put('E_c2f4_rhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2f4_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_c2f4_rhf

   subroutine test_c2h4_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2h4_twisted.xyz', 'basis_sets/cc-pvdz.json', .true., '', 1, 2, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF C2H4+ cc-pVDZ', TOL_E_HF)
      call golden_put('E_c2h4_uhf', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4_uhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_c2h4_uhf

   subroutine test_c2h4_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'c2h4_twisted.xyz', 'basis_sets/cc-pvdz.json', .true., 'b3lyp', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'B3LYP C2H4 cc-pVDZ', TOL_E_DFT)
      call golden_put('E_c2h4_b3lyp', p%rs%energy)
      if (.not. allocated(error)) call compare_anchor(error, p, 'c2h4_b3lyp', TOL_ANCHOR_DFT)
      call p%bas%release()
   end subroutine test_c2h4_b3lyp

   subroutine test_sub_frame(error)
      !! D2 -> C2: the reported group, `detected`, the frame and the SO dimensions as pyscf.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_basis_t) :: bas
      integer, allocatable :: elem(:), dims(:)
      real(dp) :: ax(3, 3), ref(2)
      integer :: im, ib, u, nao_ref, nao
      logical :: swapped, ok
      character(len=64) :: label
      do im = 1, NMOL
         call detect_mol(MOL(im), elem, geom, error, requested='C2')
         if (allocated(error)) return
         call check(error, geom%group%order, 2, "two operations")
         if (allocated(error)) return
         call pyscf_frame(MOL(im), '_C2', ax)
         ok = frame_map(geom%rot, ax, swapped)
         ! only the C2 axis matters for the labels (A, B): the z axis must be pyscf's
         call check(error, abs(abs(dot_product(geom%rot(3, :), ax(3, :))) - 1.0_dp) < 1.0e-6_dp, &
                    trim(MOL(im))//": z of the requested C2 is pyscf's axis")
         if (allocated(error)) return
         if (MOL(im) == 'c2f4_twisted') then
            call check(error, maxval(abs(abs(geom%rot) - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))) < 1.0e-6_dp, &
                       "aligned input: pyscf keeps the input axes for C2, so does terco")
            if (allocated(error)) return
         end if
         call ref_require('sym_d2.ref', u)
         do ib = 1, NBAS
            label = trim(MOL(im))//' C2 '//trim(BAS_NAME(ib))
            call make_basis(BAS_FILE(ib), elem, geom, ib == 2, bas, error)
            if (allocated(error)) then
               call ref_close(u)
               return
            end if
            call analyse(bas, geom, trim(label), error, dims)
            nao = merge(bas%nao_sph, bas%nao, bas%spherical)
            call bas%release()
            if (allocated(error)) then
               call ref_close(u)
               return
            end if
            call goto_tag(u, 'NAO_'//trim(MOL(im))//'_C2_'//trim(BAS_NAME(ib)))
            call ref_read_int(u, 'NAO_'//trim(MOL(im))//'_C2_'//trim(BAS_NAME(ib)), nao_ref)
            call ref_read_vector(u, 'SO_'//trim(MOL(im))//'_C2_'//trim(BAS_NAME(ib)), 2, ref)
            call check(error, nao == nao_ref .and. size(dims) == 2, trim(label)//": AO count and irreps against pyscf")
            if (allocated(error)) then
               call ref_close(u)
               return
            end if
            call check(error, all(dims == nint(ref)), trim(label)//": SO dimensions (A, B) against pyscf")
            if (allocated(error)) then
               call ref_close(u)
               return
            end if
         end do
         call ref_close(u)
      end do
   end subroutine test_sub_frame

   subroutine test_sub_scf(error)
      !! D2 -> C2: on the D2-symmetrised geometry the C2, D2 and C1 energies agree; irreps and energies vs pyscf.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: pc, pd
      real(dp) :: de
      call make_pair(error, 'c2h4_twisted.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, pc, sub='C2', sym_first=.true.)
      if (allocated(error)) return
      call make_pair(error, 'c2h4_twisted.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, pd, sym_first=.true.)
      if (allocated(error)) return
      call compare_pair(error, pc, 'C2 of D2, C2H4 6-31G*', TOL_E_HF)
      if (allocated(error)) return
      de = abs(pc%rs%energy - pd%rs%energy)
      print '(a,es10.2)', '  C2H4 6-31G*: |E(C2) - E(D2)| ', de
      call check(error, de <= TOL_E_HF, "E(C2) = E(D2): "//to_string(de, "es10.2"))
      if (allocated(error)) return
      call compare_anchor(error, pc, 'c2h4_c2_rhf', TOL_ANCHOR_HF)
      call pc%bas%release()
      call pd%bas%release()
      if (allocated(error)) return
      call make_pair(error, 'c2f4_twisted.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, pc, sub='C2', sym_first=.true.)
      if (allocated(error)) return
      call make_pair(error, 'c2f4_twisted.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, pd, sym_first=.true.)
      if (allocated(error)) return
      call compare_pair(error, pc, 'C2 of D2, C2F4 cc-pVDZ', TOL_E_HF)
      if (allocated(error)) return
      de = abs(pc%rs%energy - pd%rs%energy)
      print '(a,es10.2)', '  C2F4 cc-pVDZ: |E(C2) - E(D2)| ', de
      call check(error, de <= TOL_E_HF, "E(C2) = E(D2): "//to_string(de, "es10.2"))
      if (allocated(error)) return
      call compare_anchor(error, pc, 'c2f4_c2_rhf', TOL_ANCHOR_HF)
      call pc%bas%release()
      call pd%bas%release()
   end subroutine test_sub_scf

   subroutine test_gradient(error)
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, 'c2h4_twisted.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'c2f4_twisted.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'c2h4_twisted.xyz', 'basis_sets/6-31gs.json', .false., 1, 2)
      if (allocated(error)) return
      call noise_case(error, 'c2h4_twisted.xyz')
      if (allocated(error)) return
      call noise_case(error, 'c2f4_twisted.xyz')
   end subroutine test_gradient

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
      geom = sym_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_detect: "//err%get_message())
      if (allocated(error)) return
      allocate (g(3, nat))
      do a = 1, nat
         do ic = 1, 3
            g(ic, a) = sin(1.7_dp*real(3*a + ic, dp)) + 0.3_dp
         end do
      end do
      gs = g
      call sym_grad_symmetrize(geom%group, geom%perm, gs)
      moved = maxval(abs(gs - g))
      call check(error, moved > 0.1_dp, xyz//": the projection changed a non-symmetric gradient")
      if (allocated(error)) return
      g2 = gs
      call sym_grad_symmetrize(geom%group, geom%perm, g2)
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
                            sym_geom=p%geom)
      call trc_scf_gradient(p%bas, lp, p%r1%nspin, p%nalpha, p%nbeta, p%r1%cmo, p%r1%eps, 1.0e-12_dp, g1)
      call lp%release()
      call p%bas%release()
      dg = maxval(abs(gsym - g1))
      print '(a,a,a,es10.2)', '  ', xyz, ': |sym gradient - C1 gradient| ', dg
      call check_equivariance(error, p%geom, gsym, xyz//" symmetrised gradient", asym, forbidden)
      if (allocated(error)) return
      call check(error, dg <= 1.0e-9_dp, xyz//": symmetrised gradient differs from C1 by "//to_string(dg, "es10.2"))
   end subroutine grad_case

   subroutine test_choose(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_sym_request_t) :: req
      type(error_t) :: e1, e2, e3, e4
      integer, allocatable :: elem(:)
      character(len=3) :: g
      call detect_mol('c2h4_twisted', elem, geom, error)
      if (allocated(error)) return
      req%mode = SYM_AUTO
      g = sym_choose_group(req, geom, e1)
      call check(error, g == 'D2' .and. .not. e1%has_error(), 'auto picks D2')
      if (allocated(error)) return
      req%mode = SYM_EXPLICIT
      req%group = 'D2'
      g = sym_choose_group(req, geom, e2)
      call check(error, g == 'D2' .and. .not. e2%has_error(), 'explicit D2 is accepted')
      if (allocated(error)) return
      call detect_mol('c2h4_twisted', elem, geom, error, requested='C2')
      if (allocated(error)) return
      req%group = 'C2'
      g = sym_choose_group(req, geom, e3)
      call check(error, g == 'C2' .and. .not. e3%has_error(), 'explicit C2 is accepted and reported as C2: '//e3%get_message())
      if (allocated(error)) return
      req%group = 'Cs'
      g = sym_choose_group(req, geom, e4)
      call check(error, e4%has_error(), 'Cs is not a subgroup of D2')
   end subroutine test_choose

end module check_sym_d2_suite

program check_sym_d2
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_sym_d2_suite, only: collect_sym_d2
   use golden_io, only: golden_begin, golden_end
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, is
   stat = 0
   call golden_begin("sym_d2")
   testsuites = [new_testsuite("sym_d2", collect_sym_d2)]
   do is = 1, size(testsuites)
      write (error_unit, '(a)') "Testing: "//testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
   end do
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
   call golden_end()
end program check_sym_d2
