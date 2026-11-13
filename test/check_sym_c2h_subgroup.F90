module check_sym_c2h_subgroup_suite
   !!
   !! Requested subgroups of C2h: Cs, C2 and Ci requested on C2h molecules,
   !! against pyscf 2.14.0 (symmetry='Cs' / 'C2' / 'Ci'; for these pyscf keeps the
   !! C2h frame, `as_subgroup` changes the axes for C2h only) and against the C2h
   !! and C1 runs of terco.
   !!
   !! Molecules: trans-diazene (planar) and anti ClH2C-CH2Cl (nonplanar) in a
   !! general orientation, and anti ClH2C-CH2Cl with axis-aligned input twice:
   !! clh2c_ch2cl_alignx (C2 along x: sigma_h and C2 do not hold in the input axes,
   !! so pyscf turns to the C2h frame for Cs and C2, but keeps the input axes for
   !! Ci, which holds in any frame) and clh2c_ch2cl_alignz (C2 along z: every
   !! subgroup holds in the input axes). The latter two check pyscf
   !! keep-the-input-axes rule (`input_frame_holds`).
   !!
   !! test_frame:   the reported group is the requested one, `detected` is C2h, the
   !!               frame (rows x, y, z up to sign) and origin are pyscf
   !!               `_symm_axes` and `_symm_orig`, and the operations are
   !!               bijections that leave the symmetrised geometry invariant.
   !! test_so_dims: SO dimensions per irrep equal mol.symm_orb (6-31G* Cartesian,
   !!               cc-pVDZ spherical).
   !! test_scf:     on the geometry symmetrised to C2h, the energy of the subgroup
   !!               run equals the C2h run and the C1 run (<= 1e-10 Eh); on the raw input,
   !!               the electrons per irrep, the orbital energies per irrep and the
   !!               energy equal pyscf (sym_c2h_subgroup.ref).
   !! test_refusals: a request that is not a proper subgroup keeps the full group.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_vector
   implicit none
   private
   public :: collect_sym_c2h_subgroup

   real(dp), parameter :: TOL_E = 1.0e-10_dp      !! the central-oracle gate, HF
   real(dp), parameter :: TOL_ANCHOR = 1.0e-6_dp  !! terco vs pyscf, orbital energies and energy
   integer, parameter :: NMOL = 4, NSUB = 3
   character(len=20), parameter :: MOL(NMOL) = [character(len=20) :: "n2h2_trans", "clh2c_ch2cl_anti", &
                                                      "clh2c_ch2cl_alignx", "clh2c_ch2cl_alignz"]
   character(len=2), parameter :: SUB(NSUB) = [character(len=2) :: "Cs", "C2", "Ci"]

   type :: scf_case_t
      character(len=20) :: tag, mol
      character(len=2) :: sub
      character(len=32) :: basis
      logical :: spherical
      integer :: mult
      integer :: charge = 0
   end type scf_case_t

   type(scf_case_t), parameter :: CASES(10) = [ &
                                  scf_case_t("n2h2_cs_rhf", "n2h2_trans", "Cs", "basis_sets/cc-pvdz.json", .true., 1), &
                                  scf_case_t("n2h2_c2_rhf", "n2h2_trans", "C2", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("n2h2_ci_uhf", "n2h2_trans", "Ci", "basis_sets/cc-pvdz.json", .true., 2, 1), &
                                  scf_case_t("anti_cs_rhf", "clh2c_ch2cl_anti", "Cs", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("anti_c2_rhf", "clh2c_ch2cl_anti", "C2", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("anti_ci_rhf", "clh2c_ch2cl_anti", "Ci", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("alignx_cs_rhf", "clh2c_ch2cl_alignx", "Cs", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("alignx_ci_rhf", "clh2c_ch2cl_alignx", "Ci", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("alignz_cs_rhf", "clh2c_ch2cl_alignz", "Cs", "basis_sets/6-31gs.json", .false., 1), &
                                  scf_case_t("alignz_c2_rhf", "clh2c_ch2cl_alignz", "C2", "basis_sets/6-31gs.json", .false., 1)]

contains

   subroutine collect_sym_c2h_subgroup(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Cs, C2 and Ci of C2h: reported group, frame and origin as in pyscf", test_frame), &
                  new_unittest("Cs, C2 and Ci of C2h: SO dimensions vs pyscf", test_so_dims), &
                  new_unittest("Cs, C2 and Ci of C2h: energy equals C2h and C1, irreps vs pyscf", test_scf), &
                  new_unittest("requests that are not subgroups, C1, and the plain call", test_refusals) &
                  ]
   end subroutine collect_sym_c2h_subgroup

   subroutine detect(xyz, sub, geom, elem, error, r_in)
      !! `sym_geom_detect` on the geometry of `xyz.xyz` (or `r_in`), with the requested group `sub` (not given: none).
      character(len=*), intent(in) :: xyz
      character(len=*), intent(in) :: sub
      type(trc_sym_geom_t), intent(out) :: geom
      integer, allocatable, intent(out) :: elem(:)
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in), optional :: r_in(:, :)
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      integer :: nat
      call read_xyz(trim(xyz)//'.xyz', nat, elem, r)
      if (present(r_in)) r = r_in
      if (len_trim(sub) > 0) then
         geom = sym_geom_detect(r, elem, err, basis_id=elem, requested=trim(sub))
      else
         geom = sym_geom_detect(r, elem, err, basis_id=elem)
      end if
      call check(error, .not. err%has_error(), "sym_geom_detect "//trim(xyz)//" "//trim(sub)//": "//err%get_message())
   end subroutine detect

   subroutine test_frame(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)
      real(dp) :: ax(9), orig(3), g(3, 3)
      integer :: im, is, k, i, j, kk, u, nat
      logical, allocatable :: seen(:)

      call ref_require('sym_c2h_subgroup.ref', u)
      do im = 1, NMOL
         do is = 1, NSUB
            call detect(MOL(im), SUB(is), geom, elem, error)
            if (allocated(error)) return
            nat = size(elem)
            call ref_read_vector(u, 'AXES_'//trim(MOL(im))//'_'//SUB(is), 9, ax)
            call ref_read_vector(u, 'ORIG_'//trim(MOL(im))//'_'//SUB(is), 3, orig)
            call check(error, geom%detected, "C2h", "detected group of "//trim(MOL(im)))
            if (allocated(error)) return
            call check(error, geom%group%name, SUB(is), trim(MOL(im))//": reported group")
            if (allocated(error)) return
            call check(error, geom%group%order, 2, "two operations")
            if (allocated(error)) return
            call check(error, maxval(abs(geom%shift - orig)) < 1.0e-6_dp, "origin of "//trim(MOL(im))//" "//SUB(is))
            if (allocated(error)) return
            do k = 1, 3
               call check(error, abs(abs(dot_product(geom%rot(k, :), ax(3*k - 2:3*k))) - 1.0_dp) < 1.0e-7_dp, &
                          trim(MOL(im))//" "//SUB(is)//": axis "//to_string(k)//" against pyscf")
               if (allocated(error)) return
            end do
            g = matmul(geom%rot, transpose(geom%rot))
            do i = 1, 3
               g(i, i) = g(i, i) - 1.0_dp
            end do
            call check(error, maxval(abs(g)) < 1.0e-14_dp, "orthogonal frame")
            if (allocated(error)) return
            call check(error, all(geom%perm(:, 1) == [(i, i=1, nat)]), "E is the identity")
            if (allocated(error)) return
            allocate (seen(nat))
            do kk = 1, 2
               seen = .false.
               do i = 1, nat
                  seen(geom%perm(i, kk)) = .true.
               end do
               call check(error, all(seen), "perm is a bijection for "//trim(MOL(im)))
               if (allocated(error)) return
               do i = 1, nat
                  do j = 1, 3
                     call check(error, abs(merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(kk), j - 1))*geom%r_std(j, i) - &
                                           geom%r_std(j, geom%perm(i, kk))) < 1.0e-14_dp, "symmetrised coordinates are invariant")
                     if (allocated(error)) return
                  end do
               end do
            end do
            deallocate (seen)
         end do
      end do
      call ref_close(u)
   end subroutine test_frame

   subroutine test_so_dims(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_basis_t) :: bas
      type(trc_sym_ao_t) :: ao
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp) :: ref(2), ax(9), orig(3)
      integer :: im, is, ib, u, nao_ref
      character(len=8), parameter :: BN(2) = [character(len=8) :: "631gs", "ccpvdz"]
      character(len=32), parameter :: BF(2) = [character(len=32) :: "basis_sets/6-31gs.json", "basis_sets/cc-pvdz.json"]

      call ref_require('sym_c2h_subgroup.ref', u)
      do im = 1, NMOL  ! the frame records come first in the file
         do is = 1, NSUB
            call ref_read_vector(u, 'AXES_'//trim(MOL(im))//'_'//SUB(is), 9, ax)
            call ref_read_vector(u, 'ORIG_'//trim(MOL(im))//'_'//SUB(is), 3, orig)
         end do
      end do
      do im = 1, NMOL
         do is = 1, NSUB
            call detect(MOL(im), SUB(is), geom, elem, error)
            if (allocated(error)) return
            do ib = 1, 2
               call trc_basis_from_json(BF(ib), size(elem), elem, geom%r_std, bas, err, force_cart=ib == 1)
               call check(error, .not. err%has_error(), "basis: "//err%get_message())
               if (allocated(error)) return
               bas%spherical = ib == 2
               ao = sym_ao_build(bas, geom, err)
               call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
               if (allocated(error)) return
               call ref_read_int(u, 'NAO_'//trim(MOL(im))//'_'//SUB(is)//'_'//trim(BN(ib)), nao_ref)
               call ref_read_vector(u, 'SO_'//trim(MOL(im))//'_'//SUB(is)//'_'//trim(BN(ib)), 2, ref)
               call check(error, sum(ao%irrep_dim) == nao_ref .and. all(ao%irrep_dim == nint(ref)), &
                          trim(MOL(im))//" "//SUB(is)//" "//trim(BN(ib))//": SO dimensions "//to_string(ao%irrep_dim(1))// &
                          " "//to_string(ao%irrep_dim(2))//", pyscf "//to_string(nint(ref(1)))//" "//to_string(nint(ref(2))))
               call bas%release()
               if (allocated(error)) return
            end do
         end do
      end do
      call ref_close(u)
   end subroutine test_so_dims

   subroutine run_geom(error, geom, elem, c, rs, ao, with_sym)
      !! SCF for the case `c` on the standard frame of `geom`, with the irrep
      !! blocks of `geom` or (with_sym false) as plain C1.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t), intent(in) :: geom
      integer, intent(in) :: elem(:)
      type(scf_case_t), intent(in) :: c
      type(trc_scf_result_t), intent(out) :: rs
      type(trc_sym_ao_t), intent(out) :: ao
      logical, intent(in) :: with_sym
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_sym_scf_t) :: sc
      integer :: nelec, na, nb, dev

      dev = trc_bind_device(0)
      call trc_basis_from_json(c%basis, size(elem), elem, geom%r_std, bas, err, force_cart=.not. c%spherical)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = c%spherical
      call bas%to_device()
      nelec = sum(elem) - c%charge
      na = (nelec + c%mult - 1)/2
      nb = nelec - na
      opts%conv_energy = 1.0e-10_dp   ! just above the accumulation-order noise floor
      opts%conv_diis = 1.0e-9_dp
      opts%max_iteration = 100
      if (with_sym) then
         ao = sym_ao_build(bas, geom, err)
         call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
         if (allocated(error)) then
            call bas%release()
            return
         end if
         sc%ao = ao
         call trc_scf_run(bas, na, nb, opts, rs, sym=sc)
         call sc%release()
      else
         call trc_scf_run(bas, na, nb, opts, rs)
      end if
      call bas%release()
      call check(error, rs%converged, trim(c%tag)//": SCF did not converge: "//trim(rs%message))
   end subroutine run_geom

   subroutine test_scf(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ic
      do ic = 1, size(CASES)
         call one_case(error, CASES(ic))
         if (allocated(error)) return
      end do
   end subroutine test_scf

   subroutine one_case(error, c)
      type(error_type), allocatable, intent(out) :: error
      type(scf_case_t), intent(in) :: c
      type(trc_sym_geom_t) :: gfull, gsub, graw
      type(trc_scf_result_t) :: rfull, rsub, r1, rraw
      type(trc_sym_ao_t) :: afull, asub, a1, araw
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :)
      integer :: nat

      ! Oracle: the geometry symmetrised to C2h, so that the C2h, Cs/C2 and C1 runs
      ! are the same molecule.
      call detect(c%mol, '', gfull, elem, error)
      if (allocated(error)) return
      nat = size(elem)
      r = gfull%r_sym
      call detect(c%mol, c%sub, gsub, elem, error, r_in=r)
      if (allocated(error)) return
      call detect(c%mol, '', gfull, elem, error, r_in=r)
      if (allocated(error)) return
      call check(error, gfull%group%name, "C2h", trim(c%tag)//": full group")
      if (allocated(error)) return
      call check(error, gsub%group%name, c%sub, trim(c%tag)//": requested group")
      if (allocated(error)) return
      call run_geom(error, gfull, elem, c, rfull, afull, .true.)
      if (allocated(error)) return
      call run_geom(error, gsub, elem, c, rsub, asub, .true.)
      if (allocated(error)) return
      call run_geom(error, gsub, elem, c, r1, a1, .false.)
      if (allocated(error)) return
      print '(a,a,a,es10.2,a,es10.2)', '  ', trim(c%tag), ': |E(sub) - E(C2h)| ', abs(rsub%energy - rfull%energy), &
         '  |E(sub) - E(C1)| ', abs(rsub%energy - r1%energy)
      call check(error, abs(rsub%energy - rfull%energy) <= TOL_E, trim(c%tag)//": E(sub) - E(C2h) = "// &
                 to_string(rsub%energy - rfull%energy, "es10.2"))
      if (allocated(error)) return
      call check(error, abs(rsub%energy - r1%energy) <= TOL_E, trim(c%tag)//": E(sub) - E(C1) = "// &
                 to_string(rsub%energy - r1%energy, "es10.2"))
      if (allocated(error)) return

      ! Anchor: the raw input geometry, as pyscf sees it.
      call detect(c%mol, c%sub, graw, elem, error)
      if (allocated(error)) return
      call run_geom(error, graw, elem, c, rraw, araw, .true.)
      if (allocated(error)) return
      call compare_anchor(error, rraw, araw, c, sum(elem) - c%charge)
   end subroutine one_case

   subroutine compare_anchor(error, rs, ao, c, nelec)
      !! Electrons per irrep, orbital energies per irrep and energy against pyscf.
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_result_t), intent(in) :: rs
      type(trc_sym_ao_t), intent(in) :: ao
      type(scf_case_t), intent(in) :: c
      integer, intent(in) :: nelec
      integer :: u, nspin, ir, s, m, nocc(2), cnt(2, 2), k, nref
      real(dp) :: eref, nel(4), worst
      real(dp), allocatable :: ev(:), eir(:)
      character(len=512) :: line
      integer :: ios

      call ref_require('sym_c2h_subgroup.ref', u)
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) then
            print '(a)', 'check_sym_c2h_subgroup: no case '//trim(c%tag)//' in sym_c2h_subgroup.ref'
            stop 1
         end if
         if (index(line, 'NSPIN_'//trim(c%tag)//' ') == 1) exit
      end do
      read (line(len('NSPIN_'//trim(c%tag)) + 1:), *) nspin
      read (u, '(a)') line
      read (line(len('E_'//trim(c%tag)) + 1:), *) eref
      call ref_read_vector(u, 'NELEC_'//trim(c%tag), 4, nel)
      nocc(1) = (nelec + c%mult - 1)/2
      nocc(2) = nelec - nocc(1)
      cnt = 0
      do s = 1, rs%nspin
         do m = 1, nocc(s)
            ir = rs%mo_irrep(m, s)
            cnt(ir, s) = cnt(ir, s) + 1
         end do
      end do
      if (rs%nspin == 1) cnt(:, 2) = cnt(:, 1)
      call check(error, all(cnt(:, 1) == nint(nel(1:2))) .and. all(cnt(:, 2) == nint(nel(3:4))), &
                 trim(c%tag)//": electrons per irrep "//to_string(cnt(1, 1))//to_string(cnt(2, 1))//" / "// &
                 to_string(cnt(1, 2))//to_string(cnt(2, 2))//" differ from pyscf "//to_string(nint(nel(1)))// &
                 to_string(nint(nel(2)))//" / "//to_string(nint(nel(3)))//to_string(nint(nel(4))))
      if (allocated(error)) then
         call ref_close(u)
         return
      end if
      worst = 0.0_dp
      do s = 1, nspin
         do ir = 1, 2
            nref = ao%irrep_dim(ir)
            allocate (eir(nref))
            call ref_read_vector(u, 'EPS'//to_string(ir)//'_'//trim(c%tag)//'_'//to_string(s), nref, eir)
            ev = pack(rs%eps(1:ao%nao, s), rs%mo_irrep(1:ao%nao, s) == ir)
            call check(error, size(ev) == nref, trim(c%tag)//": orbitals per irrep differ from pyscf")
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
      print '(a,a,a,es10.2,a,es10.2)', '  ', trim(c%tag), ': max |eps - eps(pyscf)| ', worst, '  |E - E(pyscf)| ', &
         abs(rs%energy - eref)
      call check(error, worst <= TOL_ANCHOR, trim(c%tag)//": orbital energies differ from pyscf by "// &
                 to_string(worst, "es10.2"))
      if (allocated(error)) return
      call check(error, abs(rs%energy - eref) <= TOL_ANCHOR, trim(c%tag)//": energy differs from pyscf by "// &
                 to_string(abs(rs%energy - eref), "es10.2"))
   end subroutine compare_anchor

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:)

      ! C2v is not a subgroup of C2h: sym_geom_detect keeps the full group.
      call detect('n2h2_trans', 'C2v', geom, elem, error)
      if (allocated(error)) return
      call check(error, geom%group%name == "C2h" .and. geom%detected == "C2h", &
         "C2v requested on trans-N2H2: the full group is kept")
      if (allocated(error)) return
      ! The full group requested: the plain call.
      call detect('n2h2_trans', 'C2h', geom, elem, error)
      if (allocated(error)) return
      call check(error, geom%group%name == "C2h" .and. geom%group%order == 4, "C2h on trans-N2H2")
      if (allocated(error)) return
      ! Without a request a C2h molecule is reported as C2h, not as a subgroup.
      call detect('clh2c_ch2cl_anti', '', geom, elem, error)
      if (allocated(error)) return
      call check(error, geom%group%name == "C2h", "anti ClH2C-CH2Cl is reported as C2h")
      if (allocated(error)) return
      ! A molecule with no rule for the pair (Cs inside Cs): the group itself.
      call detect('hocl', 'Cs', geom, elem, error)
      if (allocated(error)) return
      call check(error, geom%group%name == "Cs" .and. geom%detected == "Cs", "Cs on a Cs molecule")
   end subroutine test_refusals

end module check_sym_c2h_subgroup_suite

program check_sym_c2h_subgroup
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_sym_c2h_subgroup_suite, only: collect_sym_c2h_subgroup
   implicit none
   integer :: stat
   stat = 0
   call run_testsuite(collect_sym_c2h_subgroup, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_c2h_subgroup
