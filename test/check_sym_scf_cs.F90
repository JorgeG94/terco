!
! The SO-blocked SCF (trc_sym_scf, trc_scf_run's `sym` argument) for the group
! Cs, the only one beyond C1 enabled at this stage.
!
! THE CENTRAL ORACLE. A molecule on the symmetrised, standard-frame geometry
! of trc_sym_geom is run twice, with the irrep blocks and as plain C1, and the
! two must agree: total energy (<= 1e-10 Eh for HF, <= 1e-9 for B3LYP),
! orbital energies and the density. Cases:
!   RHF   HOCl (planar: every atom is its own mirror image) cc-pVDZ, spherical
!   RHF   CH2Cl2 (nonplanar, the mirror swaps atoms) 6-31G*, Cartesian
!   UHF   HCO radical cc-pVDZ, and the CH2Cl2 radical cation 6-31G* (nonplanar)
!   B3LYP HOCl cc-pVDZ
! Anchors (pyscf 2.14.0, mol.symmetry = 'Cs', sym_scf_cs.ref): the electrons
! per irrep (`mf.get_irrep_nelec()`) and the orbital energies of each irrep.
! The gradient: the components that symmetry forces to zero are <= 1e-10 in
! the raw gradient of the symmetric run, vanish after `sym_geom` symmetrises
! it, and the symmetrised gradient equals the C1 gradient.
! Also: fixed per-irrep occupations, --soscf refused with symmetry, and the
! incremental Fock build with symmetry.
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_scf_cs_suite
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_basis_json, only: trc_basis_from_json
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_scf_grad, only: trc_scf_gradient
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_sym_scf, only: trc_sym_scf_t, sym_grad_symmetrize
   use trc_test_basis, only: read_xyz
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_vector
   implicit none
   private
   public :: collect_sym_scf_cs

   real(dp), parameter :: TOL_E_HF = 1.0e-10_dp   !! the brief's gate, HF
   real(dp), parameter :: TOL_E_DFT = 1.0e-9_dp   !! the brief's gate, DFT
   real(dp), parameter :: TOL_EPS = 1.0e-8_dp     !! orbital energies, symmetry vs C1
   real(dp), parameter :: TOL_D = 1.0e-7_dp       !! density, symmetry vs C1
   real(dp), parameter :: TOL_ANCHOR_HF = 1.0e-6_dp   !! terco vs pyscf, HF orbital energies (measured <= 2.5e-8)
   real(dp), parameter :: TOL_ANCHOR_DFT = 1.0e-5_dp  !! B3LYP: the two codes use different grids (measured 8.7e-7)
   real(dp), parameter :: TOL_GRAD = 1.0e-10_dp

   type :: pair_t
      !! One molecule run with and without symmetry.
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: ao
      type(trc_scf_result_t) :: rs, r1
      integer :: nalpha = 0, nbeta = 0
   end type pair_t

contains

   subroutine collect_sym_scf_cs(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("RHF HOCl cc-pVDZ (spherical): symmetry vs C1 and pyscf", test_hocl_rhf), &
                  new_unittest("RHF CH2Cl2 6-31G* (Cartesian, mirror moves atoms): vs C1 and pyscf", test_ch2cl2_rhf), &
                  new_unittest("UHF HCO cc-pVDZ: symmetry vs C1 and pyscf", test_hco_uhf), &
                  new_unittest("UHF CH2Cl2 cation 6-31G* (nonplanar): symmetry vs C1", test_ch2cl2_cation_uhf), &
                  new_unittest("B3LYP HOCl cc-pVDZ: symmetry vs C1 and pyscf", test_hocl_b3lyp), &
                  new_unittest("RHF CH2Cl2 without incremental Fock builds: symmetry vs C1", test_full_builds), &
                  new_unittest("gradient: forbidden components vanish, equals C1", test_gradient), &
                  new_unittest("fixed per-irrep occupations", test_fixed_occupations), &
                  new_unittest("--soscf is refused with symmetry; C1 group is the old path", test_refusals) &
                  ]
   end subroutine collect_sym_scf_cs

   subroutine make_pair(error, xyz, basis_json, spherical, functional, charge, mult, p, incremental, with_c1)
      !! Detect Cs, build the basis on the standard frame, and run the SCF
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
      call check(error, p%geom%group%name, "Cs", "group of "//xyz)
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
      call ref_require('sym_scf_cs.ref', u)
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
            print '(a)', 'check_sym_scf_cs: no case '//tag//' in sym_scf_cs.ref'
            stop 1
         end if
         if (index(line, 'NSPIN_'//tag//' ') == 1) exit
      end do
      read (line(len('NSPIN_'//tag) + 1:), *) nspin
      read (u, '(a)') line
      read (line(len('E_'//tag) + 1:), *) e
      ! positioned at the NELEC record, which the caller reads
   end subroutine find_case

   subroutine test_hocl_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'hocl.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF HOCl cc-pVDZ', TOL_E_HF)
      if (.not. allocated(error)) call compare_anchor(error, p, 'hocl_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_hocl_rhf

   subroutine test_ch2cl2_rhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF CH2Cl2 6-31G*', TOL_E_HF)
      if (.not. allocated(error)) call check(error, p%rs%n_increments > 0 .and. p%r1%n_increments > 0, &
                                             "the default runs must have used incremental Fock builds")
      if (.not. allocated(error)) call compare_anchor(error, p, 'ch2cl2_rhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_ch2cl2_rhf

   subroutine test_hco_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'hco.xyz', 'basis_sets/cc-pvdz.json', .true., '', 0, 2, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF HCO cc-pVDZ', TOL_E_HF)
      if (.not. allocated(error)) call compare_anchor(error, p, 'hco_uhf', TOL_ANCHOR_HF)
      call p%bas%release()
   end subroutine test_hco_uhf

   subroutine test_ch2cl2_cation_uhf(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., '', 1, 2, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'UHF CH2Cl2+ 6-31G*', TOL_E_HF)
      call p%bas%release()
   end subroutine test_ch2cl2_cation_uhf

   subroutine test_full_builds(error)
      !! The same oracle with every Fock build full. With the default
      !! (incremental) schedule, test_ch2cl2_rhf checks the increments are used.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, incremental=.false.)
      if (allocated(error)) return
      call compare_pair(error, p, 'RHF CH2Cl2 6-31G*, full builds', TOL_E_HF)
      if (.not. allocated(error)) call check(error, p%rs%n_increments == 0 .and. p%r1%n_increments == 0, &
                                             "incremental builds were off")
      call p%bas%release()
   end subroutine test_full_builds

   subroutine test_hocl_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      call make_pair(error, 'hocl.xyz', 'basis_sets/cc-pvdz.json', .true., 'b3lyp', 0, 1, p)
      if (allocated(error)) return
      call compare_pair(error, p, 'B3LYP HOCl cc-pVDZ', TOL_E_DFT)
      if (.not. allocated(error)) call compare_anchor(error, p, 'hocl_b3lyp', TOL_ANCHOR_DFT)
      call p%bas%release()
   end subroutine test_hocl_b3lyp

   subroutine test_gradient(error)
      !! Forbidden components of the gradient, and the symmetrised gradient
      !! against C1, for planar HOCl (RHF), nonplanar CH2Cl2 (RHF) and HCO (UHF).
      type(error_type), allocatable, intent(out) :: error
      call grad_case(error, 'hocl.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'ch2cl2.xyz', 'basis_sets/6-31gs.json', .false., 0, 1)
      if (allocated(error)) return
      call grad_case(error, 'hco.xyz', 'basis_sets/6-31gs.json', .false., 0, 2)
      if (allocated(error)) return
      call noise_case(error)
   end subroutine test_gradient

   subroutine noise_case(error)
      !! The projection itself, on a gradient that is NOT symmetric: the
      !! symmetrised result is equivariant, has no forbidden component, is
      !! unchanged by a second application, and keeps the allowed part of an
      !! already symmetric gradient.
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(error_t) :: err
      integer, allocatable :: elem(:)
      real(dp), allocatable :: r(:, :), g(:, :), gs(:, :), g2(:, :)
      real(dp) :: asym, sg, forbidden, idem
      integer :: a, k, ic, nat
      call read_xyz('ch2cl2.xyz', nat, elem, r)
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
      g2 = gs
      call sym_grad_symmetrize(geom%group, geom%perm, g2)
      idem = maxval(abs(g2 - gs))
      asym = 0.0_dp
      forbidden = 0.0_dp
      do a = 1, nat
         do k = 2, geom%group%order
            if (geom%perm(a, k) == a) forbidden = max(forbidden, abs(gs(3, a)))
            do ic = 1, 3
               sg = merge(-1.0_dp, 1.0_dp, btest(geom%group%ops(k), ic - 1))
               asym = max(asym, abs(gs(ic, geom%perm(a, k)) - sg*gs(ic, a)))
            end do
         end do
      end do
      call check(error, asym <= 1.0e-14_dp .and. forbidden == 0.0_dp .and. idem <= 1.0e-14_dp, &
                 "gradient projection: equivariance "//to_string(asym, "es10.2")//", forbidden "// &
                 to_string(forbidden, "es10.2")//", idempotence "//to_string(idem, "es10.2"))
   end subroutine noise_case

   subroutine grad_case(error, xyz, basis_json, spherical, charge, mult)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: xyz, basis_json
      logical, intent(in) :: spherical
      integer, intent(in) :: charge, mult
      type(pair_t) :: p
      type(trc_linalg_t) :: lp
      real(dp), allocatable :: graw(:, :), gsym(:, :), g1(:, :)
      real(dp) :: forbidden, dg, asym, after
      integer :: a, k, nat, ic
      real(dp) :: sg

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

      ! Cs: the mirror normal is z. An atom on the plane (its own image) has no
      ! z force; for the others, force(image) = s * force(atom), s = (1, 1, -1).
      forbidden = 0.0_dp
      after = 0.0_dp
      asym = 0.0_dp
      do a = 1, nat
         do k = 2, p%geom%group%order
            if (p%geom%perm(a, k) == a) then
               forbidden = max(forbidden, abs(graw(3, a)))
               after = max(after, abs(gsym(3, a)))
            end if
            do ic = 1, 3
               sg = merge(-1.0_dp, 1.0_dp, btest(p%geom%group%ops(k), ic - 1))
               asym = max(asym, abs(gsym(ic, p%geom%perm(a, k)) - sg*gsym(ic, a)))
            end do
         end do
      end do
      dg = maxval(abs(gsym - g1))
      print '(a,a,a,es10.2,a,es10.2,a,es10.2,a,es10.2)', '  ', xyz, ': forbidden raw ', forbidden, '  after symmetrisation ', after, &
         '  equivariance ', asym, '  |sym - C1| ', dg
      call check(error, forbidden <= TOL_GRAD, xyz//": forbidden gradient component "//to_string(forbidden, "es10.2"))
      if (allocated(error)) return
      call check(error, asym <= 1.0e-14_dp, xyz//": symmetrised gradient is not equivariant: "//to_string(asym, "es10.2"))
      if (allocated(error)) return
      call check(error, dg <= 1.0e-9_dp, xyz//": symmetrised gradient differs from C1 by "//to_string(dg, "es10.2"))
   end subroutine grad_case

   subroutine test_fixed_occupations(error)
      !! pyscf's irrep_nelec: the aufbau counts reproduce the aufbau energy; a
      !! pair moved from A' to A" gives a higher, converged state with the
      !! requested labels in the density's orbitals; a wrong sum is refused.
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      type(trc_sym_scf_t) :: sc
      type(trc_scf_result_t) :: rf
      type(trc_scf_options_t) :: opts
      integer :: m, nocc, cnt(2)

      call make_pair(error, 'hocl.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, with_c1=.false.)
      if (allocated(error)) return
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      opts%max_iteration = 200
      nocc = p%nalpha
      cnt = 0
      do m = 1, nocc
         cnt(p%rs%mo_irrep(m, 1)) = cnt(p%rs%mo_irrep(m, 1)) + 1
      end do

      ! aufbau counts, fixed
      sc%ao = p%ao
      allocate (sc%irrep_nocc(2, 1))
      sc%irrep_nocc(:, 1) = cnt
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rf, sym=sc)
      call check(error, rf%converged, "fixed aufbau occupations: not converged: "//trim(rf%message))
      if (.not. allocated(error)) call check(error, abs(rf%energy - p%rs%energy) <= TOL_E_HF, &
         "fixed aufbau occupations: energy differs from aufbau by "//to_string(abs(rf%energy - p%rs%energy), "es10.2"))
      if (allocated(error)) then
         call p%bas%release()
         return
      end if

      ! a pair from A' to A"
      sc%irrep_nocc(:, 1) = [cnt(1) - 1, cnt(2) + 1]
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
      sc%irrep_nocc(:, 1) = [cnt(1), cnt(2) + 1]
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rf, sym=sc)
      call check(error, .not. rf%converged .and. index(rf%message, "irrep_nocc") > 0, &
                 "a wrong irrep_nocc sum must be refused, got: "//trim(rf%message))
      call p%bas%release()
   end subroutine test_fixed_occupations

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(pair_t) :: p
      type(trc_sym_scf_t) :: sc
      type(trc_scf_result_t) :: rr, r0
      type(trc_scf_options_t) :: opts

      call make_pair(error, 'hocl.xyz', 'basis_sets/6-31gs.json', .false., '', 0, 1, p, with_c1=.true.)
      if (allocated(error)) return

      ! --soscf with symmetry: refused with a message, before any work
      opts%soscf = .true.
      sc%ao = p%ao
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rr, sym=sc)
      call check(error, .not. rr%converged .and. index(rr%message, "symmetry") > 0, &
                 "--soscf with symmetry must be refused, got: "//trim(rr%message))
      if (allocated(error)) then
         call p%bas%release()
         return
      end if

      ! a C1 group is the old path: no labels, and the same result as a run
      ! without `sym`. Not bitwise: under OpenMP (and on the device) the Fock
      ! build accumulates in a run-dependent order, so two plain C1 runs of
      ! this case differ by up to 4.5e-13 Eh and 3.1e-13 in D (measured at 4
      ! and 8 threads; 0 at 1 thread). The gate is ~20x that noise and far
      ! below any real symmetry effect (the Cs oracle cases above).
      opts%soscf = .false.
      opts%conv_energy = 1.0e-12_dp
      opts%conv_diis = 1.0e-10_dp
      sc%ao%nop = 1
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, rr, sym=sc)
      call trc_scf_run(p%bas, p%nalpha, p%nbeta, opts, r0)
      call check(error, rr%converged .and. r0%converged .and. .not. allocated(rr%mo_irrep), &
                 "a one-operation group must take the unchanged C1 path (no irrep labels)")
      if (allocated(error)) then
         call p%bas%release()
         return
      end if
      call check(error, abs(rr%energy - r0%energy) <= 1.0e-11_dp .and. &
                 maxval(abs(rr%dmat - r0%dmat)) <= 1.0e-11_dp, &
                 "a one-operation group must give the C1 result (to accumulation-order noise)")
      call p%bas%release()
   end subroutine test_refusals

end module check_sym_scf_cs_suite

program check_sym_scf_cs
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use check_sym_scf_cs_suite, only: collect_sym_scf_cs
   implicit none
   type(testsuite_type), allocatable :: testsuites(:)
   integer :: stat, is
   stat = 0
   testsuites = [new_testsuite("sym_scf_cs", collect_sym_scf_cs)]
   do is = 1, size(testsuites)
      write (error_unit, '(a)') "Testing: "//testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
   end do
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_scf_cs
