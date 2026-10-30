!
! The SCF driver with a polarisable continuum, end to end, against pyscf.
!
! trc_scf_run with `pcm_model` set: the cavity, K and the nuclear potential
! once, then every iteration the surface potential of the current density,
! the charges, their operator in the Fock matrix and q.v/2 in the energy.
! Each case is compared two ways against test/pyscf_refs/pcm_scf*.ref
! (generate.py's gen_pcm_scf: pyscf's own RHF/UHF/RKS/UKS with .PCM(),
! same radii, scale, Lebedev order and dielectric, converged to conv_tol
! 1e-12 and conv_tol_grad 1e-8 on terco's own DFT grid):
!
!   converged:   terco's converged total energy and E_diel against pyscf's;
!   at-density:  terco's energy and E_diel AT pyscf's converged density --
!                trc_scf_run with max_iteration=1 and that density as dguess,
!                whose energy is evaluated at dguess before the Roothaan
!                step (the "energy at this density" in src/scf/trc_scf_driver.F90),
!                i.e. pyscf's energy_tot(dm) -- so a disagreement there is
!                the functional, not the convergence.
!
! The sets, chosen by the first argument:
!
!   main  water and formamide, 6-31G* Cartesian: HF, B3LYP and wB97X-V, each
!         C-PCM and IEF-PCM (restricted); the OH and NH2 radicals, UHF
!         C-PCM and B3LYP IEF-PCM (unrestricted); water/cc-pVDZ in its own
!         spherical AO space, HF C-PCM (converged energies only: pyscf's
!         spherical AO convention is not terco's, so its density is not
!         carried over)
!
! OH's converged run starts from pyscf's converged density rather than
! from GWH. Its pi pair is degenerate in vacuum and only split by the
! Lebedev cavity (and the DFT grid), so the energy along the rotation
! between the two singly occupied choices is flat to 3e-8 Eh (UHF) and
! 4e-7 Eh (B3LYP), with a stationary point at each end: gfortran's terco
! settled in one, ifx's in the other, from the same guess. Started inside
! pyscf's basin, the converged comparison still tests that terco's own
! iterations hold pyscf's solution to 1e-8; NH2, whose singly occupied
! orbital is not degenerate, is converged from GWH like every other case.
!   f     water/cc-pVTZ Cartesian, f on oxygen: HF C-PCM and B3LYP IEF-PCM;
!         skips (test-drive's skip_test) below TERCO_LMAX=3
!
! Writes pcm_scf_probe_<set>.bin first -- the geometry, the basis file, the
! grids terco used, its converged density to start pyscf from (as scf_rks
! does) -- which is all generate.py needs to rebuild each case; a missing
! reference fails only after that is safely written.
!
! One test-drive case per (molecule, functional, PCM model) in the chosen
! set; the shared setup -- the probe write and the sequential pyscf-
! reference read -- runs once -- see check_xc_energy_suite's own
! `next_ifn` comment for why, and why this suite stays `parallel=.false.`.
!
module scf_pcm_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_tables, only: LMAX
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3
   use golden_io, only: golden_begin, golden_put, golden_end
   implicit none
   private
   public :: collect_scf_pcm, set_arg

   real(dp), parameter :: EPS_WATER = 78.3553_dp
   real(dp), parameter :: TOL_E = 1.0e-8_dp
   integer, parameter :: NMAIN = 17, NF = 2
   character(len=16), parameter :: main_xyz(NMAIN) = [character(len=16) :: &
                                                      "water", "water", "water", "water", "water", "water", &
                                                      "formamide", "formamide", "formamide", "formamide", &
                                                      "formamide", "formamide", "oh", "oh", "nh2", "nh2", "water"]
   character(len=32), parameter :: main_bas(NMAIN) = [character(len=32) :: &
                                                      "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", &
                                                      "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", &
                                                      "6-31gs", "6-31gs", "6-31gs", "6-31gs", "cc-pvdz"]
   character(len=8), parameter :: main_fn(NMAIN) = [character(len=8) :: &
                                                    "", "", "b3lyp", "b3lyp", "wb97x-v", "wb97x-v", &
                                                    "", "", "b3lyp", "b3lyp", "wb97x-v", "wb97x-v", &
                                                    "", "b3lyp", "", "b3lyp", ""]
   character(len=8), parameter :: main_pcm(NMAIN) = [character(len=8) :: &
                                                     "cpcm", "iefpcm", "cpcm", "iefpcm", "cpcm", "iefpcm", &
                                                     "cpcm", "iefpcm", "cpcm", "iefpcm", "cpcm", "iefpcm", &
                                                     "cpcm", "iefpcm", "cpcm", "iefpcm", "cpcm"]
   integer, parameter :: main_na(NMAIN) = [5, 5, 5, 5, 5, 5, 12, 12, 12, 12, 12, 12, 5, 5, 5, 5, 5]
   integer, parameter :: main_nb(NMAIN) = [5, 5, 5, 5, 5, 5, 12, 12, 12, 12, 12, 12, 4, 4, 4, 4, 5]
   logical, parameter :: main_sph(NMAIN) = [.false., .false., .false., .false., .false., .false., &
                                            .false., .false., .false., .false., .false., .false., &
                                            .false., .false., .false., .false., .true.]
   logical, parameter :: main_seed(NMAIN) = [.false., .false., .false., .false., .false., .false., &
                                             .false., .false., .false., .false., .false., .false., &
                                             .true., .true., .false., .false., .false.]
   !! converge from pyscf's density (see the header: OH's degenerate pi pair)
   character(len=16), parameter :: f_xyz(NF) = [character(len=16) :: "water", "water"]
   character(len=32), parameter :: f_bas(NF) = [character(len=32) :: "cc-pvtz", "cc-pvtz"]
   character(len=8), parameter :: f_fn(NF) = [character(len=8) :: "", "b3lyp"]
   character(len=8), parameter :: f_pcm(NF) = [character(len=8) :: "cpcm", "iefpcm"]
   integer, parameter :: f_na(NF) = [5, 5], f_nb(NF) = [5, 5]
   logical, parameter :: f_sph(NF) = [.false., .false.]
   logical, parameter :: f_seed(NF) = [.false., .false.]

   character(len=16) :: set_arg = "main"

   logical :: computed = .false.
   logical :: skip_all = .false.
   integer :: ncase_module = 0
   logical, allocatable :: ok_case(:)
   character(len=256), allocatable :: msg_case(:)

contains

   subroutine collect_scf_pcm(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      integer :: ic, ncase
      character(len=16), allocatable :: cxyz(:)
      character(len=32), allocatable :: cbas(:)
      character(len=8), allocatable :: cfn(:), cpcm(:)

      select case (trim(set_arg))
      case ("main")
         ncase = NMAIN
         cxyz = main_xyz; cbas = main_bas; cfn = main_fn; cpcm = main_pcm
      case ("f")
         ncase = NF
         cxyz = f_xyz; cbas = f_bas; cfn = f_fn; cpcm = f_pcm
      case default
         error stop 'usage: scf_pcm [main|f]'
      end select

      allocate (testsuite(ncase))
      do ic = 1, ncase
         testsuite(ic) = new_unittest(trim(cxyz(ic))//"/"//trim(cbas(ic))//"/"// &
                                      trim(merge(cfn(ic), 'hf      ', len_trim(cfn(ic)) > 0))//"/"// &
                                      trim(cpcm(ic)), test_case)
      end do
   end subroutine collect_scf_pcm

   subroutine test_case(error)
      type(error_type), allocatable, intent(out) :: error
      integer, save :: next_ic = 0
      call ensure_computed()
      next_ic = next_ic + 1
      if (skip_all) then
         call skip_test(error, "the f set needs TERCO_LMAX >= 3")
         return
      end if
      call check(error, ok_case(next_ic), trim(msg_case(next_ic)))
   end subroutine test_case

   subroutine ensure_computed()
      if (computed) return
      call run_all()
      computed = .true.
   end subroutine ensure_computed

   subroutine run_all()
      character(len=16), allocatable :: cxyz(:)
      character(len=32), allocatable :: cbas(:)
      character(len=8), allocatable :: cfn(:), cpcm(:)
      integer, allocatable :: cna(:), cnb(:)
      logical, allocatable :: csph(:), cseed(:)
      integer :: ncase, ic, unit, natm, refu, ncase_ref, nspin, nao_scf
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      type(trc_basis_t) :: bas
      type(trc_scf_options_t) :: opts, opts2
      type(trc_scf_result_t) :: res, res2
      type(error_t) :: err
      logical :: has_ref
      integer :: nspin_ref, nao_ref
      real(dp) :: e_ref, ep_ref, e_dm, ep_dm, de, dat, dp_c, dp_at
      real(dp), allocatable :: dm_ref(:, :, :)
      character(len=16) :: tag
      character(len=48) :: what
      character(len=64) :: gkey

      select case (trim(set_arg))
      case ("main")
         ncase = NMAIN
         cxyz = main_xyz; cbas = main_bas; cfn = main_fn; cpcm = main_pcm
         cna = main_na; cnb = main_nb; csph = main_sph; cseed = main_seed
      case ("f")
         ncase = NF
         cxyz = f_xyz; cbas = f_bas; cfn = f_fn; cpcm = f_pcm
         cna = f_na; cnb = f_nb; csph = f_sph; cseed = f_seed
      end select
      ncase_module = ncase
      allocate (ok_case(ncase), msg_case(ncase))
      ok_case = .false.
      msg_case = "pyscf reference missing"

      call golden_begin('scf_pcm_'//trim(set_arg))

      ! The f set needs f in the four-centre kernels; below that there is
      ! nothing to run yet.
      if (trim(set_arg) == "f" .and. LMAX < 3) then
         print '(a,i0,a)', '  SKIP: the f set needs TERCO_LMAX >= 3; this build has LMAX = ', LMAX, '.'
         skip_all = .true.
         return
      end if

      open (newunit=unit, file='pcm_scf_probe_'//trim(set_arg)//'.bin', access='stream', &
            form='unformatted', status='replace')
      write (unit) int(ncase, int64)
      call ref_open('pcm_scf_'//trim(set_arg)//'.ref', refu, has_ref)
      if (has_ref) then
         call ref_read_int(refu, 'NCASE', ncase_ref)
         if (ncase_ref /= ncase) error stop 'scf_pcm: pyscf reference case count does not match this build'
      end if

      do ic = 1, ncase
         call read_xyz(trim(cxyz(ic))//'.xyz', natm, zint, at_r)
         call trc_basis_from_json('basis_sets/'//trim(cbas(ic))//'.json', natm, zint, at_r, bas, err)
         if (err%has_error()) error stop 'scf_pcm: '//err%get_message()
         bas%spherical = csph(ic)
         call bas%to_device()
         nspin = merge(2, 1, cna(ic) /= cnb(ic))
         nao_scf = merge(bas%nao_sph, bas%nao, csph(ic))

         ! The reference first: a seeded case starts from its density.
         if (has_ref) then
            write (tag, '(i0)') ic - 1
            call ref_read_int(refu, 'NSPIN_'//trim(tag), nspin_ref)
            call ref_read_int(refu, 'NAO_'//trim(tag), nao_ref)
            call ref_read_real(refu, 'E_'//trim(tag), e_ref)
            call ref_read_real(refu, 'EPCM_'//trim(tag), ep_ref)
            if (nspin_ref /= nspin .or. nao_ref /= nao_scf) &
               error stop 'scf_pcm: pyscf reference shape does not match this build'
            if (.not. csph(ic)) then
               call ref_read_real(refu, 'EDM_'//trim(tag), e_dm)
               call ref_read_real(refu, 'EPCMDM_'//trim(tag), ep_dm)
               if (allocated(dm_ref)) deallocate (dm_ref)
               allocate (dm_ref(bas%nao, bas%nao, nspin))
               call ref_read_matrix3(refu, 'DM_'//trim(tag), bas%nao, bas%nao, nspin, dm_ref)
            end if
         end if

         opts = trc_scf_options_t()
         opts%functional = cfn(ic)
         opts%grid_level = 3
         ! VV10's grid at its coarsest: pyscf is handed the same points, so the
         ! comparison is as exact on it as on any other, and VV10's host pass
         ! is a double sum over this grid -- level 1 put sixteen seconds an
         ! iteration on formamide, level 0 under two.
         opts%nlc_grid_level = 0
         opts%conv_energy = 1.0e-11_dp
         opts%conv_diis = 1.0e-8_dp
         opts%max_iteration = 200
         opts%pcm_model = cpcm(ic)
         opts%pcm_eps = EPS_WATER
         if (cseed(ic) .and. has_ref) then
            call trc_scf_run(bas, cna(ic), cnb(ic), opts, res, dguess=dm_ref)
         else
            call trc_scf_run(bas, cna(ic), cnb(ic), opts, res)
         end if
         if (.not. res%converged) print '(a,a)', "  NOT CONVERGED: ", trim(res%message)
         print '(a,a10,a8,a8,a,i1,a,i3,a,f20.12,a,f16.12,a,i0)', "  ", trim(cxyz(ic)), trim(cfn(ic)), &
            trim(cpcm(ic)), " nspin ", res%nspin, " iter ", res%iterations, "  E = ", res%energy, &
            "  E_diel = ", res%e_pcm, "  points ", res%pcm_npoints
         write (gkey, '(i0)') ic
         gkey = trim(gkey)//'_'//trim(cxyz(ic))//'_'//trim(cbas(ic))//'_'// &
                trim(merge(cfn(ic), 'hf      ', len_trim(cfn(ic)) > 0))//'_'//trim(cpcm(ic))
         call golden_put('E_'//trim(gkey), res%energy)
         call golden_put('E_DIEL_'//trim(gkey), res%e_pcm)

         write (unit) int(bas%nao, int64), int(merge(bas%nao_sph, bas%nao, csph(ic)), int64), &
            int(natm, int64), int(res%nspin, int64), int(cna(ic), int64), int(cnb(ic), int64), &
            int(merge(1, 0, csph(ic)), int64)
         write (unit) cxyz(ic), cbas(ic), cfn(ic), cpcm(ic)
         write (unit) EPS_WATER, int(opts%pcm_lebedev, int64), opts%pcm_radii_scale
         write (unit) int(zint, int64), at_r
         if (len_trim(cfn(ic)) > 0) then
            write (unit) int(size(res%grid_w), int64)
            write (unit) res%grid_r, res%grid_w
         else
            write (unit) 0_int64
         end if
         if (allocated(res%grid_w_nlc)) then
            write (unit) int(size(res%grid_w_nlc), int64)
            write (unit) res%grid_r_nlc, res%grid_w_nlc
         else
            write (unit) 0_int64
         end if
         ! The Cartesian overlap, for the per-AO normalisation the generator
         ! carries pyscf's density across with.
         block
            type(trc_pairlist_t) :: pl
            real(dp), allocatable :: s(:, :), t(:, :), v(:, :)
            allocate (s(bas%nao, bas%nao), t(bas%nao, bas%nao), v(bas%nao, bas%nao))
            call pl%build(bas, 1.0e-14_dp)
            call pl%to_device()
            !$acc enter data create(s, t, v)
            call trc_1e(bas, pl, s, t, v)
            !$acc update self(s)
            !$acc exit data delete(s, t, v)
            call pl%release()
            write (unit) s
         end block
         write (unit) res%energy, res%e_pcm
         ! The converged density, as scf_rks writes its own: generate.py starts
         ! pyscf from it, so both codes settle in the same basin (see the
         ! header on OH).
         write (unit) res%dmat

         if (has_ref) then
            what = trim(cxyz(ic))//'/'//trim(cbas(ic))//' '//trim(merge(cfn(ic), 'hf      ', &
                   len_trim(cfn(ic)) > 0))//' '//trim(cpcm(ic))
            if (cseed(ic)) what = trim(what)//' (from pyscf''s D)'
            de = abs(res%energy - e_ref)
            dp_c = abs(res%e_pcm - ep_ref)
            ok_case(ic) = res%converged .and. de < TOL_E .and. dp_c < TOL_E
            write (msg_case(ic), '(a,es10.2)') trim(what)//': E/E_diel converged vs pyscf, worst = ', max(de, dp_c)
            if (.not. csph(ic)) then
               opts2 = opts
               opts2%max_iteration = 1
               call trc_scf_run(bas, cna(ic), cnb(ic), opts2, res2, dguess=dm_ref)
               dat = abs(res2%energy - e_dm)
               dp_at = abs(res2%e_pcm - ep_dm)
               call golden_put('E_AT_PYSCF_D_'//trim(gkey), res2%energy)
               call golden_put('E_DIEL_AT_PYSCF_D_'//trim(gkey), res2%e_pcm)
               ok_case(ic) = ok_case(ic) .and. dat < TOL_E .and. dp_at < 1.0e-9_dp
               write (msg_case(ic), '(a,es10.2)') trim(what)// &
                  ': E/E_diel converged and at-density vs pyscf, worst = ', max(max(de, dp_c), max(dat, dp_at))
            end if
         end if
         call bas%release()
      end do
      call golden_end()
      close (unit)
      if (has_ref) call ref_close(refu)
      if (.not. has_ref) error stop 'FAIL: missing reference pyscf_refs/pcm_scf_'//trim(set_arg)//'.ref -- run ' &
         //'test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)'
   end subroutine run_all

end module scf_pcm_suite

program scf_pcm
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use scf_pcm_suite, only: collect_scf_pcm, set_arg
   implicit none
   integer :: stat
   if (command_argument_count() >= 1) call get_command_argument(1, set_arg)
   stat = 0
   call run_testsuite(collect_scf_pcm, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program scf_pcm
