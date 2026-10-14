!
! The SCF driver with a polarisable continuum, end to end, against pyscf.
!
! trc_scf_run with `pcm_model` set: the cavity, K and the nuclear potential
! once, then every iteration the surface potential of the current density,
! the charges, their operator in the Fock matrix and q.v/2 in the energy.
! Each case is compared two ways against test/pyscf_refs/pcm_scf*.ref
! (generate.py's gen_pcm_scf: pyscf's own RHF/UHF/RKS/UKS with .PCM(),
! same radii, scale, Lebedev order and dielectric, converged to conv_tol
! 1e-12 and conv_tol_grad 1e-9 on terco's own DFT grid):
!
!   converged:   terco's converged total energy and E_diel against pyscf's;
!   at-density:  terco's energy and E_diel AT pyscf's converged density --
!                trc_scf_run with max_iter=1 and that density as dguess,
!                whose energy is evaluated at dguess before the Roothaan
!                step (the "energy at this density" in src/trc_scf.F90),
!                i.e. pyscf's energy_tot(dm) -- so a disagreement there is
!                the functional, not the convergence.
!
! The sets, chosen by the first argument:
!
!   main  water and formamide, 6-31G* Cartesian: HF, B3LYP and wB97X-V, each
!         C-PCM and IEF-PCM (restricted); the OH radical UHF C-PCM and
!         B3LYP IEF-PCM (unrestricted); water/cc-pVDZ in its own spherical
!         AO space, HF C-PCM (converged energies only: pyscf's spherical
!         AO convention is not terco's, so its density is not carried over)
!   f     water/cc-pVTZ Cartesian, f on oxygen: HF C-PCM and B3LYP IEF-PCM;
!         skips (SKIP, exit 0) below TERCO_LMAX=3
!
! Writes pcm_scf_probe_<set>.bin first -- the geometry, the basis file, the
! grids terco used, its converged density to start pyscf from (as scf_rks
! does) -- which is all generate.py needs to rebuild each case; a missing
! reference fails only after that is safely written.
!
program scf_pcm
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e
   use trc_tables, only: LMAX
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3, ref_report
   implicit none

   real(dp), parameter :: EPS_WATER = 78.3553_dp
   real(dp), parameter :: TOL_E = 1.0e-8_dp
   integer, parameter :: NMAIN = 15, NF = 2
   character(len=16), parameter :: main_xyz(NMAIN) = [character(len=16) :: &
                                                      "water", "water", "water", "water", "water", "water", &
                                                      "formamide", "formamide", "formamide", "formamide", &
                                                      "formamide", "formamide", "oh", "oh", "water"]
   character(len=32), parameter :: main_bas(NMAIN) = [character(len=32) :: &
                                                      "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", &
                                                      "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", "6-31gs", &
                                                      "6-31gs", "6-31gs", "cc-pvdz"]
   character(len=8), parameter :: main_fn(NMAIN) = [character(len=8) :: &
                                                    "", "", "b3lyp", "b3lyp", "wb97x-v", "wb97x-v", &
                                                    "", "", "b3lyp", "b3lyp", "wb97x-v", "wb97x-v", &
                                                    "", "b3lyp", ""]
   character(len=8), parameter :: main_pcm(NMAIN) = [character(len=8) :: &
                                                     "cpcm", "iefpcm", "cpcm", "iefpcm", "cpcm", "iefpcm", &
                                                     "cpcm", "iefpcm", "cpcm", "iefpcm", "cpcm", "iefpcm", &
                                                     "cpcm", "iefpcm", "cpcm"]
   integer, parameter :: main_na(NMAIN) = [5, 5, 5, 5, 5, 5, 12, 12, 12, 12, 12, 12, 5, 5, 5]
   integer, parameter :: main_nb(NMAIN) = [5, 5, 5, 5, 5, 5, 12, 12, 12, 12, 12, 12, 4, 4, 5]
   logical, parameter :: main_sph(NMAIN) = [.false., .false., .false., .false., .false., .false., &
                                            .false., .false., .false., .false., .false., .false., &
                                            .false., .false., .true.]
   character(len=16), parameter :: f_xyz(NF) = [character(len=16) :: "water", "water"]
   character(len=32), parameter :: f_bas(NF) = [character(len=32) :: "cc-pvtz", "cc-pvtz"]
   character(len=8), parameter :: f_fn(NF) = [character(len=8) :: "", "b3lyp"]
   character(len=8), parameter :: f_pcm(NF) = [character(len=8) :: "cpcm", "iefpcm"]
   integer, parameter :: f_na(NF) = [5, 5], f_nb(NF) = [5, 5]
   logical, parameter :: f_sph(NF) = [.false., .false.]

   character(len=16) :: set, tag
   character(len=16), allocatable :: cxyz(:)
   character(len=32), allocatable :: cbas(:)
   character(len=8), allocatable :: cfn(:), cpcm(:)
   integer, allocatable :: cna(:), cnb(:)
   logical, allocatable :: csph(:)
   integer :: ncase, ic, unit, natm, refu, nfail, ncase_ref
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   type(trc_basis_t) :: bas
   type(trc_scf_options_t) :: opts, opts2
   type(trc_scf_result_t) :: res, res2
   type(error_t) :: err
   logical :: has_ref, all_ok

   set = "main"
   if (command_argument_count() >= 1) call get_command_argument(1, set)
   select case (trim(set))
   case ("main")
      ncase = NMAIN
      cxyz = main_xyz; cbas = main_bas; cfn = main_fn; cpcm = main_pcm
      cna = main_na; cnb = main_nb; csph = main_sph
   case ("f")
      ncase = NF
      cxyz = f_xyz; cbas = f_bas; cfn = f_fn; cpcm = f_pcm
      cna = f_na; cnb = f_nb; csph = f_sph
   case default
      print '(a)', 'usage: scf_pcm [main|f]'; stop 1
   end select

   ! The f set needs f in the four-centre kernels; below that there is
   ! nothing to run yet.
   if (trim(set) == "f" .and. LMAX < 3) then
      print '(a,i0,a)', '  SKIP: the f set needs TERCO_LMAX >= 3; this build has LMAX = ', LMAX, '.'
      stop 0
   end if

   open (newunit=unit, file='pcm_scf_probe_'//trim(set)//'.bin', access='stream', &
         form='unformatted', status='replace')
   write (unit) int(ncase, kind=8)
   call ref_open('pcm_scf_'//trim(set)//'.ref', refu, has_ref)
   if (has_ref) then
      call ref_read_int(refu, 'NCASE', ncase_ref)
      if (ncase_ref /= ncase) then
         print '(a)', 'scf_pcm: pyscf reference case count does not match this program'; stop 1
      end if
   end if
   nfail = 0
   all_ok = .true.

   do ic = 1, ncase
      call read_xyz(trim(cxyz(ic))//'.xyz', natm, zint, at_r)
      call trc_basis_from_json('basis_sets/'//trim(cbas(ic))//'.json', natm, zint, at_r, bas, err)
      if (err%has_error()) then
         print '(a)', 'scf_pcm: '//err%get_message(); stop 1
      end if
      bas%spherical = csph(ic)
      call bas%to_device()

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
      opts%max_iter = 200
      opts%pcm_model = cpcm(ic)
      opts%pcm_eps = EPS_WATER
      call trc_scf_run(bas, cna(ic), cnb(ic), opts, res)
      if (.not. res%converged) then
         print '(a,a)', "  NOT CONVERGED: ", trim(res%message)
         all_ok = .false.
      end if
      print '(a,a10,a8,a8,a,i1,a,i3,a,f20.12,a,f16.12,a,i0)', "  ", trim(cxyz(ic)), trim(cfn(ic)), &
         trim(cpcm(ic)), " nspin ", res%nspin, " iter ", res%iterations, "  E = ", res%energy, &
         "  E_diel = ", res%e_pcm, "  points ", res%pcm_npoints

      write (unit) int(bas%nao, kind=8), int(merge(bas%nao_sph, bas%nao, csph(ic)), kind=8), &
         int(natm, kind=8), int(res%nspin, kind=8), int(cna(ic), kind=8), int(cnb(ic), kind=8), &
         int(merge(1, 0, csph(ic)), kind=8)
      write (unit) cxyz(ic), cbas(ic), cfn(ic), cpcm(ic)
      write (unit) EPS_WATER, int(opts%pcm_lebedev, kind=8), opts%pcm_radii_scale
      write (unit) int(zint, kind=8), at_r
      if (len_trim(cfn(ic)) > 0) then
         write (unit) int(size(res%grid_w), kind=8)
         write (unit) res%grid_r, res%grid_w
      else
         write (unit) 0_8
      end if
      if (allocated(res%grid_w_nlc)) then
         write (unit) int(size(res%grid_w_nlc), kind=8)
         write (unit) res%grid_r_nlc, res%grid_w_nlc
      else
         write (unit) 0_8
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
      ! pyscf from it, so both codes settle in the same basin. The OH
      ! radical needs that -- its pi pair is degenerate in a Lebedev cavity
      ! as in vacuum, and the rotation between the two singly occupied
      ! choices is so flat that UHF finds stationary points 3e-8 Eh apart.
      write (unit) res%dmat

      if (has_ref) then
         block
            integer :: nspin_ref, nao_ref
            real(dp) :: e_ref, ep_ref, e_dm, ep_dm, de, dat, dp_c, dp_at
            real(dp), allocatable :: dm_ref(:, :, :)
            character(len=48) :: what
            write (tag, '(i0)') ic - 1
            call ref_read_int(refu, 'NSPIN_'//trim(tag), nspin_ref)
            call ref_read_int(refu, 'NAO_'//trim(tag), nao_ref)
            call ref_read_real(refu, 'E_'//trim(tag), e_ref)
            call ref_read_real(refu, 'EPCM_'//trim(tag), ep_ref)
            if (nspin_ref /= res%nspin .or. nao_ref /= merge(bas%nao_sph, bas%nao, csph(ic))) then
               print '(a)', 'scf_pcm: pyscf reference shape does not match this build'; stop 1
            end if
            what = trim(cxyz(ic))//'/'//trim(cbas(ic))//' '//trim(merge(cfn(ic), 'hf      ', &
                   len_trim(cfn(ic)) > 0))//' '//trim(cpcm(ic))
            de = abs(res%energy - e_ref)
            dp_c = abs(res%e_pcm - ep_ref)
            call ref_report(nfail, de < TOL_E, trim(what)//': E converged', de)
            call ref_report(nfail, dp_c < TOL_E, trim(what)//': E_diel converged', dp_c)
            if (.not. csph(ic)) then
               call ref_read_real(refu, 'EDM_'//trim(tag), e_dm)
               call ref_read_real(refu, 'EPCMDM_'//trim(tag), ep_dm)
               allocate (dm_ref(bas%nao, bas%nao, res%nspin))
               call ref_read_matrix3(refu, 'DM_'//trim(tag), bas%nao, bas%nao, res%nspin, dm_ref)
               opts2 = opts
               opts2%max_iter = 1
               call trc_scf_run(bas, cna(ic), cnb(ic), opts2, res2, dguess=dm_ref)
               dat = abs(res2%energy - e_dm)
               dp_at = abs(res2%e_pcm - ep_dm)
               call ref_report(nfail, dat < TOL_E, trim(what)//': E at pyscf''s D', dat)
               call ref_report(nfail, dp_at < 1.0e-9_dp, trim(what)//': E_diel at pyscf''s D', dp_at)
            end if
         end block
      end if
      call bas%release()
   end do
   close (unit)
   if (has_ref) call ref_close(refu)
   if (.not. has_ref) then
      print '(a)', 'FAIL: missing reference pyscf_refs/pcm_scf_'//trim(set)//'.ref -- run '// &
         'test/pyscf_refs/generate.py (see test/pyscf_refs/README.md)'
      stop 1
   end if
   if (.not. all_ok .or. nfail > 0) then
      print '(a,i0,a)', 'scf_pcm: ', nfail, ' comparison(s) FAILED'; stop 1
   end if
   print '(a)', "scf_pcm: every case converged and agrees with pyscf's PCM; RESULT: PASS"

end program scf_pcm
