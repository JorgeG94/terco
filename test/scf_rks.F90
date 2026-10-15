!
! The SCF driver end to end, against pyscf on the same grid.
!
! trc_scf is the half of terco's contract that used to be the caller's:
! basis and electron count in, converged energy and density out, every
! matrix resident on the device. This runs it on two molecules in 6-31G
! UNCONTRACTED to its primitives -- every primitive its own shell with
! coefficient one, because pyscf renormalises a contracted shell to unit
! self-overlap and terco does not, and a total energy compared at 1e-8
! cannot carry the 1e-6 that leaves in the basis -- and writes what pyscf
! needs to reproduce each case exactly: the shells as built, the grid, the
! energy and the final density.
!
!   water, closed shell:  RHF, then SVWN, PBE, BLYP, B3LYP, wB97X-V, wB97X,
!                         wB97 restricted
!   OH radical, doublet:  UHF, then PBE, wB97X-V, wB97X unrestricted
!
! So the restricted and unrestricted paths, Hartree-Fock and Kohn-Sham,
! LDA, GGA, global and range-separated hybrid, with and without VV10's
! nonlocal correlation, all go through the one driver. `rks_ref.py` rebuilds
! the same basis in the same order and runs pyscf's RHF/UHF/RKS/UKS on the
! same points from terco's converged density -- VV10's own (coarser) grid
! too, for a functional that has it, set as pyscf's `nlcgrids` exactly as
! the semilocal grid is set as `grids`.
!
program scf_rks
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz, uncontract_basis
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_int, ref_read_real, ref_read_str, &
                           ref_read_matrix3, ref_report
   implicit none

   integer, parameter :: NCASE = 12
   character(len=8), parameter :: mol(NCASE) = [character(len=8) :: &
                                                "water", "water", "water", "water", "water", "water", "water", &
                                                "water", "oh", "oh", "oh", "oh"]
   character(len=8), parameter :: fn(NCASE) = [character(len=8) :: &
                                               "", "svwn", "pbe", "blyp", "b3lyp", "wb97x-v", "wb97x", "wb97", &
                                               "", "pbe", "wb97x-v", "wb97x"]
   integer, parameter :: nalpha(NCASE) = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
   integer, parameter :: nbeta(NCASE) = [5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4]
   integer :: ic, unit, natm, nsh, maxnp, ish
   integer, allocatable :: sh_l(:), sh_np(:), zint(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), at_z(:), at_r(:, :)
   type(trc_basis_t) :: bas
   type(trc_scf_options_t) :: opts, opts2
   type(trc_scf_result_t) :: res, res2
   character(len=8) :: current
   logical :: all_ok
   integer :: refu, ncase_ref, nfail
   character(len=16) :: tag
   logical :: has_ref

   open (newunit=unit, file='rks_probe.bin', access='stream', form='unformatted', status='replace')
   write (unit) int(NCASE, kind=8)
   all_ok = .true.
   current = ""

   ! Opened (not required) before any case runs: a missing reference must
   ! not stop this program before rks_probe.bin is written, or
   ! test/pyscf_refs/generate.py has nothing to bootstrap from. The FAIL
   ! happens below, after the dump is safely written.
   call ref_open('rks.ref', refu, has_ref)
   nfail = 0
   if (has_ref) then
      call ref_read_int(refu, 'NCASE', ncase_ref)
      if (ncase_ref /= NCASE) then
         print '(a)', 'scf_rks: pyscf reference case count does not match this build'; stop 1
      end if
   end if
   do ic = 1, NCASE
      if (mol(ic) /= current) then
         if (len_trim(current) > 0) call bas%release()
         current = mol(ic)
         call read_xyz(trim(current)//'.xyz', natm, zint, at_r)
         if (allocated(at_z)) deallocate (at_z)
         allocate (at_z(natm))
         at_z = real(zint, dp)
         block
            type(trc_basis_t) :: jbas
            type(error_t) :: err
            call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, jbas, err)
            if (err%has_error()) then
               print '(a)', 'scf_rks: '//err%get_message(); stop 1
            end if
            call uncontract_basis(jbas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
         end block
         call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
         call bas%to_device()
         print '(a)', ""
         print '(a,a,a,i0,a,i0,a)', "scf_rks: ", trim(current), ", 6-31G uncontracted, ", nsh, " shells, ", bas%nao, " AOs"
      end if
      opts = trc_scf_options_t()
      opts%functional = fn(ic)
      opts%grid_level = 3
      opts%conv_energy = 1.0e-11_dp
      opts%conv_density = 1.0e-8_dp
      ! conv_density no longer gates (trc_scf_driver.F90); the commutator norm does.
      ! At its default 1e-5 the slow UKS wB97X-V case could stop up to ~2e-7
      ! above its minimum (seen under ifx), past the 1e-8 gate below, while
      ! terco's energy AT pyscf's density agreed to 2.5e-11. Converge to the
      ! gate's own precision instead.
      opts%conv_diis = 1.0e-8_dp
      call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts, res)
      if (.not. res%converged) then
         print '(a,a)', "  NOT CONVERGED: ", trim(res%message)
         all_ok = .false.
      end if
      if (len_trim(fn(ic)) > 0) then
         print '(a,a8,a,i1,a,i3,a,f22.12,a,f14.10)', "  ", fn(ic), " nspin ", res%nspin, " iter ", res%iterations, &
            "  E = ", res%energy, "  N(grid) = ", res%nelec_grid
      else
         print '(a,a8,a,i1,a,i3,a,f22.12)', "  ", "hf", " nspin ", res%nspin, " iter ", res%iterations, &
            "  E = ", res%energy
      end if

      ! Everything pyscf needs to reproduce this case.
      write (unit) int(bas%nao, kind=8), int(natm, kind=8), int(nsh, kind=8), int(res%nspin, kind=8), &
         int(nalpha(ic), kind=8), int(nbeta(ic), kind=8)
      write (unit) at_z, at_r
      do ish = 1, nsh
         write (unit) int(centre_of(ish), kind=8), int(sh_l(ish), kind=8), sh_e(1, ish)
      end do
      write (unit) fn(ic)
      if (len_trim(fn(ic)) > 0) then
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
      write (unit) res%energy
      write (unit) res%dmat

      ! --- against the committed pyscf reference (test/pyscf_refs/generate.py) ---
      if (has_ref) then
      block
         character(len=16) :: name_ref
         integer :: nspin_ref, nao_ref
         real(dp) :: e_ref, de, dat, de_tol, pyscf_conv_tol, pyscf_conv_tol_grad, pyscf_grad_norm
         real(dp), allocatable :: dm_ref(:, :, :)
         write (tag, '(i0)') ic - 1
         call ref_read_str(refu, 'NAME_'//trim(tag), name_ref)
         call ref_read_int(refu, 'NSPIN_'//trim(tag), nspin_ref)
         call ref_read_int(refu, 'NAO_'//trim(tag), nao_ref)
         call ref_read_real(refu, 'E_'//trim(tag), e_ref)
         ! Provenance only (not asserted on here): the pyscf convergence
         ! this energy was committed at -- see generate.py's converge_tight,
         ! which refuses a case whose gradient norm does not clear
         ! conv_tol_grad, tighter than pyscf's own sqrt(conv_tol) default.
         call ref_read_real(refu, 'CONV_TOL_'//trim(tag), pyscf_conv_tol)
         call ref_read_real(refu, 'CONV_TOL_GRAD_'//trim(tag), pyscf_conv_tol_grad)
         call ref_read_real(refu, 'GRAD_NORM_'//trim(tag), pyscf_grad_norm)
         if (nspin_ref /= res%nspin .or. nao_ref /= bas%nao) then
            print '(a)', 'scf_rks: pyscf reference shape does not match this build'; stop 1
         end if
         allocate (dm_ref(bas%nao, bas%nao, res%nspin))
         call ref_read_matrix3(refu, 'DM_'//trim(tag), bas%nao, bas%nao, res%nspin, dm_ref)
         de = abs(res%energy - e_ref)
         ! The at-density control: terco's OWN energy at pyscf's committed,
         ! converged density -- max_iter=1 with dguess runs exactly one
         ! Roothaan iteration, whose etot (src/scf/trc_scf_driver.F90's own "energy at
         ! this density" comment) is computed from dguess BEFORE the
         ! diagonalisation updates it, i.e. exactly energy_tot(dm0)'s
         ! equivalent -- never an iterated answer.
         opts2 = opts
         opts2%max_iter = 1
         call trc_scf_run(bas, nalpha(ic), nbeta(ic), opts2, res2, dguess=dm_ref)
         dat = abs(res2%energy - e_ref)
         print '(a,a,a,es10.2,a,es10.2)', '    ', trim(name_ref), &
            ': de(converged) ', de, '  de(at pyscf D) ', dat
         ! The CONVERGED energy of an open-shell DFT case gets a looser bound
         ! than everything else. OH's 2-Pi hole can sit anywhere in the
         ! degenerate pi plane, and with a spherically averaged SAD guess
         ! nothing picks the orientation, so round-off does -- which differs
         ! with the order of an OpenMP reduction. The grid is not invariant
         ! under an arbitrary rotation about the bond, so orientations that
         ! are equivalent exactly differ by grid error, measured up to
         ! 6.1e-7 (pbe) on the OpenMP build. The at-density check is what
         ! pins correctness, and it stays at 1e-8 for every case, as does the
         ! converged energy for HF (no grid) and every closed-shell case.
         if (res%nspin == 2 .and. len_trim(fn(ic)) > 0) then
            de_tol = 2.0e-6_dp
         else
            de_tol = 1.0e-8_dp
         end if
         call ref_report(nfail, de < de_tol .and. dat < 1.0e-8_dp, &
                         trim(name_ref)//': E(converged) and E(terco at pyscf D) vs pyscf', max(de, dat))
         deallocate (dm_ref)
      end block
      end if
   end do
   close (unit)
   call bas%release()
   if (has_ref) call ref_close(refu)
   if (.not. has_ref) then
      print '(a)', 'FAIL: missing reference pyscf_refs/rks.ref -- run test/pyscf_refs/generate.py ' &
         //'(see test/pyscf_refs/README.md)'
      stop 1
   end if
   if (.not. all_ok .or. nfail > 0) stop 1
   print '(a)', ""
   print '(a)', "scf_rks: all cases converged and agree with the pyscf reference; RESULT: PASS"

contains

   integer function centre_of(ish)
      integer, intent(in) :: ish
      integer :: a
      centre_of = 0
      do a = 1, natm
         if (all(abs(sh_r(:, ish) - at_r(:, a)) < 1.0e-12_dp)) centre_of = a
      end do
      if (centre_of == 0) then
         print '(a)', "scf_rks: a shell is on no atom"; stop 1
      end if
   end function centre_of

end program scf_rks
