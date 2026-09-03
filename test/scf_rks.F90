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
!   water, closed shell:  RHF, then SVWN, PBE, BLYP, B3LYP restricted
!   OH radical, doublet:  UHF, then PBE unrestricted
!
! So the restricted and unrestricted paths, Hartree-Fock and Kohn-Sham,
! LDA, GGA and hybrid, all go through the one driver. `rks_ref.py` rebuilds
! the same basis in the same order and runs pyscf's RHF/UHF/RKS/UKS on the
! same points from terco's converged density.
!
program scf_rks
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   integer, parameter :: NCASE = 7
   character(len=8), parameter :: mol(NCASE) = [character(len=8) :: &
                                                "water", "water", "water", "water", "water", "oh", "oh"]
   character(len=8), parameter :: fn(NCASE) = [character(len=8) :: &
                                               "", "svwn", "pbe", "blyp", "b3lyp", "", "pbe"]
   integer, parameter :: nalpha(NCASE) = [5, 5, 5, 5, 5, 5, 5]
   integer, parameter :: nbeta(NCASE) = [5, 5, 5, 5, 5, 4, 4]
   integer :: ic, unit, natm, nsh, maxnp, ish
   integer, allocatable :: sh_l(:), sh_np(:), zint(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), at_z(:), at_r(:, :)
   type(trc_basis_t) :: bas
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   character(len=8) :: current
   logical :: all_ok

   open (newunit=unit, file='rks_probe.bin', access='stream', form='unformatted', status='replace')
   write (unit) int(NCASE, kind=8)
   all_ok = .true.
   current = ""
   do ic = 1, NCASE
      if (mol(ic) /= current) then
         if (len_trim(current) > 0) call bas%release()
         current = mol(ic)
         call read_xyz(trim(current)//'.xyz', natm, zint, at_r)
         if (allocated(at_z)) deallocate (at_z)
         allocate (at_z(natm))
         at_z = real(zint, dp)
         call build_631g(natm, zint, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp, uncontracted=.true.)
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
      write (unit) res%energy
      write (unit) res%dmat
   end do
   close (unit)
   call bas%release()
   if (.not. all_ok) stop 1
   print '(a)', ""
   print '(a)', "scf_rks: all cases converged; run rks_ref.py for the pyscf comparison"

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
