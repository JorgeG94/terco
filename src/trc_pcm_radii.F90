!
! Cavity radii for the polarisable continuum.
!
! A continuum model's answer depends on the cavity above almost anything
! else in it: the solvation energy goes roughly as 1/R, so a radius wrong by
! a few percent is a solvation energy wrong by a few percent, and no
! identity downstream catches it the way the electron count catches a wrong
! density. So the table is refused outside its range rather than
! extrapolated.
!
! The values are pyscf's `pcm.modified_Bondi`, element for element: Bondi's
! van der Waals radii [J. Phys. Chem. 68, 441 (1964)], with Mantina et al.
! [J. Phys. Chem. A 113, 5806 (2009)] for the elements Bondi did not cover,
! and hydrogen at 1.1 Angstrom rather than Bondi's 1.2 (the "modified"
! part). They are van der Waals radii, not the Bragg-Slater ones that size
! the DFT partition cells in src/grid: those are systematically smaller and
! would shrink every cavity.
!
! Hydrogen through calcium only. pyscf's table continues, but past calcium
! its entries are either a 2.0 Angstrom placeholder or, from copper on, one
! element out of step with Bondi's own (its copper is Bondi's nickel), and
! neither is something to reproduce on purpose.
!
! Converted to bohr with pyscf's own constant (0.52917721092 Angstrom), not
! CODATA 2018's, so a cavity here and pyscf's agree to the last bit of the
! radius rather than to the eighth digit.
!
module trc_pcm_radii
   use trc_boys, only: dp
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: pcm_cavity_radius
   public :: PCM_MAX_ELEMENT, PCM_DEFAULT_RADII_SCALE

   integer, parameter :: PCM_MAX_ELEMENT = 20
      !! Highest atomic number with a radius here
   real(dp), parameter :: PCM_DEFAULT_RADII_SCALE = 1.2_dp
      !! van der Waals radius to cavity radius; pyscf's and Q-Chem's default
   real(dp), parameter :: BOHR_ANGSTROM = 0.52917721092_dp
      !! pyscf's `radii.BOHR`, deliberately; see the header
   real(dp), parameter :: VDW_ANGSTROM(PCM_MAX_ELEMENT) = [ &
                          1.10_dp, 1.40_dp, &                                              ! H  He
                          1.82_dp, 1.53_dp, 1.92_dp, 1.70_dp, 1.55_dp, 1.52_dp, 1.47_dp, 1.54_dp, & ! Li..Ne
                          2.27_dp, 1.73_dp, 1.84_dp, 2.10_dp, 1.80_dp, 1.80_dp, 1.75_dp, 1.88_dp, & ! Na..Ar
                          2.75_dp, 2.31_dp]                                                ! K  Ca
      !! van der Waals radii in Angstrom, indexed by atomic number

contains

   subroutine pcm_cavity_radius(atomic_number, scale, radius, error)
      !! The scaled cavity radius of one element, in bohr
      integer, intent(in) :: atomic_number
      real(dp), intent(in) :: scale         !! usually PCM_DEFAULT_RADII_SCALE
      real(dp), intent(out) :: radius
      type(error_t), intent(inout) :: error
      character(len=16) :: num

      radius = 0.0_dp
      if (atomic_number < 1 .or. atomic_number > PCM_MAX_ELEMENT) then
         write (num, '(i0)') atomic_number
         call error%set(ERROR_VALIDATION, "PCM: no cavity radius for element "//trim(num)// &
                        "; the continuum carries van der Waals radii for hydrogen through "// &
                        "calcium only, and refuses rather than guesses one (a wrong radius "// &
                        "is a wrong solvation energy that nothing downstream would question)")
         return
      end if
      if (scale <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "PCM: the cavity radius scale must be positive")
         return
      end if
      radius = scale*VDW_ANGSTROM(atomic_number)/BOHR_ANGSTROM
   end subroutine pcm_cavity_radius

end module trc_pcm_radii
