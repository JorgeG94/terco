!
! The single source of truth for mathematical and physical constants and
! atomic masses: every other module takes them from here, so a unit
! conversion can't silently differ between two parts of the code.
!
! Values carry their provenance: CODATA 2018 (the exact SI values where the
! 2019 redefinition fixed them), an exact definition, or derived from other
! constants here. HARTREE_TO_KJMOL is the published CODATA value; recomputing
! it from HARTREE_TO_J*NA_SI differs in the 13th digit, which is immaterial.
!
! BOHR_TO_ANGSTROM_CODATA2010 exists on purpose: PCM converts its cavity
! radii with pyscf's 2010 value so the radii agree with pyscf bit for bit.
!
module trc_constants
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   ! --- mathematical ----------------------------------------------------------
   public :: PI, TWO_PI, SQRT_PI

   ! --- CODATA 2018 (exact SI where the 2019 redefinition fixed them) ---------
   public :: BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR
   public :: HARTREE_TO_EV, HARTREE_TO_J, HARTREE_TO_KCALMOL, HARTREE_TO_KJMOL, HARTREE_TO_CM1
   public :: AMU_TO_KG, AMU_TO_ME
   public :: KB_SI, H_SI, HBAR_SI, C_SI, NA_SI, E_CHARGE_SI
   public :: KB_HARTREE, R_CALMOLK, CAL_TO_J, ATM_TO_PA, AU_TO_DEBYE
   public :: HARTREE_TO_CALMOL, CM1_TO_KELVIN
   public :: HARTREE_PER_AMU_BOHR2_TO_CM1, ROTTEMP_AMUA2_TO_K, ROTCONST_AMUA2_TO_GHZ

   ! --- legacy, kept on purpose -------------------------------------------------
   public :: BOHR_TO_ANGSTROM_CODATA2010

   ! --- atomic masses -----------------------------------------------------------
   public :: N_ELEMENTS, trc_atomic_mass, trc_atomic_mass_average

   real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
      !! Mathematical constant, exact; the most precise literal already in
      !! the tree (trc_dft_grid.F90, trc_dft_radial.F90), kept at full
      !! precision -- the extra digits beyond a double's ~17 significant
      !! figures round away, they do not change the stored value.
   real(dp), parameter :: TWO_PI = 2.0_dp*PI
      !! Derived. The translational and rotational partition functions
      !! (2*pi*m*kT/h^2, 8*pi^2*I*kT/h^2) use this form, not PI itself.
   real(dp), parameter :: SQRT_PI = sqrt(PI)
      !! Derived (SQRT is one of the intrinsics a PARAMETER initialization
      !! expression may fold at compile time). The nonlinear-rotor
      !! partition function (sqrt(pi)/sigma * ...) uses this directly;
      !! trc_pcm.F90 already computes it inline as `sqrt(PI)`.

   real(dp), parameter :: BOHR_TO_ANGSTROM = 0.529177210903_dp
      !! CODATA 2018 Bohr radius, in Angstrom.
   real(dp), parameter :: ANGSTROM_TO_BOHR = 1.0_dp/BOHR_TO_ANGSTROM
      !! Derived; app/terco_app.F90's BOHR_PER_ANGSTROM.
   real(dp), parameter :: HARTREE_TO_EV = 27.211386245988_dp
      !! CODATA 2018 Hartree energy, in eV.
   real(dp), parameter :: HARTREE_TO_J = 4.3597447222071e-18_dp
      !! CODATA 2018 Hartree energy, in J.
   real(dp), parameter :: HARTREE_TO_KCALMOL = 627.5094740631_dp
      !! Derived: HARTREE_TO_J * NA_SI / (1000 * CAL_TO_J), i.e. the
      !! thermochemical-calorie kcal/mol a Hartree converts to.
   real(dp), parameter :: HARTREE_TO_KJMOL = 2625.4996394799_dp
      !! Derived: HARTREE_TO_J * NA_SI / 1000.
   real(dp), parameter :: HARTREE_TO_CM1 = 219474.6313632_dp
      !! CODATA 2018 hartree-inverse meter relationship, converted to cm^-1.
   real(dp), parameter :: AMU_TO_KG = 1.66053906660e-27_dp
      !! CODATA 2018 atomic mass constant.
   real(dp), parameter :: AMU_TO_ME = 1822.888486209_dp
      !! CODATA 2018 atomic mass unit-electron mass ratio.
   real(dp), parameter :: KB_SI = 1.380649e-23_dp
      !! Exact SI (2019 redefinition): Boltzmann constant, J/K.
   real(dp), parameter :: H_SI = 6.62607015e-34_dp
      !! Exact SI (2019 redefinition): Planck constant, J s.
   real(dp), parameter :: HBAR_SI = H_SI/TWO_PI
      !! Derived: the reduced Planck constant, J s.
   real(dp), parameter :: C_SI = 299792458.0_dp
      !! Exact SI: speed of light in vacuum, m/s.
   real(dp), parameter :: NA_SI = 6.02214076e23_dp
      !! Exact SI (2019 redefinition): Avogadro constant, mol^-1.
   real(dp), parameter :: E_CHARGE_SI = 1.602176634e-19_dp
      !! Exact SI (2019 redefinition): elementary charge, C.
   real(dp), parameter :: KB_HARTREE = 3.1668115634556e-6_dp
      !! Derived: KB_SI / HARTREE_TO_J, the Boltzmann constant in Hartree/K.
   real(dp), parameter :: R_CALMOLK = 1.98720425864_dp
      !! Derived: the molar gas constant R = KB_SI * NA_SI, in the
      !! thermochemical cal/(mol K) (CAL_TO_J, not the 4.1868 "IT" calorie).
   real(dp), parameter :: CAL_TO_J = 4.184_dp
      !! Exact by definition: the thermochemical calorie, J.
   real(dp), parameter :: ATM_TO_PA = 101325.0_dp
      !! Exact by definition: the standard atmosphere, Pa.
   real(dp), parameter :: AU_TO_DEBYE = 2.541746473_dp
      !! CODATA 2018: the atomic unit of electric dipole moment (e*a0),
      !! in Debye.
   real(dp), parameter :: HARTREE_TO_CALMOL = HARTREE_TO_KCALMOL*1000.0_dp
      !! Derived: HARTREE_TO_KCALMOL, in cal/mol instead of kcal/mol -- the
      !! denominator trc_thermo.F90 uses to fold an entropy term (already
      !! in cal/(mol K)) back into a Gibbs correction in Hartree.
   real(dp), parameter :: CM1_TO_KELVIN = H_SI*C_SI*1.0e2_dp/KB_SI
      !! Derived: h*c/k_B, converted from m^-1 to cm^-1 (the 1e2) -- a
      !! vibrational wavenumber's characteristic temperature, theta_v =
      !! CM1_TO_KELVIN * nu_tilde. The CODATA 2014 value is
      !! 1.4387773538277; this one, computed from the CODATA 2018
      !! constants already above, is 1.438776877, 3.3e-7 relative lower -- negligible next to the SCF/gradient
      !! noise a finite-difference Hessian carries (trc_vib.F90's header).
   real(dp), parameter :: HARTREE_PER_AMU_BOHR2_TO_CM1 = &
      sqrt(HARTREE_TO_J/(AMU_TO_KG*(BOHR_TO_ANGSTROM*1.0e-10_dp)**2))/(TWO_PI*C_SI)*1.0e-2_dp
      !! Derived: sqrt(E_h / (amu * a0^2)) in cm^-1 -- the factor
      !! trc_vib.F90 multiplies sign(lambda)*sqrt(|lambda|) by to turn a
      !! mass-weighted-Hessian eigenvalue (Hartree/(amu bohr^2)) into a
      !! wavenumber. The 1e-10 is Angstrom -> metre (folded into
      !! BOHR_TO_ANGSTROM to get a0 in metres) and the 1e-2 is m^-1 ->
      !! cm^-1; both exact unit conversions, not measured constants. Same
      !! formula pyscf's hessian/thermo.py uses (`au2hz`, CODATA-2014-ish
      !! nist values there); the two agree to ~1e-9 relative (measured),
      !! which is why a finite-difference frequency test can hold a tight
      !! tolerance despite the different CODATA years.
   real(dp), parameter :: ROTTEMP_AMUA2_TO_K = H_SI**2/(8.0_dp*PI**2*AMU_TO_KG*(1.0e-10_dp)**2*KB_SI)
      !! Derived: h^2 / (8 pi^2 * amu * Angstrom^2 * k_B) -- the
      !! rotational temperature theta_r = ROTTEMP_AMUA2_TO_K / I, I in
      !! amu*Angstrom^2. Evaluates to 24.2544 (pyscf agrees); the 24.2637
      !! sometimes quoted for it is wrong in the 4th digit, which is why
      !! this is derived rather than typed in.
   real(dp), parameter :: ROTCONST_AMUA2_TO_GHZ = H_SI/(8.0_dp*PI**2*AMU_TO_KG*(1.0e-10_dp)**2)*1.0e-9_dp
      !! Derived: h / (8 pi^2 * amu * Angstrom^2), in GHz -- a rotational
      !! constant B = ROTCONST_AMUA2_TO_GHZ / I. Evaluates to 505.379;
      !! 505379.07 is the same coefficient in MHz, an easy one to mislabel.

   real(dp), parameter :: BOHR_TO_ANGSTROM_CODATA2010 = 0.52917721092_dp
      !! CODATA 2010, NOT 2018 -- kept so trc_pcm_radii.F90's cavity table,
      !! converted with pyscf's own `radii.BOHR`, agrees with pyscf to the
      !! last bit of the radius rather than to the eighth digit. See that
      !! file's header for why a continuum model cares.

   integer, parameter :: N_ELEMENTS = 118
      !! Highest atomic number either mass table below covers.

   real(dp), parameter :: ATOMIC_WEIGHT_TABLE(N_ELEMENTS) = [ &
      1.008_dp, 4.002602_dp, 6.94_dp, 9.0121831_dp, 10.81_dp, &
      12.011_dp, 14.007_dp, 15.999_dp, 18.998403163_dp, 20.1797_dp, &
      22.98976928_dp, 24.305_dp, 26.9815385_dp, 28.085_dp, 30.973761998_dp, &
      32.06_dp, 35.45_dp, 39.948_dp, 39.0983_dp, 40.078_dp, &
      44.955908_dp, 47.867_dp, 50.9415_dp, 51.9961_dp, 54.938044_dp, &
      55.845_dp, 58.933194_dp, 58.6934_dp, 63.546_dp, 65.38_dp, &
      69.723_dp, 72.63_dp, 74.921595_dp, 78.971_dp, 79.904_dp, &
      83.798_dp, 85.4678_dp, 87.62_dp, 88.90584_dp, 91.224_dp, &
      92.90637_dp, 95.95_dp, 97.90721_dp, 101.07_dp, 102.9055_dp, &
      106.42_dp, 107.8682_dp, 112.414_dp, 114.818_dp, 118.71_dp, &
      121.76_dp, 127.6_dp, 126.90447_dp, 131.293_dp, 132.90545196_dp, &
      137.327_dp, 138.90547_dp, 140.116_dp, 140.90766_dp, 144.242_dp, &
      144.91276_dp, 150.36_dp, 151.964_dp, 157.25_dp, 158.92535_dp, &
      162.5_dp, 164.93033_dp, 167.259_dp, 168.93422_dp, 173.054_dp, &
      174.9668_dp, 178.49_dp, 180.94788_dp, 183.84_dp, 186.207_dp, &
      190.23_dp, 192.217_dp, 195.084_dp, 196.966569_dp, 200.592_dp, &
      204.38_dp, 207.2_dp, 208.9804_dp, 208.98243_dp, 209.98715_dp, &
      222.01758_dp, 223.01974_dp, 226.02541_dp, 227.02775_dp, 232.0377_dp, &
      231.03588_dp, 238.02891_dp, 237.04817_dp, 244.06421_dp, 243.06138_dp, &
      247.07035_dp, 247.07031_dp, 251.07959_dp, 252.083_dp, 257.09511_dp, &
      258.09843_dp, 259.101_dp, 262.11_dp, 267.122_dp, 268.126_dp, &
      271.134_dp, 270.133_dp, 269.1338_dp, 278.156_dp, 281.165_dp, &
      281.166_dp, 285.177_dp, 286.182_dp, 289.19_dp, 289.194_dp, &
      293.204_dp, 293.208_dp, 294.214_dp]
      !! Isotope-averaged standard atomic weights, Z = 1..118, from pyscf
      !! 2.14.0's `pyscf.data.elements.MASSES` (its own IUPAC-sourced
      !! table). For the handful of elements with no stable isotope at
      !! all (Tc, Pm, Po and everything from At on), pyscf substitutes the
      !! most long-lived isotope's mass here too, rather than leaving a
      !! NaN -- see that array's own header comment.

   real(dp), parameter :: ISOTOPE_MASS_TABLE(N_ELEMENTS) = [ &
      1.007825_dp, 4.002603_dp, 7.016004_dp, 9.012182_dp, 11.009305_dp, &
      12.0_dp, 14.003074_dp, 15.994915_dp, 18.998403_dp, 19.99244_dp, &
      22.98977_dp, 23.985042_dp, 26.981538_dp, 27.976927_dp, 30.973762_dp, &
      31.972071_dp, 34.968853_dp, 39.962383_dp, 38.963707_dp, 39.962591_dp, &
      44.95591_dp, 47.947947_dp, 50.943964_dp, 51.940512_dp, 54.93805_dp, &
      55.934942_dp, 58.9332_dp, 57.935348_dp, 62.929601_dp, 63.929147_dp, &
      68.925581_dp, 73.921178_dp, 74.921596_dp, 79.916522_dp, 78.918338_dp, &
      83.911507_dp, 84.911789_dp, 87.905614_dp, 88.905848_dp, 89.904704_dp, &
      92.906378_dp, 97.905408_dp, 98.907216_dp, 101.90435_dp, 102.905504_dp, &
      105.903483_dp, 106.905093_dp, 113.903358_dp, 114.903878_dp, 119.902197_dp, &
      120.903818_dp, 129.906223_dp, 126.904468_dp, 131.904154_dp, 132.905447_dp, &
      137.905241_dp, 138.906348_dp, 139.905435_dp, 140.907648_dp, 141.907719_dp, &
      144.912744_dp, 151.919729_dp, 152.921227_dp, 157.924101_dp, 158.925343_dp, &
      163.929171_dp, 164.930319_dp, 165.93029_dp, 168.934211_dp, 173.938858_dp, &
      174.940768_dp, 179.946549_dp, 180.947996_dp, 183.950933_dp, 186.955751_dp, &
      191.961479_dp, 192.962924_dp, 194.964774_dp, 196.966552_dp, 201.970626_dp, &
      204.974412_dp, 207.976636_dp, 208.980383_dp, 208.982416_dp, 209.987131_dp, &
      222.01757_dp, 223.019731_dp, 226.025403_dp, 227.027747_dp, 232.03805_dp, &
      231.035879_dp, 238.050783_dp, 237.048167_dp, 244.064198_dp, 243.061373_dp, &
      247.070347_dp, 247.070299_dp, 251.07958_dp, 252.082972_dp, 257.095099_dp, &
      258.098425_dp, 259.101024_dp, 262.109692_dp, 267.0_dp, 268.0_dp, &
      269.0_dp, 270.0_dp, 270.0_dp, 278.0_dp, 281.0_dp, &
      282.0_dp, 285.0_dp, 286.0_dp, 289.0_dp, 290.0_dp, &
      293.0_dp, 294.0_dp, 294.0_dp]
      !! Mass of the most abundant (or, past bismuth, the longest-lived)
      !! isotope, Z = 1..118, from pyscf 2.14.0's
      !! `pyscf.data.elements.COMMON_ISOTOPE_MASSES` (its own copy of the
      !! NC State isotopic-mass table). Past rutherfordium (Z=104) pyscf
      !! has no measured value and fills in the mass number as a float
      !! (267.0, 268.0, ...); this table carries that placeholder too --
      !! fine for mass-weighting a Hessian on an element terco's basis
      !! sets do not reach anyway.

contains

   function trc_atomic_mass(z, error) result(mass)
      !! The most abundant isotope's mass, in amu. Use this to mass-weight
      !! a Cartesian Hessian the way a vibrational analysis conventionally
      !! does (the isotope terco's SCF itself never sees).
      !!
      !! Not PURE: a pure function's dummy arguments must all be
      !! INTENT(IN) (unlike a pure subroutine's, a function's result is
      !! supposed to be the only thing it hands back), so the out-of-range
      !! refusal below -- trc_error's own convention for a user-facing
      !! lookup, the same one src/grid's grid builders use -- rules it out.
      integer, intent(in) :: z
      type(error_t), intent(inout) :: error
      real(dp) :: mass
      if (z < 1 .or. z > N_ELEMENTS) then
         call error%set(ERROR_VALIDATION, "trc_constants: atomic number "//to_string(z)// &
                        " is out of range [1, "//to_string(N_ELEMENTS)//"]")
         mass = 0.0_dp
         return
      end if
      mass = ISOTOPE_MASS_TABLE(z)
   end function trc_atomic_mass

   function trc_atomic_mass_average(z, error) result(mass)
      !! The isotope-averaged standard atomic weight, in amu. See
      !! trc_atomic_mass's header for why this is not PURE either.
      integer, intent(in) :: z
      type(error_t), intent(inout) :: error
      real(dp) :: mass
      if (z < 1 .or. z > N_ELEMENTS) then
         call error%set(ERROR_VALIDATION, "trc_constants: atomic number "//to_string(z)// &
                        " is out of range [1, "//to_string(N_ELEMENTS)//"]")
         mass = 0.0_dp
         return
      end if
      mass = ATOMIC_WEIGHT_TABLE(z)
   end function trc_atomic_mass_average

end module trc_constants
