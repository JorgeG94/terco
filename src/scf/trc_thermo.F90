!
! Ideal-gas rigid-rotor-harmonic-oscillator (RRHO) thermochemistry from
! trc_vib's own vibrational frequencies.
!
! Follows pyscf's hessian/thermo.py `thermo()` formulas and conventions
! (read, not guessed -- see the comments below for where each term maps),
! because that is this project's external reference, with three
! deliberate inputs pyscf's own convenience wrapper does not expose and
! mqc (HESSIAN_VIB_THERMO_SPEC.md, defects D1/D2) never passes at all:
!
!   * the rotational symmetry number sigma is a REQUIRED argument. pyscf's
!     own `thermo()` derives it from point-group detection
!     (`rotational_symmetry_number`); this module leaves that to the
!     caller (the CLI's own symmetry detection, or 1 if it has none),
!     exactly so a caller is never silently handed sigma=1 for a water
!     molecule (sigma=2) without choosing to be.
!   * the spin multiplicity is a REQUIRED argument, for the same reason
!     (mqc's defect D2: an open-shell reference always got S_elec=0).
!   * ONE vibrational mode set. trc_vib never returns a translation or
!     rotation mode in the first place (see that module's own header), so
!     "which frequencies enter which term" is pyscf's single rule --
!     every nu > 0 enters ZPE and the vibrational partition function,
!     nu <= 0 (imaginary) is skipped everywhere -- with no second,
!     inconsistent mode set to reconcile (mqc's defects D4/D6/D7).
!
! UNITS, deliberately pyscf's split, not SI throughout: energies (ZPE, the
! thermal corrections, the totals) in Hartree per molecule; entropies and
! heat capacities in cal/(mol K), the way HESSIAN_VIB_THERMO_SPEC.md's own
! reference numbers are reported. Converting between the two uses
! KB_HARTREE (the Boltzmann constant in Hartree/K -- pyscf's own `R_Eh` is
! numerically this SAME constant, kB/HARTREE2J, despite the name: pyscf
! reports its own S/Cv in Eh/K, per molecule, not cal/(mol K)) and
! R_CALMOLK (kB * Avogadro, in the thermochemical calorie) respectively;
! both from trc_constants, never recomputed here.
!
! NO low-frequency cutoff. pyscf's `thermo()` applies none beyond nu > 0
! (read its source: `idx = freq.real > 0`, nothing else) -- mqc's own
! "drop nu < 10 cm^-1 or u > 100" (HESSIAN_VIB_THERMO_SPEC.md section 3.2)
! is mqc-specific, not pyscf's, and is not reproduced here.
!
module trc_thermo
   use pic_types, only: dp, default_int
   use pic_lapack_interfaces, only: pic_syev
   use pic_strings, only: to_string
   use pic_logger, only: global_logger
   use trc_constants, only: PI, SQRT_PI, KB_HARTREE, R_CALMOLK, KB_SI, H_SI, AMU_TO_KG, &
                            BOHR_TO_ANGSTROM, CM1_TO_KELVIN, ROTTEMP_AMUA2_TO_K, ROTCONST_AMUA2_TO_GHZ, &
                            HARTREE_TO_CALMOL
   use trc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   implicit none
   private

   public :: trc_thermo_result_t
   public :: trc_thermochemistry
   public :: trc_thermo_log_summary

   real(dp), parameter :: DEFAULT_TEMPERATURE_K = 298.15_dp
   real(dp), parameter :: DEFAULT_PRESSURE_PA = 101325.0_dp
      !! 1 atm; trc_constants' own ATM_TO_PA, spelled out so a caller
      !! reading this module's signature sees the number without having
      !! to chase the constants module too.
   real(dp), parameter :: MOMENT_LINEAR_TOL_AMUA2 = 1.0e-6_dp
      !! Below this (amu*Angstrom^2), a principal moment of inertia is
      !! treated as zero -- mqc's own threshold (HESSIAN_VIB_THERMO_SPEC.md
      !! section 3.1), which pyscf's rotor-type test (rot_const > 1e8 GHz,
      !! the same statement inverted) agrees with to several more digits
      !! than this analysis needs.

   type :: trc_thermo_result_t
      real(dp) :: temperature_k = DEFAULT_TEMPERATURE_K
      real(dp) :: pressure_pa = DEFAULT_PRESSURE_PA
      integer :: sigma = 1
      integer :: mult = 1
      logical :: is_linear = .false.
      logical :: is_atom = .false.
      integer :: n_vib = 0
      integer :: n_imaginary = 0
      real(dp) :: molecular_mass_amu = 0.0_dp
      real(dp) :: moments_of_inertia_amua2(3) = 0.0_dp   !! ascending
      real(dp) :: rot_const_ghz(3) = 0.0_dp               !! 0 where the moment is ~0
      real(dp) :: zpe_hartree = 0.0_dp
      real(dp) :: e_elec_hartree = 0.0_dp
      real(dp) :: delta_e_hartree = 0.0_dp
      real(dp) :: delta_h_hartree = 0.0_dp
      real(dp) :: delta_g_hartree = 0.0_dp
      real(dp) :: e_zpe_hartree = 0.0_dp                  !! E_elec + ZPE
      real(dp) :: e_tot_hartree = 0.0_dp
      real(dp) :: h_tot_hartree = 0.0_dp
      real(dp) :: g_tot_hartree = 0.0_dp
      real(dp) :: s_trans_calmolk = 0.0_dp
      real(dp) :: s_rot_calmolk = 0.0_dp
      real(dp) :: s_vib_calmolk = 0.0_dp
      real(dp) :: s_elec_calmolk = 0.0_dp
      real(dp) :: s_tot_calmolk = 0.0_dp
      real(dp) :: cv_trans_calmolk = 0.0_dp
      real(dp) :: cv_rot_calmolk = 0.0_dp
      real(dp) :: cv_vib_calmolk = 0.0_dp
      real(dp) :: cv_elec_calmolk = 0.0_dp
      real(dp) :: cv_tot_calmolk = 0.0_dp
      real(dp) :: q_trans = 0.0_dp
      real(dp) :: q_rot = 0.0_dp
      real(dp) :: q_vib = 0.0_dp
   end type trc_thermo_result_t

contains

   !
   ! freq_vib_cm1: trc_vibrational_analysis's own `freq_cm1`, the true
   ! vibrational modes only (negative = imaginary). sigma and mult are
   ! mandatory (see this module's header); temperature_k/pressure_pa
   ! default to 298.15 K / 1 atm, the one pair this module DOES default,
   ! because unlike sigma/mult there is no silently-wrong value to guess.
   !
   subroutine trc_thermochemistry(mass, coords_bohr, freq_vib_cm1, e_elec_hartree, sigma, mult, &
                                   result, error, temperature_k, pressure_pa)
      real(dp), intent(in) :: mass(:)
      real(dp), intent(in) :: coords_bohr(:, :)
      real(dp), intent(in) :: freq_vib_cm1(:)
      real(dp), intent(in) :: e_elec_hartree
      integer, intent(in) :: sigma
      integer, intent(in) :: mult
      type(trc_thermo_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: temperature_k
      real(dp), intent(in), optional :: pressure_pa

      integer :: natm, k, info
      real(dp) :: temperature, pressure
      real(dp) :: mass_tot_amu, mass_tot_kg, q_trans
      real(dp) :: theta_r, theta1, theta2, theta3
      real(dp), allocatable :: theta_v(:), u(:), ev(:)
      integer :: n_vib_warm
      logical :: idx_mask_any

      natm = size(mass)
      if (natm < 1 .or. size(coords_bohr, 1) /= 3 .or. size(coords_bohr, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "trc_thermochemistry: mass/coords_bohr shape mismatch")
         return
      end if
      if (any(mass <= 0.0_dp)) then
         call error%set(ERROR_VALIDATION, "trc_thermochemistry: every mass must be positive")
         return
      end if
      if (sigma < 1) then
         call error%set(ERROR_VALIDATION, "trc_thermochemistry: sigma (rotational symmetry number) "// &
                        "must be >= 1, got "//to_string(sigma))
         return
      end if
      if (mult < 1) then
         call error%set(ERROR_VALIDATION, "trc_thermochemistry: mult (spin multiplicity) must be "// &
                        ">= 1, got "//to_string(mult))
         return
      end if

      temperature = DEFAULT_TEMPERATURE_K
      if (present(temperature_k)) temperature = temperature_k
      pressure = DEFAULT_PRESSURE_PA
      if (present(pressure_pa)) pressure = pressure_pa
      if (temperature <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "trc_thermochemistry: temperature_k must be > 0")
         return
      end if
      if (pressure <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "trc_thermochemistry: pressure_pa must be > 0")
         return
      end if

      result%temperature_k = temperature
      result%pressure_pa = pressure
      result%sigma = sigma
      result%mult = mult
      result%e_elec_hartree = e_elec_hartree
      result%n_vib = size(freq_vib_cm1)
      result%n_imaginary = count(freq_vib_cm1 <= 0.0_dp)   ! among the vibrations only -- mqc's D4

      mass_tot_amu = sum(mass)
      result%molecular_mass_amu = mass_tot_amu

      call trc_moments_of_inertia(mass, coords_bohr, result%moments_of_inertia_amua2, info)
      if (info /= 0) then
         call error%set(ERROR_GENERIC, "trc_thermochemistry: the inertia-tensor eigh did not converge")
         return
      end if
      result%is_atom = (natm == 1) .or. (result%moments_of_inertia_amua2(3) < MOMENT_LINEAR_TOL_AMUA2)
      result%is_linear = (.not. result%is_atom) .and. natm >= 2 .and. &
                          (result%moments_of_inertia_amua2(1) < MOMENT_LINEAR_TOL_AMUA2)
      do k = 1, 3
         if (result%moments_of_inertia_amua2(k) > MOMENT_LINEAR_TOL_AMUA2) then
            result%rot_const_ghz(k) = ROTCONST_AMUA2_TO_GHZ/result%moments_of_inertia_amua2(k)
         end if
      end do

      ! --- electronic --------------------------------------------------------
      result%s_elec_calmolk = R_CALMOLK*log(real(mult, dp))
      result%cv_elec_calmolk = 0.0_dp

      ! --- translational (Sackur-Tetrode) -------------------------------------
      mass_tot_kg = mass_tot_amu*AMU_TO_KG
      q_trans = (2.0_dp*PI*mass_tot_kg*KB_SI*temperature/H_SI**2)**1.5_dp*KB_SI*temperature/pressure
      result%q_trans = q_trans
      result%s_trans_calmolk = R_CALMOLK*(2.5_dp + log(q_trans))
      result%cv_trans_calmolk = 1.5_dp*R_CALMOLK

      ! --- rotational ----------------------------------------------------------
      if (result%is_atom) then
         result%q_rot = 0.0_dp
         result%s_rot_calmolk = 0.0_dp
         result%cv_rot_calmolk = 0.0_dp
      else if (result%is_linear) then
         theta_r = ROTTEMP_AMUA2_TO_K/result%moments_of_inertia_amua2(3)
         result%q_rot = temperature/(real(sigma, dp)*theta_r)
         result%s_rot_calmolk = R_CALMOLK*(1.0_dp + log(result%q_rot))
         result%cv_rot_calmolk = R_CALMOLK
      else
         theta1 = ROTTEMP_AMUA2_TO_K/result%moments_of_inertia_amua2(1)
         theta2 = ROTTEMP_AMUA2_TO_K/result%moments_of_inertia_amua2(2)
         theta3 = ROTTEMP_AMUA2_TO_K/result%moments_of_inertia_amua2(3)
         result%q_rot = SQRT_PI*temperature**1.5_dp/(real(sigma, dp)*sqrt(theta1*theta2*theta3))
         result%s_rot_calmolk = R_CALMOLK*(1.5_dp + log(result%q_rot))
         result%cv_rot_calmolk = 1.5_dp*R_CALMOLK
      end if

      ! --- vibrational (every true mode with nu > 0; no low-frequency cutoff) --
      n_vib_warm = count(freq_vib_cm1 > 0.0_dp)
      allocate (theta_v(n_vib_warm), u(n_vib_warm), ev(n_vib_warm))
      theta_v = pack(CM1_TO_KELVIN*freq_vib_cm1, freq_vib_cm1 > 0.0_dp)
      u = theta_v/temperature
      ev = exp(-u)
      idx_mask_any = n_vib_warm > 0

      result%zpe_hartree = 0.5_dp*KB_HARTREE*sum(theta_v)
      if (idx_mask_any) then
         result%s_vib_calmolk = R_CALMOLK*sum(u*ev/(1.0_dp - ev) - log(1.0_dp - ev))
         result%cv_vib_calmolk = R_CALMOLK*sum(ev*u**2/(1.0_dp - ev)**2)
         result%q_vib = product(1.0_dp/(1.0_dp - ev))
      else
         result%s_vib_calmolk = 0.0_dp
         result%cv_vib_calmolk = 0.0_dp
         result%q_vib = 1.0_dp
      end if

      ! --- totals --------------------------------------------------------------
      block
         real(dp) :: e_trans, e_rot, e_vib_thermal
         e_trans = 1.5_dp*KB_HARTREE*temperature
         e_rot = merge(0.0_dp, merge(KB_HARTREE*temperature, 1.5_dp*KB_HARTREE*temperature, result%is_linear), &
                       result%is_atom)
         e_vib_thermal = merge(KB_HARTREE*sum(theta_v*ev/(1.0_dp - ev)), 0.0_dp, idx_mask_any)

         result%delta_e_hartree = result%zpe_hartree + e_trans + e_rot + e_vib_thermal
      end block
      result%delta_h_hartree = result%delta_e_hartree + KB_HARTREE*temperature
      result%s_tot_calmolk = result%s_trans_calmolk + result%s_rot_calmolk + result%s_vib_calmolk + &
                              result%s_elec_calmolk
      result%delta_g_hartree = result%delta_h_hartree - temperature*result%s_tot_calmolk/HARTREE_TO_CALMOL
      result%cv_tot_calmolk = result%cv_trans_calmolk + result%cv_rot_calmolk + result%cv_vib_calmolk + &
                              result%cv_elec_calmolk

      result%e_zpe_hartree = e_elec_hartree + result%zpe_hartree
      result%e_tot_hartree = e_elec_hartree + result%delta_e_hartree
      result%h_tot_hartree = e_elec_hartree + result%delta_h_hartree
      result%g_tot_hartree = e_elec_hartree + result%delta_g_hartree
   end subroutine trc_thermochemistry

   !
   ! Principal moments of inertia, ascending, in amu*Angstrom^2 -- the same
   ! convention pyscf's `rotation_const` builds from before converting to
   ! GHz, computed independently of trc_vib's own TR-generator SVD (this
   ! is pyscf's `thermo()` own recipe: moments of inertia first, rotor
   ! type from THOSE, never reusing harmonic_analysis's projector).
   !
   subroutine trc_moments_of_inertia(mass, coords_bohr, moments_amua2, info)
      real(dp), intent(in) :: mass(:)
      real(dp), intent(in) :: coords_bohr(:, :)
      real(dp), intent(out) :: moments_amua2(3)
      integer, intent(out) :: info
      integer :: natm, a, r, s
      real(dp) :: com(3), rvec(3), im(3, 3)
      integer(default_int) :: linfo

      natm = size(mass)
      com = 0.0_dp
      do a = 1, natm
         com = com + mass(a)*coords_bohr(:, a)*BOHR_TO_ANGSTROM
      end do
      com = com/sum(mass)

      im = 0.0_dp
      do a = 1, natm
         rvec = coords_bohr(:, a)*BOHR_TO_ANGSTROM - com
         do s = 1, 3
            do r = 1, 3
               im(r, s) = im(r, s) - mass(a)*rvec(r)*rvec(s)
            end do
         end do
         im(1, 1) = im(1, 1) + mass(a)*sum(rvec**2)
         im(2, 2) = im(2, 2) + mass(a)*sum(rvec**2)
         im(3, 3) = im(3, 3) + mass(a)*sum(rvec**2)
      end do

      call pic_syev(im, moments_amua2, jobz="V", uplo="U", info=linfo)
      info = int(linfo)
      where (abs(moments_amua2) < 1.0e-9_dp) moments_amua2 = 0.0_dp
   end subroutine trc_moments_of_inertia

   !
   ! The thermochemistry summary through global_logger, for the CLI's
   ! later --hessian report.
   !
   subroutine trc_thermo_log_summary(r)
      type(trc_thermo_result_t), intent(in) :: r
      character(len=256) :: line

      write (line, "(a,f7.2,a,f10.1,a,i0,a,i0)") "trc_thermo: T=", r%temperature_k, " K  P=", &
         r%pressure_pa, " Pa  sigma=", r%sigma, "  mult=", r%mult
      call global_logger%info(trim(line))
      write (line, "(a,i0,a,i0,a,l1)") "  n_vib=", r%n_vib, "  n_imaginary=", r%n_imaginary, &
         "  is_linear=", r%is_linear
      call global_logger%info(trim(line))
      write (line, "(a,f16.8,a)") "  ZPE            = ", r%zpe_hartree, " Eh"
      call global_logger%info(trim(line))
      write (line, "(a,f16.8,a,f16.8,a,f16.8)") "  dE= ", r%delta_e_hartree, "  dH= ", r%delta_h_hartree, &
         "  dG= ", r%delta_g_hartree
      call global_logger%info(trim(line))
      write (line, "(a,f16.8,a,f16.8,a,f16.8)") "  E_tot= ", r%e_tot_hartree, "  H_tot= ", r%h_tot_hartree, &
         "  G_tot= ", r%g_tot_hartree
      call global_logger%info(trim(line))
      write (line, "(a,f12.6,a,f12.6,a,f12.6,a,f12.6)") "  S trans/rot/vib/elec (cal/mol/K) = ", &
         r%s_trans_calmolk, " / ", r%s_rot_calmolk, " / ", r%s_vib_calmolk, " / ", r%s_elec_calmolk
      call global_logger%info(trim(line))
   end subroutine trc_thermo_log_summary

end module trc_thermo
