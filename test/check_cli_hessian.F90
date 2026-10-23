!
! terco's --hessian, through the installed CLI binary (the SAME path
! app/main.F90 drives: trc_run_hessian/trc_hessian, then trc_vib/
! trc_thermo) against pyscf's own analytic RHF/STO-3G Hessian and
! harmonic_analysis/thermo already in test/pyscf_refs/vib_thermo_water.ref
! -- the same reference check_fd_hessian_vs_pyscf.F90 and
! check_vib_thermo_vs_pyscf.F90 use, at the SAME geometry
! (test/vib_thermo_water.xyz, Angstrom, copied verbatim from
! generate.py's gen_vib_thermo literal).
!
! --symmetry-number 2 is passed explicitly: water's rotational symmetry
! number (the .ref's own SIGMA) is 2, and the CLI's own default is 1
! (never auto-detected, per --hessian's own warning) -- without this the
! comparison is not wrong, it is a DIFFERENT, equally legitimate
! calculation (sigma=1), which is exactly the point of requiring the flag.
!
! TOLERANCES. Measured 2026-10-22 (gfortran 15 host build): the CLI's
! frequencies differ from pyscf's analytic ones by -0.061, +0.019, +0.016
! cm^-1 -- the h = 0.005 Bohr finite-difference error, the same numbers
! check_fd_hessian_vs_pyscf measures through the library; the electronic
! energy agrees to 4e-13 Eh, so reading the geometry in Angstrom adds
! nothing. FREQ_TOL is 2x the worst mode. ZPE, dE, dH, dG and G_tot carry
! that error through half the sum of the frequency errors (6.7e-8 Eh
! measured), hence ENERGY_TOL 1.5e-7.
program check_cli_hessian
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_char
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use pyscf_ref_io, only: ref_require, ref_close, ref_read_int, ref_read_real, ref_read_vector, &
                           ref_read_matrix, ref_read_matrix3, ref_report
   implicit none

   character(len=:), allocatable :: cli_exe, probe_json, probe_out
   character(len=1024) :: arg
   character(len=:), allocatable :: cmd
   integer :: rc, n_modes_cli, n_modes_ref, k, nfail
   real(dp) :: zpe_ref, de_ref, dh_ref, dg_ref, gtot_ref
   real(dp), allocatable :: freq_ref(:)
   real(dp) :: zpe_cli, de_cli, dh_cli, dg_cli, gtot_cli
   real(dp), allocatable :: freq_cli(:)
   type(json_file) :: json
   logical :: found
   real(dp), parameter :: ENERGY_TOL = 1.5e-7_dp
   real(dp), parameter :: FREQ_TOL = 0.12_dp

   nfail = 0

   call get_command_argument(1, arg)
   cli_exe = trim(arg)
   call get_command_argument(2, arg)
   probe_json = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0) then
      print '(a)', 'check_cli_hessian: usage: check_cli_hessian <terco CLI binary> <scratch .json path>'
      stop 1
   end if
   probe_out = probe_json//'.out'

   cmd = trim(cli_exe)//' vib_thermo_water.xyz -b basis_sets/sto-3g.json --hessian '// &
         '--symmetry-number 2 --silent --json '//probe_json//' > '//probe_out//' 2>&1'
   call execute_command_line(cmd, exitstat=rc, wait=.true.)
   if (rc /= 0) then
      print '(a,i0)', 'check_cli_hessian: the CLI run itself failed, exit code ', rc; stop 1
   end if

   call json%initialize()
   call json%load(filename=probe_json)
   if (json%failed()) then
      print '(a)', 'check_cli_hessian: cannot read '//probe_json; stop 1
   end if
   call json%get('frequencies.n_modes', n_modes_cli, found)
   call expect(found, 'check_cli_hessian: no frequencies.n_modes in the CLI''s JSON output')
   call json%get('frequencies.frequencies_cm1', freq_cli, found)
   call expect(found, 'check_cli_hessian: no frequencies.frequencies_cm1 in the CLI''s JSON output')
   call json%get('thermochemistry.zero_point_energy_hartree', zpe_cli, found)
   call expect(found, 'check_cli_hessian: no thermochemistry.zero_point_energy_hartree')
   call json%get('thermochemistry.thermal_corrections_hartree.to_energy', de_cli, found)
   call json%get('thermochemistry.thermal_corrections_hartree.to_enthalpy', dh_cli, found)
   call json%get('thermochemistry.thermal_corrections_hartree.to_gibbs', dg_cli, found)
   call json%get('thermochemistry.total_energies_hartree.electronic_plus_thermal_G', gtot_cli, found)
   call json%destroy()

   ! pyscf_ref_io is a strictly sequential reader (one tag after another,
   ! that module's own header): walk every tag ahead of the ones this
   ! check needs too, in the order generate.py wrote them -- the same
   ! full list check_vib_thermo_vs_pyscf.F90 reads.
   call read_water_ref()

   call expect(n_modes_cli == n_modes_ref, 'check_cli_hessian: n_modes matches the .ref')
   call expect(size(freq_cli) == n_modes_ref, 'check_cli_hessian: frequencies_cm1 has n_modes entries')
   if (size(freq_cli) == n_modes_ref) then
      do k = 1, n_modes_ref
         call ref_report(nfail, abs(freq_cli(k) - freq_ref(k)) < FREQ_TOL, &
                         'check_cli_hessian: frequency '//trim(to_string(k))//' vs pyscf, cm^-1', &
                         abs(freq_cli(k) - freq_ref(k)))
      end do
   end if
   call ref_report(nfail, abs(zpe_cli - zpe_ref) < ENERGY_TOL, 'check_cli_hessian: ZPE vs pyscf, Eh', &
                   abs(zpe_cli - zpe_ref))
   call ref_report(nfail, abs(de_cli - de_ref) < ENERGY_TOL, 'check_cli_hessian: delta E vs pyscf, Eh', &
                   abs(de_cli - de_ref))
   call ref_report(nfail, abs(dh_cli - dh_ref) < ENERGY_TOL, 'check_cli_hessian: delta H vs pyscf, Eh', &
                   abs(dh_cli - dh_ref))
   call ref_report(nfail, abs(dg_cli - dg_ref) < ENERGY_TOL, 'check_cli_hessian: delta G vs pyscf, Eh', &
                   abs(dg_cli - dg_ref))
   call ref_report(nfail, abs(gtot_cli - gtot_ref) < ENERGY_TOL, 'check_cli_hessian: G_tot vs pyscf, Eh', &
                   abs(gtot_cli - gtot_ref))

   if (nfail > 0) then
      print '(a,i0,a)', 'check_cli_hessian: ', nfail, ' check(s) FAILED'; stop 1
   end if
   print '(a)', 'check_cli_hessian: RESULT: PASS'

contains

   subroutine expect(cond, what)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         print '(a)', what; stop 1
      end if
   end subroutine expect

   !
   ! vib_thermo_water.ref's own tags, in generate.py's write order, up to
   ! and including the ones this check compares -- see
   ! check_vib_thermo_vs_pyscf.F90 for the full tag list this file is a
   ! subset of.
   !
   subroutine read_water_ref()
      integer :: uu, nat2, sigma, mult
      real(dp) :: ztmp(3), coordstmp(3, 3), masstmp(3)
      real(dp) :: temperature_k, pressure_pa, e_elec, dummy
      real(dp) :: hesstmp(9, 9)
      real(dp), allocatable :: redmasstmp(:), fctmp(:)
      real(dp) :: normmodetmp(3, 3, 3)
      integer :: is_linear_ref, n_imag_ref
      real(dp) :: momentstmp(3), rotconsttmp(3)
      call ref_require('vib_thermo_water.ref', uu)
      call ref_read_int(uu, 'NAT', nat2)
      call ref_read_vector(uu, 'Z', nat2, ztmp)
      call ref_read_matrix(uu, 'COORDS', 3, nat2, coordstmp)
      call ref_read_vector(uu, 'MASS', nat2, masstmp)
      call ref_read_int(uu, 'SIGMA', sigma)
      call ref_read_int(uu, 'MULT', mult)
      call ref_read_real(uu, 'TEMPERATURE_K', temperature_k)
      call ref_read_real(uu, 'PRESSURE_PA', pressure_pa)
      call ref_read_real(uu, 'E_ELEC', e_elec)
      call ref_read_matrix(uu, 'HESS', 9, 9, hesstmp)
      call ref_read_int(uu, 'N_MODES', n_modes_ref)
      allocate (freq_ref(n_modes_ref), redmasstmp(n_modes_ref), fctmp(n_modes_ref))
      call ref_read_vector(uu, 'FREQ_CM1', n_modes_ref, freq_ref)
      call ref_read_vector(uu, 'REDUCED_MASS_AMU', n_modes_ref, redmasstmp)
      call ref_read_vector(uu, 'FORCE_CONST_AU', n_modes_ref, fctmp)
      call ref_read_matrix3(uu, 'NORM_MODE', 3, nat2, n_modes_ref, normmodetmp)
      call ref_read_int(uu, 'IS_LINEAR', is_linear_ref)
      call ref_read_int(uu, 'N_IMAGINARY', n_imag_ref)
      call ref_read_vector(uu, 'MOMENTS_AMUA2', 3, momentstmp)
      call ref_read_vector(uu, 'ROTCONST_GHZ', 3, rotconsttmp)
      call ref_read_real(uu, 'ZPE', zpe_ref)
      call ref_read_real(uu, 'DELTA_E', de_ref)
      call ref_read_real(uu, 'DELTA_H', dh_ref)
      call ref_read_real(uu, 'DELTA_G', dg_ref)
      call ref_read_real(uu, 'E_TOT', dummy)
      call ref_read_real(uu, 'H_TOT', dummy)
      call ref_read_real(uu, 'G_TOT', gtot_ref)
      call ref_close(uu)
   end subroutine read_water_ref

end program check_cli_hessian
