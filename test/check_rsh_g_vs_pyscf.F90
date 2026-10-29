module check_rsh_g_suite
   !! Range-separated hybrids on a basis with g functions, against pyscf.
   !!
   !! cc-pVQZ puts a g shell on every heavy atom, so the long-range exchange
   !! pass (erf(omega r12)/r12, j_scale = 0) reaches every exchange engine
   !! the Fock build has: HGP/Rys for the classes up to f, Rys for the g
   !! classes up to TERCO_RYS_G_LTOT and McMurchie-Davidson exchange
   !! (trc_mdg) above it. A full SCF energy is the end-to-end statement;
   !! check_erf_fock_vs_pyscf's water_ccpvqz case is the matrix-level one.
   !!
   !! The basis is spherical, as pyscf runs it (cart=False), so the MD pass
   !! takes its spherical path. Each case dumps terco's own DFT grid (and
   !! VV10 grid) to rsh_g_probe_<case>.bin; test/pyscf_refs/generate.py
   !! converges pyscf on exactly those points and commits the energy. The
   !! grid is the one input the two codes would otherwise build differently.
   !!
   !! Needs TERCO_LMAX >= 4: below it every case skips.
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use trc_api, only: trc_basis_t
   use trc_tables, only: LMAX
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use pyscf_ref_io, only: ref_open, ref_close, ref_read_real
   implicit none
   private
   public :: collect_rsh_g

   real(dp), parameter :: TOL = 1.0e-9_dp
   !! Eh, terco against pyscf on the same grid
   character(len=*), parameter :: BASIS = "basis_sets/cc-pvqz-pyscf17.json"

contains

   subroutine collect_rsh_g(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("water_wb97x", test_water_wb97x), &
                  new_unittest("nh2_wb97x", test_nh2_wb97x), &
                  new_unittest("water_wb97xv", test_water_wb97xv) &
                  ]
   end subroutine collect_rsh_g

   subroutine test_water_wb97x(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "water_wb97x", "water.xyz", "wb97x", 0)
   end subroutine test_water_wb97x

   subroutine test_nh2_wb97x(error)
      !! NH2 (2B1): open shell with a non-degenerate ground state, so both
      !! codes land in the same UKS solution without a shared guess.
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "nh2_wb97x", "nh2.xyz", "wb97x", 1)
   end subroutine test_nh2_wb97x

   subroutine test_water_wb97xv(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "water_wb97xv", "water.xyz", "wb97x-v", 0)
   end subroutine test_water_wb97xv

   subroutine run_case(error, name, xyz, functional, nunpaired)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: name, xyz, functional
      integer, intent(in) :: nunpaired
      type(trc_basis_t) :: bas
      type(error_t) :: err
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: natm, nelec, na, nb, unit, refu
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp) :: e_ref, de
      logical :: ok
      character(len=160) :: msg
      character(len=8) :: fn8

      if (LMAX < 4) then
         call skip_test(error, "SKIP: needs TERCO_LMAX >= 4 (cc-pVQZ has g)")
         return
      end if
      call read_xyz(xyz, natm, zint, at_r)
      nelec = sum(zint)
      nb = (nelec - nunpaired)/2
      na = nelec - nb
      call trc_basis_from_json(BASIS, natm, zint, at_r, bas, err)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      call check(error, bas%spherical, "cc-pVQZ did not come back spherical")
      if (allocated(error)) return
      call bas%to_device()

      opts = trc_scf_options_t()
      opts%functional = functional
      opts%grid_level = 3
      opts%conv_energy = 1.0e-11_dp
      opts%conv_diis = 1.0e-8_dp
      call trc_scf_run(bas, na, nb, opts, res)
      call bas%release()
      call check(error, res%converged, name//": terco SCF did not converge: "//trim(res%message))
      if (allocated(error)) return

      ! The probe generate.py reads: molecule, spin, functional, both grids
      ! and terco's energy (for its own printout only).
      fn8 = functional
      open (newunit=unit, file="rsh_g_probe_"//name//".bin", access="stream", form="unformatted", &
            status="replace")
      write (unit) int(natm, int64), int(res%nspin, int64), int(na, int64), int(nb, int64)
      write (unit) int(zint, int64), at_r
      write (unit) fn8
      write (unit) int(size(res%grid_w), int64)
      write (unit) res%grid_r, res%grid_w
      if (allocated(res%grid_w_nlc)) then
         write (unit) int(size(res%grid_w_nlc), int64)
         write (unit) res%grid_r_nlc, res%grid_w_nlc
      else
         write (unit) 0_int64
      end if
      write (unit) res%energy
      close (unit)

      call ref_open("rsh_g_"//name//".ref", refu, ok)
      call check(error, ok, "missing pyscf_refs/rsh_g_"//name//".ref -- run test/pyscf_refs/generate.py rsh_g_"//name)
      if (allocated(error)) return
      call ref_read_real(refu, "E", e_ref)
      call ref_close(refu)
      de = abs(res%energy - e_ref)
      write (msg, "(a,a,f20.12,a,f20.12,a,es9.2)") name, ": E(terco) ", res%energy, "  E(pyscf) ", e_ref, &
         "  |dE| ", de
      print "(4x,a)", trim(msg)
      call check(error, de <= TOL, trim(msg))
   end subroutine run_case

end module check_rsh_g_suite

program check_rsh_g_vs_pyscf
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, run_selected
   use check_rsh_g_suite, only: collect_rsh_g
   implicit none
   integer :: stat
   character(len=64) :: name
   stat = 0
   ! One case per ctest entry (and per generate.py case): the name is the
   ! first argument; none runs them all.
   if (command_argument_count() >= 1) then
      call get_command_argument(1, name)
      call run_selected(collect_rsh_g, trim(name), error_unit, stat)
   else
      call run_testsuite(collect_rsh_g, error_unit, stat, parallel=.false.)
   end if
   if (stat > 0) then
      write (error_unit, '(i0,a)') stat, " test(s) failed"
      error stop 1
   end if
end program check_rsh_g_vs_pyscf
