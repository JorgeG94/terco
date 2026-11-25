module check_ccsd_grad_ref_suite
   !!
   !! The RI-CCSD and RI-CCSD(T) gradient references (test/pyscf_refs/ccsd_grad_water_*.ref, generate.py's
   !! gen_ccsd_grad) against terco's own RI-CCSD(T) energy, before any terco gradient exists: water distorted to C1
   !! (test/ccgrad_water.xyz), STO-3G and cc-pVDZ, frozen core 0 and 1, cc-pVDZ-RIFIT in the orbital basis's
   !! convention (Cartesian with STO-3G, spherical with cc-pVDZ).
   !!
   !! Per case:
   !!   - terco's SCF, RI-CCSD correlation energy, (T) and the CCSD(T) total against pyscf's DF-RCCSD(T) of the
   !!     reference, at the tolerances of check_capi_rccsd (same method, same molecule class, SCF to 1e-13 / 1e-11):
   !!     SCF 1e-10, E_CCSD 1e-10, E(T) 3e-12, total 3e-11;
   !!   - the reference is self-consistent: the energy rebuilt from its own densities, sum gamma f + 1/2 sum g3 b,
   !!     is the CCSD (and then the CCSD(T)) correlation energy to 1e-9, the one-particle density is
   !!     traceless over the occupied and virtual spaces together, the translational sum of the gradient vanishes
   !!     (1e-9), the Lambda overlaps and norms are positive and the gradient is nonzero in every component
   !!     (the molecule has no symmetry).
   !! The STO-3G cases are Cartesian throughout: terco's auxiliary basis follows the orbital basis's convention (sto-3g.json
   !! is untagged), so cc-pVDZ-RIFIT is Cartesian there too, and the pyscf side is built with cart=True to match.
   !! The achieved differences are printed. The gradient itself is checked against terco's energy in check_ccsd_grad_fd.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use ccsd_grad_support, only: NCASE, CASE_REF, CASE_BASIS, CASE_NFROZEN, ccsd_grad_ref_t, read_ccsd_grad_ref, &
                                terco_energies, geometry
   implicit none
   private
   public :: collect_ccsd_grad_ref

   !> About ten times the worst difference measured over the four cases on a host build: SCF 2.3e-13,
   !> E_CCSD 6.8e-13, E(T) 6.4e-14, total 9.8e-13.
   real(dp), parameter :: TOL_SCF = 3.0e-12_dp
   real(dp), parameter :: TOL_CORR = 7.0e-12_dp
   real(dp), parameter :: TOL_T = 7.0e-13_dp
   real(dp), parameter :: TOL_TOTAL = 1.0e-11_dp
   real(dp), parameter :: TOL_DENSITY = 1.0e-9_dp

contains

   subroutine collect_ccsd_grad_ref(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("the references reproduce their own energy from the densities", test_self_consistent), &
                  new_unittest("terco RI-CCSD(T) matches the references: STO-3G and cc-pVDZ, frozen core 0 and 1", &
                               test_terco_energies) &
                  ]
   end subroutine collect_ccsd_grad_ref

   subroutine test_self_consistent(error)
      type(error_type), allocatable, intent(out) :: error
      type(ccsd_grad_ref_t) :: ref
      integer :: ic, v, k, natm
      real(dp) :: e_corr(2), trace
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)

      call geometry(natm, zint, at_r)
      do ic = 1, NCASE
         call read_ccsd_grad_ref(trim(CASE_REF(ic)), ref)
         call check(error, ref%natm == natm .and. ref%nfrozen == CASE_NFROZEN(ic), trim(CASE_REF(ic))//": header")
         if (allocated(error)) return
         e_corr = [ref%e_ccsd_corr, ref%e_ccsd_corr + ref%e_t]
         do v = 1, 2
            call check(error, abs(ref%probes(10, v) - e_corr(v)) < TOL_DENSITY, trim(CASE_REF(ic))// &
                       ": the densities give "//to_string(ref%probes(10, v), "f20.12")//" for "// &
                       to_string(e_corr(v), "f20.12"))
            if (allocated(error)) return
            trace = ref%probes(6, v) + ref%probes(7, v)
            call check(error, abs(trace) < TOL_DENSITY, trim(CASE_REF(ic))//": gamma trace "//to_string(trace, "es10.2"))
            if (allocated(error)) return
            do k = 1, 3
               call check(error, abs(sum(ref%grad(:, k, v))) < TOL_DENSITY, trim(CASE_REF(ic))// &
                          ": the gradient has a net force along axis "//to_string(k))
               if (allocated(error)) return
            end do
            call check(error, all(abs(ref%grad(:, :, v)) > 1.0e-6_dp), trim(CASE_REF(ic))// &
                       ": a gradient component vanishes: the molecule is not asymmetric enough")
            if (allocated(error)) return
            call check(error, all(ref%probes(1:2, v) > 0.0_dp) .and. ref%probes(5, v) > 0.0_dp .and. &
                       ref%probes(8, v) > 0.0_dp, trim(CASE_REF(ic))//": a norm is not positive")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_self_consistent

   subroutine test_terco_energies(error)
      type(error_type), allocatable, intent(out) :: error
      type(ccsd_grad_ref_t) :: ref
      integer :: ic, natm, rc, iterations
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :)
      real(dp) :: e_scf, e_corr, e_t
      character(len=:), allocatable :: line

      call geometry(natm, zint, at_r)
      do ic = 1, NCASE
         call read_ccsd_grad_ref(trim(CASE_REF(ic)), ref)
         call terco_energies(trim(CASE_BASIS(ic)), CASE_NFROZEN(ic), zint, at_r, e_scf, e_corr, e_t, rc, iterations)
         call check(error, rc == 0, trim(CASE_REF(ic))//": terco's RI-CCSD(T) did not converge, status "//to_string(rc))
         if (allocated(error)) return
         line = "  "//trim(CASE_REF(ic))//": |SCF| "//to_string(abs(e_scf - ref%e_scf), "es9.2")//", |E_CCSD| "// &
                to_string(abs(e_corr - ref%e_ccsd_corr), "es9.2")//", |E(T)| "//to_string(abs(e_t - ref%e_t), "es9.2")// &
                ", |total| "//to_string(abs(e_scf + e_corr + e_t - ref%e_ccsd_t), "es9.2")
         print "(a)", line
         call check(error, abs(e_scf - ref%e_scf) < TOL_SCF, trim(CASE_REF(ic))//": E_SCF "// &
                    to_string(e_scf, "f20.12")//" vs pyscf "//to_string(ref%e_scf, "f20.12"))
         if (allocated(error)) return
         call check(error, abs(e_corr - ref%e_ccsd_corr) < TOL_CORR, trim(CASE_REF(ic))//": E_CCSD "// &
                    to_string(e_corr, "f20.12")//" vs pyscf "//to_string(ref%e_ccsd_corr, "f20.12"))
         if (allocated(error)) return
         call check(error, abs(e_t - ref%e_t) < TOL_T, trim(CASE_REF(ic))//": E(T) "//to_string(e_t, "f20.14")// &
                    " vs pyscf "//to_string(ref%e_t, "f20.14"))
         if (allocated(error)) return
         call check(error, abs(e_scf + e_corr + e_t - ref%e_ccsd_t) < TOL_TOTAL, trim(CASE_REF(ic))// &
                    ": the CCSD(T) total energy")
         if (allocated(error)) return
      end do
   end subroutine test_terco_energies

end module check_ccsd_grad_ref_suite

program check_ccsd_grad_ref
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_ccsd_grad_ref_suite, only: collect_ccsd_grad_ref
   use ccsd_grad_support, only: supported
   implicit none
   integer :: stat
   stat = 0
   if (.not. supported()) then
      print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT, no case can run"
      stop 0
   end if
   call run_testsuite(collect_ccsd_grad_ref, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ccsd_grad_ref
