module check_ccsd_grad_fd_suite
   !!
   !! The RI-CCSD and RI-CCSD(T) gradient references against central differences of terco's OWN energy.
   !! check_ccsd_grad_ref anchors the references' energies to terco; this closes the loop: the gradient in
   !! test/pyscf_refs/ccsd_grad_water_*.ref (generate.py's gen_ccsd_grad) is the derivative of the energy terco
   !! computes, not only of pyscf's.
   !!
   !! The energies are terco's trc_run_rccsd totals (exact-ERI RHF to 1e-13 / 1e-11, amplitudes to 1e-12 / 1e-10,
   !! cc-pVDZ-RIFIT), CCSD and CCSD(T) from the same run. The derivative is the five-point stencil
   !! (f(-2h) - 8 f(-h) + 8 f(h) - f(2h)) / 12h with h = 0.01 Bohr: truncation error h^4 f^(5) / 30 and energy
   !! noise (1e-11) / h. Every coordinate of the four cases (STO-3G and cc-pVDZ, frozen core 0 and 1): 36 coordinates,
   !! four energies each, each an SCF and a RI-CCSD(T) of a three-atom molecule, about a minute on the host, so that run
   !! is LONG. With the argument `short` only STO-3G with one frozen core orbital and three coordinates run, which every
   !! pull request does. The largest difference per case is printed.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use ccsd_grad_support, only: ccsd_grad_ref_t, read_ccsd_grad_ref, terco_energies, geometry, CASE_REF, CASE_BASIS, &
                                CASE_NFROZEN
   implicit none
   private
   public :: collect_ccsd_grad_fd, collect_ccsd_grad_fd_short

   real(dp), parameter :: STEP = 1.0e-2_dp
   real(dp), parameter :: TOL_FD = 1.0e-8_dp
   !! measured worst difference 3.2e-9 (O z), the h^4 truncation of the stencil; the energy noise (1e-11)/h is 1e-9

contains

   subroutine collect_ccsd_grad_fd(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("STO-3G frozen core 0: every coordinate", test_sto3g_fc0), &
                  new_unittest("STO-3G frozen core 1: every coordinate", test_sto3g_fc1), &
                  new_unittest("cc-pVDZ frozen core 0: every coordinate", test_ccpvdz_fc0), &
                  new_unittest("cc-pVDZ frozen core 1: every coordinate", test_ccpvdz_fc1) &
                  ]
   end subroutine collect_ccsd_grad_fd

   subroutine collect_ccsd_grad_fd_short(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("STO-3G frozen core 1: three coordinates", test_short)]
   end subroutine collect_ccsd_grad_fd_short

   subroutine test_short(error)
      !! O z, H1 y and H2 x: one coordinate per atom and one per axis
      type(error_type), allocatable, intent(out) :: error
      call fd_case(error, 2, [3, 5, 7])
   end subroutine test_short

   subroutine test_sto3g_fc0(error)
      type(error_type), allocatable, intent(out) :: error
      call fd_case(error, 1, [1, 2, 3, 4, 5, 6, 7, 8, 9])
   end subroutine test_sto3g_fc0

   subroutine test_sto3g_fc1(error)
      type(error_type), allocatable, intent(out) :: error
      call fd_case(error, 2, [1, 2, 3, 4, 5, 6, 7, 8, 9])
   end subroutine test_sto3g_fc1

   subroutine test_ccpvdz_fc0(error)
      type(error_type), allocatable, intent(out) :: error
      call fd_case(error, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9])
   end subroutine test_ccpvdz_fc0

   subroutine test_ccpvdz_fc1(error)
      type(error_type), allocatable, intent(out) :: error
      call fd_case(error, 4, [1, 2, 3, 4, 5, 6, 7, 8, 9])
   end subroutine test_ccpvdz_fc1

   subroutine fd_case(error, ic, coords)
      !! case `ic` of ccsd_grad_support, the flattened coordinates (atom - 1)*3 + axis in `coords`
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ic, coords(:)

      type(ccsd_grad_ref_t) :: ref
      integer :: natm, rc, k, atom, axis, j
      integer, allocatable :: zint(:)
      real(dp), allocatable :: at_r(:, :), moved(:, :)
      real(dp) :: e_scf, e_corr, e_t, f_ccsd(-2:2), f_t(-2:2), fd_ccsd, fd_t, worst(2), diff(2)

      call geometry(natm, zint, at_r)
      call read_ccsd_grad_ref(trim(CASE_REF(ic)), ref)
      worst = 0.0_dp
      do k = 1, size(coords)
         atom = (coords(k) - 1)/3 + 1
         axis = coords(k) - 3*(atom - 1)
         do j = -2, 2
            if (j == 0) cycle
            moved = at_r
            moved(axis, atom) = moved(axis, atom) + real(j, dp)*STEP
            call terco_energies(trim(CASE_BASIS(ic)), CASE_NFROZEN(ic), zint, moved, e_scf, e_corr, e_t, rc)
            call check(error, rc == 0, trim(CASE_REF(ic))//": terco's RI-CCSD(T) failed at a displaced geometry")
            if (allocated(error)) return
            f_ccsd(j) = e_scf + e_corr
            f_t(j) = e_scf + e_corr + e_t
         end do
         fd_ccsd = (f_ccsd(-2) - 8.0_dp*f_ccsd(-1) + 8.0_dp*f_ccsd(1) - f_ccsd(2))/(12.0_dp*STEP)
         fd_t = (f_t(-2) - 8.0_dp*f_t(-1) + 8.0_dp*f_t(1) - f_t(2))/(12.0_dp*STEP)
         diff = [abs(fd_ccsd - ref%grad(atom, axis, 1)), abs(fd_t - ref%grad(atom, axis, 2))]
         worst = max(worst, diff)
         print "(a)", "  "//trim(CASE_REF(ic))//" atom "//to_string(atom)//" axis "//to_string(axis)//": CCSD "// &
            to_string(fd_ccsd, "f14.9")//" vs "//to_string(ref%grad(atom, axis, 1), "f14.9")//", CCSD(T) "// &
            to_string(fd_t, "f14.9")//" vs "//to_string(ref%grad(atom, axis, 2), "f14.9")//", |diff| "// &
            to_string(diff(1), "es9.2")//" "//to_string(diff(2), "es9.2")
         call check(error, diff(1) < TOL_FD, trim(CASE_REF(ic))//": CCSD finite difference "//to_string(fd_ccsd, "f14.9")// &
                    " vs reference "//to_string(ref%grad(atom, axis, 1), "f14.9"))
         if (allocated(error)) return
         call check(error, diff(2) < TOL_FD, trim(CASE_REF(ic))//": CCSD(T) finite difference "//to_string(fd_t, "f14.9")// &
                    " vs reference "//to_string(ref%grad(atom, axis, 2), "f14.9"))
         if (allocated(error)) return
      end do
      print "(a)", "  "//trim(CASE_REF(ic))//": worst |finite difference - reference| CCSD "//to_string(worst(1), "es9.2")// &
         ", CCSD(T) "//to_string(worst(2), "es9.2")
   end subroutine fd_case

end module check_ccsd_grad_fd_suite

program check_ccsd_grad_fd
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_ccsd_grad_fd_suite, only: collect_ccsd_grad_fd, collect_ccsd_grad_fd_short
   use ccsd_grad_support, only: supported
   implicit none
   integer :: stat
   character(len=16) :: mode
   stat = 0
   call get_command_argument(1, mode)
   if (.not. supported()) then
      print "(a)", "SKIP: this build's density-fitting ceilings are below cc-pVDZ-RIFIT, no case can run"
      stop 0
   end if
   if (mode == "short") then
      call run_testsuite(collect_ccsd_grad_fd_short, error_unit, stat, parallel=.false.)
   else
      call run_testsuite(collect_ccsd_grad_fd, error_unit, stat, parallel=.false.)
   end if
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_ccsd_grad_fd
