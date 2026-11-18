module check_cli_tddft_opt_suite
   !!
   !! `terco formaldehyde.xyz --tddft --opt`: the geometry optimisation of the first singlet excited state (CIS /
   !! STO-3G, the n -> pi* state that lengthens C=O), end to end through DL-FIND (a -DTERCO_ENABLE_DLFIND=ON build):
   !!  1. it converges, and the energy the optimiser reports is E_SCF + omega of the followed root: the "tddft" block of
   !!     the JSON, at the final geometry, has the same target-root total energy as `optimization.final_energy`;
   !!  1b. that energy is pyscf's TDA E(S1) at the converged geometry, where pyscf's own excited-state gradient vanishes;
   !!  2. the excited-state gradient the run reports at the end is the optimiser's: maximum component below the
   !!     convergence threshold (so the optimisation reached a stationary point of the S1 surface, not of the ground
   !!     state: the ground-state gradient of that geometry is far from zero);
   !!  3. the geometry relaxed on the excited surface: the energy fell by more than 10 mEh from the vertical one, and
   !!     C=O lengthened by more than 0.05 Angstrom;
   !!  4. the state is watched: a gap window wider than the two roots' spacing makes the run warn about the other root;
   !!  5. terco's own finite-difference Hessian (--opt-init-hessian terco with a Newton step) has no excited-state
   !!     gradient to difference and is refused before the SCF.
   !! Arguments: the terco executable, a scratch .json path and a scratch directory. A test-drive suite.
   !!
   use pic_types, only: dp
   use pic_strings, only: to_string
   use json_module, only: json_file
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_cli_tddft_opt, cli_exe, probe_json, work_dir

   character(len=:), allocatable :: cli_exe, probe_json, work_dir
   real(dp), parameter :: GTOL = 1.0e-5_dp
   real(dp), parameter :: E_S1_PYSCF = -112.214126569198_dp
   !! pyscf 2.14 TDA (CIS) RHF/STO-3G, E_SCF + omega_1, at the geometry this optimisation converges to; pyscf's own
   !! analytic S1 gradient (grad.tdrhf) there is 1.1e-8 Eh/Bohr in every component, so it is the S1 minimum. terco's
   !! STO-3G energies sit about 1.2e-8 Eh above pyscf's (the basis's own rounding).
   real(dp), parameter :: E_TOL = 5.0e-8_dp

contains

   subroutine collect_cli_tddft_opt(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("S1 optimisation of formaldehyde converges to a stationary point of the S1 surface", &
                               test_opt), &
                  new_unittest("a gap window wider than the root spacing warns", test_gap_warning), &
                  new_unittest("--opt --tddft with terco's own Hessian is refused", test_hessian_refused) &
                  ]
   end subroutine collect_cli_tddft_opt

   subroutine stage(rc)
      !! The inputs, copied to the scratch directory the CLI runs in (it writes its xyz files to its working directory).
      integer, intent(out) :: rc
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f formaldehyde.xyz '//work_dir//'/formaldehyde.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', exitstat=rc, wait=.true.)
   end subroutine stage

   function command(extra, json_path, out_path) result(cmd)
      character(len=*), intent(in) :: extra, json_path, out_path
      character(len=:), allocatable :: cmd
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' formaldehyde.xyz -b basis_sets/sto-3g.json --tddft --tddft-roots 2 '// &
            '--opt --opt-gtol 1e-5 --econv 1e-12 --gconv 1e-9 '//extra//' --json '//json_path//' > '//out_path//' 2>&1'
   end function command

   subroutine test_opt(error)
      type(error_type), allocatable, intent(out) :: error
      type(json_file) :: json
      real(dp) :: e_final, e_target, gmax, gtot
      real(dp), allocatable :: traj_e(:), first(:), last(:)
      character(len=:), allocatable :: of
      logical :: converged, found
      integer :: rc, steps
      real(dp) :: co_start, co_end

      call stage(rc)
      call execute_command_line(command('', probe_json, probe_json//'.out'), exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run failed, exit code "//to_string(rc)//" (see "//probe_json//".out)")
      if (allocated(error)) return
      call json%initialize()
      call json%load(filename=probe_json)
      call check(error, .not. json%failed(), "cannot read "//probe_json)
      if (allocated(error)) return
      call json%get('optimization.converged', converged, found)
      call check(error, found .and. converged, "the optimisation did not converge")
      if (.not. allocated(error)) then
         call json%get('optimization.steps', steps, found)
         call json%get('optimization.final_energy', e_final, found)
         call json%get('tddft.target root total energy', e_target, found)
         call check(error, found, "no tddft block in the JSON")
      end if
      if (.not. allocated(error)) then
         call check(error, abs(e_final - E_S1_PYSCF) < E_TOL, "the S1 energy of the optimised geometry against pyscf's, "// &
                    to_string(abs(e_final - E_S1_PYSCF), "es10.2"))
      end if
      if (.not. allocated(error)) then
         call check(error, abs(e_final - e_target) < 1.0e-8_dp, "the optimiser's energy is the target root's "// &
                    "E_SCF + omega: "//to_string(e_final, "f16.10")//" against "//to_string(e_target, "f16.10"))
      end if
      if (.not. allocated(error)) then
         call json%get('gradient of', of, found)
         call check(error, found .and. index(of, "tddft root 1") > 0, "the final gradient is of the excited root")
      end if
      if (.not. allocated(error)) then
         call json%get('optimization.final_gradient_max', gmax, found)
         call json%get('gradient max', gtot, found)
         call check(error, gmax < GTOL, "the optimiser's final gradient is above its threshold, "//to_string(gmax, "es10.2"))
      end if
      if (.not. allocated(error)) then
         ! The gradient the run reports afterwards is the S1 one, re-evaluated at the final geometry: also stationary.
         call check(error, gtot < 3.0_dp*GTOL, "the excited-state gradient at the final geometry, per-atom norm "// &
                    to_string(gtot, "es10.2"))
      end if
      if (.not. allocated(error)) then
         call json%get('optimization.trajectory_energies', traj_e, found)
         call check(error, found, "no trajectory")
         if (.not. allocated(error)) then
            call check(error, traj_e(1) - e_final > 1.0e-2_dp, "the energy should fall by more than 10 mEh from the "// &
                       "vertical state, fell "//to_string(traj_e(1) - e_final, "es10.2"))
         end if
      end if
      if (.not. allocated(error)) then
         call json%get('optimization.final_geometry_angstrom(1)', first, found)
         call json%get('optimization.final_geometry_angstrom(2)', last, found)
         co_end = sqrt(sum((first - last)**2))
         co_start = 1.203_dp
         call check(error, co_end - co_start > 0.05_dp, "C=O should lengthen by more than 0.05 Angstrom in S1, now "// &
                    to_string(co_end, "f8.4"))
      end if
      call json%destroy()
   end subroutine test_opt

   subroutine test_gap_warning(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=512) :: line
      logical :: seen
      integer :: rc, u, ios

      call stage(rc)
      call execute_command_line(command('--tddft-gap-warn 1.0', probe_json//'.gap.json', probe_json//'.gap.out'), &
                                exitstat=rc, wait=.true.)
      call check(error, rc == 0, "the CLI run failed, exit code "//to_string(rc))
      if (allocated(error)) return
      seen = .false.
      open (newunit=u, file=probe_json//'.gap.out', status="old", action="read", iostat=ios)
      do while (ios == 0)
         read (u, "(a)", iostat=ios) line
         if (ios == 0 .and. index(line, "lies only") > 0) seen = .true.
      end do
      close (u)
      call check(error, seen, "no warning about the other root (see "//probe_json//".gap.out)")
   end subroutine test_gap_warning

   subroutine test_hessian_refused(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=512) :: line
      logical :: seen
      integer :: rc, u, ios

      call stage(rc)
      call execute_command_line(command('--opt-algorithm nr --opt-init-hessian terco', probe_json//'.nr.json', &
                                        probe_json//'.nr.out'), exitstat=rc, wait=.true.)
      call check(error, rc /= 0, "not refused (exit code 0)")
      if (allocated(error)) return
      seen = .false.
      open (newunit=u, file=probe_json//'.nr.out', status="old", action="read", iostat=ios)
      do while (ios == 0)
         read (u, "(a)", iostat=ios) line
         if (ios == 0 .and. index(line, "no excited-state gradient to difference") > 0) seen = .true.
      end do
      close (u)
      call check(error, seen, "the refusal does not say why (see "//probe_json//".nr.out)")
   end subroutine test_hessian_refused

end module check_cli_tddft_opt_suite

program check_cli_tddft_opt
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_cli_tddft_opt_suite, only: collect_cli_tddft_opt, cli_exe, probe_json, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat

   call get_command_argument(1, arg)
   cli_exe = trim(arg)
   call get_command_argument(2, arg)
   probe_json = trim(arg)
   call get_command_argument(3, arg)
   work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(probe_json) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_cli_tddft_opt: usage: check_cli_tddft_opt <terco CLI> <scratch .json> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_cli_tddft_opt, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_cli_tddft_opt
