!
! Invalid --opt-* combinations are refused by app/main.F90's own CLI
! validation BEFORE trc_dlfind_optimize/api_dl_find ever runs -- the
! point of validating up front at all (DL-FIND itself answers an invalid
! setting with `error stop`, which is unacceptable to let through; see
! app/trc_dlfind_bridge.F90's own header). Every case here must exit
! non-zero with the SPECIFIC message app/main.F90's own `die()` call
! raises for it, and NEVER reach an "SCF :" line (the marker that
! validation passed and a run started) -- asserting only "non-zero exit"
! would let a missing input file or an unrelated crash pass as "refused"
! just as easily as the intended validation.
!
! A test-drive suite.
!
module check_opt_refusals_suite
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed
   implicit none
   private
   public :: collect_opt_refusals, cli_exe, work_dir

   character(len=:), allocatable :: cli_exe
      !! The terco CLI binary, from the command line.
   character(len=:), allocatable :: work_dir
      !! A scratch directory this check owns, for the --silent output.

contains

   subroutine collect_opt_refusals(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("--opt-saddle-method dimer needs --opt-target saddle", test_dimer_needs_saddle), &
         new_unittest("--opt-connect needs --opt-target saddle", test_connect_needs_saddle), &
         new_unittest("an unknown --opt-algorithm is refused", test_bad_algorithm), &
         new_unittest("an unknown --opt-coords is refused", test_bad_coords), &
         new_unittest("an unknown --opt-target is refused", test_bad_target), &
         new_unittest("an unknown --opt-saddle-method is refused", test_bad_saddle_method), &
         new_unittest("--opt-algorithm together with --opt-target saddle is refused", test_algorithm_with_saddle), &
         new_unittest("--opt-zero-modes above 3*natm is refused", test_zero_modes_too_large), &
         new_unittest("a negative --opt-zero-modes is refused", test_zero_modes_negative), &
         new_unittest("an unknown --opt-init-hessian is refused", test_bad_init_hessian), &
         new_unittest("an unknown --opt-hessian-update is refused", test_bad_hessian_update), &
         new_unittest("an unknown --opt-linesearch is refused", test_bad_linesearch), &
         new_unittest("a --opt-print-level DL-FIND does not accept is refused", test_bad_print_level), &
         new_unittest("--opt-init-hessian with an algorithm that ignores it is refused", &
                      test_init_hessian_wrong_algorithm), &
         new_unittest("--opt-hessian-update with an algorithm that ignores it is refused", &
                      test_hessian_update_wrong_algorithm), &
         new_unittest("--opt-max-updates with an algorithm that ignores it is refused", &
                      test_max_updates_wrong_algorithm), &
         new_unittest("--opt-zero-modes with an algorithm that ignores it is refused", &
                      test_zero_modes_wrong_algorithm), &
         new_unittest("--opt-soft-modes with an algorithm that ignores it is refused", &
                      test_soft_modes_wrong_algorithm), &
         new_unittest("--opt-algorithm without --opt is refused, not silently ignored", &
                      test_opt_option_without_opt), &
         new_unittest("--opt-algorithm prfo on a single atom is refused (no internal DOF)", &
                      test_single_atom_prfo), &
         new_unittest("--opt-max-updates 0 is accepted (DL-FIND's own 'always recompute')", &
                      test_max_updates_zero_accepted), &
         new_unittest("--opt-global-pop-size without --opt-algorithm stochastic/genetic is refused", &
                      test_global_opt_without_algorithm), &
         new_unittest("--opt-seed without --opt-algorithm stochastic/genetic or --opt-saddle-method "// &
                      "dimer is refused", test_seed_without_global_or_dimer), &
         new_unittest("--opt-seed IS accepted for --opt-saddle-method dimer", test_seed_with_dimer_accepted), &
         new_unittest("--opt-global-mutation-rate needs --opt-algorithm genetic, not stochastic", &
                      test_global_mutation_wrong_algorithm), &
         new_unittest("--opt-global-contraction needs --opt-algorithm stochastic, not genetic", &
                      test_global_contraction_wrong_algorithm), &
         new_unittest("an out-of-range --opt-global-distribution is refused", test_global_distribution_range), &
         new_unittest("an out-of-range --opt-global-contraction is refused", test_global_contraction_range), &
         new_unittest("--opt-global-pop-size 0 is refused", test_global_pop_size_zero), &
         new_unittest("--opt-neb together with --opt-target saddle is refused", test_neb_with_saddle), &
         new_unittest("--opt-neb together with an algorithm it cannot use is refused", test_neb_bad_algorithm), &
         new_unittest("a NEB sub-option without --opt-neb itself is refused", test_neb_suboption_orphan), &
         new_unittest("--opt-neb-images below DL-FIND's own minimum is refused", test_neb_images_too_few), &
         new_unittest("an unknown --opt-neb-endpoints is refused", test_neb_bad_endpoints), &
         new_unittest("a missing --opt-neb endpoint file is refused", test_neb_missing_file), &
         new_unittest("an --opt-neb endpoint with mismatched atoms is refused", test_neb_mismatched_atoms), &
         new_unittest("--opt-neb-via without --opt-neb itself is refused", test_neb_via_orphan), &
         new_unittest("a missing --opt-neb-via file is refused", test_neb_via_missing_file), &
         new_unittest("an --opt-neb-via frame with mismatched atoms is refused", test_neb_via_mismatched_atoms), &
         new_unittest("two atoms landing on top of each other on the initial NEB path is refused", &
                      test_neb_initial_path_overlap), &
         new_unittest("--opt-frozen together with --opt-coords dlc is refused", test_frozen_with_dlc), &
         new_unittest("--opt-frozen atom index 0 is refused", test_frozen_index_zero), &
         new_unittest("--opt-frozen atom index out of range is refused", test_frozen_index_too_large), &
         new_unittest("a non-numeric --opt-frozen entry is refused", test_frozen_bogus), &
         new_unittest("--opt-constrain with Cartesian coordinates is refused", test_constrain_cartesian_coords), &
         new_unittest("--opt-constrain 'cartesian' with --opt-coords dlc is refused", &
                      test_constrain_cartesian_type_dlc), &
         new_unittest("an unknown --opt-constrain type is refused", test_constrain_bad_type), &
         new_unittest("an out-of-range --opt-constrain atom index is refused", test_constrain_index_too_large), &
         new_unittest("--opt-constrain atom index 0 is refused", test_constrain_index_zero), &
         new_unittest("--opt-constrain with too few atom indices is refused", test_constrain_too_few_indices), &
         new_unittest("--opt-constrain with a repeated atom in one entry is refused", &
                      test_constrain_repeated_atom), &
         new_unittest("--opt-constrain across an HDLC residue boundary is refused", &
                      test_constrain_crosses_residue), &
         new_unittest("--opt-constrain touching a --opt-frozen atom is refused", test_constrain_touches_frozen), &
         new_unittest("none of the refusal cases above started an SCF or left xyz output in test/ itself", &
                      test_no_stray) &
         ]
   end subroutine collect_opt_refusals

   subroutine test_dimer_needs_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-saddle-method dimer', 'needs --opt-target saddle', error)
   end subroutine test_dimer_needs_saddle

   subroutine test_connect_needs_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-connect', 'needs --opt-target saddle', error)
   end subroutine test_connect_needs_saddle

   subroutine test_bad_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm bogus', '--opt-algorithm must be', error)
   end subroutine test_bad_algorithm

   subroutine test_bad_coords(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-coords bogus', '--opt-coords must be', error)
   end subroutine test_bad_coords

   subroutine test_bad_target(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target bogus', '--opt-target must be', error)
   end subroutine test_bad_target

   subroutine test_bad_saddle_method(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method bogus', &
                  '--opt-saddle-method must be', error)
   end subroutine test_bad_saddle_method

   subroutine test_algorithm_with_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-algorithm lbfgs', &
                  'is not used with --opt-target', error)
   end subroutine test_algorithm_with_saddle

   subroutine test_zero_modes_too_large(error)
      type(error_type), allocatable, intent(out) :: error
      ! nzero out of range: 3*natm = 9 for water.
      call refuse('water.xyz', '--opt --opt-zero-modes 10', 'must be between 0 and 3*natm', error)
   end subroutine test_zero_modes_too_large

   subroutine test_zero_modes_negative(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-zero-modes -1', 'must be between 0 and 3*natm', error)
   end subroutine test_zero_modes_negative

   subroutine test_bad_init_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method prfo --opt-init-hessian bogus', &
                  '--opt-init-hessian must be', error)
   end subroutine test_bad_init_hessian

   subroutine test_bad_hessian_update(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-target saddle --opt-saddle-method prfo --opt-hessian-update bogus', &
                  '--opt-hessian-update must be', error)
   end subroutine test_bad_hessian_update

   subroutine test_bad_linesearch(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-linesearch bogus', '--opt-linesearch must be', error)
   end subroutine test_bad_linesearch

   subroutine test_bad_print_level(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-print-level 3', '--opt-print-level must be', error)
   end subroutine test_bad_print_level

   !
   ! The Hessian-related options are meaningful only for P-RFO/NR (or
   ! --opt-target saddle, which forces one of P-RFO/dimer) -- refused
   ! for the CLI's own default algorithm (lbfgs) rather than silently
   ! ignored (app/main.F90's own opt_hessian_relevant/opt_zero_soft_relevant).
   !
   subroutine test_init_hessian_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-init-hessian terco', '--opt-init-hessian needs --opt-algorithm', error)
   end subroutine test_init_hessian_wrong_algorithm

   subroutine test_hessian_update_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-hessian-update none', '--opt-hessian-update needs --opt-algorithm', &
                  error)
   end subroutine test_hessian_update_wrong_algorithm

   subroutine test_max_updates_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-max-updates 10', '--opt-max-updates needs --opt-algorithm', error)
   end subroutine test_max_updates_wrong_algorithm

   subroutine test_zero_modes_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      ! In range (2 <= 3*natm=9) so the bound check does not fire first:
      ! this must be refused by the algorithm-relevance check instead.
      call refuse('water.xyz', '--opt --opt-zero-modes 2', '--opt-zero-modes needs --opt-algorithm', error)
   end subroutine test_zero_modes_wrong_algorithm

   subroutine test_soft_modes_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-soft-modes 0.01', '--opt-soft-modes needs --opt-algorithm', error)
   end subroutine test_soft_modes_wrong_algorithm

   !
   ! Every --opt-* option is parsed only inside `if (r%opt) then` --
   ! given without --opt itself, it must be refused, not silently
   ! ignored (app/main.F90's own check just before that block).
   !
   subroutine test_opt_option_without_opt(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt-algorithm lbfgs', '--opt-algorithm needs --opt', error)
   end subroutine test_opt_option_without_opt

   !
   ! A single atom has no internal degrees of freedom for P-RFO/NR/the
   ! dimer method's own eigenvector-following to search over.
   !
   subroutine test_single_atom_prfo(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('o_atom.xyz', '--opt --opt-algorithm prfo', 'need more than one atom', error)
   end subroutine test_single_atom_prfo

   !
   ! DL-FIND itself accepts maxupd=0 (dlf_formstep.f90's own
   ! iupd>=maxupd test, iupd starting at 0: it recomputes the Hessian
   ! fresh every step) -- NOT refused, unlike the cases above.
   !
   subroutine test_max_updates_zero_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      call execute_command_line('mkdir -p '//work_dir, exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -rf basis_sets '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      out_path = 'opt_refusal_zero.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-target saddle --opt-saddle-method prfo --opt-max-updates 0 --opt-maxiter 1 --silent > '// &
            out_path//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      ! --opt-maxiter 1 will not converge (exit non-zero), but NOT with
      ! the "--opt-max-updates must be at least 0" validation message --
      ! this is a normal non-convergence exit, never a CLI refusal.
      if (scan_for(work_dir//'/'//out_path, '--opt-max-updates must be')) then
         call test_failed(error, "--opt-max-updates 0 was refused at the CLI, but DL-FIND itself accepts it")
      end if
   end subroutine test_max_updates_zero_accepted

   !
   ! --opt-global-*/--opt-seed without --opt-algorithm stochastic/genetic
   ! (and --opt-seed without --opt-saddle-method dimer either), an
   ! --opt-global-* parameter the OTHER search does not read, and an
   ! out-of-range one (parse_opt_global_options, app/main.F90).
   !
   subroutine test_global_opt_without_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-global-pop-size 10', 'needs --opt-algorithm stochastic or genetic', &
                  error)
   end subroutine test_global_opt_without_algorithm

   subroutine test_seed_without_global_or_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-seed 1', &
                  'needs --opt-algorithm stochastic/genetic or --opt-saddle-method dimer', error)
   end subroutine test_seed_without_global_or_dimer

   !
   ! --opt-seed is accepted together with --opt-saddle-method dimer
   ! (trc_dlfind_bridge's own dimer-seed fix, this file's header) -- a
   ! connect run converges well inside --opt-maxiter 150, so reaching an
   ! "SCF :" line here means the --opt-seed/dimer combination itself was
   ! accepted, the opposite assertion from `refuse`.
   !
   subroutine test_seed_with_dimer_accepted(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc
      logical :: saw_refusal
      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f hcn_ts_guess.xyz '//work_dir//'/hcn_ts_guess.xyz', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      out_path = 'opt_refusal_seed_dimer.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' hcn_ts_guess.xyz -b basis_sets/sto-3g.json --opt '// &
            '--opt-target saddle --opt-saddle-method dimer --opt-seed 1 --opt-maxiter 150 --opt-gtol 1e-5 '// &
            '--silent > '//out_path//' 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)
      saw_refusal = scan_for(work_dir//'/'//out_path, 'needs --opt-algorithm')
      if (saw_refusal) then
         call test_failed(error, "--opt-seed with --opt-saddle-method dimer was refused at the CLI")
      end if
   end subroutine test_seed_with_dimer_accepted

   subroutine test_global_mutation_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm stochastic --opt-global-mutation-rate 0.1', &
                  'needs --opt-algorithm genetic', error)
   end subroutine test_global_mutation_wrong_algorithm

   subroutine test_global_contraction_wrong_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm genetic --opt-global-contraction 0.5', &
                  'needs --opt-algorithm stochastic', error)
   end subroutine test_global_contraction_wrong_algorithm

   subroutine test_global_distribution_range(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm stochastic --opt-global-distribution 5', &
                  'must be 1 (box), 2 (sphere) or 3 (Gaussian)', error)
   end subroutine test_global_distribution_range

   subroutine test_global_contraction_range(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm stochastic --opt-global-contraction 1.5', &
                  'must be in (0, 1)', error)
   end subroutine test_global_contraction_range

   subroutine test_global_pop_size_zero(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-algorithm genetic --opt-global-pop-size 0', &
                  '--opt-global-pop-size must be at least 1', error)
   end subroutine test_global_pop_size_zero

   !
   ! --opt-neb: mutually exclusive with a saddle target/algorithm this
   ! bridge does not offer for a chain.
   !
   subroutine test_neb_with_saddle(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb water.xyz --opt-target saddle', 'needs --opt-target minimum', error)
   end subroutine test_neb_with_saddle

   subroutine test_neb_bad_algorithm(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb water.xyz --opt-algorithm prfo', &
                  'needs a Hessian-free, non-global --opt-algorithm', error)
   end subroutine test_neb_bad_algorithm

   subroutine test_neb_suboption_orphan(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb-images 10', 'need --opt-neb', error)
   end subroutine test_neb_suboption_orphan

   subroutine test_neb_images_too_few(error)
      type(error_type), allocatable, intent(out) :: error
      ! Below DL-FIND's own minimum (3 -- dlf_neb.f90: "More than 2
      ! images have to be used in NEB"); water.xyz as its own endpoint is
      ! fine here, this never reaches DL-FIND either way.
      call refuse('water.xyz', '--opt --opt-neb water.xyz --opt-neb-images 2', '--opt-neb-images must be at '// &
                  'least 3', error)
   end subroutine test_neb_images_too_few

   subroutine test_neb_bad_endpoints(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb water.xyz --opt-neb-endpoints bogus', &
                  '--opt-neb-endpoints must be free, perpendicular or frozen', error)
   end subroutine test_neb_bad_endpoints

   subroutine test_neb_missing_file(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb does_not_exist.xyz', 'cannot open geometry file', error)
   end subroutine test_neb_missing_file

   !
   ! hcn.xyz (H,C,N) as the endpoint of a water (O,H,H) start: SAME atom
   ! count (3) but different elements, so check_neb_endpoint's per-atom Z
   ! check fires, not its atom-count check.
   !
   subroutine test_neb_mismatched_atoms(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb hcn.xyz', 'differs between the start geometry', error)
   end subroutine test_neb_mismatched_atoms

   subroutine test_neb_via_orphan(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb-via hcn.xyz', 'need --opt-neb', error)
   end subroutine test_neb_via_orphan

   subroutine test_neb_via_missing_file(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb water.xyz --opt-neb-via does_not_exist.xyz', &
                  "--opt-neb-via: cannot open", error)
   end subroutine test_neb_via_missing_file

   subroutine test_neb_via_mismatched_atoms(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-neb water.xyz --opt-neb-via hcn.xyz', 'differs between the '// &
                  'start geometry', error)
   end subroutine test_neb_via_mismatched_atoms

   !
   ! water.xyz with its two H atoms swapped as the NEB product endpoint:
   ! linear interpolation toward it puts both H atoms on the SAME point
   ! (the midpoint of their swapped positions, water's own mirror plane)
   ! at the image halfway along the path -- a cheap, deterministic way
   ! to make app/trc_dlfind_bridge.F90%check_neb_initial_path's "two
   ! atoms closer than 0.5 Bohr" pre-check actually fire, not just reason
   ! about it. This refusal comes back from trc_dlfind_optimize itself
   ! (opt_res%failed), not a CLI `die`, but still before any SCF runs --
   ! the SAME "no SCF line" assertion `refuse` makes applies.
   !
   subroutine test_neb_initial_path_overlap(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: u, rc
      call execute_command_line('mkdir -p '//work_dir, exitstat=rc, wait=.true.)
      open (newunit=u, file=work_dir//'/water_hswap.xyz', status='replace', action='write')
      write (u, '(a)') '3'
      write (u, '(a)') 'water.xyz with its two H atoms swapped'
      write (u, '(a)') 'O   0.0000000   0.0000000   0.1173000'
      write (u, '(a)') 'H   0.0000000  -0.7572000  -0.4692000'
      write (u, '(a)') 'H   0.0000000   0.7572000  -0.4692000'
      close (u)
      call refuse('water.xyz', '--opt --opt-neb water_hswap.xyz --opt-neb-images 10', 'Bohr apart at image', error)
   end subroutine test_neb_initial_path_overlap

   !
   ! --opt-frozen: refused together with pure DLC/DLC-TC (DL-FIND itself
   ! rejects a frozen atom there -- dlf_coords.f90's own check), and a bad
   ! index.
   !
   subroutine test_frozen_with_dlc(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-frozen 1 --opt-coords dlc', &
                  'refused together with --opt-coords dlc/dlc-tc', error)
   end subroutine test_frozen_with_dlc

   subroutine test_frozen_index_zero(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-frozen 0', 'is not a valid 1-based atom index', error)
   end subroutine test_frozen_index_zero

   subroutine test_frozen_index_too_large(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-frozen 99', 'is not a valid 1-based atom index', error)
   end subroutine test_frozen_index_too_large

   subroutine test_frozen_bogus(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-frozen bogus', 'is not a valid 1-based atom index', error)
   end subroutine test_frozen_bogus

   !
   ! --opt-constrain: refused with Cartesian coordinates outright, and
   ! 'cartesian i' refused with DLC/DLC-TC specifically (dlf_hdlc_
   ! interface.f90's own hdlc%internal=(icoord>=3) check); plus bad
   ! types, bad indices, repeated atoms, residue boundaries and frozen
   ! atoms.
   !
   subroutine test_constrain_cartesian_coords(error)
      type(error_type), allocatable, intent(out) :: error
      ! default --opt-coords is cartesian
      call refuse('water.xyz', '--opt --opt-constrain "angle 2 1 3"', &
                  'needs --opt-coords hdlc, hdlc-tc, dlc or dlc-tc', error)
   end subroutine test_constrain_cartesian_coords

   subroutine test_constrain_cartesian_type_dlc(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-constrain "cartesian 1" --opt-coords dlc', &
                  "needs --opt-coords hdlc or hdlc-tc", error)
   end subroutine test_constrain_cartesian_type_dlc

   subroutine test_constrain_bad_type(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-constrain "bogus 1 2" --opt-coords hdlc', 'unknown type', error)
   end subroutine test_constrain_bad_type

   subroutine test_constrain_index_too_large(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-constrain "bond 1 99" --opt-coords hdlc', &
                  'atom index out of range', error)
   end subroutine test_constrain_index_too_large

   !
   ! A typed 0 used to pass parse_constraints' own range check (it tested
   ! idx(k) /= 0 "is this slot even used", not "is 0 a valid atom"), then
   ! reached DL-FIND's own assign_cons, which reports "crosses residue
   ! boundary" for atom 0 (it belongs to no residue) and `error stop`s.
   !
   subroutine test_constrain_index_zero(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-constrain "bond 0 2" --opt-coords hdlc', &
                  'atom index out of range', error)
   end subroutine test_constrain_index_zero

   subroutine test_constrain_too_few_indices(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-constrain "bond 1" --opt-coords hdlc', &
                  'does not have the atom indices', error)
   end subroutine test_constrain_too_few_indices

   subroutine test_constrain_repeated_atom(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-constrain "angle 1 2 1" --opt-coords hdlc', &
                  'the same atom appears twice', error)
   end subroutine test_constrain_repeated_atom

   !
   ! w3.xyz is three separate, far-apart water molecules -- three HDLC
   ! residues under --opt-coords hdlc (trc_constants%trc_covalent_radius-
   ! based connected components, build_opt_residues). Atom 1 (O of the
   ! first molecule) and atom 4 (H of the second) are in different
   ! residues, which DL-FIND's own assign_cons refuses with "crosses
   ! residue boundary" before this fix caught it here instead.
   !
   subroutine test_constrain_crosses_residue(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('w3.xyz', '--opt --opt-constrain "bond 1 4" --opt-coords hdlc', &
                  'are in different HDLC residues', error)
   end subroutine test_constrain_crosses_residue

   subroutine test_constrain_touches_frozen(error)
      type(error_type), allocatable, intent(out) :: error
      call refuse('water.xyz', '--opt --opt-frozen 1 --opt-constrain "bond 1 2" --opt-coords hdlc', &
                  'is frozen (--opt-frozen)', error)
   end subroutine test_constrain_touches_frozen

   !
   ! None of the refusal cases above ever reaches --opt's own
   ! write-the-xyz-outputs code (every one is refused before DL-FIND
   ! runs at all) -- checked anyway, for the same reason every other opt
   ! check checks it: a future refusal that stops catching an invalid
   ! combination should fail LOUDLY here, not just silently start
   ! optimising and littering test/ again.
   !
   subroutine test_no_stray(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: exists
      inquire (file='water_opt.xyz', exist=exists)
      if (exists) then
         call test_failed(error, "water_opt.xyz was left in test/ itself"); return
      end if
      inquire (file='water_traj.xyz', exist=exists)
      if (exists) call test_failed(error, "water_traj.xyz was left in test/ itself")
   end subroutine test_no_stray

   !
   ! Runs `xyz --opt <flags> --silent` (flags must include --opt unless
   ! the case is testing an --opt-* option given WITHOUT --opt) and
   ! requires: a non-zero exit, no "SCF :" line (validation refused
   ! before any run started), and `expect_substr` present somewhere in
   ! the redirected stdout/stderr -- the specific message, not just "it
   ! failed somehow".
   !
   subroutine refuse(xyz, flags, expect_substr, error)
      character(len=*), intent(in) :: xyz, flags, expect_substr
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: cmd, out_path
      integer :: rc, u, ios
      character(len=1024) :: line
      logical :: saw_scf_line, saw_message

      call execute_command_line('mkdir -p '//work_dir//'/basis_sets', exitstat=rc, wait=.true.)
      call execute_command_line('cp -f '//xyz//' '//work_dir//'/'//xyz, exitstat=rc, wait=.true.)
      call execute_command_line('cp -f basis_sets/sto-3g.json '//work_dir//'/basis_sets/sto-3g.json', &
                                exitstat=rc, wait=.true.)
      out_path = work_dir//'/opt_refusal.out'
      cmd = 'cd '//work_dir//' && '//trim(cli_exe)//' '//xyz//' -b basis_sets/sto-3g.json '//flags// &
            ' --silent > opt_refusal.out 2>&1'
      call execute_command_line(cmd, exitstat=rc, wait=.true.)

      saw_scf_line = .false.
      saw_message = .false.
      open (newunit=u, file=out_path, status='old', action='read', iostat=ios)
      if (ios == 0) then
         do
            read (u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'SCF              :') > 0) saw_scf_line = .true.
            if (index(line, expect_substr) > 0) saw_message = .true.
         end do
         close (u)
      end if

      if (rc == 0) then
         call test_failed(error, "'"//flags//"' exited 0 -- should have been refused"); return
      end if
      if (saw_scf_line) then
         call test_failed(error, "'"//flags//"' ran an SCF before being refused"); return
      end if
      if (.not. saw_message) then
         call test_failed(error, "'"//flags//"' was refused, but the output did not contain '"// &
                          expect_substr//"' -- it failed for a different reason than the one this case tests")
      end if
   end subroutine refuse

   !
   ! Whether `substr` appears on any line of `path` -- used where a case
   ! must assert a message's ABSENCE (test_max_updates_zero_accepted,
   ! test_seed_with_dimer_accepted).
   !
   function scan_for(path, substr) result(found)
      character(len=*), intent(in) :: path, substr
      logical :: found
      integer :: u, ios
      character(len=1024) :: line
      found = .false.
      open (newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, substr) > 0) found = .true.
      end do
      close (u)
   end function scan_for

end module check_opt_refusals_suite

program check_opt_refusals
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_opt_refusals_suite, only: collect_opt_refusals, cli_exe, work_dir
   implicit none
   character(len=1024) :: arg
   integer :: stat
   call get_command_argument(1, arg); cli_exe = trim(arg)
   call get_command_argument(2, arg); work_dir = trim(arg)
   if (len(cli_exe) == 0 .or. len(work_dir) == 0) then
      write (error_unit, "(a)") "check_opt_refusals: usage: check_opt_refusals <terco CLI> <scratch dir>"
      error stop 1
   end if
   stat = 0
   call run_testsuite(collect_opt_refusals, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_opt_refusals
