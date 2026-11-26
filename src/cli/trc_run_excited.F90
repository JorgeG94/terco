submodule (trc_run) trc_run_excited
   !! The CLI side of TDDFT/TDHF (`--tddft`, `--tddft-*`): the options, what the method refuses before any
   !! work, and the calls into the library (bodies of `trc_run`). The gradient and the optimisation reach the
   !! library through `--grad` and `--opt`, which follow the root `--tddft-root` while `--tddft` is on. The
   !! report and the JSON come from rank 0 as for the rest of the CLI.
   use pic_strings, only: to_string
   use trc_constants, only: HARTREE_TO_EV
   use trc_capi, only: trc_set_tddft, trc_set_tddft_convergence, trc_set_tddft_missed_root_iterations, &
                       trc_set_tddft_density_screen, trc_run_tddft, trc_tddft_nroots, &
                       trc_tddft_window_open, trc_tddft_energies, trc_tddft_character, trc_tddft_energy
   use trc_scf_settings, only: trc_scf_options_t
   use trc_tddft_run, only: trc_tddft_run_options_t, trc_tddft_run_refusal, trc_tddft_run_pair_label
   implicit none

   integer(c_int), parameter :: FOLLOWED_ROOT_BUDGET = 24_c_int
   !! The missed-root budget --grad and --opt use unless the user gives one: they follow a root by index.

   ! An option as the user met it: the flag, and the TOML key that gives the same setting.
   character(len=*), parameter :: N_ROOTS = "--tddft-roots ([tddft] roots)"
   character(len=*), parameter :: N_ROOT = "--tddft-root ([tddft] target_root)"
   character(len=*), parameter :: N_APPROX = "--tddft-approx ([tddft] approximation)"
   character(len=*), parameter :: N_SPIN = "--tddft-spin ([tddft] spin)"
   character(len=*), parameter :: N_TOL = "--tddft-tol ([tddft] energy_tol)"
   character(len=*), parameter :: N_ITER = "--tddft-iterations ([tddft] max_iterations)"
   character(len=*), parameter :: N_MISSED = "--tddft-missed-root-iterations ([tddft] missed_root_iterations)"
   character(len=*), parameter :: N_CPKS = "--tddft-cpks-tol ([tddft] cpks_tol)"
   character(len=*), parameter :: N_SWEEP = "--tddft-sweep ([tddft] roots_per_sweep)"
   character(len=*), parameter :: N_GAP = "--tddft-gap-warn ([tddft] gap_warning)"

contains

   subroutine die_options_alone(r)
      !! The refusal of `--tddft-*` options given without `--tddft`, in a routine of its own (see
      !! trc_run_cc.F90's `die_options_alone` for why a first `run_die` call lives apart).
      type(run_t), intent(in) :: r

      call run_die(r, "the --tddft-* options (the [tddft] table) need --tddft ([tddft] enabled = true)")
   end subroutine die_options_alone

   subroutine parse_values(r)
      !! The numeric options, each checked as it is read.
      type(run_t), intent(inout) :: r
      real(dp) :: window

      if (r%cli%is_set("tddft-roots")) then
         call r%cli%get("tddft-roots", r%tddft_nroots, r%err)
         call run_check(r, N_ROOTS)
         if (r%tddft_nroots < 1) call run_die(r, N_ROOTS//" must be at least 1")
      end if
      if (r%cli%is_set("tddft-root")) then
         call r%cli%get("tddft-root", r%tddft_root, r%err)
         call run_check(r, N_ROOT)
      end if
      if (r%tddft_root < 1 .or. r%tddft_root > r%tddft_nroots) then
         call run_die(r, N_ROOT//" is "//to_string(r%tddft_root)//": it counts from 1 and must not exceed the "// &
                      to_string(r%tddft_nroots)//" roots asked for (--tddft-roots)")
      end if
      if (r%cli%is_set("tddft-tol")) then
         call r%cli%get("tddft-tol", r%tddft_tol, r%err)
         call run_check(r, N_TOL)
         if (.not. (r%tddft_tol > 0.0_dp .and. r%tddft_tol <= huge(1.0_dp))) call run_die(r, N_TOL//" must be positive")
      end if
      if (r%cli%is_set("tddft-iterations")) then
         call r%cli%get("tddft-iterations", r%tddft_iterations, r%err)
         call run_check(r, N_ITER)
         if (r%tddft_iterations < 1) call run_die(r, N_ITER//" must be at least 1")
      end if
      if (r%cli%is_set("tddft-missed-root-iterations")) then
         call r%cli%get("tddft-missed-root-iterations", r%tddft_missed_root, r%err)
         call run_check(r, N_MISSED)
         if (r%tddft_missed_root < 0) call run_die(r, N_MISSED//" must be zero or more")
      end if
      if (r%cli%is_set("tddft-cpks-tol")) then
         call r%cli%get("tddft-cpks-tol", r%tddft_cpks_tol, r%err)
         call run_check(r, N_CPKS)
         if (.not. (r%tddft_cpks_tol > 0.0_dp .and. r%tddft_cpks_tol <= huge(1.0_dp))) then
            call run_die(r, N_CPKS//" must be positive")
         end if
      end if
      if (r%cli%is_set("tddft-sweep")) then
         call r%cli%get("tddft-sweep", r%tddft_sweep, r%err)
         call run_check(r, N_SWEEP)
         if (r%tddft_sweep < 1) call run_die(r, N_SWEEP//" must be at least 1")
      end if
      if (r%cli%is_set("tddft-gap-warn")) then
         call r%cli%get("tddft-gap-warn", window, r%err)
         call run_check(r, N_GAP)
         if (.not. (window >= 0.0_dp .and. window <= huge(1.0_dp))) then
            call run_die(r, N_GAP//" must be a finite number, zero or more (zero: no gap warning)")
         end if
         r%tddft_gap_warn = window
      end if
   end subroutine parse_values

   subroutine parse_words(r)
      !! `--tddft-approx` and `--tddft-spin`: the words, and the refusal of the triplet.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval

      if (r%cli%is_set("tddft-approx")) then
         call r%cli%get("tddft-approx", sval, r%err)
         call run_check(r, N_APPROX)
         sval = run_to_lower_alloc(trim(sval))
         if (sval == "tda") then
            r%tddft_tda = .true.
         else if (sval == "rpa") then
            r%tddft_tda = .false.
         else
            call run_die(r, N_APPROX//" must be tda (CIS for Hartree-Fock) or rpa (TDHF / Casida), got '"//sval//"'")
         end if
      end if
      if (r%cli%is_set("tddft-spin")) then
         call r%cli%get("tddft-spin", sval, r%err)
         call run_check(r, N_SPIN)
         sval = run_to_lower_alloc(trim(sval))
         if (sval == "triplet") then
            call run_die(r, N_SPIN//": triplet excitations are not implemented (only singlets, of a closed-shell reference)")
         else if (sval /= "singlet") then
            call run_die(r, N_SPIN//" must be singlet, got '"//sval//"'")
         end if
      end if
   end subroutine parse_words

   subroutine refuse_combinations(r)
      !! What the method cannot be combined with, from the command line alone.
      type(run_t), intent(inout) :: r
      type(trc_scf_options_t) :: sopts
      type(trc_tddft_run_options_t) :: opts
      character(len=256) :: reason

      if (r%rimp2) call run_die(r, "--tddft cannot be combined with --rimp2 or a double-hybrid functional: "// &
                                "there is no excited-state correlation energy here")
      if (r%cc) call run_die(r, "--tddft cannot be combined with --riccsd or --riccsd-t")
      if (r%uhf /= 0) then
         if (r%cli%is_set("uhf")) then
            call run_die(r, "--tddft is closed-shell only; --uhf "//to_string(r%uhf)//" asks for unpaired electrons "// &
                         "(UHF/UKS excited states are not implemented)")
         else
            call run_die(r, "--tddft is closed-shell only; the "//to_string(r%nelec)//" electrons of this molecule "// &
                         "(charge "//to_string(r%chrg)//") cannot all be paired: check --charge")
         end if
      end if
      if (r%hessian) call run_die(r, "--tddft has no Hessian: it cannot be combined with --hessian or --opt-freq")
      if (r%df) call run_die(r, "--tddft is not implemented over a --df (density-fitted) reference")
      if (len(r%pcm) > 0) call run_die(r, "--tddft has no solvent response: it cannot be combined with --pcm")
      if (r%nranks > 1) call run_die(r, "--tddft runs on one rank (the response path is not rank-split)")
      ! The library's own list, from the functional and the options alone: meta-GGA, range-separated hybrids, VV10,
      ! double hybrids, an unknown functional, a root out of range.
      sopts%functional = r%func
      opts%nroots = r%tddft_nroots
      opts%target_root = r%tddft_root
      reason = trc_tddft_run_refusal(sopts, opts)
      if (len_trim(reason) > 0) call run_die(r, trim(reason))
   end subroutine refuse_combinations

   module subroutine run_tddft_parse(r)
      !! Reads `--tddft` and the `--tddft-*` options. Every combination the method does not cover is refused here,
      !! before the SCF, with the library's own message where it has one: a triplet, an unrestricted reference, a
      !! meta-GGA, range-separated or VV10 functional, --df, --pcm.
      type(run_t), intent(inout) :: r
      logical :: any_option

      r%tddft = r%cli%is_set("tddft")
      any_option = r%cli%is_set("tddft-roots") .or. r%cli%is_set("tddft-root") .or. r%cli%is_set("tddft-approx") .or. &
                   r%cli%is_set("tddft-spin") .or. r%cli%is_set("tddft-tol") .or. r%cli%is_set("tddft-iterations") .or. &
                   r%cli%is_set("tddft-missed-root-iterations") .or. &
                   r%cli%is_set("tddft-cpks-tol") .or. r%cli%is_set("tddft-sweep") .or. &
                   r%cli%is_set("tddft-gap-warn") .or. r%cli%is_set("tddft-strict-root") .or. &
                   r%cli%is_set("no-tddft-grid-response") .or. &
                   r%cli%is_set("tddft-density-screen") .or. r%cli%is_set("no-tddft-density-screen") .or. &
                   r%cli%is_set("tddft-grid-response")
      if (.not. r%tddft) then
         if (any_option) call die_options_alone(r)
         return
      end if
      if (r%cli%is_set("tddft-grid-response") .and. r%cli%is_set("no-tddft-grid-response")) then
         call run_die(r, "--tddft-grid-response and --no-tddft-grid-response contradict each other")
      end if
      r%tddft_grid_response = .not. r%cli%is_set("no-tddft-grid-response")
      if (r%cli%is_set("tddft-density-screen") .and. r%cli%is_set("no-tddft-density-screen")) then
         call run_die(r, "--tddft-density-screen and --no-tddft-density-screen contradict each other")
      end if
      r%tddft_density_screen = .not. r%cli%is_set("no-tddft-density-screen")
      r%tddft_strict = r%cli%is_set("tddft-strict-root")
      call parse_words(r)
      call parse_values(r)
      call refuse_combinations(r)
      if (r%tddft_root > r%tddft_nroots) call run_die(r, N_ROOT//" exceeds --tddft-roots")
      if (r%rank0 .and. len(r%func) == 0 .and. .not. r%tddft_grid_response) then
         call global_logger%warning("--no-tddft-grid-response has no effect for Hartree-Fock")
      end if
   end subroutine run_tddft_parse

   module subroutine run_tddft_configure(r)
      !! The library's TDDFT settings from the parsed options. The solver settings are passed only when one was
      !! given, so the library's defaults (1e-8, 50 iterations, a Z-vector tolerance of 1e-9) stand otherwise.
      type(run_t), intent(in) :: r

      call run_api(r, trc_set_tddft(r%h, int(r%tddft_nroots, c_int), merge(1_c_int, 0_c_int, r%tddft_tda), 1_c_int, &
                                    int(r%tddft_root, c_int), merge(1_c_int, 0_c_int, r%tddft_grid_response)))
      if (r%tddft_tol > 0.0_dp .or. r%tddft_iterations > 0 .or. r%tddft_cpks_tol > 0.0_dp .or. r%tddft_sweep > 0) then
         call run_api(r, trc_set_tddft_convergence(r%h, real(r%tddft_tol, c_double), int(r%tddft_iterations, c_int), &
                                                   real(r%tddft_cpks_tol, c_double), 0_c_int, &
                                                   int(r%tddft_sweep, c_int)))
      end if
      if (r%tddft_missed_root >= 0) then
         call run_api(r, trc_set_tddft_missed_root_iterations(r%h, int(r%tddft_missed_root, c_int)))
      else if (r%grad .or. r%opt) then
         ! --grad and --opt follow a root BY INDEX: a root the solver missed shifts the index to another state, and the
         ! run differentiates or optimises the wrong surface with only a warning. So they keep the missed-root guard
         ! unless the user chose a budget; single points keep the library default (0).
         call run_api(r, trc_set_tddft_missed_root_iterations(r%h, FOLLOWED_ROOT_BUDGET))
         if (r%rank0) call global_logger%info("--tddft with --grad/--opt follows a root by index: the missed-root "// &
                                              "budget is "//to_string(int(FOLLOWED_ROOT_BUDGET))// &
                                              " (--tddft-missed-root-iterations N overrides it)")
      end if
      if (.not. r%tddft_density_screen) call run_api(r, trc_set_tddft_density_screen(r%h, 0_c_int))
   end subroutine run_tddft_configure

   module subroutine run_tddft_execute(r)
      !! Runs the excitation energies on the SCF the handle holds and reads the roots back. A refusal or a solve
      !! that did not converge stops the run with the library's message (`run_api`).
      type(run_t), intent(inout) :: r
      integer(c_int) :: nroots, nopen
      real(c_double), allocatable :: omega(:), oscillator(:), ynorm(:), dipole(:), weight(:)
      integer(c_int), allocatable :: hole(:), particle(:)

      call run_api(r, trc_run_tddft(r%h))
      call run_api(r, trc_tddft_nroots(r%h, nroots))
      allocate (omega(nroots), oscillator(nroots), ynorm(nroots), dipole(3*nroots), weight(nroots))
      allocate (hole(nroots), particle(nroots))
      call run_api(r, trc_tddft_energies(r%h, nroots, omega, oscillator, ynorm, dipole))
      call run_api(r, trc_tddft_character(r%h, nroots, hole, particle, weight))
      call run_api(r, trc_tddft_energy(r%h, r%tddft_etot))
      call run_api(r, trc_tddft_window_open(r%h, nopen))
      r%tddft_n = int(nroots)
      r%tddft_open = int(nopen)
      r%tddft_omega = real(omega, dp)
      r%tddft_osc = real(oscillator, dp)
      r%tddft_ynorm = real(ynorm, dp)
      r%tddft_weight = real(weight, dp)
      r%tddft_hole = int(hole)
      r%tddft_particle = int(particle)
      allocate (r%tddft_dipole(3, nroots))
      r%tddft_dipole = reshape(real(dipole, dp), [3, int(nroots)])
      r%tddft_ran = .true.
      ! The energy --grad and --opt differentiate is the target root's.
      r%etot_total = r%tddft_etot
   end subroutine run_tddft_execute

   module subroutine run_tddft_summary(r)
      !! The excited states, in the output's own table: the energy in hartree and eV, the oscillator strength, ||Y||^2
      !! (zero for the TDA) and the dominant orbital pair with its share of X^2 - Y^2. The root --grad and --opt
      !! follow is marked.
      type(run_t), intent(in) :: r
      character(len=1024) :: logline
      integer :: k

      call global_logger%info("")
      call global_logger%info("   excited states (singlet, "//merge("TDA", "RPA", r%tddft_tda)//", "// &
                              to_string(r%tddft_n)//" roots; the root --grad/--opt follow is marked *)")
      call global_logger%info("   ----------------------------------------------------------------------------------")
      call global_logger%info("      root   omega (Eh)      omega (eV)         f        ||Y||^2   dominant pair (weight)")
      do k = 1, r%tddft_n
         write (logline, "(a,a1,i4,f15.8,f14.6,f12.6,f12.6,3x,a,a,f4.2,a)") "   ", merge("*", " ", k == r%tddft_root), k, &
            r%tddft_omega(k), r%tddft_omega(k)*HARTREE_TO_EV, r%tddft_osc(k), r%tddft_ynorm(k), &
            trc_tddft_run_pair_label(r%tddft_hole(k), r%tddft_particle(k)), "  (", r%tddft_weight(k), ")"
         call global_logger%info(trim(logline))
      end do
      call global_logger%info("   ----------------------------------------------------------------------------------")
      write (logline, "(a,i0,a,f18.12,a)") "   total energy of root ", r%tddft_root, ": ", real(r%tddft_etot, dp), " Eh"
      call global_logger%info(trim(logline))
   end subroutine run_tddft_summary

end submodule trc_run_excited
