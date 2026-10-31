submodule (trc_run) trc_run_cli
   !! Defining and parsing the `terco` command line: every option, its
   !! validation and the input files it names (bodies of `trc_run`).
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_capi, only: TRC_GUESS_CORE, TRC_GUESS_GWH, TRC_GUESS_SAD, TRC_GUESS_SADQ, TRC_GUESS_SAP, TRC_GUESS_MINAO
   use trc_xc_functional, only: XC_AVAILABLE_FUNCTIONALS, trc_xc_functional_t, xc_functional_by_name, &
                                xc_is_double_hybrid
   use trc_app, only: app_read_xyz, app_resolve_basis, app_set_env, app_help_for_topic, app_set_log_level
   use trc_dlfind_bridge, only: trc_dlfind_available, &
                                TRC_OPT_COORDS_CARTESIAN, TRC_OPT_COORDS_HDLC, TRC_OPT_COORDS_HDLC_TC, &
                                TRC_OPT_COORDS_DLC, TRC_OPT_COORDS_DLC_TC, &
                                TRC_OPT_GUESS_DENSITY, TRC_OPT_GUESS_PROJECTED, &
                                TRC_OPT_ALG_SD, TRC_OPT_ALG_CG_AUTO, TRC_OPT_ALG_CG, TRC_OPT_ALG_LBFGS, &
                                TRC_OPT_ALG_PRFO, TRC_OPT_ALG_NR, TRC_OPT_ALG_DAMPED, &
                                TRC_OPT_TARGET_MINIMUM, TRC_OPT_TARGET_SADDLE, &
                                TRC_OPT_SADDLE_PRFO, TRC_OPT_SADDLE_DIMER, &
                                TRC_OPT_HESS_TERCO, TRC_OPT_HESS_DLFIND, &
                                TRC_OPT_HUPD_NONE, TRC_OPT_HUPD_POWELL, TRC_OPT_HUPD_BOFILL, TRC_OPT_HUPD_BFGS, &
                                TRC_OPT_LINE_AUTO, TRC_OPT_LINE_SCALING, TRC_OPT_LINE_ENERGY_TRUST, &
                                TRC_OPT_LINE_GRADIENT_TRUST, &
                                TRC_OPT_ALG_STOCHASTIC, TRC_OPT_ALG_GENETIC, &
                                TRC_OPT_NEB_FREE, TRC_OPT_NEB_PERP, TRC_OPT_NEB_FROZEN, &
                                TRC_OPT_CONS_CARTESIAN
   implicit none

contains

   module subroutine run_define_cli(r)
      !! Declare every option and flag of the `terco` command on `r%cli`.
      type(run_t), intent(inout) :: r
   call r%cli%set_program("terco", "Hartree-Fock and Kohn-Sham DFT, in standard Fortran")
   call r%cli%add_positional("geometry", "molecule in xyz format, Angstrom", required=.false.)
   call r%cli%add_option("basis", "basis set: a Basis Set Exchange JSON file, or a name looked up "// &
                         "in $TERCO_BASIS_PATH and then the installed basis sets (def2-svp)", short="b", &
                         default="def2-svp")
   call r%cli%add_option("chrg", "molecular charge (default: a .CHRG file here, else 0)", short="c")
   call r%cli%add_option("uhf", "number of unpaired electrons (default: a .UHF file here, else 0 "// &
                         "or 1 by electron count)", short="u")
   call r%cli%add_option("func", "hf (the default), or a functional: "//XC_AVAILABLE_FUNCTIONALS, &
                         short="f", default="hf")
   call r%cli%add_option("guess", "sad, sadq, sap, minao, gwh or core (default: sad; sadq for an ion). "// &
                         "minao always starts from neutral atoms, also for an ion; it is the default with --fukui", short="g")
   call r%cli%add_option("charges", "atomic charge analysis: hirshfeld (none by default)")
   call r%cli%add_option("fukui", "Hirshfeld-condensed Fukui indices, by finite difference in the "// &
                         "electron count: hirshfeld (none by default; needs a closed-shell neutral; "// &
                         "the guess defaults to minao for all three states)")
   call r%cli%add_flag("cart", "Cartesian AOs even when the basis is defined spherical")
   call r%cli%add_option("grid", "DFT integration grid level (3)", default="3")
   call r%cli%add_option("nlc-grid", "grid level for VV10 nonlocal correlation, e.g. wb97x-v (1)", default="1")
   call r%cli%add_option("iterations", "maximum SCF iterations (100)", default="100")
   call r%cli%add_option("econv", "energy convergence, Eh (1e-10)", default="1e-10")
   call r%cli%add_option("gconv", "convergence on the DIIS error norm |FDS - SDF| (1e-5)", default="1e-5")
   call r%cli%add_option("diis", "DIIS subspace size, Fock/error pairs kept (10); more helps a stalled SCF", &
                         default="10")
   call r%cli%add_flag("adiis", "ADIIS before CDIIS: for SCFs that stall with |FDS - SDF| stuck near 1e-2 "// &
                       "(small gaps, charged open shells)")
   call r%cli%add_option("level-shift", "virtual level shift (Eh) during the damped start, before DIIS (1.0)")
   call r%cli%add_option("damp", "weight of the new density while damping, (0, 1] (0.3)")
   call r%cli%add_option("damp-steps", "damped steps before DIIS (by spin: 0 closed shell, 4 open shell)")
   call r%cli%add_option("diis-start", "commutator norm below which DIIS may start (by spin)")
   call r%cli%add_flag("soscf", "second-order SCF: hand over from DIIS to a trust-region Newton step on the "// &
                       "orbital rotations once |FDS - SDF| is small enough. Restricted Hartree-Fock only; "// &
                       "refused with an open shell, a functional, --pcm, --df or more than one rank")
   call r%cli%add_option("soscf-start", "commutator norm at which --soscf hands over from DIIS (1e-2)")
   call r%cli%add_option("thresh", "integral screening threshold (1e-12)", default="1e-12")
   call r%cli%add_option("lindep", "linear-dependence threshold: eigenvalues of S normalised to a unit "// &
                         "diagonal (pyscf/ORCA's convention) below it are dropped (canonical orthogonalisation; "// &
                         "1e-8; 0 drops none, though the minao guess keeps its own 1e-8). "// &
                         "Diffuse bases on large molecules may need ~5e-6")
   call r%cli%add_flag("incremental", "build later Fock matrices from the density difference (closed shell)")
   call r%cli%add_flag("df", "density-fitted J and K instead of the four-centre engine (restricted, "// &
                       "Cartesian only; fits both, never DF-J with 4-centre K)")
   call r%cli%add_option("jkfit", "fitting basis for --df (def2-universal-jkfit)", default="def2-universal-jkfit")
   call r%cli%add_option("pcm", "polarisable continuum solvation: cpcm or iefpcm (pyscf's SWIG "// &
                         "cavity, 302 points per sphere, radii 1.2 x Bondi); needs --eps")
   call r%cli%add_option("eps", "solvent dielectric constant for --pcm, e.g. 78.3553 for water")
   call r%cli%add_flag("rimp2", "RI-MP2 correlation energy after the SCF (restricted or "// &
                       "unrestricted; unrestricted refuses more than one rank)")
   call r%cli%add_option("aux", "auxiliary basis for --rimp2 (cc-pvdz-rifit)", default="cc-pvdz-rifit")
   call r%cli%add_option("frozen", "frozen core orbitals for --rimp2; -1, the default, counts them by element", &
                         default="-1")
   call r%cli%add_option("aux-block", "auxiliary functions per three-index block for --rimp2 (default: as many "// &
                         "as keep the block under 2 GB)")
   call r%cli%add_option("save-scf", "write the SCF result (orbitals, density, energies) to this file, "// &
                         "for a later --restart")
   call r%cli%add_option("restart", "read the SCF result from a --save-scf file instead of running the "// &
                         "SCF, e.g. to go straight to --rimp2; the molecule, basis and method must match")
   call r%cli%add_option("json", "also write a JSON summary to this file")
   call r%cli%add_flag("grad", "Print the analytic nuclear gradient (Eh/Bohr) after the SCF -- "// &
                       "RHF/UHF, or KS-DFT for an LDA/GGA --func (global and range-separated "// &
                       "hybrids and VV10 included); with --df, the fitted RHF/RKS gradient; "// &
                       "refused for a meta-GGA, --rimp2, a double hybrid or --pcm; with more than "// &
                       "one rank it is computed on rank 0 alone, with a warning")
   call r%cli%add_flag("hessian", "finite-difference Cartesian Hessian (central differences of the "// &
                       "analytic gradient) after the SCF, then harmonic vibrational analysis and RRHO "// &
                       "thermochemistry; refused wherever --grad is (meta-GGA, --rimp2, a double "// &
                       "hybrid, --pcm), checked before any displacement runs; single rank only for now")
   call r%cli%add_option("step", "--hessian's displacement, Bohr (0.005)", default="0.005")
   call r%cli%add_option("temperature", "temperature for --hessian's thermochemistry, K (298.15)", &
                         default="298.15")
   call r%cli%add_option("pressure", "pressure for --hessian's thermochemistry, atm (1.0)", default="1.0")
   call r%cli%add_option("symmetry-number", "rotational symmetry number for --hessian's thermochemistry "// &
                         "(1; not determined automatically -- a warning names the quantities it affects)", &
                         default="1")
   call r%cli%add_option("masses", "isotope (trc_atomic_mass, the most-abundant isotope) or average "// &
                         "(trc_atomic_mass_average, the isotope-averaged standard atomic weight), for "// &
                         "--hessian's vibrational analysis (isotope)", default="isotope")
   call r%cli%add_flag("opt-freq", "--opt followed by --hessian at the optimised geometry (Gaussian's "// &
                       "Opt Freq); shorthand for giving both, and --opt-*/--hessian options still apply")
   call r%cli%add_flag("opt", "optimise the geometry with DL-FIND (needs a TERCO_ENABLE_DLFIND=ON build; "// &
                       "refused otherwise). Warm-starts the SCF between steps; writes <geometry>_opt.xyz "// &
                       "and <geometry>_traj.xyz into the CURRENT WORKING DIRECTORY (named from the input "// &
                       "geometry's own base name), and an 'optimization' block to --json. "// &
                       "Combined with --hessian, the Hessian runs at the OPTIMISED geometry. If a step's "// &
                       "SCF fails to converge even after one retry, DL-FIND has no recovery: terco writes "// &
                       "the files/JSON above (marked failed, with the last ACCEPTED geometry) and the "// &
                       "process then stops abnormally from inside DL-FIND, not through terco's own exit path")
   call r%cli%add_option("opt-coords", "internal coordinates for --opt: cartesian, hdlc, hdlc-tc, dlc "// &
                         "or dlc-tc (cartesian)", default="cartesian")
   call r%cli%add_option("opt-maxiter", "maximum number of --opt steps (100)", default="100")
   call r%cli%add_option("opt-gtol", "--opt convergence: this run's own verdict is the max Cartesian "// &
                         "gradient component of the last accepted step, Eh/Bohr (4.5e-4) -- NOT "// &
                         "DL-FIND's own five-criterion internal test, which in hdlc/dlc is over internal "// &
                         "coordinates and is not this number", default="4.5e-4")
   call r%cli%add_option("opt-guess", "SCF warm start between --opt steps: density (the previous step's "// &
                         "converged density) or projected (its orbitals re-orthonormalised in the new "// &
                         "overlap) (density)", default="density")
   call r%cli%add_option("opt-algorithm", "--opt's minimiser: sd, cg, cg-auto, lbfgs, nr, damped, prfo, "// &
                         "stochastic or genetic (lbfgs); refused together with --opt-target saddle, "// &
                         "which picks its own algorithm from --opt-saddle-method. stochastic/genetic "// &
                         "take --opt-seed and the --opt-global-* parameters instead of the single-state "// &
                         "controls above", default="lbfgs")
   call r%cli%add_option("opt-target", "minimum or saddle (minimum)", default="minimum")
   call r%cli%add_option("opt-saddle-method", "prfo (uses a Hessian, --opt-init-hessian) or dimer (no "// &
                         "Hessian); needs --opt-target saddle (prfo)", default="prfo")
   call r%cli%add_flag("opt-connect", "after a saddle search (--opt-target saddle), displace along the "// &
                       "imaginary mode both ways and minimise each to a downhill minimum (DL-FIND task "// &
                       "1011); needs --opt-target saddle. Writes <geometry>_ts.xyz, <geometry>_minus.xyz "// &
                       "and <geometry>_plus.xyz in addition to the usual --opt output")
   call r%cli%add_option("opt-init-hessian", "for --opt-target saddle --opt-saddle-method prfo (or "// &
                         "--opt-algorithm nr): terco (terco's own FD Hessian, trc_run_hessian, via "// &
                         "cb_get_hessian) or dlfind (DL-FIND's own two-point FD of the gradient) (terco)", &
                         default="terco")
   call r%cli%add_option("opt-hessian-update", "Hessian update for prfo/nr: none, powell, bofill or bfgs "// &
                         "(DL-FIND's own default, bofill)")
   call r%cli%add_option("opt-max-updates", "maximum Hessian updates before prfo/nr recomputes one "// &
                         "(DL-FIND's own default, 50)")
   call r%cli%add_option("opt-zero-modes", "modes assumed to have zero vibrational frequency, for "// &
                         "prfo/nr's Hessian (nzero); default 6 with --opt-coords cartesian and "// &
                         "--opt-algorithm prfo/nr (minimum or saddle) or --opt-saddle-method dimer (5 "// &
                         "if the start geometry is linear -- otherwise the search climbs a rotation), "// &
                         "DL-FIND's own default (0) otherwise")
   call r%cli%add_option("opt-soft-modes", "eigenmodes of prfo/nr's Hessian below this value (Eh/Bohr^2) "// &
                         "are treated as soft (DL-FIND's own default: none)")
   call r%cli%add_option("opt-max-step", "DL-FIND's own ceiling on a trust-radius step, Bohr -- also the "// &
                         "ceiling the trust radius may grow to across accepted steps (DL-FIND's own "// &
                         "default, 0.5, unless given); raising it lets the optimiser attempt a larger, "// &
                         "possibly-rejected trial sooner")
   call r%cli%add_option("opt-linesearch", "trust radius/line search: auto (DL-FIND's own default, which "// &
                         "depends on --opt-algorithm), scaling, energy-trust or gradient-trust (auto)", &
                         default="auto")
   call r%cli%add_option("opt-scalestep", "step scale factor for --opt-linesearch scaling (DL-FIND's "// &
                         "own default, 1.0)")
   call r%cli%add_option("opt-lbfgs-memory", "L-BFGS history length for --opt-algorithm lbfgs "// &
                         "(DL-FIND's own default, min(3*natm, 50))")
   call r%cli%add_option("opt-etol", "--opt's energy-change convergence criterion, Eh (DL-FIND's own "// &
                         "default, --opt-gtol/450)")
   call r%cli%add_option("opt-print-level", "DL-FIND's own stdout verbosity (printl): 0, 2, 4 or 6 "// &
                         "(follows this process's own logger level: 0 normally, 2 under --verbose, 4 "// &
                         "under --debug)")
   call r%cli%add_option("opt-neb", "nudged elastic band: the PRODUCT geometry (an xyz file, same atoms "// &
                         "and order as the start geometry); needs --opt-target minimum. Writes "// &
                         "<geometry>_neb.xyz, one frame per image")
   call r%cli%add_option("opt-neb-images", "number of NEB images, DL-FIND's own nimage -- the path has "// &
                         "this many minus one points, both endpoints included (DL-FIND's own default, 10)")
   call r%cli%add_option("opt-neb-spring", "NEB spring constant, Eh/Bohr^2 (DL-FIND's own default, 0.01)")
   call r%cli%add_flag("opt-neb-climb", "climbing-image NEB: once the path is close to converged, refine "// &
                       "its highest-energy image toward the saddle point. DL-FIND has no separate on/off "// &
                       "switch for this: giving the flag lowers its own neb_climb_test threshold to the "// &
                       "default that actually spawns the climbing image; leaving it out keeps that "// &
                       "threshold at a value no real run reaches")
   call r%cli%add_option("opt-neb-endpoints", "free, perpendicular (DL-FIND's own default: endpoints move "// &
                         "only perpendicular to the path) or frozen", default="perpendicular")
   call r%cli%add_option("opt-neb-via", "one or more intermediate geometries for --opt-neb's own initial "// &
                         "path (an xyz file, one or more frames concatenated -- the multi-frame "// &
                         "convention <geometry>_traj.xyz already writes; same atoms and order as the "// &
                         "start geometry, checked frame by frame). DL-FIND's own linear interpolation "// &
                         "between the start and product endpoints slides a bridging atom straight "// &
                         "through another for two COLLINEAR endpoints (HCN<->HNC is one); a bent via "// &
                         "point fixes that. Needs --opt-neb")
   call r%cli%add_option("opt-frozen", "1-based atom indices to freeze, e.g. '1-3,7' (DL-FIND's own spec "// &
                         "= -1; refused together with --opt-coords dlc/dlc-tc, which DL-FIND itself "// &
                         "rejects for a frozen atom)")
   call r%cli%add_option("opt-constrain", "internal-coordinate constraints (needs --opt-coords hdlc, "// &
                         "hdlc-tc, dlc or dlc-tc -- refused with cartesian), ';'-separated: "// &
                         "'bond i j', 'angle i j k' (j is the vertex), 'torsion i j k l' (j-k is the "// &
                         "central bond), 'cartesian i' (needs hdlc/hdlc-tc: DLC cannot constrain a "// &
                         "Cartesian component), 'bond-difference i j k' (bond(i,j) - bond(j,k)). The "// &
                         "constrained value is always the START geometry's -- DL-FIND's own public "// &
                         "interface has no way to set a different target")
   ! --opt-algorithm stochastic/genetic (iopt 51/52, dlf_parallel_opt.f90)
   ! and the population/cycle/radius/mutation settings it reads --
   ! refused together with any other --opt-algorithm (parse_opt_global_
   ! options). DL-FIND's own defaults are documented per option below;
   ! app/trc_dlfind_bridge.F90's own header has the full dl-find.f90
   ! table and which ones apply to stochastic (SS) vs genetic (GA) only.
   call r%cli%add_option("opt-seed", "integer seed for --opt-algorithm stochastic/genetic (sets "// &
                         "DLF_PO_SEED) or --opt-saddle-method dimer (seeds the Fortran intrinsic "// &
                         "generator dlf_dimer.f90's own initial direction reads); without it both "// &
                         "default to the wall clock, so two runs never agree")
   call r%cli%add_option("opt-global-pop-size", "population size (DL-FIND's own default, 25)")
   call r%cli%add_option("opt-global-maxcycle", "maximum cycles (DL-FIND's own default, 10000)")
   call r%cli%add_option("opt-global-init-pop-size", "initial population size, genetic only "// &
                         "(DL-FIND's own default, 2x --opt-global-pop-size)")
   call r%cli%add_option("opt-global-nsave", "how many of the best points to keep (DL-FIND's own "// &
                         "default, 10)")
   call r%cli%add_option("opt-global-distribution", "stochastic's initial distribution: 1 uniform "// &
                         "in a box, 2 uniform in a sphere, 3 Gaussian (DL-FIND's own default)")
   call r%cli%add_option("opt-global-reset", "cycles without improvement before a reset, genetic "// &
                         "only (DL-FIND's own default, 500)")
   call r%cli%add_option("opt-global-radius", "initial search radius, Bohr in Cartesians "// &
                         "(DL-FIND's own default, 1.0)")
   call r%cli%add_option("opt-global-contraction", "radius contraction per cycle, stochastic only "// &
                         "(DL-FIND's own default, 0.9)")
   call r%cli%add_option("opt-global-tolerance-r", "radius convergence tolerance, stochastic only "// &
                         "(DL-FIND's own default, 1e-8)")
   call r%cli%add_option("opt-global-tolerance-g", "gradient-norm convergence tolerance "// &
                         "(DL-FIND's own default, 1e-3)")
   call r%cli%add_option("opt-global-mutation-rate", "mutation rate, genetic only (DL-FIND's own "// &
                         "default, 0.15)")
   call r%cli%add_option("opt-global-death-rate", "fraction of the population replaced per cycle, "// &
                         "genetic only (DL-FIND's own default, 0.5)")
   call r%cli%add_option("opt-global-scalefac", "step scale factor, stochastic only (DL-FIND's own "// &
                         "default, 10.0)")
   call r%cli%add_option("logger", "log verbosity: knowledge, error, warning, performance, info, "// &
                         "large_info, verbose or debug (info); also sets DL-FIND's own printl", default="info")
   call r%cli%add_flag("silent", "only the summary: no SCF iterations, timings or guess details", short="s")
   call r%cli%add_flag("version", "print the version and stop")
   end subroutine run_define_cli

   module subroutine run_parse(r)
      !! Parse the command line and read every input; stop on a bad one, or after
      !! `--help`/`--version`. The banner is printed last (rank 0).
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
      logical :: log_ok
   call r%cli%parse(r%err)
   if (r%err%has_error()) call run_die(r, r%err%get_message()//new_line("a")//"try 'terco --help'")
   if (r%cli%help_requested()) then
      ! pic_cli's --help takes no value, so `--help opt` leaves "opt" in the
      ! geometry positional; read it back as the help topic.
      call r%cli%get("geometry", sval, r%err)
      if (.not. allocated(sval) .or. r%err%has_error()) sval = ""
      if (r%rank0) call global_logger%info(app_help_for_topic(r%cli%help_text(), sval))
      call run_quit(r, 0)
   end if
   if (r%cli%is_set("version")) then
      if (r%rank0) call global_logger%info("terco version "//VERSION)
      call run_quit(r, 0)
   end if
   call r%cli%get("logger", sval, r%err)
   call run_check(r, "--logger")
   call app_set_log_level(sval, log_ok)
   if (.not. log_ok) call run_die(r, "--logger must be knowledge, error, warning, performance, info, large_info, "// &
                              "verbose or debug, got '"//trim(sval)//"'")
   call r%cli%get("geometry", r%geom, r%err)
   if (r%err%has_error() .or. len_trim(r%geom) == 0) then
     call run_die(r, "no geometry given: terco molecule.xyz [options]"//new_line("a")//"try 'terco --help'")
   end if

   call run_parse_inputs(r)
   call run_parse_scf(r)
   call run_parse_features(r)
   call run_parse_opt_flags(r)
   call run_parse_paths(r)

   if (r%rank0) call run_banner(r)
   end subroutine run_parse

   module subroutine run_parse_inputs(r)
      !! The geometry, basis, charge, multiplicity, functional, guess and charge-analysis inputs.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
   ! ---- inputs ----------------------------------------------------------
   call app_read_xyz(r%geom, r%natm, r%z, r%xyz, r%title, r%err)
   if (r%err%has_error()) call run_die(r, r%err%get_message())

   call r%cli%get("basis", r%basis_arg, r%err)
   call app_resolve_basis(r%basis_arg, r%basis_path, r%err)
   if (r%err%has_error()) call run_die(r, r%err%get_message())

   r%chrg = run_int_option(r, "chrg", ".CHRG", 0)
   r%nelec = sum(r%z) - r%chrg
   if (r%nelec < 1) call run_die(r, "charge "//to_string(r%chrg)//" leaves "//to_string(r%nelec)//" electrons")
   r%uhf = run_int_option(r, "uhf", ".UHF", mod(r%nelec, 2))
   if (r%uhf < 0 .or. r%uhf > r%nelec .or. mod(r%nelec - r%uhf, 2) /= 0) then
     call run_die(r, "--uhf "//to_string(r%uhf)//" is impossible with "//to_string(r%nelec)//" electrons: the unpaired "// &
               "count must be between 0 and the electron count, with the same parity")
   end if
   r%mult = r%uhf + 1

   call r%cli%get("func", sval, r%err)
   call run_to_lower_into(trim(sval), r%func)
   if (r%func == "hf" .or. r%func == "rhf" .or. r%func == "uhf") r%func = ""

   ! A double hybrid's PT2 part is always RI-MP2 on the converged KS
   ! orbitals (no canonical, non-fitted MP2 path exists in terco at all):
   ! --func b2plyp runs it without needing --rimp2 too, with --aux's usual
   ! default (cc-pvdz-rifit) unless overridden below.
   if (len(r%func) > 0) then
      block
         type(trc_xc_functional_t) :: func_spec
         if (xc_functional_by_name(r%func, func_spec)) then
            r%dh = xc_is_double_hybrid(func_spec)
            r%vv10 = func_spec%nlc_b /= 0.0_dp
         end if
         if (r%dh) then
            r%pt2_cos = real(func_spec%pt2_cos, c_double)
            r%pt2_css = real(func_spec%pt2_css, c_double)
         end if
      end block
   end if

   r%guess = "sad"
   if (r%chrg /= 0) r%guess = "sadq"
   ! --fukui runs an anion and a cation too; minao is what converges the
   ! hard ones (an NA-040-type nitrosamine anion stalls from SAD), and the
   ! ions inherit the neutral's guess kind.
   if (r%cli%is_set("fukui")) r%guess = "minao"
   if (r%cli%is_set("guess")) then
      call r%cli%get("guess", sval, r%err)
      call run_to_lower_into(trim(sval), r%guess)
   end if
   select case (r%guess)
   case ("sad")
   r%guess_kind = TRC_GUESS_SAD
   case ("sadq")
   r%guess_kind = TRC_GUESS_SADQ
   case ("sap")
   r%guess_kind = TRC_GUESS_SAP
   case ("gwh")
   r%guess_kind = TRC_GUESS_GWH
   case ("core", "hcore")
   r%guess_kind = TRC_GUESS_CORE
   case ("minao")
   r%guess_kind = TRC_GUESS_MINAO
   case default
      call run_die(r, "unknown guess '"//r%guess//"': use sad, sadq, sap, minao, gwh or core")
   end select

   r%charges = ""
   if (r%cli%is_set("charges")) then
      call r%cli%get("charges", sval, r%err)
      call run_to_lower_into(trim(sval), r%charges)
      if (r%charges /= "hirshfeld") call run_die(r, "unknown --charges '"//r%charges//"': use hirshfeld")
   end if

   r%fukui = ""
   if (r%cli%is_set("fukui")) then
      call r%cli%get("fukui", sval, r%err)
      call run_to_lower_into(trim(sval), r%fukui)
      if (r%fukui /= "hirshfeld") call run_die(r, "unknown --fukui '"//r%fukui//"': use hirshfeld")
   end if
   end subroutine run_parse_inputs

   module subroutine run_parse_scf(r)
      !! The SCF's numerical options (grid, convergence, DIIS and its aids, screening, frozen core).
      type(run_t), intent(inout) :: r
   call r%cli%get("grid", r%grid, r%err)
   call run_check(r, "--grid")
   call r%cli%get("nlc-grid", r%nlc_grid, r%err)
   call run_check(r, "--nlc-grid")
   call r%cli%get("iterations", r%maxiter, r%err)
   call run_check(r, "--iterations")
   call r%cli%get("econv", r%econv, r%err)
   call run_check(r, "--econv")
   call r%cli%get("gconv", r%gconv, r%err)
   call run_check(r, "--gconv")
   call r%cli%get("diis", r%ndiis, r%err)
   call run_check(r, "--diis")
   if (r%ndiis < 1) call run_die(r, "--diis must be at least 1")
   if (r%cli%is_set("level-shift")) then
      call r%cli%get("level-shift", r%level_shift, r%err)
      call run_check(r, "--level-shift")
      if (.not. (r%level_shift >= 0.0_dp)) call run_die(r, "--level-shift must be >= 0 (0 turns it off)")
   end if
   r%adiis = r%cli%is_set("adiis")
   r%accel_set = r%adiis .or. r%cli%is_set("damp") .or. r%cli%is_set("damp-steps") .or. r%cli%is_set("diis-start")
   if (r%cli%is_set("damp")) then
      call r%cli%get("damp", r%damp, r%err)
      call run_check(r, "--damp")
      if (.not. (r%damp > 0.0_dp .and. r%damp <= 1.0_dp)) call run_die(r, "--damp must be in (0, 1]")
   end if
   if (r%cli%is_set("damp-steps")) then
      call r%cli%get("damp-steps", r%ndamp, r%err)
      call run_check(r, "--damp-steps")
      if (r%ndamp < 0) call run_die(r, "--damp-steps must be >= 0")
   end if
   if (r%cli%is_set("diis-start")) then
      call r%cli%get("diis-start", r%diis_start, r%err)
      call run_check(r, "--diis-start")
      if (.not. (r%diis_start > 0.0_dp)) call run_die(r, "--diis-start must be > 0")
   end if
   r%soscf = r%cli%is_set("soscf")
   if (r%cli%is_set("soscf-start")) then
      call r%cli%get("soscf-start", r%soscf_start, r%err)
      call run_check(r, "--soscf-start")
      if (.not. (r%soscf_start > 0.0_dp)) call run_die(r, "--soscf-start must be > 0")
   end if
   call r%cli%get("thresh", r%thresh, r%err)
   call run_check(r, "--thresh")
   if (r%cli%is_set("lindep")) then
      call r%cli%get("lindep", r%lindep, r%err)
      call run_check(r, "--lindep")
      ! Written so a NaN fails it too.
      if (.not. (r%lindep >= 0.0_dp .and. r%lindep < 1.0_dp)) call run_die(r, "--lindep must be in [0, 1)")
   end if
   call r%cli%get("frozen", r%nfrozen, r%err)
   call run_check(r, "--frozen")
   end subroutine run_parse_scf

   module subroutine run_parse_features(r)
      !! The feature flags: RI-MP2, DF, gradient, Hessian, `--opt`, PCM, and the basis files they need.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
   r%rimp2 = r%cli%is_set("rimp2") .or. r%dh
   r%df = r%cli%is_set("df")
   r%grad = r%cli%is_set("grad")
   r%hessian = r%cli%is_set("hessian") .or. r%cli%is_set("opt-freq")
   r%cart = r%cli%is_set("cart")
   if (r%hessian) then
      call r%cli%get("step", r%hess_step, r%err)
      call run_check(r, "--step")
      if (r%hess_step <= 0.0_dp) call run_die(r, "--step must be positive, got "//to_string(r%hess_step))
      call r%cli%get("temperature", r%temperature_k, r%err)
      call run_check(r, "--temperature")
      if (r%temperature_k <= 0.0_dp) call run_die(r, "--temperature must be positive")
      call r%cli%get("pressure", r%pressure_atm, r%err)
      call run_check(r, "--pressure")
      if (r%pressure_atm <= 0.0_dp) call run_die(r, "--pressure must be positive")
      call r%cli%get("symmetry-number", r%sym_number, r%err)
      call run_check(r, "--symmetry-number")
      if (r%sym_number < 1) call run_die(r, "--symmetry-number must be at least 1")
      if (r%rank0) call global_logger%warning("--hessian: the rotational symmetry number ("// &
         to_string(r%sym_number)//") is not determined automatically -- it affects the rotational "// &
         "entropy and so the Gibbs free energy; pass --symmetry-number if the molecule's point group "// &
         "has one above 1 (water: 2, ammonia: 3, ...)")
      call r%cli%get("masses", sval, r%err)
      call run_check(r, "--masses")
      call run_to_lower_into(trim(sval), r%masses_kind)
      if (r%masses_kind /= "isotope" .and. r%masses_kind /= "average") then
         call run_die(r, "--masses must be isotope or average, got '"//r%masses_kind//"'")
      end if
   end if
   r%opt = r%cli%is_set("opt") .or. r%cli%is_set("opt-freq")
   r%pcm = ""
   r%eps = 0.0_dp
   if (r%cli%is_set("pcm")) then
      call r%cli%get("pcm", sval, r%err)
      call run_to_lower_into(trim(sval), r%pcm)
      if (r%pcm /= "cpcm" .and. r%pcm /= "iefpcm") call run_die(r, "unknown --pcm '"//r%pcm//"': use cpcm or iefpcm")
      if (.not. r%cli%is_set("eps")) call run_die(r, "--pcm needs the solvent's dielectric constant: "// &
                                              "--eps 78.3553 for water; nothing is assumed")
      call r%cli%get("eps", r%eps, r%err)
      call run_check(r, "--eps")
   else if (r%cli%is_set("eps")) then
      call run_die(r, "--eps is the dielectric constant for --pcm, which was not given")
   end if
   if (r%rimp2) then
      ! RI-MP2 (restricted or unrestricted, UHF/UKS) now follows the
      ! basis's own AO space, spherical or Cartesian, same as a plain
      ! SCF -- trc_run_rimp2 folds the MO coefficients into the Cartesian
      ! space its three-centre tensor runs in itself, see trc_capi.F90's
      ! own header. `--cart` is still the explicit override. An open-shell
      ! --rimp2 refuses more than one rank, deeper in (trc_run_rimp2) --
      ! that needs the communicator this parses later, not knowable yet.
      call r%cli%get("aux", sval, r%err)
      call app_resolve_basis(sval, r%aux_path, r%err)
      if (r%err%has_error()) call run_die(r, "--aux: "//r%err%get_message())
   end if
   if (r%df) then
      if (r%cli%is_set("incremental")) call run_die(r, "--df with --incremental is refused -- a "// &
                                               "density-fitted build is a full build every "// &
                                               "iteration, and the increment buys nothing")
      ! `--df` follows the basis's own AO space too, same as --rimp2 above
      ! (trc_scf_df.F90's df_jk_setup only ever touches the Cartesian
      ! d_c/gmat_c trc_scf_fock.F90's scf_fock_to_cart/scf_fock_to_sph fold at the boundary;
      ! it never reads b%spherical itself). `--cart` is still the explicit
      ! override.
      call r%cli%get("jkfit", sval, r%err)
      call app_resolve_basis(sval, r%jkfit_path, r%err)
      if (r%err%has_error()) call run_die(r, "--jkfit: "//r%err%get_message())
   end if
   end subroutine run_parse_features

   module subroutine run_parse_opt_flags(r)
      !! Refuse an `--opt-*` option given without `--opt`, then parse them all when `--opt` is on.
      type(run_t), intent(inout) :: r
   ! Every --opt-* option is parsed inside `if (r%opt) then` below, so one
   ! given without --opt itself would otherwise do nothing silently --
   ! refused here instead, before that block, so a command line that sets
   ! one always means what it says.
   if (.not. r%opt) then
      block
         character(len=*), parameter :: opt_options(*) = [character(len=18) :: &
            "opt-coords", "opt-guess", "opt-maxiter", "opt-gtol", "opt-target", "opt-saddle-method", &
            "opt-connect", "opt-algorithm", "opt-init-hessian", "opt-hessian-update", "opt-max-updates", &
            "opt-zero-modes", "opt-soft-modes", "opt-max-step", "opt-linesearch", "opt-scalestep", &
            "opt-lbfgs-memory", "opt-etol", "opt-print-level"]
         integer :: k
         do k = 1, size(opt_options)
            if (r%cli%is_set(trim(opt_options(k)))) call run_die(r, "--"//trim(opt_options(k))//" needs --opt")
         end do
      end block
   end if
   ! Validated LAST among the feature flags above: opt_capability_check
   ! reads r%func/r%rimp2/r%dh/r%pcm, every one of which must already be
   ! resolved (not just r%cli%is_set) for the check to see what the run
   ! will actually do -- a PCM run that reached DL-FIND with this check
   ! run too early crashed with `error stop` instead of this message
   ! (this is exactly why opt_capability_check exists: DL-FIND itself
   ! has no graceful way to decline a run it cannot do).
   if (r%opt) then
      call run_parse_opt_basic(r)
      call run_parse_opt_target(r)
      call run_parse_opt_global_options(r)
      call run_parse_opt_neb(r)
      call run_parse_opt_hessian(r)
      call run_parse_opt_steps(r)
      call run_parse_opt_geometry(r)
   end if
   end subroutine run_parse_opt_flags

   module subroutine run_parse_opt_basic(r)
      !! The `--opt` capability checks, coordinates, warm-start guess, step limit and gradient tolerance.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
      character(len=:), allocatable :: opt_refusal
   if (.not. trc_dlfind_available()) call run_die(r, "--opt needs libdlfind: this terco was built with "// &
                                               "-DTERCO_ENABLE_DLFIND=OFF (the default). Reconfigure "// &
                                               "with -DTERCO_ENABLE_DLFIND=ON.")
   if (r%nranks > 1) call run_die(r, "--opt is single-rank only for now (DL-FIND itself is serial; it is "// &
                                "driven from the CLI, never run collectively)")
   if (len(r%restart) > 0) call run_die(r, "--opt cannot start from --restart: restart loads one SCF "// &
                                    "result, and --opt needs to move the geometry and build a new "// &
                                    "one at every step")
   opt_refusal = run_opt_capability_check(r)
   if (len(opt_refusal) > 0) call run_die(r, opt_refusal)
   call r%cli%get("opt-coords", sval, r%err)
   call run_to_lower_into(trim(sval), r%opt_coords_arg)
   select case (r%opt_coords_arg)
   case ("cartesian")
      r%opt_coords_kind = TRC_OPT_COORDS_CARTESIAN
   case ("hdlc")
      r%opt_coords_kind = TRC_OPT_COORDS_HDLC
   case ("hdlc-tc")
      r%opt_coords_kind = TRC_OPT_COORDS_HDLC_TC
   case ("dlc")
      r%opt_coords_kind = TRC_OPT_COORDS_DLC
   case ("dlc-tc")
      r%opt_coords_kind = TRC_OPT_COORDS_DLC_TC
   case default
      call run_die(r, "--opt-coords must be cartesian, hdlc, hdlc-tc, dlc or dlc-tc, got '"//r%opt_coords_arg//"'")
   end select
   call r%cli%get("opt-guess", sval, r%err)
   call run_to_lower_into(trim(sval), r%opt_guess_arg)
   select case (r%opt_guess_arg)
   case ("density")
      r%opt_guess_kind = TRC_OPT_GUESS_DENSITY
   case ("projected")
      r%opt_guess_kind = TRC_OPT_GUESS_PROJECTED
   case default
      call run_die(r, "--opt-guess must be density or projected, got '"//r%opt_guess_arg//"'")
   end select
   call r%cli%get("opt-maxiter", r%opt_maxiter, r%err)
   call run_check(r, "--opt-maxiter")
   if (r%opt_maxiter < 1) call run_die(r, "--opt-maxiter must be at least 1")
   call r%cli%get("opt-gtol", r%opt_gtol, r%err)
   call run_check(r, "--opt-gtol")
   if (r%opt_gtol <= 0.0_dp) call run_die(r, "--opt-gtol must be positive")
   end subroutine run_parse_opt_basic

   module subroutine run_parse_opt_target(r)
      !! `--opt-target`, `--opt-saddle-method`, `--opt-connect` and `--opt-algorithm`.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
   ! ---- target, saddle method, connect -------------------------------
   call r%cli%get("opt-target", sval, r%err)
   call run_to_lower_into(trim(sval), r%opt_target_arg)
   select case (r%opt_target_arg)
   case ("minimum")
      r%opt_target_kind = TRC_OPT_TARGET_MINIMUM
   case ("saddle")
      r%opt_target_kind = TRC_OPT_TARGET_SADDLE
   case default
      call run_die(r, "--opt-target must be minimum or saddle, got '"//r%opt_target_arg//"'")
   end select
   r%opt_connect = r%cli%is_set("opt-connect")
   if (r%opt_target_kind == TRC_OPT_TARGET_MINIMUM) then
      if (r%cli%is_set("opt-saddle-method")) call run_die(r, "--opt-saddle-method needs --opt-target saddle")
      if (r%opt_connect) call run_die(r, "--opt-connect needs --opt-target saddle")
      if (r%cli%is_set("opt-algorithm")) then
         call r%cli%get("opt-algorithm", sval, r%err)
         call run_to_lower_into(trim(sval), r%opt_algorithm_arg)
         select case (r%opt_algorithm_arg)
         case ("sd")
            r%opt_algorithm_kind = TRC_OPT_ALG_SD
         case ("cg-auto")
            r%opt_algorithm_kind = TRC_OPT_ALG_CG_AUTO
         case ("cg")
            r%opt_algorithm_kind = TRC_OPT_ALG_CG
         case ("lbfgs")
            r%opt_algorithm_kind = TRC_OPT_ALG_LBFGS
         case ("prfo")
            r%opt_algorithm_kind = TRC_OPT_ALG_PRFO   ! P-RFO as a plain minimiser: valid, just unusual
         case ("nr")
            r%opt_algorithm_kind = TRC_OPT_ALG_NR
         case ("damped")
            r%opt_algorithm_kind = TRC_OPT_ALG_DAMPED
         case ("stochastic")
            r%opt_algorithm_kind = TRC_OPT_ALG_STOCHASTIC
         case ("genetic")
            r%opt_algorithm_kind = TRC_OPT_ALG_GENETIC
         case default
            call run_die(r, "--opt-algorithm must be sd, cg, cg-auto, lbfgs, nr, damped, prfo, stochastic "// &
                     "or genetic, got '"//r%opt_algorithm_arg//"'")
         end select
      else
         r%opt_algorithm_kind = TRC_OPT_ALG_LBFGS
      end if
      r%opt_is_global = (r%opt_algorithm_kind == TRC_OPT_ALG_STOCHASTIC .or. &
                         r%opt_algorithm_kind == TRC_OPT_ALG_GENETIC)
   else
      ! target = saddle: the saddle method fixes DL-FIND's own
      ! algorithm/coordinate family (cb_get_params) -- --opt-algorithm
      ! is refused rather than silently ignored, so a command line
      ! that sets both never means something different from what it
      ! says.
      if (r%cli%is_set("opt-algorithm")) call run_die(r, "--opt-algorithm is not used with --opt-target "// &
                                                  "saddle: --opt-saddle-method (prfo or dimer) fixes "// &
                                                  "DL-FIND's own algorithm/coordinate family")
      call r%cli%get("opt-saddle-method", sval, r%err)
      call run_to_lower_into(trim(sval), r%opt_saddle_method_arg)
      select case (r%opt_saddle_method_arg)
      case ("prfo")
         r%opt_saddle_method_kind = TRC_OPT_SADDLE_PRFO
      case ("dimer")
         r%opt_saddle_method_kind = TRC_OPT_SADDLE_DIMER
      case default
         call run_die(r, "--opt-saddle-method must be prfo or dimer, got '"//r%opt_saddle_method_arg//"'")
      end select
      ! --opt-connect with --opt-saddle-method prfo used to be refused
      ! here (a bad array copy inside libdlfind's own dlf_formstep.f90,
      ! confirmed with a debug build): cmake/modules/Findlibdlfind.cmake's
      ! pin now carries the fix (3c27977) -- see
      ! app/trc_dlfind_bridge.F90's own header for the detail.
      ! test/check_opt_connect.F90 runs connect with both methods.
   end if
   end subroutine run_parse_opt_target

   module subroutine run_parse_opt_neb(r)
      !! `--opt-neb` and its companions: the product endpoint, images, spring, climbing image, endpoints, via points.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
   ! ---- --opt-neb: nudged elastic band / climbing-image NEB -----------
   r%opt_neb = r%cli%is_set("opt-neb")
   if (r%opt_neb) then
      if (r%opt_target_kind == TRC_OPT_TARGET_SADDLE) then
         call run_die(r, "--opt-neb needs --opt-target minimum: the saddle point is this run's OWN "// &
                  "highest-energy image (--opt-neb-climb), not a separate --opt-target saddle search")
      end if
      if (any(r%opt_algorithm_kind == [TRC_OPT_ALG_PRFO, TRC_OPT_ALG_NR, TRC_OPT_ALG_STOCHASTIC, &
                                       TRC_OPT_ALG_GENETIC])) then
         call run_die(r, "--opt-neb needs a Hessian-free, non-global --opt-algorithm (sd, cg, cg-auto, "// &
                  "lbfgs or damped), got '"//r%opt_algorithm_arg//"'")
      end if
      call r%cli%get("opt-neb", r%opt_neb_endpoint_path, r%err)
      block
         integer :: natm2
         integer, allocatable :: z2(:)
         real(dp), allocatable :: xyz2(:, :)
         character(len=:), allocatable :: title2
         call app_read_xyz(r%opt_neb_endpoint_path, natm2, z2, xyz2, title2, r%err)
         if (r%err%has_error()) call run_die(r, "--opt-neb: "//r%err%get_message())
         call run_check_neb_endpoint(r%natm, r%z, z2, natm2, r%err)
         if (r%err%has_error()) call run_die(r, r%err%get_message())
         r%opt_neb_product_xyz = xyz2
      end block
      if (r%cli%is_set("opt-neb-images")) then
         call r%cli%get("opt-neb-images", r%opt_neb_images, r%err)
         call run_check(r, "--opt-neb-images")
         if (r%opt_neb_images < 3) call run_die(r, "--opt-neb-images must be at least 3 -- DL-FIND's own "// &
                  "minimum (dlf_neb.f90: 'More than 2 images have to be used in NEB')")
      end if
      if (r%cli%is_set("opt-neb-spring")) then
         call r%cli%get("opt-neb-spring", r%opt_neb_spring, r%err)
         call run_check(r, "--opt-neb-spring")
         if (r%opt_neb_spring <= 0.0_dp) call run_die(r, "--opt-neb-spring must be positive")
      end if
      r%opt_neb_climb = r%cli%is_set("opt-neb-climb")
      call r%cli%get("opt-neb-endpoints", sval, r%err)
      call run_to_lower_into(trim(sval), r%opt_neb_endpoints_arg)
      select case (r%opt_neb_endpoints_arg)
      case ("free")
         r%opt_neb_endpoints_kind = TRC_OPT_NEB_FREE
      case ("perpendicular")
         r%opt_neb_endpoints_kind = TRC_OPT_NEB_PERP
      case ("frozen")
         r%opt_neb_endpoints_kind = TRC_OPT_NEB_FROZEN
      case default
         call run_die(r, "--opt-neb-endpoints must be free, perpendicular or frozen, got '"// &
                  r%opt_neb_endpoints_arg//"'")
      end select
      if (r%cli%is_set("opt-neb-via")) then
         call r%cli%get("opt-neb-via", r%opt_neb_via_path, r%err)
         block
            integer :: n_via
            call run_read_xyz_frames(r%opt_neb_via_path, r%natm, r%z, n_via, r%opt_neb_via_xyz, r%err)
            if (r%err%has_error()) call run_die(r, r%err%get_message())
         end block
      end if
   else if (r%cli%is_set("opt-neb-images") .or. r%cli%is_set("opt-neb-spring") .or. &
           r%cli%is_set("opt-neb-climb") .or. r%cli%is_set("opt-neb-endpoints") .or. &
           r%cli%is_set("opt-neb-via")) then
      call run_die(r, "--opt-neb-images/--opt-neb-spring/--opt-neb-climb/--opt-neb-endpoints/--opt-neb-via "// &
               "need --opt-neb (the product endpoint)")
   end if
   end subroutine run_parse_opt_neb

   module subroutine run_parse_opt_hessian(r)
      !! The Hessian controls (prfo/nr/dimer), each refused for an algorithm that never reads it.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
      logical :: opt_hessian_relevant, opt_zero_soft_relevant
   ! ---- Hessian handling (meaningful for prfo/nr; harmless, DL-FIND
   ! simply ignores them, for every other algorithm) -------------------
   call r%cli%get("opt-init-hessian", sval, r%err)
   call run_to_lower_into(trim(sval), r%opt_init_hessian_arg)
   select case (r%opt_init_hessian_arg)
   case ("terco")
      r%opt_hessian_source_kind = TRC_OPT_HESS_TERCO
   case ("dlfind")
      r%opt_hessian_source_kind = TRC_OPT_HESS_DLFIND
   case default
      call run_die(r, "--opt-init-hessian must be terco or dlfind, got '"//r%opt_init_hessian_arg//"'")
   end select
   if (r%cli%is_set("opt-hessian-update")) then
      call r%cli%get("opt-hessian-update", sval, r%err)
      call run_to_lower_into(trim(sval), r%opt_hessian_update_arg)
      select case (r%opt_hessian_update_arg)
      case ("none")
         r%opt_hessian_update_kind = TRC_OPT_HUPD_NONE
      case ("powell")
         r%opt_hessian_update_kind = TRC_OPT_HUPD_POWELL
      case ("bofill")
         r%opt_hessian_update_kind = TRC_OPT_HUPD_BOFILL
      case ("bfgs")
         r%opt_hessian_update_kind = TRC_OPT_HUPD_BFGS
      case default
         call run_die(r, "--opt-hessian-update must be none, powell, bofill or bfgs, got '"// &
                  r%opt_hessian_update_arg//"'")
      end select
   end if
   if (r%cli%is_set("opt-max-updates")) then
      call r%cli%get("opt-max-updates", r%opt_max_updates, r%err)
      call run_check(r, "--opt-max-updates")
      ! 0 is accepted by DL-FIND itself (dlf_formstep.f90's own
      ! iupd>=maxupd test, iupd starting at 0): it recomputes the
      ! Hessian fresh every step rather than ever updating it, not a
      ! degenerate "no Hessian" case.
      if (r%opt_max_updates < 0) call run_die(r, "--opt-max-updates must be at least 0")
   end if
   if (r%cli%is_set("opt-zero-modes")) then
      call r%cli%get("opt-zero-modes", r%opt_zero_modes, r%err)
      call run_check(r, "--opt-zero-modes")
      if (r%opt_zero_modes < 0 .or. r%opt_zero_modes > 3*r%natm) then
         call run_die(r, "--opt-zero-modes must be between 0 and 3*natm ("//to_string(3*r%natm)// &
                  "), got "//to_string(r%opt_zero_modes))
      end if
   end if
   if (r%cli%is_set("opt-soft-modes")) then
      call r%cli%get("opt-soft-modes", r%opt_soft_modes, r%err)
      call run_check(r, "--opt-soft-modes")
      if (r%opt_soft_modes <= 0.0_dp) call run_die(r, "--opt-soft-modes must be positive")
   end if

   ! --opt-init-hessian/--opt-hessian-update/--opt-max-updates:
   ! meaningful only for P-RFO or NR (cb_get_params' own
   ! is_prfo_or_nr branch; DL-FIND never reads inithessian/update/
   ! maxupd for any other algorithm) -- refused rather than silently
   ! ignored, so a command line that sets one means what it says.
   opt_hessian_relevant = (r%opt_target_kind == TRC_OPT_TARGET_SADDLE .and. &
                           r%opt_saddle_method_kind == TRC_OPT_SADDLE_PRFO) .or. &
                          (r%opt_target_kind == TRC_OPT_TARGET_MINIMUM .and. &
                           (r%opt_algorithm_kind == TRC_OPT_ALG_PRFO .or. r%opt_algorithm_kind == TRC_OPT_ALG_NR))
   if (.not. opt_hessian_relevant) then
      if (r%cli%is_set("opt-init-hessian")) call run_die(r, "--opt-init-hessian needs --opt-algorithm prfo/nr "// &
               "or --opt-target saddle --opt-saddle-method prfo")
      if (r%cli%is_set("opt-hessian-update")) call run_die(r, "--opt-hessian-update needs --opt-algorithm "// &
               "prfo/nr or --opt-target saddle --opt-saddle-method prfo")
      if (r%cli%is_set("opt-max-updates")) call run_die(r, "--opt-max-updates needs --opt-algorithm prfo/nr "// &
               "or --opt-target saddle --opt-saddle-method prfo")
   end if

   ! --opt-zero-modes/--opt-soft-modes: meaningful for P-RFO, NR or
   ! the dimer method (cb_get_params' own nzero/soft handling) --
   ! refused the same way for every algorithm that never reads them.
   opt_zero_soft_relevant = (r%opt_target_kind == TRC_OPT_TARGET_SADDLE) .or. &
                            (r%opt_target_kind == TRC_OPT_TARGET_MINIMUM .and. &
                             (r%opt_algorithm_kind == TRC_OPT_ALG_PRFO .or. r%opt_algorithm_kind == TRC_OPT_ALG_NR))
   if (.not. opt_zero_soft_relevant) then
      if (r%cli%is_set("opt-zero-modes")) call run_die(r, "--opt-zero-modes needs --opt-algorithm prfo/nr "// &
               "or --opt-target saddle")
      if (r%cli%is_set("opt-soft-modes")) call run_die(r, "--opt-soft-modes needs --opt-algorithm prfo/nr "// &
               "or --opt-target saddle")
   end if

   ! A single atom has no internal degrees of freedom: P-RFO/NR/the
   ! dimer method in Cartesian coordinates default nzero to 6 (5 if
   ! linear, never reached with one atom), which already exceeds
   ! nvar=3 before any eigenvector-following step could even be
   ! attempted -- refused here rather than left to DL-FIND, which has
   ! no graceful way to decline a run it cannot do.
   if (r%natm == 1 .and. r%opt_coords_kind == TRC_OPT_COORDS_CARTESIAN .and. opt_zero_soft_relevant) then
      call run_die(r, "--opt-algorithm prfo/nr and --opt-target saddle need more than one atom in Cartesian "// &
               "coordinates: a single atom has no internal degrees of freedom to search over")
   end if
   if (r%cli%is_set("opt-zero-modes") .and. r%opt_zero_modes > 3*r%natm - 1) then
      call run_die(r, "--opt-zero-modes "//to_string(r%opt_zero_modes)//" leaves fewer than one true "// &
               "degree of freedom for a "//to_string(r%natm)//"-atom molecule (3*natm = "// &
               to_string(3*r%natm)//")")
   end if
   end subroutine run_parse_opt_hessian

   module subroutine run_parse_opt_steps(r)
      !! The step and trust controls: `--opt-max-step`, line search, scale step, L-BFGS memory, energy tolerance, print level.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
   ! ---- step/trust controls -------------------------------------------
   if (r%cli%is_set("opt-max-step")) then
      call r%cli%get("opt-max-step", r%opt_max_step, r%err)
      call run_check(r, "--opt-max-step")
      if (r%opt_max_step <= 0.0_dp) call run_die(r, "--opt-max-step must be positive")
   end if
   call r%cli%get("opt-linesearch", sval, r%err)
   call run_to_lower_into(trim(sval), r%opt_linesearch_arg)
   select case (r%opt_linesearch_arg)
   case ("auto")
      r%opt_linesearch_kind = TRC_OPT_LINE_AUTO
   case ("scaling")
      r%opt_linesearch_kind = TRC_OPT_LINE_SCALING
   case ("energy-trust")
      r%opt_linesearch_kind = TRC_OPT_LINE_ENERGY_TRUST
   case ("gradient-trust")
      r%opt_linesearch_kind = TRC_OPT_LINE_GRADIENT_TRUST
   case default
      call run_die(r, "--opt-linesearch must be auto, scaling, energy-trust or gradient-trust, got '"// &
               r%opt_linesearch_arg//"'")
   end select
   if (r%cli%is_set("opt-scalestep")) then
      call r%cli%get("opt-scalestep", r%opt_scalestep, r%err)
      call run_check(r, "--opt-scalestep")
      if (r%opt_scalestep <= 0.0_dp) call run_die(r, "--opt-scalestep must be positive")
   end if
   if (r%cli%is_set("opt-lbfgs-memory")) then
      call r%cli%get("opt-lbfgs-memory", r%opt_lbfgs_memory, r%err)
      call run_check(r, "--opt-lbfgs-memory")
      if (r%opt_lbfgs_memory < 1) call run_die(r, "--opt-lbfgs-memory must be at least 1")
   end if
   if (r%cli%is_set("opt-etol")) then
      call r%cli%get("opt-etol", r%opt_etol, r%err)
      call run_check(r, "--opt-etol")
      if (r%opt_etol <= 0.0_dp) call run_die(r, "--opt-etol must be positive")
   end if
   if (r%cli%is_set("opt-print-level")) then
      call r%cli%get("opt-print-level", r%opt_print_level, r%err)
      call run_check(r, "--opt-print-level")
      if (all(r%opt_print_level /= [0, 2, 4, 6])) call run_die(r, "--opt-print-level must be 0, 2, 4 or 6")
   end if
   end subroutine run_parse_opt_steps

   module subroutine run_parse_opt_geometry(r)
      !! `--opt-frozen` and `--opt-constrain`.
      type(run_t), intent(inout) :: r
   ! ---- --opt-frozen: atoms = -1 in the residue/spec block -------------
   if (r%cli%is_set("opt-frozen")) then
      if (r%opt_coords_kind == TRC_OPT_COORDS_DLC .or. r%opt_coords_kind == TRC_OPT_COORDS_DLC_TC) then
         call run_die(r, "--opt-frozen is refused together with --opt-coords dlc/dlc-tc: DL-FIND itself "// &
                  "rejects a frozen atom in pure internals (dlf_coords.f90: 'Frozen atoms are not "// &
                  "pure internals')")
      end if
      call r%cli%get("opt-frozen", r%opt_frozen_arg, r%err)
      call run_parse_frozen_list(r%opt_frozen_arg, r%natm, r%opt_frozen_atoms, r%err)
      if (r%err%has_error()) call run_die(r, r%err%get_message())
   end if

   ! ---- --opt-constrain: internal-coordinate constraints ---------------
   if (r%cli%is_set("opt-constrain")) then
      if (r%opt_coords_kind == TRC_OPT_COORDS_CARTESIAN) then
         call run_die(r, "--opt-constrain needs --opt-coords hdlc, hdlc-tc, dlc or dlc-tc: DL-FIND's own "// &
                  "internal-coordinate machinery is what enforces a constraint; Cartesian coordinates "// &
                  "would just ignore it")
      end if
      call r%cli%get("opt-constrain", r%opt_constrain_arg, r%err)
      call run_parse_constraints(r%opt_constrain_arg, r%natm, r%opt_cons_kind, r%opt_cons_a, r%opt_cons_b, &
                             r%opt_cons_c, r%opt_cons_d, r%err)
      if (r%err%has_error()) call run_die(r, r%err%get_message())
      if ((r%opt_coords_kind == TRC_OPT_COORDS_DLC .or. r%opt_coords_kind == TRC_OPT_COORDS_DLC_TC) .and. &
          any(r%opt_cons_kind == TRC_OPT_CONS_CARTESIAN)) then
         call run_die(r, "--opt-constrain 'cartesian': needs --opt-coords hdlc or hdlc-tc -- DL-FIND's own "// &
                  "DLC cannot constrain a Cartesian component (dlf_hdlc_interface.f90: 'no Cartesian "// &
                  "components can be constrained with DLC')")
      end if
   end if
   end subroutine run_parse_opt_geometry

   module subroutine run_parse_paths(r)
      !! Output paths (`--json`, `--save-scf`, `--restart`) and the guess/Hirshfeld data files the library finds by environment.
      type(run_t), intent(inout) :: r
      character(len=:), allocatable :: sval
   r%json = ""
   if (r%cli%is_set("json")) call r%cli%get("json", r%json, r%err)
   r%save_scf = ""
   if (r%cli%is_set("save-scf")) call r%cli%get("save-scf", r%save_scf, r%err)
   r%restart = ""
   if (r%cli%is_set("restart")) call r%cli%get("restart", r%restart, r%err)

   ! SAP's fitted atomic potentials ship as a basis-format file; find it the
   ! way a basis is found, unless the caller has already said where it is.
   if (r%guess_kind == TRC_GUESS_SAP) then
      block
         character(len=1024) :: e
         integer :: n, s
         call get_environment_variable("TRC_SAP_BASIS", e, n, s)
         if (s /= 0 .or. n == 0) then
            call app_resolve_basis("sap_helfem_large", sval, r%err)
            if (r%err%has_error()) call run_die(r, "--guess sap needs its fitted potentials: "//r%err%get_message())
            call app_set_env("TRC_SAP_BASIS", sval)
         end if
      end block
   end if

   ! minao's own minimal-ANO basis ships the same way SAP's fitted
   ! potentials do: found like a basis unless the caller has already set
   ! TRC_MINAO_BASIS.
   if (r%guess_kind == TRC_GUESS_MINAO) then
      block
         character(len=1024) :: e
         integer :: n, s
         call get_environment_variable("TRC_MINAO_BASIS", e, n, s)
         if (s /= 0 .or. n == 0) then
            call app_resolve_basis("minao", sval, r%err)
            if (r%err%has_error()) call run_die(r, "--guess minao needs its basis: "//r%err%get_message())
            call app_set_env("TRC_MINAO_BASIS", sval)
         end if
      end block
   end if

   ! The Hirshfeld free atoms ship next to the basis sets too: without this
   ! the library falls back to a path relative to the working directory,
   ! which only exists inside the test tree.
   if (r%charges == "hirshfeld" .or. r%fukui == "hirshfeld") then
      block
         character(len=1024) :: e
         integer :: n, s
         call get_environment_variable("TRC_HIRSHFELD_DATA", e, n, s)
         if (s /= 0 .or. n == 0) then
            call app_resolve_basis("hirshfeld_proatoms", sval, r%err)
            if (r%err%has_error()) call run_die(r, "Hirshfeld charges need the free-atom table: "//r%err%get_message())
            call app_set_env("TRC_HIRSHFELD_DATA", sval)
         end if
      end block
   end if
   end subroutine run_parse_paths

end submodule trc_run_cli
