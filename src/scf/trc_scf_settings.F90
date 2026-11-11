!
! trc_scf_options_t (the public knobs), trc_scf_result_t (what a run hands
! back) and scf_settings_t (those knobs turned into what the loop actually
! reads: a functional, a spin count, damping/DIIS thresholds, env overrides).
! `resolve` runs once, on the host, before anything is allocated or put on
! the device, and is the only place the run's rejections live (unknown or
! incomplete functional, --df's own restrictions) -- every other setup
! failure comes later, once real device state exists to release
! (trc_scf_context).
!
! WHY A SHARED SPIN RULE. `trc_run_scf`/`trc_load_scf` (trc_capi.F90) need
! nspin before trc_scf_run is even reachable -- to size the guess density and
! to validate a density handed in from outside -- so they used to carry their
! own copy of the nalpha/nbeta -> nspin rule, missing the `frac_occ` case
! this module's rule already resolves. `scf_settings_spin_rule` is that rule, pulled
! out so there is exactly one place it can drift from this module's own use
! of it.
!
module trc_scf_settings
   use pic_types, only: dp
   use pic_mpi_lib, only: comm_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_xc, only: XC_RHO_TOL
   use trc_pcm, only: pcm_model_from_name, PCM_MODEL_NONE, PCM_DEFAULT_LEBEDEV
   use trc_pcm_radii, only: PCM_DEFAULT_RADII_SCALE
   use trc_decontract, only: trc_decontract_options_t
   implicit none
   integer, parameter :: ROTATION_SIZE = 9
   !! A 3x3 rotation stored in column order (see the rotation components below).
   private

   public :: trc_scf_options_t, trc_scf_guess_t, trc_scf_result_t, trc_scf_timings_t, scf_settings_t, scf_settings_spin_rule
   public :: scf_energies_t, scf_settings_eri_thresh_env, ROTATION_SIZE
   public :: scf_settings_homo_lumo_gap
   public :: scf_settings_fock_densities
   ! incr_env/incr_reset_env/incr_real_env stay
   ! private: nothing outside this module uses them.
   public :: DEFAULT_NLC_GRID_LEVEL
   public :: DEFAULT_SOSCF_START

   real(dp), parameter :: DEFAULT_SOSCF_START = 1.0e-2_dp
   !! `trc_scf_options_t%soscf_start`'s own backend default, read here by
   !! `resolve` and by `trc_capi.F90`'s `trc_set_soscf` so the sentinel (0)
   !! means the same number everywhere it is spelled.

   integer, parameter :: DEFAULT_NLC_GRID_LEVEL = 1
   !! VV10's own default grid level, coarser than the semilocal grid's
   !! own default (`trc_dft_grid`'s `DEFAULT_GRID_LEVEL = 3`) by
   !! convention (SG-1-like): VV10's double sum over grid points costs
   !! the square of this grid's size. Both `trc_scf_options_t%nlc_grid_level`
   !! below and `trc_scf_grad%trc_scf_gradient`'s own fallback (its
   !! `nlc_grid_level` argument absent) read this SAME constant, so the
   !! two cannot drift apart the way an SCF energy on one grid and its
   !! own gradient on another once did (feat/grad-vv10's own fix, see
   !! check_vv10_grad_fixedd_fd.F90).

   type :: trc_scf_options_t
      character(len=16) :: functional = ""   !! empty is Hartree-Fock
      character(len=8) :: guess = "gwh"
      !! Starting orbitals when no density is given: "gwh" (Wolfsberg-Helmholz
      !! over the overlap and the core Hamiltonian, the default) or "core".
      !! The core guess is hopeless past a few atoms -- cholesterol starts
      !! 660 hartree high and never recovers -- and GWH costs nothing.
      logical :: unrestricted = .false.      !! forced; an open shell is unrestricted regardless
      real(dp) :: conv_energy = 1.0e-10_dp
      real(dp) :: conv_density = 1.0e-7_dp   !! RMS change per element, printed; no longer gates
      real(dp) :: conv_diis = 1.0e-5_dp
      !! Frobenius norm of the orthogonalised commutator X^T (FDS - SDF) X, the DIIS
      !! error, below which the SCF is converged once the energy change is too. The
      !! density change used to gate and on cholesterol/PBE sat at the quadrature and
      !! screening noise floor: energy converged by iteration 40, then 37 more
      !! iterations waiting on an RMS change wandering between 1e-7 and 1e-6. PySCF
      !! stops on the gradient norm at sqrt(conv_tol) for the same reason.
      logical :: frac_occ = .false.
      !! FRACTIONAL OCCUPATIONS, for the free-atom runs behind the SAD guess. Each
      !! iteration fills a spin channel by aufbau, the degenerate frontier sharing
      !! what is left evenly, so the atom comes out spherical. Restricted fills
      !! `nelec_frac` electrons two per orbital; unrestricted fills `nalpha` and
      !! `nbeta` one per orbital in their own channels (carbon: two alpha thirds in
      !! each p, no beta p). `nalpha`/`nbeta` are ignored for the count only in the
      !! restricted case.
      real(dp) :: nelec_frac = 0.0_dp
      integer :: max_iteration = 100
      integer :: batch_res = 1
      !! Shell-pair batch resolution: buckets per decade of Schwarz bound
      !! when the pairs are binned. 1 is decades, which is what this always
      !! did; higher values make the pre-launch screen finer, at more and
      !! smaller batches. Settable through trc_set_batching.
      type(trc_decontract_options_t) :: decontract
      !! How the four-centre engine splits the basis into its working basis
      !! (trc_decontract). The integrals are exact in any working basis; the
      !! policy sets the cost and, at the screening threshold, what is
      !! dropped (trc_decontract's DECON_SHELL_COST_DEFAULT says how much).
      !! Gradients, the response solves and the RI-MP2 gradient build their
      !! own engines with it too. TRC_DECON_CUTOFF, TRC_DECON_DIFFUSE,
      !! TRC_DECON_SHELL_COST and TRC_DECON_SCREEN_WHOLE override it; with
      !! more than one rank, rank 0's environment decides.
      logical :: incremental = .true.
      !! Rebuild the two-electron matrix from the DENSITY DIFFERENCE, accumulating
      !! onto the last one. G is linear in D, so this is exact but for the
      !! screening, and the screening is the point: the difference shrinks as the
      !! SCF converges while D does not. trc_scf_incremental decides which builds
      !! are increments and carries the error argument; in short, the first build
      !! is full, every later one is an increment, and convergence is declared only
      !! on a full build, so the reported energy and density are a full build's.
      !! The exchange-correlation potential and the continuum are rebuilt in full
      !! from the density every iteration either way. Off for SAD's free atoms and
      !! under `df` (a full build every iteration, nothing to screen);
      !! TRC_INCREMENTAL=0 or 1 overrides it without a rebuild.
      !!
      !! Measured (V100, RHF cc-pVTZ, `eri_thresh` 3.16e-11, `conv_energy` 1e-8),
      !! full builds / these defaults: gly10 14 / 14 iterations, 69.9 / 45.4 s;
      !! gly20 14 / 14, 216.6 / 139.9 s; water20 12 / 12, 22.3 / 13.6 s; final
      !! energies within 2e-12 Eh. The defaults before this schedule (screen 1e-4
      !! tighter, from |FDS - SDF| < 1e-2, a full build every 20) took 69.8 s on
      !! gly10: the first increment came at iteration 9 and cost more than a full
      !! build. `incr_reset` and the three below are that schedule's knobs, each
      !! with an environment hook (TRC_INCR_RESET, TRC_INCR_THRESH_SCALE,
      !! TRC_INCR_EARLY_SCALE, TRC_INCR_EARLY_UNTIL, TRC_INCR_START).
      real(dp) :: incr_predict = 1.0_dp
      !! The build after one whose commutator, times its ratio to the one
      !! before, falls below `incr_predict * conv_diis` is full: the
      !! iteration predicted to converge (trc_scf_incremental, rule 2). 0
      !! turns the prediction off, so convergence is first met on an
      !! increment and retaken on the next, full, build (rule 1);
      !! TRC_INCR_PREDICT overrides it.
      logical :: incremental_requested = .false.
      !! The caller set `incremental` on itself (trc_set_incremental(1),
      !! TRC_INCREMENTAL=1) rather than leaving the default. Read only to
      !! decide how loudly `df` turning increments off is reported: a warning
      !! when they were asked for, an info line when they were the default.
      integer :: incr_reset = 30
      !! At most this many increments in a row before a full build (0: no
      !! limit). A loose cap only: the stall rule rebuilds as soon as the
      !! accumulated error shows (trc_scf_incremental, rule 3), and none of the
      !! measured runs reached 30 increments without one.
      real(dp) :: incr_thresh_scale = 1.0_dp
      !! Multiplies `eri_thresh` for an increment once |FDS - SDF| is below
      !! `incr_early_until`: the full build's own screen. It was 1e-4, chosen on
      !! a 450-atom water cluster in 6-31G where the energy bounced by 1e-4 to
      !! 1e-6 Eh at 1 -- before the Schwarz bounds were fixed (they were 4 pi too
      !! small). Re-measured on w172/6-31G (516 atoms, 2236 functions, 1e-12) at
      !! 1: 12 iterations as with full builds, the energy within 6e-11 Eh, no
      !! bounce. A tighter screen only costs time (gly10 cc-pVTZ, increments from
      !! iteration 2, before the exchange-only density screen: 57.0 s at 1,
      !! 72.5 s at 1e-2, 86.6 s at 1e-4, the same 14 iterations).
      real(dp) :: incr_early_scale = 100.0_dp
      !! Multiplies `eri_thresh` for an increment while the previous
      !! |FDS - SDF| is at least `incr_early_until`; leaving that phase takes a
      !! full build, which discards what the loose increments left in G. At 100
      !! the increments' accumulated error holds the commutator at 2e-5 (water20)
      !! to 5e-5 (gly10, gly20) when they run to convergence, so above 1e-2 it is
      !! invisible: the iterates match full builds to three digits in
      !! |FDS - SDF|. At 1000 that floor reaches ~5e-4 and gly20 took 16 iterations
      !! instead of 14.
      real(dp) :: incr_early_until = 1.0e-2_dp
      !! The commutator norm that ends the loose phase (above). 1e-3 kept gly20
      !! loose too long (the switch rebuild moved the commutator 10%, one more
      !! iteration); 1e-2 matched the full-build iteration counts everywhere.
      real(dp) :: incr_start = huge(1.0_dp)
      !! Do not increment until the previous commutator norm is below this; the
      !! default increments from iteration 2. The old 1e-2 held the increments
      !! back to iterations 7 to 9, so the first five builds were all full.
      integer :: ndiis = 10
      integer :: ndamp = -1                  !! damped steps before DIIS, at least
      !! The damped start, -1 meaning "by spin": a closed shell goes to DIIS from
      !! the first iteration (damping took cholesterol from 29 iterations to 51 and
      !! water from 10 to 19), while an open shell keeps four damped steps and a 0.5
      !! start, since UHF on OH from GWH oscillates for a hundred iterations without
      !! them. Any value >= 0 here, or a finite `diis_start`, is taken as given.
      real(dp) :: diis_start = -1.0_dp     !! commutator norm below which DIIS may start; -1 by spin
      !! DIIS also waits until the largest commutator element is below this;
      !! extrapolating from a wild guess wanders, damping does not.
      real(dp) :: damp = 0.3_dp              !! weight of the new density while damping
      logical :: adiis = .false.
      !! ADIIS (Hu & Yang, JCP 132, 054109) while damping has ended but the
      !! commutator is still too large for CDIIS (see `diis_start`): minimise an
      !! energy functional over the convex hull of the Fock/density history instead
      !! of the level-shifted Roothaan step. Off by default: `ndamp`/`level_shift`
      !! already reach the CDIIS window in the same iteration count on the existing
      !! suite. It earns its keep on hard cases (stretched bonds, small gaps,
      !! charged open shells); see the robustness test set.
      logical :: soscf = .false.
      !! Second-order SCF (`--soscf`): once Pulay DIIS brings the commutator
      !! below `soscf_start`, hand over to a trust-region Newton step on the
      !! orbital rotations (`trc_scf_soscf`), ported from metalquicha's
      !! `mqc_czt_soscf`. Opt in, because each Newton step costs a Krylov
      !! subspace of extra Fock builds (one per Hessian-vector product) on
      !! top of one more to try each backtracked step -- cheaper in
      !! iterations, not always cheaper in wall time. Restricted
      !! Hartree-Fock only, like the backend it was ported from: refused
      !! with an open shell, a functional, `--pcm` or `--df` (`resolve`,
      !! below), never run silently as plain DIIS.
      real(dp) :: soscf_start = 0.0_dp
      !! The commutator at which DIIS hands over to Newton; 0 means the
      !! backend's own default (1e-2, metalquicha's `DEFAULT_SOSCF_START`,
      !! the same scale as `ACCEL_SWITCH` there).
      real(dp) :: level_shift = 1.0_dp
      !! Level shift on the virtual space while DIIS is off, in hartree.
      !! From a Wolfsberg-Helmholz guess cholesterol sits 20 hartree above
      !! its answer and a damped Roothaan step alone throws it 300 hartree
      !! away; raising the virtuals keeps the occupied space from collapsing.
      real(dp) :: lindep_thresh = 1.0e-8_dp
      !! Canonical orthogonalisation: an eigenvalue of S, normalised to a unit
      !! diagonal first (Sn = D S D, D = diag(1/sqrt(S_ii)) -- pyscf's
      !! `remove_linear_dep_`/ORCA's own convention; S itself is not
      !! unit-diagonal, every spherical d function has S_ii = 5/(4*pi)), below
      !! this is linear dependence, projected out of X rather than scaled by
      !! Sn^-1/4 (which would divide by nearly zero). Diffuse bases
      !! (def2-TZVPD, aug-cc-pVDZ) put several Sn eigenvalues below 1e-8 on
      !! anything but a single atom. 0 restores the old `error stop` when S
      !! is not numerically positive definite.
      real(dp) :: eri_thresh = 1.0e-12_dp
      logical :: df = .false.
      !! Density-fitted J and K (`--df`): ALWAYS fits both, never DF-J with
      !! 4c K. Restricted (closed-shell) reference and a Cartesian orbital
      !! basis only, for now; a range-separated functional needs the erf fit
      !! this milestone does not have. `--incremental` is refused with it --
      !! DF is a full build every iteration, and the increment buys nothing.
      !! A caller that sets this must also hand `trc_scf_run` a `dfjk`
      !! object already fitted on the same geometry (trc_capi.F90 owns one
      !! per context handle; see trc_scf_df.F90's own header on why it
      !! outlives a single run).
      integer :: grid_level = 3
      real(dp) :: grid_rotation(ROTATION_SIZE) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
      !! Rotation of the Lebedev points of the DFT and VV10 grids
      !! (`dft_grid_options_t%angular_rotation`); the identity unless a
      !! symmetric run in a standard frame wants the grid of the user's axes.
      !! Rank 1, the 3x3 matrix in column order: ifx 2025.0 rejects an empty
      !! structure constructor of a type with a rank-2 component default
      !! (error #6595).
      integer :: nlc_grid_level = DEFAULT_NLC_GRID_LEVEL
      !! Grid level for VV10's nonlocal correlation, when the functional has
      !! it. Coarser than the semilocal grid by convention (SG-1-like): VV10
      !! is a double sum over grid points, so its cost goes as the square of
      !! this grid's size.
      integer :: grid_max_pts = 512
      real(dp) :: rho_tol = XC_RHO_TOL
      logical :: verbose = .false.
      logical :: force_cartesian = .false.
      !! Run the Cartesian AO space even when the basis says spherical
      !! (`trc_basis_t%spherical`, from a JSON basis's "gto_spherical"), for
      !! comparison against pyscf's cart=True or terco's Cartesian-only references.
      !! `trc_scf_run` does not read this: it reads `b%spherical`, already resolved
      !! by the caller (trc_capi.F90's `basis_json` applies this option to the
      !! loaded basis).
      character(len=8) :: pcm_model = ""
      !! Polarisable continuum: "" (off, the default), "cpcm" or "iefpcm"
      !! (trc_pcm). The surface, K and the nuclear potential are built once
      !! per run; each iteration solves for the charges at the current
      !! density and adds their full potential to the Fock matrix (both
      !! spins), and q.v/2 to the energy.
      real(dp) :: pcm_eps = 0.0_dp
      !! Solvent dielectric constant; required, > 1, when pcm_model is set
      !! (water is 78.3553). No solvent-name table, so nothing is assumed.
      integer :: pcm_lebedev = PCM_DEFAULT_LEBEDEV
      !! Lebedev points per atomic sphere; one of the orders with a fitted
      !! SWIG exponent (302, pyscf's default, here too)
      real(dp) :: pcm_radii_scale = PCM_DEFAULT_RADII_SCALE
      !! van der Waals radius to cavity radius (1.2, as pyscf and Q-Chem)
      real(dp), allocatable :: pcm_radii_override(:)
      !! Final, already-scaled per-atom cavity radii in bohr, trc_set_molecule's
      !! own atom order -- unallocated (the default) reads the table and
      !! pcm_radii_scale above instead. trc_set_pcm_radii's own route in;
      !! see trc_capi.F90 for why a caller would want this (matching another
      !! backend's own radii convention, rather than terco's).
   end type trc_scf_options_t

   type :: trc_scf_guess_t
      !! The two rarer guesses a caller can hand `trc_scf_run` besides a
      !! density (`dguess`), bundled as one optional argument (`guess`) so the
      !! run's argument list stays within the limit. Plain values (copies);
      !! an unallocated component is a guess not given. `dguess` beats
      !! `hamiltonian` beats `orbitals`, each beating the plain
      !! `trc_scf_options_t%guess` string.
      real(dp), allocatable :: hamiltonian(:, :)
      !! A guess HAMILTONIAN to occupy instead of hcore or the GWH matrix --
      !! how SAP gets in (see trc_sap), which cannot come through `dguess`
      !! because what it produces is a one-electron operator and not a
      !! density. Ignored if `dguess` is also given: a density the caller
      !! already has is more specific than a matrix to make one from.
      real(dp), allocatable :: orbitals(:, :, :)
      !! A PREVIOUS geometry's converged orbitals (SAME molecule/basis, so
      !! the same nao), re-orthonormalised into a guess density at THIS
      !! geometry's own overlap once `scf_context%setup` has built it
      !! (`trc_scf_context%scf_context_project_density`) -- trc_capi.F90's
      !! TRC_GUESS_PROJECTED, lowest priority of the three guess-density
      !! arguments (checked after `dguess` and `hamiltonian`, same as every
      !! other guess kind the plain `trc_scf_options_t%guess` string selects).
   end type trc_scf_guess_t

   type :: scf_energies_t
      !! What the SCF loop computes at each density, gathered so `trc_scf_run`
      !! hands it to the result in one `update` call. All zero until set.
      real(dp) :: total = 0.0_dp
      real(dp) :: one = 0.0_dp         !! Tr(D H)
      real(dp) :: two = 0.0_dp         !! Coulomb and exact exchange
      real(dp) :: xc = 0.0_dp
      real(dp) :: nelec_grid = 0.0_dp  !! what the grid integrated the density to
      real(dp) :: pcm = 0.0_dp         !! PCM dielectric energy q.v/2, included in `total`
      real(dp) :: pcm_qtotal = 0.0_dp  !! total apparent surface charge
   end type scf_energies_t

   type :: trc_scf_timings_t
      !! The SCF's own stage times, seconds, read off the loop's pic timers
      !! once the run ends (the Fock timer stops after a synchronous device
      !! kernel; the XC build syncs explicitly). `has_*` say which stages the run had; a stage it did not
      !! have stays 0 and is not reported.
      real(dp) :: wall = 0.0_dp        !! `trc_scf_run` from entry to the last iteration
      real(dp) :: setup = 0.0_dp       !! once per geometry: pairs, 1e, ERI context, linalg
      real(dp) :: fock_total = 0.0_dp  !! all Fock (two-electron) builds of the loop
      real(dp) :: fock_first = 0.0_dp  !! build 1, which carries one-off device and kernel setup
      integer :: n_fock = 0            !! Fock builds taken by the main loop
      real(dp) :: grid = 0.0_dp, xc = 0.0_dp, vv10 = 0.0_dp   !! xc excludes vv10
      real(dp) :: pcm_setup = 0.0_dp, pcm_total = 0.0_dp      !! cavity and K once; the loop's solves summed
      logical :: has_grid = .false., has_vv10 = .false., has_pcm = .false., has_df = .false.
   end type trc_scf_timings_t

   type :: trc_scf_result_t
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: nspin = 1
      real(dp) :: energy = 0.0_dp      !! total
      real(dp) :: e_nuc = 0.0_dp
      real(dp) :: e_one = 0.0_dp       !! Tr(D H)
      real(dp) :: e_two = 0.0_dp       !! Coulomb and exact exchange
      real(dp) :: e_xc = 0.0_dp
      real(dp) :: nelec_grid = 0.0_dp  !! what the grid integrated the density to
      real(dp) :: e_pcm = 0.0_dp       !! PCM dielectric energy q.v/2, included in `energy`
      real(dp) :: pcm_qtotal = 0.0_dp  !! total apparent surface charge
      integer :: pcm_npoints = 0       !! surface points the cavity kept
      type(trc_scf_timings_t) :: timings
      integer :: n_increments = 0
      !! Fock builds taken from the density difference (trc_scf_incremental);
      !! the rest of `iterations` were full builds, the last one always.
      logical :: sym_declined = .false.
      !! The symmetric run was not made: `sym%decline_lindep` and functions were
      !! dropped (`message` says so, `iterations` is 0). The caller runs in C1.
      integer :: n_lindep = 0
      !! Functions dropped by canonical orthogonalisation: S eigenvalues
      !! below `opts%lindep_thresh`. 0 for every compact basis; see
      !! `lindep_thresh` for why a diffuse one is different.
      real(dp), allocatable :: dmat(:, :, :)   !! (nao, nao, nspin); restricted holds the total
      real(dp), allocatable :: cmo(:, :, :)    !! (nao, nao, nspin) orbitals
      real(dp), allocatable :: eps(:, :)       !! (nao, nspin) orbital energies
      integer, allocatable :: mo_irrep(:, :)   !! (nao, nspin) irrep of each orbital (trc_sym_scf); unallocated without symmetry
      real(dp), allocatable :: grid_r(:, :)    !! (3, npts) the grid used, for a same-grid comparison
      real(dp), allocatable :: grid_w(:)
      real(dp), allocatable :: grid_r_nlc(:, :)  !! VV10's own (coarser) grid, unallocated when the
      real(dp), allocatable :: grid_w_nlc(:)     !! functional has no nonlocal correlation.
      character(len=256) :: message = ""
   contains
      procedure :: update => scf_result_update
   end type trc_scf_result_t

   type :: scf_settings_t
      !! Everything `trc_scf_options_t` resolves to before the loop starts:
      !! immutable once `resolve` returns `ok`. Two kinds of field --
      !! policy read straight off `opts` with a functional lookup or a
      !! by-spin default folded in, and the env-hook overrides
      !! (`incr_env` and friends) read once here rather than once per
      !! iteration asking the same question.
      type(trc_xc_functional_t) :: func
      logical :: dft = .false.
      logical :: nlc = .false.
      real(dp) :: exx = 1.0_dp
      real(dp) :: exx_lr = 0.0_dp
      real(dp) :: omega = 0.0_dp
      integer :: nspin = 1
      integer :: nocc(2) = 0
      real(dp) :: occ = 2.0_dp
      integer :: ndamp = 0
      real(dp) :: diis_start = huge(1.0_dp)
      logical :: soscf = .false.
      real(dp) :: soscf_start = DEFAULT_SOSCF_START
      real(dp) :: conv_energy = 0.0_dp
      real(dp) :: conv_diis = 0.0_dp
      integer :: max_iteration = 0
      !! `opts%conv_energy`/`%conv_diis`/`%max_iteration`, resolved through
      !! so `trc_scf_soscf%trc_scf_soscf_run` reads its convergence test and
      !! remaining-iteration budget off `scf_config` alone, the same bundle
      !! every other per-iteration step takes, rather than `scf_options` as
      !! an extra argument.
      logical :: pcm_on = .false.
      logical :: incr_on = .false.
      logical :: incr_off_df = .false.
      !! Increments were on (by default or asked for) but `df` turned them off.
      logical :: incr_df_requested = .false.
      !! ... and they had been asked for explicitly (`incremental_requested`).
      logical :: df = .false.
      logical :: talk = .false.
      logical :: root = .true.
      real(dp) :: eri_thresh_run = 0.0_dp
      real(dp) :: incr_start_run = 0.0_dp
      real(dp) :: incr_scale_run = 0.0_dp
      real(dp) :: incr_early_scale_run = 0.0_dp
      real(dp) :: incr_early_until_run = 0.0_dp
      real(dp) :: incr_predict_run = 1.0_dp
      !! `opts%incr_predict`, or TRC_INCR_PREDICT: the next build is full once the
      !! commutator it is predicted to reach is below this times `conv_diis`.
      real(dp) :: incr_endgame_run = 0.0_dp
      !! TRC_INCR_ENDGAME, an experiment hook only: every build is full once
      !! |FDS - SDF| < this times `conv_diis` (0, the default: the schedule's
      !! own prediction alone, trc_scf_incremental's rule 2).
      integer :: incr_reset_run = 0
      real(dp) :: rho_tol = 0.0_dp   !! the XC batch screen, `opts%rho_tol`
   contains
      procedure :: resolve => scf_settings_resolve
   end type scf_settings_t

contains

   !
   ! The end-of-run copy `trc_scf_run` used to spell out field by field: the
   ! converged density and the energies. Orbitals, grids and counts are set
   ! where they are produced.
   !
   subroutine scf_result_update(this, energies, dmat)
      class(trc_scf_result_t), intent(inout) :: this
      type(scf_energies_t), intent(in) :: energies
      real(dp), intent(in) :: dmat(:, :, :)
      this%dmat = dmat
      this%energy = energies%total
      this%e_one = energies%one
      this%e_two = energies%two
      this%e_xc = energies%xc
      this%nelec_grid = energies%nelec_grid
      this%e_pcm = energies%pcm
      this%pcm_qtotal = energies%pcm_qtotal
   end subroutine scf_result_update

   !
   ! Everything `trc_scf_run` used to work out between its declarations and
   ! its first allocation, in the same order: the functional (and its
   ! rejections), the spin rule, the by-spin DIIS/damping defaults, whether
   ! PCM or the incremental build are on, and the env-hook overrides. `ok` is
   ! false only for the functional and --df rejections below; every other
   ! setup failure (grid, PCM) needs device state built first and is
   ! `trc_scf_context`'s to report.
   !
   function scf_settings_resolve(this, opts, nalpha, nbeta, comm, msg) result(ok)
      class(scf_settings_t), intent(out) :: this
      type(trc_scf_options_t), intent(in) :: opts
      integer, intent(in) :: nalpha, nbeta
      type(comm_t), intent(in), optional :: comm
      character(len=*), intent(out) :: msg
      logical :: ok
      ok = .true.
      msg = ""
      this%talk = opts%verbose
      this%rho_tol = opts%rho_tol
      this%root = .true.
      if (present(comm)) this%root = comm%rank() == 0
      this%talk = this%talk .and. this%root
      this%dft = len_trim(opts%functional) > 0
      this%nlc = .false.
      this%exx = 1.0_dp
      this%exx_lr = 0.0_dp
      this%omega = 0.0_dp
      if (this%dft) then
         if (.not. xc_functional_by_name(opts%functional, this%func)) then
            msg = "trc_scf: unknown functional "//trim(opts%functional)
            ok = .false.
            return
         end if
         ! A registered functional terco cannot evaluate yet (meta-GGA,
         ! dispersion, ...) comes back with `missing` set rather than as an
         ! unknown name; running it would converge on whatever subset of its
         ! pieces is present.
         if (len_trim(this%func%missing) > 0) then
            msg = "trc_scf: functional "//trim(opts%functional)// &
                  " is not available: "//trim(this%func%missing)
            ok = .false.
            return
         end if
         this%exx = this%func%exx
         this%nlc = this%func%nlc_b /= 0.0_dp
         ! VV10 over MPI: trc_vv10's vv10_drive splits pass 2's outer grid
         ! points into a contiguous range per rank and allreduces the
         ! result -- no rejection needed here.
         this%exx_lr = this%func%exx_lr
         this%omega = this%func%omega
      end if
      call scf_settings_spin_rule(nalpha, nbeta, opts%unrestricted, opts%frac_occ, opts%nelec_frac, &
                         this%nspin, this%nocc)
      this%occ = merge(2.0_dp, 1.0_dp, this%nspin == 1)
      this%ndamp = opts%ndamp
      if (this%ndamp < 0) this%ndamp = merge(4, 0, this%nspin == 2)
      this%diis_start = opts%diis_start
      if (this%diis_start < 0.0_dp) this%diis_start = merge(0.5_dp, huge(1.0_dp), this%nspin == 2)
      this%pcm_on = pcm_model_from_name(opts%pcm_model) /= PCM_MODEL_NONE
      ! Not for SAD's free-atom runs (`frac_occ`): a few functions, nothing
      ! to screen, and the guess stays the same density whatever the
      ! molecule's own setting.
      ! Nor for --df, a full build every iteration with no screening for a
      ! difference density to buy; the CLI refuses an explicit --incremental
      ! with it.
      this%incr_on = incr_env(opts%incremental) .and. .not. opts%frac_occ .and. .not. opts%df
      this%incr_off_df = incr_env(opts%incremental) .and. .not. opts%frac_occ .and. opts%df
      this%incr_df_requested = this%incr_off_df .and. (opts%incremental_requested .or. incr_env_set())
      this%df = opts%df
      if (this%df .and. this%nspin == 2) then
         msg = "trc_scf: --df supports a restricted (closed-shell) reference only, for now"
         ok = .false.
         return
      end if
      if (this%df .and. this%exx_lr /= 0.0_dp) then
         msg = "trc_scf: --df with a range-separated functional is refused -- the erf-attenuated "// &
               "fit does not exist yet, and mixing a fitted K with a 4-centre long-range one is "// &
               "never done silently"
         ok = .false.
         return
      end if
      this%soscf = opts%soscf
      this%soscf_start = opts%soscf_start
      if (this%soscf_start <= 0.0_dp) this%soscf_start = DEFAULT_SOSCF_START
      this%conv_energy = opts%conv_energy
      this%conv_diis = opts%conv_diis
      this%max_iteration = opts%max_iteration
      if (this%soscf) then
         ! Restricted Hartree-Fock only, like the backend `trc_scf_soscf`
         ! was ported from (metalquicha's `mqc_czt_soscf`/`mqc_czt_rhf`):
         ! refused here, before anything is built, rather than silently
         ! falling back to plain DIIS -- the same reason `--df` with an
         ! unsupported combination is refused above instead of ignored.
         if (this%nspin == 2) then
            msg = "trc_scf: --soscf supports a restricted (closed-shell) reference only, for now"
            ok = .false.
            return
         end if
         if (this%dft) then
            msg = "trc_scf: --soscf supports Hartree-Fock only, for now -- "// &
                  "the orbital Hessian has no exchange-correlation kernel term here"
            ok = .false.
            return
         end if
         if (this%pcm_on) then
            msg = "trc_scf: --soscf refuses a continuum solvent -- the orbital Hessian has no "// &
                  "PCM response term here"
            ok = .false.
            return
         end if
         if (this%df) then
            msg = "trc_scf: --soscf refuses --df -- its Hessian-vector product runs the four-centre "// &
                  "engine directly and has no density-fitted form"
            ok = .false.
            return
         end if
         if (present(comm)) then
            if (comm%size() > 1) then
               msg = "trc_scf: --soscf runs on one rank only, for now"
               ok = .false.
               return
            end if
         end if
      end if
      this%eri_thresh_run = scf_settings_eri_thresh_env(opts%eri_thresh)
      this%incr_start_run = incr_real_env("TRC_INCR_START", opts%incr_start)
      this%incr_reset_run = incr_reset_env(opts%incr_reset)
      this%incr_scale_run = incr_real_env("TRC_INCR_THRESH_SCALE", opts%incr_thresh_scale)
      this%incr_early_scale_run = incr_real_env("TRC_INCR_EARLY_SCALE", opts%incr_early_scale)
      this%incr_early_until_run = incr_real_env("TRC_INCR_EARLY_UNTIL", opts%incr_early_until)
      this%incr_endgame_run = incr_real_env("TRC_INCR_ENDGAME", 0.0_dp)
      this%incr_predict_run = incr_real_env("TRC_INCR_PREDICT", opts%incr_predict)
   end function scf_settings_resolve

   !
   ! nalpha/nbeta (+ whether the case is forced unrestricted, and whether
   ! this is one of SAD's fractional-occupation atomic runs) -> nspin and
   ! the per-spin occupation count the guess and the loop both fill to.
   ! `trc_scf_run` and `trc_capi.F90`'s `trc_run_scf`/`trc_load_scf` all need
   ! this before a density exists to look at; one rule, called from both.
   !
   pure subroutine scf_settings_spin_rule(nalpha, nbeta, unrestricted, frac_occ, nelec_frac, nspin, nocc)
      integer, intent(in) :: nalpha, nbeta
      logical, intent(in) :: unrestricted, frac_occ
      real(dp), intent(in) :: nelec_frac
      integer, intent(out) :: nspin
      integer, intent(out) :: nocc(2)
      nspin = 1
      if ((nalpha /= nbeta .or. unrestricted) .and. .not. (frac_occ .and. .not. unrestricted)) nspin = 2
      nocc = [nalpha, nbeta]
      ! The guess fills whole orbitals at `nocc`; the loop refills by aufbau
      ! from `nelec_frac` every iteration (`trc_scf_fock`'s orbital step).
      if (frac_occ .and. nspin == 1) nocc = int(nelec_frac/2.0_dp)
   end subroutine scf_settings_spin_rule

   logical function incr_env(dflt) result(res)
      !! Experiment hook for the incremental build, so its effect can be
      !! measured through a driver that does not expose the option:
      !! TRC_INCREMENTAL=0 turns it off, 1 on, anything else leaves `dflt`.
      logical, intent(in) :: dflt
      character(len=8) :: v
      call get_environment_variable("TRC_INCREMENTAL", v)
      res = dflt
      if (trim(v) == "1") res = .true.
      if (trim(v) == "0") res = .false.
   end function incr_env

   logical function incr_env_set() result(res)
      !! TRC_INCREMENTAL=1: increments asked for through the environment.
      character(len=8) :: v
      call get_environment_variable("TRC_INCREMENTAL", v)
      res = (trim(v) == "1")
   end function incr_env_set

   integer function incr_reset_env(dflt) result(res)
      !! Full rebuild interval, and an integral threshold override. A
      !! difference build carries the same absolute screening error as a full
      !! one, so the error accumulates across increments and DIIS stalls
      !! before the energy converges; the cure is a tighter integral
      !! threshold on the difference builds, and resetting often enough that
      !! what does accumulate stays below it. Both are measured here rather
      !! than assumed.
      integer, intent(in) :: dflt
      character(len=8) :: v
      integer :: ios, k
      call get_environment_variable("TRC_INCR_RESET", v)
      res = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) k
         if (ios == 0 .and. k >= 0) res = k
      end if
   end function incr_reset_env

   real(dp) function incr_real_env(name, dflt) result(res)
      !! Experiment hook, like `incr_reset_env`: TRC_INCR_THRESH_SCALE,
      !! TRC_INCR_EARLY_SCALE, TRC_INCR_EARLY_UNTIL and TRC_INCR_START sweep
      !! their option without a rebuild. A positive number replaces `dflt`.
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable(name, v)
      res = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) res = t
      end if
   end function incr_real_env

   real(dp) function scf_settings_eri_thresh_env(dflt) result(res)
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_ERI_THRESH", v)
      res = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) res = t
      end if
   end function scf_settings_eri_thresh_env

   pure function scf_settings_homo_lumo_gap(nao, nspin, nmo, nocc, eps) result(gap)
      !! The smallest HOMO-LUMO gap over the spins, each spin's own LUMO
      !! minus its own HOMO, so a radical's beta LUMO below its alpha HOMO is
      !! not a negative gap. Only the first `nmo` levels are real: with
      !! functions dropped by canonical orthogonalisation the trailing
      !! nao - nmo are the lifted null states. huge() when no spin has both a
      !! HOMO and a LUMO.
      integer, intent(in) :: nao, nspin, nmo
      integer, intent(in) :: nocc(2)
      real(dp), intent(in) :: eps(nao, nspin)
      real(dp) :: gap
      integer :: s
      gap = huge(1.0_dp)
      do s = 1, nspin
         if (nocc(s) < 1 .or. nocc(s) >= nmo) cycle
         gap = min(gap, eps(nocc(s) + 1, s) - eps(nocc(s), s))
      end do
   end function scf_settings_homo_lumo_gap

   pure function scf_settings_fock_densities(config) result(ndens)
      !! Densities in the widest four-centre pass of one Fock build: 2 for
      !! an unrestricted run with exact exchange (both spins' K together,
      !! trc_eri's `fock_resident_many_lent`), 1 otherwise. What sizes the
      !! scratch the SCF lends that build (trc_eri's `fock_work_size`).
      type(scf_settings_t), intent(in) :: config
      integer :: ndens
      ndens = 1
      if (config%nspin == 2 .and. (config%exx /= 0.0_dp .or. config%exx_lr /= 0.0_dp)) ndens = 2
   end function scf_settings_fock_densities

end module trc_scf_settings
