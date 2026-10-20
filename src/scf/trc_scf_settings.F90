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
! this module's rule already resolves. `scf_spin_rule` is that rule, pulled
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
   implicit none
   private

   public :: trc_scf_options_t, trc_scf_result_t, scf_settings_t, scf_spin_rule, scf_energies_t
   public :: incr_env, incr_reset_env, incr_thresh_scale_env, incr_start_env, eri_thresh_env
   public :: DEFAULT_NLC_GRID_LEVEL

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
      integer :: max_iter = 100
      integer :: batch_res = 1
      !! Shell-pair batch resolution: buckets per decade of Schwarz bound
      !! when the pairs are binned. 1 is decades, which is what this always
      !! did; higher values make the pre-launch screen finer, at more and
      !! smaller batches. Settable through trc_set_batching.
      logical :: incremental = .false.
      !! Rebuild the two-electron matrix from the DENSITY DIFFERENCE, accumulating
      !! onto the last one. G is linear in D so this is exact; it buys screening,
      !! because the difference shrinks as the SCF converges while D does not. On a
      !! 123-atom silica slice the late iterations dominate the time and are the ones
      !! the difference screens hardest. A full build is taken every `incr_reset`
      !! iterations as a drift backstop. TRC_INCREMENTAL=1 turns it on without a
      !! rebuild.
      !!
      !! It was off because the difference density was given the same absolute
      !! screening threshold as a full build: the dropped quartets are bounded
      !! absolutely, which is no bound once the density is heading to zero. On
      !! w172/6-31G (450 atoms, 1950 AOs) the energy bounced by 1e-4 to 1e-6 Eh per
      !! iteration while RMS(D) reached 1e-10, since each difference build drops a
      !! different set and those errors do not cancel. `incr_thresh_scale` and
      !! `incr_start` are the fix: 200 unconverged iterations became 20, at the same
      !! energy as incremental-off (-11395.98363707), in 23.3 s against 24.4 s.
      integer :: incr_reset = 20
      real(dp) :: incr_thresh_scale = 1.0e-4_dp
      !! Multiplies `eri_thresh` for an incremental build only, a fixed tighter
      !! screen. Scaling by ||ddif||/||dref|| was tried first and is backwards: it
      !! loosens as ddif shrinks, losing the speed an increment exists for (25.6 s
      !! on w172). A fixed floor keeps per-quartet error small enough that unlucky
      !! non-cancelling drops cannot add up past `conv_energy`. 1e-4 reproduces the
      !! incremental-off energy to 2e-10 on w172 while costing less than a full build.
      real(dp) :: incr_start = 1.0e-2_dp
      !! Do not increment until the previous iteration's commutator norm is below
      !! this. Early on D swings by a large fraction of itself (on w172 the energy
      !! jumps 1267 Eh between iterations 3 and 4), so ddif is no smaller than D and
      !! the fixed floor makes an increment survey more quartets than a full build.
      !! Gating on `errmax < incr_start` turned a 34 s run into 23.3 s (reference
      !! 24.4 s) at the same energy.
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
      real(dp) :: level_shift = 1.0_dp
      !! Level shift on the virtual space while DIIS is off, in hartree.
      !! From a Wolfsberg-Helmholz guess cholesterol sits 20 hartree above
      !! its answer and a damped Roothaan step alone throws it 300 hartree
      !! away; raising the virtuals keeps the occupied space from collapsing.
      real(dp) :: lindep_thresh = 1.0e-8_dp
      !! Canonical orthogonalisation: an eigenvalue of S below this is linear
      !! dependence, projected out of X rather than scaled by S^-1/4 (which would
      !! divide by nearly zero). Diffuse bases (def2-TZVPD, aug-cc-pVDZ) put several
      !! S eigenvalues below 1e-8 on anything but a single atom. 0 restores the old
      !! `error stop` when S is not numerically positive definite.
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
      integer :: n_lindep = 0
      !! Functions dropped by canonical orthogonalisation: S eigenvalues
      !! below `opts%lindep_thresh`. 0 for every compact basis; see
      !! `lindep_thresh` for why a diffuse one is different.
      real(dp), allocatable :: dmat(:, :, :)   !! (nao, nao, nspin); restricted holds the total
      real(dp), allocatable :: cmo(:, :, :)    !! (nao, nao, nspin) orbitals
      real(dp), allocatable :: eps(:, :)       !! (nao, nspin) orbital energies
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
      logical :: pcm_on = .false.
      logical :: incr_on = .false.
      logical :: df = .false.
      logical :: talk = .false.
      logical :: root = .true.
      real(dp) :: eri_thresh_run = 0.0_dp
      real(dp) :: incr_start_run = 0.0_dp
      real(dp) :: incr_scale_run = 0.0_dp
      integer :: incr_reset_run = 0
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
      call scf_spin_rule(nalpha, nbeta, opts%unrestricted, opts%frac_occ, opts%nelec_frac, &
                         this%nspin, this%nocc)
      this%occ = merge(2.0_dp, 1.0_dp, this%nspin == 1)
      this%ndamp = opts%ndamp
      if (this%ndamp < 0) this%ndamp = merge(4, 0, this%nspin == 2)
      this%diis_start = opts%diis_start
      if (this%diis_start < 0.0_dp) this%diis_start = merge(0.5_dp, huge(1.0_dp), this%nspin == 2)
      this%pcm_on = pcm_model_from_name(opts%pcm_model) /= PCM_MODEL_NONE
      this%incr_on = (opts%incremental .or. incr_env()) .and. this%nspin == 1
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
      if (this%df .and. (opts%incremental .or. incr_env())) then
         msg = "trc_scf: --df with --incremental is refused -- a density-fitted build is a full "// &
               "build every iteration, and the increment buys nothing"
         ok = .false.
         return
      end if
      this%eri_thresh_run = eri_thresh_env(opts%eri_thresh)
      this%incr_start_run = incr_start_env(opts%incr_start)
      this%incr_reset_run = incr_reset_env(opts%incr_reset)
      this%incr_scale_run = incr_thresh_scale_env(opts%incr_thresh_scale)
   end function scf_settings_resolve

   !
   ! nalpha/nbeta (+ whether the case is forced unrestricted, and whether
   ! this is one of SAD's fractional-occupation atomic runs) -> nspin and
   ! the per-spin occupation count the guess and the loop both fill to.
   ! `trc_scf_run` and `trc_capi.F90`'s `trc_run_scf`/`trc_load_scf` all need
   ! this before a density exists to look at; one rule, called from both.
   !
   pure subroutine scf_spin_rule(nalpha, nbeta, unrestricted, frac_occ, nelec_frac, nspin, nocc)
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
   end subroutine scf_spin_rule

   logical function incr_env() result(res)
      !! Experiment hook for the incremental build, so its effect can be
      !! measured through a driver that does not expose the option yet.
      character(len=8) :: v
      call get_environment_variable("TRC_INCREMENTAL", v)
      res = (trim(v) == "1")
   end function incr_env

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
         if (ios == 0 .and. k > 0) res = k
      end if
   end function incr_reset_env

   real(dp) function incr_thresh_scale_env(dflt) result(res)
      !! Experiment hook, like `incr_reset_env`: TRC_INCR_THRESH_SCALE sweeps
      !! `incr_thresh_scale` without a rebuild.
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_INCR_THRESH_SCALE", v)
      res = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) res = t
      end if
   end function incr_thresh_scale_env

   real(dp) function incr_start_env(dflt) result(res)
      !! Experiment hook, like `incr_reset_env`: TRC_INCR_START sweeps
      !! `incr_start` without a rebuild.
      real(dp), intent(in) :: dflt
      character(len=16) :: v
      integer :: ios
      real(dp) :: t
      call get_environment_variable("TRC_INCR_START", v)
      res = dflt
      if (len_trim(v) > 0) then
         read (v, *, iostat=ios) t
         if (ios == 0 .and. t > 0.0_dp) res = t
      end if
   end function incr_start_env

   real(dp) function eri_thresh_env(dflt) result(res)
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
   end function eri_thresh_env

end module trc_scf_settings
