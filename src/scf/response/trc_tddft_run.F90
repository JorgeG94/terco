module trc_tddft_run
   !!
   !! TDDFT / TDHF excitation energies and the excited-state gradient of one
   !! root, from a converged closed-shell SCF: the entry the C ABI and the
   !! `terco` command line call, and the one to use from Fortran without
   !! assembling the pair list, the integral engine, the XC grid and the
   !! reference by hand.
   !!
   !!     type(trc_tddft_run_options_t) :: opts
   !!     type(trc_tddft_run_result_t) :: res
   !!     call trc_scf_run(bas, nocc, nocc, sopts, scf)
   !!     opts%nroots = 4                ! TDA by default; opts%tda = .false. for the RPA
   !!     call trc_tddft_run_solve(bas, la, scf, nocc, sopts, opts, res)
   !!     opts%target_root = 2           ! 1-based
   !!     call trc_tddft_run_gradient(bas, la, scf, nocc, sopts, opts, res)
   !!     ! res%omega, res%osc, res%ynorm, res%grad (3, natm) of root 2
   !!
   !! WHAT IT BUILDS. A pair list and a four-centre engine at the SCF's
   !! screening threshold and decontraction policy (the one trc_run_polarizability
   !! uses: the response runs on the SCF's working basis), and for Kohn-Sham the
   !! XC grid. The grid is the SCF's own: a caller that kept the SCF's
   !! `trc_grid_cache_t` passes it and nothing is rebuilt; otherwise it is
   !! regenerated from the same settings (level, rotation, batch size), and the
   !! point count is checked against the SCF's.
   !!
   !! GRID RESPONSE IS ON BY DEFAULT. The SCF's grid moves with the nuclei, so the
   !! exact derivative of terco's own excitation energy includes the motion of
   !! the points and of the partition weights. `grid_response` (default true)
   !! needs the partition that built the grid, which is why the grid is rebuilt
   !! with `dft_grid_build` (or taken from the cache) rather than from the SCF's
   !! bare points. Turn it off to get the Pulay-only gradient, which is what
   !! pyscf's tdrks computes (about 1e-5 Eh/bohr away on water PBE).
   !!
   !! WHAT IS REFUSED, in `res%err` and before any work (trc_tddft_run_refusal):
   !! triplets, UHF/UKS, meta-GGAs, range-separated hybrids, VV10, double
   !! hybrids, a spherical basis, a reference that is not (nao, nao) (functions
   !! dropped for linear dependence), a root outside 1..nroots.
   !!
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name
   use trc_dft_grid, only: dft_grid_t, dft_grid_build, dft_grid_options_t
   use trc_decontract, only: trc_decontract_options_t, decontract_options_for_exchange
   use trc_scf_settings, only: trc_scf_options_t, trc_scf_result_t, scf_settings_eri_thresh_env
   use trc_scf_context, only: trc_grid_cache_t
   use trc_tddft, only: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_reference_t, trc_tddft_solve_rks, &
                        trc_tddft_refuse_reason
   use trc_tddft_grad, only: trc_tddft_grad_options_t, trc_tddft_gradient_rks
   use trc_tddft_rpa_grad, only: trc_tddft_rpa_gradient_rks
   implicit none
   private

   public :: trc_tddft_run_options_t, trc_tddft_run_result_t
   public :: trc_tddft_run_solve, trc_tddft_run_gradient, trc_tddft_run_refusal, trc_tddft_run_gap_note
   public :: trc_tddft_run_pair_label

   real(dp), parameter :: SCREEN_NEAR_DEGENERATE = 1.0e-6_dp
   !! Roots closer than this (hartree) have no unique gradient.

   type :: trc_tddft_run_options_t
      !! What to compute. Every default is the library's recommended one.
      integer :: nroots = 3
      !! Lowest singlet roots wanted (at least 1).
      logical :: tda = .true.
      !! True: Tamm-Dancoff (CIS for Hartree-Fock). False: the full RPA.
      logical :: singlet = .true.
      !! Only singlets exist; false is refused.
      integer :: target_root = 1
      !! The root `trc_tddft_run_gradient` differentiates, 1-based, at most `nroots`.
      real(dp) :: conv_tol = 1.0e-8_dp
      !! The Davidson residual norm, tighter than trc_tddft's own 1e-6 because
      !! the gradient is first order in the amplitude error.
      integer :: max_iteration = 50
      !! Davidson iterations.
      integer :: missed_root_iterations = 24
      !! trc_tddft_options_t%missed_root_iterations: iterations after the probe in
      !! which candidate Ritz pairs beyond the buffer are still iterated. 0 iterates
      !! none (they are reported open); a negative value is refused.
      logical :: grid_response = .true.
      !! Kohn-Sham only: differentiate the grid's points and partition weights
      !! too, the exact derivative of terco's own energy. False: Pulay only.
      integer :: roots_per_sweep = 0
      !! trc_tddft_grad_options_t%roots_per_sweep; 0 takes its default. The
      !! gradient is of one root, so this changes nothing but the pass layout.
      real(dp) :: cpks_tol = 1.0e-9_dp
      !! The Z-vector's residual norm.
      integer :: cpks_max_iteration = 100
      logical :: verbose = .false.
   end type trc_tddft_run_options_t

   type :: trc_tddft_run_result_t
      !! The excitation energies of a solve, then the gradient of its target root.
      logical :: converged = .false.
      type(error_t) :: err
      !! The refusal or failure of the last call; no error when it succeeded.
      integer :: nroots = 0
      integer :: iterations = 0
      logical :: tda = .true.
      !! The approximation of the solve this result holds; the gradient needs the same one.
      integer :: target_root = 0
      !! The root `grad` belongs to (0: no gradient yet).
      real(dp) :: e_scf = 0.0_dp
      !! The ground-state SCF energy the roots sit on.
      real(dp), allocatable :: omega(:)
      !! (nroots) excitation energies, hartree.
      real(dp), allocatable :: osc(:)
      !! (nroots) oscillator strengths.
      real(dp), allocatable :: ynorm(:)
      !! (nroots) sum Y^2 (X.X - Y.Y = 1): zero for the TDA, the RPA's de-excitation weight.
      real(dp), allocatable :: tdip(:, :)
      !! (3, nroots) transition dipoles, atomic units.
      integer, allocatable :: hole(:), particle(:)
      !! (nroots) the dominant orbital pair of each root: the occupied orbital
      !! counted down from the HOMO (1 = HOMO) and the virtual counted up from
      !! the LUMO (1 = LUMO). Orbital phases and the order of a degenerate
      !! pair do not enter, so it is comparable between geometries.
      real(dp), allocatable :: weight(:)
      !! (nroots) the share of that pair, (X^2 - Y^2) of the pair over the sum.
      real(dp), allocatable :: grad(:, :)
      !! (3, natm) hartree/bohr: the full gradient of E_SCF + omega(target_root).
      integer :: cpks_iterations = 0
      type(trc_tddft_result_t) :: td
      !! The solver's own result, amplitudes included, for the gradient.
   end type trc_tddft_run_result_t

   type :: run_state_t
      !! What a solve or a gradient builds and tears down.
      type(trc_pairlist_t) :: pl
      type(trc_eri_t) :: eri
      type(trc_xc_grid_t) :: xg
      type(dft_grid_t) :: grid
      type(trc_xc_functional_t) :: func
      type(trc_tddft_reference_t) :: ref
      logical :: dft = .false.
      logical :: own_grid = .false.
      !! `xg` belongs to this state (to release); false when the cache's grid is used instead.
      logical :: from_cache = .false.
      !! The SCF's grid cache supplies the XC grid and its partition; `xg` is unused.
      logical :: built = .false.
   end type run_state_t

contains

   function trc_tddft_run_refusal(sopts, opts, b) result(msg)
      !! Why this run cannot happen, or an empty string: the solver's own
      !! reasons (functional, spin), an unrestricted reference, a spherical
      !! basis, an out-of-range root. Answered from the settings alone.
      type(trc_scf_options_t), intent(in) :: sopts
      type(trc_tddft_run_options_t), intent(in) :: opts
      type(trc_basis_t), intent(in), optional :: b
      character(len=256) :: msg
      type(trc_tddft_options_t) :: probe
      type(trc_xc_functional_t) :: func
      logical :: found

      probe%nroots = opts%nroots
      probe%singlet = opts%singlet
      probe%missed_root_iterations = opts%missed_root_iterations
      msg = ""
      if (len_trim(sopts%functional) > 0) then
         found = xc_functional_by_name(trim(sopts%functional), func)
         if (.not. found) then
            msg = "trc_tddft: unknown functional '"//trim(sopts%functional)//"'"
         else
            msg = trc_tddft_refuse_reason(probe, func)
         end if
      else
         msg = trc_tddft_refuse_reason(probe)
      end if
      if (len_trim(msg) > 0) return
      if (opts%target_root < 1 .or. opts%target_root > opts%nroots) then
         msg = "trc_tddft: the target root "//to_string(opts%target_root)//" is outside 1.."//to_string(opts%nroots)
      else if (sopts%df) then
         msg = "trc_tddft: a density-fitted (--df) reference is not implemented"
      else if (len_trim(sopts%pcm_model) > 0) then
         msg = "trc_tddft: a continuum solvent has no response here"
      end if
      if (len_trim(msg) > 0) return
      if (present(b)) then
         if (b%spherical) msg = "trc_tddft: a spherical basis is not implemented (use Cartesian functions: --cart)"
      end if
   end function trc_tddft_run_refusal

   pure function trc_tddft_run_pair_label(hole, particle) result(text)
      !! "HOMO-1 -> LUMO+2": an orbital pair counted from the frontier orbitals (1 = HOMO, 1 = LUMO).
      integer, intent(in) :: hole, particle
      character(len=:), allocatable :: text
      text = "HOMO"
      if (hole > 1) text = text//"-"//to_string(hole - 1)
      text = text//" -> LUMO"
      if (particle > 1) text = text//"+"//to_string(particle - 1)
   end function trc_tddft_run_pair_label

   pure function trc_tddft_run_gap_note(omega, root, window) result(near)
      !! The nearest root of the set to `root` closer than `window` hartree
      !! (the index, 0 if none): the optimiser's check that the followed state
      !! is not about to cross another.
      real(dp), intent(in) :: omega(:)
      integer, intent(in) :: root
      real(dp), intent(in) :: window
      integer :: near
      integer :: k
      real(dp) :: best

      near = 0
      best = window
      do k = 1, size(omega)
         if (k == root) cycle
         if (abs(omega(k) - omega(root)) < best) then
            best = abs(omega(k) - omega(root))
            near = k
         end if
      end do
   end function trc_tddft_run_gap_note

   subroutine run_state_build(b, scf, nocc, sopts, grid_response, state, err, cache)
      !! The pair list, the engine, the grid and the reference of `scf`.
      type(trc_basis_t), intent(in) :: b
      type(trc_scf_result_t), intent(in) :: scf
      integer, intent(in) :: nocc
      type(trc_scf_options_t), intent(in) :: sopts
      logical, intent(in) :: grid_response
      type(run_state_t), intent(inout) :: state
      type(error_t), intent(inout) :: err
      type(trc_grid_cache_t), intent(inout), optional :: cache
      type(trc_decontract_options_t) :: policy
      type(dft_grid_options_t) :: gopts
      real(dp) :: thresh
      logical :: found, use_cache
      integer :: npoint

      state%dft = len_trim(sopts%functional) > 0
      thresh = scf_settings_eri_thresh_env(sopts%eri_thresh)
      policy = decontract_options_for_exchange(sopts%decontract, 1.0_dp, 0.0_dp)
      if (state%dft) then
         found = xc_functional_by_name(trim(sopts%functional), state%func)
         if (.not. found) then
            call err%set(ERROR_VALIDATION, "trc_tddft: unknown functional '"//trim(sopts%functional)//"'")
            return
         end if
         policy = decontract_options_for_exchange(sopts%decontract, state%func%exx, state%func%exx_lr)
      end if
      call state%pl%build(b, thresh)
      call state%pl%to_device()
      call state%eri%build(b, thresh, decontract=policy)
      state%built = .true.
      if (.not. state%dft) then
         call state%ref%set(nocc, scf%cmo(:, :, 1), scf%eps(:, 1), scf%dmat(:, :, 1))
         return
      end if
      npoint = 0
      if (allocated(scf%grid_w)) npoint = size(scf%grid_w)
      use_cache = .false.
      if (present(cache)) then
         use_cache = cache%grid_built .and. (cache%grid%n_points == npoint .or. npoint == 0)
      end if
      if (use_cache) then
         state%grid%partition = cache%grid%partition
         state%from_cache = .true.
      else if (grid_response .or. npoint == 0) then
         ! A result read back from a file has no grid of its own: regenerate it from the settings.
         gopts%angular_rotation = sopts%grid_rotation
         call dft_grid_build(b%at_r, b%at_elem, state%grid, err, level=sopts%grid_level, options=gopts)
         if (err%has_error()) return
         if (npoint > 0 .and. state%grid%n_points /= npoint) then
            call err%set(ERROR_VALIDATION, "trc_tddft: the rebuilt grid has "//to_string(state%grid%n_points)// &
                         " points, the SCF's has "//to_string(npoint)//"; the SCF grid settings changed")
            return
         end if
         call state%xg%build_from_grid(state%grid, b, max_pts=sopts%grid_max_pts)
         call state%xg%to_device()
         state%own_grid = .true.
      else
         call state%xg%build(npoint, scf%grid_r, scf%grid_w, b, max_pts=sopts%grid_max_pts)
         call state%xg%to_device()
         state%own_grid = .true.
      end if
      call state%ref%set(nocc, scf%cmo(:, :, 1), scf%eps(:, 1), scf%dmat(:, :, 1), state%func)
   end subroutine run_state_build

   subroutine run_state_release(state)
      type(run_state_t), intent(inout) :: state
      if (.not. state%built) return
      if (state%own_grid) call state%xg%release()
      call state%eri%release()
      call state%pl%release()
      state%built = .false.
   end subroutine run_state_release

   subroutine solve_on_state(b, la, state, topts, td, cache)
      !! The solver on the state's grid, or the cache's.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(run_state_t), intent(inout) :: state
      type(trc_tddft_options_t), intent(in) :: topts
      type(trc_tddft_result_t), intent(out) :: td
      type(trc_grid_cache_t), intent(inout), optional :: cache
      if (state%from_cache) then
         call trc_tddft_solve_rks(b, state%pl, state%eri, la, state%ref, topts, td, cache%xg)
      else
         call trc_tddft_solve_rks(b, state%pl, state%eri, la, state%ref, topts, td, state%xg)
      end if
   end subroutine solve_on_state

   subroutine gradient_on_state(b, la, state, tda, one, gopts, cache)
      !! The TDA or RPA gradient of a one-root result on the state's grid, or the cache's.
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(run_state_t), intent(inout) :: state
      logical, intent(in) :: tda
      type(trc_tddft_result_t), intent(in) :: one
      type(trc_tddft_grad_options_t), intent(inout) :: gopts
      type(trc_grid_cache_t), intent(inout), optional :: cache
      if (state%from_cache) then
         if (tda) then
            call trc_tddft_gradient_rks(b, state%pl, state%eri, la, state%ref, one, gopts, cache%xg)
         else
            call trc_tddft_rpa_gradient_rks(b, state%pl, state%eri, la, state%ref, one, gopts, cache%xg)
         end if
      else
         if (tda) then
            call trc_tddft_gradient_rks(b, state%pl, state%eri, la, state%ref, one, gopts, state%xg)
         else
            call trc_tddft_rpa_gradient_rks(b, state%pl, state%eri, la, state%ref, one, gopts, state%xg)
         end if
      end if
   end subroutine gradient_on_state

   subroutine characterise(td, nocc, res)
      !! The dominant orbital pair, its share, ||Y||^2 and the oscillator
      !! strengths of every root, into `res`.
      type(trc_tddft_result_t), intent(in) :: td
      integer, intent(in) :: nocc
      type(trc_tddft_run_result_t), intent(inout) :: res
      real(dp), allocatable :: share(:, :)
      integer :: k, nr
      integer :: best(2)

      nr = size(td%omega)
      allocate (res%ynorm(nr), res%hole(nr), res%particle(nr), res%weight(nr))
      allocate (share(size(td%xamp, 1), size(td%xamp, 2)))
      do k = 1, nr
         res%ynorm(k) = sum(td%yamp(:, :, k)**2)
         share = td%xamp(:, :, k)**2 - td%yamp(:, :, k)**2
         best = maxloc(abs(share))
         res%hole(k) = nocc - best(1) + 1
         res%particle(k) = best(2)
         res%weight(k) = share(best(1), best(2))/sum(share)
      end do
   end subroutine characterise

   subroutine result_reset(res)
      !! No roots, no gradient, no error: what a failed or refused solve leaves, so a stale result cannot be read.
      type(trc_tddft_run_result_t), intent(inout) :: res
      res%converged = .false.
      res%err = error_t()
      res%nroots = 0
      res%iterations = 0
      res%target_root = 0
      if (allocated(res%grad)) deallocate (res%grad)
      if (allocated(res%omega)) deallocate (res%omega, res%osc, res%tdip)
      if (allocated(res%ynorm)) deallocate (res%ynorm, res%hole, res%particle, res%weight)
   end subroutine result_reset

   subroutine trc_tddft_run_solve(b, la, scf, nocc, sopts, opts, res, cache)
      !!
      !! The lowest `opts%nroots` singlet roots of the converged restricted SCF
      !! `scf` (orbitals in the Cartesian basis `b`, `nocc` doubly occupied,
      !! `sopts` the settings the SCF ran with): `res%omega`, `res%osc`,
      !! `res%ynorm`, `res%tdip` and the dominant pairs. A refusal or a failed
      !! solve is in `res%err`. `cache` is the SCF's own grid cache, to save
      !! rebuilding the grid; the Kohn-Sham solve needs no partition.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_scf_result_t), intent(in) :: scf
      integer, intent(in) :: nocc
      type(trc_scf_options_t), intent(in) :: sopts
      type(trc_tddft_run_options_t), intent(in) :: opts
      type(trc_tddft_run_result_t), intent(inout) :: res
      type(trc_grid_cache_t), intent(inout), optional :: cache
      type(run_state_t) :: state
      type(trc_tddft_options_t) :: topts
      character(len=256) :: reason

      call result_reset(res)
      res%tda = opts%tda
      reason = trc_tddft_run_refusal(sopts, opts, b)
      if (len_trim(reason) == 0 .and. (scf%nspin /= 1 .or. size(scf%cmo, 1) /= b%nao .or. size(scf%cmo, 2) /= b%nao)) then
         reason = "trc_tddft: needs a restricted reference whose orbitals are (nao, nao) in this basis"
      end if
      if (len_trim(reason) > 0) then
         call res%err%set(ERROR_VALIDATION, trim(reason))
         return
      end if
      call run_state_build(b, scf, nocc, sopts, .false., state, res%err, cache)
      if (res%err%has_error()) then
         call run_state_release(state)
         return
      end if
      topts%nroots = opts%nroots
      topts%tda = opts%tda
      topts%singlet = opts%singlet
      topts%conv_tol = opts%conv_tol
      topts%max_iteration = opts%max_iteration
      topts%missed_root_iterations = opts%missed_root_iterations
      topts%verbose = opts%verbose
      call solve_on_state(b, la, state, topts, res%td, cache)
      call run_state_release(state)
      if (.not. res%td%converged) then
         res%err = res%td%err
         return
      end if
      res%nroots = opts%nroots
      res%iterations = res%td%iterations
      res%e_scf = scf%energy
      res%omega = res%td%omega
      res%osc = res%td%osc
      res%tdip = res%td%tdip
      if (allocated(res%ynorm)) deallocate (res%ynorm, res%hole, res%particle, res%weight)
      call characterise(res%td, nocc, res)
      res%converged = .true.
   end subroutine trc_tddft_run_solve

   subroutine single_root(td, k, one)
      !! The solver result of root k alone (amplitudes and energy), so the
      !! gradient solves one Z-vector, not all of them.
      type(trc_tddft_result_t), intent(in) :: td
      integer, intent(in) :: k
      type(trc_tddft_result_t), intent(out) :: one
      one%converged = .true.
      one%nroots = 1
      allocate (one%omega(1), source=td%omega(k:k))
      allocate (one%xamp(size(td%xamp, 1), size(td%xamp, 2), 1), source=td%xamp(:, :, k:k))
      allocate (one%yamp(size(td%xamp, 1), size(td%xamp, 2), 1), source=td%yamp(:, :, k:k))
   end subroutine single_root

   subroutine trc_tddft_run_gradient(b, la, scf, nocc, sopts, opts, res, cache)
      !!
      !! The nuclear gradient of E_SCF + omega of root `opts%target_root` of the
      !! solve in `res`, into `res%grad` (3, natm), hartree/bohr: the TDA
      !! gradient for a TDA solve, the RPA one otherwise. Kohn-Sham runs the
      !! grid response unless `opts%grid_response` is false (module header).
      !! A refusal or failure is in `res%err`.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_scf_result_t), intent(in) :: scf
      integer, intent(in) :: nocc
      type(trc_scf_options_t), intent(in) :: sopts
      type(trc_tddft_run_options_t), intent(in) :: opts
      type(trc_tddft_run_result_t), intent(inout) :: res
      type(trc_grid_cache_t), intent(inout), optional :: cache
      type(run_state_t) :: state
      type(trc_tddft_run_options_t) :: solved
      type(trc_tddft_result_t) :: one
      type(trc_tddft_grad_options_t) :: gopts
      character(len=256) :: reason

      res%err = error_t()
      if (allocated(res%grad)) deallocate (res%grad)
      res%target_root = 0
      if (.not. res%converged) then
         reason = "trc_tddft: solve the excitation energies first"
      else
         ! The root is checked against the roots the solve produced, not against the options' count.
         solved = opts
         solved%nroots = res%nroots
         reason = trc_tddft_run_refusal(sopts, solved, b)
         if (len_trim(reason) == 0 .and. opts%tda .neqv. res%tda) then
            reason = "trc_tddft: the gradient's approximation ("//merge("TDA", "RPA", opts%tda)// &
                     ") is not that of the solve ("//merge("TDA", "RPA", res%tda)//")"
         end if
      end if
      if (len_trim(reason) > 0) then
         call res%err%set(ERROR_VALIDATION, trim(reason))
         return
      end if
      call warn_degenerate(res%omega, opts%target_root)
      call run_state_build(b, scf, nocc, sopts, opts%grid_response, state, res%err, cache)
      if (res%err%has_error()) then
         call run_state_release(state)
         return
      end if
      call single_root(res%td, opts%target_root, one)
      gopts%roots_per_sweep = opts%roots_per_sweep
      gopts%cpks_tol = opts%cpks_tol
      gopts%cpks_max_iteration = opts%cpks_max_iteration
      gopts%verbose = opts%verbose
      gopts%grid_response = opts%grid_response .and. state%dft
      if (gopts%grid_response) gopts%partition = state%grid%partition
      call gradient_on_state(b, la, state, opts%tda, one, gopts, cache)
      call run_state_release(state)
      if (gopts%err%has_error()) then
         res%err = gopts%err
         return
      end if
      res%grad = gopts%grad(:, :, 1)
      res%target_root = opts%target_root
      res%cpks_iterations = maxval(gopts%cpks_iterations)
   end subroutine trc_tddft_run_gradient

   subroutine warn_degenerate(omega, root)
      !! A root inside a degenerate multiplet has no unique gradient: say so.
      real(dp), intent(in) :: omega(:)
      integer, intent(in) :: root
      integer :: near
      near = trc_tddft_run_gap_note(omega, root, SCREEN_NEAR_DEGENERATE)
      if (near > 0) then
         call global_logger%warning("tddft: root "//to_string(root)//" is degenerate with root "//to_string(near)// &
                                    " (|d omega| < 1e-6): its gradient is not unique")
      end if
   end subroutine warn_degenerate

end module trc_tddft_run
