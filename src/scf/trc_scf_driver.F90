!
!
! The SCF: a basis and an electron count in, a converged energy and
! density out, with every matrix resident on the device. Hartree-Fock or
! Kohn-Sham (any functional in trc_xc_functional), restricted or unrestricted.
!
! WHAT CROSSES THE BUS
! --------------------
! Per iteration: two energies, a norm, and, once DIIS runs, three packed
! triangles (the new Fock matrix and commutator down, the extrapolated Fock
! matrix up): the DIIS history lives on the host, where it costs no device
! memory (trc_scf_diis's own header). The density, Fock matrices, orbitals
! and orthogonaliser stay on the device; the Fock build is the resident
! driver, the XC integrator finds its density already present, and the
! products, commutator and eigenproblem are trc_linalg (cuBLAS/cuSOLVER
! under OpenACC). This follows
! metalquicha's device SCF, which measured host LAPACK plus copies at three
! times the cost of the GPU Fock build. A host build runs the same source with
! LAPACK, and the numbers agree to rounding, so the host suite tests this module.
!
! THE ITERATION
! -------------
! A few damped steps from the guess, then Pulay's DIIS on the commutator
! FDS - SDF in the orthogonal basis, per spin with one set of coefficients.
! Converged when the energy change and the RMS density change are both below
! threshold. The DIIS system (at most 11 by 11) is solved by pivoted
! elimination as in metalquicha's mqc_diis, with the error block scaled to
! order one against the constraint border first: near convergence the overlaps
! are 1e-16 and an unscaled elimination makes the weights noise, showing as a
! density that drifts from the fixed point while the energy sits flat.
!
! UNRESTRICTED
! ------------
! F_s = H + J[D_a + D_b] - a_x K[D_s] + V_xc,s. The Fock driver returns
! j J[D] - k K[D]/2 for a density it takes to be a total, so J is one call on
! D_a + D_b with k = 0 and -a_x K[D_s] one call on 2 D_s with j = 0, k = a_x:
! three builds per iteration against one restricted. A batched build over the
! three densities is the obvious improvement and is not done.
!
! SPHERICAL BASES
! ---------------
! `b%spherical` (trc_api.F90, from a JSON basis's function_type) changes the
! size of every matrix this module owns (S, H, D, F, C, eps, X, the DIIS
! history) from `b%nao` to `b%nao_sph`; everything operating on them is
! unchanged, just at a smaller `nao`. Everything touching a kernel stays
! Cartesian: `trc_1e`, `trc_eri_t%fock_resident`, and the XC integrator.
! Cartesian-sized scratch matrices (suffix `_c`: `s_c`, `h_c`, `d_c`,
! `gmat_c`, `vxc_c`, ...) sit between the kernels and the SCF's own matrices,
! built once per geometry (`s_c`/`h_c`, now `trc_scf_context`'s `smat`/`hcore`)
! or once per iteration (`d_c` from the current density; `gmat_c`/`vxc_c` from
! the Fock/XC calls, now `trc_scf_state`'s), and folded by `sph_cart_to_sph`.
! `sph` guards every call into trc_sph.F90, so the Cartesian path is the same
! code path as before spherical support existed.
!
! The XC integrator gets `d_c`, the Cartesian density rebuilt from `res%dmat`
! each iteration, not spherical AO values on the grid: the transform is shared
! with the Coulomb/exchange build, and `collocation_shell` evaluates Cartesian
! components straight off a flat table, so a second AO space there would
! duplicate trc_sph.F90's block structure.
!
! SAD and GWH reach this module already in the right space (SAD's atomic
! sub-SCFs run on a `subset` basis carrying `spherical`; GWH is built here from
! `hcore`). SAP stays Cartesian and is transformed once at the trc_capi.F90
! call site, as it has no basis handle to carry `spherical`.
!
! POLARISABLE CONTINUUM
! ---------------------
! With `pcm_model` set, trc_pcm builds the cavity, K and the nuclear potential
! once. Every iteration, at the density the energy is evaluated at, the surface
! potential of the total Cartesian density gives the charges, q.v/2 into the
! energy and the charges' full operator into each spin's Fock matrix (as
! pyscf.solvent.pcm's SCFWithSolvent, which also adds it before DIIS). It is
! rebuilt in full each iteration, independent of the incremental-ERI path.
!
! MODULES
! -------
! `trc_scf_settings` resolves `trc_scf_options_t` into what the loop reads
! (`scf_settings_t`, "scf_config" below); `trc_scf_context` builds everything once
! per geometry (`scf_ctx_t`, "scf_context": the pair list, S/H, the ERI engine, the
! linear-algebra handles, the orthogonaliser, and -- when asked for -- the
! PCM cavity and the DFT grid); `trc_scf_state` owns the per-iteration arrays
! (`scf_state_t`, "scf_state") and the one enter/exit pair that puts them, `scf_context`'s
! own matrices and `scf_result` on the device together. `scf_fock_build_g`/`scf_fock_assemble`
! and the DIIS push/combine helpers below stay in THIS file for now (they
! are not yet split into `trc_scf_fock`/`trc_scf_diis`): every one of them
! already takes its arrays as explicit-shape dummies rather than touching
! `scf_context`/`scf_state` inside a `do concurrent` body, which is what lets this compile
! correctly once those arrays are struct components instead of plain locals.
!
module trc_scf_driver
   use pic_types, only: dp, int64, int_index
   use trc_basis_types, only: trc_basis_t
   use trc_ecp_types, only: trc_ecp_t
   use trc_error, only: error_t
   use trc_eri, only: trc_eri_t
   use trc_linalg, only: trc_linalg_t
   use trc_scf_settings, only: trc_scf_options_t, trc_scf_guess_t, trc_scf_result_t, trc_scf_timings_t, scf_settings_t, &
                               scf_energies_t, &
                               scf_settings_homo_lumo_gap
   use trc_scf_context, only: scf_ctx_t, scf_context_extras_t, scf_context_bcast_flat, trc_grid_cache_t, &
                              trc_pair_eri_cache_t, scf_context_project_density
   use trc_scf_state, only: scf_state_t, scf_timers_t
   use trc_scf_fock, only: scf_fock_build_g, scf_fock_build_vxc, scf_fock_build_pcm, scf_fock_energy, scf_fock_assemble, &
                           scf_fock_commutator_error, scf_fock_to_cart, scf_fock_fold_g
   use trc_scf_df, only: df_jk_t
   use trc_scf_diis, only: cdiis_t, adiis_t, scf_diis_level_shift, scf_diis_damp
   use trc_scf_soscf, only: trc_scf_soscf_run
   use trc_scf_incremental, only: scf_incremental_t
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_atom_orbits
   use trc_sym_run, only: sym_run_skeleton_mode, sym_run_skeleton_size_rule, SYM_SKEL_AUTO, SYM_SKEL_OFF
   use pic_types, only: default_int
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_sorting, only: sort_index
   use pic_mpi_lib, only: comm_t
   implicit none
   private

   public :: trc_scf_options_t, trc_scf_guess_t, trc_scf_result_t, trc_scf_run, trc_grid_cache_t, trc_pair_eri_cache_t
   ! check_scf_stabilise unit-tests the MOM reordering on a synthetic problem.
   public :: scf_driver_mom_weights, scf_driver_mom_permute
   ! ... and the stall guard of the shift that level_shift_until keeps on.
   public :: scf_shift_guard_t, scf_driver_shift_stalled, scf_driver_shift_now, SHIFT_STALL_WINDOW, SHIFT_STALL_FACTOR
   public :: scf_driver_mom_active, MOM_START

   real(dp), parameter :: MOM_START = 0.1_dp
   !! Commutator norm below which MOM starts following the density of a run
   !! that began from a guess it built itself. Below this the occupied space
   !! is settled enough to follow; the stalls MOM is for sit at 1e-2..1e-1,
   !! while a poor guess has a far larger commutator and MOM would lock it in.

   integer, parameter :: SHIFT_STALL_WINDOW = 10
   !! Iterations, shifted only because of `level_shift_until`, over which the
   !! stall guard asks for progress. The plain O3+ UHF/6-31G* run (36
   !! iterations from the SAD guess, 35 from SADQ; including a plateau of seven
   !! iterations at 4e-3) never has a window of ten whose best commutator is
   !! above 0.28 of the best before it,
   !! while the shifted run with `level_shift_until` = 1e-3 sits at 3e-3 to
   !! 4e-3 for twenty iterations. Shorter windows (5, 6, 8) saved two or three
   !! iterations there and fire on slower but steady convergence (a window of
   !! 5 fires at a rate above 0.87 per iteration, one of 10 only above 0.93).
   real(dp), parameter :: SHIFT_STALL_FACTOR = 0.5_dp
   !! The best commutator of the last `SHIFT_STALL_WINDOW` iterations must be
   !! below this fraction of the best before them, or the guard fires. On
   !! O3+ the factor hardly matters (0.3, 0.5 and 0.7 end in 48, 48 and 49
   !! iterations): the shifted run is either parked or falling fast.

   type :: scf_shift_guard_t
      !! The memory of the stall guard of `level_shift_until`: the last
      !! `SHIFT_STALL_WINDOW` commutator norms seen on iterations shifted
      !! only because of that option, and the best of all the earlier ones.
      !! It reads nothing but the commutator (broadcast, so identical on every
      !! rank) and a count, so every rank drops the shift on the same
      !! iteration. Once `dropped` it stays so for the run.
      integer :: n_seen = 0
      real(dp) :: recent(SHIFT_STALL_WINDOW) = huge(1.0_dp)   !! a ring buffer, the oldest overwritten first
      real(dp) :: older_best = huge(1.0_dp)   !! the smallest commutator that has left the window
      logical :: dropped = .false.
   contains
      procedure :: observe => scf_driver_shift_observe
   end type scf_shift_guard_t

   integer, parameter :: ORB_AUFBAU = 0, ORB_FRAC = 1, ORB_MOM = 2
   !! How a diagonalisation's eigenvectors are occupied: lowest first, the
   !! fractional occupations of the free-atom runs, or maximum overlap with
   !! the previous density.

   type :: scf_skeleton_t
      !! The verdict on the skeleton Fock build of a symmetric run, and what
      !! the run needs of it: the decision, the reason it came out so (for the
      !! run's own log), the orbit count the size rule saw and the group
      !! geometry the engine is built with.
      logical :: on = .false.   !! the engine holds the group's orbits and G is a skeleton sum
      integer :: norbit = 0     !! atoms inequivalent under the group, 0 when the size rule was not consulted
      character(len=48) :: why = "asked for"   !! why `on` came out as it did
      type(trc_sym_geom_t) :: geom   !! the group and atom permutations the engine is built with
   end type scf_skeleton_t

contains

   subroutine dump_density(nao, d)
      !! The density at a chosen iteration, written where a benchmark can read
      !! it back. Off unless TRC_DUMP_DENSITY names a file. Stream format:
      !! an 8-byte nao, then the matrix, which is what the test readers expect.
      !! The SCF writes the density leaving iteration 2. The model density
      !! bench_gc invents is small and banded, so it screens far harder than
      !! anything an SCF actually meets, and a Fock build timed on it flatters
      !! the code by about three times. The density lives on the DEVICE (every
      !! update since the guess was a gemm up there), so the caller pulls it
      !! back first; dumping without that wrote the block-diagonal atomic SAD
      !! guess and called it iteration 2.
      integer, intent(in) :: nao
      real(dp), intent(in) :: d(nao, nao)
      character(len=512) :: path
      integer :: u, ios
      integer(int64) :: prev
      call get_environment_variable("TRC_DUMP_DENSITY", path)
      if (len_trim(path) == 0) return
      ! The atomic SCFs of a SAD guess come through here too, and would
      ! overwrite the molecular density with a 5x5 one. Keep the largest.
      open (newunit=u, file=trim(path), access="stream", form="unformatted", &
            status="old", action="read", iostat=ios)
      if (ios == 0) then
         read (u, iostat=ios) prev
         close (u)
         if (ios == 0) then
            if (prev >= int(nao, int64)) return
         end if
      end if
      open (newunit=u, file=trim(path), access="stream", form="unformatted", &
            status="replace", action="write", iostat=ios)
      if (ios /= 0) return
      write (u) int(nao, int64)
      write (u) d
      close (u)
      call global_logger%info("density dumped to "//trim(path)//", nao "//to_string(nao))
   end subroutine dump_density

   !
   ! Run the SCF. `nalpha` and `nbeta` are the electron counts; the case is
   ! restricted when they are equal and `unrestricted` is not forced. The
   ! basis must already be on the device. `dguess` is (nao, nao, nspin) in
   ! the result's convention; without it the guess is the core Hamiltonian.
   !
   recursive subroutine trc_scf_run(scf_basis, nalpha, nbeta, scf_options, scf_result, dguess, comm, guess, la, &
                                    scf_dfjk, grid_cache, pair_eri_cache, sym, ecp)
      type(trc_basis_t), intent(in) :: scf_basis
      integer, intent(in) :: nalpha, nbeta
      type(trc_scf_options_t), intent(in) :: scf_options
      type(trc_scf_result_t), intent(out) :: scf_result
      real(dp), intent(in), optional :: dguess(:, :, :)
      type(trc_scf_guess_t), intent(in), optional :: guess
      !! The guess HAMILTONIAN (SAP) or a PREVIOUS geometry's orbitals to
      !! project, when `dguess` is not given: see the type. The orbitals are
      !! re-orthonormalised at THIS geometry's overlap once
      !! `scf_context%setup` below has built it.
      type(comm_t), intent(in), optional :: comm
      !! Ranks. Every rank runs the same iteration on the same matrices;
      !! the quartets are split across them and the Fock matrix summed
      !! back, so the density stays identical everywhere. Bind devices
      !! first (trc_bind_device). The XC batches are split across ranks too.
      type(trc_linalg_t), intent(inout), optional :: la
      !! cuBLAS/cuSOLVER handles made once by the caller and carried through
      !! the calculation (trc_linalg's init is create-once and grows the
      !! workspace to nao). Absent, the run makes its own and releases them.
      type(df_jk_t), intent(inout), optional :: scf_dfjk
      !! The density-fitted J/K object for `scf_options%df` (`--df`), already fitted
      !! by the caller -- required when `scf_options%df` is true (trc_scf_context's
      !! `setup` refuses otherwise), ignored when it is false. Owned by the
      !! caller across more than one run (trc_capi.F90's context handle), the
      !! same reason `la` is a caller-owned, optional argument rather than
      !! something this routine builds for itself.
      type(trc_grid_cache_t), intent(inout), optional :: grid_cache
      !! A caller-owned DFT grid cache (`scf_context%setup` builds into it).
      !! Present, the converged grid stays alive and device-resident after
      !! this call, for a gradient right after (trc_capi.F90's
      !! `trc_run_gradient`); absent, the run makes its own and releases it.
      type(trc_pair_eri_cache_t), intent(inout), optional :: pair_eri_cache
      !! A caller-owned pair-list/ERI cache, the same way: present, the pair
      !! list and (outside `--df`) the four-centre engine stay alive for a
      !! gradient at the SAME geometry and threshold; absent, the run's own.
      type(trc_sym_scf_t), intent(inout), optional :: sym
      !! SYMMETRY (see the header of trc_sym_scf): present, with a group of more than
      !! one operation in `sym%ao`, the Fock matrix is diagonalised per irrep
      !! block, orbitals carry irrep labels (`scf_result%mo_irrep`), and D and
      !! F are symmetrised each iteration. `scf_basis` must sit on the
      !! standard frame `sym%ao` was built for. Absent, or a one-operation
      !! group, the run is the C1 code path unchanged. The Fock build is C1
      !! either way. Refused (`scf_result%message`): --soscf, whose Newton
      !! rotations ignore the blocks, and fractional occupations.
      type(trc_ecp_t), intent(in), optional :: ecp   !! passed through to `scf_context%setup`

      integer :: nao, nao_c, s, it, n2
      logical :: ok, sph
      type(scf_settings_t) :: scf_config
      type(scf_ctx_t) :: scf_context
      type(scf_context_extras_t) :: extras   !! the point group (when the skeleton build is on) and the ECP, for `setup`
      type(scf_state_t) :: scf_state
      type(scf_timers_t) :: scf_timers
      !! This iteration's share of the two totals (tf_it/t_fock, tx_it/t_xc,
      !! tr_it/t_rest), so the table says where the time went while it is
      !! going rather than only at the end. A Fock build is not the same
      !! cost twice: the SAD guess is atomic, so density screening throws
      !! nearly all of iteration 1 away, and the per-iteration cost then
      !! climbs to a plateau. Comparing two runs by total/iterations is
      !! therefore only valid at equal iteration counts.
      character(len=256) :: msg
      character(len=160) :: line
      logical :: incr_now
      type(scf_incremental_t) :: incremental   !! which builds are increments, trc_scf_incremental
      type(cdiis_t) :: cdiis
      type(adiis_t) :: adiis
      !! ADIIS's history (fock, density) is kept apart from CDIIS's own (fock,
      !! commutator): ADIIS runs exactly in the window CDIIS has not opened
      !! yet, so the two never need the same entries, and keeping them
      !! separate means `scf_options%adiis = .false.` leaves the CDIIS path
      !! untouched, to the bit.
      logical :: adiis_active
      logical :: shifted, mom_active
      !! F carried the level shift this iteration; MOM follows the density yet (`scf_driver_mom_active`).
      type(scf_shift_guard_t) :: shift_guard
      integer :: orb_mode, nswap_s
      integer :: mom_swaps(2)   !! orbitals MOM occupied above aufbau, per spin
      type(error_t) :: err
      type(scf_energies_t) :: scf_energies
      !! Every energy and integrated electron count the loop produces; all zero
      !! until the first iteration, so a run with max_iteration=0 still reports zeros.
      logical :: diis_on
      logical :: sym_on          !! the SO-blocked path (a symmetric group was handed in)
      logical :: have_hamiltonian, have_orbitals   !! which of `guess`'s parts were given
      type(scf_skeleton_t) :: skeleton   !! the skeleton Fock build's verdict (`on` only under `sym_on`)
      real(dp) :: ctl(2)
      real(dp) :: eold, drms, errmax
      integer :: soscf_switch_iteration   !! DIIS iteration --soscf took over after, 0 if never
      logical :: diag_failed      !! this iteration's Fock diagonalisation (every rank's own copy)
      integer :: diag_info        !! the FIRST failing spin's own eigensolver info, kept across the spin loop
      integer :: info_s           !! this spin's own raw info from the call just made
      real(dp) :: diag_flag(4)    !! [0.0/1.0 failed, info, mom_swaps(1:2)], broadcast from root the same way `ctl` is

      ! The handles and caches the run reads: the caller's when given, else
      ! the run's own, made here and released when the run (the recursive
      ! call) returns. Below these three blocks all three are present.
      if (.not. present(la)) then
         block
            type(trc_linalg_t) :: la_run
            call trc_scf_run(scf_basis, nalpha, nbeta, scf_options, scf_result, dguess, comm, guess, la_run, &
                             scf_dfjk, grid_cache, pair_eri_cache, sym, ecp)
            call la_run%release()
         end block
         return
      end if
      if (.not. present(grid_cache)) then
         block
            type(trc_grid_cache_t) :: grid_run
            call trc_scf_run(scf_basis, nalpha, nbeta, scf_options, scf_result, dguess, comm, guess, la, &
                             scf_dfjk, grid_run, pair_eri_cache, sym, ecp)
            call grid_run%release()
         end block
         return
      end if
      if (.not. present(pair_eri_cache)) then
         block
            type(trc_pair_eri_cache_t) :: pair_eri_run
            call trc_scf_run(scf_basis, nalpha, nbeta, scf_options, scf_result, dguess, comm, guess, la, &
                             scf_dfjk, grid_cache, pair_eri_run, sym, ecp)
            call pair_eri_run%release()
         end block
         return
      end if

      ! Everything `scf_options` resolves to before the loop starts: the functional
      ! lookup and its rejections, the spin rule, the by-spin damping/DIIS
      ! defaults, and the env-hook overrides (trc_scf_settings). `scf_config` is
      ! immutable from here on, and the rest of this routine reads it
      ! directly.
      call scf_timers%run%start()
      ok = scf_config%resolve(scf_options, nalpha, nbeta, msg=msg, comm=comm)
      if (.not. ok) then
         scf_result%message = msg
         return
      end if
      scf_result%nspin = scf_config%nspin
      sym_on = .false.
      if (present(sym)) sym_on = sym%is_active()
      if (sym_on) then
         ! Second-order SCF rotates orbitals with its own Newton step and
         ! rebuilds the density from the global orthogonaliser, so it would
         ! neither keep the irrep labels nor the symmetrisation: refused
         ! rather than run without them. The fractional occupations of the
         ! free-atom runs have no meaning for the irreps of a molecule.
         if (scf_config%soscf) then
            scf_result%message = "trc_scf: --soscf is not supported with symmetry (its Newton rotations "// &
                                 "ignore the irrep blocks); run with symmetry off"
            return
         end if
         if (scf_options%frac_occ) then
            scf_result%message = "trc_scf: fractional occupations are not supported with symmetry"
            return
         end if
         ! MOM picks the occupied orbitals by overlap across irreps, which the
         ! per-irrep occupations of the symmetric run fix instead.
         if (scf_options%mom) then
            scf_result%message = "trc_scf: mom is not supported with symmetry; run with symmetry off"
            return
         end if
      end if
      ! --- everything once per geometry: pairs, 1e, S/H fold, ERI, linalg,
      ! orthogonaliser, and (when asked for) the PCM cavity and DFT grid ---
      ! The skeleton Fock build is decided first: the
      ! engine is built with the point group or not. The cavity, K and its
      ! LU are built and solved whole on every rank (replicated: cheap next
      ! to the integrals); pcm%operator splits the two integral passes over
      ! the ranks itself, given `comm`.
      call scf_driver_skeleton_decide(scf_basis, sym, sym_on, scf_config, comm, skeleton, err)
      if (err%has_error()) then
         scf_result%message = "trc_scf: "//err%get_message()
         return
      end if
      scf_result%sym_skeleton = skeleton%on
      scf_result%sym_atom_orbits = skeleton%norbit
      if (skeleton%on) then
         extras%has_sym_geom = .true.
         extras%sym_geom = skeleton%geom
      end if
      if (present(ecp)) extras%ecp = ecp
      call scf_context%setup(scf_basis, scf_options, scf_config, la, grid_cache, pair_eri_cache, err, comm=comm, &
                             dfjk=scf_dfjk, extras=extras)
      if (err%has_error()) then
         scf_result%message = "trc_scf: "//err%get_message()
         return
      end if
      nao = scf_context%nao
      nao_c = scf_context%nao_c
      sph = scf_context%sph
      n2 = nao*nao
      scf_state%skeleton = skeleton%on
      scf_result%e_nuc = scf_context%e_nuc
      scf_result%n_lindep = scf_context%n_lindep
      if (sym_on) then
         call scf_driver_sym_prepare(sym, scf_context, scf_config, scf_options, scf_result, ok)
         if (.not. ok) return
      end if
      if (scf_config%pcm_on) scf_result%pcm_npoints = scf_context%pcm%surface%n_points
      if (scf_config%dft) then
         scf_result%grid_r = grid_cache%grid%coords
         scf_result%grid_w = grid_cache%grid%weights
         if (scf_config%nlc) then
            scf_result%grid_r_nlc = grid_cache%grid_nlc%coords
            scf_result%grid_w_nlc = grid_cache%grid_nlc%weights
         end if
         ! The DFT pipeline's own step table (grid build, semilocal XC,
         ! VV10), reported alongside the existing fock(s)/xc(s) columns
         ! below rather than instead of them: `report`'s per-step rate
         ! columns are new here (grid build has no flop count to show
         ! them against), but the three totals are not new numbers, only
         ! a second, uniform place they land -- see trc_scf_state's own
         ! header on `xcst`. "dft", not "scf": its total is grid + XC only.
         call scf_timers%xcst%init("dft", [character(len=8) :: "grid", "xc", "vv10"])
         call scf_timers%xcst%add_time(1, scf_context%t_grid)
      end if

      ! ONLY WHAT THIS RUN USES IS FULL SIZE (trc_scf_state%alloc owns the
      ! array sizes, 1 x 1 for an array the run never touches). One that a
      ! run never touches is nao^2 of device memory for nothing: on ionic05
      ! in aug-cc-pVTZ (4255 functions, 4950 Cartesian) an HF run carried
      ! ~2.4 GB of XC, unrestricted, range-separated, incremental and PCM
      ! matrices it had no use for.
      call scf_state%alloc(scf_context, scf_config)
      call cdiis%alloc(nao, scf_config%nspin, scf_options%ndiis)
      allocate (scf_result%dmat(nao, nao, scf_config%nspin), scf_result%cmo(nao, nao, scf_config%nspin), &
             scf_result%eps(nao, scf_config%nspin))
      ! Only when asked for: a second ndiis-deep history the size of
      ! cdiis's own. Both live on the host, packed (trc_scf_diis's own
      ! header), but at gly100's 16,258 functions one is still 10 GB.
      if (scf_options%adiis) call adiis%alloc(nao, scf_config%nspin, scf_options%ndiis)

      ! --- the guess, on the host, then everything goes up -------------------
      diag_failed = .false.
      diag_info = 0
      have_hamiltonian = .false.
      have_orbitals = .false.
      if (present(guess)) then
         have_hamiltonian = allocated(guess%hamiltonian)
         have_orbitals = allocated(guess%orbitals)
      end if
      if (present(dguess)) then
         scf_result%dmat = dguess
      else if (have_hamiltonian) then
         do s = 1, scf_config%nspin
            call scf_driver_guess_spin(la, scf_context, scf_config%occ, guess%hamiltonian, scf_config%nocc(s), &
                                       scf_result%dmat(:, :, s), info_s)
            call scf_driver_note_failure(info_s, diag_failed, diag_info)
         end do
      else if (have_orbitals) then
         call scf_context_project_density(la, nao, scf_config%nspin, scf_config%nocc, guess%orbitals, &
                                          scf_context%smat, scf_result%dmat)
      else if (scf_options%guess == "core") then
         do s = 1, scf_config%nspin
            call scf_driver_guess_spin(la, scf_context, scf_config%occ, scf_context%hcore, scf_config%nocc(s), &
                                       scf_result%dmat(:, :, s), info_s)
            call scf_driver_note_failure(info_s, diag_failed, diag_info)
         end do
      else
         ! GWH: F_ij = 1.75/2 S_ij (H_ii + H_jj), diagonal H_ii.
         call scf_driver_gwh_matrix(nao, scf_context%hcore, scf_context%smat, scf_state%w1)
         do s = 1, scf_config%nspin
            call scf_driver_guess_spin(la, scf_context, scf_config%occ, scf_state%w1, scf_config%nocc(s), &
                                       scf_result%dmat(:, :, s), info_s)
            call scf_driver_note_failure(info_s, diag_failed, diag_info)
         end do
      end if
      ! A diagonalisation failure in the INITIAL guess (core/GWH/a supplied
      ! Hamiltonian -- SAD/SADQ build their own density without this
      ! routine, through `dguess`): report it the same way `scf_ctx_setup`
      ! itself does on its own refusals (iterations is still 0 here, so
      ! trc_capi.F90's trc_run_scf reports TRC_ERR_UNSUPPORTED, not
      ! TRC_ERR_NOCONV -- this never started, rather than failed to finish).
      if (diag_failed) then
         scf_result%message = "trc_scf: the initial guess's own diagonalisation did not converge "// &
                              "(eigensolver info="//to_string(diag_info)//")"
         if (sym_on) call sym%release()
         call scf_context%release()
         return
      end if
      scf_result%cmo = 0.0_dp
      scf_result%eps = 0.0_dp

      call scf_state%to_device(scf_context, scf_result)
      if (sym_on) then
         ! The guess density is a C1 diagonalisation (or SAD, or a restart):
         ! average it over the group so the first Fock build sees a symmetric D.
         call scf_driver_symmetrize(sym, scf_result%dmat, scf_state%w1)
         if (scf_config%talk) call scf_driver_log_symmetry(sym, skeleton)
      end if

      eold = 0.0_dp
      ! Read by the first iteration's `incr_now` before any commutator exists.
      errmax = huge(1.0_dp)
      call incremental%setup(scf_config)
      diis_on = .false.
      mom_swaps = 0
      shifted = .false.
      mom_active = .false.
      soscf_switch_iteration = 0
      ! What the SCF is actually running on. Two callers can hand over the
      ! same molecule and basis name and different shells -- one had a
      ! reader that kept dead primitives -- and that only shows here.
      if (scf_config%talk) call global_logger%info("   basis: "//to_string(nao)//" functions, " &
         //to_string(scf_basis%nshell)//" shells, "//to_string(sum(scf_basis%sh_np))//" primitives, max " &
         //to_string(maxval(scf_basis%sh_np)))
      if (scf_config%talk) call log_incremental_mode(scf_config)
      if (scf_config%talk) then
         if (scf_config%dft) then
            call global_logger%info("   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s)     xc(s) linalg(s)")
         else
            call global_logger%info("   it        E(total)            dE          RMS(D)      |FDS-SDF|" &
               //"   fock(s) linalg(s)")
         end if
      end if
      do it = 1, scf_options%max_iteration
         call scf_timers%fock%start()
         incr_now = .false.
         ! --- D_cart := T D_sph T^T, the one transform every Cartesian
         ! quantity below (the Fock build, the XC density) is built from.
         ! Not taken at all when `sph` is false -- see this module's own
         ! header on why that keeps the Cartesian path byte-for-byte
         ! unchanged rather than merely numerically equivalent.
         do s = 1, scf_config%nspin
            call scf_fock_to_cart(scf_basis, sph, nao, nao_c, scf_result%dmat(:, :, s), scf_state%d_c(:, :, s))
         end do
         ! --- two-electron part, resident, entirely in Cartesian space -------
         !
         ! INCREMENTAL. G(D) = G(D_ref) + G(D - D_ref), exactly, because G is
         ! linear in D. The point is the screening: the kernel screens on the
         ! density it is handed, and the difference is small and getting
         ! smaller while the density is neither. See `scf_fock_build_g`'s own header.
         !
         incr_now = incremental%next_build(it, errmax, it == scf_options%max_iteration)
         call scf_fock_build_g(scf_basis, scf_context, scf_config, scf_state, incremental%scale, pair_eri_cache%eri, la, &
                               scf_dfjk)
         ! --- G_sph := T^T G_cart T, back up to the space the rest of the
         ! iteration runs in, into `fock` until the assembly below.
         call scf_fock_fold_g(scf_basis, scf_context, scf_config, scf_state)
         ! The skeleton build returns G summed over the orbit representatives
         ! only. It is not symmetrized here: the energy below is Tr(D G) with
         ! D invariant, which the average over the group leaves unchanged, and
         ! the symmetrization of F after the assembly (below, whenever `sym`
         ! is active) makes the matrix itself invariant before DIIS sees it.
         scf_timers%tf_it = scf_timers%fock%get_elapsed_time()
         scf_timers%t_fock = scf_timers%t_fock + scf_timers%tf_it
         scf_timers%n_fock = scf_timers%n_fock + 1
         if (scf_timers%n_fock == 1) scf_timers%tf_first = scf_timers%tf_it
         call scf_timers%xc%start()
         ! --- exchange-correlation: finds its (Cartesian) density present ----
         if (scf_config%dft) then
            call scf_fock_build_vxc(scf_basis, grid_cache, la, scf_config, scf_state, scf_energies, scf_timers, comm)
         end if
         scf_timers%tx_it = scf_timers%xc%get_elapsed_time()
         scf_timers%t_xc = scf_timers%t_xc + scf_timers%tx_it
         call scf_timers%rest%start()
         ! --- the continuum, at the TOTAL density (pyscf's choice for UKS
         ! too: one set of charges, one operator added to both spins) -----
         if (scf_config%pcm_on) then
            call scf_timers%pcm%start()
            call scf_fock_build_pcm(scf_basis, scf_context, scf_config, scf_state, scf_energies, la, &
                                    pair_eri_cache%pl, comm)
            scf_timers%t_pcm = scf_timers%t_pcm + scf_timers%pcm%get_elapsed_time()
            ! the PCM step is its own column, not part of the rest
            call scf_timers%rest%start()
         end if
         ! --- energy at this density -----------------------------------------
         call scf_fock_energy(scf_context, scf_config, scf_state, scf_result, scf_energies, la)
         ! --- Fock and its commutator, per spin ------------------------------
         call scf_fock_assemble(scf_context, scf_config, scf_state)
         ! F is symmetric up to the round-off of the build and, for DFT, the
         ! grid noise of the quadrature; averaging over the group removes both
         ! before DIIS and the commutator see it.
         if (sym_on) call scf_driver_symmetrize(sym, scf_state%fock, scf_state%w1)
         call scf_fock_commutator_error(scf_context, scf_config, scf_state, scf_result, la)
         ! `dold` must hold `scf_result%dmat` from BEFORE DIIS/diagonalisation touch
         ! it -- see trc_scf_state's own header.
         call scf_state%save_dold(nao, scf_config%nspin, scf_result%dmat)
         ! --- DIIS, once the commutator is small enough to extrapolate from --
         errmax = sqrt(la%dot_3d(scf_state%errv, scf_state%errv))
         ! Over ranks the Fock matrices agree only to the XC atomics' 1e-13,
         ! so the numbers every branch below and the convergence test read
         ! are rank 0's, on every rank: one rank leaving the loop while
         ! another iterates on would hang the next collective.
         if (present(comm)) then
            if (comm%size() > 1) then
               ctl = [scf_energies%total, errmax]
               call scf_context_bcast_flat(comm, ctl, 2)
               scf_energies%total = ctl(1)
               errmax = ctl(2)
            end if
         end if
         if (.not. diis_on .and. it > scf_config%ndamp) diis_on = errmax < scf_config%diis_start
         !
         ! A/EDIIS while damping has ended but the commutator is still too large for
         ! CDIIS. CDIIS has no constraint on its coefficients, hence the
         ! `diis_start` gate: far from convergence it can overshoot. ADIIS (Hu &
         ! Yang, JCP 132, 054109) minimises an energy functional over the convex
         ! hull of the history, so it only moves towards a visited density: reliable
         ! where CDIIS is not and slow where CDIIS is fast, hence the two switch
         ! rather than blend (gmshpc: a blended step turned a stretched-water case
         ! converging in 21 CDIIS iterations into one that never converged).
         !
         adiis_active = scf_options%adiis .and. .not. diis_on .and. it > scf_config%ndamp
         mom_active = scf_driver_mom_active(mom_active, scf_options, diis_on, errmax)
         call scf_driver_shift_now(scf_options, diis_on, adiis_active, errmax, shift_guard, cdiis, scf_config%root, shifted)
         ! The new density is rank 0's alone: extrapolation, eigensolve and density
         ! run on one device and the density is broadcast, since dsyevd need not
         ! return the same vectors on two devices. Orbitals and energies come over
         ! once, after the loop.
         diag_failed = .false.
         diag_info = 0
         if (scf_config%root) then
            if (adiis_active) then
               call adiis%push(nao, scf_config%nspin, scf_state%fock, scf_state%dold, scf_state%w1)
               call adiis%extrapolate(nao, scf_config%nspin, scf_state%fock, scf_state%w1)
            end if
            if (diis_on) then
               call cdiis%push(nao, scf_config%nspin, scf_state%fock, scf_state%errv, scf_state%w1)
               call cdiis%extrapolate(nao, scf_config%nspin, scf_state%fock, scf_state%w1, ok)
            end if
            ! --- level shift, F += mu (S - S D S / occ), onto the extrapolated F
            if (shifted) then
               call scf_diis_level_shift(la, nao, scf_config%nspin, scf_options%level_shift, scf_config%occ, &
                      scf_context%smat, scf_result%dmat, scf_state%w1, scf_state%w2, scf_state%fock)
            end if
            orb_mode = scf_driver_occupation_mode(scf_options, mom_active)
            ! --- diagonalise, per spin, and the new density ---------------------
            mom_swaps = 0
            do s = 1, scf_config%nspin
               if (sym_on) then
                  ! MOM is refused with symmetry at the top of this routine.
                  call scf_driver_sym_spin_step(la, scf_config, scf_state, scf_result, sym, s, info_s)
                  nswap_s = 0
               else
                  call scf_driver_spin_step(la, scf_context, scf_config, scf_state, scf_result, scf_options, s, &
                                            orb_mode == ORB_MOM, nswap_s, info_s)
               end if
               mom_swaps(s) = nswap_s
               ! The FIRST failing spin's own info -- scf_driver_note_failure's own header.
               call scf_driver_note_failure(info_s, diag_failed, diag_info)
            end do
            if (.not. diag_failed .and. .not. diis_on .and. .not. adiis_active) then
               call scf_diis_damp(nao, scf_config%nspin, scf_options%damp, scf_state%dold, scf_result%dmat)
            end if
            if (sym_on .and. .not. diag_failed) call scf_driver_symmetrize(sym, scf_result%dmat, scf_state%w1)
         end if
         ! Every rank must learn of a diagonalisation failure the same
         ! cycle it happened (root's own `diag_failed` above is otherwise
         ! invisible to the others), or a non-root rank iterates on into
         ! the next collective while root has already exited below --
         ! the same hang `ctl`'s own broadcast, right above, exists to
         ! avoid for the energy/commutator-norm pair.
         if (present(comm)) then
            if (comm%size() > 1) then
               diag_flag = [merge(1.0_dp, 0.0_dp, diag_failed), real(diag_info, dp), real(mom_swaps, dp)]
               call scf_context_bcast_flat(comm, diag_flag, 4)
               diag_failed = diag_flag(1) > 0.5_dp
               diag_info = nint(diag_flag(2))
               mom_swaps = nint(diag_flag(3:4))
            end if
         end if
         if (diag_failed) then
            ! A malformed Fock matrix (NOT the ordinary "has not converged
            ! yet" case -- dsyev/cuSOLVER could not even diagonalise it):
            ! stop iterating rather than build the next Fock from
            ! eps_s/cmo_s/d_s that la%syev never actually
            ! wrote. scf_result%converged stays .false.; trc_capi.F90's
            ! trc_run_scf already turns that into TRC_ERR_NOCONV, and
            ! scf_result%iterations > 0 keeps it off the "setup failed"
            ! (TRC_ERR_UNSUPPORTED) path -- the SAME two fields, and the
            ! SAME mapping, an ordinary non-converged run already uses.
            scf_result%message = "trc_scf: the Fock-matrix diagonalisation did not converge "// &
                                 "(eigensolver info="//to_string(diag_info)//")"
            scf_result%iterations = it
            exit
         end if
         ! Every rank must hold the same density to the bit. The Fock sum is
         ! broadcast already; the XC potential is accumulated with atomics
         ! whose order differs per rank, 1e-13 in the energy on water, and a
         ! difference that small compounds through a diverging iteration.
         if (present(comm)) then
            if (comm%size() > 1) then
               !$acc update self(scf_result%dmat)
               call scf_context_bcast_flat(comm, scf_result%dmat, n2*scf_config%nspin)
               !$acc update device(scf_result%dmat)
            end if
         end if
         call scf_state%update_errv(nao, scf_config%nspin, scf_result%dmat)
         drms = sqrt(la%dot_3d(scf_state%errv, scf_state%errv))/real(nao, dp)
         ! The linear algebra between the XC build and the next Fock build:
         ! energy, Fock assembly, commutator, DIIS, eigensolve, density and
         ! its broadcast (PCM excluded, it has its own timer).
         scf_timers%tr_it = scf_timers%rest%get_elapsed_time()
         scf_timers%t_rest = scf_timers%t_rest + scf_timers%tr_it
         if (scf_config%talk) then
            if (scf_config%dft) then
               write (line, "(i5,f22.12,es14.4,es14.4,es14.4,3f10.3,a)") &
                  it, scf_energies%total, scf_energies%total - eold, drms, errmax, scf_timers%tf_it, scf_timers%tx_it, &
                  scf_timers%tr_it, merge(" i", "  ", incr_now)
            else
               write (line, "(i5,f22.12,es14.4,es14.4,es14.4,2f10.3,a)") &
                  it, scf_energies%total, scf_energies%total - eold, drms, errmax, scf_timers%tf_it, &
                  scf_timers%tr_it, merge(" i", "  ", incr_now)
            end if
            call global_logger%info(trim(line))
            if (any(mom_swaps > 0)) call global_logger%info("         mom: kept "//to_string(sum(mom_swaps))// &
                                                            " orbital(s) occupied above the aufbau set")
         end if
         scf_result%iterations = it
         ! A real density for the benchmarks (dump_density's comment says why).
         if (it == 2) then
            !$acc update self(scf_result%dmat)
            call dump_density(nao, scf_result%dmat(:, :, 1))
         end if
         ! Converged only on a full build (trc_scf_incremental's rule 1).
         if (incremental%accept(it > 1 .and. abs(scf_energies%total - eold) < scf_options%conv_energy &
                                .and. errmax < scf_options%conv_diis, errmax)) then
            scf_result%converged = .true.
            exit
         end if
         ! --soscf: hand over to trc_scf_soscf_run, its own header, once the commutator is small enough.
         if (scf_config%soscf .and. errmax < scf_config%soscf_start) then
            soscf_switch_iteration = it
            exit
         end if
         eold = scf_energies%total
      end do
      scf_result%n_increments = incremental%n_increments
      scf_result%shift_dropped = shift_guard%dropped
      if (soscf_switch_iteration > 0 .and. .not. scf_result%converged) then
         call trc_scf_soscf_run(scf_basis, scf_context, scf_config, scf_state, scf_result, scf_energies, la, &
                               pair_eri_cache%eri)
         if (scf_state%eps_rewritten) shifted = .false.
      end if
      ! The density the energy was evaluated at, and the orbitals that made it.
      !$acc update self(scf_state%dold, scf_result%cmo, scf_result%eps)
      if (present(comm)) then
         if (comm%size() > 1) then
            call scf_context_bcast_flat(comm, scf_result%cmo, n2*scf_config%nspin)
            call scf_context_bcast_flat(comm, scf_result%eps, nao*scf_config%nspin)
         end if
      end if
      call scf_state%release_device(scf_context, scf_result)
      if (sym_on) then
         scf_result%mo_irrep = sym%mo_irrep
         call sym%release()
      end if
      scf_result%dmat = scf_state%dold
      scf_energies%pcm_qtotal = scf_context%pcm%q_total
      call scf_result%update(scf_energies, scf_state%dold)
      ! The not-converged message (`diag_failed` set a sharper one) and the
      ! HOMO-LUMO gap check, which `mom` and a shifted last iteration skip.
      call scf_driver_eps_notes(scf_config, scf_options, shifted, nao - scf_result%n_lindep, scf_result)
      call scf_driver_report(scf_context, scf_config, scf_timers, scf_energies, scf_result, pair_eri_cache%eri, &
                             scf_dfjk, comm)
      call scf_context%release()

   end subroutine trc_scf_run

   subroutine scf_driver_skeleton_decide(scf_basis, sym, sym_on, scf_config, comm, skeleton, err)
      !! THE SKELETON FOCK BUILD under a symmetric
      !! group: the four-centre engine is built with the point group, launches
      !! only the orbit representatives, and F is symmetrized after the assembly
      !! (as for the full build). Not with density fitting (its own path) nor over several ranks
      !! (the plan is unmeasured there); `sym%skeleton = SYM_SKEL_OFF` keeps the
      !! full build and the symmetrization of F alone, and SYM_SKEL_AUTO (the
      !! default) leaves it to the size rule below. Either way the group, the
      !! frame, the SO basis, the irrep labels and every reported quantity are
      !! the same: only which quartets are launched changes.
      type(trc_basis_t), intent(in) :: scf_basis
      type(trc_sym_scf_t), intent(inout), optional :: sym
      logical, intent(in) :: sym_on   !! `sym` is present with a group of more than one operation
      type(scf_settings_t), intent(in) :: scf_config
      type(comm_t), intent(in), optional :: comm
      type(scf_skeleton_t), intent(out) :: skeleton
      type(error_t), intent(out) :: err
      integer :: mode   !! what was asked for: SYM_SKEL_AUTO/ON/OFF, TRC_SYM_SKELETON applied
      mode = SYM_SKEL_AUTO
      skeleton%norbit = 0
      skeleton%why = "asked for"
      skeleton%on = sym_on
      if (skeleton%on .and. scf_config%df) then
         skeleton%on = .false.
         skeleton%why = "density fitting"
      end if
      if (skeleton%on) then
         mode = sym_run_skeleton_mode(sym%skeleton)
         if (mode == SYM_SKEL_OFF) then
            skeleton%on = .false.
            skeleton%why = "asked for"
         end if
      end if
      if (skeleton%on .and. present(comm)) then
         if (comm%size() > 1) then
            skeleton%on = .false.
            skeleton%why = "several ranks"
         end if
      end if
      if (skeleton%on) then
         call sym%geometry(scf_basis, skeleton%geom, err)
         if (err%has_error()) return
         ! THE SIZE RULE. The skeleton build is a loss on a molecule whose
         ! atom orbits are small, measured at 0.44x the symmetry-off ERI time
         ! on SF6 and 1.41x-2.78x on benzene, naphthalene and pentacene
         ! (sym_run_skeleton_size_rule carries the numbers). The atom
         ! permutations are in hand either way, so the orbit count is free.
         if (mode == SYM_SKEL_AUTO) then
            skeleton%norbit = sym_geom_atom_orbits(skeleton%geom%perm)
            skeleton%on = sym_run_skeleton_size_rule(scf_basis%natm, skeleton%norbit)
            skeleton%why = "size rule: "//to_string(scf_basis%natm)//" atoms in "// &
                           to_string(skeleton%norbit)//" orbits"
         end if
      end if
   end subroutine scf_driver_skeleton_decide

   subroutine scf_driver_sym_prepare(sym, scf_context, scf_config, scf_options, scf_result, ok)
      !! Per-irrep orthogonalisers from the overlap just built, and the
      !! device data of the symmetrisation (see the header of trc_sym_scf).
      !! `ok` false: the run must stop, with `scf_result%message` set and
      !! `sym` and `scf_context` already released.
      type(trc_sym_scf_t), intent(inout) :: sym
      type(scf_ctx_t), intent(inout) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(trc_scf_options_t), intent(in) :: scf_options
      type(trc_scf_result_t), intent(inout) :: scf_result
      logical, intent(out) :: ok
      type(error_t) :: err
      ok = .false.
      call sym%prepare(scf_context%smat, scf_config%nspin, scf_options%lindep_thresh, err)
      if (.not. err%has_error()) call sym%check_occupations(scf_config%nocc(1:scf_config%nspin), err)
      if (err%has_error()) then
         scf_result%message = "trc_scf: "//err%get_message()
         call sym%release()
         call scf_context%release()
         return
      end if
      scf_result%n_lindep = scf_context%nao - sym%nmo
      if (sym%decline_lindep .and. (scf_context%n_lindep > 0 .or. sym%nmo < scf_context%nao)) then
         scf_result%sym_declined = .true.
         scf_result%message = "trc_scf: canonical orthogonalisation drops functions, which the per-irrep "// &
                              "blocks treat differently from C1; symmetry declined"
         call sym%release()
         call scf_context%release()
         return
      end if
      ok = .true.
   end subroutine scf_driver_sym_prepare

   subroutine scf_driver_symmetrize(sym, a, work)
      !! Average each spin's matrix of `a` (a density or a Fock matrix) over
      !! the group, in place; `work` is `sym%symmetrize`'s scratch.
      type(trc_sym_scf_t), intent(in) :: sym
      real(dp), contiguous, intent(inout) :: a(:, :, :), work(:, :)
      integer :: s
      do s = 1, size(a, 3)
         call sym%symmetrize(size(a, 1), a(:, :, s), work)
      end do
   end subroutine scf_driver_symmetrize

   subroutine scf_driver_log_symmetry(sym, skeleton)
      !! The symmetric run's group and which Fock build it uses, in the log.
      type(trc_sym_scf_t), intent(in) :: sym
      type(scf_skeleton_t), intent(in) :: skeleton
      call global_logger%info("   symmetry: "//trim(sym%ao%group%name)//", "// &
         to_string(sym%ao%nirrep)//" irreps, SO dimensions "//dims_string(sym%ao%irrep_dim))
      ! Which Fock build a symmetric run is using is not visible in any
      ! reported number -- both paths give the same F -- so it is said
      ! here, with the size rule's own counts when the rule decided.
      if (skeleton%on) then
         call global_logger%info("   symmetry: skeleton Fock build ("//trim(skeleton%why)//")")
      else
         call global_logger%info("   symmetry: full Fock build ("//trim(skeleton%why)//")")
      end if
   end subroutine scf_driver_log_symmetry

   subroutine scf_driver_guess_spin(la, scf_context, occ, h, nocc_s, dout, info)
      !! Guess density for one spin from a Fock-like `h`, into `dout`: the
      !! eigenproblem and the products on the device, the result copied back
      !! for the resident copy made afterwards. `h` is a host array, copied
      !! up for the region; the loop runs the same sequence on resident ones.
      type(trc_linalg_t), intent(inout) :: la
      type(scf_ctx_t), intent(in) :: scf_context
      real(dp), intent(in) :: occ   !! occupation per orbital (`scf_config%occ`)
      real(dp), intent(in) :: h(scf_context%nao, scf_context%nao)
      integer, intent(in) :: nocc_s
      real(dp), intent(out) :: dout(scf_context%nao, scf_context%nao)
      integer, intent(out) :: info
      real(dp), allocatable :: fp(:, :), c(:, :), t(:, :), e(:), d(:, :)
      integer :: nao
      nao = scf_context%nao
      allocate (fp(nao, nao), c(nao, nao), t(nao, nao), e(nao), d(nao, nao))
      !$acc data copyin(h, scf_context%x) create(fp, c, t, e) copyout(d)
      call scf_driver_diagonalize(la, scf_context, h, t, fp, e, c, info)
      if (info == 0) call la%density(c, nocc_s, occ, d)
      !$acc end data
      dout = d
   end subroutine scf_driver_guess_spin

   subroutine scf_driver_diagonalize(la, scf_context, f, t, fp, eps_s, cmo_s, info)
      !! One spin's eigenproblem from a Fock-like matrix: F' = X^T F X, its
      !! eigenvectors, C = X C'. Everything is on the device already (the
      !! loop's resident arrays, or the data region of `scf_driver_guess_spin`);
      !! `t` and `fp` are the caller's n x n workspaces, which hold nothing on
      !! return that matters.
      !!
      !! `info` comes back non-zero when F' failed to diagonalise (a
      !! malformed Fock matrix -- a density read past the end of an
      !! undersized buffer was one real cause, found and fixed on this
      !! branch): the caller then stops using `cmo_s` for this spin and the
      !! main iteration reports non-convergence instead of continuing on a
      !! matrix la%syev never actually filled in.
      type(trc_linalg_t), intent(inout) :: la
      type(scf_ctx_t), intent(in) :: scf_context
      real(dp), intent(in) :: f(scf_context%nao, scf_context%nao)
      real(dp), intent(inout) :: t(scf_context%nao, scf_context%nao), fp(scf_context%nao, scf_context%nao)
      real(dp), intent(inout) :: eps_s(scf_context%nao), cmo_s(scf_context%nao, scf_context%nao)
      integer, intent(out) :: info
      call la%gemm_mat(f, scf_context%x, t)
      call la%gemm_mat(scf_context%x, t, fp, transa="T")
      call scf_driver_null_shift(la, scf_context, scf_context%nao, fp)
      call la%syev(fp, eps_s, info)
      if (info /= 0) return
      call la%gemm_mat(scf_context%x, fp, cmo_s)
   end subroutine scf_driver_diagonalize

   subroutine scf_driver_scale_columns(n, cmo, ofrac, t)
      !! t(:, j) = cmo(:, j)*sqrt(ofrac(j)): the column-scaled copy of C that
      !! a density with fractional or per-irrep occupations is one GEMM of.
      !! Its own procedure so the kernel names plain arrays, never a
      !! derived-type component.
      integer, intent(in) :: n
      real(dp), intent(in) :: cmo(n, n), ofrac(n)
      real(dp), intent(inout) :: t(n, n)
      integer :: i, j
      do concurrent(i=1:n, j=1:n)
         t(i, j) = cmo(i, j)*sqrt(ofrac(j))
      end do
   end subroutine scf_driver_scale_columns

   subroutine scf_driver_fractional_density(la, eps_s, cmo_s, nelec, omax, t, ofrac, d_s)
      !! D = sum_i o_i c_i c_i^T with the fractional occupations of the
      !! free-atom runs: Aufbau over `eps_s` (`nelec` electrons, `omax` per
      !! orbital), through a column-scaled copy of C in `t`.
      type(trc_linalg_t), intent(inout) :: la
      real(dp), contiguous, intent(inout) :: eps_s(:), cmo_s(:, :), t(:, :), ofrac(:), d_s(:, :)
      real(dp), intent(in) :: nelec, omax
      integer :: nfrac
      !$acc update self(eps_s)
      call frac_occupations(size(eps_s), eps_s, nelec, omax, ofrac, nfrac)
      !$acc update device(ofrac)
      call scf_driver_scale_columns(size(eps_s), cmo_s, ofrac, t)
      ! D = occ * C(:,1:k) C(:,1:k)^T with occ = 1, k = nfrac, on the
      ! column-scaled copy -- the same GEMM la%density does on t with
      ! those arguments (A(1:nao,1:nfrac) is the same block either way).
      call la%density(t, nfrac, 1.0_dp, d_s)
   end subroutine scf_driver_fractional_density

   subroutine scf_driver_spin_step(la, scf_context, scf_config, scf_state, scf_result, scf_options, spin, mom, nswap, &
                                   info)
      !! One spin's step from its Fock matrix to its density: the eigenproblem
      !! (`scf_driver_diagonalize`), then D = occ C(:,1:nocc) C(:,1:nocc)^T,
      !! or, with fractional occupations, `scf_driver_fractional_density`.
      !! With `mom` the occupied columns are first chosen by overlap with the
      !! density `scf_result%dmat` still holds, the previous one
      !! (`scf_driver_mom_weights`); `nswap` counts those aufbau would skip.
      !! Everything is on the device already. `info` as in the former.
      type(trc_linalg_t), intent(inout) :: la
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(trc_scf_options_t), intent(in) :: scf_options
      integer, intent(in) :: spin
      logical, intent(in) :: mom
      integer, intent(out) :: nswap
      integer, intent(out) :: info
      integer :: nao, nocc
      nswap = 0
      nao = scf_context%nao
      nocc = scf_config%nocc(spin)
      call scf_driver_diagonalize(la, scf_context, scf_state%fock(:, :, spin), scf_state%w1, scf_state%w2, &
                                  scf_result%eps(:, spin), scf_result%cmo(:, :, spin), info)
      if (info /= 0) return
      if (mom .and. nocc > 0 .and. nocc < nao) then
         call scf_driver_mom_weights(la, nao, scf_context%smat, scf_result%dmat(:, :, spin), scf_result%cmo(:, :, spin), &
                                     scf_state%w1, scf_state%w2, scf_state%ofrac)
         call scf_driver_mom_permute(nao, nocc, scf_state%ofrac, scf_result%cmo(:, :, spin), scf_result%eps(:, spin), &
                                     scf_state%w1, nswap)
      end if
      if (scf_options%frac_occ) then
         ! `ofrac` is an explicit argument (not read off `scf_state` by the
         ! kernel) so the `do concurrent` below it never names a
         ! derived-type component.
         call scf_driver_fractional_density(la, scf_result%eps(:, spin), scf_result%cmo(:, :, spin), &
                                            merge(scf_options%nelec_frac, real(scf_config%nocc(spin), dp), &
                                                  scf_config%nspin == 1), &
                                            merge(2.0_dp, 1.0_dp, scf_config%nspin == 1), scf_state%w1, &
                                            scf_state%ofrac, scf_result%dmat(:, :, spin))
      else
         ! k <= 0 is la%density's own zeroing branch (do concurrent).
         call la%density(scf_result%cmo(:, :, spin), scf_config%nocc(spin), scf_config%occ, scf_result%dmat(:, :, spin))
      end if
   end subroutine scf_driver_spin_step

   subroutine scf_driver_sym_spin_step(la, scf_config, scf_state, scf_result, sym, spin, info)
      !! `scf_driver_spin_step` for the SO-blocked path: per-irrep
      !! eigenproblems merged by energy (trc_sym_scf), then the density
      !! from the lowest `nocc` orbitals, or, with fixed per-irrep
      !! occupations, from the occupied ones of each irrep. Everything is
      !! on the device already, as there.
      type(trc_linalg_t), intent(inout) :: la
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_state_t), intent(inout) :: scf_state
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(trc_sym_scf_t), intent(inout) :: sym
      integer, intent(in) :: spin
      integer, intent(out) :: info
      call sym%diagonalize(scf_state%fock(:, :, spin), spin, scf_config%occ, scf_result%eps(:, spin), &
                           scf_result%cmo(:, :, spin), info, scf_state%ofrac)
      if (info /= 0) return
      if (sym%fixed) then
         !$acc update device(scf_state%ofrac)
         call scf_driver_scale_columns(size(scf_state%ofrac), scf_result%cmo(:, :, spin), scf_state%ofrac, scf_state%w1)
         call la%density(scf_state%w1, sym%nlast, 1.0_dp, scf_result%dmat(:, :, spin))
      else
         call la%density(scf_result%cmo(:, :, spin), scf_config%nocc(spin), scf_config%occ, scf_result%dmat(:, :, spin))
      end if
   end subroutine scf_driver_sym_spin_step

   subroutine scf_driver_report(scf_context, scf_config, scf_timers, scf_energies, scf_result, eri, scf_dfjk, comm)
      !! The end of a run's accounting: the DFT step table's totals, then,
      !! when the run talks, the timing tables, the continuum's line and the
      !! density fit's own report.
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_settings_t), intent(in) :: scf_config
      type(scf_timers_t), intent(inout) :: scf_timers
      type(scf_energies_t), intent(in) :: scf_energies
      type(trc_scf_result_t), intent(inout) :: scf_result
      type(trc_eri_t), intent(in) :: eri   !! the four-centre engine, for `report_timing`
      type(df_jk_t), intent(inout), optional :: scf_dfjk
      type(comm_t), intent(in), optional :: comm
      character(len=160) :: line
      integer :: rank, nranks
      if (scf_config%dft) then
         ! The run's own totals, read once here: grid build was added after
         ! setup, above. xc is all of scf_fock_build_vxc but VV10 (t_xc includes it),
         ! not trc_xc's points/pairs timers, which miss the functional and,
         ! over ranks, the collectives on V (3.3 of 10.3 s on 16 ranks), so
         ! the rows add up to the total. Points/pairs are in the summary.
         call scf_timers%xcst%add_time(2, scf_timers%t_xc - scf_timers%t_vv10)
         call scf_timers%xcst%add_time(3, scf_timers%t_vv10)
      end if
      ! Before the quiet return: --json and the regression bench read these
      ! whether or not the run talks.
      call scf_driver_fill_timings(scf_context, scf_timers, scf_config, scf_result%timings)
      if (.not. scf_config%talk) return
      rank = 0
      nranks = 1
      if (present(comm)) then
         rank = comm%rank()
         nranks = comm%size()
      end if
      call report_timing(scf_context, scf_timers, scf_config%dft, eri)
      if (scf_config%dft) call scf_timers%xcst%report(scf_context%t_grid + scf_timers%t_xc, rank, nranks)
      if (scf_config%pcm_on) then
         write (line, "(a,f20.12,a,f10.6,a,f10.3,a,f10.4,a)") "   pcm: E_diel = ", scf_energies%pcm, &
            "  q_total = ", scf_context%pcm%q_total, "  setup ", scf_context%t_pcm_setup, " s, ", &
            scf_timers%t_pcm/real(max(1, scf_result%iterations), dp), " s per iteration"
         call global_logger%info(trim(line))
      end if
      if (scf_config%df) call scf_dfjk%report(scf_timers%t_fock, rank, nranks)
   end subroutine scf_driver_report

   function dims_string(dims) result(str)
      !! "[a, b, ...]" for the one-line symmetry summary.
      integer, intent(in) :: dims(:)
      character(len=:), allocatable :: str
      integer :: i
      str = "["
      do i = 1, size(dims)
         if (i > 1) str = str//", "
         str = str//to_string(dims(i))
      end do
      str = str//"]"
   end function dims_string

   subroutine scf_driver_fill_timings(scf_context, scf_timers, scf_config, tm)
      !! The run's stage times into the result, for `--json` and the
      !! regression bench: the totals `report_timing` prints, plus the
      !! first build kept apart.
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_timers_t), intent(in) :: scf_timers
      type(scf_settings_t), intent(in) :: scf_config
      type(trc_scf_timings_t), intent(out) :: tm
      tm%wall = scf_timers%run%get_elapsed_time()
      tm%setup = scf_context%t_setup
      tm%fock_total = scf_timers%t_fock
      tm%fock_first = scf_timers%tf_first
      tm%n_fock = scf_timers%n_fock
      tm%has_grid = scf_config%dft
      tm%has_vv10 = scf_config%dft .and. scf_config%nlc
      tm%has_pcm = scf_config%pcm_on
      tm%has_df = scf_config%df
      tm%grid = scf_context%t_grid
      tm%xc = scf_timers%t_xc - scf_timers%t_vv10
      tm%vv10 = scf_timers%t_vv10
      tm%pcm_setup = scf_context%t_pcm_setup
      tm%pcm_total = scf_timers%t_pcm
   end subroutine scf_driver_fill_timings

   subroutine report_timing(scf_context, scf_timers, dft, eri)
      !! Where the run's time went, through pic's logger at PERFORMANCE,
      !! which its default INFO level shows: the phases, then setup and
      !! the ERI engine's build broken into their steps, then what parity
      !! skipping did with the last Fock build's work list.
      type(scf_ctx_t), intent(in) :: scf_context
      type(scf_timers_t), intent(in) :: scf_timers
      logical, intent(in) :: dft
      type(trc_eri_t), intent(in) :: eri   !! the ERI engine: its build timings and the last build's parity counters
      character(len=160) :: line
      ! `total` is the run's own wall time, so it also holds what no column
      ! times (PCM, the guess density's upload, the final transfers).
      write (line, "(a,8f10.3)") "          ", scf_context%t_setup, scf_context%t_grid, scf_timers%t_fock, scf_timers%t_xc, &
         scf_timers%tx_pts, scf_timers%tx_prs, scf_timers%t_rest, scf_timers%run%get_elapsed_time()
      call global_logger%performance( &
         "   time (s)   setup      grid      fock        xc   (points     pairs)      linalg     total")
      call global_logger%performance(trim(line))
      call global_logger%performance( &
         "   setup (s)  pairs        1e   sph fold       eri    linalg    orthog")
      write (line, "(a,6f10.3)") "          ", scf_context%ts
      call global_logger%performance(trim(line))
      call global_logger%performance( &
         "   eri (s) hgp pairs   schwarz      bins  quartets    upload")
      write (line, "(a,5f10.3)") "          ", eri%t_build
      call global_logger%performance(trim(line))
      ! Parity skipping, from the engine's own
      ! counters: a dead axis of 0 means no mask was taken at all, and a
      ! masked share with no variant launches means the masks were taken but
      ! no class that reached a kernel has a masked variant in this build.
      ! `variant cost` is the same share as variant/masked work but weighted
      ! by the modelled cost of a work item of each class, which is the figure
      ! that tracks the speed-up: counting work items alone reads 70% where
      ! the cost-weighted share is 50%, because the uncovered classes are the
      ! expensive ones (gen_perclass.py's PARITY_RYS_TOP comment).
      call global_logger%performance( &
         "   parity  dead axis         nwork   masked work  variant work  variant launches  variant cost")
      write (line, "(a,i10,3i14,i18,f14.3)") "          ", eri%parity_axis, eri%nwork, eri%nwork_masked, eri%nwork_variant, &
         eri%nlaunch_variant, eri%parity_cost_variant/max(eri%parity_cost_masked, 1.0_dp)
      call global_logger%performance(trim(line))
      if (dft) then
         call global_logger%performance( &
            "   grid (s)   points  partition  batching    upload")
         write (line, "(a,4f10.3)") "          ", scf_context%t_grid_sub
         call global_logger%performance(trim(line))
      end if
   end subroutine report_timing

   pure subroutine frac_occupations(n, eps, nelec, omax, o, nocc_out)
      !! Aufbau over the eigenvalues, `omax` per orbital (2 restricted, 1 per
      !! spin channel), the frontier's degenerate set (within FRAC_TOL)
      !! sharing the remainder evenly. `nocc_out` is the last orbital with
      !! any occupation, for the GEMM.
      integer, intent(in) :: n
      real(dp), intent(in) :: nelec, omax
      real(dp), intent(in) :: eps(n)
      real(dp), intent(out) :: o(n)
      integer, intent(out) :: nocc_out
      real(dp), parameter :: FRAC_TOL = 1.0e-3_dp
      integer :: nfull, i, g0, g1, ninside
      real(dp) :: rem, ef, share
      nfull = int(nelec/omax)
      rem = nelec - omax*real(nfull, dp)
      o = 0.0_dp
      if (nfull > 0) o(1:min(nfull, n)) = omax
      nocc_out = min(nfull, n)
      if (nfull >= n) return
      ! the frontier: the first not-full orbital when electrons are left,
      ! else the last full one (an open shell of even count, carbon 2p^2)
      if (rem > 0.0_dp .or. nfull == 0) then
         ef = eps(nfull + 1)
      else
         ef = eps(nfull)
      end if
      g0 = n + 1
      g1 = 0
      do i = 1, n
         if (abs(eps(i) - ef) < FRAC_TOL) then
            g0 = min(g0, i)
            g1 = max(g1, i)
         end if
      end do
      if (g1 < g0) return
      ninside = count([(i <= nfull, i=g0, g1)])
      share = (omax*real(ninside, dp) + rem)/real(g1 - g0 + 1, dp)
      o(g0:g1) = share
      nocc_out = g1
   end subroutine frac_occupations

   subroutine scf_driver_null_shift(la, scf_context, n, fp)
      !! With functions dropped by canonical orthogonalisation, F' = X^T F X
      !! has a zero eigenvalue per dropped direction, a zero orbital that
      !! aufbau occupied below a positive HOMO. Add shift*(I - X^T S X) to
      !! lift them above every real level (2|F'|_F + 1 bounds the spectrum):
      !! never occupied, always last, at eps(n - n_lindep + 1:n). A no-op
      !! when nothing was dropped.
      type(trc_linalg_t), intent(in) :: la
      type(scf_ctx_t), intent(in) :: scf_context
      integer, intent(in) :: n
      real(dp), intent(inout) :: fp(n, n)
      if (scf_context%n_lindep == 0) return
      call null_shift_add(n, 2.0_dp*sqrt(la%dot_2d(fp, fp)) + 1.0_dp, scf_context%pnull, fp)
   end subroutine scf_driver_null_shift

   subroutine null_shift_add(n, shift, p, a)
      !! a = a + shift*p, on the device under the GPU build.
      integer, intent(in) :: n
      real(dp), intent(in) :: shift
      real(dp), intent(in) :: p(n, n)
      real(dp), intent(inout) :: a(n, n)
      integer :: i, j
      do concurrent(i=1:n, j=1:n)
         a(i, j) = a(i, j) + shift*p(i, j)
      end do
   end subroutine null_shift_add

   subroutine scf_driver_gap_check(scf_config, skip, nao, nmo, eps)
      !! A warning, never a refusal, when the final orbitals' HOMO-LUMO gap
      !! is tiny. A closed-shell molecule run at the wrong charge stalled for
      !! 60 iterations with |FDS - SDF| stuck near 1e-2 before anyone looked
      !! at the charge; a gap this small says to look there first. Rank 0
      !! only; `skip` for fractional occupations (no aufbau gap) or no
      !! iteration run.
      type(scf_settings_t), intent(in) :: scf_config
      logical, intent(in) :: skip
      integer, intent(in) :: nao, nmo
      real(dp), intent(in) :: eps(nao, scf_config%nspin)
      real(dp), parameter :: GAP_WARN = 1.0e-2_dp
      !! 0.27 eV. Closed-shell organics sit well above it even with a pure
      !! functional (singlet O2 at PBE/6-31G*, a degenerate pi* pair half
      !! filled, is 0.039); only a genuinely degenerate or wrongly charged
      !! system gets this close.
      real(dp) :: gap
      character(len=160) :: line
      if (skip .or. .not. scf_config%root) return
      gap = scf_settings_homo_lumo_gap(nao, scf_config%nspin, nmo, scf_config%nocc, eps)
      if (gap < GAP_WARN) then
         write (line, "(a,es10.3,a)") "trc_scf: HOMO-LUMO gap is only ", gap, &
            " Eh; check the charge and multiplicity"
         call global_logger%warning(trim(line))
      end if
   end subroutine scf_driver_gap_check

   subroutine scf_driver_shift_now(opts, diis_on, adiis_active, errmax, guard, cdiis, root, shift)
      !! Whether this iteration's Fock matrix carries the level shift
      !! F += mu (S - S D S / occ). While DIIS and A/EDIIS are off it always
      !! does. With `level_shift_until` it stays on through them as well, onto
      !! the extrapolated F so the DIIS history stays unshifted, while the
      !! commutator is at least `level_shift_until` and the stall guard
      !! (`scf_shift_guard_t`) has not dropped it. In every case an iteration
      !! that can converge (commutator below `conv_diis`) diagonalises an
      !! unshifted F, or the reported virtual eps would carry the shift (a
      !! restart from a converged density inside the damped start did this).
      !! Only the iterations shifted by `level_shift_until` alone feed the
      !! guard; the shift of the damped start is never watched.
      type(trc_scf_options_t), intent(in) :: opts
      logical, intent(in) :: diis_on, adiis_active, root
      real(dp), intent(in) :: errmax
      type(scf_shift_guard_t), intent(inout) :: guard
      type(cdiis_t), intent(inout) :: cdiis
      logical, intent(out) :: shift
      logical :: fired
      shift = .false.
      if (opts%level_shift <= 0.0_dp .or. errmax < opts%conv_diis) return
      shift = .not. diis_on .and. .not. adiis_active
      if (shift .or. opts%level_shift_until <= 0.0_dp .or. guard%dropped .or. errmax < opts%level_shift_until) return
      call guard%observe(errmax, fired)
      shift = .not. fired
      if (.not. fired) return
      ! The history holds ten iterations of a shifted, stalled run. Plain DIIS
      ! extrapolating over it stays parked until they have been pushed out
      ! (O3+ 6-31G*: 83 iterations with the history kept), so start it afresh.
      if (diis_on) call cdiis%restart()
      if (root) call global_logger%info("trc_scf: level_shift_until dropped the level shift after "// &
                                        to_string(guard%n_seen)//" iterations shifted through DIIS: the commutator is not "// &
                                        "below "//to_string(SHIFT_STALL_FACTOR, "f4.2")//" of its earlier best over the last "// &
                                        to_string(SHIFT_STALL_WINDOW)//" ("//to_string(minval(guard%recent), "es9.2")// &
                                        " against "//to_string(guard%older_best, "es9.2")//")")
   end subroutine scf_driver_shift_now

   pure function scf_driver_shift_stalled(recent_best, older_best) result(stalled)
      !! The criterion of the stall guard: the best commutator of the window has
      !! not fallen below `SHIFT_STALL_FACTOR` times the best before it.
      real(dp), intent(in) :: recent_best, older_best
      logical :: stalled
      stalled = recent_best > SHIFT_STALL_FACTOR*older_best
   end function scf_driver_shift_stalled

   pure subroutine scf_driver_shift_observe(this, errmax, fired)
      !! Record the commutator norm of one more iteration shifted only because
      !! of `level_shift_until`. `fired` is true on the call that drops the
      !! shift: from the (`SHIFT_STALL_WINDOW` + 1)th observation on, once the
      !! best of the window fails `scf_driver_shift_stalled` against the best
      !! before it. Nothing is observed after that.
      class(scf_shift_guard_t), intent(inout) :: this
      real(dp), intent(in) :: errmax
      logical, intent(out) :: fired
      integer :: slot
      fired = .false.
      if (this%dropped) return
      slot = mod(this%n_seen, SHIFT_STALL_WINDOW) + 1
      if (this%n_seen >= SHIFT_STALL_WINDOW) this%older_best = min(this%older_best, this%recent(slot))
      this%recent(slot) = errmax
      this%n_seen = this%n_seen + 1
      if (this%n_seen > SHIFT_STALL_WINDOW) fired = scf_driver_shift_stalled(minval(this%recent), this%older_best)
      this%dropped = fired
   end subroutine scf_driver_shift_observe

   pure function scf_driver_mom_active(was_active, opts, diis_on, errmax) result(active)
      !! Whether MOM follows the density this iteration; once true it stays
      !! true (a latch). MOM locks the occupation of the density it starts
      !! from, so it must start from one worth following: DIIS is on AND either
      !! the caller chose the starting density (`opts%mom_follow_guess`: a
      !! density handed in, a restart, a projection) or the commutator is below
      !! `MOM_START`. A density terco builds itself, SAD included, has no claim
      !! to be followed: a restricted run has DIIS on from iteration 1, and
      !! following a poor core guess there converged water/6-31G* to a
      !! determinant 7.9 Eh above the ground state.
      logical, intent(in) :: was_active, diis_on
      type(trc_scf_options_t), intent(in) :: opts
      real(dp), intent(in) :: errmax
      logical :: active
      active = was_active
      if (opts%mom .and. diis_on .and. (opts%mom_follow_guess .or. errmax < MOM_START)) active = .true.
   end function scf_driver_mom_active

   pure function scf_driver_occupation_mode(opts, mom_active) result(mode)
      !! `ORB_FRAC` for the free-atom runs, else `ORB_MOM` once MOM is active
      !! (`scf_driver_mom_active`), else aufbau.
      type(trc_scf_options_t), intent(in) :: opts
      logical, intent(in) :: mom_active
      integer :: mode
      mode = ORB_AUFBAU
      if (opts%frac_occ) then
         mode = ORB_FRAC
      else if (mom_active) then
         mode = ORB_MOM
      end if
   end function scf_driver_occupation_mode

   subroutine scf_driver_eps_notes(scf_config, scf_options, shifted, nmo, scf_result)
      !! What the returned orbital energies are not, then the HOMO-LUMO gap
      !! warning where they can be read as one. With `mom` they need not ascend
      !! across the occupied/virtual boundary: a spin whose LUMO lies below its
      !! HOMO (a negative gap, read off the final `eps` itself, so it holds
      !! after the second-order phase too, which keeps the occupied space but
      !! rewrites `eps`) holds a non-aufbau occupation. After a non-converged
      !! run that ended on a shifted iteration whose `eps` the main loop
      !! produced (`shifted`), the virtuals carry +level_shift. Either way the
      !! gap check would say something false, so it is skipped and the run says
      !! why. Also skipped for fractional occupations and a run with no
      !! iteration.
      type(scf_settings_t), intent(in) :: scf_config
      type(trc_scf_options_t), intent(in) :: scf_options
      logical, intent(in) :: shifted
      integer, intent(in) :: nmo
      type(trc_scf_result_t), intent(inout) :: scf_result
      logical :: skip_gap
      skip_gap = scf_options%frac_occ .or. scf_result%iterations == 0
      if (scf_options%mom .and. .not. skip_gap) then
         if (scf_settings_homo_lumo_gap(size(scf_result%eps, 1), scf_config%nspin, nmo, scf_config%nocc, &
                                        scf_result%eps) < 0.0_dp) then
            skip_gap = .true.
            if (scf_config%root) call global_logger%info("trc_scf: mom kept a non-aufbau occupation; the orbital "// &
                                                         "energies are not in aufbau order and the HOMO-LUMO gap is not checked")
         end if
      end if
      if (.not. scf_result%converged) then
         if (shifted) skip_gap = .true.
         if (len_trim(scf_result%message) == 0) then
            scf_result%message = "trc_scf: not converged"
            if (shifted) scf_result%message = trim(scf_result%message)//" (the last iteration was level shifted: "// &
                                              "the virtual orbital energies include the shift)"
         end if
      end if
      call scf_driver_gap_check(scf_config, skip_gap, size(scf_result%eps, 1), nmo, scf_result%eps)
   end subroutine scf_driver_eps_notes

   subroutine scf_driver_mom_weights(la, nao, smat, d_old, cmo, t, sc, p)
      !! The maximum overlap method's occupation weights (Gilbert, Besley & Gill,
      !! JPCA 112, 13164), here as a projection onto the previous occupied
      !! space,
      !!
      !!   p_j = c_j^T S D_old S c_j,
      !!
      !! the sum of SQUARED overlaps of orbital j with the previous occupied
      !! orbitals times occ, i.e. occ times the squared length of the
      !! projection of c_j onto the previous occupied space. The weight of the paper
      !! is the linear sum of the overlaps themselves, which depends on the
      !! arbitrary sign of each previous orbital (and on how the occupied space
      !! is rotated into orbitals); this one does not. The two agree where it
      !! matters, an orbital inside the previous occupied space (full weight)
      !! or orthogonal to it (zero), and can order partly overlapping orbitals
      !! differently. A damped D_old (not idempotent) only blurs the weights.
      !! Device arrays in; `t` and `sc` (n x n) are workspace and `p` (n)
      !! comes back on the host as well.
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nao
      real(dp), intent(in) :: smat(nao, nao), d_old(nao, nao), cmo(nao, nao)
      real(dp), intent(inout) :: t(nao, nao), sc(nao, nao), p(nao)
      real(dp) :: acc
      integer :: j, k
      ! sc = S C, t = D_old S C, p_j = sum_k sc(k,j) t(k,j)
      call la%gemm_mat(smat, cmo, sc)
      call la%gemm_mat(d_old, sc, t)
      do concurrent(j=1:nao) local(acc, k)
         acc = 0.0_dp
         do k = 1, nao
            acc = acc + sc(k, j)*t(k, j)
         end do
         p(j) = acc
      end do
      !$acc update self(p)
   end subroutine scf_driver_mom_weights

   subroutine scf_driver_mom_permute(nao, nocc, p, cmo, eps, t, nswap)
      !! Reorder the eigenvectors `cmo` (and `eps`) so the `nocc` columns with
      !! the largest weight `p` (host, `scf_driver_mom_weights`) come first;
      !! each block keeps energy order. `cmo`, `eps` and the n x n workspace
      !! `t` are device arrays. `nswap` counts occupied columns that aufbau
      !! would have left empty; nothing moves when it is 0.
      integer, intent(in) :: nao, nocc
      real(dp), intent(in) :: p(nao)
      real(dp), intent(inout) :: cmo(nao, nao), eps(nao), t(nao, nao)
      integer, intent(out) :: nswap
      integer(int_index), allocatable :: order(:)
      integer, allocatable :: perm(:)
      logical, allocatable :: occupied(:)
      real(dp), allocatable :: e_old(:)
      integer :: i, j, k, ij
      nswap = 0
      if (nocc <= 0 .or. nocc >= nao) return
      allocate (order(nao), occupied(nao), perm(nao))
      ! sort_index sorts its argument in place, so sort a copy.
      e_old = p
      call sort_index(e_old, order, reverse=.true.)
      occupied = .false.
      do i = 1, nocc
         occupied(order(i)) = .true.
      end do
      nswap = count(occupied(nocc + 1:nao))
      if (nswap == 0) return
      k = 0
      do i = 1, nao
         if (occupied(i)) then
            k = k + 1
            perm(k) = i
         end if
      end do
      do i = 1, nao
         if (.not. occupied(i)) then
            k = k + 1
            perm(k) = i
         end if
      end do
      !$acc update self(eps)
      e_old = eps
      eps = e_old(perm)
      !$acc update device(eps)
      !$acc enter data copyin(perm)
      ! One flattened index (ij = i + (j - 1) nao), not a two-index loop.
      do concurrent(ij=1:nao*nao) local(i, j)
         i = mod(ij - 1, nao) + 1
         j = (ij - 1)/nao + 1
         t(i, j) = cmo(i, perm(j))
      end do
      do concurrent(ij=1:nao*nao) local(i, j)
         i = mod(ij - 1, nao) + 1
         j = (ij - 1)/nao + 1
         cmo(i, j) = t(i, j)
      end do
      !$acc exit data delete(perm)
   end subroutine scf_driver_mom_permute

   subroutine scf_driver_gwh_matrix(nao, hcore, smat, w1)
      !! The Wolfsberg-Helmholz guess Fock-like matrix, F_ij = 1.75/2 S_ij
      !! (H_ii + H_jj), diagonal H_ii. Host only, plain loops, not `do
      !! concurrent`: `w1` is not yet mapped to the device at this point
      !! in `trc_scf_run` (`scf_state%to_device` runs after the guess is
      !! built), and nvfortran's `-stdpar=gpu` would otherwise try to
      !! launch a kernel on an array `-gpu=mem:separate` has never put
      !! there.
      integer, intent(in) :: nao
      real(dp), intent(in) :: hcore(nao, nao), smat(nao, nao)
      real(dp), intent(inout) :: w1(nao, nao)
      integer :: i, j
      do j = 1, nao
         do i = 1, nao
            w1(i, j) = 0.875_dp*smat(i, j)*(hcore(i, i) + hcore(j, j))
         end do
         w1(j, j) = hcore(j, j)
      end do
   end subroutine scf_driver_gwh_matrix

   subroutine scf_driver_note_failure(info_s, diag_failed, diag_info)
      !! The FIRST failing spin's own diagonalisation info, not whichever
      !! spin ran last -- see this file's own call sites for why that
      !! matters for an unrestricted run.
      integer, intent(in) :: info_s
      logical, intent(inout) :: diag_failed
      integer, intent(inout) :: diag_info
      if (info_s /= 0 .and. .not. diag_failed) diag_info = info_s
      if (info_s /= 0) diag_failed = .true.
   end subroutine scf_driver_note_failure

   subroutine log_incremental_mode(scf_config)
      !! Say so when increments are on: they change the timings. An `i` after
      !! an iteration's line marks a Fock matrix built from the density difference.
      !! Under --df they are off: a warning when the caller asked for them
      !! (trc_set_incremental(1), TRC_INCREMENTAL=1), one info line when they
      !! were only the default, so every --df run is not a warning.
      type(scf_settings_t), intent(in) :: scf_config
      if (scf_config%incr_on) then
         call global_logger%info("   incremental Fock (marked i); converged only on a full build")
      else if (scf_config%incr_df_requested) then
         call global_logger%warning("   incremental Fock builds are not supported with density fitting; "// &
                                    "every iteration is a full build")
      else if (scf_config%incr_off_df) then
         call global_logger%info("   incremental Fock: off with density fitting (always full builds)")
      end if
   end subroutine log_incremental_mode

end module trc_scf_driver
