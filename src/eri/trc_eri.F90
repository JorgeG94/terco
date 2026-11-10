!
! The four-centre Fock build, behind one object.
!
! `trc_eri_t` is built from a basis and a threshold, owns every four-centre
! structure (pair bins, HGP primitive pairs, Schwarz bounds), stays
! device-resident, and replaces three differently-conventioned build routines:
!
!     call eri%build(basis, thresh)
!     do while (.not. converged)
!        call eri%fock(basis, dmat, vhf)
!        ! F = Hcore + vhf; diagonalise, DIIS, new density -- the caller's
!     end do
!     call eri%release()
!
! `fock` returns the two-electron part only, vhf = J - K/2 for the TOTAL
! density D (trace = nelec). NOT 2J - K, which is twice this and is what the
! internal kernel and benchmarks produce: an SCF written against 2J - K
! converged smoothly to -98.50 for water instead of -75.99, since a factor of
! two in the two-electron energy does not diverge. The kernel accumulates a
! folded half-matrix at twice scale, so vhf = (raw + raw^T)/8; both factors
! are applied here.
!
module trc_eri
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use trc_boys, only: boys_init
   use trc_rys_select, only: rys_select_tables_init
   use trc_tables, only: tables_init, LMAX
   use pic_strings, only: to_string
   use trc_cart, only: cart_init
   use trc_hgp, only: hgp_build_pairs, trc_pair_prune_t
   use trc_screen, only: screen_schwarz_bounds_rys, screen_pair_index
   use trc_bins, only: pair_bins_t, bins_build_binned_pairs, bins_build_pair_quartets, bins_pair_image_source, &
                       bins_rep_weights, trc_bin_shells_t, trc_bin_symmetry_t
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t, trc_bin_launch_t
   use trc_error, only: error_t
   use trc_sym_geom, only: trc_sym_geom_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_pairs, only: sym_pairs_shell_map, sym_pairs_symmetrize_q, sym_pairs_symmetrize_dsh, sym_pairs_wanted
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   use pic_timer, only: timer_type
   use pic_profiler, only: profiler_start, profiler_stop
   use trc_binkernel, only: binkernel_fock_bins, binkernel_have_variants
   use trc_parity, only: parity_choose_axis, parity_axperm, parity_ao_perm, parity_mode_env, &
                         PARITY_OFF, PARITY_FORCE_FRAME
   use trc_bin_tiles, only: bin_tiles_build
   use trc_grad, only: grad_bins
   use trc_rys_grad_kernels, only: RYSG_LMAX
   use trc_mdj, only: trc_mdj_t
   use trc_mdg, only: trc_mdg_t, trc_mdg_scales_t
   use trc_binkernel, only: TRC_PRIM_MARGIN
   use trc_basis_types, only: trc_basis_t
   use trc_decontract, only: decon_t, decontract_basis, decontract_expand, decontract_fold, decontract_release, &
                             trc_decontract_options_t, decontract_whole_screen, &
                             decontract_expand_host, decontract_fold_host
   implicit none
   private
   public :: trc_eri_t, eri_lmax_refusal, eri_grad_lmax_refusal, trc_fock_options_t
   public :: eri_mdj_size_rule, TRC_J_ENGINE_AUTO, TRC_J_ENGINE_BINS, TRC_J_ENGINE_MDJ

   integer, parameter :: TRC_J_ENGINE_AUTO = 0
   integer, parameter :: TRC_J_ENGINE_BINS = 1
   integer, parameter :: TRC_J_ENGINE_MDJ = 2

   ! eri_mdj_size_rule's coefficients: its log time ratio, MD J over the bins'
   ! shared J+K pass, fitted to the measurements in its own comment.
   real(dp), parameter :: RULE_C0 = 3.97_dp, RULE_CFILL = 0.464_dp, RULE_CPPP = 1.144_dp
   real(dp), parameter :: RULE_CFPA = -2.227_dp, RULE_CNPAIR = -0.254_dp
   real(dp), parameter :: WHOLE_BINS_CREDIT = 0.157_dp
   !! ln(1.17): the least the bins have measured faster over the split
   !! policy's working basis than over the all-split one the size rule was
   !! fitted on (V100, RHF, 1e-11, cc-pVDZ/aug-cc-pVDZ, C/N/O s kept whole:
   !! water40 1.17x, NA-040 1.19x, gly20 1.31x, aug 1.34-1.44x). eri_build
   !! gives J to MD J over the policy's bins only when the rule predicts
   !! it faster than the all-split bins by more than this.

   integer, parameter :: NBUILD_STAGES = 5  !! timed stages of eri%build, see t_build

   type :: trc_fock_options_t
      !! What `fock_resident_lent`/`fock_resident_many_lent` build, the
      !! optional arguments of `fock_resident` as one value: gmat = j_scale
      !! J[D] - k_scale K[D]/2 of the operator erf(omega r12)/r12 (omega 0:
      !! the Coulomb operator). `fock` takes the same value. Plain values only (see trc_gemm_batch_t's
      !! own comment on nvfortran's structure constructors).
      real(dp) :: j_scale = 1.0_dp
      !! Fraction of the Coulomb term. Zero is what the second pass of a
      !! range-separated hybrid wants: long-range exchange and nothing else.
      real(dp) :: k_scale = 1.0_dp
      !! Fraction of exact exchange: one is Hartree-Fock, zero pure DFT
      !! exchange, a hybrid between. The coupled-perturbed equations pass it
      !! explicitly.
      real(dp) :: omega = 0.0_dp
      !! Range-separation parameter, bohr^-1. Zero builds J and K with the
      !! ordinary Coulomb operator; omega > 0 with erf(omega r12)/r12 instead
      !! (wB97X-V's long-range exchange pass wants k_scale = 0.833 and
      !! usually j_scale = 0). omega < 0 (erfc short-range) is not
      !! implemented. See scripts/gen_perclass.py for the substitution.
      logical :: density_screen = .true.
      !! Weight the Schwarz bound by the density it multiplies before
      !! screening: a quartet whose density elements are all negligible
      !! contributes nothing however large its integral.
      !!
      !! A COUPLED-PERTURBED SOLVE MUST SET FALSE. Its density is a trial
      !! vector the solver drives towards zero, so a screen keyed on the
      !! density's magnitude tightens as the solve proceeds -- the operator
      !! stops being the same linear map from one matvec to the next and the
      !! iteration cannot converge. False recovers the plain Schwarz screen,
      !! which depends on the basis alone and leaves the operator fixed.
      logical :: count_survivors = .false.
      !! `fock_resident` only: also count the quartets that survive
      !! screening, into this%nkept. A full extra pass, so it is opt-in.
      integer :: j_engine = TRC_J_ENGINE_AUTO
      !! Which engine forms J: TRC_J_ENGINE_AUTO (mdj_requested's rule),
      !! TRC_J_ENGINE_BINS or TRC_J_ENGINE_MDJ. TRC_MDJ=0/1 in the
      !! environment overrides it either way. A g basis always takes MD J.
      logical :: skeleton = .false.
      !! Return the SKELETON matrix: the bins launch only the orbit
      !! representatives of one side of each bin pair, weighted (see
      !! `fock`'s skeleton), and the result must be symmetrized by the caller
      !! (trc_sym_scf's symmetrize) in the caller's AO space. Needs `build`
      !! to have been given the point group, and an invariant density.
   end type trc_fock_options_t

   type :: trc_eri_t
      real(dp) :: thresh = 1.0e-10_dp
      integer :: nbas = 0, nao = 0
      type(pair_bins_t) :: bins
      type(trc_shell_pairs_t) :: pairs
      !! The working basis's shells (l, AO offset, atom) and their HGP
      !! primitive pairs: what the kernels actually read. The i-side
      !! exponent `pairs%a` is read by the Rys gradient kernels only; it is
      !! built anyway, a one-time O(npp) cost at geometry-build time.
      type(trc_pair_quartets_t) :: quartets
      !! The same pairs, bin-interleaved (bins_build_pair_quartets): a type
      !! of their own, not components of `bins` -- see that routine's own
      !! comment for why.
      real(dp), allocatable :: dsh(:, :)
      !! Shell-block density bound, refreshed per iteration from D.
      type(trc_basis_t) :: pb
      !! The basis the integrals are actually evaluated over, and the
      !! transform back to the one the caller speaks. On a generally
      !! contracted basis `dec%active` is true, `pb` is the split basis and
      !! `nao_w` is its (larger) function count; the public `nao` stays the
      !! contracted one so no caller sees the difference. On a segmented
      !! basis nothing here is allocated and `nao_w == nao`.
      type(decon_t) :: dec
      integer :: nao_w = 0
      type(trc_mdj_t) :: md
      type(trc_mdg_t) :: mg
      !! The MD exchange pass over the same pairs, for the g classes above TERCO_RYS_G_LTOT.
      logical :: upair = .false.
      !! The unique-pair build (docs/SYMMETRY_UNIQUE_PAIRS.md): `pairs` holds
      !! primitive pairs for the orbit representatives only, the Schwarz bounds
      !! were computed for them alone, and the bin-interleaved tables of the
      !! images were rebuilt from them. Set by `build` when it is given a point
      !! group (TRC_UNIQUE_PAIRS=0 turns it off). The paths that read `pairs`
      !! directly (the non-symmetric Fock build, the gradient) call
      !! `complete_pairs` first.
      integer :: upair_nshell_pairs = 0
      !! Ordered shell pairs whose primitive pairs `build` formed (all of
      !! nshell**2 without the unique-pair build).
      integer(int64) :: upair_npp_stored = 0
      !! Primitive pairs held by `pairs` after `build`.
      integer(int64) :: upair_npp_bins = 0
      !! Primitive pairs of the pairs in the bins: what the full list holds, and
      !! what the bin-interleaved tables are made of.
      logical :: needs_md = .false.
      !! A shell above f: HGP and Rys have no kernel for its classes, so J
      !! comes from trc_mdj and K from trc_mdg whatever the environment says.
      integer :: mdj_force = -1
      !! TRC_MDJ at build time (`mdj_force_env`), rank 0's on every rank.
      real(dp) :: pcut_scale = 1.0_dp
      !! Multiplies every kernel's primitive-quartet cut: below 1 when a
      !! general contraction kept whole needs it (trc_decontract's
      !! `decontract_whole_screen`).
      logical :: has_f = .false.
      !! A shell of f or higher: the J-engine rule of a coupled-perturbed
      !! build and of several densities (see mdj_requested).
      logical :: mdj_extended = .false.
      real(dp) :: mdj_score = 0.0_dp
      !! eri_mdj_size_score on the working basis `build` judged: the
      !! predicted ln(t_MDJ / t_bins); mdj_extended is its sign.
      !! eri_mdj_size_rule's verdict on this basis and geometry, set by `build`:
      !! the bins' shared J+K pass is slower than MD J plus exchange-only bins.
      !! The McMurchie-Davidson J engine over the same bins, built on the
      !! first `coulomb` call: an HGP/Rys-only run never pays for it.
      logical :: on_device = .false.
      integer :: nlaunch = 0
      integer(int64) :: nkept = 0
      !! Quartets that survive both screens in the last build. On the FOCK
      !! path, counted only when fock_resident is asked for it, because there
      !! it costs a full extra enumeration pass. On the GRADIENT path it is
      !! always set: the ordered-survivor compaction already copies its screen
      !! pass's atomic counter back to the host once per batch, so summing it
      !! is free, and `gradient` had no counter of its own at all before S16.
      integer(int64) :: nwork = 0
      !! Work items (quartets before the in-kernel screens) of the last
      !! build, Fock or gradient. For a skeleton build this is the restricted
      !! count; comparing it against the full build's is how a test sees that
      !! the restriction engaged rather than silently falling back.
      integer(int64) :: nwork_masked = 0, nwork_variant = 0
      integer :: nlaunch_variant = 0
      real(dp) :: parity_ops_base = 0.0_dp, parity_ops_masked = 0.0_dp
      !! Parity skipping, of the last bin launch (trc_bin_launch_t): work items
      !! of the segments with the dead z axis of the kernels' frame, the part
      !! of them that ran in a masked kernel variant, and those launches. The
      !! counters a test reads to see that the masks were used and not
      !! silently dropped.
      real(dp) :: parity_cost_masked = 0.0_dp, parity_cost_variant = 0.0_dp
      !! Of the same launch: the modelled unmasked cost of every masked launch
      !! and of the part of it a variant took (trc_bin_launch_t's
      !! `cost_masked`/`cost_variant`). Their ratio is how much of the masked
      !! work the variant set reaches; `nwork_variant / nwork_masked` reads
      !! much higher because it counts an (ss|ss) work item and a (dd|dd) one
      !! the same.
      integer :: parity_mode = -1
      !! Set before `build`: PARITY_OFF, PARITY_AUTO or PARITY_FORCE_FRAME
      !! (trc_parity). -1, the default, takes TRC_PARITY from the environment
      !! (PARITY_AUTO without it). Only an engine built with a point group
      !! has a dead axis at all.
      integer :: parity_axis = 0
      !! The dead axis of the last `build` (1 to 3 of the standard frame), 0
      !! for none: the engine's mask and frame.
      type(comm_t) :: comm
      !! Ranks. Every rank builds the same bins and takes every nranks-th
      !! item of the sorted work list from its rank; the Fock matrix is then
      !! summed across ranks. One rank is rank 0 of 1 and no communicator.
      integer :: rank = 0, nranks = 1
      logical :: distributed = .false.
      real(dp) :: t_build(NBUILD_STAGES) = 0.0_dp
      !! Wall seconds of the last build's stages: HGP pairs, Schwarz bounds,
      !! bins, pair quartets, upload to the device. trc_scf reports them with
      !! the rest of its setup.
      !!
      !! Used to be six: a separate MMD primitive-pair pass fed Schwarz
      !! before the HGP pairs existed. Schwarz is now computed from the HGP
      !! pairs themselves (trc_screen.F90's screen_schwarz_bounds_hgp), so that pass
      !! and its stage are both gone -- see build_structures.
   contains
      procedure :: build   => eri_build
      procedure :: fock      => eri_fock
      procedure :: fock_many  => eri_fock_many
      procedure :: fock_nosym => eri_fock_nosym
      procedure :: fock_resident => eri_fock_resident
      procedure :: fock_resident_many => eri_fock_resident_many
      procedure :: fock_resident_lent => eri_fock_resident_lent
      procedure :: fock_resident_many_lent => eri_fock_resident_many_lent
      procedure :: fock_work_size => eri_fock_work_size
      procedure :: gradient => eri_gradient
      procedure :: coulomb => eri_coulomb
      procedure :: release => eri_release
   end type trc_eri_t

contains

   function eri_lmax_refusal(b) result(msg)
      !! "" when this build's four-centre integrals cover every shell of
      !! `b`; otherwise why not, naming the rebuild that would.
      type(trc_basis_t), intent(in) :: b
      character(len=:), allocatable :: msg
      character(len=*), parameter :: LNAME = "spdfghik"
      integer, parameter :: LMAX_BUILDABLE = 4
      !! The highest TERCO_LMAX the generators support.
      integer :: lb
      msg = ""
      if (b%nshell == 0) return
      lb = maxval(b%sh_l)
      if (lb <= LMAX) return
      msg = "the basis has "//LNAME(lb + 1:lb + 1)//" functions (l = "//to_string(lb)//"), but this build's " &
            //"four-centre integrals stop at "//LNAME(LMAX + 1:LMAX + 1)//" (l = "//to_string(LMAX)//"); "
      if (lb <= LMAX_BUILDABLE) then
         msg = msg//"rebuild with -DTERCO_LMAX="//to_string(lb)
      else
         msg = msg//"terco's integrals go up to g (TERCO_LMAX="//to_string(LMAX_BUILDABLE)//")"
      end if
   end function eri_lmax_refusal

   function eri_grad_lmax_refusal(b) result(msg)
      !! "" when this build has four-centre GRADIENT kernels for every shell
      !! of `b`; otherwise why not. Stricter than eri_lmax_refusal only on a
      !! TERCO_LMAX=4 build: the gradient kernels stop at f (RYSG_LMAX), since
      !! no g gradient has been validated.
      type(trc_basis_t), intent(in) :: b
      character(len=:), allocatable :: msg
      integer :: lb
      msg = eri_lmax_refusal(b)
      if (len(msg) > 0) return
      if (b%nshell == 0) return
      lb = maxval(b%sh_l)
      if (lb <= RYSG_LMAX) return
      msg = "the basis has functions of l = "//to_string(lb)//", but the analytic gradient stops at l = " &
            //to_string(RYSG_LMAX)//" (no gradient kernels above it)"
   end function eri_grad_lmax_refusal

   !
   ! Everything that depends on the geometry and not on the density.
   !
   ! Two structures, both of which the four-centre path needs and neither of
   ! which a caller should have to assemble: the HGP primitive pairs the
   ! kernels read, and the Schwarz vector, now formed from THOSE pairs
   ! directly (trc_screen.F90's screen_schwarz_bounds_hgp) rather than from a second,
   ! MMD-only primitive-pair set built and discarded just to feed it -- see
   ! that routine's header for why this is not only faster but also closes a
   ! scale mismatch the MMD route had.
   !
   subroutine eri_build(this, b, thresh, comm, batch_res, decontract, sym)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      type(comm_t), intent(in), optional :: comm
      !! Ranks to split the quartets over. Bind each rank to its device
      !! (trc_bind_device) BEFORE any object goes to the device.
      integer, intent(in), optional :: batch_res
      !! Buckets per decade of Schwarz bound in the pair bins; 1 is decades.
      type(trc_decontract_options_t), intent(in), optional :: decontract
      !! How the working basis is split (trc_decontract); the defaults if
      !! absent. The integrals over any working basis are exact; the policy
      !! changes the cost, and what the screens drop at `thresh` (see
      !! trc_decontract's DECON_SHELL_COST_DEFAULT for what is guaranteed).
      type(trc_sym_geom_t), intent(in), optional :: sym
      !! The point group of the molecule (trc_sym_geom's sym_geom_detect), with `b`
      !! on its standard frame. Present, and of order above one: the bins are
      !! built closed under the group, with the orbit tables of the skeleton
      !! Fock build (`fock`'s skeleton option). The full Fock builds are
      !! unchanged by it, up to the order of summation. Absent: nothing differs
      !! from a build without symmetry.

      integer :: bres
      type(trc_decontract_options_t) :: decon, split
      logical :: no_decontract
      character(len=8) :: bt_env
      character(len=:), allocatable :: refusal

      ! A shell above the generated kernels has no kernel at all: the
      ! integrals came out as NaN and the run died much later, in the
      ! eigensolver, with nothing pointing at the basis. The SCF and the C
      ! API refuse this with a status first; this is the last line for a
      ! caller that builds the engine directly, which has no error path.
      refusal = eri_lmax_refusal(b)
      if (len(refusal) > 0) error stop "trc_eri: "//refusal

      this%rank = 0
      this%nranks = 1
      this%distributed = .false.
      if (present(comm)) then
         this%comm = comm
         this%rank = comm%rank()
         this%nranks = comm%size()
         this%distributed = this%nranks > 1
      end if
      bres = 1
      if (present(batch_res)) bres = batch_res

      call this%release()
      this%thresh = thresh
      this%nao = b%nao
      this%needs_md = maxval(b%sh_l) > 3
      this%has_f = maxval(b%sh_l) >= 3

      !
      ! SPLIT A GENERAL CONTRACTION BEFORE ANYTHING IS BUILT.
      !
      ! Everything below -- pairs, Schwarz, the HGP primitive pairs, the bins,
      ! the density bound -- is built over `wb`, which is the split basis when
      ! there is one and the caller's basis otherwise. Only the transform of D
      ! on the way in and G on the way out knows the difference, so every
      ! entry point and every caller gets this for free. See trc_decontract
      ! for why it is worth 2.4x on cc-pVDZ.
      !
      ! TRC_NO_DECONTRACT=1 forces the old path, for bisecting a discrepancy.
      !
      bt_env = " "
      call get_environment_variable("TRC_NO_DECONTRACT", bt_env)
      no_decontract = len_trim(bt_env) > 0
      ! The policy is the caller's, or the defaults, and the environment can
      ! override either (decon_options).
      if (present(decontract)) then
         decon = decon_options(decontract)
      else
         decon = decon_options(trc_decontract_options_t())
      end if
      this%mdj_force = mdj_force_env()
      ! Every rank must build the SAME working basis: G is summed over
      ! ranks element by element, and two working bases of the same size
      ! would add up to a wrong matrix with nothing to show for it. The
      ! environment is per process, so rank 0's reading of it is the one.
      if (this%distributed) call share_policy(this%comm, no_decontract, decon, this%mdj_force)
      if (no_decontract) then
         call build_over(this, b, thresh, bres, sym)
         return
      end if
      !
      ! WHICH ENGINE FORMS J DECIDES THE SPLIT, AND THE SPLIT FEEDS THE RULE
      ! THAT DECIDES THE ENGINE.
      !
      ! MD J wants every general contraction split (`md_coulomb`: its cost
      ! goes with primitive pairs), the bins want the policy's own working
      ! basis. eri_mdj_size_rule picks the engine from the working basis's
      ! counts, and it was fitted on the all-split working basis. So the
      ! all-split basis is built first and the rule judged on exactly what
      ! it was fitted on. Only when it then gives J to the bins, and the
      ! policy would keep something whole, is the engine rebuilt over the
      ! policy's basis, with the verdict kept: the rule is not re-asked on
      ! counts it was never fitted to. A g basis (MD exchange), TRC_MDJ=1
      ! and a caller's md_coulomb (no exact exchange) are MD J already, and
      ! TRC_MDJ=0 is the bins already: one pass each.
      !
      if (this%needs_md .or. this%mdj_force == 1) decon%md_coulomb = .true.
      if (decon%md_coulomb .or. this%mdj_force == 0) then
         call decontract_basis(b, this%pb, this%dec, decon)
         call build_over(this, b, thresh, bres, sym)
         return
      end if
      split = decon
      split%md_coulomb = .true.
      call decontract_basis(b, this%pb, this%dec, split)
      call build_over(this, b, thresh, bres, sym)
      block
         type(trc_basis_t) :: pb_policy
         type(decon_t) :: dec_policy
         logical :: same
         call decontract_basis(b, pb_policy, dec_policy, decon)
         same = (dec_policy%active .eqv. this%dec%active)
         if (same .and. dec_policy%active) same = pb_policy%nshell == this%pb%nshell .and. &
                                                  dec_policy%nao_p == this%dec%nao_p
         call decontract_release(dec_policy)
         call pb_policy%release()
         if (same) return
      end block
      ! The rule compares MD J with the bins over the all-split basis; the
      ! bins over the policy's basis are faster than that by at least
      ! WHOLE_BINS_CREDIT, so MD J keeps J only if it wins by more.
      call log_engine_choice(this%mdj_score, this%mdj_score < -WHOLE_BINS_CREDIT)
      if (this%mdj_score < -WHOLE_BINS_CREDIT) return
      call this%release()
      this%thresh = thresh
      this%nao = b%nao
      this%needs_md = maxval(b%sh_l) > 3
      this%has_f = maxval(b%sh_l) >= 3
      call decontract_basis(b, this%pb, this%dec, decon)
      call build_over(this, b, thresh, bres, sym)
      this%mdj_extended = .false.
   end subroutine eri_build

   subroutine log_engine_choice(score, md_j)
      !! Under TRC_BUILD_TIMING, the size rule's score on the all-split
      !! basis and what eri_build made of it.
      real(dp), intent(in) :: score
      logical, intent(in) :: md_j
      character(len=8) :: e
      e = " "
      call get_environment_variable("TRC_BUILD_TIMING", e)
      if (len_trim(e) == 0) return
      call global_logger%info("  [J engine] size score "//to_string(score)//": "// &
                              merge("MD J over the all-split basis  ", "bins over the policy's basis   ", md_j))
   end subroutine log_engine_choice

   subroutine build_over(this, b, thresh, bres, sym)
      !! Everything eri_build does once `this%dec` holds the working basis:
      !! the transform to the device, then the structures over the working
      !! basis (this%pb when split, `b` otherwise).
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bres
      type(trc_sym_geom_t), intent(in), optional :: sym
      !! eri_build's `sym`, passed through to build_structures.
      if (this%dec%active) then
         ! The transform runs inside `fock_resident`, where everything is
         ! device-resident; the maps have to be up there with it.
         ! The containing type first, then its arrays. Mapping a component
         ! without the type that holds it is what nvfortran reports as
         ! "partially present", and the same rule is why this%bins is in the
         ! clause below alongside this%bins%sp_i.
         !$acc enter data copyin(this%dec)
         !$acc enter data copyin(this%dec%f_off, this%dec%f_c, this%dec%f_a, &
         !$acc                   this%dec%t_off, this%dec%t_p, this%dec%t_a)
         this%dec%on_device = .true.
      end if
      if (this%dec%active) then
         call build_structures(this, this%pb, thresh, bres, sym)
      else
         call build_structures(this, b, thresh, bres, sym)
      end if
   end subroutine build_over

   !
   ! Everything eri_build does once the basis to build over has been chosen.
   !
   function decon_options(base) result(opts)
      !! The split policy the build uses: `base` (the caller's, or the
      !! defaults), with TRC_DECON_CUTOFF, TRC_DECON_DIFFUSE (0 or 1),
      !! TRC_DECON_SHELL_COST and TRC_DECON_SCREEN_WHOLE (0 or 1) overriding
      !! it, so a policy can be swept without a rebuild. An unreadable value leaves the field as it was.
      type(trc_decontract_options_t), intent(in) :: base
      type(trc_decontract_options_t) :: opts
      character(len=32) :: e
      integer :: ios, flag
      real(dp) :: x
      opts = base
      e = " "
      call get_environment_variable("TRC_DECON_CUTOFF", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) x
         if (ios == 0) opts%cutoff = x
      end if
      e = " "
      call get_environment_variable("TRC_DECON_DIFFUSE", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) flag
         if (ios == 0) opts%diffuse = flag /= 0
      end if
      e = " "
      call get_environment_variable("TRC_DECON_SHELL_COST", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) x
         if (ios == 0) opts%shell_cost = x
      end if
      e = " "
      call get_environment_variable("TRC_DECON_SCREEN_WHOLE", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) flag
         if (ios == 0) opts%screen_whole = flag /= 0
      end if
   end function decon_options

   subroutine share_policy(comm, no_decontract, decon, mdj_force)
      !! Rank 0's split policy, TRC_NO_DECONTRACT and TRC_MDJ, on every rank.
      type(comm_t), intent(in) :: comm
      logical, intent(inout) :: no_decontract
      type(trc_decontract_options_t), intent(inout) :: decon
      integer, intent(inout) :: mdj_force
      integer, parameter :: NPOLICY = 7
      real(dp) :: buffer(NPOLICY)
      buffer = [merge(1.0_dp, 0.0_dp, no_decontract), decon%cutoff, merge(1.0_dp, 0.0_dp, decon%diffuse), &
             decon%shell_cost, merge(1.0_dp, 0.0_dp, decon%screen_whole), merge(1.0_dp, 0.0_dp, decon%md_coulomb), &
             real(mdj_force, dp)]
      call bcast(comm, buffer, NPOLICY, 0)
      no_decontract = buffer(1) /= 0.0_dp
      decon%cutoff = buffer(2)
      decon%diffuse = buffer(3) /= 0.0_dp
      decon%shell_cost = buffer(4)
      decon%screen_whole = buffer(5) /= 0.0_dp
      decon%md_coulomb = buffer(6) /= 0.0_dp
      mdj_force = nint(buffer(7))
   end subroutine share_policy

   subroutine build_structures(this, b, thresh, bres, sym)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh
      integer, intent(in) :: bres
      type(trc_sym_geom_t), intent(in), optional :: sym

      integer :: i, j
      integer, allocatable :: shmap(:, :), shmask(:)
      type(trc_bin_shells_t) :: bsh
      type(trc_bin_symmetry_t) :: bsym
      integer :: pmode
      logical :: orbits
      logical, allocatable :: want(:, :)
      real(dp) :: uneven
      real(dp), allocatable :: qs(:), one(:)
      type(trc_pair_prune_t) :: prune
      ! Stage timing, kept in this%t_build: `build` is five routines with
      ! very different costs and no way to tell from the outside which one
      ! is the slow one.
      type(timer_type) :: tm

      this%t_build = 0.0_dp
      this%parity_axis = 0
      pmode = this%parity_mode
      if (pmode < 0) pmode = parity_mode_env()
      call tm%start()

      ! All three, not just Boys: the HGP primitive pairs and the Schwarz
      ! bounds now computed from them both read the Cartesian index tables.
      ! Without them the bounds come back identically zero, every pair is
      ! pre-screened away and the Fock matrix is silently empty. That is what
      ! happened the first time this ran.
      call boys_init()
      call tables_init()
      call rys_select_tables_init()
      call cart_init()
      this%nbas = b%nshell
      this%nao_w = b%nao

      this%pairs%nbas = b%nshell
      this%pairs%natm = b%natm
      allocate (this%pairs%sh_l(b%nshell), this%pairs%ao_off(b%nshell), this%pairs%sh_at(b%nshell))
      this%pairs%sh_l   = b%sh_l
      this%pairs%ao_off = b%sh_ao
      this%pairs%sh_at  = b%sh_at

      ! The basis already folded common_fac_sp into its coefficients, so the
      ! separate per-shell factor hgp_build_pairs expects is unity. Passing
      ! the factor twice is the classic way to be wrong by a smooth number
      ! here -- see screen_schwarz_bounds_hgp's header for the MMD route that used
      ! to get this wrong.
      allocate (one(b%nshell))
      one = 1.0_dp

      ! The kernel multiplies each pair by sh_c, so sh_c is the amplitude
      ! the prune judges by -- except on a general contraction kept whole,
      ! which is judged by the coefficients its split would carry, and
      ! whose unevenness tightens the kernels' primitive cut to match
      ! (trc_decontract's `decontract_whole_screen`). The cutoff sits three
      ! orders below the screening threshold.
      allocate (prune%camp(b%maxnp, b%nshell))
      prune%camp = b%sh_c
      prune%cut = thresh*1.0e-3_dp
      this%pcut_scale = 1.0_dp
      if (this%dec%screen_whole) then
         call decontract_whole_screen(b, prune%camp, uneven)
         this%pcut_scale = uneven**(-4)
      end if
      orbits = .false.
      if (present(sym)) orbits = sym%group%order > 1
      this%upair = .false.
      if (orbits) then
         ! The point group first: the unique-pair build forms primitive pairs
         ! and Schwarz bounds for the orbit representatives only.
         call orbit_shell_map(b, sym, shmap)
         this%upair = unique_pairs_requested()
      end if
      if (this%upair) then
         want = sym_pairs_wanted(shmap, b%sh_l)
         call hgp_build_pairs(b%sh_np, b%sh_e, b%sh_c, b%sh_r, one, prune, this%pairs, want)
         this%upair_nshell_pairs = count(want)
         deallocate (want)
      else
         call hgp_build_pairs(b%sh_np, b%sh_e, b%sh_c, b%sh_r, one, prune, this%pairs)
         this%upair_nshell_pairs = b%nshell*b%nshell
      end if
      deallocate (prune%camp)
      this%upair_npp_stored = this%pairs%npp
      this%t_build(1) = tm%get_elapsed_time()
      call tm%start()

      !
      ! THE HGP PAIRS GO ON THE DEVICE ONCE, HERE -- NOT IMPLICITLY, PER
      ! SCHWARZ LAUNCH.
      !
      ! They stay resident for the Fock kernels regardless of what Schwarz
      ! does with them, so entering them before computing the bound (instead
      ! of after, alongside everything else below) costs nothing extra and
      ! is what lets screen_schwarz_bounds_rys read them without a transfer of its
      ! own. The exit pairing these before `eri_release` is unchanged.
      !
      !$acc enter data copyin(this%pairs%off, this%pairs%n, this%pairs%p, this%pairs%r, &
      !$acc                   this%pairs%ra, this%pairs%rb, this%pairs%c)

      allocate (qs(b%nshell*(b%nshell + 1)/2))
      ! screen_schwarz_bounds_rys: the Rys-quadrature diagonal kernel (one thread
      ! per (pair, component pair), trc_schwarz_rys.F90), which replaces
      ! screen_schwarz_bounds_hgp's one-thread-per-pair whole-block scan. Same Q,
      ! see trc_screen.F90's own header; screen_schwarz_bounds_hgp stays for
      ! check_schwarz_rys_vs_hgp's comparison.
      call screen_schwarz_bounds_rys(this%pairs, qs)
      this%t_build(2) = tm%get_elapsed_time()
      call tm%start()

      !
      ! ON A SPLIT BASIS THE BOUND CARRIES THE TRANSFORM'S AMPLITUDE.
      !
      ! `qs` bounds a primitive integral, but what has to stay below
      ! `thresh` is the error in the contracted matrix the caller gets
      ! back, and one contracted quartet is the sum of every primitive
      ! quartet its shells can form. `dec%amp` is the column sum of |C|
      ! that turns one into the other; see build_maps in trc_decontract
      ! for the measurement that says it is worth 1.04e-05 Ha on water.
      !
      ! Scaled here rather than inside screen_schwarz_bounds_rys so that routine
      ! stays a statement about integrals, and because both consumers --
      ! the pair prune in bins_build_binned_pairs and the exact per-quartet
      ! test in the kernel, which reads sp_q -- then get it for free.
      !
      ! `amp` is also there, over the caller's own shells, when nothing was
      ! split but a general contraction was kept whole: it then screens as
      ! its split would have (trc_decontract's `whole_amplitude`).
      if (allocated(this%dec%amp)) then
         do i = 1, b%nshell
            do j = 1, i
               qs(screen_pair_index(i, j)) = qs(screen_pair_index(i, j)) &
                                      *this%dec%amp(i)*this%dec%amp(j)
            end do
         end do
      end if

      bsh%sh_l = b%sh_l
      bsh%sh_np = b%sh_np
      bsh%sh_r = b%sh_r
      if (orbits) then
         ! The bins are closed under the group: Q is made invariant first (the
         ! equivalent bounds differ in the last bits, and the cut on Q and the
         ! size bucket must not tell the members of an orbit apart). With the
         ! unique-pair build the images have no bound of their own (Q = 0
         ! there), and the maximum over the orbit is the representative's.
         call sym_pairs_symmetrize_q(shmap, qs)
         ! Parity skipping: the dead axis of the point group, and the shells
         ! in its plane, which bin apart.
         if (pmode /= PARITY_OFF .and. allocated(sym%perm)) then
            call parity_choose_axis(sym%group, sym%perm, b%sh_at, b%sh_r, this%parity_axis, shmask)
         end if
         ! A mask is only meaningful in a frame whose z is the dead axis: z
         ! itself, or the permuted frame, which is applied only when a masked
         ! kernel exists to profit from it (or a test forces it).
         if (this%parity_axis /= 3 .and. .not. (binkernel_have_variants() .or. pmode == PARITY_FORCE_FRAME)) then
            this%parity_axis = 0
         end if
         bsym%shmap = shmap
         if (this%parity_axis > 0) bsym%shmask = shmask
         call bins_build_binned_pairs(bsh, qs, thresh, this%bins, pp_n=this%pairs%n, res=bres, sym=bsym, orb=this%quartets%orb)
      else
         call bins_build_binned_pairs(bsh, qs, thresh, this%bins, pp_n=this%pairs%n, res=bres)
      end if
      this%mdj_score = eri_mdj_size_score(b%nshell, this%bins%npair, prim_pairs(b, this%bins), b%nao)
      this%mdj_extended = this%mdj_score < 0.0_dp
      this%t_build(3) = tm%get_elapsed_time()
      call tm%start()

      if (orbits) then
         this%quartets%orb%ops = sym%group%ops
         call bins_build_pair_quartets(this%bins, this%pairs, this%quartets)
         call bin_tiles_build(this%quartets%tiles, this%bins, b%nshell, this%quartets%orb)
      else
         call bins_build_pair_quartets(this%bins, this%pairs, this%quartets)
         call bin_tiles_build(this%quartets%tiles, this%bins, b%nshell)
      end if
      ! The engine's frame: the dead axis is the kernels' z (a permutation of
      ! D and G per build, hence only when the axis is not z already). Independent
      ! of the pair-quartet build above (sets quartets%axperm/ao_perm*, which
      ! bins_build_pair_quartets and bin_tiles_build neither read nor write).
      if (this%parity_axis > 0 .and. this%parity_axis /= 3) then
         this%quartets%axperm = parity_axperm(this%parity_axis)
         allocate (this%quartets%ao_perm(b%nao), this%quartets%ao_perm_inv(b%nao))
         call parity_ao_perm(b%nshell, b%sh_l, b%sh_ao, b%nao, this%quartets%axperm, &
                             this%quartets%ao_perm, this%quartets%ao_perm_inv)
      end if
      this%upair_npp_bins = sum(int(this%quartets%n, int64))
      this%t_build(4) = tm%get_elapsed_time()
      call tm%start()

      allocate (this%dsh(b%nshell, b%nshell))
      this%dsh = huge(1.0_dp)*1.0e-30_dp

      deallocate (qs, one)

      ! The pair arrays are already device-resident -- entered above,
      ! before screen_schwarz_bounds_rys ran.
      !$acc enter data copyin(this%pairs%sh_l, this%pairs%ao_off, this%pairs%sh_at, this%dsh, &
      !$acc                   this%bins, this%bins%sp_i, this%bins%sp_j, &
      !$acc                   this%bins%sp_q, &
      !$acc                   this%quartets%z, this%quartets%c, this%quartets%a, &
      !$acc                   this%quartets%px, this%quartets%py, this%quartets%pz, &
      !$acc                   this%quartets%ax, this%quartets%ay, this%quartets%az, &
      !$acc                   this%quartets%abx, this%quartets%aby, this%quartets%abz, &
      !$acc                   this%quartets%n)
      if (allocated(this%quartets%ao_perm)) then
         !$acc enter data copyin(this%quartets%ao_perm, this%quartets%ao_perm_inv)
      end if
      this%on_device = .true.
      this%t_build(5) = tm%get_elapsed_time()
   end subroutine build_structures

   subroutine orbit_shell_map(b, geom, shmap)
      !! The shell map of the point group on the working basis `b`: the AO
      !! maps are built on a Cartesian copy, the space the four-centre build
      !! works in (trc_sym_pairs).
      type(trc_basis_t), intent(in) :: b
      type(trc_sym_geom_t), intent(in) :: geom
      integer, allocatable, intent(out) :: shmap(:, :)
      type(trc_basis_t) :: wbc
      type(trc_sym_ao_t) :: ao
      type(error_t) :: err
      wbc = b
      wbc%spherical = .false.
      ao = sym_ao_build(wbc, geom, err)
      if (err%has_error()) error stop "terco: eri%build with a point group: "//err%get_message()
      call sym_pairs_shell_map(wbc, ao, geom, shmap, err)
      if (err%has_error()) error stop "terco: eri%build with a point group: "//err%get_message()
   end subroutine orbit_shell_map

   logical function unique_pairs_requested() result(res)
      !! The unique-pair build is on for every engine built with a point group;
      !! TRC_UNIQUE_PAIRS=0 forms the primitive pairs and Schwarz bounds of
      !! every shell pair instead (the merged scheme's build), to compare.
      character(len=8) :: e
      e = " "
      call get_environment_variable("TRC_UNIQUE_PAIRS", e)
      res = trim(adjustl(e)) /= "0"
   end function unique_pairs_requested

   subroutine complete_pairs(this)
      !! Give every pair of the bins its own primitive pairs in `pairs`, rebuilt
      !! from its representative's (the centres flipped by the operation's
      !! mask, as the bin-interleaved tables were): what the paths that read
      !! `pairs` through the bins' (sp_i, sp_j) need, the non-symmetric Fock
      !! build and the gradient. Host work, once, appended after the
      !! representatives' entries; the device copy is replaced. A no-op
      !! without the unique-pair build.
      class(trc_eri_t), intent(inout) :: this
      integer :: p, key, ksrc, nbas, npp0, npp1, kk, d, o_src, o_dst, n_src
      integer, allocatable :: n_old(:)
      real(dp) :: sg(3)
      real(dp), allocatable :: pp(:), pc(:), pa(:), pr(:, :), pra(:, :), prb(:, :)

      if (.not. this%upair) return
      nbas = this%pairs%nbas
      npp0 = this%pairs%npp
      n_old = this%pairs%n
      npp1 = npp0
      do p = 1, this%bins%npair
         key = (this%bins%sp_i(p) - 1)*nbas + this%bins%sp_j(p)
         if (n_old(key) > 0) cycle
         call bins_pair_image_source(this%bins%sp_i(p), this%bins%sp_j(p), nbas, n_old, this%quartets%orb, ksrc, sg)
         npp1 = npp1 + n_old(ksrc)
      end do
      allocate (pp(npp1), pc(npp1), pa(npp1), pr(npp1, 3), pra(npp1, 3), prb(npp1, 3))
      pp(1:npp0) = this%pairs%p
      pc(1:npp0) = this%pairs%c
      pa(1:npp0) = this%pairs%a
      pr(1:npp0, :) = this%pairs%r
      pra(1:npp0, :) = this%pairs%ra
      prb(1:npp0, :) = this%pairs%rb
      npp1 = npp0
      do p = 1, this%bins%npair
         key = (this%bins%sp_i(p) - 1)*nbas + this%bins%sp_j(p)
         if (n_old(key) > 0) cycle
         call bins_pair_image_source(this%bins%sp_i(p), this%bins%sp_j(p), nbas, n_old, this%quartets%orb, ksrc, sg)
         n_src = n_old(ksrc)
         o_src = this%pairs%off(ksrc)
         o_dst = npp1
         do kk = 1, n_src
            pp(o_dst + kk) = this%pairs%p(o_src + kk)
            pc(o_dst + kk) = this%pairs%c(o_src + kk)
            pa(o_dst + kk) = this%pairs%a(o_src + kk)
            do d = 1, 3
               pr(o_dst + kk, d) = sg(d)*this%pairs%r(o_src + kk, d)
               pra(o_dst + kk, d) = sg(d)*this%pairs%ra(o_src + kk, d)
               prb(o_dst + kk, d) = sg(d)*this%pairs%rb(o_src + kk, d)
            end do
         end do
         this%pairs%off(key) = o_dst
         this%pairs%n(key) = n_src
         npp1 = npp1 + n_src
      end do
      if (this%on_device) then
         !$acc exit data delete(this%pairs%off, this%pairs%n, this%pairs%p, this%pairs%r, &
         !$acc                  this%pairs%ra, this%pairs%rb, this%pairs%c)
      end if
      call move_alloc(pp, this%pairs%p)
      call move_alloc(pc, this%pairs%c)
      call move_alloc(pa, this%pairs%a)
      call move_alloc(pr, this%pairs%r)
      call move_alloc(pra, this%pairs%ra)
      call move_alloc(prb, this%pairs%rb)
      this%pairs%npp = npp1
      if (this%on_device) then
         !$acc enter data copyin(this%pairs%off, this%pairs%n, this%pairs%p, this%pairs%r, &
         !$acc                   this%pairs%ra, this%pairs%rb, this%pairs%c)
      end if
      this%upair = .false.
   end subroutine complete_pairs

   subroutine skeleton_dsh(this)
      !! Make the density bound dsh invariant under the group (its maximum
      !! over the orbit of each shell pair), which the skeleton build's
      !! screens need: the surviving quartets must be a closed set. Host
      !! code on the device-resident bound (one copy each way).
      class(trc_eri_t), intent(inout) :: this
      if (.not. this%quartets%orb%built) error stop "terco: a skeleton Fock build needs eri%build to be given the point group"
      !$acc update self(this%dsh)
      call sym_pairs_symmetrize_dsh(this%quartets%orb%shmap, this%dsh)
      !$acc update device(this%dsh)
   end subroutine skeleton_dsh

   !
   ! G = 2J - K for a given density.
   !
   ! `dmat` and `gmat` are host arrays; the density goes up and the result
   ! comes back. A device-resident variant belongs here for a caller whose
   ! density already lives on the GPU, but it is a separate entry point rather
   ! than a flag -- confusing the two is a segfault, not a wrong number.
   !
   subroutine eri_fock(this, b, dmat, gmat, opts, hcore)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      !! j*J - k*K/2 for the total density, plus Hcore if given.
      type(trc_fock_options_t), intent(in), optional :: opts
      !! What to build: the scalings j_scale and k_scale, the range-separation
      !! parameter omega, the density screen, the skeleton matrix; see
      !! trc_fock_options_t for each. Absent: Hartree-Fock, Coulomb operator,
      !! density screen on, full Fock matrix. Its `j_engine` and
      !! `count_survivors` are not read here.
      real(dp), intent(in), optional :: hcore(this%nao, this%nao)
      !! Core Hamiltonian. If present the result is the full Fock matrix
      !! rather than its two-electron part -- the caller almost always wants
      !! that, and making them add it themselves is friction for nothing.

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j
      real(dp) :: jfac, kfac, ow, jbins
      type(trc_bin_launch_t) :: launch
      type(trc_fock_options_t) :: o
      logical :: md_j, skel

      if (present(opts)) o = opts
      skel = o%skeleton
      jfac = o%j_scale
      kfac = o%k_scale
      ow = o%omega
      if (.not. (ow >= 0.0_dp)) error stop "terco: omega < 0 or NaN (erfc/short-range) is not implemented"

      n = this%nao
      nw = this%nao_w
      ! On a split basis the kernels work in the primitive basis: the density
      ! goes up through C D C^T and the result comes back through C^T G_p C,
      ! a few tens of megaflops each. See trc_decontract.
      allocate (jmat(1, nw, nw), kmat(1, nw, nw), dwork(1, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      if (this%dec%active) then
         call decontract_expand_host(this%dec, dmat, dp_w)
      else
         dp_w = dmat
      end if
      dwork(1, :, :) = dp_w
      jmat = 0.0_dp
      kmat = 0.0_dp

      if (o%density_screen) then
         call density_blocks(this, dwork, 1)
      else
         this%dsh = 1.0_dp
      end if
      if (skel) then
         if (.not. this%quartets%orb%built) error stop "terco: a skeleton Fock build needs eri%build to be given the point group"
         call sym_pairs_symmetrize_dsh(this%quartets%orb%shmap, this%dsh)
      end if
      ! A g basis: J from the MD J engine below, the bins and MD exchange for
      ! exchange only. MD J has no erf operator, so with omega > 0 the bins and
      ! MD exchange form J as well (an RSH's long-range pass has j_scale = 0).
      md_j = this%needs_md .and. ow == 0.0_dp
      jbins = jfac
      if (md_j) jbins = 0.0_dp

      !$acc enter data copyin(dwork, jmat, kmat)
      !$acc update device(this%dsh)
      launch = bin_launch(this, 1, jbins, kfac, ow, skel)
      call binkernel_fock_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, dwork, jmat, kmat)
      this%nlaunch = launch%nlaunch
      this%nwork = launch%nwork
      this%nwork_masked = launch%nwork_masked
      this%nwork_variant = launch%nwork_variant
      this%nlaunch_variant = launch%nlaunch_variant
      this%parity_ops_base = launch%ops_base
      this%parity_ops_masked = launch%ops_masked
      this%parity_cost_masked = launch%cost_masked
      this%parity_cost_variant = launch%cost_variant
      call mdg_pass(this, b, dwork, jmat, jbins, kfac, ow, skel)
      if (md_j .and. jfac /= 0.0_dp) call add_md_j(this, b, dp_w, jmat, jfac, skel)
      !$acc wait
      !$acc update self(jmat)
      !$acc exit data delete(dwork, jmat, kmat)

      ! Symmetrise in the basis the integrals were formed in, THEN fold. The
      ! other order is a different matrix -- C is not square.
      do j = 1, nw
         do i = 1, nw
            gp_w(i, j) = 0.5_dp*(jmat(1, i, j) + jmat(1, j, i))
         end do
      end do
      if (this%dec%active) then
         call decontract_fold_host(this%dec, gp_w, gmat)
      else
         do j = 1, n
            do i = 1, n
               gmat(i, j) = gp_w(i, j)
            end do
         end do
      end if

      ! vhf = (raw + raw^T)/2.
      !
      ! Determined by measurement, not derivation, and the two disagreed: the
      ! benchmark path reports 2J - K as (raw + raw^T)/4, which would make
      ! J - K/2 equal (raw + raw^T)/8. Against gpu4pyscf on the same density
      ! that is uniformly a factor of four too small -- elementwise ratio
      ! 0.2500 at every element, min = median = max, so a pure scale and not
      ! a structural error. The benchmark and this path must differ by two
      ! somewhere in how the folded halves are counted; the constant here is
      ! the one checked against an external code on more than one system.
      if (present(hcore)) then
         do j = 1, n
            do i = 1, n
               gmat(i, j) = gmat(i, j) + hcore(i, j)
            end do
         end do
      end if
      ! hcore, if it was added, is on every rank: sum the two-electron part only.
      if (this%distributed) then
         if (present(hcore)) gmat = gmat - hcore
         call sum_ranks(this, gmat, n*n)
         if (present(hcore)) gmat = gmat + hcore
      end if

      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock

   type(trc_bin_launch_t) function bin_launch(this, ndens, jfac, kfac, omega, skeleton) result(launch)
      !! One launch over this object's bins, for `ndens` densities in the
      !! working basis: the threshold, the ranks and the operand shape are
      !! the object's own; the operator is the caller's.
      class(trc_eri_t), intent(in) :: this
      integer,  intent(in) :: ndens
      real(dp), intent(in) :: jfac, kfac, omega
      logical, intent(in), optional :: skeleton
      if (present(skeleton)) launch%skeleton = skeleton
      launch%thresh = this%thresh
      launch%jfac = jfac
      launch%kfac = kfac
      launch%omega = omega
      launch%ndens = ndens
      launch%nao = this%nao_w
      launch%rank = this%rank
      launch%nranks = this%nranks
      launch%pcut_scale = this%pcut_scale
   end function bin_launch

   logical function mdj_requested(force, pure_functional, size_rule, j_engine) result(res)
      !! Whether the Coulomb part of a resident Fock build comes from the
      !! McMurchie-Davidson J engine instead of the bins. TRC_MDJ=1 forces it
      !! on and TRC_MDJ=0 off (`force`, read once at build time by
      !! `mdj_force_env` and made the same on every rank), so both stay
      !! measurable; then `j_engine` (trc_fock_options_t's), if it names an
      !! engine; otherwise on when the functional is pure (no exact exchange,
      !! so the bins are skipped altogether) or when `size_rule` says so. The
      !! callers pass eri_mdj_size_rule's verdict (eri%mdj_extended) for a
      !! density-screened build; a coupled-perturbed build (density_screen =
      !! .false.) keeps the rule that predates it, MD J for an f basis
      !! (eri%has_f), so its numbers do not move. See eri_mdj_size_rule for
      !! the measurements.
      integer, intent(in) :: force
      logical, intent(in) :: pure_functional, size_rule
      integer, intent(in) :: j_engine
      select case (force)
      case (1)
         res = .true.
      case (0)
         res = .false.
      case default
         select case (j_engine)
         case (TRC_J_ENGINE_MDJ)
            res = .true.
         case (TRC_J_ENGINE_BINS)
            res = .false.
         case default
            res = pure_functional .or. size_rule
         end select
      end select
   end function mdj_requested

   integer function mdj_force_env() result(force)
      !! TRC_MDJ as `mdj_requested` takes it: 1, 0, or -1 when unset.
      character(len=8) :: e
      e = " "
      call get_environment_variable("TRC_MDJ", e)
      select case (trim(e))
      case ("1")
         force = 1
      case ("0")
         force = 0
      case default
         force = -1
      end select
   end function mdj_force_env

   pure logical function eri_mdj_size_rule(nshell, npair, nprim_pairs, nao) result(res)
      !! Whether MD J plus exchange-only bins beats the bins' shared J+K
      !! pass for a hybrid or HF build, from four counts `build` has anyway:
      !! the working basis's shells and functions, its Schwarz-significant
      !! shell pairs and their primitive pairs (before the primitive prune).
      !!
      !! In the bins J and K share one quartet list, so J costs only the
      !! quartets that are Coulomb- but not exchange-significant. On a
      !! compact cluster that is few of them, and the bins win; along an
      !! extended chain it is most, and MD J, whose Hermite stage is cheap
      !! per primitive quartet, wins. The fit is a power law in
      !!   fill = npair / (nshell (nshell + 1) / 2)   compactness: K shares J's quartets
      !!   ppp  = nprim_pairs / npair                 contraction: MD J pays per primitive pair
      !!   fpa  = nao / nshell                        functions per shell: bins pay per Cartesian quartet
      !!   npair                                      size
      !! and MD J is chosen when its predicted log time ratio is negative:
      !!   3.97 + 0.464 ln(fill) + 1.144 ln(ppp) - 2.227 ln(fpa) - 0.254 ln(npair) < 0
      !! Every sign is the one the cost argument predicts. The coefficients
      !! are a linear programme's (hinge loss weighted by |ln ratio|) over the
      !! RHF and B3LYP rows below without f; the f rows were not fitted.
      !!
      !! Measured, Tesla V100-DGXS-32GB, RHF, --thresh 1e-11, full builds
      !! (--no-incremental), mean Fock s/build of builds 2-4, TRC_MDJ=0 / 1;
      !! `*` marks the engine the rule picks:
      !!
      !!   cc-pVDZ:
      !!     water10 240 0.027* / 0.037;  water20 480 0.137* / 0.169
      !!     water30 720 0.426* / 0.478;  water40 960 0.991* / 1.033
      !!     gly5 379 0.096* / 0.110;  gly10 734 0.437 / 0.405*
      !!     gly15 1089 1.032 / 0.834*;  gly20 1444 1.946 / 1.349*
      !!     gly50 3574 14.640 / 5.477*;  NA-040 886 2.106 / 1.968*
      !!   def2-SVP:
      !!     water10 240 0.017* / 0.030;  water20 480 0.089* / 0.130
      !!     water30 720 0.264* / 0.360;  water40 960 0.606* / 0.772
      !!     gly5 379 0.052* / 0.075;  gly10 734 0.213* / 0.259
      !!     gly15 1089 0.508* / 0.543;  gly20 1444 0.931 / 0.881*
      !!     gly50 3574 6.448 / 3.720*;  NA-040 886 1.078* / 1.263
      !!   6-31G*:
      !!     water10 190 0.012* / 0.025;  water20 380 0.047* / 0.085
      !!     water30 570 0.124* / 0.213;  water40 760 0.273* / 0.444
      !!     gly5 349 0.039* / 0.069;  gly10 679 0.151* / 0.232
      !!     gly15 1009 0.362* / 0.489;  gly20 1339 0.665* / 0.814
      !!     gly50 3319 4.449 / 3.559*
      !!   def2-SVPD:
      !!     water10 390 0.149* / 0.167;  water20 780 1.659 / 1.650*
      !!     water30 1170 6.779 / 6.162*;  water40 1560 18.502 / 15.415*
      !!     gly5 574 0.418* / 0.431;  gly10 1109 1.963 / 1.801*
      !!     gly15 1644 4.724 / 3.892*;  gly20 2179 8.848 / 6.463*
      !!     gly50 5389 72.466 / 27.377*
      !!   aug-cc-pVDZ:
      !!     water10 410 0.253* / 0.260;  water20 820 2.851 / 2.588*
      !!     water30 1230 12.225 / 9.714*;  water40 1640 36.399 / 24.815*
      !!     gly5 636 1.065 / 0.970*;  gly10 1231 5.445 / 4.321*
      !!     gly15 1826 13.995 / 9.690*;  gly20 2421 27.715 / 16.773*
      !!     gly50 5991 280.655 / 87.308*
      !!   def2-TZVP:
      !!     water10 430 0.144* / 0.186;  water20 860 0.862 / 0.899*
      !!     gly5 753 0.988 / 0.993*;  gly10 1463 4.470 / 3.939*
      !!     gly20 2883 19.396 / 12.828*
      !!   cc-pVTZ:
      !!     water10 580 0.359 / 0.350*;  water20 1160 2.420 / 2.042*
      !!     gly5 868 1.726 / 1.458*;  gly10 1678 7.644 / 5.561*
      !!     gly20 3298 34.297 / 17.174*
      !!   cc-pVDZ B3LYP:
      !!     water20 480 0.141* / 0.173;  gly10 734 0.441 / 0.427*
      !!     gly20 1444 1.954 / 1.391*;  gly50 3574 14.694 / 5.533*
      !!   def2-SVP B3LYP:
      !!     water40 960 0.618* / 0.806;  gly20 1444 0.936 / 0.910*
      !!   cc-pVDZ PBE:
      !!     water20 480 0.085 / 0.067*;  gly10 734 0.246 / 0.125*
      !!     gly20 1444 1.109 / 0.474*;  gly50 3574 9.040 / 2.901*
      !!   def2-SVP PBE:
      !!     water40 960 0.356 / 0.338*;  gly20 1444 0.531 / 0.374*
      !!
      !! The table's nao is the caller's basis; the rule's counts are the
      !! working basis's (larger on a split general contraction). The rule
      !! picks the faster engine in 67 of 69 rows; it misses water20
      !! def2-TZVP (4.4% slower) and gly5 def2-TZVP (0.5%). The rule it
      !! replaces (MD J for an f basis or a pure functional, bins otherwise)
      !! lost 3.2x at gly50 aug-cc-pVDZ, 2.7x at gly50 cc-pVDZ and 29% at
      !! water10 def2-TZVP. B3LYP's exchange is HF's, and its rows follow the
      !! HF ones; a pure functional (PBE) always takes MD J, which wins every
      !! PBE row. The two engines' iteration-4 energies agree to 3.5e-9
      !! (gly50 6-31G*; most below 1e-9), and converged (--econv 1e-12)
      !! NA-040 cc-pVDZ gives -2270.555778831719 by MD J (this rule) against
      !! -2270.555778831765 by the bins, with 1e-14 references of
      !! -2270.555778834471 (bins) and -2270.555778834457 (MD J): the
      !! threshold's own 2.7e-9 is the same with either engine.
      integer, intent(in) :: nshell, npair, nao
      real(dp), intent(in) :: nprim_pairs
      res = eri_mdj_size_score(nshell, npair, nprim_pairs, nao) < 0.0_dp
   end function eri_mdj_size_rule

   pure real(dp) function eri_mdj_size_score(nshell, npair, nprim_pairs, nao) result(score)
      !! eri_mdj_size_rule's predicted ln(t_MDJ / t_bins); the rule is its
      !! sign. Zero (no verdict) for an empty basis.
      integer, intent(in) :: nshell, npair, nao
      real(dp), intent(in) :: nprim_pairs
      real(dp) :: fill, ppp, fpa
      score = 0.0_dp
      if (nshell < 1 .or. npair < 1) return
      fill = real(npair, dp)/(0.5_dp*real(nshell, dp)*real(nshell + 1, dp))
      ppp = nprim_pairs/real(npair, dp)
      fpa = real(nao, dp)/real(nshell, dp)
      score = RULE_C0 + RULE_CFILL*log(fill) + RULE_CPPP*log(ppp) + RULE_CFPA*log(fpa) &
              + RULE_CNPAIR*log(real(npair, dp))
   end function eri_mdj_size_score

   real(dp) function prim_pairs(b, bins) result(n)
      !! Primitive pairs of `bins`' shell pairs, before the primitive prune.
      type(trc_basis_t), intent(in) :: b
      type(pair_bins_t), intent(in) :: bins
      integer :: s
      n = 0.0_dp
      do s = 1, bins%npair
         n = n + real(b%sh_np(bins%sp_i(s))*b%sh_np(bins%sp_j(s)), dp)
      end do
   end function prim_pairs

   subroutine mdg_pass(this, b, dwork, jmat, jfac, kfac, omega, skeleton)
      !! The g classes binkernel_fock_bins leaves to MD exchange (those above
      !! TERCO_RYS_G_LTOT), into the same folded matrix (both on the device).
      !! Builds the MD pair list and its sorted groups on first use.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)    :: dwork(:, :, :)
      real(dp), intent(inout) :: jmat(:, :, :)
      real(dp), intent(in)    :: jfac, kfac
      real(dp), intent(in)    :: omega
      !! the operator's erf range parameter, as binkernel_fock_bins takes it (0: Coulomb)
      logical, intent(in), optional :: skeleton
      !! Only the orbit representatives as bra, weighted by their orbit size: the
      !! skeleton matrix, as in the bins. Needs the object built with the point group.
      type(trc_mdg_scales_t) :: sc
      logical :: skel
      skel = .false.
      if (present(skeleton)) skel = skeleton
      if (.not. this%needs_md) return
      call ensure_md(this, b)
      ! The representative groups are formed when first asked for; an object
      ! first used for full builds is rebuilt then.
      if (this%mg%built .and. skel .and. .not. this%mg%has_reps) call this%mg%release()
      if (.not. this%mg%built) then
         ! Exchange in spherical components when the caller's basis is
         ! spherical; TRC_MDG_CART=1 keeps it Cartesian, for comparing.
         block
            character(len=8) :: ec
            logical :: sph
            ec = " "
            call get_environment_variable("TRC_MDG_CART", ec)
            sph = b%spherical .and. trim(ec) /= "1"
            if (this%dec%active) then
               call this%mg%build(this%md, this%pb, sph, skel)
            else
               call this%mg%build(this%md, b, sph, skel)
            end if
         end block
      end if
      sc%thresh = this%thresh
      sc%pcut = this%thresh*TRC_PRIM_MARGIN*this%pcut_scale
      sc%omega = omega
      sc%jfac = jfac
      sc%kfac = kfac
      sc%rank = this%rank
      sc%nranks = this%nranks
      sc%skeleton = skel
      call this%mg%fock(this%md, this%dsh, dwork, jmat, sc)
   end subroutine mdg_pass

   subroutine ensure_md(this, b)
      !! The MD pair list, built on first use. An object built with the point
      !! group gives it the pair orbit weights, for the skeleton bra.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      if (this%md%built) return
      if (this%quartets%orb%built) then
         if (this%dec%active) then
            call this%md%build(this%pb, this%bins, this%thresh, bins_rep_weights(this%bins, this%quartets%orb))
         else
            call this%md%build(b, this%bins, this%thresh, bins_rep_weights(this%bins, this%quartets%orb))
         end if
      else if (this%dec%active) then
         call this%md%build(this%pb, this%bins, this%thresh)
      else
         call this%md%build(b, this%bins, this%thresh)
      end if
   end subroutine ensure_md

   subroutine add_md_j(this, b, dp_w, jmat, jfac, skeleton)
      !! jmat(1, :, :) += jfac J[dp_w] from the MD engine; dp_w is a host
      !! array (eri_fock's), jmat is on the device. This rank's share of it,
      !! as in eri_fock_resident: the caller sums the ranks.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)    :: dp_w(:, :)
      real(dp), intent(inout) :: jmat(:, :, :)
      real(dp), intent(in)    :: jfac
      logical, intent(in), optional :: skeleton
      !! the skeleton J: the orbit representatives' bra only (trc_mdj)
      real(dp), allocatable :: jw(:, :)
      integer :: nw, i, j
      call ensure_md(this, b)
      nw = size(dp_w, 1)
      allocate (jw(nw, nw))
      !$acc enter data copyin(dp_w) create(jw)
      call this%md%coulomb(nw, dp_w, jw, this%rank, this%nranks, skeleton)
      !$acc update self(jw)
      !$acc exit data delete(dp_w, jw)
      ! jmat on the device; the fold symmetrises (jmat + jmat^T)/2, and J is
      ! symmetric, so adding it once is adding it once.
      !$acc update self(jmat)
      do j = 1, nw
         do i = 1, nw
            jmat(1, i, j) = jmat(1, i, j) + jfac*jw(i, j)
         end do
      end do
      !$acc update device(jmat)
      deallocate (jw)
   end subroutine add_md_j

   subroutine add_scaled_device(n, f, x, y)
      !! y = y + f x on the device.
      integer,  intent(in)    :: n
      real(dp), intent(in)    :: f
      real(dp), intent(in)    :: x(n)
      real(dp), intent(inout) :: y(n)
      integer :: i
      do concurrent(i = 1:n)
         y(i) = y(i) + f*x(i)
      end do
   end subroutine add_scaled_device

   subroutine eri_coulomb(this, b, dmat, jmat)
      !! J for a density by the McMurchie-Davidson J engine, which never forms
      !! an integral quartet (trc_mdj). Same conventions as `fock` at
      !! k_scale = 0: host arrays, the caller's basis, J of the total density.
      !!
      !! Every rank computes the whole of J; the engine does not split work
      !! over ranks yet, so there is nothing to sum.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: jmat(this%nao, this%nao)

      real(dp), allocatable :: dp_w(:, :), jp_w(:, :)
      integer :: nw

      call ensure_md(this, b)
      nw = this%nao_w
      allocate (dp_w(nw, nw), jp_w(nw, nw))
      if (this%dec%active) then
         call decontract_expand_host(this%dec, dmat, dp_w)
      else
         dp_w = dmat
      end if
      !$acc enter data copyin(dp_w) create(jp_w)
      call this%md%coulomb(nw, dp_w, jp_w)
      !$acc update self(jp_w)
      !$acc exit data delete(dp_w, jp_w)
      if (this%dec%active) then
         call decontract_fold_host(this%dec, jp_w, jmat)
      else
         jmat = jp_w
      end if
      deallocate (dp_w, jp_w)
   end subroutine eri_coulomb

   !
   ! The two-electron part of the nuclear gradient, Rys only: there is no
   ! HGP gradient. `dens` is (nao, nao, ndens): ndens=1 is RHF with the
   ! total density; ndens=2 is UHF with Da and Db. `grad` is (3, natm), the
   ! two-electron contribution ONLY -- the one-electron
   ! derivative terms, the overlap (energy-weighted density) term, nuclear
   ! repulsion, XC and PCM are the caller's, same scope line eri_fock draws
   ! around Hcore.
   !
   subroutine eri_gradient(this, b, density, grad, k_scale, j_scale, omega, skeleton)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: density(:, :, :)
      !! (nao, nao, ndens): ndens = size(density, 3) is read from this, not
      !! passed separately, so the call stays eri%gradient(bas, density, grad)
      !! with no extra integer for the caller to keep in sync.
      real(dp), intent(out) :: grad(3, b%natm)
      real(dp), intent(in), optional :: k_scale, j_scale
      !! Same meaning as eri_fock's: a hybrid functional's exchange fraction
      !! multiplies kfac; default one is Hartree-Fock.
      real(dp), intent(in), optional :: omega
      !! Range-separation parameter, bohr^-1; same meaning as eri_fock's
      !! own omega. Absent or zero (the default) is the ordinary Coulomb
      !! operator; omega > 0 switches to erf(omega*r12)/r12 through the
      !! same rho -> rho_theta substitution eri_fock's Fock kernels make
      !! (scripts/gen_perclass.py), threaded here into the Rys gradient
      !! kernels' own per-root coefficients (scripts/gen_rys_grad.py's
      !! `_RANGE_SEP_BLOCK`). omega == 0 is bit-identical to the
      !! unattenuated path -- that substitution collapses to an exact
      !! multiply by one at every step. omega < 0 (erfc/short-range) is
      !! not implemented, same restriction as eri_fock's own.
      logical, intent(in), optional :: skeleton
      !! Return the SKELETON gradient (docs/SYMMETRY_GRADIENT.md): only the
      !! orbit representatives of one side of each bin pair are launched,
      !! weighted by their orbit size, so `grad` is NOT a gradient until the
      !! caller projects it onto the symmetric subspace
      !! (`sym_scf_grad_symmetrize` of trc_sym_scf, which trc_scf_gradient
      !! already applies). Unlike `fock`'s skeleton that projection is not a
      !! correction but the whole of the symmetry: measured, the raw skeleton
      !! gradient differs from the full one by 0.58 to 1.84 against
      !! components of 0.27 to 2.14. Requires `build` to have been given the
      !! point group and `density` to be invariant under it (the density
      !! bound is then symmetrized too). Default false.
      !!
      !! NOT ALWAYS CHEAPER, and there is no size rule for it (the skeleton
      !! Fock build has none either). When every pair orbit is a singleton --
      !! a PLANAR Cs molecule, where the only operation is the molecular
      !! plane and every atom is its own image -- the restriction drops
      !! nothing, while a diagonal bin pair still becomes the rectangle
      !! reps(A) x A at weight 1/2 instead of the triangle, so the plan grows
      !! slightly: measured on HOCl/6-31G*, 9355 work items against the full
      !! plan's 9142 (+2.3%), with the raw skeleton gradient already equal to
      !! the full one. check_sym_grad's test_cs_planar pins that.

      real(dp), allocatable :: dwork(:, :, :)
      integer :: ndens, nw, id
      real(dp) :: jfac, kfac, ow
      logical :: skel
      type(trc_bin_launch_t) :: launch
      if (this%needs_md) error stop "terco: g functions reach eri%fock and eri%fock_resident only, for now"

      ow = 0.0_dp
      if (present(omega)) ow = omega
      if (.not. (ow >= 0.0_dp)) error stop "terco: eri%gradient: omega < 0 or NaN (erfc/short-range) is not implemented"
      skel = .false.
      if (present(skeleton)) skel = skeleton
      if (skel .and. .not. this%quartets%orb%built) error stop &
         "terco: a skeleton gradient needs eri%build to be given the point group"

      ndens = size(density, 3)

      ! A shell whose centre matched no atom exactly (basis_build leaves
      ! sh_at = 0) has nowhere to put its force: grad_bins would drop or
      ! misindex its contributions. A ghost function is fine for an energy,
      ! not for a gradient; a recomputed centre is a caller bug either way.
      if (any(this%pairs%sh_at == 0)) error stop &
         "terco: eri%gradient: a shell centre matches no atom coordinate exactly"
      ! The gradient kernels read the primitive pairs through the bins' pair list.
      call complete_pairs(this)

      jfac = 1.0_dp
      kfac = 1.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale

      ! On a general contraction the kernels run over the split (working)
      ! basis, same as eri_fock's dp_w. A shell -> atom map is unaffected by
      ! the split -- a split shell keeps its parent's centre, trc_decontract
      ! never moves one -- so `grad` needs no fold-back the way a two-index
      ! G does; only the density going IN has to be expanded.
      nw = this%nao_w
      allocate (dwork(ndens, nw, nw))
      do id = 1, ndens
         if (this%dec%active) then
            call decontract_expand_host(this%dec, density(:, :, id), dwork(id, :, :))
         else
            dwork(id, :, :) = density(:, :, id)
         end if
      end do

      grad = 0.0_dp
      call density_blocks(this, dwork, ndens)
      ! The skeleton's screens must drop a representative and all its images
      ! together, so the bound they test has to be invariant: the surviving
      ! quartet set is then closed under the group, which is what makes the
      ! orbit weights exact. The host call, as eri_fock's skeleton makes it:
      ! density_blocks has just written dsh on the HOST and the copy up is
      ! below, so skeleton_dsh (which starts by copying dsh back DOWN) would
      ! overwrite it with the previous build's bound.
      if (skel) call sym_pairs_symmetrize_dsh(this%quartets%orb%shmap, this%dsh)

      !$acc enter data copyin(dwork) copyin(grad)
      !$acc update device(this%dsh)
      launch = bin_launch(this, ndens, jfac, kfac, ow, skel)
      call grad_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, dwork, grad)
      !$acc wait
      !$acc update self(grad)
      !$acc exit data delete(dwork, grad)
      deallocate (dwork)
      ! The gradient's own work and survivor counts. The Fock entries have
      ! reported these since S15; without them a skeleton gradient that
      ! silently fell back to the full plan would look identical from
      ! outside, which is the failure mode check_sym_grad's "the skeleton is
      ! smaller" assertion exists to catch.
      this%nwork = launch%nwork
      this%nkept = launch%nkept

      if (this%distributed) call sum_ranks(this, grad, 3*b%natm)
   end subroutine eri_gradient

   subroutine sum_ranks(this, g, n)
      !! Sum a matrix over the ranks, in place, every rank receiving the total.
      !!
      !! Then rank 0's copy is broadcast. An allreduce is not obliged to leave
      !! the same bits on every rank -- 1e-13 apart in practice -- and an SCF
      !! iterated from bit-different matrices on different ranks drifts apart,
      !! slowly when it converges and catastrophically when it does not. The
      !! broadcast costs one more nao^2 message and makes the ranks identical.
      class(trc_eri_t), intent(in) :: this
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(this%comm, g, op=MPI_SUM)
      call bcast(this%comm, g, n, 0)
   end subroutine sum_ranks

   !
   ! Largest |D| in each shell-pair block, for the in-kernel density screen.
   !
   ! Over a BATCH the bound is the maximum across every density, so all of them
   ! see exactly the same quartets. That is deliberate and matches mqc: a batch
   ! then gives bit-comparable results to the same densities passed one at a
   ! time, which is what makes the single-density wrapper safe. Screening each
   ! set on its own bound would be tighter and would make the batch disagree
   ! with the loop.
   !
   !
   ! J - K/2 for a density that is NOT symmetric.
   !
   ! Routes through the enumerated kernel, which reads D_mn and D_nm as
   ! distinct and writes J and K separately, so no fold is applied and the
   ! antisymmetric part of K survives. This is what the coupled-perturbed
   ! Hessian needs; mqc calls it `build_fock_direct_nosym`.
   !
   subroutine eri_fock_nosym(this, b, dmat, gmat, k_scale, j_scale, &
                             density_screen, omega)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      real(dp), intent(in), optional :: k_scale, j_scale
      real(dp), intent(in), optional :: omega
      !! Range-separation parameter, bohr^-1; see eri%fock's omega.
      logical, intent(in), optional :: density_screen
      !! Weight the Schwarz bound by the density it multiplies before
      !! screening. Default true, which is what an SCF wants: a quartet whose
      !! density elements are all negligible contributes nothing however large
      !! its integral.
      !!
      !! A COUPLED-PERTURBED SOLVE MUST PASS FALSE. Its density is a trial
      !! vector the solver drives towards zero, so a screen keyed on the
      !! density's magnitude tightens as the solve proceeds -- the operator
      !! stops being the same linear map from one matvec to the next and the
      !! iteration cannot converge. False recovers the plain Schwarz screen,
      !! which depends on the basis alone and leaves the operator fixed.

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j
      real(dp) :: jfac, kfac, ow
      type(trc_bin_launch_t) :: launch
      if (this%needs_md) error stop "terco: g functions reach eri%fock and eri%fock_resident only, for now"

      jfac = 1.0_dp
      kfac = 1.0_dp
      ow = 0.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale
      if (present(omega)) ow = omega
      if (.not. (ow >= 0.0_dp)) error stop "terco: omega < 0 or NaN (erfc/short-range) is not implemented"

      n = this%nao
      nw = this%nao_w
      allocate (jmat(1, nw, nw), kmat(1, nw, nw), dwork(1, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      if (this%dec%active) then
         call decontract_expand_host(this%dec, dmat, dp_w)
      else
         dp_w = dmat
      end if
      dwork(1, :, :) = dp_w
      jmat = 0.0_dp
      kmat = 0.0_dp
      if (screen_on(density_screen)) then
         call density_blocks(this, dwork, 1)
      else
         this%dsh = 1.0_dp
      end if

      !$acc enter data copyin(dwork, jmat, kmat)
      !$acc update device(this%dsh)
      launch = bin_launch(this, 1, jfac, kfac, ow)
      launch%nosym = .true.
      ! The non-symmetric kernel reads the primitive pairs through the bins.
      call complete_pairs(this)
      call binkernel_fock_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, dwork, jmat, kmat)
      this%nlaunch = launch%nlaunch
      this%nwork = launch%nwork
      this%nwork_masked = launch%nwork_masked
      this%nwork_variant = launch%nwork_variant
      this%nlaunch_variant = launch%nlaunch_variant
      this%parity_ops_base = launch%ops_base
      this%parity_ops_masked = launch%ops_masked
      this%parity_cost_masked = launch%cost_masked
      this%parity_cost_variant = launch%cost_variant
      !$acc wait
      !$acc update self(jmat, kmat)
      !$acc exit data delete(dwork, jmat, kmat)

      ! J and K come out in full and unfolded, so this is the definition.
      do j = 1, nw
         do i = 1, nw
            gp_w(i, j) = jfac*jmat(1, i, j) - 0.5_dp*kfac*kmat(1, i, j)
         end do
      end do
      if (this%dec%active) then
         call decontract_fold_host(this%dec, gp_w, gmat)
      else
         do j = 1, n
            do i = 1, n
               gmat(i, j) = gp_w(i, j)
            end do
         end do
      end if
      if (this%distributed) call sum_ranks(this, gmat, n*n)

      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock_nosym

   !
   ! Fock build with the density and result ALREADY ON THE DEVICE.
   !
   ! `fock` copies the density up and the result back every call. At 832
   ! functions that is 11 MB each way, and the coupled-perturbed equations run
   ! a hundred of them for matrices that never need to be on the host at all.
   !
   ! Here the caller owns the residency:
   !
   !     !$acc enter data copyin(dmat) create(gmat)
   !     do ...
   !        call eri%fock_resident(b, dmat, gmat)   ! nothing crosses
   !     end do
   !     !$acc update self(gmat)
   !     !$acc exit data delete(dmat, gmat)
   !
   ! A SEPARATE ENTRY POINT rather than a flag on `fock`. Passing a host array
   ! to this routine is a segfault, not a wrong number, and a flag makes that
   ! one typo away; two names make it a compile-time choice.
   !
   ! The scratch it needs -- three working-basis matrices, see
   ! `fock_resident_lent` -- is allocated per call here. A caller with device
   ! memory to spare between builds lends it to `fock_resident_lent` instead;
   ! the SCF lends its linear-algebra workspace.
   !
   subroutine eri_fock_resident(this, b, dmat, gmat, opts)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      !! Both must already be present on the device.
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      type(trc_fock_options_t), intent(in), optional :: opts
      !! What to build; see trc_fock_options_t and `eri_fock`. Absent:
      !! Hartree-Fock, Coulomb operator, density screen on. A COUPLED-PERTURBED
      !! SOLVE MUST SET density_screen FALSE. `count_survivors` also counts the
      !! quartets that survive screening, into this%nkept (a full extra pass
      !! over the work list, so it is opt-in); `j_engine` is read as in
      !! `fock_resident_lent`.
      type(trc_fock_options_t) :: o
      real(dp), allocatable :: work(:)
      if (present(opts)) o = opts
      allocate (work(this%fock_work_size(1)))
      !$acc enter data create(work)
      call this%fock_resident_lent(b, dmat, gmat, o, work)
      !$acc exit data delete(work)
      deallocate (work)
   end subroutine eri_fock_resident

   pure function eri_fock_work_size(this, ndens) result(nelem)
      !! Doubles of scratch `fock_resident_lent` (ndens = 1) and
      !! `fock_resident_many_lent` need: three working-basis matrices for one
      !! density, 2 ndens + 1 for several.
      class(trc_eri_t), intent(in) :: this
      integer, intent(in) :: ndens
      integer(int64) :: nelem
      nelem = int(this%nao_w, int64)**2*int(merge(3, 2*ndens + 1, ndens == 1), int64)
   end function eri_fock_work_size

   subroutine eri_fock_resident_lent(this, b, dmat, gmat, opts, work)
      !! `fock_resident` on scratch the caller lends: `work`, device-resident,
      !! at least `fock_work_size(1)` doubles, contents undefined on entry
      !! and on return. The SCF lends its linear-algebra workspace
      !! (trc_linalg's `reserve`), which is dead during a Fock build, so the
      !! build costs no device memory of its own.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      !! Both must already be present on the device.
      real(dp), intent(in)  :: dmat(this%nao, this%nao)
      real(dp), intent(out) :: gmat(this%nao, this%nao)
      type(trc_fock_options_t), intent(in) :: opts
      real(dp), intent(inout) :: work(1, this%nao_w, this%nao_w, 3)
      !! THREE working-basis matrices (the old build allocated five, plus
      !! one for MD J): the working-basis density, which is also the
      !! kernel's (1, nw, nw) operand (the same storage as an (nw, nw)
      !! matrix, so it is passed as one where that is wanted); the folded
      !! accumulator, symmetrised in place, which also stands in, unread,
      !! for the exchange accumulator only the non-symmetric kernel
      !! (eri_fock_nosym) writes; and MD J's own result. Each is 1.66 nao^2
      !! on gly20/cc-pVTZ (4252 decontracted functions against 3298).

      integer :: n, nw, i, j
      real(dp) :: jfac, kfac, ow, jbins
      type(trc_bin_launch_t) :: launch
      logical :: md_j

      jfac = opts%j_scale
      kfac = opts%k_scale
      ow = opts%omega
      if (.not. (ow >= 0.0_dp)) error stop "terco: omega < 0 or NaN (erfc/short-range) is not implemented"
      n = this%nao
      nw = this%nao_w
      ! See mdj_requested: the Coulomb part from the McMurchie-Davidson J engine, and
      ! the bins only for exchange (or not at all without it). Ordinary
      ! Coulomb operator only; an erf pass keeps the bins.
      md_j = (mdj_requested(this%mdj_force, kfac == 0.0_dp, merge(this%mdj_extended, this%has_f, &
                            opts%density_screen), opts%j_engine) .or. this%needs_md) &
             .and. jfac /= 0.0_dp .and. ow == 0.0_dp
      jbins = jfac
      if (md_j) jbins = 0.0_dp

      ! The density into the basis the integrals are formed in, on the device.
      ! Nothing crosses the bus here -- that is the whole point of `resident`.
      if (this%dec%active) then
         call decontract_expand(this%dec, dmat, work(:, :, :, 1))
      else
         do concurrent(j=1:n, i=1:n)
            work(1, i, j, 1) = dmat(i, j)
         end do
      end if

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(j=1:nw, i=1:nw)
         work(1, i, j, 2) = 0.0_dp
      end do

      ! The screening bound is computed ON THE DEVICE. Doing it host-side, as
      ! `fock` does, would mean pulling the whole density down and would undo
      ! the entire point of this routine.
      if (opts%density_screen) then
         call density_blocks_device(this%nbas, nw, this%pairs%sh_l, this%pairs%ao_off, work(:, :, :, 1), this%dsh)
      else
         ! density_screen=.false. (a coupled-perturbed solve): fill_ones_device,
         ! not a `do concurrent` reading `this%dsh`/`this%nbas` directly --
         ! trc_scf_fock's own header: never read a derived-type component
         ! inside `do concurrent` on the GPU build.
         call fill_ones_device(this%nbas, this%dsh)
      end if
      if (opts%skeleton) call skeleton_dsh(this)

      !
      ! TRC_KEPT_STATS=1 reports what fraction of the launched quartets
      ! actually survive the in-kernel screen, which is the only honest
      ! bound on what any reordering of the work could buy. It costs a
      ! SECOND pass over the whole work list, so a build timed with it on
      ! is not a build time.
      !
      block
         character(len=8) :: ke
         logical :: want
         ke = " "
         call get_environment_variable("TRC_KEPT_STATS", ke)
         want = len_trim(ke) > 0
         want = want .or. opts%count_survivors
         if (want) then
            ! jbins, not jfac: count what the bins will screen, and with the
            ! MD J engine forming J that is exchange alone.
            launch = bin_launch(this, 1, jbins, kfac, ow, opts%skeleton)
            launch%count_kept = .true.
            call binkernel_fock_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, work(:, :, :, 1), &
                                     work(:, :, :, 2), work(:, :, :, 2))
            this%nlaunch = launch%nlaunch
            this%nwork = launch%nwork
            this%nwork_masked = launch%nwork_masked
            this%nwork_variant = launch%nwork_variant
            this%nlaunch_variant = launch%nlaunch_variant
      this%parity_ops_base = launch%ops_base
      this%parity_ops_masked = launch%ops_masked
      this%parity_cost_masked = launch%cost_masked
      this%parity_cost_variant = launch%cost_variant
            this%nkept = launch%nkept
            if (len_trim(ke) > 0) then
               block
                  character(len=1024) :: logline
                  write (logline, "(a,i0,a,i0,a,f6.2,a)") &
                     "  [kept] nwork ", this%nwork, "  nkept ", this%nkept, "  survive ", &
                     100.0_dp*real(this%nkept, dp)/max(real(this%nwork, dp), 1.0_dp), "%"
                  call global_logger%info(trim(logline))
               end block
            end if
         end if
      end block
      if (.not. (md_j .and. kfac == 0.0_dp)) then
         call profiler_start("eri: bins (hgp/rys)", nvtx_only=.true.)
         launch = bin_launch(this, 1, jbins, kfac, ow, opts%skeleton)
         call binkernel_fock_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, work(:, :, :, 1), &
                                     work(:, :, :, 2), work(:, :, :, 2))
         this%nlaunch = launch%nlaunch
         this%nwork = launch%nwork
         this%nwork_masked = launch%nwork_masked
         this%nwork_variant = launch%nwork_variant
         this%nlaunch_variant = launch%nlaunch_variant
      this%parity_ops_base = launch%ops_base
      this%parity_ops_masked = launch%ops_masked
      this%parity_cost_masked = launch%cost_masked
      this%parity_cost_variant = launch%cost_variant
         !$acc wait
         call profiler_stop("eri: bins (hgp/rys)")
         call profiler_start("eri: md exchange", nvtx_only=.true.)
         call mdg_pass(this, b, work(:, :, :, 1), work(:, :, :, 2), jbins, kfac, ow, opts%skeleton)
         !$acc wait
         call profiler_stop("eri: md exchange")
      end if
      !$acc wait

      ! Fold on the device, straight into the caller's resident result.
      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      !
      ! Symmetrise in the basis the integrals were formed in, THEN fold back.
      ! The other order is a different matrix, because C is not square. In
      ! place, one thread per (i, j) pair below the diagonal; the diagonal is
      ! already 0.5 (x + x) = x.
      call symmetrise_in_place_device(nw, work(:, :, :, 2))
      if (md_j) then
         call ensure_md(this, b)
         call profiler_start("eri: md J", nvtx_only=.true.)
         call this%md%coulomb(nw, work(:, :, :, 1), work(:, :, :, 3), this%rank, this%nranks, opts%skeleton)
         call profiler_stop("eri: md J")
         ! Each rank forms its share of J (trc_mdj splits the Hermite stage's
         ! ket chunks); sum_ranks totals them with the rest of the matrix.
         call add_scaled_device(nw*nw, jfac, work(:, :, :, 3), work(:, :, :, 2))
      end if
      if (this%dec%active) then
         call decontract_fold(this%dec, work(:, :, :, 2), gmat)
      else
         do concurrent(j=1:n, i=1:n)
            gmat(i, j) = work(1, i, j, 2)
         end do
      end if

      ! Resident on every rank; the sum goes through the host for now.
      if (this%distributed) then
         !$acc update self(gmat)
         call sum_ranks(this, gmat, n*n)
         !$acc update device(gmat)
      end if
   end subroutine eri_fock_resident_lent


   subroutine eri_fock_resident_many(this, b, ndens, dmats, gmats, k_scale, j_scale, omega)
      !! `fock_resident_many_lent` on scratch of its own; see that for what
      !! is built.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ndens
      real(dp), intent(in)  :: dmats(this%nao, this%nao, ndens)
      !! Both must already be present on the device.
      real(dp), intent(out) :: gmats(this%nao, this%nao, ndens)
      real(dp), intent(in), optional :: k_scale, j_scale
      real(dp), intent(in), optional :: omega
      type(trc_fock_options_t) :: opts
      real(dp), allocatable :: work(:)
      if (present(j_scale)) opts%j_scale = j_scale
      if (present(k_scale)) opts%k_scale = k_scale
      if (present(omega)) opts%omega = omega
      allocate (work(this%fock_work_size(ndens)))
      !$acc enter data create(work)
      call this%fock_resident_many_lent(b, ndens, dmats, gmats, opts, work)
      !$acc exit data delete(work)
      deallocate (work)
   end subroutine eri_fock_resident_many

   subroutine eri_fock_resident_many_lent(this, b, ndens, dmats, gmats, opts, work)
      !! eri_fock_resident for `ndens` densities in ONE pass over the
      !! integrals: gmats(:, :, d) = j_scale J[D_d] - k_scale K[D_d]/2, every
      !! quartet evaluated once and digested into each density's own G (the
      !! bin kernels' `ndens`). An unrestricted SCF's two exchange builds
      !! (trc_scf_fock's g_uks) are the case it is for: the integrals are most
      !! of a K build, the per-density digestion the rest.
      !!
      !! The MD paths (MD J; MD exchange for the g classes) take one density
      !! at a time, so when either would run -- a basis with classes routed to
      !! MD exchange, or a Coulomb part mdj_requested would give to MD J --
      !! this is eri_fock_resident once per density instead: the same
      !! numbers, without the shared pass. eri_mdj_size_rule decides it as for one
      !! density, though it was fitted to RHF: unrestricted, MD J did better
      !! still (V100, cation doublets, cc-pVDZ, 1e-11, s/build bins / MD J:
      !! gly10 0.652 / 0.535, gly15 1.406 / 1.105, NA-040 3.122 / 2.480; and
      !! where the rule keeps the bins, water40 1.449 / 1.266, gly5 0.164 /
      !! 0.142), so the rule leaves some UHF time on the table, never loses it.
      !!
      !! Always density-screened, the bound being the maximum over the
      !! densities per shell pair (each sees at least the quartets it would
      !! alone). A coupled-perturbed solve, which must not screen on its
      !! trial densities, uses eri_fock_resident.
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ndens
      real(dp), intent(in)  :: dmats(this%nao, this%nao, ndens)
      !! Both must already be present on the device.
      real(dp), intent(out) :: gmats(this%nao, this%nao, ndens)
      type(trc_fock_options_t), intent(in) :: opts
      !! `density_screen` and `count_survivors` are not read: this build is
      !! always density-screened. `skeleton` is as for `fock_resident_lent`
      !! (every density must be invariant under the group).
      real(dp), intent(inout) :: work(this%nao_w, this%nao_w, 2*ndens + 1)
      !! At least `fock_work_size(ndens)` doubles, device-resident, as for
      !! `fock_resident_lent`: the kernel's (ndens, nw, nw) operand and
      !! accumulator (`kmat`, which only the non-symmetric kernel writes, is
      !! the accumulator again, unread), then one working-basis matrix for
      !! each density's expansion and, once the operands are packed, its
      !! symmetrised result: 2 ndens + 1 matrices, where the old build
      !! allocated 3 ndens + 2.

      integer :: n, nw, d, eng
      real(dp) :: jfac, kfac, ow
      type(trc_bin_launch_t) :: launch
      logical :: use_md

      jfac = opts%j_scale
      kfac = opts%k_scale
      ow = opts%omega
      if (.not. (ow >= 0.0_dp)) error stop "terco: omega < 0 or NaN (erfc/short-range) is not implemented"

      ! MD J by eri_mdj_size_rule, or for an f basis as before it: an f basis
      ! the rule would give the bins (a small compact one) keeps MD J here
      ! (unmeasured with two densities), and it is passed on as forced so
      ! each density's build makes the same choice.
      use_md = jfac /= 0.0_dp .and. ow == 0.0_dp .and. &
               mdj_requested(this%mdj_force, kfac == 0.0_dp, this%mdj_extended .or. this%has_f, &
                             opts%j_engine)
      if (ndens == 1 .or. this%needs_md .or. use_md) then
         eng = opts%j_engine
         if (use_md) eng = TRC_J_ENGINE_MDJ
         ! One density at a time on the same scratch: fock_work_size(1) is
         ! never more than fock_work_size(ndens).
         do d = 1, ndens
            call this%fock_resident_lent(b, dmats(:, :, d), gmats(:, :, d), &
                                         trc_fock_options_t(j_scale=jfac, k_scale=kfac, omega=ow, &
                                                            j_engine=eng, skeleton=opts%skeleton), work)
         end do
         return
      end if

      n = this%nao
      nw = this%nao_w

      ! Each density into the integrals' basis and into the kernel's
      ! (ndens, nao, nao) layout -- density index fastest, so one (mu,nu)'s
      ! ndens updates coalesce -- all on the device, with the screening bound
      ! the maximum over the densities.
      do d = 1, ndens
         if (this%dec%active) then
            call decontract_expand(this%dec, dmats(:, :, d), work(:, :, 2*ndens + 1))
         else
            call copy_square_device(n, dmats(:, :, d), work(:, :, 2*ndens + 1))
         end if
         call pack_density_device(nw, ndens, d, work(:, :, 2*ndens + 1), work(:, :, 1:ndens), &
                                  work(:, :, ndens + 1:2*ndens))
         if (d == 1) then
            call density_blocks_device(this%nbas, nw, this%pairs%sh_l, this%pairs%ao_off, &
                                       work(:, :, 2*ndens + 1), this%dsh)
         else
            call density_blocks_device_max(this%nbas, nw, this%pairs%sh_l, this%pairs%ao_off, &
                                           work(:, :, 2*ndens + 1), this%dsh)
         end if
      end do

      ! The skeleton plan needs the bound invariant under the group: the
      ! maximum over the densities is symmetrized like a single density's.
      if (opts%skeleton) call skeleton_dsh(this)

      call profiler_start("eri: bins (hgp/rys), several densities", nvtx_only=.true.)
      launch = bin_launch(this, ndens, jfac, kfac, ow, opts%skeleton)
      call binkernel_fock_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, work(:, :, 1:ndens), &
                               work(:, :, ndens + 1:2*ndens), work(:, :, ndens + 1:2*ndens))
      this%nlaunch = launch%nlaunch
      this%nwork = launch%nwork
      this%nwork_masked = launch%nwork_masked
      this%nwork_variant = launch%nwork_variant
      this%nlaunch_variant = launch%nlaunch_variant
      this%parity_ops_base = launch%ops_base
      this%parity_ops_masked = launch%ops_masked
      this%parity_cost_masked = launch%cost_masked
      this%parity_cost_variant = launch%cost_variant
      !$acc wait
      call profiler_stop("eri: bins (hgp/rys), several densities")

      ! Symmetrise each in the integrals' basis, THEN fold back (the other
      ! order is a different matrix: the decontraction C is not square).
      do d = 1, ndens
         call unpack_symmetrise_device(nw, ndens, d, work(:, :, ndens + 1:2*ndens), work(:, :, 2*ndens + 1))
         if (this%dec%active) then
            call decontract_fold(this%dec, work(:, :, 2*ndens + 1), gmats(:, :, d))
         else
            call copy_square_device(n, work(:, :, 2*ndens + 1), gmats(:, :, d))
         end if
      end do

      ! Resident on every rank; the sum goes through the host for now.
      if (this%distributed) then
         !$acc update self(gmats)
         call sum_ranks(this, gmats, n*n*ndens)
         !$acc update device(gmats)
      end if
   end subroutine eri_fock_resident_many_lent

   !
   ! Shell-block density bound, computed on the device.
   !
   ! The host version reads a density that has to be there anyway. This one
   ! exists so `fock_resident` never pulls the density down: one thread per
   ! shell pair, each scanning its own block, writing `dsh` in place while it
   ! is device-resident.
   !
   pure logical function screen_on(flag) result(res)
      !! `density_screen` with its default applied, in one place.
      !!
      !! Written as a function rather than repeated three times because the
      !! default is the DANGEROUS one for a response solve: a caller that forgets
      !! the argument gets a screen that changes with the trial vector, and a
      !! CPHF that will not converge for a reason nothing in the output points
      !! at. One place to read the default is one place to get it wrong.
      logical, intent(in), optional :: flag
      res = .true.
      if (present(flag)) res = flag
   end function screen_on

   subroutine copy_square_device(n, a, c)
      !! c := a for two resident (n, n) matrices.
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(out) :: c(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         c(i, j) = a(i, j)
      end do
   end subroutine copy_square_device

   subroutine pack_density_device(nw, ndens, d, dp_w, dwork, jmat)
      !! Density `d` into slot d of the kernel's (ndens, nw, nw) layout, and
      !! that slot's accumulator zeroed. Device to device.
      integer, intent(in) :: nw, ndens, d
      real(dp), intent(in) :: dp_w(nw, nw)
      real(dp), intent(inout) :: dwork(ndens, nw, nw), jmat(ndens, nw, nw)
      integer :: i, j
      do concurrent(j=1:nw, i=1:nw)
         dwork(d, i, j) = dp_w(i, j)
         jmat(d, i, j) = 0.0_dp
      end do
   end subroutine pack_density_device

   subroutine symmetrise_in_place_device(n, a)
      !! a := (a + a^T)/2 on a resident (n, n) matrix, in place: the
      !! iteration at (i, j), i > j, owns both (i, j) and (j, i), so no two
      !! iterations touch the same element. Bit for bit the out-of-place
      !! 0.5 (a(i,j) + a(j,i)), which is the same sum either way round.
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         if (i > j) then
            a(i, j) = 0.5_dp*(a(i, j) + a(j, i))
            a(j, i) = a(i, j)
         end if
      end do
   end subroutine symmetrise_in_place_device

   subroutine unpack_symmetrise_device(nw, ndens, d, jmat, gp_w)
      !! Slot d of the kernel's result, symmetrised.
      integer, intent(in) :: nw, ndens, d
      real(dp), intent(in) :: jmat(ndens, nw, nw)
      real(dp), intent(out) :: gp_w(nw, nw)
      integer :: i, j
      do concurrent(j=1:nw, i=1:nw)
         gp_w(i, j) = 0.5_dp*(jmat(d, i, j) + jmat(d, j, i))
      end do
   end subroutine unpack_symmetrise_device

   subroutine fill_ones_device(nshell, dsh)
      !! `dsh := 1`, as an explicit-shape dummy `do concurrent` (no derived-
      !! type component in its body) -- see eri_fock_resident's own comment
      !! on why this exists.
      integer, intent(in) :: nshell
      real(dp), intent(inout) :: dsh(nshell, nshell)
      integer :: i, j
      do concurrent(j=1:nshell, i=1:nshell)
         dsh(i, j) = 1.0_dp
      end do
   end subroutine fill_ones_device

   subroutine density_blocks_device(nshell, nao, sh_l, sh_ao, dmat, dsh)
      integer,  intent(in)    :: nshell, nao
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: dsh(nshell, nshell)

      integer :: ij, i, j, a, c, ni, nj
      real(dp) :: m

      ! Every per-iteration scalar named. This one writes the screening bound
      ! that decides which quartets are computed at all, so stating the
      ! locality rather than leaving it to be inferred is worth the line.
      do concurrent(j=1:nshell, i=1:nshell) local(ni, nj, m, a, c)
         ni = (sh_l(i) + 1)*(sh_l(i) + 2)/2
         nj = (sh_l(j) + 1)*(sh_l(j) + 2)/2
         m = 0.0_dp
         do c = 0, nj - 1
            do a = 0, ni - 1
               m = max(m, abs(dmat(sh_ao(i) + a, sh_ao(j) + c)))
            end do
         end do
         dsh(i, j) = m
      end do
   end subroutine density_blocks_device

   subroutine density_blocks_device_max(nshell, nao, sh_l, sh_ao, dmat, dsh)
      !! density_blocks_device for a further density: dsh becomes the larger
      !! of what it holds and this density's own bound, per shell pair.
      integer,  intent(in)    :: nshell, nao
      integer,  intent(in)    :: sh_l(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: dmat(nao, nao)
      real(dp), intent(inout) :: dsh(nshell, nshell)

      integer :: i, j, a, c, ni, nj
      real(dp) :: m

      do concurrent(j=1:nshell, i=1:nshell) local(ni, nj, m, a, c)
         ni = (sh_l(i) + 1)*(sh_l(i) + 2)/2
         nj = (sh_l(j) + 1)*(sh_l(j) + 2)/2
         m = dsh(i, j)
         do c = 0, nj - 1
            do a = 0, ni - 1
               m = max(m, abs(dmat(sh_ao(i) + a, sh_ao(j) + c)))
            end do
         end do
         dsh(i, j) = m
      end do
   end subroutine density_blocks_device_max

   !
   ! Shell-block density bound, over the basis the integrals are built in.
   !
   ! That is `this%pairs%sh_l`/`this%pairs%ao_off`, not the caller's basis: on a split
   ! general contraction the two differ, and bounding the wrong one indexes
   ! past the end of the density it was handed.
   !
   subroutine density_blocks(this, dmat, ndens)
      class(trc_eri_t), intent(inout) :: this
      integer, intent(in) :: ndens
      real(dp), intent(in) :: dmat(ndens, this%nao_w, this%nao_w)
      integer :: i, j, ni, nj, d
      this%dsh = 0.0_dp
      do d = 1, ndens
         do j = 1, this%nbas
            nj = (this%pairs%sh_l(j) + 1)*(this%pairs%sh_l(j) + 2)/2
            do i = 1, this%nbas
               ni = (this%pairs%sh_l(i) + 1)*(this%pairs%sh_l(i) + 2)/2
               this%dsh(i, j) = max(this%dsh(i, j), maxval(abs( &
                  dmat(d, this%pairs%ao_off(i):this%pairs%ao_off(i) + ni - 1, &
                       this%pairs%ao_off(j):this%pairs%ao_off(j) + nj - 1))))
            end do
         end do
      end do
   end subroutine density_blocks

   !
   ! Many densities, one pass over the integrals.
   !
   ! The point of the whole exercise: in a direct scheme the integral
   ! evaluation dominates and the contractions against it are nearly free, so
   ! contracting one quartet against every density in hand is the difference
   ! between one integral pass and N of them. The coupled-perturbed equations
   ! for dynamic polarizabilities need about a hundred right-hand sides.
   !
   ! The ceiling is one over the digestion fraction -- measured at 24 to 37%
   ! here, so between 2.7x and 4.2x, approached slowly.
   !
   subroutine eri_fock_many(this, b, ndens, dmats, gmats, k_scale, j_scale, &
                            density_screen, omega)
      class(trc_eri_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ndens
      real(dp), intent(in)  :: dmats(ndens, this%nao, this%nao)
      !! (ndens, nao, nao): the density index is FASTEST on purpose. With it
      !! slowest the N atomic updates for one (mu,nu) are nao^2 apart and the
      !! batch degrades past N=4; adjacent, they coalesce.
      real(dp), intent(out) :: gmats(ndens, this%nao, this%nao)
      real(dp), intent(in), optional :: k_scale, j_scale
      real(dp), intent(in), optional :: omega
      !! Range-separation parameter, bohr^-1; see eri%fock's omega.
      logical, intent(in), optional :: density_screen
      !! Weight the Schwarz bound by the density before screening (default
      !! true). A COUPLED-PERTURBED SOLVE MUST PASS FALSE -- see `eri_fock`.
      !!
      !! Over a batch the bound is the maximum across every density, so all of
      !! them see the same quartets whichever way this is set; that is what
      !! keeps a batch bit-comparable with the same densities passed one at a
      !! time.
      !!
      !! No skeleton option: this host-array entry point (the C API's
      !! trc_fock_many) always runs the full plan, also on an engine built with
      !! the point group. The SCF's several-density build is
      !! `fock_resident_many_lent`, which has it.

      real(dp), allocatable :: jmat(:, :, :), kmat(:, :, :), dwork(:, :, :)
      real(dp), allocatable :: dp_w(:, :), gp_w(:, :)
      integer :: n, nw, i, j, d
      real(dp) :: jfac, kfac, ow
      type(trc_bin_launch_t) :: launch
      if (this%needs_md) error stop "terco: g functions reach eri%fock and eri%fock_resident only, for now"

      jfac = 1.0_dp
      kfac = 1.0_dp
      ow = 0.0_dp
      if (present(j_scale)) jfac = j_scale
      if (present(k_scale)) kfac = k_scale
      if (present(omega)) ow = omega
      if (.not. (ow >= 0.0_dp)) error stop "terco: omega < 0 or NaN (erfc/short-range) is not implemented"

      n = this%nao
      nw = this%nao_w
      allocate (jmat(ndens, nw, nw), kmat(ndens, nw, nw), dwork(ndens, nw, nw))
      allocate (dp_w(nw, nw), gp_w(nw, nw))
      ! Every density transformed on the way in, every result folded on the
      ! way out. Both are per-density, and both are negligible against the
      ! one pass over the integrals they share.
      do d = 1, ndens
         if (this%dec%active) then
            call decontract_expand_host(this%dec, dmats(d, :, :), dp_w)
         else
            dp_w = dmats(d, :, :)
         end if
         dwork(d, :, :) = dp_w
      end do
      jmat = 0.0_dp
      kmat = 0.0_dp

      if (screen_on(density_screen)) then
         call density_blocks(this, dwork, ndens)
      else
         this%dsh = 1.0_dp
      end if

      !$acc enter data copyin(dwork, jmat, kmat)
      !$acc update device(this%dsh)
      launch = bin_launch(this, ndens, jfac, kfac, ow)
      call binkernel_fock_bins(this%bins, this%pairs, this%quartets, launch, this%dsh, dwork, jmat, kmat)
      this%nlaunch = launch%nlaunch
      this%nwork = launch%nwork
      this%nwork_masked = launch%nwork_masked
      this%nwork_variant = launch%nwork_variant
      this%nlaunch_variant = launch%nlaunch_variant
      this%parity_ops_base = launch%ops_base
      this%parity_ops_masked = launch%ops_masked
      this%parity_cost_masked = launch%cost_masked
      this%parity_cost_variant = launch%cost_variant
      !$acc wait
      !$acc update self(jmat)
      !$acc exit data delete(dwork, jmat, kmat)

      do d = 1, ndens
         do j = 1, nw
            do i = 1, nw
               gp_w(i, j) = 0.5_dp*(jmat(d, i, j) + jmat(d, j, i))
            end do
         end do
         if (this%dec%active) then
            call decontract_fold_host(this%dec, gp_w, gmats(d, :, :))
         else
            do j = 1, n
               do i = 1, n
                  gmats(d, i, j) = gp_w(i, j)
               end do
            end do
         end if
      end do

      if (this%distributed) call sum_ranks(this, gmats, ndens*n*n)
      deallocate (jmat, kmat, dwork, dp_w, gp_w)
   end subroutine eri_fock_many

   subroutine eri_release(this)
      class(trc_eri_t), intent(inout) :: this
      call this%mg%release()
      call this%md%release()
      if (this%on_device) then
         if (allocated(this%quartets%ao_perm)) then
            !$acc exit data delete(this%quartets%ao_perm, this%quartets%ao_perm_inv)
         end if
         !$acc exit data delete(this%pairs%sh_l, this%pairs%ao_off, this%pairs%sh_at, this%dsh, &
         !$acc                  this%bins%sp_i, this%bins%sp_j, &
         !$acc                  this%bins%sp_q, this%bins, &
         !$acc                  this%quartets%z, this%quartets%c, this%quartets%a, &
         !$acc                  this%quartets%px, this%quartets%py, this%quartets%pz, &
         !$acc                  this%quartets%ax, this%quartets%ay, this%quartets%az, &
         !$acc                  this%quartets%abx, this%quartets%aby, this%quartets%abz, &
         !$acc                  this%quartets%n, &
         !$acc                  this%pairs%off, this%pairs%n, this%pairs%p, this%pairs%r, &
         !$acc                  this%pairs%ra, this%pairs%rb, this%pairs%c)
         this%on_device = .false.
      end if
      ! Default initialisation frees every allocatable component.
      this%pairs = trc_shell_pairs_t()
      this%quartets = trc_pair_quartets_t()
      if (allocated(this%dsh))    deallocate (this%dsh)
      call decontract_release(this%dec)
      call this%pb%release()
      this%nbas = 0
      this%nao = 0
      this%nao_w = 0
   end subroutine eri_release

end module trc_eri
