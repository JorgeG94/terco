module trc_grad
   !!
   !! Two-electron gradient, over the same binned shell-pair plan the Fock
   !! build uses: one launch per angular-momentum class, Rys only -- there is
   !! no HGP gradient (out of scope).
   !!
   !! RANGE SEPARATION: `omega` switches the operator to erf(omega*r12)/r12,
   !! same meaning and same rho -> rho_theta substitution as eri_fock's own
   !! (scripts/gen_perclass.py) -- threaded only into the COMPACT pass
   !! (grad_classes/rysg_compact_dispatch), never the screen pass: the
   !! Schwarz bound stays the plain Coulomb one, which bounds the attenuated
   !! integral too (eri_fock's own comment). omega == 0 is bit-identical to
   !! the unattenuated path -- see scripts/gen_rys_grad.py's `_RANGE_SEP_BLOCK`.
   !!
   !! SCREENING: the same threshold as the energy/Fock build (`thresh`), not a
   !! separate gradient-specific tolerance. Gamma already carries the density
   !! weighting per quartet exactly as the Fock digestion's jab/jcd/kac/kad/kbc/
   !! kbd blocks do (see scripts/gen_rys_grad.py). The screen assumes, as the
   !! energy build does, that the derivative integral is of the size of the
   !! energy Schwarz product; it is not a rigorous bound on it (derivative
   !! integrals are larger), so what `thresh` buys is measured, not proved.
   !! A modest tightening factor was considered and rejected
   !! for lack of any measurement showing the plain bound insufficient for the
   !! ~1e-7 relative finite-difference tolerance -- see
   !! test/check_rys_grad_fd.F90. Measured on the 450-atom water cluster,
   !! 6-31G*, same converged 1e-13 density: the gradient at thresh = 1e-10
   !! differs from the one at 1e-13 by at most 3.9e-7 Eh/Bohr (rms 8.5e-8,
   !! against components up to 92), about 40x below a tight geometry
   !! optimisation's max-force criterion (1.5e-5).
   !!
   !! The density factor of that bound is the largest density block the
   !! quartet's Gamma multiplies: 4 D_ab, 4 D_cd for the Coulomb term, D_ac,
   !! D_ad, D_bc, D_bd for the exchange term (D_xy the shell-block maximum
   !! `dsh`, which for UHF and a pair table is the maximum over the spins or
   !! the pool). Only the blocks of the terms a LAUNCH HAS are counted: the
   !! Coulomb ones when jfac /= 0 (a pair table: when some pair of the sweep
   !! has a nonzero Coulomb weight), the exchange ones when kfac /= 0 (some
   !! nonzero exchange weight), decided once per launch in `grad_classes` and
   !! passed to the screen pass as two logicals, so no warp diverges on them.
   !! A launch with both terms counts all six blocks, as before, bit for bit;
   !! one with neither has no gradient to form and does not launch. This is
   !! sound, and changes no accuracy statement, because with jfac = 0 the
   !! kernel's Gamma has no D_ab D_cd term at all (likewise D_ac D_bd, D_ad
   !! D_bc for kfac = 0): |Gamma| is bounded by the blocks of the terms
   !! present, so a Coulomb block that is large cannot make an exchange-only
   !! quartet matter, and the bound for the terms present is the one that
   !! already held. The meaning of `thresh` for those terms is unchanged. It
   !! matters for the MD Coulomb gradient (`eri%grad_j_engine`): the Rys pass
   !! behind it is exchange-only, and the J blocks of a converged density are
   !! large where its K blocks are small.
   !!
   !! WHICH DEFAULT LAUNCHES CHANGE (no switch is involved): a launch with both
   !! terms counts all six blocks exactly as before, so the gradients of HF,
   !! global hybrids, RI-MP2 and double hybrids with exact exchange and hybrid
   !! TDDFT pair tables are arithmetically unchanged. A launch with one term
   !! screens on that term's blocks only: a pure functional's gradient
   !! (Coulomb only), the long-range pass of a range-separated hybrid
   !! (exchange only), a full-range pass whose exchange factor is zero (wB97,
   !! say), a pure-functional TDDFT/RPA pair table. Their gradients move at
   !! screening level, by quartets the density factor of the terms present
   !! already allowed to be dropped: the block maximum over those terms is
   !! the same single-block maximum, so the same Q*D <= thresh test, just
   !! without the blocks of a term the launch does not have.
   !!
   !! THE EXCHANGE-ONLY RYS PASS BEHIND THE MD COULOMB GRADIENT screens on a
   !! PRODUCT of blocks (`eri%grad_j_engine`, the default route where it
   !! exists; the launch flag `grad_k_product`, set by eri_gradient when it
   !! took the MD route). The
   !! exchange gradient's density factor is Gamma_K = -(kfac/4)(D_ac D_bd +
   !! D_ad D_bc) for RHF (the total density), so |Gamma_K| <= (|kfac|/2) P with
   !! P = max(D_ac D_bd, D_ad D_bc); for UHF it is
   !! -(kfac/2) sum_s (D^s_ac D^s_bd + D^s_ad D^s_bc), and with `dsh` the
   !! maximum over the spins |Gamma_K| <= 2 |kfac| P, four times the RHF
   !! constant (the same physical system has a UHF `dsh` about half the RHF
   !! one, so its P is about a quarter). That launch therefore drops a quartet when
   !!     Q * max(dsh(a,c) dsh(b,d), dsh(a,d) dsh(b,c)) <= thresh * m,
   !! m = `grad_k_margin_of(ndens)` = GRAD_K_PRODUCT_MARGIN for RHF and a
   !! quarter of it for UHF, so that RHF and UHF drop the same Gamma bound
   !! (halving a density is exact, so UHF with D/2 per spin keeps exactly the
   !! RHF survivors). The test does not see kfac, as the single-block one never
   !! did. Everything else is untouched and keeps the single-block test bit
   !! for bit: a launch with a Coulomb term (fused or Coulomb only), a pair
   !! table (see `grad_screen_kmargin`), and every exchange-only launch that
   !! is not behind the MD route -- notably the long-range pass of a
   !! range-separated hybrid and a K-only call from any other caller. The
   !! margin is needed because the screen multiplies the bound on |Gamma| by
   !! the ENERGY Schwarz product while the derivative integrals exceed it; see
   !! GRAD_K_PRODUCT_MARGIN for the measurements and what they do not cover.
   !!
   !! PLANNING. That pass is also planned on the exchange tiles (trc_bin_tiles,
   !! `bin_tiles_reach` with the margin): a tile pair is not launched when
   !! Q_A Q_B C^2 <= thresh*m, with C the largest density block between the
   !! two tiles' shells (all four factors of the products are such blocks),
   !! which bounds every quartet's product above, so the tiles drop only
   !! quartets the screen pass drops anyway; the survivor list and the gradient
   !! are unchanged and only the screen pass's visits (`nwork`) shrink. The
   !! fused, the single-block, the pair-table and the skeleton launches keep
   !! the whole-bin plan. TRC_EXCHANGE_TILES=0 turns it off. Measured, RHF
   !! def2-SVP, thresh 1e-12, one V100, the two-electron gradient with the MD
   !! route on, whole-bin against tiled: glycine 20-mer 1.27 -> 1.23 s, 50-mer
   !! 5.28 -> 4.67 s, 100-mer 18.6 -> 14.3 s, 40 waters 1.28 -> 1.29 s
   !! (gradients equal to 1e-12): without the tiles the screen pass kept 9% of
   !! what it visited on the 50-mer, 33% on the water cluster, whose tiles'
   !! hulls cover most of the molecule.
   !!
   !! WHICH DEFAULT LAUNCHES CHANGE with the product test: the exchange pass
   !! behind the MD Coulomb gradient, which is the default route for HF and
   !! hybrids from trc_eri's MDJ_GRAD_MIN_NAO functions (GRAD_K_PRODUCT_MARGIN
   !! carries what was measured there); the fused route below that size, and
   !! every pure functional (no exchange pass), are untouched. For the
   !! long-range pass of a range-separated hybrid, the attenuated integral is
   !! bounded by the Coulomb Schwarz quantities (Cauchy-Schwarz with the
   !! positive-definite erf kernel), which is what makes the same Schwarz
   !! screen valid there; the margin has not been measured on it (see the
   !! parameter). Skeleton gradients stay sound: the
   !! symmetrised `dsh` is invariant under the group, so each factor of the
   !! product, and hence the product, is the same for every image of a
   !! quartet (the MD route is not taken for a skeleton gradient anyway). A
   !! decontracted basis screens on the density expanded to the working
   !! basis, as the single-block test does (not measured).
   !!
   use pic_types, only: dp, int32, int64
   use trc_bins, only: pair_bins_t, trc_bin_plan_t
   use trc_bin_tiles, only: bin_tiles_requested, bin_tiles_reach, bin_tiles_plan
   use trc_eri_data, only: trc_shell_pairs_t, trc_pair_quartets_t, trc_bin_launch_t
   use trc_rys_grad_kernels, only: rysg_has_class, rysg_screen_dispatch, &
                                   rysg_compact_dispatch, rysg_pairs_dispatch, RYSG_PAIRS_NCB
   implicit none
   private
   public :: grad_bins, trc_grad_pair_t, trc_grad_pairtab_t, RYSG_PAIRS_NCB
   public :: GRAD_K_PRODUCT_MARGIN, grad_k_margin_of

   type :: trc_grad_pair_t
      !! One separable contribution to a column's gradient: the derivative, at
      !! fixed density, of  Tr(D_left (coulomb J - exchange K)[D_right]), with
      !! J and K built from the integrals as they move. A cross pair
      !! (left /= right) counts BOTH orderings, so it needs no mirror entry; a
      !! pair with left == right counts once.
      integer(int32) :: left = 0
      !! Pool index of the left density.
      integer(int32) :: right = 0
      !! Pool index of the right density.
      integer(int32) :: column = 0
      !! Output column (1..ncol) the pair adds to.
      real(dp) :: coulomb = 0.0_dp
      !! Weight of the J term; the closed-shell energy 1/2 D.J[D] is 0.5.
      real(dp) :: exchange = 0.0_dp
      !! Weight of the K term; the closed-shell energy -1/4 D.K[D] is 0.25.
   end type trc_grad_pair_t

   type :: trc_grad_pairtab_t
      !! The pair table one kernel sweep reads, in the kernels' own terms: the
      !! columns of the sweep, each owning a run of pairs, with the weights
      !! already folded (both orderings, antisymmetry). Built by
      !! eri_gradient_pairs; `pstart` has ncb + 1 entries. The default table
      !! is the one column of eri%gradient: no pairs (`pl` unallocated), the
      !! RHF/UHF contraction of the density.
      integer :: ncb = 1
      !! Gradient columns the sweep fills, 1..RYSG_PAIRS_NCB.
      integer, allocatable :: pstart(:), pl(:), pr(:)
      real(dp), allocatable :: pj(:), pk(:)
      real(dp), allocatable :: gp(:, :, :, :)
      !! The compact pass's lane-private gradient copies (3, natm, GRAD_NPRIV,
      !! ncb), allocated and released by grad_bins. Carried here, not as a
      !! component of trc_eri_t or a further argument: nvfortran 26.5 has
      !! hung or ICE'd on both.
   end type trc_grad_pairtab_t

   real(dp), parameter :: GRAD_K_PRODUCT_MARGIN = 1.0e-2_dp
   !! The margin m of the exchange-only product screen behind the MD Coulomb
   !! gradient (the module's paragraph on it): a quartet is dropped when
   !! Q*max(D_ac D_bd, D_ad D_bc) <= thresh*m, RHF. One source of truth: the
   !! host passes it to the screen pass as a launch scalar (`kmarg`) and
   !! check_grad_screen pins its value with a literal.
   !!
   !! Why not m = 1. The product bounds |Gamma| but the screen multiplies it
   !! by the ENERGY Schwarz product, and the derivative integrals are larger
   !! than that, so a bare product bound drops far too many quartets (an
   !! earlier attempt at the bare bound for the fused kernel raised the error
   !! about 100x). MEASURED, one V100, NVHPC 26.9, def2-SVP, thresh 1e-12, RHF,
   !! MD Coulomb gradient on so the Rys pass is exchange only; the two-electron
   !! gradient stage, and the max gradient error against a fused reference at
   !! thresh 1e-14:
   !!
   !!   | system | fused (the default) | single-block exchange screen | m = 1 | m = 1e-1 | m = 1e-2 | m = 1e-3 |
   !!   |---|---|---|---|---|---|---|
   !!   | gly20 time (s) | 2.23 | 1.86 | 1.02 | 1.12 | 1.23 | 1.34 |
   !!   | gly20 error | 2.9e-11 | 2.8e-11 | 3.8e-10 | 4.8e-11 | 2.7e-11 | 2.7e-11 |
   !!   | 40 waters time (s) | 1.34 | 1.75 | 0.97 | 1.11 | 1.27 | 1.45 |
   !!   | 40 waters error | 1.2e-10 | (not measured) | 5.7e-10 | 1.0e-10 | 9.9e-11 | 1.0e-10 |
   !!   | gly50 time (s) | 13.14 | 7.14 | 4.71 | 5.00 | 5.32 | 5.65 |
   !!
   !! (gly50's errors are all 1.96e-7, set by the SCF-convergence difference
   !! between the reference and the runs, not by the screen.) m = 1 costs a
   !! factor 13 in error on gly20; m = 1e-2 matched the single-block screen's
   !! error on these and was clearly cheaper than it (gly20 1.86 -> 1.23 s, 40
   !! waters 1.75 -> 1.27 s, gly50 7.14 -> 5.32 s).
   !!
   !! A fourth, small system is less favourable. Host measurement, gly3
   !! 6-31G*, converged RHF density, the exchange-only pass against Schwarz
   !! only: at thresh 1e-12 the product screen at m = 1e-2 gave an error of
   !! 8.5e-12 against 4.1e-12 for the single-block screen (2.1x worse; 1.5x
   !! worse with omega = 0.3) while keeping only 8% fewer quartets; at 1e-10
   !! 3.5e-10 against 5.4e-10; at 1e-8 5.2e-8 against 9.6e-8. Both errors at
   !! 1e-12 are below 1e-11, but "as accurate as the single-block screen" is
   !! therefore not established in general: the margin is calibrated on three
   !! systems at def2-SVP, RHF, and is within about 2x of the single-block
   !! screen's error on a fourth.
   !!
   !! The MD route is what AUTO takes for HF and hybrids from
   !! trc_eri's MDJ_GRAD_MIN_NAO functions, so the product test is on a default
   !! route, and it was measured beyond that calibration. One V100, NVHPC 26.9,
   !! thresh 1e-12, whole SCF gradients (UHF takes a quarter of the margin, by
   !! the constants in the module's paragraph): the largest component of MD
   !! route minus fused route, and each route's error against a fused run at
   !! thresh 1e-14 where the three SCFs took the same number of iterations
   !! (elsewhere the reference differs by its SCF, not by a screen):
   !!
   !!   | case | MD - fused | MD error | fused error |
   !!   |---|---|---|---|
   !!   | gly20 triplet, UHF, def2-SVP | 2.2e-11 | 2.9e-11 | 3.6e-11 |
   !!   | 40 waters cation, UHF, def2-SVP | 8.1e-11 | | |
   !!   | 40 waters cation, UPBE0, def2-SVP | 8.5e-11 | | |
   !!   | 30 waters triplet, UPBE0, def2-SVP | 2.7e-11 | 4.9e-10 | 5.0e-10 |
   !!   | 20 waters, HF, def2-TZVP | 1.4e-11 | | |
   !!   | 20 waters, PBE0, def2-TZVP | 1.1e-11 | 1.3e-11 | 1.6e-11 |
   !!   | 20 waters triplet, UHF, def2-TZVP | 3.9e-11 | 6.2e-11 | 7.3e-11 |
   !!   | gly20, HF, def2-TZVP | 2.2e-11 | 5.5e-11 | 4.2e-11 |
   !!   | 20 waters, HF, aug-cc-pVDZ | 8.6e-12 | 2.6e-11 | 2.8e-11 |
   !!   | 20 waters, PBE0, aug-cc-pVDZ | 3.5e-12 | 7.0e-12 | 7.1e-12 |
   !!   | gly20, HF, aug-cc-pVDZ | 3.1e-11 | 2.4e-10 | 2.5e-10 |
   !!   | gly20, wB97X, def2-SVP | 9.2e-12 | 1.5e-11 | 1.8e-11 |
   !!   | 30 waters, RI-MP2 gradient, cc-pVDZ | 9.8e-12 | | |
   !!
   !! The routes differ by under 1e-10 everywhere and the MD route's error is
   !! at most 1.3x the fused one (gly20 def2-TZVP). The RI-MP2 row is the
   !! gradient that forms Q(D + P) - Q(P) from two calls on different
   !! densities. NOT measured: finite-difference Hessians, thresholds other
   !! than 1e-12 on these cases, and elements beyond the second row.
   !!
   !! NOT SHIPPED, for a later decision: with the product test also on the
   !! long-range pass of a range-separated hybrid (the default route), an
   !! experimental build measured, one V100, def2-SVP, wB97X, two-electron
   !! gradient stage against main: gly20 4.54 -> 2.99 s, error against each
   !! build's own thresh-1e-14 run 1.8e-11 -> 1.6e-11; 40 waters 2.74 ->
   !! 2.11 s, 6.3e-11 -> 5.6e-11 (part of that gain is the exchange-only block
   !! selection of an earlier commit, already default).
   !!
   !! Where the margin is NOT applied: pair tables (`gradient_pairs`) keep the
   !! single-block test. Their weights pk are folded host-side and can exceed 1
   !! (both orderings of a cross pair), the bound for a column is
   !! |Gamma_K| <= 4 (sum |pk|) P with P from the pool maximum, so a fixed
   !! margin would have to scale with the largest column sum, and nothing here
   !! measures it.

   integer, parameter :: GCLASS_RADIX = 11
   !! Same fixed radix as trc_rys_select's RYS_CLASS_RADIX -- duplicated as a
   !! literal rather than imported, for the same reason: importing would
   !! pull in a module this one has no other use for, for one integer.

   integer(int64), parameter :: GRAD_COMPACT_CAP = 2_int64**25
   !! Memory bound for ORDERED SURVIVOR COMPACTION's per-batch survivor-list
   !! scratch (out_seg/out_iab/out_icd: three 4-byte integer arrays, plus an
   !! 8-byte atomic counter). A class run is walked in segment-aligned
   !! batches of at most this many work items (`grad_classes` below),
   !! and a single segment larger than this is itself walked in windows of
   !! this many, so that scratch never exceeds 3*4*GRAD_COMPACT_CAP bytes
   !! (~384 MiB at the default) however large a class run or a segment is.
   !! 2**25 is a round number comfortably inside a GPU's memory for three
   !! such arrays; no class run on w150 6-31G* is anywhere near large enough
   !! to need more than one batch -- this exists so a larger system does
   !! not exhaust device memory.

   integer, parameter :: GRAD_NPRIV = 32
   !! LANE-PRIVATE GRADIENT COPIES: how many copies gp(3, natm, GRAD_NPRIV)
   !! of the gradient the compact pass adds into (survivor i into copy
   !! mod(i-1, GRAD_NPRIV)+1; see scripts/gen_rys_grad.py's compact
   !! driver). One per lane of a warp. Measured on (ss|ss)+(ps|ss), w150
   !! 6-31G* 1e-13, V100, compact pass: 3.11 s with one grad, 1.93 s with
   !! 8 copies, 1.91 s with 32, 1.91 s with 64. 32 copies are 345 KB at
   !! 450 atoms: L2-resident, zeroed and summed once per eri%gradient
   !! (the two tiny kernels are microseconds). Whole eri%gradient on w150:
   !! 17.9 -> 16.4 s at 1e-13, 9.29 -> 8.80 s at 1e-10, with classes of
   !! L > 3 left on copy 1 (gen_rys_grad.py).

contains

   subroutine grad_bins(b, pairs, quartets, launch, dsh, dmat, grad, ptab)
      !!
      !! Walk admitted bin pairs on the host (same plan binkernel_fock_bins builds), sort
      !! by class, run `grad_classes` over the contiguous class runs. No
      !! density-weighted bin screen: Gamma is not guaranteed to admit the same
      !! bin-level bound a symmetric Fock density does (same reasoning as
      !! `eri_fock_nosym`), so every launch walks the plain Schwarz-screened plan,
      !! except the MD route's exchange-only product pass, which is planned on
      !! the exchange tiles (see PLANNING in the module header).
      !!
      !! SKELETON (launch%skeleton): the same
      !! restriction and the same per-segment weight the skeleton FOCK build uses
      !! -- only the orbit representatives of the cheaper side of each admitted
      !! bin pair, weighted by their orbit size. 1.47x (C2v water) to 3.85x
      !! (D2h pentacene) fewer work items, and the result is not a gradient until
      !! the caller projects it over the group; trc_scf_gradient's
      !! two_electron_4c does that. Only the launch COUNT is unchanged: the
      !! segments are sub-ranges of their bins, which is why the compact pass
      !! takes the bin strides sSA/sSB (see scripts/gen_rys_grad.py) and not just
      !! the segment extents sA/sNB.
      !!
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      !! The shells (ao_off, sh_at, natm) and their pairs, device-resident.
      type(trc_pair_quartets_t), intent(inout) :: quartets
      !! The bin-interleaved pairs, `a` included, device-resident. Its exchange
      !! tiles (`tiles`, host only) take this build's reach when the plan is
      !! tiled.
      type(trc_bin_launch_t), intent(inout) :: launch
      !! thresh, jfac, kfac, omega, ranks, the density's shape and
      !! `skeleton`; nwork and nkept come back. omega is
      !! the range-separation parameter, bohr^-1, as eri_fock's; zero is the
      !! ordinary Coulomb operator the screen is already built for -- see
      !! grad_classes for why the SCREEN pass takes no omega of its
      !! own.
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      type(trc_grad_pairtab_t), intent(inout) :: ptab
      !! Carries the column count of `grad`. Without pairs (`pl` unallocated,
      !! the default table): the one-column RHF/UHF contraction of `dmat`.
      !! With them: each column contracted from the same integral derivatives,
      !! `dmat` being the pool the table indexes.
      real(dp), intent(inout) :: grad(3, pairs%natm, ptab%ncb)

      type(trc_bin_plan_t) :: plan
      integer :: ga, gc, gs, gcol, natm, ncol
      real(dp) :: gsum
      logical :: cnt_j, cnt_k, tiled
      real(dp) :: kmarg

      natm = pairs%natm
      ncol = ptab%ncb
      ! The pair kernels hold RYSG_PAIRS_NCB columns of private scratch;
      ! a larger block would index past it.
      if (ncol < 1 .or. ncol > RYSG_PAIRS_NCB) error stop &
         "terco: grad_bins: a pair table must carry between 1 and RYSG_PAIRS_NCB columns"
      if (ncol > 1 .and. .not. allocated(ptab%pl)) error stop "terco: grad_bins: a multi-column sweep needs pairs"
      ! Neither term: Gamma is zero, there is nothing to contract, and no plan
      ! is built or launched, so the work counters read what ran (zero).
      call grad_terms(launch, ptab, cnt_j, cnt_k)
      if (.not. (cnt_j .or. cnt_k)) then
         launch%nwork = 0
         launch%nkept = 0
         return
      end if
      if (launch%skeleton) then
         if (.not. quartets%orb%built) error stop "terco: a skeleton gradient needs the bins' pair orbits"
         ! CANONICAL emission: the
         ! cheaper side is still the one that is cut, so the whole 1.47x-3.85x
         ! work reduction is kept, but the emitted descriptor is ordered by
         ! angular momentum so that (sLA, sLB) >= (sLC, sLD). Every class the
         ! plan can then present is inside the gradient filter's own built set
         ! by construction rather than by the bin ordering, so no new kernel
         ! and no new class run -- and so no extra launch -- comes with it.
         call plan%build_skeleton(b, quartets%orb, launch%thresh, quartets%bin_pbase, canonical=.true.)
      else
         ! The exchange-only product pass is planned on the exchange tiles
         ! (trc_bin_tiles, GRADIENT PRODUCT TEST): tile pairs whose every
         ! quartet the screen kernel's own test would reject are never
         ! launched, so the survivors and the gradient are unchanged. Only
         ! with that test (kmarg > 0, which excludes a pair table, a Coulomb
         ! term and the single-block passes) and not for a skeleton plan.
         ! `dsh` is device-resident here, as bin_tiles_reach needs.
         kmarg = grad_screen_kmargin(launch, ptab, cnt_j, cnt_k)
         tiled = kmarg > 0.0_dp .and. quartets%tiles%built .and. bin_tiles_requested()
         if (tiled) tiled = quartets%tiles%nbas == pairs%nbas
         if (tiled) then
            call bin_tiles_reach(quartets%tiles, dsh, kmarg)
            call bin_tiles_plan(quartets%tiles, plan, b, launch%thresh, quartets%bin_pbase)
            ! The gradient has kernels for canonical classes only, and the
            ! tile planner writes the walk order.
            call plan%canonicalise()
         else
            call plan%build(b, launch%thresh, quartets%bin_pbase, .false., canonical=.true.)
         end if
      end if
      launch%nwork = plan%nwork
      launch%nkept = 0
      if (plan%nwork == 0) then
         call plan%release()
         return
      end if
      call plan%sort_by_class(GCLASS_RADIX)

      ! The kernels read the segment's orbit weight themselves, so one launch
      ! still covers a whole class run. A plan with no weights -- every
      ! non-skeleton build -- is handed an array of ones rather than a branch
      ! in the kernels' per-segment prologue, exactly as the Fock launcher
      ! does it (binkernel_fock_bins). AFTER the sort, for the same reason:
      ! sort_by_class folds nint(2*sW) into its key whenever sW is allocated,
      ! and an unweighted plan has to keep the segment order it has today.
      if (.not. allocated(plan%sW)) then
         allocate (plan%sW(plan%nseg))
         plan%sW = 1.0_dp
      end if

      !$acc enter data copyin(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sPA, plan%sPB, &
      !$acc                   plan%sSA, plan%sSB, plan%sD, plan%sW, plan%sOff)
      allocate (ptab%gp(3, natm, GRAD_NPRIV, ncol))
      !$acc enter data create(ptab%gp)
      do concurrent(gcol=1:ncol, gs=1:GRAD_NPRIV, ga=1:natm, gc=1:3)
         ptab%gp(gc, ga, gs, gcol) = 0.0_dp
      end do
      if (allocated(ptab%pl)) then
         !$acc enter data copyin(ptab%pstart, ptab%pl, ptab%pr, ptab%pj, ptab%pk)
      end if
      call grad_classes(plan, b, pairs, quartets, launch, dsh, dmat, ptab)
      if (allocated(ptab%pl)) then
         !$acc exit data delete(ptab%pstart, ptab%pl, ptab%pr, ptab%pj, ptab%pk)
      end if
      ! Sum the copies into grad, on the device (grad is already there:
      ! eri_gradient's enter data).
      do concurrent(gcol=1:ncol, ga=1:natm, gc=1:3) local(gsum, gs)
         gsum = 0.0_dp
         do gs = 1, GRAD_NPRIV
            gsum = gsum + ptab%gp(gc, ga, gs, gcol)
         end do
         grad(gc, ga, gcol) = grad(gc, ga, gcol) + gsum
      end do
      !$acc exit data delete(ptab%gp)
      deallocate (ptab%gp)
      !$acc exit data delete(plan%sA, plan%sNB, plan%sOA, plan%sOB, plan%sPA, plan%sPB, &
      !$acc                  plan%sSA, plan%sSB, plan%sD, plan%sW, plan%sOff)
      call plan%release()
   end subroutine grad_bins

   subroutine grad_terms(launch, ptab, cnt_j, cnt_k)
      !! Which terms a launch has, so the screen counts the density blocks of
      !! those only (the module's SCREENING paragraph). A pair table has a
      !! term when any pair of the sweep gives it a weight: one survivor list
      !! serves every column and pair, so a term absent from all of them is
      !! the only one that may be dropped. `launch%grad_all_blocks` (a test
      !! hook) forces both.
      type(trc_bin_launch_t), intent(in) :: launch
      type(trc_grad_pairtab_t), intent(in) :: ptab
      logical, intent(out) :: cnt_j, cnt_k
      if (allocated(ptab%pl)) then
         cnt_j = any(ptab%pj /= 0.0_dp)
         cnt_k = any(ptab%pk /= 0.0_dp)
      else
         cnt_j = launch%jfac /= 0.0_dp
         cnt_k = launch%kfac /= 0.0_dp
      end if
      if (launch%grad_all_blocks) then
         cnt_j = .true.
         cnt_k = .true.
      end if
   end subroutine grad_terms

   pure function grad_k_margin_of(ndens) result(m)
      !! The margin of the exchange-only product screen for a density count:
      !! GRAD_K_PRODUCT_MARGIN for RHF, a quarter of it for UHF (the module's
      !! paragraph on the MD-route exchange pass: UHF's Gamma_K constant in P is four times
      !! RHF's, so the same Gamma bound needs a quarter of the margin).
      integer, intent(in) :: ndens
      real(dp) :: m
      m = GRAD_K_PRODUCT_MARGIN
      if (ndens == 2) m = 0.25_dp*GRAD_K_PRODUCT_MARGIN
   end function grad_k_margin_of

   function grad_screen_kmargin(launch, ptab, cnt_j, cnt_k) result(kmarg)
      !! The `kmarg` the screen pass is given: the product test's margin for
      !! the exchange-only launch behind the MD Coulomb gradient
      !! (`launch%grad_k_product`, set by eri_gradient; not inferred from
      !! jfac = 0, which the long-range pass of a range-separated hybrid and a
      !! K-only caller have too), zero (the single-block test) otherwise.
      !! Zero too for a pair table (the
      !! parameter's comment says why) and under the test hook
      !! `launch%grad_k_single`, which keeps the single-block exchange screen
      !! as the comparison the product test is judged against.
      type(trc_bin_launch_t), intent(in) :: launch
      type(trc_grad_pairtab_t), intent(in) :: ptab
      logical, intent(in) :: cnt_j, cnt_k
      real(dp) :: kmarg
      kmarg = 0.0_dp
      if (.not. launch%grad_k_product) return
      if (allocated(ptab%pl) .or. launch%grad_k_single .or. launch%grad_all_blocks) return
      if (cnt_k .and. .not. cnt_j) kmarg = grad_k_margin_of(launch%ndens)
   end function grad_screen_kmargin

   subroutine grad_classes(plan, b, pairs, quartets, launch, dsh, dmat, ptab)
      !!
      !! ORDERED SURVIVOR COMPACTION glue, one class run [c0, c1] of the
      !! class-sorted plan at a time. Screens the run into a dense (seg, iab, icd) survivor list with
      !! `rysg_screen_dispatch`, then runs `rysg_compact_dispatch` over exactly
      !! that many survivors -- see scripts/gen_rys_grad.py's
      !! screen_and_compact_text for the technique the per-class kernels
      !! implement and scripts/gen_rys_grad_dispatch.py for the two
      !! dispatchers called here.
      !!
      !! RANKS: an MPI rank takes a CONTIGUOUS slice of this class's own
      !! segment range [c0, c1] rather than an interleaved work-item stride.
      !! Every
      !! segment is assigned to exactly one rank, so every work item still
      !! runs exactly once across all ranks and the summed gradient is
      !! unchanged; load is balanced by segment count rather than by work
      !! item, and a rank's own slice stays segment-contiguous so the screen
      !! pass's RUNLEN=32 chunks keep sharing a ket within it. A tiled plan
      !! (PLANNING, module header) depends on the density, so the ranks must
      !! hold identical `dsh` for their slices to partition the work, as the
      !! Fock tile path already assumes. Its balance has
      !! not been measured -- MPI is outside today's validation
      !! (TERCO_ENABLE_MPI=OFF); nranks=1 reduces this to the whole [c0, c1]
      !! range.
      !!
      !! BATCHING: a rank's own slice is then walked forward in segment-
      !! aligned batches of at most GRAD_COMPACT_CAP work items; a segment
      !! larger than that on its own is walked in windows of GRAD_COMPACT_CAP
      !! work items, which the screen pass takes as an explicit [w0, w1]. This
      !! bounds the survivor-list scratch's memory however large a class run
      !! or a single segment is -- see GRAD_COMPACT_CAP's own comment.
      !!
      !! OMEGA reaches only `rysg_compact_dispatch` below, never
      !! `rysg_screen_dispatch`: the screen's Schwarz bound is the plain Coulomb
      !! one and stays valid at omega > 0 (pointwise larger -- see eri_fock's
      !! own comment), so the screen pass takes no omega of its own.
      !!
      type(trc_bin_plan_t), intent(in) :: plan
      !! Sorted by class, its arrays already on the device.
      type(pair_bins_t), intent(in) :: b
      type(trc_shell_pairs_t), intent(in) :: pairs
      type(trc_pair_quartets_t), intent(in) :: quartets
      type(trc_bin_launch_t), intent(inout) :: launch
      !! `intent(inout)`: the screen pass's own survivor count is summed into
      !! launch%nkept here, which is what the oracle's "the skeleton plan is
      !! smaller than the full one" assertion reads (eri%nkept). It is the
      !! atomic counter this pass already copies back per batch, so it costs
      !! nothing -- unlike the Fock path's count_kept, which is a whole extra
      !! pass and stays behind TRC_KEPT_STATS.
      real(dp), intent(in)    :: dsh(pairs%nbas, pairs%nbas)
      real(dp), intent(in)    :: dmat(launch%ndens, launch%nao, launch%nao)
      type(trc_grad_pairtab_t), intent(inout) :: ptab
      !! `gp` holds the lane-private gradient copies grad_bins sums into grad,
      !! one set per column. With pairs (`pl` allocated) the pair-table
      !! kernels run instead of the RHF/UHF ones.

      integer :: c0, c1, key, nseg, nseg_run, seg_lo, seg_hi, b0, b1
      integer, allocatable :: out_seg(:), out_iab(:), out_icd(:)
      integer :: nmax
      logical :: cnt_j, cnt_k
      real(dp) :: kmarg
      integer(int64) :: w0, w1, gend, nout
      integer(int64) :: ctr(1)
      !! A 1-element array, not a bare scalar: nvfortran's mem:separate
      !! present-table lookup for an atomic-capture target threaded through
      !! two subroutine layers hangs on a bare scalar here (same finding as
      !! the (ps|ss) prototype's own ctr).

      call grad_terms(launch, ptab, cnt_j, cnt_k)
      if (.not. (cnt_j .or. cnt_k)) return
      kmarg = grad_screen_kmargin(launch, ptab, cnt_j, cnt_k)

      nseg = plan%nseg
      c0 = 1
      do while (c0 <= nseg)
         key = plan%class_key(c0, GCLASS_RADIX)
         c1 = c0
         do while (c1 < nseg)
            if (plan%class_key(c1 + 1, GCLASS_RADIX) /= key) exit
            c1 = c1 + 1
         end do
         if (.not. rysg_has_class(key)) error stop &
            "terco: grad_bins met a class the build did not generate a gradient kernel for"
         nseg_run = c1 - c0 + 1
         seg_lo = c0 + int((int(launch%rank, int64)*int(nseg_run, int64))/int(launch%nranks, int64))
         seg_hi = c0 + int((int(launch%rank + 1, int64)*int(nseg_run, int64))/int(launch%nranks, int64)) - 1

         b0 = seg_lo
         do while (b0 <= seg_hi)
            b1 = b0
            do while (b1 < seg_hi)
               if (plan%sOff(b1 + 2) - plan%sOff(b0) > GRAD_COMPACT_CAP) exit
               b1 = b1 + 1
            end do

            ! Batches of more than one segment are at most GRAD_COMPACT_CAP
            ! work items by construction, so this inner loop runs more than
            ! once only for a single segment larger than that.
            w0 = plan%sOff(b0) + 1
            gend = plan%sOff(b1 + 1)
            do while (w0 <= gend)
               w1 = min(gend, w0 + GRAD_COMPACT_CAP - 1)
               nmax = int(w1 - w0 + 1)
               allocate (out_seg(nmax), out_iab(nmax), out_icd(nmax))
               ctr(1) = 0_int64
               !$acc enter data create(out_seg, out_iab, out_icd) copyin(ctr)
               call rysg_screen_dispatch(key, b0, b1, w0, w1, nseg, plan%sOff, plan%sA, plan%sNB, &
                                         plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sD, &
                                         b%npair, b%sp_i, b%sp_j, b%sp_q, launch%thresh, cnt_j, cnt_k, kmarg, dsh, pairs%nbas, &
                                         pairs%sh_at, out_seg, out_iab, out_icd, nmax, ctr)
               !$acc update host(ctr)
               nout = min(ctr(1), int(nmax, int64))
               launch%nkept = launch%nkept + nout  ! this rank's survivors (trc_eri_t%nkept)
               if (allocated(ptab%pl)) then
                  call rysg_pairs_dispatch(key, nout, out_seg, out_iab, out_icd, nseg, plan%sA, plan%sNB, &
                                           plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sSA, plan%sSB, &
                                           plan%sD, plan%sW, &
                                           b%npair, b%sp_i, b%sp_j, b%sp_q, &
                                           launch%thresh*1.0e-3_dp*launch%pcut_scale, launch%omega, &
                                           pairs%nbas, launch%nao, pairs%ao_off, pairs%sh_at, pairs%natm, &
                                           quartets%nqp, quartets%z, quartets%c, quartets%a, &
                                           quartets%px, quartets%py, quartets%pz, &
                                           quartets%ax, quartets%ay, quartets%az, &
                                           quartets%abx, quartets%aby, quartets%abz, quartets%n, &
                                           launch%ndens, dmat, ptab%ncb, ptab%pstart, ptab%pl, ptab%pr, &
                                           ptab%pj, ptab%pk, GRAD_NPRIV, ptab%gp)
               else
                  call rysg_compact_dispatch(key, nout, out_seg, out_iab, out_icd, nseg, plan%sA, plan%sNB, &
                                             plan%sOA, plan%sOB, plan%sPA, plan%sPB, plan%sSA, plan%sSB, &
                                             plan%sD, plan%sW, &
                                             b%npair, b%sp_i, b%sp_j, b%sp_q, launch%thresh*1.0e-3_dp*launch%pcut_scale, &
                                             launch%jfac, launch%kfac, launch%omega, pairs%nbas, launch%nao, &
                                             pairs%ao_off, pairs%sh_at, pairs%natm, &
                                             quartets%nqp, quartets%z, quartets%c, quartets%a, &
                                             quartets%px, quartets%py, quartets%pz, &
                                             quartets%ax, quartets%ay, quartets%az, &
                                             quartets%abx, quartets%aby, quartets%abz, quartets%n, &
                                             launch%ndens, dmat, GRAD_NPRIV, ptab%gp)
               end if
               !$acc exit data delete(out_seg, out_iab, out_icd, ctr)
               deallocate (out_seg, out_iab, out_icd)
               w0 = w1 + 1
            end do
            b0 = b1 + 1
         end do
         c0 = c1 + 1
      end do
   end subroutine grad_classes

end module trc_grad
