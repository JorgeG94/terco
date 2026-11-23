module trc_tddft
   !!
   !! Excitation energies of a closed-shell RHF/RKS reference, in the
   !! Tamm-Dancoff approximation (TDA; CIS for Hartree-Fock) or the full random
   !! phase approximation (RPA: TDHF for Hartree-Fock, the Casida equations for
   !! Kohn-Sham): the lowest `nroots` singlet roots by a block Davidson
   !! iteration, with the transition dipoles and oscillator strengths.
   !!
   !! THE EQUATION (singlet, restricted)
   !! ----------------------------------
   !! With X_ia the amplitude of the i -> a excitation (i occupied, a virtual)
   !! and T = C_occ X C_vir^T its AO transition density (NON-symmetric):
   !!
   !!     [A X]_ia = (e_a - e_i) X_ia
   !!                + [C_occ^T (2 J[T] - c_x K[T] + 2 f_xc[T]) C_vir]_ia
   !!
   !! J[T]_mn = (mn|ls) T_ls and K[T]_mn = (ml|ns) T_ls, the pyscf get_jk(hermi=0)
   !! convention; c_x is the functional's exact-exchange fraction (`func%exx`,
   !! 1 for Hartree-Fock, where f_xc is absent). The exchange of a
   !! non-symmetric density is not the exchange of its symmetric part, which is
   !! why this uses `eri%fock_split_many_lent` (the folded kernels on the
   !! symmetric and antisymmetric parts of T, the device-resident form of
   !! `fock_split_many`; `fock_nosym_many` is its reference): it
   !! returns j J - k K/2, so the call passes j_scale = 2 and k_scale = 2 c_x. f_xc[T] is
   !! `trc_xc_fxc_rks_apply_many_lent` (it symmetrises T itself, the only part the
   !! density response sees), against a cache built once from the ground-state
   !! density. A is real symmetric and, for a stable ground state, positive
   !! definite; the eigenvalues are the excitation energies in hartree.
   !!
   !! THE RPA EQUATION
   !! ----------------
   !! With B_ia,jb = 2 (ia|jb) - c_x (ib|ja) + 2 f_xc the singlet RPA problem is
   !!
   !!     [  A   B ] [X]       [X]
   !!     [ -B  -A ] [Y]  = w  [Y]
   !!
   !! and, in p = X + Y and m = X - Y,  (A-B) m = w p,  (A+B) p = w m, so
   !! w^2 is an eigenvalue of (A-B)(A+B). The two products are
   !!
   !!     (A+B) p = (e_a - e_i) p + P[ 4 J[Ts] - 2 c_x K[Ts] + 4 f_xc[Ts] ]
   !!     (A-B) m = (e_a - e_i) m + P[           - 2 c_x K[Ta]          ]
   !!
   !! with P[G] = C_occ^T G C_vir, T = C_occ x C_vir^T the transition density of
   !! a trial vector x, Ts = (T + T^T)/2 its symmetric and Ta = (T - T^T)/2 its
   !! ANTISYMMETRIC part. (K[T] + K[T^T] = 2 K[Ts] and K[T] - K[T^T] = 2 K[Ta];
   !! J[T^T] = J[T], so J sees only Ts, and f_xc, a functional of the density,
   !! only the symmetric part.) A-B therefore has no J and no f_xc term: for a
   !! pure functional it is the diagonal (e_a - e_i), and the antisymmetric
   !! half of the batch is not built. Ts is exactly symmetric and Ta exactly
   !! antisymmetric (`tddft_split_density`), the two cases the folded kernels
   !! are exact for, so the symmetric batch is one `fock_folded_many_lent`
   !! (j_scale = 4, k_scale = 4 c_x) and the antisymmetric one a second, K
   !! only (j_scale = 0, k_scale = 4 c_x, `antisymmetric`), which keeps the
   !! antisymmetric part of K that a symmetric fold would drop. The enumerated
   !! `fock_nosym_many_lent` they replace costs several times as much per
   !! density. Measured on one V100 (NVHPC 26.9; 40 waters, def2-SVP
   !! Cartesian, PBE0, RPA, 5 roots, conv_tol 1e-5, missed_root_iterations 0;
   !! total wall with the SCF): 2306 s through the enumerated kernel, 159 s
   !! through the two folded batches, the five roots equal to the eight
   !! printed decimals (Eh). f_xc[Ts] is `fxc_apply_many` as before.
   !!
   !! The Davidson iteration holds ONE orthonormal basis V for both p and m
   !! (Stratmann, Scuseria, Frisch, J. Chem. Phys. 109, 8218 (1998)) and the
   !! products (A+B)V and (A-B)V. Its subspace problem is Am v = w u,
   !! Ap u = w v, i.e. Am Ap u = w^2 u with the small PD matrices Ap = V^T
   !! (A+B) V and Am = V^T (A-B) V: the symmetric Am^{1/2} Ap Am^{1/2} z = w^2 z
   !! gives w and, with u = Am^{1/2} z/sqrt(w), v = Ap u/w, the Ritz vectors
   !! p = V u, m = V v with p.m = 1. The residuals are r1 = (A+B)p - w m and
   !! r2 = (A-B)m - w p; the convergence measure is sqrt((|r1|^2 + |r2|^2)/2),
   !! the norm of the residual of the (X, Y) equation above. The corrections
   !! are the diagonally preconditioned residuals of X (divided by e_a - e_i - w)
   !! and of Y (divided by e_a - e_i + w), both added to V. Each iteration
   !! sends the new trial vectors through one folded symmetric batch (their
   !! Ts) and, when c_x /= 0, one folded antisymmetric batch (their Ta), and
   !! ONE `fxc_apply_many` call.
   !! A subspace Am or Ap that is not positive definite, or a non-positive
   !! root w^2, is an error: the reference may be unstable (or the roots
   !! complex). LIMITATION: this is tested on the Davidson subspace only,
   !! V^T (A-B) V and V^T (A+B) V, never on the full A-B and A+B. A matrix that
   !! is indefinite only along directions the subspace has not reached goes
   !! unnoticed, so a returned result is not a stability certificate of the
   !! reference; the converged roots themselves are valid for the RPA problem
   !! the subspace did resolve.
   !!
   !! CONVENTION. X and Y are normalised X.X - Y.Y = 1 (the TDA's Y = 0, X.X = 1
   !! included). pyscf's tdscf (X, Y) are normalised 2 (X.X - Y.Y) = 1, so
   !! X(terco) = sqrt(2) X(pyscf) and Y(terco) = sqrt(2) Y(pyscf); the sign of
   !! each root is arbitrary, a degenerate multiplet is an arbitrary
   !! rotation. The transition dipole is sqrt(2) sum_ia (X + Y)_ia <i|r|a> (the
   !! length gauge, pyscf's `transition_dipole`) and f = (2/3) w |mu|^2.
   !!
   !! THE SOLVER (TDA)
   !! ----------------
   !! Every new trial vector of an iteration goes through ONE `fock_split_many_lent`
   !! call (two folded integral passes: symmetric and antisymmetric parts) and
   !! ONE `fxc_apply_many` call, whatever the number of roots. The RPA matvec
   !! below uses `fock_folded_many_lent` (one batch per symmetry). Both route the batch through the
   !! kernels in chunks of a few densities (`eri_lent_chunk_size`), and both
   !! are density-screened by default (`trc_tddft_options_t%density_screen`).
   !! A density-keyed screen on a trial vector that shrinks would change the
   !! operator between matvecs (see `eri%fock_nosym_many`), which is why the
   !! coupled-perturbed solves keep it off; here it is safe because the Davidson
   !! trial vectors do not shrink: every one is a correction orthogonalised
   !! against the subspace and NORMALISED to unit length
   !! (`tddft_orthonormalise`; the guesses are unit vectors, the Ritz restarts
   !! orthonormal), so the transition densities are O(1) throughout and the
   !! screen's bound does not tighten as the solve converges. The residual
   !! tolerance, not the density, drives convergence.
   !!
   !! Guesses are the unit vectors of the lowest orbital-energy differences (a
   !! tie at the cut is kept whole). Each
   !! iteration diagonalises the subspace matrix, forms the residuals of the
   !! lowest `nroots` Ritz pairs and, for every unconverged one, the
   !! diagonally preconditioned correction -R/(d - w). The corrections are
   !! orthonormalised against the subspace and each other (classical
   !! Gram-Schmidt twice) and a correction that is numerically inside the span
   !! is dropped. When the subspace would exceed `max_space` it collapses to
   !! the current Ritz vectors. All the roots iterate together as a block, so a
   !! degenerate multiplet is converged as one: nothing here assumes the roots
   !! are non-degenerate (the amplitudes within a degenerate multiplet are an
   !! arbitrary orthonormal rotation of each other).
   !!
   !! MISSED ROOTS. Davidson follows the lowest `nroots` RITZ values, and a
   !! preconditioned correction of a Ritz vector stays in the symmetry of
   !! that vector. A root of a symmetry the guesses represent poorly has a
   !! Ritz value far above its eigenvalue, so with few guesses it can sit
   !! above the tracked window and never be iterated. Measured twice:
   !! formamide HF (Cs, A'/A''), nroots = 2 with 2*nroots guesses returned the
   !! 1st and 3rd roots and skipped the 2nd (0.3464 Eh); formaldehyde 6-31+G
   !! HF (C2v), nroots = 3, 4 and 5 all skipped the root at 0.348172 Eh. That
   !! root is made of the HOMO-2 -> LUMO and HOMO-2 -> LUMO+4 pairs (94% of its
   !! weight), ranks 9 and 13 of the orbital-energy differences, in a symmetry
   !! sector in which neither the 7 guesses nor the first probe's unit vectors
   !! (the 7 next differences, ranks 8 to 14, which do hold those two pairs)
   !! could help: the probe put the sector's Ritz value at 0.46 Eh, above the
   !! 5 Ritz vectors tracked (the roots at 0.157, 0.328, 0.353 and the buffer's
   !! 0.362, 0.371), so no correction was ever made in it, and the returned
   !! roots passed every residual test. Four measures: (1) a BUFFER of
   !! max(2, nroots/2) extra Ritz vectors is iterated along with the wanted ones,
   !! to a looser residual (BUFFER_TOL); (2) generous guesses, a tie at the cut
   !! kept whole; (3) one PROBE before convergence is accepted: the unit vectors
   !! of the next-lowest differences and two random vectors, which have weight in
   !! every sector (a root whose lowest pairs rank beyond the probe is reached
   !! only through them); (4) after the probe, `tddft_window` also iterates
   !! CANDIDATE Ritz pairs beyond the buffer: those whose interval
   !! [theta - |r|, theta + |r|] reaches below the nroots-th Ritz value (some
   !! eigenvalue lies within |r| of theta, so it may be a root that belongs in
   !! the window), that are not converged (|r| >= BUFFER_TOL) and have
   !! |r| < theta/2; the candidates are iterated, the lowest first.
   !! At most max(2, nroots) candidates are iterated at a time (the Ritz vectors
   !! kept are nt + that many, nt = nroots + max(2, nroots/2), further limited
   !! by max_space/2, or max_space/4 for the RPA: nroots = 1: 3 + 2 = 5, nroots
   !! = 3: 5 + 3 = 8, nroots = 6: 9 + 6 = 15, for the default max_space), and
   !! for `missed_root_iterations` iterations after the probe. Measure (4) is OPT-IN: the
   !! default is 0, so the probe runs and the candidates are reported open
   !! (`window_open`, the warning) but not iterated. The trade-off, measured: 40
   !! waters, def2-SVP Cartesian, PBE0, TDA, 5 roots, conv_tol 1e-5 (V100, total
   !! wall) take 233 s at 24 and 108 s at 0 with the same five roots to eight
   !! decimals (candidates left open at both), while formaldehyde 6-31+G
   !! Hartree-Fock at 0 returns the third root at 0.353367 Eh (TDA) and
   !! 0.346886 Eh (RPA) instead of 0.348172 and 0.338715, and reports one
   !! candidate open; setting the option to 24 recovers the root. Without (4)
   !! the probe only RECORDED the sector. In the formaldehyde case the sector's pair then
   !! descends from 0.464 to 0.3515 in one correction and converges to
   !! 0.348172. A pair with |r| >= theta/2 is neither iterated nor counted
   !! (see `tddft_window`). Cost, measured on the 51 Hartree-Fock, B3LYP and PBE
   !! cases of the basis-set sweep (water, formaldehyde, STO-3G to aug-cc-pVDZ,
   !! nroots = 3) against the solver without (3), (4) and the larger default
   !! max_space (24*max(nroots, 3) instead of 12*nroots): TDA matvecs 1842 ->
   !! 1928 in total (per case -3 to +4, median +2: the two random probe
   !! vectors), RPA 3054 -> 2939 (-22 to +6, median 0: the larger subspace
   !! collapses less); formaldehyde/6-31+G/HF, where a root is found, 56 -> 67
   !! (TDA) and 106 -> 124 (RPA).
   !!
   !! WHAT THE CHECK DOES NOT COVER. The check is a heuristic, not a proof: no
   !! Davidson variant guarantees the lowest roots without symmetry
   !! information. `res%window_open` / `res%nwindow_open` say when it ended
   !! with candidates open: candidates not iterated (the cap full or the
   !! budget used up), still unconverged when the budget ended, or abandoned
   !! because their correction lay in the space. Ritz pairs with |r| >= theta/2
   !! (random probe vectors' pairs in the high spectrum) are not candidates: a
   !! root only they could reach is NOT reported. `converged` stays
   !! true then (the returned roots are eigenpairs and the lowest found), but
   !! one of the roots may be missing, and the solver also logs a warning; the
   !! flag, not the log, is what a caller must test. The tests compare the
   !! lowest roots with the dense eigenvalues of A (A-B, A+B), not only with
   !! pyscf's.
   !!
   !! A breakdown (no usable new direction while a root is unconverged), a NaN
   !! or Inf anywhere, a non-positive excitation energy or running out of
   !! iterations is an error (`res%err`, `res%converged = .false.`), never a
   !! silently returned approximate answer. A NaN residual reads as not
   !! converged (`.not. (r < tol)`), as in trc_cpks.
   !!
   !! WHAT THIS DOES NOT DO
   !! ---------------------
   !! Triplets, UHF/UKS references, meta-GGA, VV10, range-separated and
   !! double hybrids: `trc_tddft_refuse_reason` names
   !! each refusal, and the solver reports it as an error before any work.
   !! A basis with a shell above MULT_LMAX (f functions: the multipole integrals
   !! of the transition dipole stop at d at every TERCO_LMAX) is refused the same
   !! way, up front and before the SCF-sized work, with a message that names the
   !! limit. The energies alone could be computed, but a TDDFT result without its
   !! oscillator strengths (NaN or a flag the caller must remember to test) is a
   !! surprise in a spectroscopy code, and a refusal is the one outcome nobody
   !! can miss; a caller that wants energies of an f basis has to ask for that
   !! explicitly once the dipoles reach f.
   !!
   use pic_types, only: dp, int64, default_int
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   use pic_sorting, only: sort_index
   use pic_rng, only: splitmix64_t, next_real_dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_multipoles, TRC_NMULT
   use trc_mult_kernels, only: MULT_LMAX
   use trc_eri, only: trc_eri_t, trc_fock_options_t
   use trc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, XC_FAMILY_LDA, XC_FAMILY_GGA
   use trc_xc, only: trc_xc_fxc_cache_t, trc_xc_fxc_rks_apply_many_lent
   implicit none
   private

   public :: trc_tddft_options_t, trc_tddft_result_t, trc_tddft_op_t, trc_tddft_reference_t
   public :: trc_tddft_solve_rks, trc_tddft_refuse_reason
   public :: tddft_window, tddft_window_t, tddft_space_t
   !! Public for the unit test of the window rule (a synthetic Ritz set).

   real(dp), parameter :: OSC_PREFACTOR = 2.0_dp/3.0_dp
   !! The 2/3 of the oscillator strength f = (2/3) w |mu|^2, atomic units.
   real(dp), parameter :: BUFFER_TOL = 1.0e-3_dp
   !! Residual the buffer Ritz vectors (those past `nroots`) are iterated to; at
   !! least 100 conv_tol. Their Ritz values are then good to ~1e-6.
   real(dp), parameter :: TIE_TOL = 1.0e-4_dp
   !! Orbital-energy differences closer than this (hartree) count as tied
   !! when the guesses are cut.
   real(dp), parameter :: PRECOND_FLOOR = 1.0e-4_dp
   !! The preconditioner's denominator (d - w) is kept at least this far from
   !! zero, sign preserved.
   integer, parameter :: NRAND_PROBE = 2
   !! Random vectors in the probe, beside the unit vectors of the next-lowest
   !! orbital-energy differences.
   character(len=*), parameter :: SCREEN_OFF_MSG = &
      "trc_tddft: the density screen of the excited-state matvec is OFF (density_screen = false): every matvec "// &
      "evaluates all the integrals above the Schwarz threshold, so the run is slower than the default"
   !! Logged once per solve at warning level when `density_screen` is false: a parameter, so the message is
   !! one exact-length constant (NVHPC 26.9 mis-sized a concatenation built at run time).
   integer(int64), parameter :: PROBE_SEED = 20241118_int64
   !! Seed of the probe's random vectors (a date, nothing tuned: any value works).
   real(dp), parameter :: DROP_REL = 1.0e-6_dp
   !! A correction whose norm falls below this fraction of its norm before
   !! orthogonalisation lies inside the subspace and is dropped.

   type :: trc_tddft_options_t
      integer :: nroots = 3
      !! Lowest singlet roots wanted.
      logical :: tda = .true.
      !! True: the Tamm-Dancoff approximation (CIS for Hartree-Fock). False:
      !! the full RPA (TDHF / Casida), see the module header.
      logical :: singlet = .true.
      !! Only singlets are implemented; false is refused.
      real(dp) :: conv_tol = 1.0e-6_dp
      !! Converged when ||A X - w X|| (Frobenius, X normalised) is below this
      !! for every root; for the RPA the residual of the (X, Y) equation,
      !! sqrt((|r1|^2 + |r2|^2)/2) with X.X - Y.Y = 1 (module header). The
      !! excitation energy error is of order conv_tol^2/gap, the amplitude
      !! error of order conv_tol/gap.
      integer :: max_space = 0
      !! Largest subspace before the collapse; 0 takes 24*max(nroots, 3), that is
      !! 72 vectors for up to three roots (arrays of nocc*nvir x 72, ten of them
      !! for the RPA). Below three roots 12*nroots left no room for the Ritz pairs
      !! `tddft_window` iterates beside the wanted ones (nroots = 1 on
      !! HNC/6-31+G/HF: the RPA restarted every second iteration and did not
      !! converge in 50), and at 36 the extra pairs of formaldehyde/6-31+G/HF
      !! kept being collapsed away (RPA: 202 matvecs against 124 at 72). Raised to
      !! 3*(nroots + buffer) if smaller (4*(nroots + buffer) for the RPA, whose
      !! subspace holds p and m), and capped at the dimension of the problem.
      integer :: max_iteration = 50
      !! Most Davidson iterations (matvec passes).
      integer :: missed_root_iterations = 0
      !! Iterations after the probe in which Ritz pairs beyond the buffer are still
      !! iterated (`tddft_window`, module header MISSED ROOTS). The default 0 iterates
      !! no candidate after the probe: they are reported open (`window_open`,
      !! `nwindow_open` and the warning, which names this option), so the run stays
      !! honest. A positive value (24 was the previous default) iterates them and
      !! trades matvecs against the missed-root guarantee: the sector the probe reaches
      !! descends into the window in the first few iterations, so it is a safety bound,
      !! not a typical count. Measured on 40 waters, def2-SVP Cartesian, PBE0, TDA,
      !! 5 roots, conv_tol 1e-5 (V100, total wall): 233 s at 24 against 108 s at 0,
      !! the same five roots to eight decimals, candidates left open at both.
      !! Formaldehyde 6-31+G Hartree-Fock needs it: at 0 the third root is 0.353367 Eh
      !! (TDA) and 0.346886 Eh (RPA) instead of 0.348172 and 0.338715, with one
      !! candidate reported open. A negative value is refused.
      logical :: density_screen = .true.
      !! Screen the matvec's two-electron integrals by the trial densities
      !! themselves, as the SCF Fock build does. Per launch of the kernels
      !! (a chunk of a few trial densities, `eri_lent_chunk_size`) the bound
      !! is the largest density block over the densities in the launch, per
      !! shell pair (for the exchange-only antisymmetric pass of the TDA, the
      !! antisymmetric part's); a quartet is dropped when its Schwarz bound
      !! times that block (weighted by the kernel's J and K factors) is below
      !! the integral threshold, the criterion of the SCF Fock build, so each
      !! dropped quartet contributes below it.
      !! The trial vectors are orthonormal (see the module header), their
      !! densities O(1), so the bound does not tighten as the solve converges:
      !! the operator is the same from one matvec to the next up to the
      !! threshold. The coupled-perturbed solves (CPHF, CPKS, the Z-vector of
      !! the excited-state gradient) do NOT use this: their vectors shrink to
      !! zero and they always run with the screen off.
      !! Measured on one V100 (NVHPC 26.9; def2-SVP Cartesian, TDA, 5 roots,
      !! conv_tol 1e-5, total wall including the SCF; the roots are identical to
      !! all eight printed decimals in Eh with the screen on and off), shipped
      !! device route. 40 waters (1000 functions), PBE0, missed_root_iterations
      !! 0: 193 s before the chunked matvec, 146 s chunked and unscreened, 108 s
      !! chunked and screened (the default). 40 waters, PBE0, missed_root_iterations 24:
      !! 384 s before, 233 s now; PBE, missed_root_iterations 24: 43 s before, 32 s now.
      !! gly20 (1525 functions), PBE0, budget 0: 162 s before, 88 s now. Water
      !! hexamer (150 functions), PBE0, missed_root_iterations 24: 3.14 s before, 3.16 s
      !! screened and 3.49 s unscreened now (the whole batch is one launch at
      !! this size, and the unscreened path is not faster on small systems
      !! either). False runs the unscreened matvec, which is slower, and the run
      !! says so in the log.
      integer :: nguess = 0
      !! Unit-vector guesses; 0 takes max(2*nroots, nroots + 4), capped at the
      !! dimension and extended over a tie at the cut. A probe with the
      !! next-lowest differences runs before convergence is accepted.
      logical :: verbose = .false.
   end type trc_tddft_options_t

   type :: trc_tddft_result_t
      logical :: converged = .false.
      type(error_t) :: err
      !! The refusal, breakdown, NaN or non-convergence that stopped the solve;
      !! no error when `converged`.
      integer :: nroots = 0
      integer :: iterations = 0
      !! Davidson iterations taken.
      integer :: nmatvec = 0
      !! Trial vectors multiplied by A in total.
      integer :: nscreened_apply = 0
      !! Matvec passes (calls of `apply`/`apply_rpa`, one per Davidson iteration that had new
      !! vectors) that ran with the density screen on: the operator's own count, read back when
      !! the solve ends, so it shows what the matvec did and not what the options asked for.
      !! `iterations` or fewer with the default, 0 with `density_screen = .false.`.
      integer :: space = 0
      !! Subspace dimension at the end.
      integer :: ncollapse = 0
      !! Times the subspace collapsed onto the Ritz vectors.
      logical :: window_open = .false.
      !! True when the check for a missed root (header, MISSED ROOTS) ended with
      !! candidate Ritz pairs left open: `converged` is then still true and the
      !! returned roots are the ones found, but a root below the highest of them,
      !! in a symmetry the guesses do not reach, is not excluded. False on every
      !! other path (a refusal and a failure included).
      integer :: nwindow_open = 0
      !! How many candidates were left open.
      real(dp), allocatable :: omega(:)
      !! (nroots) excitation energies, hartree, ascending.
      real(dp), allocatable :: xamp(:, :, :)
      !! (nocc, nvir, nroots) amplitudes X_ia: orthonormal for the TDA, and
      !! for the RPA normalised X.X - Y.Y = 1 (module header).
      real(dp), allocatable :: yamp(:, :, :)
      !! (nocc, nvir, nroots) de-excitation amplitudes Y_ia; zero for the TDA.
      real(dp), allocatable :: tdip(:, :)
      !! (3, nroots) transition dipoles <0|r|n>, atomic units
      !! (sqrt(2) sum_ia (X + Y)_ia <i|r|a>).
      real(dp), allocatable :: osc(:)
      !! (nroots) oscillator strengths (2/3) w |mu|^2.
      real(dp), allocatable :: resid(:)
      !! (nroots) the residual norms at the end (||A X - w X||; for the RPA
      !! the (X, Y) residual of `conv_tol`), from the solver's own
      !! bookkeeping; a caller who wants an independent figure multiplies by
      !! `trc_tddft_op_t%apply` / `apply_rpa`.
   end type trc_tddft_result_t

   type :: trc_tddft_reference_t
      !! The closed-shell SCF reference a solve starts from: the orbitals,
      !! their energies and the total density (copied, host), the number of
      !! occupied orbitals, and for Kohn-Sham the functional and an optional XC
      !! scratch budget. Filled by `set`. The grid the SCF converged on is not
      !! held here (it is device-mapped): the solve and the operator take it as
      !! an argument.
      integer :: nocc = 0
      real(dp), allocatable :: cmo(:, :), eps(:), dmat(:, :)
      !! (nao, nmo) orbitals, occupied first; (nmo); (nao, nao) total density.
      logical :: dft = .false.
      type(trc_xc_functional_t) :: func
      logical :: has_budget = .false.
      integer(int64) :: budget = 0
      !! trc_xc_fxc_rks_apply_many_lent's cap in doubles; unused unless has_budget.
   contains
      procedure :: set => tddft_reference_set
   end type trc_tddft_reference_t

   type :: tddft_window_t
      !! What `tddft_window` decides: act marks the iterated Ritz pairs,
      !! kept_index(1:nact) lists them in order, nopen counts the candidates the
      !! check leaves open. The caller allocates act and kept_index (one slot per
      !! Ritz pair).
      logical, allocatable :: act(:)
      integer, allocatable :: kept_index(:)
      integer :: nact = 0
      integer :: nopen = 0
   end type tddft_window_t

   type :: tddft_space_t
      !! The Davidson subspace: its orthonormal basis vs and A times it, avs
      !! (nov, mmax) with m columns in use, and the guess bookkeeping the probe
      !! needs: the ordering `ord` of the orbital-energy differences, the
      !! nguess guessed, nprobe unit plus nrand random probe vectors tp; dgn is
      !! the diagonal e_a - e_i in the amplitudes' flat order.
      integer :: m = 0, mmax = 0, nguess = 0, nprobe = 0
      integer :: nrand = 0
      !! Random vectors that follow the nprobe unit vectors of the probe.
      integer :: nt = 0
      !! Ritz vectors iterated: the wanted roots plus the buffer.
      integer :: nkcap = 0
      !! Most Ritz vectors iterated once the probe has run (the lowest nt always,
      !! plus those whose Ritz interval reaches the wanted window, see
      !! `tddft_window`).
      integer :: ncollapse = 0
      !! Collapses onto the Ritz vectors so far.
      real(dp), allocatable :: vs(:, :), avs(:, :), tp(:, :), dgn(:)
      real(dp), allocatable :: ams(:, :)
      !! (nov, mmax) (A-B) times vs, the RPA's second product; avs then holds (A+B) vs.
      integer, allocatable :: ord(:)
   end type tddft_space_t

   type :: trc_tddft_op_t
      !! The matvec X -> A X of the module header, for any number of vectors
      !! in one integral pass. Holds the orbitals and the XC cache; the
      !! orbitals are device-resident for the life of the object (mapped on
      !! their own, never the type), the trial densities per call.
      integer :: nao = 0, nocc = 0, nvir = 0
      logical :: ready = .false.
      logical :: dft = .false.
      real(dp) :: cx = 1.0_dp
      !! Exact-exchange fraction applied to K.
      logical :: density_screen = .false.
      !! Screen the matvec's integrals by the trial densities
      !! (`trc_tddft_options_t%density_screen`, which the solve copies here
      !! after `setup`). OFF by default, unlike the option: the screen is sound
      !! only for trial vectors kept at unit norm, which the solve guarantees
      !! and a caller that builds and applies the operator itself does not, so
      !! that caller opts in.
      logical :: count_survivors = .false.
      !! Test hook: `apply` and `apply_rpa` pass `count_survivors` to the folded builds, so `eri%nkept`
      !! reads back the quartets that survived both screens in the last call (for `apply_rpa`, summed over
      !! its symmetric and antisymmetric batches).
      integer :: napply_screened = 0
      !! Calls of `apply`/`apply_rpa` since `setup` that ran with `density_screen` true.
      logical :: force_antisym = .false.
      !! `apply_rpa` carries the antisymmetric half of its batch even when
      !! c_x = 0, where A-B is the diagonal: Ta is formed, its Fock results are
      !! zeroed on the device without a launch (K is not wanted, and the folded
      !! kernels' Coulomb half is not valid for an antisymmetric density) and
      !! the A-B transforms run on them. A test hook for that code path; it no
      !! longer measures that J[Ta] vanishes (`check_tddft_api`'s `rpa_route`
      !! does, against the enumerated kernel).
      logical :: has_budget = .false.
      integer(int64) :: budget = 0
      !! Passed to the XC response only when `has_budget`; absent it takes the
      !! routine's own default, which keeps all vectors in one pass.
      real(dp), allocatable :: cocc(:, :), cvir(:, :), de(:, :)
      !! (nao, nocc), (nao, nvir), (nocc, nvir) orbital-energy differences.
      type(trc_xc_functional_t) :: func
      type(trc_xc_fxc_cache_t) :: cache
   contains
      procedure :: setup => tddft_op_setup
      procedure :: apply => tddft_op_apply
      procedure :: apply_rpa => tddft_op_apply_rpa
      procedure :: transition_dipoles => tddft_op_dipoles
      procedure :: release => tddft_op_release
   end type trc_tddft_op_t

contains

   pure function trc_tddft_refuse_reason(opts, func, nspin, lmax) result(msg)
      !!
      !! Why this calculation cannot run through trc_tddft, or an empty string if
      !! it can. `func` absent is Hartree-Fock (CIS); `nspin` is the reference's
      !! spin count (1 restricted; default 1); `lmax` is the highest angular
      !! momentum of a shell of the basis (`maxval(b%sh_l)`), refused above the
      !! dipole integrals' ceiling MULT_LMAX.
      !!
      type(trc_tddft_options_t), intent(in) :: opts
      type(trc_xc_functional_t), intent(in), optional :: func
      integer, intent(in), optional :: nspin, lmax
      character(len=256) :: msg
      msg = ""
      if (present(lmax)) then
         if (lmax > MULT_LMAX) then
            msg = "trc_tddft: the basis has l = "//to_string(lmax)//" shells and the transition dipole needs "// &
                  "multipole integrals, which stop at l = "//to_string(MULT_LMAX)//" at every TERCO_LMAX (MULT_LMAX)"
            return
         end if
      end if
      if (present(nspin)) then
         if (nspin /= 1) then
            msg = "trc_tddft: only a closed-shell restricted reference is implemented (UHF/UKS is refused)"
            return
         end if
      end if
      if (.not. opts%singlet) then
         msg = "trc_tddft: only singlet excitations are implemented (triplets are refused)"
      else if (opts%nroots < 1) then
         msg = "trc_tddft: nroots must be at least 1"
      else if (opts%missed_root_iterations < 0) then
         msg = "trc_tddft: missed_root_iterations must be zero or more"
      else if (.not. present(func)) then
         msg = ""
      else if (len_trim(func%missing) > 0) then
         msg = "trc_tddft: functional '"//trim(func%name)//"' is not available: "//trim(func%missing)
      else if (func%omega /= 0.0_dp .or. func%exx_lr /= 0.0_dp) then
         msg = "trc_tddft: range-separated hybrids ('"//trim(func%name)//"') are refused -- the response needs " &
               //"a second, short-range-screened exchange term this module does not build"
      else if (func%nlc_b /= 0.0_dp .or. func%nlc_c /= 0.0_dp) then
         msg = "trc_tddft: VV10 nonlocal correlation ('"//trim(func%name)//"') has no response term here"
      else if (func%pt2_cos /= 0.0_dp .or. func%pt2_css /= 0.0_dp) then
         msg = "trc_tddft: double hybrids ('"//trim(func%name)//"') are refused -- the PT2 part of the " &
               //"excitation energy is not implemented"
      else if (func%family /= XC_FAMILY_LDA .and. func%family /= XC_FAMILY_GGA) then
         msg = "trc_tddft: only LDA and GGA kernels are implemented ('"//trim(func%name)//"' is a meta-GGA)"
      end if
   end function trc_tddft_refuse_reason

   subroutine tddft_reference_set(this, nocc, cmo, eps, dmat, func, budget)
      !! Copy a restricted SCF result in: `cmo` (nao, nmo, occupied columns
      !! first), `eps` (nmo), `dmat` (nao, nao total density), all host, as
      !! `trc_scf_run` returns them. `func` is given for Kohn-Sham and absent
      !! for Hartree-Fock; `budget` caps the XC scratch (doubles) and, left
      !! out, leaves the routine's default.
      class(trc_tddft_reference_t), intent(inout) :: this
      integer, intent(in) :: nocc
      real(dp), intent(in) :: cmo(:, :), eps(:), dmat(:, :)
      type(trc_xc_functional_t), intent(in), optional :: func
      integer(int64), intent(in), optional :: budget
      this%nocc = nocc
      this%cmo = cmo
      this%eps = eps
      this%dmat = dmat
      this%dft = present(func)
      if (present(func)) this%func = func
      this%has_budget = present(budget)
      if (present(budget)) this%budget = budget
   end subroutine tddft_reference_set

   subroutine tddft_op_setup(this, b, la, ref, err, xg)
      !!
      !! Keep the orbitals of `ref` on the device and, for a functional, build
      !! the f_xc cache from its ground-state density. `err` is set, and
      !! nothing is mapped, if the functional is refused or its grid is
      !! missing.
      !!
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(error_t), intent(out) :: err
      type(trc_xc_grid_t), intent(in), optional :: xg
      !! The grid the SCF converged on, built and on the device; needed for a
      !! functional.
      type(trc_tddft_options_t) :: probe
      character(len=256) :: reason
      integer :: n, nocc, nvir, i, a

      call this%release()
      reason = ""
      if (ref%dft) reason = trc_tddft_refuse_reason(probe, ref%func)
      if (len_trim(reason) > 0) then
         call err%set(ERROR_VALIDATION, trim(reason))
         return
      end if
      if (ref%dft .and. .not. present(xg)) then
         call err%set(ERROR_VALIDATION, "trc_tddft: a functional needs the SCF's xc grid")
         return
      end if
      nocc = ref%nocc
      n = size(ref%cmo, 1)
      nvir = size(ref%cmo, 2) - nocc
      if (nocc < 1 .or. nvir < 1 .or. n /= b%nao .or. size(ref%eps) < size(ref%cmo, 2)) then
         call err%set(ERROR_VALIDATION, "trc_tddft: needs at least one occupied and one virtual orbital, " &
                      //"and orbitals of the basis's size")
         return
      end if
      this%nao = n
      this%nocc = nocc
      this%nvir = nvir
      this%has_budget = ref%has_budget
      this%budget = ref%budget
      allocate (this%cocc(n, nocc), this%cvir(n, nvir), this%de(nocc, nvir))
      this%cocc = ref%cmo(:, 1:nocc)
      this%cvir = ref%cmo(:, nocc + 1:nocc + nvir)
      do a = 1, nvir
         do i = 1, nocc
            this%de(i, a) = ref%eps(nocc + a) - ref%eps(i)
         end do
      end do
      this%cx = 1.0_dp
      this%dft = ref%dft
      if (this%dft) then
         this%func = ref%func
         this%cx = ref%func%exx
      end if
      !$acc enter data copyin(this%cocc, this%cvir, this%de)
      if (this%dft) call this%cache%setup(b, xg, ref%func, ref%dmat, la)
      this%napply_screened = 0
      this%ready = .true.
   end subroutine tddft_op_setup

   subroutine tddft_add_diag(nov, nvec, de, x, ax)
      !! ax += (e_a - e_i) x for every vector, over the amplitudes flattened
      !! to one index (the callers pass their (nocc, nvir, nvec) arrays by
      !! sequence association): a loop over small run-time extents is one
      !! index here, not a multi-index do concurrent.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: de(nov), x(nov, nvec)
      real(dp), intent(inout) :: ax(nov, nvec)
      integer :: pair, k
      do concurrent(k=1:nvec, pair=1:nov)
         ax(pair, k) = ax(pair, k) + de(pair)*x(pair, k)
      end do
   end subroutine tddft_add_diag

   subroutine tddft_op_apply(this, b, eri, la, nvec, x, ax, xg)
      !!
      !! ax(:, :, k) := A x(:, :, k) for k = 1..nvec, x and ax HOST arrays
      !! (nocc, nvir, nvec): ONE `fock_split_many_lent` and ONE
      !! `fxc_apply_many` call whatever nvec is. tmat and gmat stay on the
      !! device between the two cuBLAS transforms, so what crosses per call is
      !! x up and ax down (nocc nvir nvec doubles each), not nao^2 nvec; the
      !! host still holds tmat, gmat and (functional) vxc, three
      !! (nao, nao, nvec) arrays, as the mapped arrays' host side, and the
      !! lent build adds device scratch of two (nvec, nao_w, nao_w) arrays and
      !! two (nao_w, nao_w). `xg` is needed exactly when `setup` was given a
      !! functional.
      !!
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nvec
      real(dp), intent(in) :: x(this%nocc, this%nvir, nvec)
      real(dp), intent(out) :: ax(this%nocc, this%nvir, nvec)
      type(trc_xc_grid_t), intent(in), optional :: xg
      !! Needed exactly when `setup` was given a functional.
      real(dp), allocatable :: tmat(:, :, :), gmat(:, :, :), vxc(:, :, :), half(:, :), half_occ(:, :)
      integer :: n, nocc, nvir, k

      if (.not. this%ready) error stop "trc_tddft: operator used before setup"
      if (this%dft .and. .not. present(xg)) error stop "trc_tddft: a functional's matvec needs the xc grid"
      if (nvec < 1) return
      n = this%nao
      nocc = this%nocc
      nvir = this%nvir
      allocate (tmat(n, n, nvec), gmat(n, n, nvec), half(n, nvir), half_occ(nocc, n))
      if (this%dft) then
         allocate (vxc(n, n, nvec))
      else
         allocate (vxc(1, 1, 1))
      end if
      !$acc enter data copyin(x) create(ax, tmat, gmat, vxc, half, half_occ)
      do k = 1, nvec
         call la%gemm_mat(this%cocc, x(:, :, k), half)
         call la%gemm_mat(half, this%cvir, tmat(:, :, k), transb="T")
      end do
      call eri%fock_split_many_lent(b, nvec, tmat, gmat, &
                                    trc_fock_options_t(j_scale=2.0_dp, k_scale=2.0_dp*this%cx, &
                                                      density_screen=this%density_screen, &
                                                      count_survivors=this%count_survivors))
      if (this%density_screen) this%napply_screened = this%napply_screened + 1
      if (this%dft) then
         ! Without a budget the routine's default keeps every vector in one pass.
         if (this%has_budget) then
            call trc_xc_fxc_rks_apply_many_lent(b, xg, this%cache, nvec, tmat, vxc, this%budget, la)
         else
            call trc_xc_fxc_rks_apply_many_lent(b, xg, this%cache, nvec, tmat, vxc, la=la)
         end if
      end if
      do k = 1, nvec
         call la%gemm_mat(this%cocc, gmat(:, :, k), half_occ, transa="T")
         call la%gemm_mat(half_occ, this%cvir, ax(:, :, k))
         if (this%dft) then
            call la%gemm_mat(this%cocc, vxc(:, :, k), half_occ, transa="T")
            call la%gemm_mat(half_occ, this%cvir, ax(:, :, k), alpha=2.0_dp, beta=1.0_dp)
         end if
      end do
      call tddft_add_diag(nocc*nvir, nvec, this%de, x, ax)
      !$acc update self(ax)
      !$acc exit data delete(x, ax, tmat, gmat, vxc, half, half_occ)
      deallocate (tmat, gmat, vxc, half, half_occ)
   end subroutine tddft_op_apply

   subroutine tddft_add_diag_rpa(nov, nvec, de, x, ax)
      !! ax(:, :, 1) += (e_a - e_i) x and ax(:, :, 2) += (e_a - e_i) x, the
      !! diagonal of both A+B and A-B, as `tddft_add_diag` over the flattened
      !! amplitudes.
      integer, intent(in) :: nov, nvec
      real(dp), intent(in) :: de(nov), x(nov, nvec)
      real(dp), intent(inout) :: ax(nov, nvec, 2)
      integer :: pair, k
      do concurrent(k=1:nvec, pair=1:nov)
         ax(pair, k, 1) = ax(pair, k, 1) + de(pair)*x(pair, k)
         ax(pair, k, 2) = ax(pair, k, 2) + de(pair)*x(pair, k)
      end do
   end subroutine tddft_add_diag_rpa

   subroutine tddft_split_density(n, nvec, antisym, tmat)
      !! tmat(:, :, k) (k = 1..nvec) holds the transition densities T; on
      !! return tmat(:, :, k) is the symmetric part (T + T^T)/2 and, when
      !! `antisym`, tmat(:, :, nvec + k) the antisymmetric part (T - T^T)/2.
      !! Device-resident (one density at a time), so `tmat` never visits the
      !! host.
      integer, intent(in) :: n, nvec
      logical, intent(in) :: antisym
      real(dp), intent(inout) :: tmat(n, n, *)
      integer :: k
      do k = 1, nvec
         if (antisym) call tddft_antisym_part(n, tmat(:, :, k), tmat(:, :, nvec + k))
         call tddft_sym_in_place(n, tmat(:, :, k))
      end do
   end subroutine tddft_split_density

   subroutine tddft_antisym_part(n, t, a)
      !! a := (t - t^T)/2.
      integer, intent(in) :: n
      real(dp), intent(in) :: t(n, n)
      real(dp), intent(out) :: a(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         a(i, j) = 0.5_dp*(t(i, j) - t(j, i))
      end do
   end subroutine tddft_antisym_part

   subroutine tddft_sym_in_place(n, t)
      !! t := (t + t^T)/2 in place: the iteration at (i, j), i > j, owns both
      !! (i, j) and (j, i), and the diagonal (t + t)/2 is t exactly.
      integer, intent(in) :: n
      real(dp), intent(inout) :: t(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         if (i > j) then
            t(i, j) = 0.5_dp*(t(i, j) + t(j, i))
            t(j, i) = t(i, j)
         end if
      end do
   end subroutine tddft_sym_in_place

   subroutine tddft_zero_matrix(n, g)
      !! g := 0, on the device.
      integer, intent(in) :: n
      real(dp), intent(out) :: g(n, n)
      integer :: i, j
      do concurrent(j=1:n, i=1:n)
         g(i, j) = 0.0_dp
      end do
   end subroutine tddft_zero_matrix

   subroutine tddft_op_apply_rpa(this, b, eri, la, nvec, x, ax, xg)
      !!
      !! ax(:, :, k, 1) := (A+B) x(:, :, k) and ax(:, :, k, 2) := (A-B)
      !! x(:, :, k) for k = 1..nvec, x and ax HOST arrays (nocc, nvir, nvec
      !! [, 2]): ONE `fock_folded_many_lent` call on the symmetric parts of the
      !! transition densities (j_scale = 4, k_scale = 4 c_x), a second on the
      !! antisymmetric parts, K only (j_scale = 0, `antisymmetric`), only when
      !! c_x /= 0 or `force_antisym`, and ONE `fxc_apply_many` call (the
      !! symmetric parts), whatever nvec is. With `force_antisym` and c_x = 0
      !! the antisymmetric half of the Fock result is zero and is set without
      !! a launch. The module header derives the products. Scratch is
      !! (nao, nao, 2 nvec) twice, plus (nao, nao, nvec) for a functional.
      !!
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nvec
      real(dp), intent(in) :: x(this%nocc, this%nvir, nvec)
      real(dp), intent(out) :: ax(this%nocc, this%nvir, nvec, 2)
      type(trc_xc_grid_t), intent(in), optional :: xg
      !! Needed exactly when `setup` was given a functional.
      real(dp), allocatable :: tmat(:, :, :), gmat(:, :, :), vxc(:, :, :), half(:, :), half_occ(:, :)
      integer :: n, nocc, nvir, k, ndens
      logical :: antisym
      integer(int64) :: nkept_sym

      if (.not. this%ready) error stop "trc_tddft: operator used before setup"
      if (this%dft .and. .not. present(xg)) error stop "trc_tddft: a functional's matvec needs the xc grid"
      if (nvec < 1) return
      n = this%nao
      nocc = this%nocc
      nvir = this%nvir
      antisym = this%cx /= 0.0_dp .or. this%force_antisym
      ndens = nvec
      if (antisym) ndens = 2*nvec
      allocate (tmat(n, n, ndens), gmat(n, n, ndens), half(n, nvir), half_occ(nocc, n))
      if (this%dft) then
         allocate (vxc(n, n, nvec))
      else
         allocate (vxc(1, 1, 1))
      end if
      !$acc enter data copyin(x) create(ax, tmat, gmat, vxc, half, half_occ)
      do k = 1, nvec
         call la%gemm_mat(this%cocc, x(:, :, k), half)
         call la%gemm_mat(half, this%cvir, tmat(:, :, k), transb="T")
      end do
      call tddft_split_density(n, nvec, antisym, tmat)
      call eri%fock_folded_many_lent(b, nvec, tmat(:, :, 1:nvec), gmat(:, :, 1:nvec), &
                                     trc_fock_options_t(j_scale=4.0_dp, k_scale=4.0_dp*this%cx, &
                                                       density_screen=this%density_screen, &
                                                       count_survivors=this%count_survivors))
      nkept_sym = eri%nkept
      if (antisym) then
         if (this%cx /= 0.0_dp) then
            call eri%fock_folded_many_lent(b, nvec, tmat(:, :, nvec + 1:ndens), gmat(:, :, nvec + 1:ndens), &
                                           trc_fock_options_t(j_scale=0.0_dp, k_scale=4.0_dp*this%cx, &
                                                             density_screen=this%density_screen, &
                                                             count_survivors=this%count_survivors), &
                                           antisymmetric=.true.)
            eri%nkept = eri%nkept + nkept_sym
         else
            do k = 1, nvec
               call tddft_zero_matrix(n, gmat(:, :, nvec + k))
            end do
         end if
      end if
      if (this%density_screen) this%napply_screened = this%napply_screened + 1
      if (this%dft) then
         if (this%has_budget) then
            call trc_xc_fxc_rks_apply_many_lent(b, xg, this%cache, nvec, tmat, vxc, this%budget, la)
         else
            call trc_xc_fxc_rks_apply_many_lent(b, xg, this%cache, nvec, tmat, vxc, la=la)
         end if
      end if
      do k = 1, nvec
         call la%gemm_mat(this%cocc, gmat(:, :, k), half_occ, transa="T")
         call la%gemm_mat(half_occ, this%cvir, ax(:, :, k, 1))
         if (this%dft) then
            call la%gemm_mat(this%cocc, vxc(:, :, k), half_occ, transa="T")
            call la%gemm_mat(half_occ, this%cvir, ax(:, :, k, 1), alpha=4.0_dp, beta=1.0_dp)
         end if
         if (antisym) then
            call la%gemm_mat(this%cocc, gmat(:, :, nvec + k), half_occ, transa="T")
            call la%gemm_mat(half_occ, this%cvir, ax(:, :, k, 2))
         else
            ! A-B is the diagonal: a zero (alpha = 0 never reads the operands'
            ! values) for the diagonal to be added to, on the device as well.
            call la%gemm_mat(half_occ, this%cvir, ax(:, :, k, 2), alpha=0.0_dp)
         end if
      end do
      call tddft_add_diag_rpa(nocc*nvir, nvec, this%de, x, ax)
      !$acc update self(ax)
      !$acc exit data delete(x, ax, tmat, gmat, vxc, half, half_occ)
      deallocate (tmat, gmat, vxc, half, half_occ)
   end subroutine tddft_op_apply_rpa

   subroutine tddft_op_dipoles(this, b, pl, la, nvec, x, tdip)
      !! tdip(:, k) = sqrt(2) sum_ia x(i, a, k) <i|r|a>, the singlet transition
      !! dipole of an excitation with unit-norm amplitudes, origin at zero
      !! (immaterial: <i|a> = 0). x host (nocc, nvir, nvec).
      class(trc_tddft_op_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_linalg_t), intent(inout) :: la
      integer, intent(in) :: nvec
      real(dp), intent(in) :: x(this%nocc, this%nvir, nvec)
      real(dp), intent(out) :: tdip(3, nvec)
      real(dp), allocatable :: mmat(:, :, :), dip(:, :, :), dmo(:, :, :), half(:, :)
      integer :: n, c, k

      if (.not. this%ready) error stop "trc_tddft: operator used before setup"
      n = this%nao
      ! All TRC_NMULT components are evaluated for the 3 dipoles: terco has no
      ! dipole-only multipole routine (trc_cphf/trc_cpks do the same), and this
      ! runs once per solve.
      allocate (mmat(n, n, TRC_NMULT))
      call trc_multipoles(b, pl, [0.0_dp, 0.0_dp, 0.0_dp], mmat)
      allocate (dip(n, n, 3), dmo(this%nocc, this%nvir, 3), half(n, this%nvir))
      dip = mmat(:, :, 1:3)
      deallocate (mmat)
      !$acc enter data copyin(dip) create(dmo, half)
      do c = 1, 3
         call la%gemm_mat(dip(:, :, c), this%cvir, half)
         call la%gemm_mat(this%cocc, half, dmo(:, :, c), transa="T")
      end do
      !$acc update self(dmo)
      !$acc exit data delete(dip, dmo, half)
      do k = 1, nvec
         do c = 1, 3
            tdip(c, k) = sqrt(2.0_dp)*sum(x(:, :, k)*dmo(:, :, c))
         end do
      end do
      deallocate (dip, dmo, half)
   end subroutine tddft_op_dipoles

   subroutine tddft_op_release(this)
      class(trc_tddft_op_t), intent(inout) :: this
      if (this%ready) then
         !$acc exit data delete(this%cocc, this%cvir, this%de)
      end if
      call this%cache%release()
      if (allocated(this%cocc)) deallocate (this%cocc)
      if (allocated(this%cvir)) deallocate (this%cvir)
      if (allocated(this%de)) deallocate (this%de)
      this%ready = .false.
      this%dft = .false.
      this%nao = 0
      this%nocc = 0
      this%nvir = 0
      this%budget = 0
   end subroutine tddft_op_release

   subroutine tddft_fail(res, code, msg)
      type(trc_tddft_result_t), intent(inout) :: res
      integer, intent(in) :: code
      character(len=*), intent(in) :: msg
      res%converged = .false.
      call res%err%set(code, msg)
   end subroutine tddft_fail

   subroutine tddft_report_window(res, nopen)
      !! Record the candidates the missed-root check left open (`tddft_window`):
      !! `res%window_open` / `res%nwindow_open`, and one warning in the log, which names
      !! the option that iterates them.
      type(trc_tddft_result_t), intent(inout) :: res
      integer, intent(in) :: nopen
      character(len=*), parameter :: HEAD = "trc_tddft: the missed-root check left "
      character(len=*), parameter :: TAIL = " candidate Ritz pair(s) not converged; the returned roots are converged, "// &
                                            "but a root below the highest of them, in a symmetry the guesses do "// &
                                            "not reach, is not excluded, and then a returned root may be a "// &
                                            "different state than the one wanted; raise missed_root_iterations "// &
                                            "(--tddft-missed-root-iterations) to iterate the candidates"
      character(len=:), allocatable :: num, msg
      res%nwindow_open = nopen
      res%window_open = nopen > 0
      if (nopen > 0) then
         ! Joined into an exact-length local: one nvfortran release mispasses the length of a message built in the call.
         num = to_string(nopen)
         allocate (character(len=len(HEAD) + len(num) + len(TAIL)) :: msg)
         msg = HEAD//num//TAIL
         call global_logger%warning(msg)
      end if
   end subroutine tddft_report_window

   subroutine tddft_orthonormalise(nov, m, mmax, vs, nc, tc, nnew)
      !! Orthonormalise the columns of tc(:, 1:nc) against vs(:, 1:m) and one
      !! another (classical Gram-Schmidt, twice), keeping those that stay
      !! above DROP_REL of their starting norm. The survivors are compacted into
      !! the first nnew columns of tc and appended to vs after column m. At most
      !! mmax - m are kept: vs has mmax columns, and when mmax is the problem
      !! dimension roundoff must not let a dependent vector past it.
      integer, intent(in) :: nov, m, mmax, nc
      real(dp), intent(inout) :: vs(nov, *), tc(nov, nc)
      integer, intent(out) :: nnew
      real(dp) :: n0, n1
      real(dp) :: proj(max(m + nc, 1))
      integer :: j, p, q, mm
      nnew = 0
      mm = m
      do j = 1, nc
         n0 = sqrt(sum(tc(:, j)**2))
         if (mm >= mmax) exit
         if (.not. (n0 > 0.0_dp)) cycle
         do p = 1, 2
            do q = 1, mm
               proj(q) = sum(vs(:, q)*tc(:, j))
            end do
            do q = 1, mm
               tc(:, j) = tc(:, j) - proj(q)*vs(:, q)
            end do
         end do
         n1 = sqrt(sum(tc(:, j)**2))
         if (.not. (n1 > DROP_REL*n0)) cycle
         mm = mm + 1
         nnew = nnew + 1
         vs(:, mm) = tc(:, j)/n1
      end do
      do j = 1, nnew
         tc(:, j) = vs(:, m + j)
      end do
   end subroutine tddft_orthonormalise

   subroutine tddft_probe_vectors(sp)
      !! Fill sp%tp with the probe: the unit vectors of the sp%nprobe
      !! next-lowest orbital-energy differences (none of which the guesses
      !! held) and sp%nrand random vectors. A random vector has weight in every
      !! symmetry sector, which no unit vector of the lowest differences can
      !! promise: a root whose sector has no pair among them is reached only
      !! through it. Each is weighted by 1/(e_a - e_i), which favours the
      !! low-lying pairs the low roots are made of. The numbers come from pic's
      !! SplitMix64 with a fixed seed (`pic_rng`: bit-identical on every
      !! compiler, thread count and rank), so a run reproduces.
      type(tddft_space_t), intent(inout) :: sp
      type(splitmix64_t) :: gen
      integer :: j, k, nov
      nov = size(sp%tp, 1)
      sp%tp = 0.0_dp
      do j = 1, sp%nprobe
         sp%tp(sp%ord(sp%nguess + j), j) = 1.0_dp
      end do
      call gen%seed(PROBE_SEED)
      do k = 1, sp%nrand
         do j = 1, nov
            sp%tp(j, sp%nprobe + k) = (2.0_dp*next_real_dp(gen) - 1.0_dp)/max(sp%dgn(j), PRECOND_FLOOR)
         end do
      end do
   end subroutine tddft_probe_vectors

   subroutine tddft_window(sp, nr, theta, rnorm, nkc, track, buffer_tol, win)
      !! Which of the nkc Ritz pairs (theta ascending, residual norms rnorm) are
      !! iterated, and how many candidates the check leaves open.
      !!
      !! The lowest sp%nt are always iterated. A further pair is a CANDIDATE
      !! when it may hold a root that belongs in the wanted window and is not
      !! converged: its Ritz interval [theta - |r|, theta + |r|] reaches below
      !! the nr-th Ritz value (some eigenvalue lies within |r| of theta), its
      !! residual is at least `buffer_tol`, and |r| < theta/2. A pair with
      !! |r| >= theta/2 is too poorly resolved to say anything and is neither
      !! iterated nor counted: the random probe's Ritz pairs in the high part of
      !! the spectrum are of that kind (|r| = 0.5 to 2 theta), and counting
      !! them as open (|r| < theta) flagged most of the 52 sweep cases, and
      !! iterating them too cost 20 to 55% more matvecs on HNC/6-31+G/HF and left
      !! them open anyway. So this check does not speak for such pairs: a root that only
      !! one of them could reach is outside what `window_open` reports (header,
      !! WHAT THE CHECK DOES NOT COVER).
      !! With `track` (after the probe, for at most `missed_root_iterations` iterations) the
      !! candidates are iterated, the lowest first, until sp%nkcap - sp%nt of
      !! them are taken. ONE COUNT per candidate, so a candidate at index 25
      !! costs one slot, not 25 (win%act marks the kept indices,
      !! win%kept_index(1:win%nact) lists them in order, and a collapse keeps
      !! exactly those). win%nopen counts every
      !! candidate, iterated or not; evaluated at the end of a solve, where every
      !! iterated one has converged to `buffer_tol` except after an early exit,
      !! it is the number of candidates left open (not iterated, cut short by the
      !! cap or the budget, or abandoned because the correction lay in the
      !! space), and the solver reports it as `window_open`.
      type(tddft_space_t), intent(in) :: sp
      integer, intent(in) :: nr, nkc
      real(dp), intent(in) :: theta(:), rnorm(:)
      real(dp), intent(in) :: buffer_tol
      logical, intent(in) :: track
      type(tddft_window_t), intent(inout) :: win
      integer :: k, n
      win%nact = sp%nt
      win%act = .false.
      win%act(1:sp%nt) = .true.
      win%nopen = 0
      do k = sp%nt + 1, nkc
         if (.not. (theta(k) - rnorm(k) < theta(nr) .and. rnorm(k) < 0.5_dp*theta(k) .and. rnorm(k) >= buffer_tol)) cycle
         win%nopen = win%nopen + 1
         if (track .and. win%nact < sp%nkcap) then
            win%act(k) = .true.
            win%nact = win%nact + 1
         end if
      end do
      n = 0
      do k = 1, nkc
         if (win%act(k)) then
            n = n + 1
            win%kept_index(n) = k
         end if
      end do
   end subroutine tddft_window

   subroutine tddft_probe(sp, xr, axr, nnew)
      !! Append the probe (`tddft_probe_vectors`) to the subspace, collapsing
      !! onto the Ritz vectors xr/axr first if there is no room. nnew is the
      !! number kept (0: all lay inside the space); sp%m advances by it.
      type(tddft_space_t), intent(inout) :: sp
      real(dp), intent(in) :: xr(:, :), axr(:, :)
      integer, intent(out) :: nnew
      integer :: nt, np
      nt = size(xr, 2)
      np = sp%nprobe + sp%nrand
      call tddft_probe_vectors(sp)
      if (sp%m + np > sp%mmax) then
         sp%vs(:, 1:nt) = xr
         sp%avs(:, 1:nt) = axr
         sp%m = nt
         sp%ncollapse = sp%ncollapse + 1
      end if
      call tddft_orthonormalise(size(xr, 1), sp%m, sp%mmax, sp%vs, np, sp%tp, nnew)
      sp%m = sp%m + nnew
   end subroutine tddft_probe

   subroutine tddft_rpa_collapse(sp, uc, vc)
      !! Collapse the RPA subspace onto the span of the current Ritz vectors p
      !! and m, whose coefficients in the subspace basis are the columns of
      !! uc and vc (m, nt): the orthonormalised span becomes the new basis and
      !! both products (A+B) V and (A-B) V are transformed with it, so no
      !! matvec is spent. A pair numerically inside the span is dropped.
      type(tddft_space_t), intent(inout) :: sp
      real(dp), intent(in) :: uc(:, :), vc(:, :)
      real(dp), allocatable :: cc(:, :), qm(:, :), work(:, :)
      integer :: m, nt, nq, nov
      m = sp%m
      nt = size(uc, 2)
      nov = size(sp%vs, 1)
      allocate (cc(m, 2*nt), qm(m, 2*nt), work(nov, 2*nt))
      cc(:, 1:nt) = uc(1:m, :)
      cc(:, nt + 1:2*nt) = vc(1:m, :)
      call tddft_orthonormalise(m, 0, 2*nt, qm, 2*nt, cc, nq)
      call pic_gemm(sp%vs(:, 1:m), qm(:, 1:nq), work(:, 1:nq), "N", "N", 1.0_dp, 0.0_dp)
      sp%vs(:, 1:nq) = work(:, 1:nq)
      call pic_gemm(sp%avs(:, 1:m), qm(:, 1:nq), work(:, 1:nq), "N", "N", 1.0_dp, 0.0_dp)
      sp%avs(:, 1:nq) = work(:, 1:nq)
      call pic_gemm(sp%ams(:, 1:m), qm(:, 1:nq), work(:, 1:nq), "N", "N", 1.0_dp, 0.0_dp)
      sp%ams(:, 1:nq) = work(:, 1:nq)
      sp%m = nq
      sp%ncollapse = sp%ncollapse + 1
   end subroutine tddft_rpa_collapse

   subroutine tddft_rpa_probe(sp, uc, vc, nnew, collapsed)
      !! As `tddft_probe` for the RPA subspace, after a collapse
      !! onto the Ritz vectors (uc, vc: `tddft_rpa_collapse`) if there is no room.
      !! uc and vc are coefficients in the basis the iteration started with, so
      !! once this iteration has collapsed (`collapsed`) they no longer describe
      !! sp%vs and a second collapse would keep arbitrary combinations; the
      !! space then holds the Ritz span already and the probe is only appended
      !! (what does not fit is dropped by the orthonormalisation).
      type(tddft_space_t), intent(inout) :: sp
      real(dp), intent(in) :: uc(:, :), vc(:, :)
      integer, intent(out) :: nnew
      logical, intent(in) :: collapsed
      integer :: np
      np = sp%nprobe + sp%nrand
      call tddft_probe_vectors(sp)
      if (sp%m + np > sp%mmax .and. .not. collapsed) call tddft_rpa_collapse(sp, uc, vc)
      call tddft_orthonormalise(size(sp%vs, 1), sp%m, sp%mmax, sp%vs, np, sp%tp, nnew)
      sp%m = sp%m + nnew
   end subroutine tddft_rpa_probe

   subroutine tddft_rpa_subspace(sp, ap, am, w2, zc, err)
      !! The RPA subspace problem (module header): Ap = V^T (A+B) V and
      !! Am = V^T (A-B) V of the m columns in use, symmetrised, then the
      !! eigenvalues w2 (ascending) of Am^{1/2} Ap Am^{1/2} and its
      !! eigenvectors zc. On return `am` holds Am^{1/2}. `err` is set when Am
      !! or the product is not positive definite or LAPACK fails. Only the
      !! subspace is tested: instability along a direction V has not reached is
      !! not detected.
      type(tddft_space_t), intent(in) :: sp
      real(dp), intent(inout) :: ap(:, :), am(:, :), zc(:, :)
      real(dp), intent(out) :: w2(:)
      character(len=*), intent(out) :: err
      real(dp), allocatable :: sam(:), work(:, :)
      integer :: m, j, k
      integer(default_int) :: linfo
      m = sp%m
      err = ""
      allocate (sam(m), work(m, m))
      call pic_gemm(sp%vs(:, 1:m), sp%avs(:, 1:m), ap(1:m, 1:m), "T", "N", 1.0_dp, 0.0_dp)
      call pic_gemm(sp%vs(:, 1:m), sp%ams(:, 1:m), am(1:m, 1:m), "T", "N", 1.0_dp, 0.0_dp)
      do j = 1, m
         do k = 1, j - 1
            ap(k, j) = 0.5_dp*(ap(k, j) + ap(j, k))
            ap(j, k) = ap(k, j)
            am(k, j) = 0.5_dp*(am(k, j) + am(j, k))
            am(j, k) = am(k, j)
         end do
      end do
      zc(1:m, 1:m) = am(1:m, 1:m)
      call pic_syev(zc(1:m, 1:m), sam, "V", "U", linfo)
      if (linfo /= 0 .or. .not. all(ieee_is_finite(sam))) then
         err = "trc_tddft: the subspace eigenproblem broke down (NaN or LAPACK failure)"
         return
      end if
      if (.not. (sam(1) > 0.0_dp)) then
         err = "trc_tddft: A-B is not positive definite in the Davidson subspace ("//to_string(sam(1), "es10.2") &
               //"): the reference may be unstable and the RPA roots complex"
         return
      end if
      do j = 1, m
         work(:, j) = zc(1:m, j)*sqrt(sam(j))
      end do
      call pic_gemm(work, zc(1:m, 1:m), am(1:m, 1:m), "N", "T", 1.0_dp, 0.0_dp)
      call pic_gemm(ap(1:m, 1:m), am(1:m, 1:m), work, "N", "N", 1.0_dp, 0.0_dp)
      call pic_gemm(am(1:m, 1:m), work, zc(1:m, 1:m), "N", "N", 1.0_dp, 0.0_dp)
      do j = 1, m
         do k = 1, j - 1
            zc(k, j) = 0.5_dp*(zc(k, j) + zc(j, k))
         end do
      end do
      call pic_syev(zc(1:m, 1:m), w2(1:m), "V", "U", linfo)
      if (linfo /= 0 .or. .not. all(ieee_is_finite(w2(1:m)))) then
         err = "trc_tddft: the subspace eigenproblem broke down (NaN or LAPACK failure)"
      else if (.not. (w2(1) > 0.0_dp)) then
         err = "trc_tddft: a non-positive RPA root w^2 ("//to_string(w2(1), "es10.2")//" Eh^2): A+B is not " &
               //"positive definite in the Davidson subspace, the reference may be unstable"
      end if
   end subroutine tddft_rpa_subspace

   subroutine tddft_rpa_iterate(b, eri, la, op, xg, opts, sp, res)
      !!
      !! The RPA block Davidson iteration of the module header, from the
      !! guesses in sp%vs(:, 1:sp%m) (set by `trc_tddft_solve_rks`, which also
      !! sizes the subspace). Fills res%omega, xamp, yamp, resid and space and
      !! sets res%converged, or sets res%err; the dipoles are the caller's.
      !! The control flow (buffer vectors, one probe before acceptance,
      !! breakdown, NaN, non-convergence) is that of the TDA loop.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_op_t), intent(inout) :: op
      type(trc_xc_grid_t), intent(in), optional :: xg
      type(trc_tddft_options_t), intent(in) :: opts
      type(tddft_space_t), intent(inout) :: sp
      type(trc_tddft_result_t), intent(inout) :: res

      real(dp), allocatable :: ap(:, :), am(:, :), zc(:, :), w2(:), om(:), uc(:, :), vc(:, :)
      real(dp), allocatable :: pr(:, :), mr(:, :), r1(:, :), r2(:, :), tc(:, :), rnorm(:), dgn(:), prod(:, :, :, :)
      real(dp) :: buffer_threshold, dx, dy, nrm
      integer :: nov, nocc, nvir, nr, nt, nnew, nunc, nmain, m, j, k, ioff, c0, nkc, xiter
      type(tddft_window_t) :: win
      logical :: collapsed
      character(len=256) :: reason
      character(len=1024) :: logline
      logical :: done, probed

      nocc = op%nocc
      nvir = op%nvir
      nov = nocc*nvir
      nr = opts%nroots
      nt = sp%nt
      buffer_threshold = max(BUFFER_TOL, 100.0_dp*opts%conv_tol)
      allocate (ap(sp%mmax, sp%mmax), am(sp%mmax, sp%mmax), zc(sp%mmax, sp%mmax), w2(sp%mmax))
      ! Room for every Ritz pair of the subspace once the probe has run.
      allocate (om(sp%mmax), uc(sp%mmax, sp%mmax), vc(sp%mmax, sp%mmax), pr(nov, sp%mmax), mr(nov, sp%mmax))
      allocate (r1(nov, sp%mmax), r2(nov, sp%mmax), tc(nov, 2*sp%mmax), rnorm(sp%mmax), dgn(nov))
      allocate (win%act(sp%mmax), win%kept_index(sp%mmax))
      dgn = reshape(op%de, [nov])
      probed = .false.
      done = .false.
      xiter = 0
      win%nopen = 0
      nnew = sp%m
      res%iterations = 0
      res%nmatvec = 0

      do
         ! (A+B) and (A-B) on the vectors added since the last pass: one integral pass.
         if (nnew > 0) then
            c0 = sp%m - nnew + 1
            allocate (prod(nocc, nvir, nnew, 2))
            call op%apply_rpa(b, eri, la, nnew, sp%vs(:, c0:sp%m), prod, xg)
            sp%avs(:, c0:sp%m) = reshape(prod(:, :, :, 1), [nov, nnew])
            sp%ams(:, c0:sp%m) = reshape(prod(:, :, :, 2), [nov, nnew])
            deallocate (prod)
            res%nmatvec = res%nmatvec + nnew
            if (.not. all(ieee_is_finite(sp%avs(:, c0:sp%m))) .or. .not. all(ieee_is_finite(sp%ams(:, c0:sp%m)))) then
               call tddft_fail(res, ERROR_GENERIC, "trc_tddft: the matvec returned NaN or Inf")
               exit
            end if
         end if
         res%iterations = res%iterations + 1
         collapsed = .false.
         m = sp%m

         call tddft_rpa_subspace(sp, ap, am, w2, zc, reason)
         if (len_trim(reason) > 0) then
            call tddft_fail(res, ERROR_GENERIC, trim(reason))
            exit
         end if

         ! Ritz vectors p = V u, m = V v (p.m = 1), the products on them and the
         ! residuals r1 = (A+B)p - w m, r2 = (A-B)m - w p of the lowest nt.
         ! The lowest nt, or all of them once the probe has run (`tddft_window`).
         nkc = merge(m, nt, probed)
         om(1:nkc) = sqrt(w2(1:nkc))
         call pic_gemm(am(1:m, 1:m), zc(1:m, 1:nkc), uc(1:m, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         do k = 1, nkc
            uc(1:m, k) = uc(1:m, k)/sqrt(om(k))
         end do
         call pic_gemm(ap(1:m, 1:m), uc(1:m, 1:nkc), vc(1:m, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         do k = 1, nkc
            vc(1:m, k) = vc(1:m, k)/om(k)
         end do
         call pic_gemm(sp%vs(:, 1:m), uc(1:m, 1:nkc), pr(:, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         call pic_gemm(sp%vs(:, 1:m), vc(1:m, 1:nkc), mr(:, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         call pic_gemm(sp%avs(:, 1:m), uc(1:m, 1:nkc), r1(:, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         call pic_gemm(sp%ams(:, 1:m), vc(1:m, 1:nkc), r2(:, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         do k = 1, nkc
            r1(:, k) = r1(:, k) - om(k)*mr(:, k)
            r2(:, k) = r2(:, k) - om(k)*pr(:, k)
            rnorm(k) = sqrt(0.5_dp*(sum(r1(:, k)**2) + sum(r2(:, k)**2)))
         end do
         if (all(ieee_is_finite(rnorm(1:nkc)))) then
            if (probed) xiter = xiter + 1
            call tddft_window(sp, nr, om, rnorm, nkc, probed .and. xiter <= opts%missed_root_iterations, buffer_threshold, win)
         else
            win%act(1:nkc) = .true.
            win%nact = nt
            win%kept_index(1:nt) = [(k, k=1, nt)]
            win%nopen = 0
         end if
         nunc = 0
         nmain = 0
         do k = 1, nkc
            if (.not. win%act(k)) cycle
            if (.not. (rnorm(k) < merge(opts%conv_tol, buffer_threshold, k <= nr))) then
               nunc = nunc + 1
               if (k <= nr) nmain = nmain + 1
            end if
         end do
         if (opts%verbose) then
            write (logline, "(a,i0,a,i0,a,i0,a,es9.2,a,f14.8)") "   tddft rpa it ", res%iterations, ": space ", m, &
               ", unconverged ", nunc, ", max resid ", maxval(rnorm(1:nr)), ", w1 ", om(1)
            call global_logger%info(trim(logline))
         end if
         if (.not. all(ieee_is_finite(rnorm(1:nkc)))) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: a residual is NaN or Inf")
            exit
         end if
         if (nunc == 0) then
            if (.not. probed .and. sp%nprobe + sp%nrand > 0) then
               probed = .true.
               call tddft_rpa_probe(sp, uc(:, 1:nt), vc(:, 1:nt), nnew, collapsed)
               if (nnew > 0) cycle
            end if
            done = .true.
            exit
         end if
         if (res%iterations >= opts%max_iteration) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: not converged in "//to_string(opts%max_iteration)// &
                            " iterations, residual "//to_string(maxval(rnorm(1:nr)), "es9.2")//" (conv_tol "// &
                            to_string(opts%conv_tol, "es9.2")//")")
            exit
         end if

         ! Preconditioned corrections of X and of Y for each unconverged root, packed first.
         ioff = 0
         do k = 1, nkc
            if (.not. win%act(k)) cycle
            if (rnorm(k) < merge(opts%conv_tol, buffer_threshold, k <= nr)) cycle
            ioff = ioff + 1
            do j = 1, nov
               dx = dgn(j) - om(k)
               if (abs(dx) < PRECOND_FLOOR) dx = sign(PRECOND_FLOOR, dx)
               dy = dgn(j) + om(k)
               tc(j, 2*ioff - 1) = -0.5_dp*(r1(j, k) + r2(j, k))/dx
               tc(j, 2*ioff) = 0.5_dp*(r2(j, k) - r1(j, k))/dy
            end do
         end do

         collapsed = sp%m + 2*nunc > sp%mmax
         if (collapsed) then
            call tddft_rpa_collapse(sp, uc(:, win%kept_index(1:win%nact)), vc(:, win%kept_index(1:win%nact)))
         end if
         call tddft_orthonormalise(nov, sp%m, sp%mmax, sp%vs, 2*nunc, tc, nnew)
         if (nnew == 0 .and. nmain == 0) then
            if (.not. probed .and. sp%nprobe + sp%nrand > 0) then
               probed = .true.
               call tddft_rpa_probe(sp, uc(:, 1:nt), vc(:, 1:nt), nnew, collapsed)
               if (nnew > 0) cycle
            end if
            done = .true.
            exit
         end if
         if (nnew == 0) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: breakdown, no usable new trial vector while "// &
                            to_string(nunc)//" root(s) are unconverged (residual "// &
                            to_string(maxval(rnorm(1:nr)), "es9.2")//")")
            exit
         end if
         sp%m = sp%m + nnew
      end do

      if (done) then
         call tddft_report_window(res, win%nopen)
         res%space = sp%m
         res%ncollapse = sp%ncollapse
         allocate (res%omega(nr), res%xamp(nocc, nvir, nr), res%yamp(nocc, nvir, nr), res%resid(nr))
         res%omega = om(1:nr)
         res%resid = rnorm(1:nr)
         do k = 1, nr
            nrm = sum(pr(:, k)*mr(:, k))
            if (.not. (nrm > 0.0_dp)) then
               call tddft_fail(res, ERROR_GENERIC, "trc_tddft: a Ritz vector has non-positive X.X - Y.Y")
               exit
            end if
            res%xamp(:, :, k) = reshape(0.5_dp*(pr(:, k) + mr(:, k))/sqrt(nrm), [nocc, nvir])
            res%yamp(:, :, k) = reshape(0.5_dp*(pr(:, k) - mr(:, k))/sqrt(nrm), [nocc, nvir])
         end do
         if (.not. res%err%has_error()) res%converged = .true.
      end if
   end subroutine tddft_rpa_iterate

   subroutine trc_tddft_solve_rks(b, pl, eri, la, ref, opts, res, xg)
      !!
      !! The lowest `opts%nroots` singlet roots (TDA, or the RPA when
      !! `opts%tda` is false) of the closed-shell
      !! reference `ref` (module header; `trc_tddft_reference_t%set` builds it
      !! from a `trc_scf_run` result, with the functional and its grid for
      !! Kohn-Sham and without for Hartree-Fock) and, for Kohn-Sham, `xg`, the
      !! grid the SCF converged on (built, on the device). `res%err` carries a
      !! refusal, a breakdown, a NaN or a failure to converge; `res%converged`
      !! is then false and the other fields are not to be used.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      type(trc_eri_t), intent(inout) :: eri
      type(trc_linalg_t), intent(inout) :: la
      type(trc_tddft_reference_t), intent(in) :: ref
      type(trc_tddft_options_t), intent(in) :: opts
      type(trc_tddft_result_t), intent(out) :: res
      type(trc_xc_grid_t), intent(in), optional :: xg

      type(trc_tddft_op_t) :: op
      real(dp), allocatable :: hsub(:, :), csub(:, :), wsub(:)
      real(dp), allocatable :: xr(:, :), axr(:, :), tc(:, :), rnorm(:), dgn(:), dsort(:)
      type(tddft_space_t) :: sp
      real(dp) :: den
      integer :: nov, nvir, nr, nt, nnew, nunc, nmain, nocc, i, a, j, k, ioff, nkc
      type(tddft_window_t) :: win
      integer :: xiter
      integer(default_int) :: linfo
      character(len=256) :: reason
      character(len=1024) :: logline
      logical :: done, probed
      real(dp) :: buffer_threshold

      res%converged = .false.
      res%nroots = opts%nroots
      reason = ""
      if (ref%dft) then
         reason = trc_tddft_refuse_reason(opts, ref%func, lmax=maxval(b%sh_l))
      else
         reason = trc_tddft_refuse_reason(opts, lmax=maxval(b%sh_l))
      end if
      if (len_trim(reason) > 0) then
         call tddft_fail(res, ERROR_VALIDATION, trim(reason))
         return
      end if
      call op%setup(b, la, ref, res%err, xg)
      if (res%err%has_error()) then
         res%converged = .false.
         call op%release()
         return
      end if
      op%density_screen = opts%density_screen
      if (.not. opts%density_screen) call global_logger%warning(SCREEN_OFF_MSG)

      nocc = ref%nocc
      nvir = op%nvir
      nov = nocc*nvir
      nr = opts%nroots
      if (nr > nov) then
         call tddft_fail(res, ERROR_VALIDATION, "trc_tddft: nroots = "//to_string(nr)// &
                         " exceeds the "//to_string(nov)//" occupied-virtual pairs")
         call op%release()
         return
      end if
      sp%mmax = opts%max_space
      if (sp%mmax <= 0) sp%mmax = 24*max(nr, 3)
      ! The block iterated is the nroots wanted plus a buffer of further Ritz
      ! vectors (see the header, MISSED ROOTS).
      nt = min(nr + max(2, nr/2), nov)
      buffer_threshold = max(BUFFER_TOL, 100.0_dp*opts%conv_tol)
      sp%nt = nt
      sp%mmax = min(max(sp%mmax, merge(3, 4, opts%tda)*nt), nov)
      sp%nguess = opts%nguess
      if (sp%nguess <= 0) sp%nguess = max(2*nr, nr + 4)
      sp%nguess = min(max(sp%nguess, nt), nov, sp%mmax)

      allocate (sp%vs(nov, sp%mmax), sp%avs(nov, sp%mmax), hsub(sp%mmax, sp%mmax), csub(sp%mmax, sp%mmax), wsub(sp%mmax))
      ! Room for every Ritz vector of the subspace once the probe has run (see
      ! `tddft_window`); until then only the lowest nt are formed.
      allocate (xr(nov, sp%mmax), axr(nov, sp%mmax), tc(nov, sp%mmax), rnorm(sp%mmax))
      allocate (win%act(sp%mmax), win%kept_index(sp%mmax))
      allocate (dgn(nov), dsort(nov))
      allocate (sp%ord(nov))
      ! The Ritz vectors kept across a collapse (nkcap of them, plus at most as many
      ! corrections) must fit: 2 nkcap <= mmax; the RPA keeps two vectors per pair
      ! and adds two per correction: 4 nkcap <= mmax.
      sp%nkcap = max(nt, min(nt + max(2, nr), sp%mmax/merge(2, 4, opts%tda)))
      if (.not. opts%tda) allocate (sp%ams(nov, sp%mmax))

      ! Diagonal (e_a - e_i) in the flat (i fastest) ordering of the amplitudes;
      ! the guesses are its lowest entries.
      do a = 1, nvir
         do i = 1, nocc
            dgn(i + (a - 1)*nocc) = op%de(i, a)
         end do
      end do
      ! sort_index sorts its array in place; the preconditioner needs dgn as is.
      dsort = dgn
      call sort_index(dsort, sp%ord)
      ! Do not cut a cluster of (near-)degenerate differences in two: the
      ! partner left out is a direction the solve may never generate.
      do while (sp%nguess < min(nov, sp%mmax))
         if (dsort(sp%nguess + 1) - dsort(sp%nguess) >= TIE_TOL) exit
         sp%nguess = sp%nguess + 1
      end do
      ! The probe must fit beside the Ritz vectors it may collapse onto.
      sp%nprobe = max(min(sp%nguess, nov - sp%nguess, sp%mmax - merge(1, 2, opts%tda)*nt), 0)
      ! Two random vectors where there is room beside the unit probe and the
      ! Ritz vectors a collapse keeps, and anything is left outside the space.
      sp%nrand = max(min(NRAND_PROBE, nov - sp%nguess - sp%nprobe, sp%mmax - merge(1, 2, opts%tda)*sp%nkcap - sp%nprobe), 0)
      sp%dgn = dgn
      allocate (sp%tp(nov, max(sp%nprobe + sp%nrand, 1)))
      probed = .false.
      sp%vs(:, 1:sp%nguess) = 0.0_dp
      do j = 1, sp%nguess
         sp%vs(sp%ord(j), j) = 1.0_dp
      end do
      sp%m = sp%nguess
      nnew = sp%nguess
      res%iterations = 0
      res%nmatvec = 0
      done = .false.
      xiter = 0
      win%nopen = 0

      if (.not. opts%tda) then
         call tddft_rpa_iterate(b, eri, la, op, xg, opts, sp, res)
         if (res%converged) then
            allocate (res%tdip(3, nr), res%osc(nr))
            call op%transition_dipoles(b, pl, la, nr, res%xamp + res%yamp, res%tdip)
            do k = 1, nr
               res%osc(k) = OSC_PREFACTOR*res%omega(k)*sum(res%tdip(:, k)**2)
            end do
         end if
         res%nscreened_apply = op%napply_screened
         call op%release()
         deallocate (sp%vs, sp%avs, sp%ams, hsub, csub, wsub, xr, axr, tc, rnorm, win%act, win%kept_index, dgn, dsort, &
                    sp%tp, sp%ord, sp%dgn)
         return
      end if

      do
         ! A on the vectors added since the last pass: one integral pass.
         if (nnew > 0) then
            call op%apply(b, eri, la, nnew, sp%vs(:, sp%m - nnew + 1:sp%m), sp%avs(:, sp%m - nnew + 1:sp%m), xg)
            res%nmatvec = res%nmatvec + nnew
            if (.not. all(ieee_is_finite(sp%avs(:, sp%m - nnew + 1:sp%m)))) then
               call tddft_fail(res, ERROR_GENERIC, "trc_tddft: the matvec returned NaN or Inf")
               exit
            end if
         end if
         res%iterations = res%iterations + 1

         ! Subspace problem.
         call pic_gemm(sp%vs(:, 1:sp%m), sp%avs(:, 1:sp%m), hsub(1:sp%m, 1:sp%m), "T", "N", 1.0_dp, 0.0_dp)
         do j = 1, sp%m
            do k = 1, j - 1
               hsub(k, j) = 0.5_dp*(hsub(k, j) + hsub(j, k))
            end do
         end do
         csub(1:sp%m, 1:sp%m) = hsub(1:sp%m, 1:sp%m)
         call pic_syev(csub(1:sp%m, 1:sp%m), wsub(1:sp%m), "V", "U", linfo)
         if (linfo /= 0 .or. .not. all(ieee_is_finite(wsub(1:sp%m)))) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: the subspace eigenproblem broke down (NaN or " &
                            //"LAPACK failure)")
            exit
         end if

         ! Ritz vectors, A times them, and the residuals of the lowest nt, or of
         ! all of them once the probe has run (`tddft_window` picks which count).
         nkc = merge(sp%m, nt, probed)
         call pic_gemm(sp%vs(:, 1:sp%m), csub(1:sp%m, 1:nkc), xr(:, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         call pic_gemm(sp%avs(:, 1:sp%m), csub(1:sp%m, 1:nkc), axr(:, 1:nkc), "N", "N", 1.0_dp, 0.0_dp)
         do k = 1, nkc
            tc(:, k) = axr(:, k) - wsub(k)*xr(:, k)
            rnorm(k) = sqrt(sum(tc(:, k)**2))
         end do
         if (all(ieee_is_finite(rnorm(1:nkc)))) then
            if (probed) xiter = xiter + 1
            call tddft_window(sp, nr, wsub, rnorm, nkc, probed .and. xiter <= opts%missed_root_iterations, buffer_threshold, win)
         else
            win%act(1:nkc) = .true.
            win%nact = nt
            win%kept_index(1:nt) = [(k, k=1, nt)]
            win%nopen = 0
         end if
         nunc = 0
         nmain = 0
         do k = 1, nkc
            if (.not. win%act(k)) cycle
            if (.not. (rnorm(k) < merge(opts%conv_tol, buffer_threshold, k <= nr))) then
               nunc = nunc + 1
               if (k <= nr) nmain = nmain + 1
            end if
         end do
         if (opts%verbose) then
            write (logline, "(a,i0,a,i0,a,i0,a,es9.2,a,f14.8)") "   tddft it ", res%iterations, ": space ", sp%m, &
               ", unconverged ", nunc, ", max resid ", maxval(rnorm(1:nr)), ", w1 ", wsub(1)
            call global_logger%info(trim(logline))
         end if
         if (.not. all(ieee_is_finite(rnorm(1:nkc)))) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: a residual is NaN or Inf")
            exit
         end if
         if (nunc == 0) then
            ! One probe before accepting: unit vectors of the next-lowest
            ! differences, none of which the guesses held, and random vectors.
            ! A root whose character the guesses (and so the corrections, which
            ! stay in their span) never reached enters the subspace through
            ! it; its Ritz pair then starts above the window, so from here on
            ! `tddft_window` iterates it until it is inside or provably above.
            ! One more matvec pass at least.
            if (.not. probed .and. sp%nprobe + sp%nrand > 0) then
               probed = .true.
               call tddft_probe(sp, xr(:, 1:nt), axr(:, 1:nt), nnew)
               if (nnew > 0) cycle
            end if
            done = .true.
            exit
         end if
         if (res%iterations >= opts%max_iteration) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: not converged in "//to_string(opts%max_iteration)// &
                            " iterations, residual "//to_string(maxval(rnorm(1:nr)), "es9.2")//" (conv_tol "// &
                            to_string(opts%conv_tol, "es9.2")//")")
            exit
         end if

         ! Preconditioned corrections of the unconverged roots, packed first.
         ioff = 0
         do k = 1, nkc
            if (.not. win%act(k)) cycle
            if (rnorm(k) < merge(opts%conv_tol, buffer_threshold, k <= nr)) cycle
            ioff = ioff + 1
            do j = 1, nov
               den = dgn(j) - wsub(k)
               if (abs(den) < PRECOND_FLOOR) den = sign(PRECOND_FLOOR, den)
               tc(j, ioff) = -(axr(j, k) - wsub(k)*xr(j, k))/den
            end do
         end do

         ! Collapse onto the current Ritz vectors when there is no room.
         if (sp%m + nunc > sp%mmax) then
            sp%vs(:, 1:win%nact) = xr(:, win%kept_index(1:win%nact))
            sp%avs(:, 1:win%nact) = axr(:, win%kept_index(1:win%nact))
            sp%m = win%nact
            sp%ncollapse = sp%ncollapse + 1
         end if
         call tddft_orthonormalise(nov, sp%m, sp%mmax, sp%vs, nunc, tc, nnew)
         if (nnew == 0 .and. nmain == 0) then
            ! Only buffer vectors were still unconverged and their corrections
            ! lie inside the space: nothing more to learn from them. The probe
            ! still runs once, as on the ordinary converged path.
            if (.not. probed .and. sp%nprobe + sp%nrand > 0) then
               probed = .true.
               call tddft_probe(sp, xr(:, 1:nt), axr(:, 1:nt), nnew)
               if (nnew > 0) cycle
            end if
            done = .true.
            exit
         end if
         if (nnew == 0) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: breakdown, no usable new trial vector while "// &
                            to_string(nunc)//" root(s) are unconverged (residual "// &
                            to_string(maxval(rnorm(1:nr)), "es9.2")//")")
            exit
         end if
         sp%m = sp%m + nnew
      end do

      if (done) then
         call tddft_report_window(res, win%nopen)
         if (wsub(1) <= 0.0_dp) then
            call tddft_fail(res, ERROR_GENERIC, "trc_tddft: a non-positive excitation energy ("// &
                            to_string(wsub(1), "es10.2")//" Eh): the reference may be unstable (or a triplet-like " &
                            //"instability)")
         else
            res%space = sp%m
            res%ncollapse = sp%ncollapse
            allocate (res%omega(nr), res%xamp(nocc, nvir, nr), res%yamp(nocc, nvir, nr), res%tdip(3, nr), &
                      res%osc(nr), res%resid(nr))
            res%yamp = 0.0_dp
            res%omega = wsub(1:nr)
            res%resid = rnorm(1:nr)
            do k = 1, nr
               xr(:, k) = xr(:, k)/sqrt(sum(xr(:, k)**2))
               res%xamp(:, :, k) = reshape(xr(:, k), [nocc, nvir])
            end do
            call op%transition_dipoles(b, pl, la, nr, res%xamp, res%tdip)
            do k = 1, nr
               res%osc(k) = OSC_PREFACTOR*res%omega(k)*sum(res%tdip(:, k)**2)
            end do
            res%converged = .true.
         end if
      end if
      res%nscreened_apply = op%napply_screened
      call op%release()
      deallocate (sp%vs, sp%avs, hsub, csub, wsub, xr, axr, tc, rnorm, win%act, win%kept_index, dgn, dsort, &
                 sp%tp, sp%ord, sp%dgn)
   end subroutine trc_tddft_solve_rks

end module trc_tddft
