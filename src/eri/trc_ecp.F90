!
! The scalar ECP one-electron matrix, U = U_type1 (local) + U_type2
! (semi-local), over a whole basis.
!
! Ported from libfint's src/cint_ecp.f90 and src/cint_ecp_drv.f90
! (Apache-2.0, from PySCF's pyscf/lib/gto/nr_ecp.c -- the spin-orbit term,
! `ECPtype_so_cart`, is not ported: nothing here claims it).
!
! AS `do concurrent` KERNELS, so the host and device builds run the same
! code (no separate host port, per Jorge): K1 precomputes each ECP
! CHANNEL's own radial potential on the whole grid once (plain U(r), no
! quadrature weight -- the weight is level-dependent and is folded in by
! K2 instead, once per level, which is what `trc_ecp_radial.F90`'s own
! header means by "letting every shell pair's own adaptive level loop
! slice a strided view of the SAME cached array"); K1b precomputes, once
! per (orbital shell, ECP channel group) rather than once per shell pair,
! the semi-local radial factors on the whole grid (`fac2_cache`), the
! semi-local angular factors (`ang2_cache`), and the local term's static
! centre-shift factors (`fac1_cache`) -- the three things the ORIGINAL,
! per-shell-pair port used to rebuild from scratch inside K2's own item,
! each scaling with ECP_NRS or with this module's angular-momentum
! ceilings: moving them to K1b's global, mapped-once arrays is what keeps
! K2's own item down to a handful of scalars per thread (measured:
! ecp_type2_item's `radi`/`radj` alone were ECP_NRS*LC1MAX reals each,
! ~147 KB, the dominant piece of a stack overrun that failed the kernel
! launch outright -- see this file's git history). K2 is one item per
! (orbital shell pair, ECP atom), summing that atom's own channels and
! writing its own exclusive output slot (no atomics); K3 reduces K2's
! per-atom slots into the Cartesian matrix.
!
! ADAPTIVE CONVERGENCE WITHOUT PER-LEVEL STATE. The original per-level
! call (halve the running value, add the next level's points, compare to
! the last level) needed its running value to persist ACROSS levels for
! every primitive pair (type 1) or lambda value (type 2) at once --
! `NPMAX*NPMAX` pairs' worth, in the worst case. `trc_ecp_radial.F90`'s
! `ecp_radial_type1_rad_part_levels`/`ecp_radial_type2_lab_levels`
! replace that with a single call per pair/lambda that materialises every
! level's cumulative value (same halve-then-add arithmetic, just written
! into separate columns instead of overwritten in place), so this file's
! own items process pairs/lambdas ONE AT A TIME and `pick_converged_level`
! below picks, from that small (level-count)-sized snapshot, the exact
! level the original's "two consecutive close levels" test would have
! stopped at -- the same convergence semantics, the same CPU numbers, with
! nothing here scaling with the radial grid.
!
! FIXED-SIZE, NO GENERAL CONTRACTIONS. Every scratch array a `routine seq`
! item owns is sized to this module's own ceilings (`ECP_ORB_LMAX` from
! trc_ecp_angular.F90, `ECP_LMAX`/`ECP_NRS` from trc_ecp_quadrature.F90,
! `NPMAX`, all through inc/trc_ecp_block_sizes.inc), never to a runtime shape
! (gpu-pitfalls.md's "no automatic arrays in routine seq"). The block kernels
! whose scratch grows with l are in inc/trc_ecp_block.inc, instantiated here
! (f, f) and in trc_ecp_gbra (g bra, f ket), which only the gradient of an f
! basis reaches. `trc_basis_t` never holds a general contraction
! (trc_ecp_radial.F90's own header), so every "nc" PySCF's C and libfint's
! port carry is 1 here and the corresponding loop is simply absent.
!
! NORMALISATION: common_fac_sp is folded into `trc_basis_t%sh_c` exactly
! once, in basis construction (AGENTS.md). libfint's own `scale_coeff`
! would fold it a SECOND time (it reads libcint's `env`, which does not
! carry it), so this port drops that factor here, keeping only the
! 4*pi / 16*pi^2 PySCF's own formula also carries -- ECP_DESIGN.md's note.
!
! 3-CENTRE SCREEN. K2 runs only the (shell pair, ECP atom) tasks whose
! bound can reach ECP_SCREEN_THRESH. A task is one K2 launch item, and the
! unscreened task count grows as nshell**2 * natm_ecp / 2: 24 CH3I in
! def2-SVP (600 shells, 24 ECP atoms) is 4,327,200 tasks, whose single K2
! launch took 30.3 s of a 43 s V100 run while 98.8 % of the blocks had a
! largest element below 1e-20 (2,513,255 exactly zero). The bound is built
! from every primitive of the shell and every row of the channel, so it
! assumes no exponent ordering in either. Per (orbital shell i, channel
! group g), `ecp_screen_item` returns
!     q(i,g) = sqrt( 4*pi * sum over radial levels L of  step_L *
!                    sum over the level's points k of  w_k |U_g(r_k)| F_i(r_k)**2 ),
!     F_i(r) = sum_p |c_p| (r + d)**l exp(-a_p (r - d)**2),
! with d = |R_C - R_i| and step_L the weight scale of level L (the stride of
! its points in the finest grid), the weight the `*_levels` routines of
! trc_ecp_radial.F90 give a level's points. The screen's level-L set holds
! every point of the finest grid whose index is a multiple of step_L, which
! is a superset of the points those routines add at level L (they add
! 2**L - 1 new points to the previous level's set); the extra terms are
! non-negative, so they only enlarge the bound. On the sphere of radius r about the ECP
! centre C a point is at least |r - d| and at most r + d from R_i, so F_i
! bounds the shell's largest value there; the angular integral of a product
! (type 1) and the projection onto a channel's Y_lm (type 2, norm-
! decreasing) are then at most 4*pi*F_i*F_j per radial point. Each level's
! rule has positive weights, so Cauchy-Schwarz over its points and then
! over the levels gives, for the sum Q of level rules of the product,
!     Q(F_i F_j) <= sqrt( Q(F_i**2) Q(F_j**2) ),
! and a block element, which is a combination of level rules each taken at
! whichever level converged for it (per primitive pair or per lambda), is
! bounded by  sum over the atom's groups of q(i,g) q(j,g), local and
! semi-local alike. The sum over levels, not the largest level, is what
! keeps the bound valid whatever level each component picked. q is formed
! with a running scale so q(i,g)*q(j,g) is not lost to underflow of the
! squares. U_g is the potential the quadrature actually uses (K1's cache),
! so nothing here assumes it is short-ranged: a long-ranged one gives a
! large q and nothing is dropped. A task whose bound is below
! ECP_SCREEN_THRESH is dropped: its slot is never computed and its
! convergence flag is not counted. The kept tasks are compacted into a
! list K2 runs over, and `block_buffer`/`conv_fail` hold only that list.
! check_ecp_screen.F90 asserts bound >= largest element on every task of
! several systems, through `trc_ecp_screen_stats_t`.
!
! THE GRADIENT (`trc_ecp_gradient`) reuses these items in one pass over all
! ECP atoms. A Cartesian primitive obeys
!     d/dR_c  x^l exp(-a r^2) = 2a x^(l+1) exp(-a r^2) - l_c x^(l-1) exp(-a r^2),
! so the derivative of a contracted shell is a contracted shell of l + 1
! (coefficients 2a c) and one of l - 1 (coefficients c). The "derivative
! basis" is the real shells followed by those, and a task is a (derivative
! bra shell, real ket shell, ECP atom) triple that survives the screen. The
! fold back onto the parent functions is linear, so it is applied to the
! density instead: for a derivative function mu' and component c,
! `parent(c, mu')` and `coef(c, mu')` (1, -l_c, or 0) are built once on the
! host, and the task's item reduces its block to
!     w_c = 2 sum_{mu' nu} block(mu', nu) coef(c, mu') D(parent(c, mu'), nu),
! three numbers in the task's own slot (no atomics, no nao**2 array). Only
! the bra is differentiated (D is symmetric: the ket gives the same again,
! hence the 2), and by translational invariance the ECP's own centre takes
! minus the sum, so the host adds w to the bra shell's atom and -w to the
! ECP atom. A task whose bra shell sits ON the ECP atom would add and take
! the same w from one atom, and is not run.
!
! The K1b caches are about 180 KB per (shell, channel group), so the
! gradient builds them only for the (shell, ECP atom) pairs some kept task
! reads (a fifth of them on 24 CH3I, 5 GB where the energy matrix's own
! caches take 11 GB), and splits the ECP atoms into chunks only when the
! device cannot hold them all (`ecp_grad_chunk_end`).
!
! f FUNCTIONS. The raised shell of an f function is g. Those derivative
! shells, and only those, have K1b columns of their own size (245 KB) and a
! second task kernel over their tasks, both in trc_ecp_gbra; the task list
! carries them after all the others. Every other column and task is the same
! as without f functions, the ket of a g task included, so a basis without f
! launches exactly what it did. g and higher orbital functions are refused,
! and so is an f function with a projected h (l = 5) ECP channel, whose
! angular factors would need l = 9 tables.
#include "inc/trc_dc_index.inc"

module trc_ecp
   use pic_types, only: dp, int64
   use trc_basis_types, only: trc_basis_t
   use trc_ecp_types, only: trc_ecp_t, trc_ecp_without_atoms
   use trc_ecp_quadrature, only: ECP_LMAX, ECP_LEVEL0, ECP_LEVEL_MAX, ECP_NRS, &
                                 ecp_quadrature_gauss_chebyshev
   use trc_ecp_angular, only: ECP_ORB_LMAX, ECP_BRA_LMAX, CART_POW_Y, CART_POW_Z, &
                              ecp_type1_rad_ang, ecp_type2_facs_ang
   use trc_ecp_tables, only: ECP_ANG_LMAX
   use trc_ecp_grad_types, only: trc_ecp_channels_t, trc_ecp_deriv_basis_t, trc_ecp_grad_tasks_t, &
                                 trc_ecp_k1b_cache_t, trc_ecp_chunk_t
   use trc_ecp_gbra, only: ECP_GBRA_FAC1, ECP_GBRA_FAC2, ECP_GBRA_ANG2, ECP_GBRA_COL_BYTES, &
                           trc_ecp_gbra_caches, trc_ecp_gbra_tasks
   use trc_ecp_radial, only: ecp_radial_rad_part, ecp_radial_type1_static_facs, &
                             ecp_radial_type1_rad_part_levels, ecp_radial_type2_facs_rad, &
                             ecp_radial_type2_lab_levels, ECP_NLEVEL
   use trc_sph, only: sph_cart_to_sph
   use trc_cart_perm, only: cart_perm_l2b
   use trc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use trc_device_memory, only: device_memory_free_bytes, device_memory_fits
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   implicit none
   private
   public :: trc_ecp_matrix, trc_ecp_matrix_cart, trc_ecp_screen_stats_t
   public :: trc_ecp_deriv_cart, trc_ecp_gradient, trc_ecp_gradient_refusal

   type :: trc_ecp_screen_stats_t
      !! What the 3-centre screen did, for a test. `worst_ratio` and
      !! `nzero_bound` cover the tasks that ran (every task for `screen=.false.`).
      !! The gradient fills the task counts and `nchunk` only, its tasks being
      !! (derivative shell, shell, ECP atom) ones.
      integer(int64) :: ntask = 0              !! (shell pair, ECP atom) tasks in all
      integer(int64) :: nkept = 0              !! tasks that ran
      integer(int64) :: nkept_g = 0            !! of those, the gradient's with a g derivative bra (trc_ecp_gbra)
      real(dp) :: worst_ratio = 0.0_dp         !! max over tasks of (largest |element|)/bound, bound > 0
      integer(int64) :: nzero_bound = 0        !! tasks with bound == 0 and a nonzero element
      integer :: nchunk = 0                    !! the gradient's chunks of ECP atoms (the matrix leaves 0)
   end type trc_ecp_screen_stats_t

   integer, parameter :: ECP_BLOCK_LBRA = ECP_ORB_LMAX, ECP_BLOCK_LKET = ECP_ORB_LMAX
   !! This module's instance of the block kernels (trc_ecp_block_sizes.inc).
#include "inc/trc_ecp_block_sizes.inc"
   integer, parameter :: NFMAX = NFKET                                       !! 10
   integer, parameter :: SLOT_SIZE = NFMAX*NFMAX                             !! 100: one (shell pair, atom) slot
   integer(int64), parameter :: GRAD_MIN_FREE_BYTES = 1024_int64**3
   !! Below one ECP atom's caches plus this, the gradient refuses.
   integer(int64), parameter :: GRAD_RESERVE_BYTES = 4*1024_int64**3
   !! A chunk takes further ECP atoms only while this stays free: the task
   !! kernel's local memory (20.5 KB a thread) is 3.5 GB on a V100 (nsys),
   !! reserved at launch unless an earlier kernel's reservation is still there.
   integer(int64), parameter :: GRAD_RESERVE_BRA_G_BYTES = 6*1024_int64**3
   !! The same when some task has a g derivative bra (an f basis), whose
   !! task kernel (trc_ecp_gbra) has the larger per-thread stack: 26.6 KB,
   !! 4.5 GB of local memory on a V100.

   real(dp), parameter :: ECP_SCREEN_THRESH = 1.0e-20_dp
   !! A (shell pair, ECP atom) task whose bound (see this file's header) is
   !! below this is not run: ECP_NEGLIGIBLE's scale, five orders below the
   !! 1e-15 absolute agreement check_ecp_screen.F90 requires. Each matrix
   !! element sums at most one dropped block per ECP atom, so the worst-case
   !! error is natm_ecp times this.

contains

#include "inc/trc_ecp_block.inc"

   pure subroutine ecp_k1_item(ur, rs, row0, row1, nrow, maxnp, ecp_rpow, ecp_np, ecp_e, ecp_c)
      !! K1's item: one channel group's U(r) on the whole grid. A same-module
      !! wrapper so the do concurrent body calls only this module's own
      !! items (AGENTS.md); the work is trc_ecp_radial's.
      !$acc routine seq
      integer, intent(in) :: row0, row1, nrow, maxnp
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      integer, intent(in) :: ecp_rpow(nrow), ecp_np(nrow)
      real(dp), intent(in) :: ecp_e(maxnp, nrow), ecp_c(maxnp, nrow)
      real(dp), intent(out) :: ur(0:ECP_NRS - 1)
      call ecp_radial_rad_part(ur, rs, row0, row1, nrow, maxnp, ecp_rpow, ecp_np, ecp_e, ecp_c)
   end subroutine ecp_k1_item

   subroutine trc_ecp_matrix(b, ecp, u, error)
      !!
      !! The whole ECP one-electron matrix, Cartesian internally and folded to
      !! spherical if `b%spherical`. `u` comes back allocated to b%nao_sph (or
      !! b%nao). `ecp%has_ecp` false gives back a zero matrix of the right size.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), allocatable, intent(out) :: u(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: mcart(:, :)

      call trc_ecp_matrix_cart(b, ecp, mcart, error)
      if (error%has_error()) return

      if (b%spherical) then
         allocate (u(b%nao_sph, b%nao_sph))
         ! sph_cart_to_sph is a do concurrent kernel: under mem:separate
         ! everything it reads must be mapped, as trc_scf_context's own call
         ! does (and b's shell arrays are not left mapped by the build above).
         !$acc enter data copyin(b, b%sh_l, b%sh_ao, b%sh_ao_sph, mcart) create(u)
         call sph_cart_to_sph(b, mcart, u)
         !$acc update self(u)
         !$acc exit data delete(b%sh_l, b%sh_ao, b%sh_ao_sph, mcart, u, b)
      else
         allocate (u(b%nao, b%nao))
         u = mcart
      end if
      deallocate (mcart)
   end subroutine trc_ecp_matrix

   subroutine trc_ecp_matrix_cart(b, ecp, mcart, error, screen, stats, split)
      !!
      !! The same matrix, always Cartesian (never folded), `mcart(b%nao, b%nao)`
      !! coming back freshly allocated and zeroed. `trc_ecp_matrix` is this plus
      !! the spherical fold; `trc_scf_context` calls this directly so its own
      !! Cartesian H_core (T + V_ne) picks up U BEFORE the ONE fold that already
      !! covers both (ECP_DESIGN.md: "added to H_core ... BEFORE the spherical
      !! fold" -- folding is linear, so U folded with T+V_ne together or on its
      !! own gives the same number to round-off, but doing it together is one
      !! fewer sph_cart_to_sph call per SCF geometry).
      !!
      !! The 3-centre screen drops the (shell pair, ECP atom) tasks whose bound
      !! is negligible (this file's header). `screen=.false.` runs every task
      !! instead, so a test can compare the two; `stats` returns the task
      !! counts and the bound's worst ratio to the largest element it bounds.
      !!
      !! `split` present (0 < split < b%nshell): only the shell pairs (i, j)
      !! with i > split >= j are built, the rectangle the derivative integrals
      !! need (`trc_ecp_deriv_cart`), and ONLY that rectangle is allocated:
      !! with nket = b%sh_ao(split + 1) - 1, `mcart(b%nao - nket, nket)`, row
      !! mu - nket holding AO mu. The shells' AOs must be in shell order. The
      !! square matrix over real and derivative shells together was 25 times
      !! the rectangle's size on CH3I.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), allocatable, intent(out) :: mcart(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: screen
      type(trc_ecp_screen_stats_t), intent(out), optional :: stats
      integer, intent(in), optional :: split
      logical :: do_screen
      integer :: nsplit, nrow, ncol
      type(trc_ecp_screen_stats_t) :: st

      do_screen = .true.
      if (present(screen)) do_screen = screen

      if (b%nshell > 0) then
         if (maxval(b%sh_l) > ECP_ORB_LMAX) then
            call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an orbital shell carries l > "// &
                           "ECP_ORB_LMAX (f); the ECP path is not built for g and above")
            return
         end if
         if (maxval(b%sh_np) > NPMAX) then
            call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an orbital shell has more "// &
                           "primitives than this port's fixed NPMAX ceiling")
            return
         end if
      end if

      nsplit = 0
      if (present(split)) nsplit = split
      nrow = b%nao
      ncol = b%nao
      if (nsplit /= 0) then
         if (nsplit < 0 .or. nsplit >= b%nshell) then
            call error%set(ERROR_VALIDATION, "trc_ecp_matrix: split must lie inside the shell list")
            return
         end if
         ncol = b%sh_ao(nsplit + 1) - 1
         nrow = b%nao - ncol
      end if
      allocate (mcart(nrow, ncol))
      mcart = 0.0_dp

      if (ecp%has_ecp) then
         call ecp_check_rows(b%natm, ecp, error)
         if (error%has_error()) return
         call ecp_matrix_cart(b, ecp, mcart, error, do_screen, present(stats), st, nsplit)
         if (error%has_error()) return
      end if
      if (present(stats)) stats = st
   end subroutine trc_ecp_matrix_cart

   subroutine ecp_check_rows(natm, ecp, error)
      !! A malformed or foreign ECP file must fail here, not read past the
      !! K1b caches (sized to ECP_LMAX) or the atom coordinates the rows name.
      integer, intent(in) :: natm
      type(trc_ecp_t), intent(in) :: ecp
      type(error_t), intent(inout) :: error

      if (ecp%natm /= natm) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: the ECP describes "// &
                        to_string(ecp%natm)//" atoms, the basis "//to_string(natm))
         return
      end if
      if (ecp%nrow == 0) return
      if (maxval(ecp%ecp_np) > NPMAX) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row has more "// &
                        "primitives than this port's fixed NPMAX ceiling")
         return
      end if
      if (minval(ecp%ecp_np) < 1) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row has no primitives")
         return
      end if
      if (minval(ecp%ecp_at) < 1 .or. maxval(ecp%ecp_at) > natm) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row names an atom "// &
                        "outside 1.."//to_string(natm))
         return
      end if
      if (minval(ecp%ecp_l) < -1 .or. maxval(ecp%ecp_l) > ECP_LMAX) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP channel has l outside "// &
                        "-1 (local) .. "//to_string(ECP_LMAX)//", this port's ECP_LMAX")
      end if
   end subroutine ecp_check_rows

   subroutine trc_ecp_deriv_cart(b, ecp, ia, d, error)
      !!
      !! `d(mu, nu, c) = < d mu / d R_mu,c | U_ia | nu >`, the Cartesian ECP matrix
      !! of ECP atom `ia` ALONE differentiated with respect to the centre of the
      !! bra function, c = x, y, z; `d` comes back `(b%nao, b%nao, 3)`, zeroed
      !! where nothing contributes. Both bra and ket are `b`'s Cartesian
      !! functions, whatever `b%spherical` says.
      !!
      !! THE TESTS' REFERENCE for dU/dR, not a production path: one ECP atom a
      !! call, a dense (nao, nao, 3) result, folded on the host.
      !! `trc_ecp_gradient` never calls it.
      !!
      !! Only the bra is differentiated: the ket's derivative is the transpose,
      !! and the derivative with respect to the ECP centre follows from the
      !! integral's invariance under translating all three centres together,
      !! which is why one atom's potential at a time.
      !!
      !! A Cartesian primitive obeys
      !!     d/dR_c  x^l exp(-a r^2) = 2a x^(l+1) exp(-a r^2) - l_c x^(l-1) exp(-a r^2),
      !! so the derivative of a shell is a shell of l + 1 and one of l - 1 per
      !! primitive, and no radial or angular kernel is new: the energy matrix
      !! over a basis of those uncontracted fake shells (their coefficients
      !! carry 2a c and c), against the real shells as ket, is folded back
      !! (libfint's cint_ecp_deriv, from PySCF's nr_ecp_deriv.c, does the same).
      !! The fake shells follow the real ones in the basis, so `split` builds and
      !! stores only the fake-bra, real-ket rectangle. That rectangle is the
      !! energy matrix's, whose kernels stop at f, so this reference takes no f
      !! function: l + 1 <= ECP_ORB_LMAX. `trc_ecp_gradient` does, through
      !! trc_ecp_gbra.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      integer, intent(in) :: ia
      real(dp), allocatable, intent(out) :: d(:, :, :)
      type(error_t), intent(inout) :: error

      type(trc_basis_t) :: bf
      type(trc_ecp_t) :: one
      real(dp), allocatable :: mf(:, :)
      logical, allocatable :: drop(:)
      integer, allocatable :: up(:, :), dn(:, :)
      integer :: nsh, nfk, ish, ip, np, l, nper, k, nao, ifunc, ly, lz, lx, c, row, rowd

      nao = b%nao
      nsh = b%nshell
      ! Before `drop(ia)` and the rows' atoms index anything.
      if (ecp%has_ecp) then
         if (ecp%natm /= b%natm) then
            call error%set(ERROR_VALIDATION, "trc_ecp_deriv_cart: the ECP describes "// &
                           to_string(ecp%natm)//" atoms, the basis "//to_string(b%natm))
            return
         end if
         if (ia < 1 .or. ia > b%natm) then
            call error%set(ERROR_VALIDATION, "trc_ecp_deriv_cart: atom "//to_string(ia)//" is outside 1.."// &
                           to_string(b%natm))
            return
         end if
      end if
      allocate (d(nao, nao, 3))
      d = 0.0_dp
      if (.not. ecp%has_ecp) return
      if (nsh > 0) then
         if (maxval(b%sh_l) + 1 > ECP_ORB_LMAX) then
            call error%set(ERROR_VALIDATION, "trc_ecp_deriv_cart: this dense reference raises each orbital "// &
                           "shell by one angular momentum and builds it with the energy matrix's kernels, which "// &
                           "stop at f; not available for an f function")
            return
         end if
      end if
      allocate (drop(ecp%natm))
      drop = .true.
      drop(ia) = .false.
      call trc_ecp_without_atoms(ecp, drop, one)
      if (.not. one%has_ecp) return

      ! The fake shells, after the real ones: per primitive the raised shell,
      ! then (l > 0) the lowered one.
      allocate (up(b%maxnp, nsh), dn(b%maxnp, nsh))
      up = 0
      dn = 0
      nfk = 0
      do ish = 1, nsh
         nper = merge(2, 1, b%sh_l(ish) > 0)
         do ip = 1, b%sh_np(ish)
            up(ip, ish) = nsh + nfk + 1
            if (nper == 2) dn(ip, ish) = nsh + nfk + 2
            nfk = nfk + nper
         end do
      end do
      bf%nshell = nsh + nfk
      bf%natm = b%natm
      bf%maxnp = b%maxnp
      allocate (bf%sh_l(bf%nshell), bf%sh_np(bf%nshell), bf%sh_ao(bf%nshell))
      allocate (bf%sh_e(bf%maxnp, bf%nshell), bf%sh_c(bf%maxnp, bf%nshell), bf%sh_r(3, bf%nshell))
      bf%sh_l(1:nsh) = b%sh_l
      bf%sh_np(1:nsh) = b%sh_np
      bf%sh_ao(1:nsh) = b%sh_ao
      bf%sh_e = 0.0_dp
      bf%sh_c = 0.0_dp
      bf%sh_e(:, 1:nsh) = b%sh_e
      bf%sh_c(:, 1:nsh) = b%sh_c
      bf%sh_r(:, 1:nsh) = b%sh_r
      k = nao
      do ish = 1, nsh
         l = b%sh_l(ish)
         do ip = 1, b%sh_np(ish)
            bf%sh_l(up(ip, ish)) = l + 1
            bf%sh_np(up(ip, ish)) = 1
            bf%sh_ao(up(ip, ish)) = k + 1
            bf%sh_e(1, up(ip, ish)) = b%sh_e(ip, ish)
            bf%sh_c(1, up(ip, ish)) = 2.0_dp*b%sh_e(ip, ish)*b%sh_c(ip, ish)
            bf%sh_r(:, up(ip, ish)) = b%sh_r(:, ish)
            k = k + ncart_l(l + 1)
            if (l == 0) cycle
            bf%sh_l(dn(ip, ish)) = l - 1
            bf%sh_np(dn(ip, ish)) = 1
            bf%sh_ao(dn(ip, ish)) = k + 1
            bf%sh_e(1, dn(ip, ish)) = b%sh_e(ip, ish)
            bf%sh_c(1, dn(ip, ish)) = b%sh_c(ip, ish)
            bf%sh_r(:, dn(ip, ish)) = b%sh_r(:, ish)
            k = k + ncart_l(l - 1)
         end do
      end do
      bf%nao = k
      allocate (bf%at_r(3, b%natm))
      bf%at_r = b%at_r

      call trc_ecp_matrix_cart(bf, one, mf, error, split=nsh)
      if (error%has_error()) then
         call bf%release()
         call one%release()
         return
      end if

      ! Fold: function (lx, ly, lz) of a real shell is slot ifunc of its shell in
      ! this module's lexicographic order, ly and lz from CART_POW_*, and the
      ! slot of a shell one l up or down is (ly + lz)(ly + lz + 1)/2 + lz.
      do ish = 1, nsh
         l = b%sh_l(ish)
         do ifunc = 0, ncart_l(l) - 1
            ly = CART_POW_Y(ifunc)
            lz = CART_POW_Z(ifunc)
            lx = l - ly - lz
            row = b%sh_ao(ish) + cart_perm_l2b(l, ifunc)
            do ip = 1, b%sh_np(ish)
               do c = 1, 3
                  select case (c)
                  case (1)
                     k = lexi(ly, lz)
                  case (2)
                     k = lexi(ly + 1, lz)
                  case default
                     k = lexi(ly, lz + 1)
                  end select
                  ! The rectangle's row 1 is the first AO after the real ones.
                  d(row, :, c) = d(row, :, c) + mf(bf%sh_ao(up(ip, ish)) + cart_perm_l2b(l + 1, k) - nao, :)
                  if (l == 0) cycle
                  select case (c)
                  case (1)
                     if (lx == 0) cycle
                     k = lexi(ly, lz)
                     rowd = lx
                  case (2)
                     if (ly == 0) cycle
                     k = lexi(ly - 1, lz)
                     rowd = ly
                  case default
                     if (lz == 0) cycle
                     k = lexi(ly, lz - 1)
                     rowd = lz
                  end select
                  d(row, :, c) = d(row, :, c) - real(rowd, dp)* &
                                 mf(bf%sh_ao(dn(ip, ish)) + cart_perm_l2b(l - 1, k) - nao, :)
               end do
            end do
         end do
      end do
      call bf%release()
      call one%release()
      deallocate (mf)
   end subroutine trc_ecp_deriv_cart

   pure integer function lexi(ly, lz) result(n)
      !! The slot of the function with powers (ly, lz) in a shell's lexicographic
      !! order, whatever the shell's l (CART_POW_Y/Z's order).
      integer, intent(in) :: ly, lz
      n = (ly + lz)*(ly + lz + 1)/2 + lz
   end function lexi

   subroutine trc_ecp_gradient(b, ecp, density, grad, error, screen, stats)
      !!
      !! ADDS d/dR_X sum_{mu nu} density(mu, nu) < mu | U | nu > into `grad(3, natm)`,
      !! the ECP's own term of the analytic gradient, `density` the symmetric
      !! Cartesian density the one-electron operators contract (D, or D + P for
      !! RI-MP2). One pass over the kept (derivative bra shell, ket shell, ECP
      !! atom) tasks, the density contracted inside each task (this file's
      !! header, THE GRADIENT).
      !!
      !! Orbital shells up to f; `trc_ecp_gradient_refusal` says what is not taken.
      !!
      !! `density` is used on the device: if the caller has it mapped already,
      !! that copy is read and must be current; otherwise it is copied in for
      !! the call. `grad` is host memory. `screen=.false.` runs every task and
      !! `stats` returns the task and chunk counts, for the tests.
      !! TERCO_ECP_GRAD_ATOM_CHUNK=<k> caps a chunk at k ECP atoms (TESTING:
      !! forces the chunked path on a small molecule).
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), intent(in) :: density(b%nao, b%nao)
      real(dp), intent(inout) :: grad(3, b%natm)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: screen
      type(trc_ecp_screen_stats_t), intent(out), optional :: stats

      type(trc_ecp_deriv_basis_t) :: db
      type(trc_ecp_screen_stats_t) :: st
      character(len=:), allocatable :: refusal
      logical :: do_screen

      do_screen = .true.
      if (present(screen)) do_screen = screen
      if (.not. ecp%has_ecp .or. ecp%nrow == 0 .or. b%nshell == 0) return
      ! Checked before anything indexes by atom: a foreign ECP, or a shell on
      ! no atom (a probe function), would write outside `grad` below.
      call ecp_check_rows(b%natm, ecp, error)
      if (error%has_error()) return
      if (any(b%sh_at(1:b%nshell) < 1) .or. any(b%sh_at(1:b%nshell) > b%natm)) then
         call error%set(ERROR_VALIDATION, "trc_ecp_gradient: every shell must sit on an atom of the molecule")
         return
      end if
      refusal = trc_ecp_gradient_refusal(b, ecp)
      if (len(refusal) > 0) then
         call error%set(ERROR_VALIDATION, "trc_ecp_gradient: "//refusal)
         return
      end if
      if (maxval(b%sh_np(1:b%nshell)) > NPMAX) then
         call error%set(ERROR_VALIDATION, "trc_ecp_gradient: an orbital shell has more "// &
                        "primitives than this port's fixed NPMAX ceiling")
         return
      end if

      call ecp_derivative_basis(b, db)
      call ecp_gradient_run(db, ecp, density, grad, do_screen, st, error)
      if (present(stats)) stats = st
   end subroutine trc_ecp_gradient

   function trc_ecp_gradient_refusal(b, ecp) result(refusal)
      !! Why `trc_ecp_gradient` cannot take this basis and ECP; zero-length
      !! when it can. The C API asks before any SCF runs.
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      character(len=:), allocatable :: refusal
      integer :: lmax

      refusal = ""
      if (.not. ecp%has_ecp .or. ecp%nrow == 0 .or. b%nshell == 0) return
      lmax = maxval(b%sh_l(1:b%nshell))
      if (lmax > ECP_ORB_LMAX) then
         refusal = "the ECP gradient does not cover g or higher basis functions (the ECP derivative raises "// &
                   "each shell by one angular momentum, past the ECP kernels' g ceiling)"
      else if (lmax + 1 + maxval(ecp%ecp_l(1:ecp%nrow)) > ECP_ANG_LMAX) then
         refusal = "the ECP gradient does not cover these basis functions with a projected ECP channel of l = "// &
                   to_string(maxval(ecp%ecp_l(1:ecp%nrow)))//": the raised shell and the channel together "// &
                   "pass the angular tables' l = "//to_string(ECP_ANG_LMAX)
      end if
   end function trc_ecp_gradient_refusal

   subroutine ecp_derivative_basis(b, db)
      !! The derivative basis of `b` and its fold arrays (`trc_ecp_deriv_basis_t`).
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_deriv_basis_t), intent(out) :: db

      integer :: ish, s, f, l, lf, m, fx, fy, fz, c, kind
      integer :: p(3)

      db%natm = b%natm
      db%nao = b%nao
      db%maxnp = b%maxnp
      allocate (db%at_r(3, b%natm))
      db%at_r = b%at_r(:, 1:b%natm)
      db%nreal = b%nshell
      db%nshell = b%nshell
      do ish = 1, b%nshell
         db%nshell = db%nshell + merge(2, 1, b%sh_l(ish) > 0)
      end do
      allocate (db%l(db%nshell), db%np(db%nshell))
      allocate (db%e(b%maxnp, db%nshell), db%c(b%maxnp, db%nshell), db%r(3, db%nshell))
      allocate (db%atom(db%nshell - db%nreal))
      db%nfold = ncart_l(max(maxval(b%sh_l(1:b%nshell)) + 1, ECP_ORB_LMAX))
      allocate (db%parent(3, db%nfold, db%nshell - db%nreal), db%coef(3, db%nfold, db%nshell - db%nreal))
      allocate (db%ket_ao(NFMAX, db%nreal))
      db%l(1:db%nreal) = b%sh_l(1:db%nreal)
      db%np(1:db%nreal) = b%sh_np(1:db%nreal)
      db%e(:, 1:db%nreal) = b%sh_e(:, 1:db%nreal)
      db%c(:, 1:db%nreal) = b%sh_c(:, 1:db%nreal)
      db%r(:, 1:db%nreal) = b%sh_r(:, 1:db%nreal)
      db%parent = 1
      db%coef = 0.0_dp
      db%ket_ao = 1

      s = db%nreal
      do ish = 1, db%nreal
         l = b%sh_l(ish)
         do m = 0, ncart_l(l) - 1
            db%ket_ao(m + 1, ish) = b%sh_ao(ish) + cart_perm_l2b(l, m)
         end do
         ! kind +1: the raised shell; -1: the lowered one.
         do kind = 1, -1, -2
            lf = l + kind
            if (lf < 0) cycle
            s = s + 1
            f = s - db%nreal
            db%l(s) = lf
            db%np(s) = b%sh_np(ish)
            db%e(:, s) = b%sh_e(:, ish)
            db%r(:, s) = b%sh_r(:, ish)
            db%atom(f) = b%sh_at(ish)
            if (kind == 1) then
               db%c(:, s) = 2.0_dp*b%sh_e(:, ish)*b%sh_c(:, ish)
            else
               db%c(:, s) = b%sh_c(:, ish)
            end if
            ! Function m of this shell, powers (fx, fy, fz), is the derivative
            ! term of the parent function with one c power less (raised) or
            ! more (lowered).
            do m = 0, ncart_l(lf) - 1
               fy = CART_POW_Y(m)
               fz = CART_POW_Z(m)
               fx = lf - fy - fz
               do c = 1, 3
                  p = [fx, fy, fz]
                  p(c) = p(c) - kind
                  if (p(c) < 0) cycle
                  db%parent(c, m + 1, f) = b%sh_ao(ish) + cart_perm_l2b(l, lexi(p(2), p(3)))
                  db%coef(c, m + 1, f) = merge(1.0_dp, -real(p(c), dp), kind == 1)
               end do
            end do
         end do
      end do
   end subroutine ecp_derivative_basis

   subroutine ecp_gradient_run(db, ecp, density, grad, do_screen, st, error)
      !! The device side of `trc_ecp_gradient`: K1, the screen, the task list,
      !! then K1b and the task kernel for each chunk of ECP atoms (one chunk
      !! unless the caches do not fit), and the scatter into `grad`.
      type(trc_ecp_deriv_basis_t), intent(in) :: db
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), intent(in) :: density(db%nao, db%nao)
      real(dp), intent(inout) :: grad(3, db%natm)
      logical, intent(in) :: do_screen
      type(trc_ecp_screen_stats_t), intent(inout) :: st
      type(error_t), intent(inout) :: error

      type(trc_ecp_channels_t) :: ch
      type(trc_ecp_grad_tasks_t) :: tk
      real(dp), allocatable :: sval(:, :)
      integer :: atoms(2)
      integer :: ns, ngroup, maxnp, natm, i, g, t, nfail
      integer(int64) :: bytes

      call build_groups(ecp, ch%grp_atom, ch%grp_l, ch%grp_row0, ch%grp_row1, ch%ngroup)
      if (ch%ngroup == 0) return
      call build_atom_list(ch%grp_atom, ch%ngroup, ch%atm_ecp, ch%atm_grp0, ch%atm_grp1, ch%natm_ecp)
      ns = db%nshell
      ngroup = ch%ngroup
      maxnp = db%maxnp
      natm = db%natm

      call ecp_quadrature_gauss_chebyshev(ch%rs, ch%ws, ECP_NRS)
      allocate (ch%ur(0:ECP_NRS - 1, ngroup), sval(ns, ngroup))
      !$acc enter data copyin(ch, ch%rs, ch%ws) create(ch%ur)
      call ecp_potential_cache(ecp, ngroup, ch%grp_row0, ch%grp_row1, ch%rs, ch%ur)
      !$acc enter data copyin(ch%grp_atom, ch%grp_l, ch%atm_ecp, ch%atm_grp0, ch%atm_grp1) &
      !$acc            copyin(db, db%l, db%np, db%e, db%c, db%r, db%at_r, db%parent, db%coef, db%ket_ao) &
      !$acc            copyin(density)

      sval = 0.0_dp
      if (do_screen) then
         !$acc enter data create(sval)
         do concurrent(t=1:ns*ngroup) local(i, g)
            g = (t - 1)/ns + 1
            i = t - (g - 1)*ns
            call ecp_screen_item(TRC_DC_INDEX(i), TRC_DC_INDEX(g), ns, db%l, db%np, maxnp, db%e, db%c, &
                                 db%r, natm, db%at_r, ngroup, ch%grp_atom, ch%rs, ch%ws, ch%ur(:, g), sval(i, g))
         end do
         !$acc update self(sval)
         !$acc exit data delete(sval)
      end if

      call ecp_grad_task_list(db, ch, sval, do_screen, tk, error)

      if (.not. error%has_error()) then
         st%ntask = int(ns - db%nreal, int64)*int(db%nreal, int64)*int(ch%natm_ecp, int64)
         st%nkept = tk%ntask
         st%nkept_g = tk%ng
         nfail = 0
         if (tk%ntask > 0) then
            !$acc enter data copyin(tk, tk%bra, tk%ket, tk%atm)
            atoms = [1, 1]
            do while (atoms(1) <= ch%natm_ecp)
               call ecp_grad_chunk_end(tk, ch%natm_ecp, atoms(1), atoms(2), error)
               if (error%has_error()) exit
               call ecp_grad_chunk(db, ch, tk, atoms, density, grad, nfail)
               st%nchunk = st%nchunk + 1
               atoms(1) = atoms(2) + 1
            end do
            !$acc exit data delete(tk%bra, tk%ket, tk%atm, tk)
         end if
      end if

      !$acc exit data delete(ch%grp_atom, ch%grp_l, ch%atm_ecp, ch%atm_grp0, ch%atm_grp1, density) &
      !$acc           delete(db%l, db%np, db%e, db%c, db%r, db%at_r, db%parent, db%coef, db%ket_ao, db) &
      !$acc           delete(ch%rs, ch%ws, ch%ur, ch)
      if (error%has_error()) return

      bytes = sum(int(tk%ncol, int64))*BRA_COL_BYTES + sum(int(tk%ncol_g, int64))*ECP_GBRA_COL_BYTES
      call global_logger%performance("   ECP gradient: kept "//to_string(st%nkept)//" of "//to_string(st%ntask)// &
                                     " (derivative shell, shell, ECP atom) tasks ("//to_string(st%nkept_g)// &
                                     " with a g bra), "//to_string(st%nchunk)//" chunk(s), "// &
                                     to_string(bytes/1000000_int64)//" MB of K1b caches")
      if (nfail > 0) then
         call global_logger%info("   ECP radial quadrature: "//to_string(nfail)// &
                                 " gradient tasks kept the finest level")
      end if
   end subroutine ecp_gradient_run

   subroutine ecp_grad_task_list(db, ch, sval, do_screen, tk, error)
      !! The kept tasks, the f-bra ones before the g-bra ones, each ECP atom by
      !! ECP atom and, within an atom, grouped by (bra l, ket l) so that
      !! neighbouring items do the same work; and the cache columns each atom's
      !! tasks read. The screen is the 3-centre bound of this file's header,
      !! over derivative bra shells and real ket shells.
      type(trc_ecp_deriv_basis_t), intent(in) :: db
      type(trc_ecp_channels_t), intent(in) :: ch
      real(dp), intent(in) :: sval(db%nshell, ch%ngroup)
      logical, intent(in) :: do_screen
      type(trc_ecp_grad_tasks_t), intent(out) :: tk
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: qmax(:)
      integer(int64) :: n, nlow
      integer :: pass, class, a, i, j, g, g0, g1, li, lj, li0, li1, lmax, ngroup, natm_ecp, ngrp

      ngroup = ch%ngroup
      natm_ecp = ch%natm_ecp
      allocate (qmax(ngroup))
      do g = 1, ngroup
         qmax(g) = maxval(sval(1:db%nreal, g))
      end do
      lmax = maxval(db%l)
      allocate (tk%first(natm_ecp + 1), tk%gfirst(natm_ecp + 1), tk%ncol(natm_ecp), tk%ncol_g(natm_ecp))
      allocate (tk%col(db%nshell, natm_ecp))
      tk%col = -1
      nlow = 0

      ! Pass 1 counts and marks the shells each atom reads; pass 2 fills.
      do pass = 1, 2
         n = 0
         ! Class 1: a bra of f or below; class 2: a g bra.
         do class = 1, 2
            li0 = merge(0, ECP_BRA_LMAX, class == 1)
            li1 = merge(min(lmax, ECP_ORB_LMAX), lmax, class == 1)
            do a = 1, natm_ecp
               g0 = ch%atm_grp0(a)
               g1 = ch%atm_grp1(a)
               if (class == 1) then
                  tk%first(a) = int(min(n, int(huge(tk%ntask), int64)))
               else
                  tk%gfirst(a) = int(min(n, int(huge(tk%ntask), int64)))
               end if
               do li = li0, li1
                  do lj = 0, lmax
                     do i = db%nreal + 1, db%nshell
                        if (db%l(i) /= li) cycle
                        ! On the ECP atom itself: w to the atom and -w to it again.
                        if (db%atom(i - db%nreal) == ch%atm_ecp(a)) cycle
                        ! No ket shell can pass if the largest cannot.
                        if (do_screen) then
                           if (sum(sval(i, g0:g1)*qmax(g0:g1)) < ECP_SCREEN_THRESH) cycle
                        end if
                        do j = 1, db%nreal
                           if (db%l(j) /= lj) cycle
                           if (do_screen) then
                              if (pair_bound(sval, db%nshell, ngroup, i, j, g0, g1) < ECP_SCREEN_THRESH) cycle
                           end if
                           n = n + 1
                           if (pass == 1) then
                              tk%col(i, a) = 0
                              tk%col(j, a) = 0
                           else
                              tk%bra(n) = i
                              tk%ket(n) = j
                              tk%atm(n) = a
                           end if
                        end do
                     end do
                  end do
               end do
            end do
            if (class == 1) nlow = n
         end do
         if (pass == 1) then
            if (n > int(huge(tk%ntask), int64)) then
               call error%set(ERROR_VALIDATION, "trc_ecp_gradient: "//to_string(n)//" (derivative shell, shell, "// &
                              "ECP atom) tasks survive the screen, more than a default integer can index")
               return
            end if
            tk%ntask = int(n)
            tk%ng = int(n - nlow)
            allocate (tk%bra(tk%ntask), tk%ket(tk%ntask), tk%atm(tk%ntask))
         end if
      end do
      tk%first(natm_ecp + 1) = tk%ntask - tk%ng
      tk%gfirst(natm_ecp + 1) = tk%ntask

      do a = 1, natm_ecp
         ngrp = ch%atm_grp1(a) - ch%atm_grp0(a) + 1
         tk%ncol(a) = 0
         tk%ncol_g(a) = 0
         do i = 1, db%nshell
            if (tk%col(i, a) < 0) cycle
            if (db%l(i) > ECP_ORB_LMAX) then
               tk%col(i, a) = tk%ncol_g(a)
               tk%ncol_g(a) = tk%ncol_g(a) + ngrp
            else
               tk%col(i, a) = tk%ncol(a)
               tk%ncol(a) = tk%ncol(a) + ngrp
            end if
         end do
      end do
   end subroutine ecp_grad_task_list

   subroutine ecp_grad_chunk_end(tk, natm_ecp, a0, a1, error)
      !! The last ECP atom of the chunk starting at `a0`: as many atoms as the
      !! device's free memory (less GRAD_RESERVE_BYTES) holds the K1b caches
      !! and task slots of, and at least one; all of them on the host build
      !! and in the usual case. TERCO_ECP_GRAD_ATOM_CHUNK=<k> caps the chunk at k atoms.
      type(trc_ecp_grad_tasks_t), intent(in) :: tk
      integer, intent(in) :: natm_ecp, a0
      integer, intent(out) :: a1
      type(error_t), intent(inout) :: error

      integer(int64) :: free, need, more, reserve
      character(len=16) :: v
      integer :: st, cap

      cap = 0
      call get_environment_variable("TERCO_ECP_GRAD_ATOM_CHUNK", v, status=st)
      if (st == 0) then
         read (v, *, iostat=st) cap
         if (st /= 0) cap = 0
      end if

      a1 = a0
      free = device_memory_free_bytes()
      reserve = merge(GRAD_RESERVE_BRA_G_BYTES, GRAD_RESERVE_BYTES, tk%ng > 0)
      need = ecp_grad_atom_bytes(tk, a0)
      if (.not. device_memory_fits(error, "trc_ecp_gradient: the caches of one ECP atom, with the reserve,", &
                                   need + GRAD_MIN_FREE_BYTES, "use a smaller basis or a device with more memory")) return
      do while (a1 < natm_ecp)
         if (cap >= 1 .and. a1 - a0 + 1 >= cap) exit
         more = ecp_grad_atom_bytes(tk, a1 + 1)
         if (free >= 0 .and. need + more + reserve > free) exit
         need = need + more
         a1 = a1 + 1
      end do
   end subroutine ecp_grad_chunk_end

   pure integer(int64) function ecp_grad_atom_bytes(tk, a) result(bytes)
      !! Device memory of ECP atom `a`'s share of a chunk: its cache columns
      !! (and their two index words) and its tasks' slots (w, two cache
      !! columns, the convergence flag).
      type(trc_ecp_grad_tasks_t), intent(in) :: tk
      integer, intent(in) :: a

      bytes = int(tk%ncol(a), int64)*(BRA_COL_BYTES + 8_int64) &
              + int(tk%ncol_g(a), int64)*(ECP_GBRA_COL_BYTES + 8_int64) &
              + int(tk%first(a + 1) - tk%first(a) + tk%gfirst(a + 1) - tk%gfirst(a), int64) &
              *(3*8_int64 + 3*4_int64)
   end function ecp_grad_atom_bytes

   subroutine ecp_grad_chunk(db, ch, tk, atoms, density, grad, nfail)
      !! ECP atoms atoms(1)..atoms(2): K1b over the (shell, group) columns
      !! their tasks read, the task kernel, and the scatter of its slots into
      !! `grad`; both kernels a second time, from trc_ecp_gbra, for the g
      !! shells and their tasks when there are any. `db`, `ch`, `tk`'s task
      !! arrays and `density` are mapped by `ecp_gradient_run`.
      type(trc_ecp_deriv_basis_t), intent(in) :: db
      type(trc_ecp_channels_t), intent(in) :: ch
      type(trc_ecp_grad_tasks_t), intent(in) :: tk
      integer, intent(in) :: atoms(2)
      real(dp), intent(in) :: density(db%nao, db%nao)
      real(dp), intent(inout) :: grad(3, db%natm)
      integer, intent(inout) :: nfail

      type(trc_ecp_k1b_cache_t) :: kc, gc
      type(trc_ecp_chunk_t) :: cn
      integer, allocatable :: col_sh(:), col_g(:)
      integer :: a, i, j, g, k, t, nt, ntg, t0, tg0, ncol, ncolg, base, gbase, iatom
      integer :: nreal, ns, maxnp, natm, nao, ngroup, nfold

      nreal = db%nreal
      ns = db%nshell
      maxnp = db%maxnp
      natm = db%natm
      nao = db%nao
      nfold = db%nfold
      ngroup = ch%ngroup
      ! The chunk's slots: `trc_ecp_chunk_t`.
      t0 = tk%first(atoms(1))
      nt = tk%first(atoms(2) + 1) - t0
      tg0 = tk%gfirst(atoms(1))
      ntg = tk%gfirst(atoms(2) + 1) - tg0
      if (nt + ntg == 0) return
      ncol = sum(tk%ncol(atoms(1):atoms(2)))
      ncolg = sum(tk%ncol_g(atoms(1):atoms(2)))
      cn%nt = nt
      cn%ntg = ntg
      cn%t0 = t0
      cn%tg0 = tg0
      kc%ncol = ncol
      gc%ncol = ncolg

      ! An atom's columns follow those of the chunk's earlier atoms; a shell's
      ! channel groups are consecutive.
      allocate (col_sh(ncol), col_g(ncol), cn%gcol_sh(ncolg), cn%gcol_g(ncolg), cn%ci(nt + ntg), cn%cj(nt + ntg))
      base = 0
      gbase = 0
      do a = atoms(1), atoms(2)
         do i = 1, ns
            if (tk%col(i, a) < 0) cycle
            do g = ch%atm_grp0(a), ch%atm_grp1(a)
               if (db%l(i) > ECP_ORB_LMAX) then
                  k = gbase + tk%col(i, a) + g - ch%atm_grp0(a) + 1
                  cn%gcol_sh(k) = i
                  cn%gcol_g(k) = g
               else
                  k = base + tk%col(i, a) + g - ch%atm_grp0(a) + 1
                  col_sh(k) = i
                  col_g(k) = g
               end if
            end do
         end do
         do t = tk%first(a) + 1, tk%first(a + 1)
            cn%ci(t - t0) = base + tk%col(tk%bra(t), a) + 1
            cn%cj(t - t0) = base + tk%col(tk%ket(t), a) + 1
         end do
         do t = tk%gfirst(a) + 1, tk%gfirst(a + 1)
            cn%ci(nt + t - tg0) = gbase + tk%col(tk%bra(t), a) + 1
            cn%cj(nt + t - tg0) = base + tk%col(tk%ket(t), a) + 1
         end do
         base = base + tk%ncol(a)
         gbase = gbase + tk%ncol_g(a)
      end do

      allocate (kc%fac1(0:FAC1SIZE_BRA - 1, ncol), kc%fac2(0:FAC2SIZE_BRA - 1, ncol), &
                kc%ang2(0:ANG2SIZE_BRA - 1, ncol), cn%w(3, nt + ntg), cn%conv_fail(nt + ntg))
      !$acc enter data copyin(col_sh, col_g, kc, cn)
      !$acc enter data copyin(cn%ci, cn%cj) create(kc%fac1, kc%fac2, kc%ang2, cn%w, cn%conv_fail)

      ! K1b, over the columns that are read only.
      do concurrent(k=1:ncol) local(i, g)
         i = col_sh(k)
         g = col_g(k)
         call ecp_k1b_item(TRC_DC_INDEX(i), TRC_DC_INDEX(g), ns, db%l, db%np, maxnp, db%e, db%c, db%r, &
                           natm, db%at_r, ngroup, ch%grp_atom, ch%grp_l, ch%rs, &
                           kc%fac1(:, k), kc%fac2(:, k), kc%ang2(:, k))
      end do

      ! The task kernel: one task, one block, one contraction, three numbers.
      do concurrent(t=1:nt) local(k, i, j, a)
         k = t0 + t
         i = tk%bra(k)
         j = tk%ket(k)
         a = tk%atm(k)
         call ecp_grad_item(cn%w(:, t), cn%conv_fail(t), TRC_DC_INDEX(i), TRC_DC_INDEX(j), nreal, ch%atm_ecp(a), &
                            ch%atm_grp0(a), ch%atm_grp1(a), cn%ci(t), cn%cj(t), ncol, ncol, ngroup, ch%grp_l, ns, &
                            db%l, db%np, maxnp, db%e, db%c, db%r, natm, db%at_r, ch%rs, ch%ws, ch%ur, &
                            kc%fac1, kc%fac2, kc%ang2, kc%fac1, kc%fac2, kc%ang2, &
                            nao, density, nfold, db%parent, db%coef, NFMAX, db%ket_ao)
      end do

      if (ntg > 0) then
         allocate (gc%fac1(0:ECP_GBRA_FAC1 - 1, ncolg), gc%fac2(0:ECP_GBRA_FAC2 - 1, ncolg), &
                   gc%ang2(0:ECP_GBRA_ANG2 - 1, ncolg))
         !$acc enter data copyin(gc)
         !$acc enter data copyin(cn%gcol_sh, cn%gcol_g) create(gc%fac1, gc%fac2, gc%ang2)
         call trc_ecp_gbra_caches(db, ch, cn, gc)
         call trc_ecp_gbra_tasks(db, ch, tk, density, gc, kc, cn)
         !$acc exit data delete(cn%gcol_sh, cn%gcol_g, gc%fac1, gc%fac2, gc%ang2)
         !$acc exit data delete(gc)
      end if

      !$acc update self(cn%w, cn%conv_fail)
      !$acc exit data delete(cn%ci, cn%cj, kc%fac1, kc%fac2, kc%ang2, cn%w, cn%conv_fail)
      !$acc exit data delete(col_sh, col_g, kc, cn)

      ! The bra shell's atom takes w, the ECP's atom -w (translational
      ! invariance).
      do t = 1, nt + ntg
         k = merge(t0 + t, tg0 + t - nt, t <= nt)
         iatom = db%atom(tk%bra(k) - nreal)
         grad(:, iatom) = grad(:, iatom) + cn%w(:, t)
         iatom = ch%atm_ecp(tk%atm(k))
         grad(:, iatom) = grad(:, iatom) - cn%w(:, t)
      end do
      nfail = nfail + sum(cn%conv_fail)
   end subroutine ecp_grad_chunk

   subroutine ecp_matrix_cart(b, ecp, mcart, error, do_screen, want_stats, st, split)
      !!
      !! The Cartesian build: groups the ECP rows into (atom, l) channels, then
      !! K1/K2/K3 as this module's own header describes.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), intent(inout) :: mcart(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in) :: do_screen, want_stats
      type(trc_ecp_screen_stats_t), intent(out) :: st
      integer, intent(in) :: split
      !! 0: every pair i >= j into the square matrix; otherwise only
      !! i > split >= j, into the rectangle `trc_ecp_matrix_cart` describes

      real(dp) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), allocatable :: sval(:, :)
      integer, allocatable :: kt_p(:), kt_a(:)
      real(dp) :: bound, amax
      integer(int64) :: nk64, ntot64
      integer :: nkept, kk
      real(dp), allocatable :: ur_cache(:, :)
      real(dp), allocatable :: fac1_cache(:, :, :), fac2_cache(:, :, :), ang2_cache(:, :, :)
      integer, allocatable :: grp_atom(:), grp_l(:), grp_row0(:), grp_row1(:)
      integer, allocatable :: atm_ecp(:), atm_grp0(:), atm_grp1(:)
      integer, allocatable :: pr_i(:), pr_j(:)
      real(dp), allocatable :: block_buffer(:, :)
      integer, allocatable :: conv_fail(:)
      integer :: ngroup, natm_ecp, npair, nshell
      integer :: i, j, p, g, a, nfail, t, ntask1b

      call build_groups(ecp, grp_atom, grp_l, grp_row0, grp_row1, ngroup)
      if (ngroup == 0) return
      call build_atom_list(grp_atom, ngroup, atm_ecp, atm_grp0, atm_grp1, natm_ecp)
      nshell = b%nshell

      call ecp_quadrature_gauss_chebyshev(rs, ws, ECP_NRS)
      allocate (ur_cache(0:ECP_NRS - 1, ngroup))

      !$acc enter data copyin(rs) create(ur_cache)
      call ecp_potential_cache(ecp, ngroup, grp_row0, grp_row1, rs, ur_cache)

      npair = nshell*(nshell + 1)/2
      if (split > 0) npair = (nshell - split)*split
      allocate (pr_i(npair), pr_j(npair))
      p = 0
      do i = 1, nshell
         do j = 1, i
            if (split > 0 .and. (i <= split .or. j > split)) cycle
            p = p + 1
            pr_i(p) = i
            pr_j(p) = j
         end do
      end do

      ! K1b: once per (orbital shell, ECP channel group) -- not once per
      ! shell pair -- the semi-local radial/angular factors (fac2_cache,
      ! ang2_cache) and the local term's static centre-shift factors
      ! (fac1_cache). Flattened for the same reason K2 below is: `ngroup`
      ! is as small a run-time extent as `natm_ecp`.
      allocate (fac1_cache(0:FAC1SIZE_BRA - 1, nshell, ngroup))
      allocate (fac2_cache(0:FAC2SIZE_BRA - 1, nshell, ngroup))
      allocate (ang2_cache(0:ANG2SIZE_BRA - 1, nshell, ngroup))
      allocate (sval(nshell, ngroup))

      !$acc enter data copyin(ws, b, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%at_r, grp_atom, grp_l) &
      !$acc            create(fac1_cache, fac2_cache, ang2_cache, sval)
      ntask1b = nshell*ngroup
      do concurrent(t=1:ntask1b) local(i, g)
         g = (t - 1)/nshell + 1
         i = t - (g - 1)*nshell
         call ecp_k1b_item(TRC_DC_INDEX(i), TRC_DC_INDEX(g), nshell, b%sh_l, b%sh_np, b%maxnp, b%sh_e, b%sh_c, b%sh_r, &
                           b%natm, b%at_r, ngroup, grp_atom, grp_l, rs, &
                           fac1_cache(:, i, g), fac2_cache(:, i, g), ang2_cache(:, i, g))
      end do

      ! Screen: one scalar per (orbital shell, channel group), then the
      ! per-task bound on the host from those (a multiply per group), kept
      ! tasks compacted a-major and p-minor, the order the atoms are
      ! summed in. The kernel runs only if the screen or the statistics
      ! read the bound.
      if (do_screen .or. want_stats) then
         do concurrent(t=1:ntask1b) local(i, g)
            g = (t - 1)/nshell + 1
            i = t - (g - 1)*nshell
            call ecp_screen_item(TRC_DC_INDEX(i), TRC_DC_INDEX(g), nshell, b%sh_l, b%sh_np, b%maxnp, b%sh_e, b%sh_c, &
                                 b%sh_r, b%natm, b%at_r, ngroup, grp_atom, rs, ws, ur_cache(:, g), sval(i, g))
         end do
         !$acc update self(sval)
      end if
      !$acc exit data delete(sval)

      ntot64 = int(npair, int64)*int(natm_ecp, int64)
      nk64 = 0
      if (do_screen) then
         do a = 1, natm_ecp
            do p = 1, npair
               bound = pair_bound(sval, nshell, ngroup, pr_i(p), pr_j(p), atm_grp0(a), atm_grp1(a))
               ! `.not. (<)`, so a NaN bound keeps its task.
               if (.not. (bound < ECP_SCREEN_THRESH)) nk64 = nk64 + 1
            end do
         end do
      else
         nk64 = ntot64
      end if
      ! Not for the derivative's rectangle: that is one call per ECP atom on
      ! every gradient.
      if (split == 0) then
         call global_logger%performance("   ECP 3-centre screen: kept "//to_string(nk64)//" of "// &
                                        to_string(ntot64)//" (shell pair, ECP atom) tasks")
      end if
      if (nk64 > int(huge(nkept), int64)) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: "//to_string(nk64)//" (shell pair, ECP atom) "// &
                        "tasks survive the screen, more than a default integer can index")
         !$acc exit data delete(ws, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%at_r, b, &
         !$acc                  grp_l, grp_atom, rs, ur_cache) &
         !$acc            delete(fac1_cache, fac2_cache, ang2_cache)
         return
      end if
      nkept = int(nk64)
      st%ntask = ntot64
      st%nkept = nk64

      allocate (kt_p(nkept), kt_a(nkept))
      kk = 0
      do a = 1, natm_ecp
         do p = 1, npair
            if (do_screen) then
               bound = pair_bound(sval, nshell, ngroup, pr_i(p), pr_j(p), atm_grp0(a), atm_grp1(a))
               if (bound < ECP_SCREEN_THRESH) cycle
            end if
            kk = kk + 1
            if (kk <= nkept) then
               kt_p(kk) = p
               kt_a(kk) = a
            end if
         end do
      end do
      if (kk /= nkept) then
         call error%set(ERROR_GENERIC, "trc_ecp_matrix: the screen's count and fill passes disagree ("// &
                        to_string(nkept)//" counted, "//to_string(kk)//" found); internal error")
         !$acc exit data delete(ws, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%at_r, b, &
         !$acc                  grp_l, grp_atom, rs, ur_cache) &
         !$acc            delete(fac1_cache, fac2_cache, ang2_cache)
         return
      end if

      allocate (block_buffer(SLOT_SIZE, nkept), conv_fail(nkept))
      block_buffer = 0.0_dp
      conv_fail = 0

      ! With nothing kept no kernel writes the buffers, so nothing is
      ! mapped, launched or copied back for K2.
      if (nkept > 0) then
         !$acc enter data copyin(pr_i, pr_j, kt_p, kt_a) copyin(atm_ecp, atm_grp0, atm_grp1) &
         !$acc            create(block_buffer, conv_fail)
         ! One index over the kept list, (kt_p(t), kt_a(t)) giving the pair
         ! and the ECP atom of task t, not do concurrent(p=1:npair,
         ! a=1:natm_ecp): a two-index do concurrent whose outer extent is
         ! small and set at run time (natm_ecp is 1 or 2 for most cases) has
         ! failed the kernel launch outright on the GPU build (nvfortran,
         ! cc70) rather than merely dropping a slice -- gpu-pitfalls.md's
         ! own warning about this shape.
         do concurrent(t=1:nkept) local(p, a)
            p = kt_p(t)
            a = kt_a(t)
            call ecp_pair_atom_item(pr_i(p), pr_j(p), atm_ecp(a), atm_grp0(a), atm_grp1(a), &
                                    ngroup, grp_l, &
                                    nshell, b%sh_l, b%sh_np, b%maxnp, b%sh_e, b%sh_c, b%sh_r, &
                                    b%natm, b%at_r, &
                                    rs, ws, ur_cache, fac1_cache, fac2_cache, ang2_cache, &
                                    block_buffer(:, t), conv_fail(t))
         end do
         !$acc update self(block_buffer, conv_fail)
         !$acc exit data delete(pr_i, pr_j, kt_p, kt_a, atm_ecp, atm_grp0, atm_grp1, block_buffer, conv_fail)
      end if
      !$acc exit data delete(ws, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%at_r, b, &
      !$acc                  grp_l, grp_atom, rs, ur_cache) &
      !$acc            delete(fac1_cache, fac2_cache, ang2_cache)

      if (want_stats) then
         do t = 1, nkept
            amax = maxval(abs(block_buffer(:, t)))
            p = kt_p(t)
            a = kt_a(t)
            bound = pair_bound(sval, nshell, ngroup, pr_i(p), pr_j(p), atm_grp0(a), atm_grp1(a))
            if (bound > 0.0_dp) then
               st%worst_ratio = max(st%worst_ratio, amax/bound)
            else if (amax > 0.0_dp) then
               st%nzero_bound = st%nzero_bound + 1
            end if
         end do
      end if

      ! A block whose levels never agree twice in a row keeps the finest
      ! level, as pyscf's nr_ecp.c does. It is not an error: on I2, Xe2,
      ! CsI, IBr, IO and CH3I in def2-SVP, a shell on one heavy atom against
      ! the other's channel hits this, and the finest level matches pyscf to
      ! 1e-13 (#194 review). Reported, so a change of grid can be traced.
      nfail = sum(conv_fail)
      if (nfail > 0) then
         call global_logger%info("   ECP radial quadrature: "//to_string(nfail)// &
                                 " (shell pair, ECP atom) blocks kept the finest level")
      end if

      ! K3: write each kept slot into the Cartesian block (and
      ! its symmetric mirror -- U is symmetric, scalar ECPs are real, and
      ! only i>=j pairs were built). The derivative's rectangle has no
      ! mirror: its rows start after its ncol ket functions.
      do t = 1, nkept
         if (split == 0) then
            call deposit_slot(b, mcart, pr_i(kt_p(t)), pr_j(kt_p(t)), 0, .true., block_buffer(:, t))
         else
            call deposit_slot(b, mcart, pr_i(kt_p(t)), pr_j(kt_p(t)), size(mcart, 2), .false., block_buffer(:, t))
         end if
      end do

      deallocate (ur_cache, fac1_cache, fac2_cache, ang2_cache, grp_atom, grp_l, grp_row0, grp_row1)
      deallocate (atm_ecp, atm_grp0, atm_grp1, pr_i, pr_j, block_buffer, conv_fail)
      deallocate (kt_p, kt_a, sval)
   end subroutine ecp_matrix_cart

   subroutine ecp_potential_cache(ecp, ngroup, grp_row0, grp_row1, rs, ur_cache)
      !! K1: every channel group's U(r) on the whole grid. `rs` and
      !! `ur_cache` are present on the device already; the ECP rows are
      !! mapped for the kernel only.
      type(trc_ecp_t), intent(in) :: ecp
      integer, intent(in) :: ngroup
      integer, intent(in) :: grp_row0(*), grp_row1(*)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      real(dp), intent(inout) :: ur_cache(0:ECP_NRS - 1, ngroup)
      integer :: g

      !$acc enter data copyin(ecp, ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c) &
      !$acc            copyin(grp_row0(1:ngroup), grp_row1(1:ngroup))
      do concurrent(g=1:ngroup)
         call ecp_k1_item(ur_cache(:, g), rs, grp_row0(g), grp_row1(g), &
                          ecp%nrow, ecp%maxnp, ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c)
      end do
      !$acc exit data delete(ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c, &
      !$acc                  grp_row0(1:ngroup), grp_row1(1:ngroup), ecp)
   end subroutine ecp_potential_cache

   pure real(dp) function pair_bound(q, nshell, ngroup, i, j, g0, g1) result(bound)
      !! The screen's bound on shell pair (i, j) against the ECP atom whose
      !! channel groups are g0..g1: the sum of q(i,g)*q(j,g), `q` being
      !! `ecp_screen_item`'s output.
      integer, intent(in) :: nshell, ngroup, i, j, g0, g1
      real(dp), intent(in) :: q(nshell, ngroup)
      integer :: g

      bound = 0.0_dp
      do g = g0, g1
         bound = bound + q(i, g)*q(j, g)
      end do
   end function pair_bound

   subroutine deposit_slot(b, mcart, i, j, row0, mirror, slot)
      !! Adds one kept (shell pair, ECP atom) slot into the Cartesian matrix:
      !! AO (mu, nu) goes to mcart(mu - row0, nu), and with `mirror` (the
      !! square matrix, row0 = 0) to its transpose too.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(inout) :: mcart(:, :)
      integer, intent(in) :: i, j, row0
      logical, intent(in) :: mirror
      real(dp), intent(in) :: slot(SLOT_SIZE)
      integer :: li, lj, nfi, nfj, ai, aj, mi, mj, pi, pj

      li = b%sh_l(i)
      lj = b%sh_l(j)
      nfi = ncart_l(li)
      nfj = ncart_l(lj)
      ai = b%sh_ao(i)
      aj = b%sh_ao(j)
      do mj = 1, nfj
         do mi = 1, nfi
            ! The slots are lexicographic (this module's angular tables);
            ! the matrix is in this build's order, like S, T and V.
            pi = ai + cart_perm_l2b(li, mi - 1)
            pj = aj + cart_perm_l2b(lj, mj - 1)
            mcart(pi - row0, pj) = mcart(pi - row0, pj) + slot((mj - 1)*NFMAX + mi)
            if (mirror .and. i /= j) mcart(pj, pi) = mcart(pj, pi) + slot((mj - 1)*NFMAX + mi)
         end do
      end do
   end subroutine deposit_slot

   subroutine build_groups(ecp, grp_atom, grp_l, grp_row0, grp_row1, ngroup)
      !!
      !! Group the ECP rows into (atom, l) channels: consecutive rows sharing
      !! both -- see trc_ecp_types.F90 for why the reader always emits a
      !! channel's rows consecutively.
      !!
      type(trc_ecp_t), intent(in) :: ecp
      integer, allocatable, intent(out) :: grp_atom(:), grp_l(:), grp_row0(:), grp_row1(:)
      integer, intent(out) :: ngroup
      integer :: i

      allocate (grp_atom(max(ecp%nrow, 1)), grp_l(max(ecp%nrow, 1)))
      allocate (grp_row0(max(ecp%nrow, 1)), grp_row1(max(ecp%nrow, 1)))
      ngroup = 0
      do i = 1, ecp%nrow
         if (ngroup > 0) then
            if (ecp%ecp_at(i) == grp_atom(ngroup) .and. ecp%ecp_l(i) == grp_l(ngroup)) then
               grp_row1(ngroup) = i
               cycle
            end if
         end if
         ngroup = ngroup + 1
         grp_atom(ngroup) = ecp%ecp_at(i)
         grp_l(ngroup) = ecp%ecp_l(i)
         grp_row0(ngroup) = i
         grp_row1(ngroup) = i
      end do
   end subroutine build_groups

   subroutine build_atom_list(grp_atom, ngroup, atm_ecp, atm_grp0, atm_grp1, natm_ecp)
      integer, intent(in) :: ngroup
      integer, intent(in) :: grp_atom(:)
      integer, allocatable, intent(out) :: atm_ecp(:), atm_grp0(:), atm_grp1(:)
      integer, intent(out) :: natm_ecp
      integer :: g

      allocate (atm_ecp(max(ngroup, 1)), atm_grp0(max(ngroup, 1)), atm_grp1(max(ngroup, 1)))
      natm_ecp = 0
      do g = 1, ngroup
         if (natm_ecp > 0) then
            if (grp_atom(g) == atm_ecp(natm_ecp)) then
               atm_grp1(natm_ecp) = g
               cycle
            end if
         end if
         natm_ecp = natm_ecp + 1
         atm_ecp(natm_ecp) = grp_atom(g)
         atm_grp0(natm_ecp) = g
         atm_grp1(natm_ecp) = g
      end do
   end subroutine build_atom_list

   pure subroutine ecp_screen_item(ish, g, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                   natm, at_r, ngroup, grp_atom, rs, ws, ur, qval)
      !!
      !! The screen's per-(orbital shell, channel group) scalar q(i,g) of this
      !! file's header: sqrt(4*pi * sum over levels of step_L * the level's
      !! quadrature sum of w_k |U_g(r_k)| F_i(r_k)**2). Level `lvl`'s points
      !! are every `step`-th point of the finest grid starting at `step-1`,
      !! each weighted `step*ws`, the weight the level loops of
      !! `ecp_radial_type1_rad_part_levels` and `ecp_radial_type2_lab_levels`
      !! give a level's points. That set is a superset of the points those
      !! loops add at the level (2**level - 1 new points per level); the extra
      !! terms are non-negative, so they only enlarge the bound. The terms are accumulated in units of the largest
      !! sqrt(w|U|)*F seen so far, so a shell whose F is far below 1e-154
      !! still gives a nonzero q instead of squares that underflow.
      !!
      !$acc routine seq
      integer, intent(in) :: ish, g, nshell, maxnp, natm, ngroup
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: at_r(3, natm)
      integer, intent(in) :: grp_atom(ngroup)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(out) :: qval

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: d, r, f, dr, gk, scale, total, acc
      real(dp) :: part(0:ECP_NLEVEL - 1)
      integer :: lvl, k, kk, m, ip, l

      l = sh_l(ish)
      d = sqrt(sum((at_r(:, grp_atom(g)) - sh_r(:, ish))**2))
      ! One pass over the finest grid. Point k first belongs to the level
      ! whose stride is the largest power of two dividing k+1 (the coarsest
      ! level takes every point whose stride is at least its own), and a
      ! level also holds every point of the levels below it.
      part = 0.0_dp
      scale = 0.0_dp
      do k = 0, ECP_NRS - 1
         if (ur(k) == 0.0_dp) cycle
         kk = k + 1
         m = 0
         do while (mod(kk, 2) == 0)
            kk = kk/2
            m = m + 1
         end do
         lvl = max(0, ECP_LEVEL_MAX - ECP_LEVEL0 - m)
         r = rs(k)
         dr = r - d
         f = 0.0_dp
         do ip = 1, sh_np(ish)
            f = f + abs(sh_c(ip, ish))*exp(-sh_e(ip, ish)*dr*dr)
         end do
         gk = sqrt(ws(k)*abs(ur(k)))*f*(r + d)**l
         if (gk > scale .or. gk /= gk) then
            if (scale > 0.0_dp) part = part*(scale/gk)**2
            scale = gk
         end if
         if (scale > 0.0_dp) part(lvl) = part(lvl) + (gk/scale)**2
      end do
      total = 0.0_dp
      acc = 0.0_dp
      do lvl = 0, ECP_NLEVEL - 1
         ! The cumulative level sum, times the level's own weight scale.
         acc = acc + part(lvl)
         total = total + acc*real(2**(ECP_LEVEL_MAX - ECP_LEVEL0 - lvl), dp)
      end do
      qval = scale*sqrt(4.0_dp*PI*total)
   end subroutine ecp_screen_item

   pure subroutine ecp_pair_atom_item(ish, jsh, atm_id, g0, g1, &
                                      ngroup, grp_l, &
                                      nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                      natm, at_r, &
                                      rs, ws, ur_cache, fac1_cache, fac2_cache, ang2_cache, &
                                      slot, conv_fail)
      !!
      !! K2's item: one (orbital shell pair, ECP atom). Sums every one of that
      !! atom's channels (local and projected) into its OWN output slot,
      !! `slot` -- SLOT_SIZE = NFMAX*NFMAX reals, column-major (mi, mj) with mi the
      !! FASTER index, zero-padded past the pair's actual (nfi, nfj). No
      !! atomics: nothing else ever writes this slot.
      !!
      !$acc routine seq
      integer, intent(in) :: ish, jsh, atm_id, g0, g1
      integer, intent(in) :: ngroup
      integer, intent(in) :: grp_l(ngroup)
      integer, intent(in) :: nshell, maxnp
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      !! `maxnp` is `b%maxnp`, the basis's OWN leading dimension for
      !! sh_e/sh_c -- NOT this module's fixed NPMAX. Declaring these two
      !! with NPMAX as their leading dimension while the actual argument
      !! (b%sh_e/b%sh_c) is shaped (b%maxnp, nshell) was a real bug here:
      !! an explicit-shape dummy addresses its argument by the dummy's OWN
      !! declared extents, so a smaller actual leading dimension (every
      !! basis this port's tests use has maxnp < NPMAX) silently reads
      !! every shell past the first at the WRONG stride -- found via
      !! shell 2 of Xe/def2-SVP reading an exponent of exactly 0.0.
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      integer, intent(in) :: natm
      real(dp), intent(in) :: at_r(3, natm)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), intent(in) :: ur_cache(0:ECP_NRS - 1, ngroup)
      !! K1b's own global caches -- see its own header. Column
      !! (g - 1)*nshell + shell: the pair's two shells against the SAME
      !! channel group, never recomputed.
      real(dp), intent(in) :: fac1_cache(0:FAC1SIZE_BRA - 1, nshell*ngroup)
      real(dp), intent(in) :: fac2_cache(0:FAC2SIZE_BRA - 1, nshell*ngroup)
      real(dp), intent(in) :: ang2_cache(0:ANG2SIZE_BRA - 1, nshell*ngroup)
      real(dp), intent(out) :: slot(SLOT_SIZE)
      integer, intent(out) :: conv_fail

      real(dp) :: gctr(NFMAX, NFMAX)
      integer :: mi, mj

      ! Bra and ket have one size here: the same caches serve both.
      call ecp_atom_block(gctr, conv_fail, ish, jsh, atm_id, g0, g1, (g0 - 1)*nshell + ish, (g0 - 1)*nshell + jsh, &
                          nshell, nshell*ngroup, nshell*ngroup, ngroup, grp_l, nshell, sh_l, sh_np, maxnp, &
                          sh_e, sh_c, sh_r, natm, at_r, rs, ws, ur_cache, fac1_cache, fac2_cache, ang2_cache, &
                          fac1_cache, fac2_cache, ang2_cache)

      do mj = 1, NFMAX
         do mi = 1, NFMAX
            slot((mj - 1)*NFMAX + mi) = gctr(mi, mj)
         end do
      end do
   end subroutine ecp_pair_atom_item

end module trc_ecp
