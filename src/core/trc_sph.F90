!
! The spherical (pure) AO space, on top of a basis that stays Cartesian
! underneath.
!
! All integrals are computed in Cartesian AOs (see docs/GENERATORS.md), but
! def2 and most modern bases are defined spherical. This module converts
! (nao, nao) Cartesian matrices to (nao_sph, nao_sph) ones and back, so the
! SCF (trc_scf_driver.F90) runs S, H, D, F, C, eps at nao_sph while the
! integrals underneath are built at nao.
!
! T (nao, nao_sph) is block diagonal: one (ncart(l), 2l+1) block per shell,
! at most (15, 9) through g, read from SPH_T (trc_sph_tables.F90, generated
! by scripts/gen_sph.py). Rather than dense GEMMs that multiply by zero,
! `cart_to_sph`/`sph_to_cart` do one small contraction per shell pair,
! O(nshell^2 * ncart^3), negligible next to the O(nao^4) Fock build.
!
! The transforms run on the device as `do concurrent` over a flattened
! shell-pair index; a 2D header gave wrong results under nvfortran (see
! `expand_kernel` in trc_decontract.F90). SPH_T is a parameter, so it needs
! no `!$acc declare`.
!
! No GAMESS-style rescale is folded in: terco's Cartesian convention already
! matches the solid-harmonic closed form (scripts/gen_sph.py). It does not
! compose with trc_ao_norm.F90's UNIT scale; trc_capi.F90 refuses that combo.
!
! XC: the grid integrator stays Cartesian. The SCF driver builds D_cart from
! D_sph with `sph_to_cart` once per iteration and gives it to both the Fock
! build and XC, so one transform serves both.
!
! TWO ONE-SIDED FOLDS, for a correlated step whose kernels also stay
! Cartesian (RI-MP2, --df): `sph_mo_to_cart` is the same per-shell T, but
! contracted on ONE AO index only -- an MO coefficient matrix's row, its
! column being an orbital label rather than a second AO -- and
! `sph_aux3c_to_sph` folds the auxiliary (P) index of a three-centre
! (mu nu|P) block the other way, Cartesian to spherical, leaving the
! orbital-basis pair (mu, nu) alone. `sph_metric_cart_to_sph` is the
! square two-sided fold `sph_cart_to_sph` already is, specialised to the
! auxiliary basis's own metric (P|Q).
!
! sph_mo_to_cart/sph_cart_to_sph/sph_sph_to_cart run ONCE per geometry (the
! SCF's S/H, or RI-MP2's own MO fold) or once per iteration on matrices no
! larger than nao -- the do-concurrent, one-contraction-per-shell-pair
! cost (O(nshell^2 ncart^3), this module's own header above) is negligible
! next to the Fock build there. `sph_aux3c_to_sph`/`sph_metric_cart_to_sph`
! fold something O(nao^2) PER AUXILIARY FUNCTION (a three-centre block) or
! O(naux^2) (the metric) instead: the SAME per-(mu,nu)-pair / per-row
! shell loop that is cheap at nao^2 pairs is now O(nao^2 * naux) total,
! and measured expensive (trc_rimp2_driver.F90's own RI-MP2 step table,
! w20/cc-pVDZ: 0.88 s on top of a 0.07 s three-centre build). Both fold
! the auxiliary index with a batched, strided GEMM through trc_linalg
! (cuBLAS on the device) instead: shells of the SAME angular momentum l
! share one (ncart(l), 2l+1) T block (SPH_T(:,:,l)) and a RUN of
! consecutive same-l shells (common within one atom's own shell list)
! becomes ONE `gemm_batched` call with that many batch items -- a single
! shell is just a run of length one, so every shell is covered without a
! special case. The small SPH_T table is staged into a device-resident
! scratch ONCE per call (`!$acc update device`), not re-read as a
! compile-time PARAMETER inside `host_data use_device` -- SPH_T itself was
! never `!$acc enter data`-registered (this module's own header: "needs
! no !$acc declare"), so cuBLAS's own device-pointer lookup on it directly
! would not resolve.
!
module trc_sph
   use pic_types, only: dp, int64
   use trc_basis_types, only: trc_basis_t
   use trc_sph_tables, only: SPH_T, SPH_LMAX, SPH_NCART_MAX, SPH_NSPH_MAX
   use trc_linalg, only: trc_linalg_t, trc_gemm_batch_t, trc_gemm_batch
   implicit none
   private

   public :: sph_cart_to_sph, sph_sph_to_cart, sph_mo_to_cart, sph_aux3c_to_sph, sph_metric_cart_to_sph
   public :: sph_cols_sph_to_cart

contains

   !
   ! M_sph := T^T M_cart T. Wrapper over the explicit-shape driver, for a
   ! caller holding a `trc_basis_t` rather than its raw shell arrays --
   ! trc_api.F90's own `trc_1e`/`trc_multipoles` are the same shape of
   ! wrapper for the same reason (explicit extents, declared first, so
   ! nvfortran is not walking an assumed-shape descriptor per launch).
   !
   subroutine sph_cart_to_sph(b, mcart, msph)
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: mcart(b%nao, b%nao)
      real(dp), intent(out) :: msph(b%nao_sph, b%nao_sph)
      call cart_to_sph_driver(b%nshell, b%sh_l, b%sh_ao, b%sh_ao_sph, &
                              b%nao, b%nao_sph, mcart, msph)
   end subroutine sph_cart_to_sph

   subroutine sph_sph_to_cart(b, msph, mcart)
      !! M_cart := T M_sph T^T.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in)  :: msph(b%nao_sph, b%nao_sph)
      real(dp), intent(out) :: mcart(b%nao, b%nao)
      call sph_to_cart_driver(b%nshell, b%sh_l, b%sh_ao, b%sh_ao_sph, &
                              b%nao, b%nao_sph, msph, mcart)
   end subroutine sph_sph_to_cart

   !
   ! M_sph(bi:bi+nsi-1, bj:bj+nsj-1) := Ti^T * M_cart(ai:ai+nci-1, aj:aj+ncj-1) * Tj
   ! for every shell pair (i, j), Ti/Tj the (ncart, nsph) blocks SPH_T(:,:,l)
   ! carries for that shell's own angular momentum. General rather than
   ! symmetric-only: it is run on S, H, D and F alike, and only D and F
   ! happen to be symmetric -- a caller handing it a nonsymmetric matrix
   ! (there is none today) would still get the right answer.
   !
   subroutine cart_to_sph_driver(nshell, sh_l, sh_ao, sh_ao_sph, nao_c, nao_s, mcart, msph)
      integer,  intent(in)  :: nshell, nao_c, nao_s
      integer,  intent(in)  :: sh_l(nshell), sh_ao(nshell), sh_ao_sph(nshell)
      real(dp), intent(in)  :: mcart(nao_c, nao_c)
      real(dp), intent(out) :: msph(nao_s, nao_s)
      integer :: ij, ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj
      integer :: p, q, r, s
      real(dp) :: acc, inner

      ! Flattened shell-pair index -- see this module's own header on why
      ! `do concurrent(jsh=..., ish=...)` is not used here.
      do concurrent(ij=0:nshell*nshell - 1) local(ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj, p, q, r, s, acc, inner)
         jsh = ij/nshell + 1
         ish = ij - (jsh - 1)*nshell + 1
         li = sh_l(ish)
         lj = sh_l(jsh)
         nci = (li + 1)*(li + 2)/2
         ncj = (lj + 1)*(lj + 2)/2
         nsi = 2*li + 1
         nsj = 2*lj + 1
         ai = sh_ao(ish)
         aj = sh_ao(jsh)
         bi = sh_ao_sph(ish)
         bj = sh_ao_sph(jsh)
         do s = 1, nsj
            do r = 1, nsi
               acc = 0.0_dp
               do q = 1, ncj
                  inner = 0.0_dp
                  do p = 1, nci
                     inner = inner + SPH_T(p, r, li)*mcart(ai + p - 1, aj + q - 1)
                  end do
                  acc = acc + inner*SPH_T(q, s, lj)
               end do
               msph(bi + r - 1, bj + s - 1) = acc
            end do
         end do
      end do
   end subroutine cart_to_sph_driver

   !
   ! M_cart(ai:ai+nci-1, aj:aj+ncj-1) := Ti * M_sph(bi:bi+nsi-1, bj:bj+nsj-1) * Tj^T,
   ! the reverse contraction, same shell-pair loop.
   !
   subroutine sph_to_cart_driver(nshell, sh_l, sh_ao, sh_ao_sph, nao_c, nao_s, msph, mcart)
      integer,  intent(in)  :: nshell, nao_c, nao_s
      integer,  intent(in)  :: sh_l(nshell), sh_ao(nshell), sh_ao_sph(nshell)
      real(dp), intent(in)  :: msph(nao_s, nao_s)
      real(dp), intent(out) :: mcart(nao_c, nao_c)
      integer :: ij, ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj
      integer :: p, q, r, s
      real(dp) :: acc, inner

      do concurrent(ij=0:nshell*nshell - 1) local(ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj, p, q, r, s, acc, inner)
         jsh = ij/nshell + 1
         ish = ij - (jsh - 1)*nshell + 1
         li = sh_l(ish)
         lj = sh_l(jsh)
         nci = (li + 1)*(li + 2)/2
         ncj = (lj + 1)*(lj + 2)/2
         nsi = 2*li + 1
         nsj = 2*lj + 1
         ai = sh_ao(ish)
         aj = sh_ao(jsh)
         bi = sh_ao_sph(ish)
         bj = sh_ao_sph(jsh)
         do q = 1, ncj
            do p = 1, nci
               acc = 0.0_dp
               do s = 1, nsj
                  inner = 0.0_dp
                  do r = 1, nsi
                     inner = inner + SPH_T(p, r, li)*msph(bi + r - 1, bj + s - 1)
                  end do
                  acc = acc + inner*SPH_T(q, s, lj)
               end do
               mcart(ai + p - 1, aj + q - 1) = acc
            end do
         end do
      end do
   end subroutine sph_to_cart_driver

   !
   ! Ccart := T Csph, the ONE-SIDED fold a rectangular, AO-row/label-column
   ! matrix needs: only the AO row index crosses the spherical/Cartesian
   ! boundary, same per-shell T block `sph_sph_to_cart` contracts on BOTH
   ! sides of a square two-AO-index matrix with; the column index here is
   ! some other label entirely (an MO number for RI-MP2's own use at
   ! trc_rimp2_driver.F90's `rimp2_driver_mo_to_cart`) and is left alone. `ncol`
   ! need not be `b%nao_sph` -- a caller folding only a slice of columns
   ! passes however many it has.
   !
   subroutine sph_mo_to_cart(b, ncol, csph, ccart)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ncol
      real(dp), intent(in)  :: csph(b%nao_sph, ncol)
      real(dp), intent(out) :: ccart(b%nao, ncol)
      call mo_sph_to_cart_driver(b%nshell, b%sh_l, b%sh_ao, b%sh_ao_sph, b%nao, b%nao_sph, ncol, csph, ccart)
   end subroutine sph_mo_to_cart

   subroutine mo_sph_to_cart_driver(nshell, sh_l, sh_ao, sh_ao_sph, nao_c, nao_s, ncol, csph, ccart)
      integer,  intent(in)  :: nshell, nao_c, nao_s, ncol
      integer,  intent(in)  :: sh_l(nshell), sh_ao(nshell), sh_ao_sph(nshell)
      real(dp), intent(in)  :: csph(nao_s, ncol)
      real(dp), intent(out) :: ccart(nao_c, ncol)
      integer :: ic, ish, col, li, nci, nsi, ai, bi, p, r
      real(dp) :: acc

      ! Flattened (shell, column) index -- this module's own header on why
      ! a 2D `do concurrent` header is not used for a shell-indexed bound.
      do concurrent(ic=0:nshell*ncol - 1) local(ish, col, li, nci, nsi, ai, bi, p, r, acc)
         col = ic/nshell + 1
         ish = ic - (col - 1)*nshell + 1
         li = sh_l(ish)
         nci = (li + 1)*(li + 2)/2
         nsi = 2*li + 1
         ai = sh_ao(ish)
         bi = sh_ao_sph(ish)
         do p = 1, nci
            acc = 0.0_dp
            do r = 1, nsi
               acc = acc + SPH_T(p, r, li)*csph(bi + r - 1, col)
            end do
            ccart(ai + p - 1, col) = acc
         end do
      end do
   end subroutine mo_sph_to_cart_driver

   !
   ! Tsph(:, :, q) := sum_p SPH_T(p, r, l) Tcart(:, :, p) -- the P-index
   ! fold a three-centre density-fitting block needs when the auxiliary
   ! basis is spherical (trc_scf_df.F90's `df_jk_setup`,
   ! trc_rimp2_driver.F90's own 3c build): the first two indices (mu, nu),
   ! an orbital-basis AO pair, are untouched; only the auxiliary index
   ! crosses the boundary, on the LAST axis instead of a row or column.
   ! `aux` is the block's own basis handle -- normally a shell-boundary
   ! subset (trc_api.F90's `basis_subset`), whose `nao`/`nao_sph`/`sh_ao`/
   ! `sh_ao_sph` already describe that block's own, LOCAL P numbering --
   ! but the whole auxiliary basis works just as well for a caller with
   ! no blocking at all. `lp` is the caller's own linear-algebra handles
   ! (cuBLAS under OpenACC): see this module's own header for why this is
   ! a batched GEMM, not a do-concurrent loop, and `lp` is now required.
   !
   subroutine sph_aux3c_to_sph(aux, n, tcart, tsph, lp)
      type(trc_basis_t), intent(in) :: aux
      integer, intent(in) :: n
      real(dp), intent(in)  :: tcart(n, n, aux%nao)
      real(dp), intent(out) :: tsph(n, n, aux%nao_sph)
      type(trc_linalg_t), intent(inout) :: lp
      call sph_fold_cols_cart_to_sph(aux, lp, n*n, tcart, tsph)
   end subroutine sph_aux3c_to_sph

   !
   ! Msph := T^T Mcart T, for the auxiliary basis's own metric (P|Q) --
   ! `sph_cart_to_sph`'s own two-sided contraction, but as two batched-GEMM
   ! one-sided folds (columns, then rows) through `lp` instead of a
   ! do-concurrent shell-pair loop: the metric is O(naux^2), and this
   ! module's own header measured that loop expensive at the auxiliary
   ! sizes RI-MP2/`--df` actually run at. `temporary` (naux_c, naux_sph) is
   ! the column-folded intermediate between the two passes.
   !
   subroutine sph_metric_cart_to_sph(aux, lp, mcart, msph)
      type(trc_basis_t), intent(in) :: aux
      type(trc_linalg_t), intent(inout) :: lp
      real(dp), intent(in)  :: mcart(aux%nao, aux%nao)
      real(dp), intent(out) :: msph(aux%nao_sph, aux%nao_sph)
      real(dp), allocatable :: temporary(:, :)
      allocate (temporary(aux%nao, aux%nao_sph))
      !$acc enter data create(temporary)
      call sph_fold_cols_cart_to_sph(aux, lp, aux%nao, mcart, temporary)
      call sph_fold_rows_cart_to_sph(aux, lp, aux%nao_sph, temporary, msph)
      !$acc exit data delete(temporary)
      deallocate (temporary)
   end subroutine sph_metric_cart_to_sph

   !
   ! Msph(:, Psph) := sum_Pc Mcart(:, Pc) SPH_T(Pc, Psph, l) -- the
   ! one-sided, COLUMN-wise Cartesian-to-spherical fold `sph_aux3c_to_sph`
   ! (nrow = n^2, the flattened orbital-AO-pair axis) and
   ! `sph_metric_cart_to_sph`'s own first pass (nrow = aux%nao) both are.
   ! `mcart`/`msph` are passed by their first element and read via `nrow`
   ! alone (BLAS's own `a(lda,*)` convention, trc_linalg.F90's own idiom
   ! for a sub-block of a larger array) so the actual argument's own shape
   ! -- (n, n, naux) for a 3c block, (naux_c, naux_c) for the metric --
   ! never has to match this routine's flattened (nrow, naux) view of it.
   !
   ! One `gemm_batched` call per RUN of consecutive shells sharing an
   ! angular momentum l (a single shell is a run of length one): every
   ! item of the run shares the SAME (ncart(l), 2l+1) block of `tl`
   ! (`sb = 0`), so strided-batched cuBLAS does every shell of the run in
   ! one launch, at m = nrow, k = ncart(l), n = 2l+1 -- a tall, thin,
   ! memory-bound GEMM, not the O(nrow * naux_c * naux_sph) a single dense
   ! GEMM against a block-diagonal T would spend most of its flops on.
   !
   subroutine sph_fold_cols_cart_to_sph(aux, lp, nrow, mcart, msph)
      type(trc_basis_t), intent(in) :: aux
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nrow
      real(dp), intent(in)  :: mcart(*)
      real(dp), intent(out) :: msph(*)
      real(dp) :: tl(SPH_NCART_MAX, SPH_NSPH_MAX, 0:SPH_LMAX)
      integer :: ish0, ish, l0, nci, nsi, run, ac, bc
      integer(int64) :: nrow8
      type(trc_gemm_batch_t) :: gemm_shape

      nrow8 = int(nrow, int64)
      !$acc enter data create(tl)
      tl = SPH_T
      !$acc update device(tl)

      ish0 = 1
      do while (ish0 <= aux%nshell)
         l0 = aux%sh_l(ish0)
         run = 1
         do ish = ish0 + 1, aux%nshell
            if (aux%sh_l(ish) /= l0) exit
            run = run + 1
         end do
         nci = (l0 + 1)*(l0 + 2)/2
         nsi = 2*l0 + 1
         ac = aux%sh_ao(ish0) - 1
         bc = aux%sh_ao_sph(ish0) - 1
         ! T is shared by every item (sb = 0) and read by its declared
         ! leading dimension; the M tiles are tight.
         gemm_shape = trc_gemm_batch(nrow, nsi, nci, run, ldb=SPH_NCART_MAX, sb=0_int64)
         call lp%gemm_batched(gemm_shape, mcart(1_int64 + int(ac, int64)*nrow8), tl(1, 1, l0), &
                              msph(1_int64 + int(bc, int64)*nrow8))
         ish0 = ish0 + run
      end do
      !$acc exit data delete(tl)
   end subroutine sph_fold_cols_cart_to_sph

   !
   ! Mcart(:, Pc) := sum_Psph Msph(:, Psph) SPH_T(Pc, Psph, l) -- the
   ! adjoint of `sph_fold_cols_cart_to_sph`: a column-indexed quantity whose
   ! columns are spherical AUXILIARY functions, expanded onto the Cartesian
   ! ones the derivative kernels run over. The RI-MP2 gradient
   ! (trc_rimp2_relax.F90) needs it for its three-index density Gamma^P
   ! and, applied twice around a transpose, for its two-index gamma_PQ:
   ! sum_Ps G^Ps (mu nu|Ps) = sum_Pc (T G)^Pc (mu nu|Pc), so contracting
   ! the expanded density with the Cartesian derivative integrals is exact.
   ! Same one-batched-GEMM-per-run-of-equal-l shape as the fold above, with
   ! T now transposed (`transb = "T"`).
   !
   subroutine sph_cols_sph_to_cart(aux, lp, nrow, msph, mcart)
      type(trc_basis_t), intent(in) :: aux
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: nrow
      real(dp), intent(in)  :: msph(*)
      real(dp), intent(out) :: mcart(*)
      real(dp) :: tl(SPH_NCART_MAX, SPH_NSPH_MAX, 0:SPH_LMAX)
      integer :: ish0, ish, l0, nci, nsi, run, ac, bc
      integer(int64) :: nrow8
      type(trc_gemm_batch_t) :: gemm_shape

      nrow8 = int(nrow, int64)
      !$acc enter data create(tl)
      tl = SPH_T
      !$acc update device(tl)

      ish0 = 1
      do while (ish0 <= aux%nshell)
         l0 = aux%sh_l(ish0)
         run = 1
         do ish = ish0 + 1, aux%nshell
            if (aux%sh_l(ish) /= l0) exit
            run = run + 1
         end do
         nci = (l0 + 1)*(l0 + 2)/2
         nsi = 2*l0 + 1
         ac = aux%sh_ao(ish0) - 1
         bc = aux%sh_ao_sph(ish0) - 1
         ! C (nrow x nci) = A (nrow x nsi) T^T, T stored (nci, nsi) in a
         ! (SPH_NCART_MAX, SPH_NSPH_MAX) slot and shared by every item.
         gemm_shape = trc_gemm_batch(nrow, nci, nsi, run, ldb=SPH_NCART_MAX, sb=0_int64, transb="T")
         call lp%gemm_batched(gemm_shape, msph(1_int64 + int(bc, int64)*nrow8), tl(1, 1, l0), &
                              mcart(1_int64 + int(ac, int64)*nrow8))
         ish0 = ish0 + run
      end do
      !$acc exit data delete(tl)
   end subroutine sph_cols_sph_to_cart

   !
   ! Msph(Psph, :) := sum_Pc SPH_T(Pc, Psph, l)^T Mcart(Pc, :) -- the
   ! one-sided, ROW-wise fold `sph_metric_cart_to_sph`'s own second pass
   ! needs (folding the metric's remaining Cartesian index after the
   ! column pass above already folded the other one). `mcart` is
   ! (aux%nao, ncol), `msph` is (aux%nao_sph, ncol), both read/written by
   ! their own leading dimension (`aux%nao`/`aux%nao_sph`) rather than the
   ! `ncol` this routine sees -- a shell's own ROW block is `nci`
   ! contiguous rows at a fixed column offset, not a contiguous run of
   ! elements, so (unlike the column fold) the stride between a run's
   ! items is in ROWS (`nci`/`nsi` elements), and T itself is the "A"
   ! operand here, transposed (`ta = "T"`), not "B".
   !
   subroutine sph_fold_rows_cart_to_sph(aux, lp, ncol, mcart, msph)
      type(trc_basis_t), intent(in) :: aux
      type(trc_linalg_t), intent(inout) :: lp
      integer, intent(in) :: ncol
      real(dp), intent(in)  :: mcart(*)
      real(dp), intent(out) :: msph(*)
      real(dp) :: tl(SPH_NCART_MAX, SPH_NSPH_MAX, 0:SPH_LMAX)
      integer :: ish0, ish, l0, nci, nsi, run, ac, bc, naux_c, naux_s
      type(trc_gemm_batch_t) :: gemm_shape

      naux_c = aux%nao
      naux_s = aux%nao_sph
      !$acc enter data create(tl)
      tl = SPH_T
      !$acc update device(tl)

      ish0 = 1
      do while (ish0 <= aux%nshell)
         l0 = aux%sh_l(ish0)
         run = 1
         do ish = ish0 + 1, aux%nshell
            if (aux%sh_l(ish) /= l0) exit
            run = run + 1
         end do
         nci = (l0 + 1)*(l0 + 2)/2
         nsi = 2*l0 + 1
         ac = aux%sh_ao(ish0) - 1
         bc = aux%sh_ao_sph(ish0) - 1
         ! Items step down the rows of M (nci / nsi elements), each a block
         ! of the full-height M; T is shared.
         gemm_shape = trc_gemm_batch(nsi, ncol, nci, run, lda=SPH_NCART_MAX, ldb=naux_c, ldc=naux_s, &
                                     sa=0_int64, sb=int(nci, int64), sc=int(nsi, int64), transa="T")
         call lp%gemm_batched(gemm_shape, tl(1, 1, l0), mcart(1 + ac), msph(1 + bc))
         ish0 = ish0 + run
      end do
      !$acc exit data delete(tl)
   end subroutine sph_fold_rows_cart_to_sph

end module trc_sph
