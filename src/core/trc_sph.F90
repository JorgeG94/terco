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
! orbital-basis pair (mu, nu) alone. Both are `b%nao`/`b%nao_sph` square
! `sph_sph_to_cart`/`sph_cart_to_sph` with one side of the contraction
! dropped, not a new kind of transform.
!
module trc_sph
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_sph_tables, only: SPH_T
   implicit none
   private

   public :: sph_cart_to_sph, sph_sph_to_cart, sph_mo_to_cart, sph_aux3c_to_sph

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
   ! trc_rimp2_driver.F90's `rimp2_mo_to_cart`) and is left alone. `ncol`
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
   ! Tsph(:, :, q) := sum_p SPH_T(p, r, l) Tcart(:, :, p), per auxiliary
   ! SHELL -- the P-index fold a three-centre density-fitting block needs
   ! when the auxiliary basis is spherical (trc_scf_df.F90's `df_jk_setup`,
   ! trc_rimp2_driver.F90's own metric/3c build): the first two indices
   ! (mu, nu), an orbital-basis AO pair, are untouched; only the auxiliary
   ! index crosses the boundary, on the LAST axis instead of a row or
   ! column. `aux` is the block's own basis handle -- normally a
   ! shell-boundary subset (trc_api.F90's `basis_subset`), whose `nao`/
   ! `nao_sph`/`sh_ao`/`sh_ao_sph` already describe that block's own,
   ! LOCAL P numbering -- but the whole auxiliary basis works just as well
   ! for a caller with no blocking at all.
   !
   subroutine sph_aux3c_to_sph(aux, n, tcart, tsph)
      type(trc_basis_t), intent(in) :: aux
      integer, intent(in) :: n
      real(dp), intent(in)  :: tcart(n, n, aux%nao)
      real(dp), intent(out) :: tsph(n, n, aux%nao_sph)
      call aux3c_to_sph_driver(aux%nshell, aux%sh_l, aux%sh_ao, aux%sh_ao_sph, aux%nao, aux%nao_sph, n, &
                               tcart, tsph)
   end subroutine sph_aux3c_to_sph

   subroutine aux3c_to_sph_driver(nshell, sh_l, sh_ao, sh_ao_sph, naux_c, naux_s, n, tcart, tsph)
      integer,  intent(in)  :: nshell, naux_c, naux_s, n
      integer,  intent(in)  :: sh_l(nshell), sh_ao(nshell), sh_ao_sph(nshell)
      real(dp), intent(in)  :: tcart(n, n, naux_c)
      real(dp), intent(out) :: tsph(n, n, naux_s)
      integer :: ij, mu, nu, ish, li, nci, nsi, ai, bi, p, r
      real(dp) :: acc

      ! Flattened (mu, nu) index: every aux SHELL is folded for each (mu,
      ! nu) pair by the inner, sequential loop below -- the bound that
      ! would need a shell-indexed `do concurrent` axis (the one this
      ! module's own header flags) never appears in the parallel header
      ! itself here.
      do concurrent(ij=0:n*n - 1) local(mu, nu, ish, li, nci, nsi, ai, bi, p, r, acc)
         nu = ij/n + 1
         mu = ij - (nu - 1)*n + 1
         do ish = 1, nshell
            li = sh_l(ish)
            nci = (li + 1)*(li + 2)/2
            nsi = 2*li + 1
            ai = sh_ao(ish)
            bi = sh_ao_sph(ish)
            do r = 1, nsi
               acc = 0.0_dp
               do p = 1, nci
                  acc = acc + SPH_T(p, r, li)*tcart(mu, nu, ai + p - 1)
               end do
               tsph(mu, nu, bi + r - 1) = acc
            end do
         end do
      end do
   end subroutine aux3c_to_sph_driver

end module trc_sph
