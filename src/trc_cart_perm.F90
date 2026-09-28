!
! The AO-index permutation between LEXICOGRAPHIC order (libcint's own, and
! what every source of AO-indexed data from OUTSIDE this build is authored
! in, regardless of TRC_CART_ORDER -- libfint's reference*.dat, pyscf, mqc's
! SAD JSON guesses) and this build's own Cartesian order.
!
! WHY THIS EXISTS
! ----------------
! A gamess build's own physics agrees with a lexicographic one exactly (the
! SiH4 RHF energy comparison in docs/GENERATORS.md-adjacent test notes) --
! the two differ only by which slot each d/f component lands in. So every
! place that reads AO-indexed numbers authored OUTSIDE this build has to
! know that those numbers are lexicographic and land them on the right slot,
! or it is comparing two different bases and calling the difference a bug.
! This module is the ONE place that permutation is computed, from the same
! generated tables `trc_cart_comp.inc` uses, so a reader cannot reimplement
! it slightly wrong.
!
! TWO USES
! --------
!   * A LIBRARY reader of externally-authored data: trc_guess_json.F90's SAD
!     guess is the one in production; permute on the way in, once, and every
!     caller downstream sees an ordinary build-order density.
!   * A TEST comparing against a stored oracle or a live pyscf array: the
!     oracle stays lexicographic (it is not this build's business to know
!     TRC_CART_ORDER when it was written -- it wasn't), and the test permutes
!     either the oracle's index on the way in (Fortran) or its own array on
!     the way in/out (Python, via scripts/cart_perm.py, the same idea against
!     the same cartorder.py tables).
!
! WHAT THIS IS NOT
! ----------------
! Not a normalisation fix. GAMESS's sqrt3/sqrt5 per-component scale is a
! separate, later decision at an API boundary (see cartorder.py's
! `gamess_norm`/`gamess_norm_sq`) and nothing here applies it.
!
module trc_cart_perm
   use trc_boys, only: dp
   implicit none
   private

   public :: trc_cart_order_id, trc_cart_order_name
   public :: cart_perm_l2b, cart_perm_b2l, build_ao_perm

contains

   !> 0 (lexicographic) or 1 (gamess), the same encoding TRC_CART_ORDER uses.
   pure integer function trc_cart_order_id()
#if TRC_CART_ORDER == 1
      trc_cart_order_id = 1
#else
      trc_cart_order_id = 0
#endif
   end function trc_cart_order_id

   !> "lexicographic" or "gamess", for anything that wants to print or
   !> record which build produced its output (the C API query function,
   !> a probe file header).
   pure function trc_cart_order_name() result(nm)
      character(len=:), allocatable :: nm
#if TRC_CART_ORDER == 1
      nm = "gamess"
#else
      nm = "lexicographic"
#endif
   end function trc_cart_order_name

   !> The component at LEXICOGRAPHIC local offset `ic` (0-based) of a shell
   !> of angular momentum `l` sits at THIS BUILD's local offset
   !> `cart_perm_l2b(l, ic)` (0-based). Identity at lexicographic order.
   pure integer function cart_perm_l2b(l, ic) result(r)
      !$acc routine seq
      integer, intent(in) :: l, ic
#include "trc_cart_perm_l2b.inc"
   end function cart_perm_l2b

   !> The inverse: the component at THIS BUILD's local offset `ic` (0-based)
   !> sits at LEXICOGRAPHIC local offset `cart_perm_b2l(l, ic)` (0-based).
   pure integer function cart_perm_b2l(l, ic) result(r)
      !$acc routine seq
      integer, intent(in) :: l, ic
#include "trc_cart_perm_b2l.inc"
   end function cart_perm_b2l

   !> Full-basis AO permutation, host code (a reader of external data is
   !> never device code): for a basis of `nshell` shells with angular
   !> momenta `sh_l` and 1-based AO offsets `sh_ao` (`trc_basis_t`'s own
   !> fields), build the two `nao`-length permutations.
   !>
   !>   perm_l2b(i): the BUILD-order AO the LEXICOGRAPHIC-order AO `i`
   !>                (1-based) corresponds to.
   !>   perm_b2l(j): the LEXICOGRAPHIC-order AO the BUILD-order AO `j`
   !>                (1-based) corresponds to -- perm_b2l's inverse.
   !>
   !> A reader converting an externally-authored (lexicographic) array
   !> `lex(nao)` into this build's own order does
   !> `build(j) = lex(perm_b2l(j))` for j = 1, nao (a gather); a matrix does
   !> `build(i,j) = lex(perm_b2l(i), perm_b2l(j))`.
   subroutine build_ao_perm(nshell, sh_l, sh_ao, nao, perm_l2b, perm_b2l)
      integer, intent(in)  :: nshell, nao
      integer, intent(in)  :: sh_l(nshell), sh_ao(nshell)
      integer, intent(out) :: perm_l2b(nao), perm_b2l(nao)
      integer :: ish, ic, base, ncart

      do ish = 1, nshell
         base = sh_ao(ish) - 1
         ncart = (sh_l(ish) + 1)*(sh_l(ish) + 2)/2
         do ic = 0, ncart - 1
            perm_l2b(base + ic + 1) = base + cart_perm_l2b(sh_l(ish), ic) + 1
            perm_b2l(base + ic + 1) = base + cart_perm_b2l(sh_l(ish), ic) + 1
         end do
      end do
   end subroutine build_ao_perm

end module trc_cart_perm
