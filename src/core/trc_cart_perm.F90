!
! AO-index permutation between LEXICOGRAPHIC order (libcint's; what all
! externally authored AO data is in, regardless of TRC_CART_ORDER: libfint
! references, pyscf, SAD JSON guesses) and this build's own Cartesian order.
!
! The two orders give identical physics and differ only in which slot each
! d/f component lands in, so any reader of external AO-indexed data must
! permute it. This module is the one place that permutation is computed,
! from the same generated tables as `trc_cart_comp.inc`.
!
! Users: trc_guess_json.F90 (permutes the SAD guess on the way in) and tests
! (permute the lexicographic oracle, or via scripts/cart_perm.py in Python).
!
! Not a normalisation fix: see trc_ao_norm.F90 for the GAMESS scale.
!
module trc_cart_perm
   use pic_types, only: dp
   implicit none
   private

   public :: trc_cart_order_id, trc_cart_order_name
   public :: cart_perm_l2b, cart_perm_b2l, build_ao_perm

contains

   pure integer function trc_cart_order_id()
      !! 0 (lexicographic) or 1 (gamess), the same encoding TRC_CART_ORDER uses.
#if TRC_CART_ORDER == 1
      trc_cart_order_id = 1
#else
      trc_cart_order_id = 0
#endif
   end function trc_cart_order_id

   pure function trc_cart_order_name() result(nm)
      !! "lexicographic" or "gamess", for anything that wants to print or
      !! record which build produced its output (the C API query function,
      !! a probe file header).
      character(len=:), allocatable :: nm
#if TRC_CART_ORDER == 1
      nm = "gamess"
#else
      nm = "lexicographic"
#endif
   end function trc_cart_order_name

   pure integer function cart_perm_l2b(l, ic) result(r)
      !! The component at LEXICOGRAPHIC local offset `ic` (0-based) of a shell
      !! of angular momentum `l` sits at THIS BUILD's local offset
      !! `cart_perm_l2b(l, ic)` (0-based). Identity at lexicographic order.
      !$acc routine seq
      integer, intent(in) :: l, ic
#include "trc_cart_perm_l2b.inc"
   end function cart_perm_l2b

   pure integer function cart_perm_b2l(l, ic) result(r)
      !! The inverse: the component at THIS BUILD's local offset `ic` (0-based)
      !! sits at LEXICOGRAPHIC local offset `cart_perm_b2l(l, ic)` (0-based).
      !$acc routine seq
      integer, intent(in) :: l, ic
#include "trc_cart_perm_b2l.inc"
   end function cart_perm_b2l

   subroutine build_ao_perm(nshell, sh_l, sh_ao, nao, perm_l2b, perm_b2l)
      !! Full-basis AO permutation, host code (a reader of external data is
      !! never device code): for a basis of `nshell` shells with angular
      !! momenta `sh_l` and 1-based AO offsets `sh_ao` (`trc_basis_t`'s own
      !! fields), build the two `nao`-length permutations.
      !!
      !!   perm_l2b(i): the BUILD-order AO the LEXICOGRAPHIC-order AO `i`
      !!                (1-based) corresponds to.
      !!   perm_b2l(j): the LEXICOGRAPHIC-order AO the BUILD-order AO `j`
      !!                (1-based) corresponds to -- perm_b2l's inverse.
      !!
      !! A reader converting an externally-authored (lexicographic) array
      !! `lex(nao)` into this build's own order does
      !! `build(j) = lex(perm_b2l(j))` for j = 1, nao (a gather); a matrix does
      !! `build(i,j) = lex(perm_b2l(i), perm_b2l(j))`.
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
