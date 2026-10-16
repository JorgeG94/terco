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
module trc_sph
   use pic_types, only: dp
   use trc_api, only: trc_basis_t
   use trc_sph_tables, only: SPH_T
   implicit none
   private

   public :: sph_cart_to_sph, sph_sph_to_cart

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

end module trc_sph
