!
! AO-index normalisation at the C API boundary.
!
! terco (like libcint) normalises per ANGULAR MOMENTUM: every Cartesian
! component of a d shell shares one constant, so <XY|XY> = 1/3. GAMESS
! normalises each component to unity. For x^a y^b z^c with l = a+b+c:
!
!     chi_GAMESS = n * chi_terco,   n^2 = (2l-1)!! / ((2a-1)!!(2b-1)!!(2c-1)!!)
!
! with (-1)!! = 1, so n = 1 for s and p. `gamess_norm` alone is the whole
! scale for every l: the `common_fac_sp` folded into s/p coefficients by
! `basis_build` is cancelled by the callers' own normalisation, leaving raw
! <s|s> = <p|p> = 1 (checked against pyscf/libcint). An earlier version also
! divided by `common_fac_sp`, cancelling a matching double fold in
! trc_basis_json; both were fixed together. Reference: scripts/cartorder.py
! `gamess_norm_sq`, checked in test/check_ao_norm.F90.
!
! ORDER AND NORMALISATION ARE SEPARATE: slot order comes from TRC_CART_ORDER
! (trc_cart_perm.F90), scaling from here. With N = diag(n) over one basis:
!   COVARIANT (both indices AO: S, T, H, F, J/K, V_xc, multipoles):
!     M_GAMESS = N M N (`ao_scale_cov`); also used for GAMESS guesses coming in.
!   CONTRAVARIANT (D, C): D_GAMESS = N^-1 D N^-1 (`ao_scale_contra`),
!     C_GAMESS = N^-1 C (`ao_scale_contra_row`; the column is an orbital label).
!   INVARIANT: energies, orbital energies, occupations, Tr(DS), flags.
!
! `ao_norm_vector` takes one basis's shell list and the `kind` from trc_capi.F90
! (TRC_AO_NORM_SHELL/UNIT) and is recomputed per call (O(nao); no stale cache).
! See trc_capi.F90's header for how every AO-indexed C API entry is classified.
!
module trc_ao_norm
   use pic_types, only: dp
   use trc_batch, only: cart_pow
   implicit none
   private

   public :: TRC_AO_NORM_SHELL, TRC_AO_NORM_UNIT
   public :: ao_norm_vector, ao_scale_cov, ao_scale_contra, ao_scale_contra_row
   ! Exported for test/check_ao_norm.F90's (e), so that check can compare
   ! against scripts/cartorder.py's `gamess_norm_sq` directly -- this IS
   ! the whole of `ao_norm_vector`'s UNIT scale now, not a part of it.
   public :: gamess_norm

   integer, parameter :: TRC_AO_NORM_SHELL = 0, TRC_AO_NORM_UNIT = 1
   !! Must match trc_capi.F90's TRC_AO_NORM_SHELL/TRC_AO_NORM_UNIT exactly;
   !! duplicated here (rather than `use`d from trc_capi) because trc_capi
   !! uses THIS module, and a module cannot use its own user the other way.

contains

   pure subroutine ao_norm_vector(nshell, sh_l, sh_ao, nao, kind, nrm)
      !! The per-AO scale n(nao) for a basis of `nshell` shells with angular
      !! momenta `sh_l` and 1-based first-AO offsets `sh_ao` (trc_basis_t's own
      !! fields, in THIS BUILD's Cartesian order -- `cart_pow` already reads
      !! TRC_CART_ORDER, so the (a,b,c) exponents handed to `gamess_norm` land
      !! on the right component without this routine knowing which order is in
      !! effect). `kind = TRC_AO_NORM_SHELL` gives back all ones -- terco's own
      !! convention needs no rescaling, and every covariant/contravariant call
      !! below is then a no-op multiply/divide by 1, which is why callers do
      !! not need to special-case SHELL mode themselves.
      integer, intent(in)  :: nshell, sh_l(nshell), sh_ao(nshell), nao, kind
      real(dp), intent(out) :: nrm(nao)
      integer :: ish, ic, ncart, base, px, py, pz

      if (kind /= TRC_AO_NORM_UNIT) then
         nrm = 1.0_dp
         return
      end if
      do ish = 1, nshell
         base = sh_ao(ish) - 1
         ncart = (sh_l(ish) + 1)*(sh_l(ish) + 2)/2
         do ic = 0, ncart - 1
            call cart_pow(sh_l(ish), ic, px, py, pz)
            ! `gamess_norm` alone, no extra `common_fac_sp` divide -- see the
            ! module header for why an earlier version of this line carried
            ! one, and why it was wrong.
            nrm(base + ic + 1) = gamess_norm(px, py, pz)
         end do
      end do
   end subroutine ao_norm_vector

   pure real(dp) function gamess_norm(px, py, pz) result(n)
      !! sqrt of GAMESS's per-component normalisation ratio, exactly
      !! scripts/cartorder.py's `gamess_norm` -- see that function's own
      !! docstring for the closed form and why d XY is sqrt3, f XYZ is sqrt15.
      integer, intent(in) :: px, py, pz
      n = sqrt(real(dfact(px + py + pz), dp)/ &
              real(dfact(px)*dfact(py)*dfact(pz), dp))
   end function gamess_norm

   pure integer function dfact(n) result(r)
      !! (2n-1)!! = 1*3*5*...*(2n-1), with (-1)!! = 1 by convention (an
      !! exponent of 0 contributes no factor, which is what makes an s or p
      !! component's n come out 1). NOT n!! (the double factorial of n
      !! itself) -- (2n-1)!! is the one the Cartesian normalisation ratio
      !! needs, and scripts/cartorder.py's `_dfact` is this same product.
      integer, intent(in) :: n
      integer :: k
      r = 1
      do k = 1, n
         r = r*(2*k - 1)
      end do
   end function dfact

   pure subroutine ao_scale_cov(nrm, mat)
      !! M := diag(nrm) M diag(nrm), in place. Both indices of `mat` are AO
      !! indices of the basis `nrm` was built from -- overlap, kinetic, core
      !! Hamiltonian, one spin slice of a Fock/J/K matrix, one multipole
      !! component. Also the direction a GAMESS-side density or Hamiltonian
      !! GUESS converts on the way IN (see the module header): the same
      !! congruence, run forward instead of back.
      real(dp), intent(in) :: nrm(:)
      real(dp), intent(inout) :: mat(:, :)
      integer :: i, j, n
      n = size(nrm)
      do j = 1, n
         do i = 1, n
            mat(i, j) = mat(i, j)*nrm(i)*nrm(j)
         end do
      end do
   end subroutine ao_scale_cov

   pure subroutine ao_scale_contra(nrm, mat)
      !! M := diag(1/nrm) M diag(1/nrm), in place. Both indices of `mat` are AO
      !! indices: one spin slice of a density matrix.
      real(dp), intent(in) :: nrm(:)
      real(dp), intent(inout) :: mat(:, :)
      integer :: i, j, n
      n = size(nrm)
      do j = 1, n
         do i = 1, n
            mat(i, j) = mat(i, j)/(nrm(i)*nrm(j))
         end do
      end do
   end subroutine ao_scale_contra

   pure subroutine ao_scale_contra_row(nrm, mat)
      !! M(:, j) := M(:, j) / nrm, every column j, in place. The ROW of `mat`
      !! is an AO index; the COLUMN is an orbital label, not an AO, and is left
      !! alone -- one spin slice of an MO-coefficient matrix, (nao, nmo).
      real(dp), intent(in) :: nrm(:)
      real(dp), intent(inout) :: mat(:, :)
      integer :: i, j, n, ncol
      n = size(nrm)
      ncol = size(mat, 2)
      do j = 1, ncol
         do i = 1, n
            mat(i, j) = mat(i, j)/nrm(i)
         end do
      end do
   end subroutine ao_scale_contra_row

end module trc_ao_norm
