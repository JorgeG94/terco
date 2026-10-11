!
! AO-index normalisation at the C API boundary.
!
! terco normalises every Cartesian component of a shell the way libcint
! does: per ANGULAR MOMENTUM, not per component. A d shell's six Cartesian
! functions (XX XY XZ YY YZ ZZ, in whichever order TRC_CART_ORDER picks) all
! carry the SAME normalisation constant, so <XY|XY> is not 1 -- it is 1/3.
!
! GAMESS normalises every component INDIVIDUALLY to unity: <XY|XY> = 1 there,
! always, for every component of every shell. The two conventions describe
! the identical physical orbitals; they differ only in how much of the
! prefactor sits in the basis function versus in the label "this is a d
! function", and the ratio between them is a FIXED, closed-form number per
! component:
!
!     chi_GAMESS = n * chi_terco,   n^2 = (2l-1)!! / ((2a-1)!!(2b-1)!!(2c-1)!!)
!
! for a component x^a y^b z^c of a shell with l = a+b+c, and (-1)!! = 1 by
! convention (so an s or p component, where every exponent pair (2a-1)!! is
! at most (2*1-1)!! = 1, gives a ratio of 1: THIS PART of the scale is 1 for
! s and p, and only d and higher ever need it). scripts/cartorder.py's
! `gamess_norm_sq` is the reference this closed form is checked against
! (test/check_ao_norm.F90, l = 0..4); this module reimplements the same
! formula in Fortran rather than importing it, because nothing on this side
! of the C boundary can call out to Python.
!
! THAT RATIO IS NOT THE WHOLE SCALE FOR s AND p, though -- terco (matching
! libcint) folds a SECOND, unrelated constant into an s or p shell's
! coefficients at build time (`common_fac_sp` below, the real-solid-harmonic
! scale libcint's own convention carries for historical reasons that have
! nothing to do with GAMESS or Cartesian normalisation). It makes terco's
! raw <s|s> come out 1/(4*pi) and <p_x|p_x> 3/(4*pi), not 1, so reaching
! GAMESS's true unit convention needs `gamess_norm` DIVIDED by this second
! factor too. It is exactly 1 for l >= 2, which is why d, f and g never
! showed the problem and s and p quietly did -- see `common_fac_sp`'s own
! comment for how this was found.
!
! ORDER AND NORMALISATION ARE SEPARATE, ON PURPOSE. Which slot a component
! sits in comes from TRC_CART_ORDER (trc_cart_perm.F90); how much that slot
! is scaled by comes from here. A caller in a `TRC_CART_ORDER=gamess` build
! still asks for this module's UNIT convention to get GAMESS's own overlap
! diagonal of 1 -- the two are orthogonal choices, and conflating them is
! exactly the "permutation bug that is really a scale bug" trc_cart_perm.F90
! already warns about for the order half of this same boundary.
!
! THE TRANSFORMATION. Write N = diag(n) over the AOs of one basis (an AO
! carries one n, from its own shell's l and component, regardless of which
! OTHER basis -- orbital or auxiliary -- appears on another index of the
! same array).
!
!   COVARIANT (both matrix indices are AO indices; scale absorbs the same
!   way "this row/column belongs to a d function" does): overlap, kinetic,
!   core Hamiltonian, Fock and its J/K pieces, V_xc, one multipole matrix at
!   a time. M_GAMESS = N M_terco N: `ao_scale_cov`.
!
!   CONTRAVARIANT (a density or an MO-coefficient index runs opposite to an
!   AO's own scale, because a density is built by CONTRACTING two AO
!   indices against an expansion, not by pairing two AO functions the way an
!   integral does): D_GAMESS = N^-1 D_terco N^-1 (`ao_scale_contra`, both
!   indices are AO); C_GAMESS = N^-1 C_terco (`ao_scale_contra_row`, only
!   the AO row is -- the column is an orbital LABEL, not an AO, and is
!   untouched). A density or Hamiltonian GUESS coming IN from a GAMESS-side
!   caller converts the other way, D_terco = N D_GAMESS N: still
!   `ao_scale_cov`, because it is the same diagonal congruence run forward
!   instead of back.
!
!   INVARIANT: energies (total, one-electron, two-electron, XC, RI-MP2),
!   orbital energies/occupations, Tr(DS), convergence flags, iteration
!   counts. These are basis-representation-independent by construction (an
!   eigenvalue of F C = S C eps stays an eigenvalue under ANY congruence
!   that transforms F, S and C consistently) and NOTHING here touches them.
!
! WHO CALLS THIS. `ao_norm_vector` needs one basis's own shell list
! (nshell, sh_l, sh_ao, nao -- trc_basis_t's own fields) and the `kind` a
! context or basis_box carries (trc_capi.F90's TRC_AO_NORM_SHELL/UNIT); it
! is recomputed on every call rather than cached, because it is O(nao) and a
! context's basis can change between calls, and a stale cache is a worse bug
! than a few hundred flops. See trc_capi.F90's own module header for the
! full classification of every AO-indexed C API entry against this scheme.
!
module trc_ao_norm
   use trc_boys, only: dp
   use trc_batch, only: cart_pow
   implicit none
   private

   public :: TRC_AO_NORM_SHELL, TRC_AO_NORM_UNIT
   public :: ao_norm_vector, ao_scale_cov, ao_scale_contra, ao_scale_contra_row
   !> Exported ONLY for test/check_ao_norm.F90's (e): the pure GAMESS ratio,
   !> decoupled from `ao_norm_vector`'s ADDITIONAL common_fac_sp correction
   !> below, so that check can compare against scripts/cartorder.py's
   !> `gamess_norm_sq` exactly, which is this ratio and nothing else.
   public :: gamess_norm

   !> Must match trc_capi.F90's TRC_AO_NORM_SHELL/TRC_AO_NORM_UNIT exactly;
   !> duplicated here (rather than `use`d from trc_capi) because trc_capi
   !> uses THIS module, and a module cannot use its own user the other way.
   integer, parameter :: TRC_AO_NORM_SHELL = 0, TRC_AO_NORM_UNIT = 1

contains

   !> The per-AO scale n(nao) for a basis of `nshell` shells with angular
   !> momenta `sh_l` and 1-based first-AO offsets `sh_ao` (trc_basis_t's own
   !> fields, in THIS BUILD's Cartesian order -- `cart_pow` already reads
   !> TRC_CART_ORDER, so the (a,b,c) exponents handed to `gamess_norm` land
   !> on the right component without this routine knowing which order is in
   !> effect). `kind = TRC_AO_NORM_SHELL` gives back all ones -- terco's own
   !> convention needs no rescaling, and every covariant/contravariant call
   !> below is then a no-op multiply/divide by 1, which is why callers do
   !> not need to special-case SHELL mode themselves.
   pure subroutine ao_norm_vector(nshell, sh_l, sh_ao, nao, kind, nrm)
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
            nrm(base + ic + 1) = gamess_norm(px, py, pz)/common_fac_sp(sh_l(ish))
         end do
      end do
   end subroutine ao_norm_vector

   !> sqrt of GAMESS's per-component normalisation ratio, exactly
   !> scripts/cartorder.py's `gamess_norm` -- see that function's own
   !> docstring for the closed form and why d XY is sqrt3, f XYZ is sqrt15.
   pure real(dp) function gamess_norm(px, py, pz) result(n)
      integer, intent(in) :: px, py, pz
      n = sqrt(real(dfact(px + py + pz), dp)/ &
              real(dfact(px)*dfact(py)*dfact(pz), dp))
   end function gamess_norm

   !> The OTHER baseline shift `ao_norm_vector` has to undo, on top of
   !> `gamess_norm`: terco folds libcint's own real-solid-harmonic scale
   !> into an s or p shell's coefficients at BUILD time (trc_api.F90's
   !> `common_fac_sp`, this is the exact same two constants, duplicated
   !> rather than exported because it is three lines and the alternative is
   !> widening trc_api's public surface for one boundary conversion). It
   !> makes terco's OWN raw <s|s> come out 1/(4*pi) and <p_x|p_x> 3/(4*pi),
   !> not 1 -- discovered the hard way, by this module's own test failing
   !> against an assumed baseline of 1 for s and p. It is 1 for l >= 2
   !> (trc_api.F90's `case default`), which is why d, f and g never needed
   !> this: `gamess_norm` alone was already the right answer for them, and
   !> ONLY s and p were silently wrong.
   pure real(dp) function common_fac_sp(l)
      integer, intent(in) :: l
      select case (l)
      case (0); common_fac_sp = 0.282094791773878143_dp
      case (1); common_fac_sp = 0.488602511902919921_dp
      case default; common_fac_sp = 1.0_dp
      end select
   end function common_fac_sp

   !> (2n-1)!! = 1*3*5*...*(2n-1), with (-1)!! = 1 by convention (an
   !> exponent of 0 contributes no factor, which is what makes an s or p
   !> component's n come out 1). NOT n!! (the double factorial of n
   !> itself) -- (2n-1)!! is the one the Cartesian normalisation ratio
   !> needs, and scripts/cartorder.py's `_dfact` is this same product.
   pure integer function dfact(n) result(r)
      integer, intent(in) :: n
      integer :: k
      r = 1
      do k = 1, n
         r = r*(2*k - 1)
      end do
   end function dfact

   !> M := diag(nrm) M diag(nrm), in place. Both indices of `mat` are AO
   !> indices of the basis `nrm` was built from -- overlap, kinetic, core
   !> Hamiltonian, one spin slice of a Fock/J/K matrix, one multipole
   !> component. Also the direction a GAMESS-side density or Hamiltonian
   !> GUESS converts on the way IN (see the module header): the same
   !> congruence, run forward instead of back.
   pure subroutine ao_scale_cov(nrm, mat)
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

   !> M := diag(1/nrm) M diag(1/nrm), in place. Both indices of `mat` are AO
   !> indices: one spin slice of a density matrix.
   pure subroutine ao_scale_contra(nrm, mat)
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

   !> M(:, j) := M(:, j) / nrm, every column j, in place. The ROW of `mat`
   !> is an AO index; the COLUMN is an orbital label, not an AO, and is left
   !> alone -- one spin slice of an MO-coefficient matrix, (nao, nmo).
   pure subroutine ao_scale_contra_row(nrm, mat)
      real(dp), intent(in) :: nrm(:)
      real(dp), intent(inout) :: mat(:, :)
      integer :: i, j, n, ncol
      n = size(nrm); ncol = size(mat, 2)
      do j = 1, ncol
         do i = 1, n
            mat(i, j) = mat(i, j)/nrm(i)
         end do
      end do
   end subroutine ao_scale_contra_row

end module trc_ao_norm
