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
! `NPMAX` below), never to a runtime shape (gpu-pitfalls.md's "no automatic
! arrays in routine seq"). `trc_basis_t` never holds a general contraction
! (trc_ecp_radial.F90's own header), so every "nc" PySCF's C and libfint's
! port carry is 1 here and the corresponding loop is simply absent.
!
! NORMALISATION: common_fac_sp is folded into `trc_basis_t%sh_c` exactly
! once, in basis construction (AGENTS.md). libfint's own `scale_coeff`
! would fold it a SECOND time (it reads libcint's `env`, which does not
! carry it), so this port drops that factor here, keeping only the
! 4*pi / 16*pi^2 PySCF's own formula also carries -- ECP_DESIGN.md's note.
!
! NO 3-CENTRE OVERLAP SCREEN (yet). libfint's `ecp_check_3c_overlap` skips
! a (shell pair, ECP atom) that cannot contribute; this port always
! evaluates every one. That costs time on a molecule with several ECP
! centres far from most shell pairs and changes no number -- the actual
! quadrature, run honestly on a far-separated pair, already comes out
! negligibly small (the `ecp_ints_ih_far` case: I-H at 25 bohr) -- so it is
! deferred rather than risking a screening bound built on an unverified
! "smallest exponent last" ordering of a BSE-JSON-read ECP channel or
! orbital shell (see git history/PR description for the measurement that
! would justify adding it back).
#include "inc/trc_dc_index.inc"

module trc_ecp
   use pic_types, only: dp
   use trc_basis_types, only: trc_basis_t
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_quadrature, only: ECP_LMAX, ECP_LEVEL0, ECP_LEVEL_MAX, ECP_NRS, &
                                 ecp_quadrature_gauss_chebyshev
   use trc_ecp_angular, only: ECP_ORB_LMAX, CART_POW_Y, CART_POW_Z, &
                              ecp_type1_rad_ang, ecp_type2_facs_ang
   use trc_ecp_radial, only: ecp_radial_rad_part, ecp_radial_type1_static_facs, &
                             ecp_radial_type1_rad_part_levels, ecp_radial_type2_facs_rad, &
                             ecp_radial_type2_lab_levels, ECP_NLEVEL
   use trc_sph, only: sph_cart_to_sph
   use trc_cart_perm, only: cart_perm_l2b
   use trc_error, only: error_t, ERROR_VALIDATION
   use pic_logger, only: global_logger
   use pic_strings, only: to_string
   implicit none
   private
   public :: trc_ecp_matrix, trc_ecp_matrix_cart

   integer, parameter :: NFMAX = (ECP_ORB_LMAX + 1)*(ECP_ORB_LMAX + 2)/2     !! 10
   integer, parameter :: NPMAX = 16       !! primitives per orbital shell or ECP row
   integer, parameter :: D1MAX = 2*ECP_ORB_LMAX + 1                          !! 7
   integer, parameter :: D2MAX = D1MAX*D1MAX                                 !! 49
   integer, parameter :: D3MAX = D2MAX*D1MAX                                 !! 343
   integer, parameter :: LC1MAX = ECP_ORB_LMAX + ECP_LMAX + 1                 !! 9
   integer, parameter :: D2MAX2 = LC1MAX*LC1MAX                              !! 81
   integer, parameter :: D3MAX2 = D1MAX*D2MAX2                               !! 567
   integer, parameter :: SLOT_SIZE = NFMAX*NFMAX                                  !! 100: one (shell pair, atom) slot
   !! K1b's per-(shell, channel) global caches, sized to this module's own
   !! ceilings and read at the caller's actual, smaller extent -- exactly
   !! `ecp_radial_type1_static_facs`'s and `ecp_type2_facs_ang`'s own
   !! output bounds (not the oversized ifac/jfac/angi/angj this replaces:
   !! those were declared NFMAX*D1MAX**3 and (ECP_ORB_LMAX+1)*NFMAX*
   !! (2*ECP_LMAX+1)*LC1MAX respectively in an earlier version of this
   !! file, but only the first FAC1SIZE/ANG2SIZE entries were ever read).
   integer, parameter :: FAC1SIZE = NFMAX*(ECP_ORB_LMAX + 1)**3                  !! 640
   integer, parameter :: ANG2SIZE = (ECP_ORB_LMAX + 1)*NFMAX*(2*ECP_LMAX + 1)*LC1MAX  !! 3960

   real(dp), parameter :: ECP_NEGLIGIBLE = 1.0e-20_dp
   !! A contribution below this, to a matrix whose real entries run from
   !! ~1e-3 to a few hundred (measured: Xe/def2-SVP's U diagonal reaches
   !! 206), can never matter past the 1e-9 tolerance anything here is
   !! validated to. Used to stop the adaptive radial loop from chasing a
   !! RELATIVE 1e-12 match on a term that is already ABSOLUTELY negligible
   !! -- which a large separation (ecp_ints_ih_far's 25 bohr) or a tight
   !! core primitive pair reaches often, and which the relative test alone
   !! can fail to confirm within this port's ECP_LEVEL_MAX even though the
   !! term's actual contribution to U is already settled at zero.

contains

   pure integer function ncart_l(l) result(n)
      !$acc routine seq
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)/2
   end function ncart_l

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

   subroutine trc_ecp_matrix_cart(b, ecp, mcart, error)
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
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), allocatable, intent(out) :: mcart(:, :)
      type(error_t), intent(inout) :: error

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

      allocate (mcart(b%nao, b%nao))
      mcart = 0.0_dp

      if (ecp%has_ecp) then
         ! A malformed or foreign ECP file must fail here, not read past the
         ! caches below: they are sized to ECP_LMAX, and the rows index
         ! b%at_r by atom.
         if (ecp%natm /= b%natm) then
            call error%set(ERROR_VALIDATION, "trc_ecp_matrix: the ECP describes "// &
                           to_string(ecp%natm)//" atoms, the basis "//to_string(b%natm))
            return
         end if
         if (ecp%nrow > 0) then
            if (maxval(ecp%ecp_np) > NPMAX) then
               call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row has more "// &
                              "primitives than this port's fixed NPMAX ceiling")
               return
            end if
            if (minval(ecp%ecp_np) < 1) then
               call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row has no primitives")
               return
            end if
            if (minval(ecp%ecp_at) < 1 .or. maxval(ecp%ecp_at) > b%natm) then
               call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row names an atom "// &
                              "outside 1.."//to_string(b%natm))
               return
            end if
            if (minval(ecp%ecp_l) < -1 .or. maxval(ecp%ecp_l) > ECP_LMAX) then
               call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP channel has l outside "// &
                              "-1 (local) .. "//to_string(ECP_LMAX)//", this port's ECP_LMAX")
               return
            end if
         end if
         call ecp_matrix_cart(b, ecp, mcart, error)
         if (error%has_error()) return
      end if
   end subroutine trc_ecp_matrix_cart

   subroutine ecp_matrix_cart(b, ecp, mcart, error)
      !!
      !! The Cartesian build: groups the ECP rows into (atom, l) channels, then
      !! K1/K2/K3 as this module's own header describes.
      !!
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), intent(inout) :: mcart(b%nao, b%nao)
      type(error_t), intent(inout) :: error

      real(dp) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), allocatable :: ur_cache(:, :)
      real(dp), allocatable :: fac1_cache(:, :, :), fac2_cache(:, :, :), ang2_cache(:, :, :)
      integer, allocatable :: grp_atom(:), grp_l(:), grp_row0(:), grp_row1(:)
      integer, allocatable :: atm_ecp(:), atm_grp0(:), atm_grp1(:)
      integer, allocatable :: pr_i(:), pr_j(:)
      real(dp), allocatable :: block_buffer(:, :, :)
      integer, allocatable :: conv_fail(:, :)
      integer :: ngroup, natm_ecp, npair, nshell
      integer :: i, j, p, g, a, nfail, t, ntask, ntask1b

      call build_groups(ecp, grp_atom, grp_l, grp_row0, grp_row1, ngroup)
      if (ngroup == 0) return
      call build_atom_list(grp_atom, ngroup, atm_ecp, atm_grp0, atm_grp1, natm_ecp)
      nshell = b%nshell

      call ecp_quadrature_gauss_chebyshev(rs, ws, ECP_NRS)
      allocate (ur_cache(0:ECP_NRS - 1, ngroup))

      !$acc enter data copyin(rs, ecp, ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c) &
      !$acc            copyin(grp_row0, grp_row1) create(ur_cache)
      do concurrent(g=1:ngroup)
         call ecp_k1_item(ur_cache(:, g), rs, grp_row0(g), grp_row1(g), &
                          ecp%nrow, ecp%maxnp, ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c)
      end do
      !$acc exit data delete(ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c, grp_row0, grp_row1, ecp)

      npair = nshell*(nshell + 1)/2
      allocate (pr_i(npair), pr_j(npair))
      p = 0
      do i = 1, nshell
         do j = 1, i
            p = p + 1
            pr_i(p) = i
            pr_j(p) = j
         end do
      end do

      allocate (block_buffer(SLOT_SIZE, npair, max(natm_ecp, 1)), conv_fail(npair, max(natm_ecp, 1)))
      block_buffer = 0.0_dp
      conv_fail = 0

      ! K1b: once per (orbital shell, ECP channel group) -- not once per
      ! shell pair -- the semi-local radial/angular factors (fac2_cache,
      ! ang2_cache) and the local term's static centre-shift factors
      ! (fac1_cache). Flattened for the same reason K2 below is: `ngroup`
      ! is as small a run-time extent as `natm_ecp`.
      allocate (fac1_cache(0:FAC1SIZE - 1, nshell, ngroup))
      allocate (fac2_cache(0:ECP_NRS*LC1MAX - 1, nshell, ngroup))
      allocate (ang2_cache(0:ANG2SIZE - 1, nshell, ngroup))

      !$acc enter data copyin(ws, b, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%at_r, grp_atom, grp_l) &
      !$acc            create(fac1_cache, fac2_cache, ang2_cache)
      ntask1b = nshell*ngroup
      do concurrent(t=1:ntask1b) local(i, g)
         g = (t - 1)/nshell + 1
         i = t - (g - 1)*nshell
         call ecp_k1b_item(TRC_DC_INDEX(i), TRC_DC_INDEX(g), nshell, b%sh_l, b%sh_np, b%maxnp, b%sh_e, b%sh_c, b%sh_r, &
                           b%natm, b%at_r, ngroup, grp_atom, grp_l, rs, &
                           fac1_cache(:, i, g), fac2_cache(:, i, g), ang2_cache(:, i, g))
      end do

      !$acc enter data copyin(pr_i, pr_j, b%sh_ao) copyin(atm_ecp, atm_grp0, atm_grp1) &
      !$acc            create(block_buffer, conv_fail)
      if (natm_ecp > 0) then
         ! Flattened, not do concurrent(p=1:npair, a=1:natm_ecp): a
         ! two-index do concurrent whose outer extent is small and set at
         ! run time has, in testing here, failed the kernel launch outright
         ! on the GPU build (nvfortran, cc70) rather than merely dropping a
         ! slice -- gpu-pitfalls.md's own warning about this shape, and
         ! natm_ecp (the number of distinct ECP atoms) is exactly the kind
         ! of small run-time extent it means: 1 or 2 for every case
         ! validated so far.
         ntask = npair*natm_ecp
         do concurrent(t=1:ntask) local(p, a)
            a = (t - 1)/npair + 1
            p = t - (a - 1)*npair
            call ecp_pair_atom_item(pr_i(p), pr_j(p), atm_ecp(a), atm_grp0(a), atm_grp1(a), &
                                    ngroup, grp_l, &
                                    nshell, b%sh_l, b%sh_np, b%maxnp, b%sh_e, b%sh_c, b%sh_r, &
                                    b%natm, b%at_r, &
                                    rs, ws, ur_cache, fac1_cache, fac2_cache, ang2_cache, &
                                    block_buffer(:, p, a), conv_fail(p, a))
         end do
      end if
      !$acc update self(block_buffer, conv_fail)
      !$acc exit data delete(ws, pr_i, pr_j, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%sh_ao, b%at_r, b, &
      !$acc                  grp_l, grp_atom, atm_ecp, atm_grp0, atm_grp1, rs, ur_cache, block_buffer, conv_fail) &
      !$acc            delete(fac1_cache, fac2_cache, ang2_cache)

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

      ! K3: reduce the per-atom slots and write the Cartesian block (and
      ! its symmetric mirror -- U is symmetric, scalar ECPs are real, and
      ! only i>=j pairs were built).
      do p = 1, npair
         call deposit_block(b, mcart, pr_i(p), pr_j(p), block_buffer(:, p, :), natm_ecp)
      end do

      deallocate (ur_cache, fac1_cache, fac2_cache, ang2_cache, grp_atom, grp_l, grp_row0, grp_row1)
      deallocate (atm_ecp, atm_grp0, atm_grp1, pr_i, pr_j, block_buffer, conv_fail)
   end subroutine ecp_matrix_cart

   subroutine deposit_block(b, mcart, i, j, slots, natm_ecp)
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(inout) :: mcart(b%nao, b%nao)
      integer, intent(in) :: i, j, natm_ecp
      real(dp), intent(in) :: slots(SLOT_SIZE, max(natm_ecp, 1))
      integer :: li, lj, nfi, nfj, ai, aj, mi, mj, a, pi, pj
      real(dp) :: acc

      li = b%sh_l(i)
      lj = b%sh_l(j)
      nfi = ncart_l(li)
      nfj = ncart_l(lj)
      ai = b%sh_ao(i)
      aj = b%sh_ao(j)
      do mj = 1, nfj
         do mi = 1, nfi
            acc = 0.0_dp
            do a = 1, natm_ecp
               acc = acc + slots((mj - 1)*NFMAX + mi, a)
            end do
            ! The slots are lexicographic (this module's angular tables);
            ! the matrix is in this build's order, like S, T and V.
            pi = ai + cart_perm_l2b(li, mi - 1)
            pj = aj + cart_perm_l2b(lj, mj - 1)
            mcart(pi, pj) = mcart(pi, pj) + acc
            if (i /= j) mcart(pj, pi) = mcart(pj, pi) + acc
         end do
      end do
   end subroutine deposit_block

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

   pure subroutine ecp_k1b_item(ish, g, nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                natm, at_r, ngroup, grp_atom, grp_l, rs, &
                                fac1_slot, fac2_slot, ang2_slot)
      !!
      !! K1b's item: one (orbital shell, ECP channel group). Writes ONLY the
      !! cache this group's own type (local or semi-local) needs -- the other
      !! two output arrays are left untouched, which is fine: K2 never reads a
      !! local group's fac2/ang2 slot or a semi-local group's fac1 slot.
      !!
      !$acc routine seq
      integer, intent(in) :: ish, g, nshell, maxnp, natm, ngroup
      integer, intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell), sh_r(3, nshell)
      real(dp), intent(in) :: at_r(3, natm)
      integer, intent(in) :: grp_atom(ngroup), grp_l(ngroup)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1)
      real(dp), intent(out) :: fac1_slot(0:FAC1SIZE - 1)
      real(dp), intent(out) :: fac2_slot(0:ECP_NRS*LC1MAX - 1)
      real(dp), intent(out) :: ang2_slot(0:ANG2SIZE - 1)

      integer :: li, np, lc
      real(dp) :: rca(3)
      real(dp) :: dca

      li = sh_l(ish)
      np = sh_np(ish)
      lc = grp_l(g)
      rca = at_r(:, grp_atom(g)) - sh_r(:, ish)

      if (lc == -1) then
         call ecp_radial_type1_static_facs(fac1_slot, li, rca)
      else
         dca = sqrt(sum(rca**2))
         call ecp_radial_type2_facs_rad(fac2_slot, li, np, sh_e(1:np, ish), sh_c(1:np, ish), &
                                        lc, dca, rs, 0, 1, ECP_NRS)
         call ecp_type2_facs_ang(ang2_slot, li, lc, rca)
      end if
   end subroutine ecp_k1b_item

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
      !! K1b's own global caches -- see its own header. Read here at
      !! (:,ish,g) and (:,jsh,g): the pair's two shells against the SAME
      !! channel group, never recomputed.
      real(dp), intent(in) :: fac1_cache(0:FAC1SIZE - 1, nshell, ngroup)
      real(dp), intent(in) :: fac2_cache(0:ECP_NRS*LC1MAX - 1, nshell, ngroup)
      real(dp), intent(in) :: ang2_cache(0:ANG2SIZE - 1, nshell, ngroup)
      real(dp), intent(out) :: slot(SLOT_SIZE)
      integer, intent(out) :: conv_fail

      real(dp) :: gctr(NFMAX, NFMAX)
      real(dp) :: rc(3), ri(3), rj(3)
      integer :: li, lj, npi, npj, mi, mj, g, local_fail

      li = sh_l(ish)
      lj = sh_l(jsh)
      npi = sh_np(ish)
      npj = sh_np(jsh)
      ri = sh_r(:, ish)
      rj = sh_r(:, jsh)
      rc = at_r(:, atm_id)
      gctr = 0.0_dp
      conv_fail = 0

      do g = g0, g1
         if (grp_l(g) == -1) then
            call ecp_type1_item(gctr, li, lj, npi, npj, &
                                sh_e(1:npi, ish), sh_c(1:npi, ish), sh_e(1:npj, jsh), sh_c(1:npj, jsh), &
                                ri, rj, rc, rs, ws, ur_cache(:, g), &
                                fac1_cache(:, ish, g), fac1_cache(:, jsh, g), local_fail)
         else
            call ecp_type2_item(gctr, li, lj, grp_l(g), &
                                rs, ws, ur_cache(:, g), &
                                fac2_cache(:, ish, g), fac2_cache(:, jsh, g), &
                                ang2_cache(:, ish, g), ang2_cache(:, jsh, g), local_fail)
         end if
         conv_fail = max(conv_fail, local_fail)
      end do

      do mj = 1, NFMAX
         do mi = 1, NFMAX
            slot((mj - 1)*NFMAX + mi) = gctr(mi, mj)
         end do
      end do
   end subroutine ecp_pair_atom_item

   pure logical function close_enough(x, y) result(c)
      !$acc routine seq
      real(dp), intent(in) :: x, y
      c = abs(x - y) <= 1.0e-12_dp*max(abs(x), abs(y))
   end function close_enough

   pure subroutine pick_converged_level(ld, d2, rad_levels, picked, failed)
      !! Picks the level a (primitive pair, type 1) or (lambda, type 2)
      !! block converged at, from its ECP_NLEVEL-column snapshot: the same
      !! "two consecutive close levels" test the original per-level loop
      !! ran while the block's own state persisted across levels, now run
      !! once over an already-materialised snapshot. `ld`, the
      !! dummy-bound explicit shape (like `sh_e(maxnp, nshell)` elsewhere
      !! in this file), is the CALLER's own fixed ceiling (D2MAX or
      !! D2MAX2), not the actual `d2`: passing the caller's array
      !! un-sliced, with its own real leading dimension, is what keeps
      !! this a reference rather than a copy-in/copy-out temporary sized
      !! by the runtime `d2` -- see `ecp_radial_type2_lab_levels`'s own
      !! header for the same trap this avoids. `picked` defaults to the
      !! last level and `failed` to true so a caller that never calls this
      !! for an exhausted block still gets the same "use the finest
      !! level, flag non-convergence" fallback the original fell through
      !! to.
      !$acc routine seq
      integer, intent(in) :: ld, d2
      real(dp), intent(in) :: rad_levels(0:ld - 1, 0:ECP_NLEVEL - 1)
      integer, intent(out) :: picked
      logical, intent(out) :: failed

      integer :: lvl, n, converged
      logical :: is_close

      converged = 0
      picked = ECP_NLEVEL - 1
      failed = .true.
      do lvl = 0, ECP_NLEVEL - 1
         is_close = .true.
         do n = 0, d2 - 1
            if (lvl == 0) then
               if (.not. close_enough(0.0_dp, rad_levels(n, lvl))) then
                  is_close = .false.
                  exit
               end if
            else
               if (.not. close_enough(rad_levels(n, lvl - 1), rad_levels(n, lvl))) then
                  is_close = .false.
                  exit
               end if
            end if
         end do
         if (is_close) then
            converged = converged + 1
            if (converged >= 2) then
               picked = lvl
               failed = .false.
               exit
            end if
         else
            converged = 0
         end if
      end do
   end subroutine pick_converged_level

   pure subroutine pick_converged_level_t2(ld, d2, rad_levels, picked, failed)
      !! Same as `pick_converged_level`, for a type-2 lambda block: a
      !! block whose every entry is already below ECP_NEGLIGIBLE is
      !! accepted at whatever level it first reaches that floor, without
      !! the relative test -- see ecp_type2_item's own header for why (a
      !! term already this small can fail a RELATIVE test indefinitely
      !! without ever being wrong to within this port's tolerance,
      !! `ecp_ints_ih_far`'s 25 bohr I-H).
      !$acc routine seq
      integer, intent(in) :: ld, d2
      real(dp), intent(in) :: rad_levels(0:ld - 1, 0:ECP_NLEVEL - 1)
      integer, intent(out) :: picked
      logical, intent(out) :: failed

      integer :: lvl, n, converged
      logical :: is_close

      converged = 0
      picked = ECP_NLEVEL - 1
      failed = .true.
      do lvl = 0, ECP_NLEVEL - 1
         if (maxval(abs(rad_levels(0:d2 - 1, lvl))) < ECP_NEGLIGIBLE) then
            is_close = .true.
         else
            is_close = .true.
            do n = 0, d2 - 1
               if (lvl == 0) then
                  if (.not. close_enough(0.0_dp, rad_levels(n, lvl))) then
                     is_close = .false.
                     exit
                  end if
               else
                  if (.not. close_enough(rad_levels(n, lvl - 1), rad_levels(n, lvl))) then
                     is_close = .false.
                     exit
                  end if
               end if
            end do
         end if
         if (is_close) then
            converged = converged + 1
            if (converged >= 2) then
               picked = lvl
               failed = .false.
               exit
            end if
         else
            converged = 0
         end if
      end do
   end subroutine pick_converged_level_t2

   pure subroutine ecp_type1_item(gctr, li, lj, npi, npj, ei, ci, ej, cj, ri, rj, rc, &
                                  rs, ws, ur, ifac, jfac, conv_fail)
      !!
      !! The type-1 (local) term for one (shell pair, ECP channel). Per-primitive-
      !! pair convergence, exactly as libfint's `ecp_type1_cart`: each (ip,jp)'s
      !! own radial integral is refined to its OWN converged level (picked by
      !! `pick_converged_level` from `ecp_radial_type1_rad_part_levels`'s whole
      !! snapshot, computed and consumed one pair at a time) and the angular
      !! contraction runs only once, after every pair has been picked (or the
      !! level ceiling is hit, which sets `conv_fail`). `ifac`/`jfac` are
      !! K1b's global cache, not rebuilt here.
      !!
      !$acc routine seq
      integer, intent(in) :: li, lj, npi, npj
      real(dp), intent(in) :: ei(npi), ci(npi), ej(npj), cj(npj)
      real(dp), intent(in) :: ri(3), rj(3), rc(3)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(in) :: ifac(0:FAC1SIZE - 1), jfac(0:FAC1SIZE - 1)
      real(dp), intent(inout) :: gctr(NFMAX, NFMAX)
      integer, intent(out) :: conv_fail

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: rca(3), rcb(3), rij(3)
      real(dp) :: cei(NPMAX), cej(NPMAX)
      real(dp) :: rad_levels(0:D2MAX - 1, 0:ECP_NLEVEL - 1)
      real(dp) :: rad_ang(0:D3MAX - 1), rad_ang_all(0:D3MAX - 1)
      integer :: lilj1, d1, d2, d3, di1, dj1
      integer :: ip, jp, mi, mj, n, i1, i2, i3, j1, j2, j3, ix, iy, iz, jx, jy, jz
      integer :: picked
      real(dp) :: r2ca, r2cb, fac
      logical :: failed

      lilj1 = li + lj + 1
      d1 = lilj1
      d2 = d1*d1
      d3 = d2*d1
      di1 = li + 1
      dj1 = lj + 1

      rca = rc - ri
      rcb = rc - rj
      r2ca = sum(rca**2)
      r2cb = sum(rcb**2)
      do ip = 1, npi
         cei(ip) = ci(ip)*exp(-ei(ip)*r2ca)*4.0_dp*PI
      end do
      do jp = 1, npj
         cej(jp) = cj(jp)*exp(-ej(jp)*r2cb)*4.0_dp*PI
      end do

      conv_fail = 0
      rad_ang_all = 0.0_dp
      do ip = 1, npi
         do jp = 1, npj
            ! A primitive pair whose centre-shift prefactor has already
            ! decayed past this floor contributes negligibly to gctr
            ! (fac*rad_ang below) however rad_levels' own magnitude looks:
            ! `ecp_radial_type1_rad_part_levels`'s exponent-shift (complete
            ! the square around kaij, rather than evaluate exp(-aij*rn^2)
            ! directly) trades a huge-but-finite intermediate for one that
            ! would otherwise underflow outright, and that intermediate's
            ! own growth tracks the SAME exponents this prefactor decays
            ! by, so it never overtakes a floor this far below any value
            ! that could matter (measured: a prefactor of 8.4e-251 next to
            ! a "val" of 8.6e+187, net contribution ~1e-63) -- so this pair
            ! never needs the quadrature run on it at all. A tight core
            ! primitive far from the ECP centre (ecp_ints_ih_far's 25 bohr
            ! I-H, ecp_ints_csh's Cs 1s next to the same-atom local channel)
            ! reaches this routinely.
            if (abs(cei(ip)*cej(jp)) < 1.0e-150_dp) cycle
            rij = ei(ip)*rca + ej(jp)*rcb
            call ecp_radial_type1_rad_part_levels(rad_levels, li + lj, &
                                      sqrt(sum(rij**2))*2.0_dp, ei(ip) + ej(jp), ur, rs, ws)
            call pick_converged_level(D2MAX, d2, rad_levels, picked, failed)
            if (failed) conv_fail = 1
            call ecp_type1_rad_ang(rad_ang, li + lj, rij, rad_levels(0:d2 - 1, picked))
            fac = cei(ip)*cej(jp)
            do n = 0, d3 - 1
               rad_ang_all(n) = rad_ang_all(n) + fac*rad_ang(n)
            end do
         end do
      end do

      do mi = 0, ncart_l(li) - 1
         iy = CART_POW_Y(mi)
         iz = CART_POW_Z(mi)
         ix = li - iy - iz
         do mj = 0, ncart_l(lj) - 1
            jy = CART_POW_Y(mj)
            jz = CART_POW_Z(mj)
            jx = lj - jy - jz
            do i1 = 0, ix
               do i2 = 0, iy
                  do i3 = 0, iz
                     do j1 = 0, jx
                        do j2 = 0, jy
                           do j3 = 0, jz
                              gctr(mi + 1, mj + 1) = gctr(mi + 1, mj + 1) &
                                 + ifac(mi*di1*di1*di1 + i1*di1*di1 + i2*di1 + i3) &
                                 *jfac(mj*dj1*dj1*dj1 + j1*dj1*dj1 + j2*dj1 + j3) &
                                 *rad_ang_all((i1 + j1)*d2 + (i2 + j2)*d1 + i3 + j3)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine ecp_type1_item

   pure subroutine ecp_type2_item(gctr, li, lj, lc, &
                                  rs, ws, ur, radi, radj, angi, angj, conv_fail)
      !!
      !! The type-2 (semi-local) term for one (shell pair, ECP channel of
      !! angular momentum lc). No nc/ic loop (trc_basis_t never holds a general
      !! contraction), so the convergence state is indexed by `lab` alone.
      !! `radi`/`radj`/`angi`/`angj` are K1b's global cache -- nothing here
      !! rebuilds a radial or angular factor, or holds one at ECP_NRS size:
      !! each `lab` (0..li+lj) is picked to its OWN converged level one at a
      !! time (`ecp_radial_type2_lab_levels` + `pick_converged_level_t2`), and
      !! only the picked, SLOT_SIZE-scale `rad_all` persists across lab values
      !! -- it has to, since the final contraction below indexes it by i+j,
      !! which repeats across different (i,j) pairs.
      !!
      !$acc routine seq
      integer, intent(in) :: li, lj, lc
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(in) :: radi(0:ECP_NRS*LC1MAX - 1), radj(0:ECP_NRS*LC1MAX - 1)
      real(dp), intent(in) :: angi(0:ANG2SIZE - 1), angj(0:ANG2SIZE - 1)
      real(dp), intent(inout) :: gctr(NFMAX, NFMAX)
      integer, intent(out) :: conv_fail

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: common_fac
      real(dp) :: rad_levels(0:D2MAX2 - 1, 0:ECP_NLEVEL - 1)
      real(dp) :: rad_all(0:D3MAX2 - 1)
      real(dp) :: buffer(0:NFMAX*(2*ECP_LMAX + 1)*LC1MAX - 1)
      !! Written at n*ljlc1 + mj, n < ncart(li)*(2lc+1), mj < lj+lc+1: pyscf's
      !! nfi*(2*ECP_LMAX+1)*(lj+ECP_LMAX+1). Without the LC1MAX factor a d-d
      !! pair against a d channel already overran it (#194 review).
      integer :: lilc1, ljlc1, d2, im, mq, dlc
      integer :: lab, i, j, n, mi, mj, picked
      real(dp) :: acc
      logical :: failed

      lilc1 = li + lc + 1
      ljlc1 = lj + lc + 1
      d2 = lilc1*ljlc1
      dlc = 2*lc + 1
      im = ncart_l(li)*dlc
      mq = dlc*ljlc1
      common_fac = 16.0_dp*PI*PI

      conv_fail = 0
      rad_all(0:(li + lj + 1)*d2 - 1) = 0.0_dp
      do lab = 0, li + lj
         call ecp_radial_type2_lab_levels(rad_levels, D2MAX2, lab, lilc1, ljlc1, ur, ws, rs, radi, radj)
         call pick_converged_level_t2(D2MAX2, d2, rad_levels, picked, failed)
         if (failed) conv_fail = 1
         do n = 0, d2 - 1
            rad_all(lab*d2 + n) = rad_levels(n, picked)
         end do
      end do

      do i = 0, li
         do j = 0, lj
            do n = 0, im - 1
               do mj = 0, ljlc1 - 1
                  acc = 0.0_dp
                  do mi = 0, lilc1 - 1
                     acc = acc + rad_all((i + j)*d2 + mi*ljlc1 + mj) &
                           *angi(i*ncart_l(li)*dlc*lilc1 + n*lilc1 + mi)
                  end do
                  buffer(n*ljlc1 + mj) = acc
               end do
            end do
            do mj = 0, ncart_l(lj) - 1
               do mi = 0, ncart_l(li) - 1
                  acc = 0.0_dp
                  do n = 0, mq - 1
                     acc = acc + buffer(mi*mq + n) &
                           *angj(j*ncart_l(lj)*dlc*ljlc1 + mj*mq + n)
                  end do
                  gctr(mi + 1, mj + 1) = gctr(mi + 1, mj + 1) + common_fac*acc
               end do
            end do
         end do
      end do
   end subroutine ecp_type2_item

end module trc_ecp
