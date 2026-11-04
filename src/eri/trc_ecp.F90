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
! slice a strided view of the SAME cached array"); K2 is one item per
! (orbital shell pair, ECP atom), summing that atom's own channels and
! writing its own exclusive output slot (no atomics); K3 reduces K2's
! per-atom slots into the Cartesian matrix.
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
!
module trc_ecp
   use pic_types, only: dp
   use trc_basis_types, only: trc_basis_t
   use trc_ecp_types, only: trc_ecp_t
   use trc_ecp_quadrature, only: ECP_LMAX, ECP_LEVEL0, ECP_LEVEL_MAX, ECP_NRS, &
                                 ecp_quadrature_gauss_chebyshev
   use trc_ecp_angular, only: ECP_ORB_LMAX, CART_POW_Y, CART_POW_Z, &
                              ecp_type1_rad_ang, ecp_type2_facs_ang
   use trc_ecp_radial, only: ecp_radial_rad_part, ecp_radial_type1_static_facs, &
                             ecp_radial_type1_rad_part, ecp_radial_type2_facs_rad
   use trc_sph, only: sph_cart_to_sph
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private
   public :: trc_ecp_matrix

   integer, parameter :: NFMAX = (ECP_ORB_LMAX + 1)*(ECP_ORB_LMAX + 2)/2     !! 10
   integer, parameter :: NPMAX = 16       !! primitives per orbital shell or ECP row
   integer, parameter :: D1MAX = 2*ECP_ORB_LMAX + 1                          !! 7
   integer, parameter :: D2MAX = D1MAX*D1MAX                                 !! 49
   integer, parameter :: D3MAX = D2MAX*D1MAX                                 !! 343
   integer, parameter :: LC1MAX = ECP_ORB_LMAX + ECP_LMAX + 1                 !! 9
   integer, parameter :: D2MAX2 = LC1MAX*LC1MAX                              !! 81
   integer, parameter :: D3MAX2 = D1MAX*D2MAX2                               !! 567
   integer, parameter :: SLOT_SIZE = NFMAX*NFMAX                                  !! 100: one (shell pair, atom) slot

   !> A contribution below this, to a matrix whose real entries run from
   !! ~1e-3 to a few hundred (measured: Xe/def2-SVP's U diagonal reaches
   !! 206), can never matter past the 1e-9 tolerance anything here is
   !! validated to. Used to stop the adaptive radial loop from chasing a
   !! RELATIVE 1e-12 match on a term that is already ABSOLUTELY negligible
   !! -- which a large separation (ecp_ints_ih_far's 25 bohr) or a tight
   !! core primitive pair reaches often, and which the relative test alone
   !! can fail to confirm within this port's ECP_LEVEL_MAX even though the
   !! term's actual contribution to U is already settled at zero.
   real(dp), parameter :: ECP_NEGLIGIBLE = 1.0e-20_dp

contains

   pure integer function ncart_l(l) result(n)
      !$acc routine seq
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)/2
   end function ncart_l

   !
   ! The whole ECP one-electron matrix, Cartesian internally and folded to
   ! spherical if `b%spherical`. `u` comes back allocated to b%nao_sph (or
   ! b%nao). `ecp%has_ecp` false gives back a zero matrix of the right size.
   !
   subroutine trc_ecp_matrix(b, ecp, u, error)
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), allocatable, intent(out) :: u(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: mcart(:, :)

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
         if (ecp%nrow > 0) then
            if (maxval(ecp%ecp_np) > NPMAX) then
               call error%set(ERROR_VALIDATION, "trc_ecp_matrix: an ECP row has more "// &
                              "primitives than this port's fixed NPMAX ceiling")
               return
            end if
         end if
         call ecp_matrix_cart(b, ecp, mcart, error)
         if (error%has_error()) return
      end if

      if (b%spherical) then
         allocate (u(b%nao_sph, b%nao_sph))
         call sph_cart_to_sph(b, mcart, u)
      else
         allocate (u(b%nao, b%nao))
         u = mcart
      end if
      deallocate (mcart)
   end subroutine trc_ecp_matrix

   !
   ! The Cartesian build: groups the ECP rows into (atom, l) channels, then
   ! K1/K2/K3 as this module's own header describes.
   !
   subroutine ecp_matrix_cart(b, ecp, mcart, error)
      type(trc_basis_t), intent(in) :: b
      type(trc_ecp_t), intent(in) :: ecp
      real(dp), intent(inout) :: mcart(b%nao, b%nao)
      type(error_t), intent(inout) :: error

      real(dp) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1)
      real(dp), allocatable :: ur_cache(:, :)
      integer, allocatable :: grp_atom(:), grp_l(:), grp_row0(:), grp_row1(:)
      integer, allocatable :: atm_ecp(:), atm_grp0(:), atm_grp1(:)
      integer, allocatable :: pr_i(:), pr_j(:)
      real(dp), allocatable :: block_buffer(:, :, :)
      integer, allocatable :: conv_fail(:, :)
      integer :: ngroup, natm_ecp, npair
      integer :: i, j, p, g, a, nfail

      call build_groups(ecp, grp_atom, grp_l, grp_row0, grp_row1, ngroup)
      if (ngroup == 0) return
      call build_atom_list(grp_atom, ngroup, atm_ecp, atm_grp0, atm_grp1, natm_ecp)

      call ecp_quadrature_gauss_chebyshev(rs, ws, ECP_NRS)
      allocate (ur_cache(0:ECP_NRS - 1, ngroup))

      !$acc enter data copyin(rs, ecp, ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c) &
      !$acc            copyin(grp_row0, grp_row1) create(ur_cache)
      do concurrent(g=1:ngroup)
         call ecp_radial_rad_part(ur_cache(:, g), rs, grp_row0(g), grp_row1(g), &
                           ecp%nrow, ecp%maxnp, ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c)
      end do
      !$acc exit data delete(ecp%ecp_rpow, ecp%ecp_np, ecp%ecp_e, ecp%ecp_c, grp_row0, grp_row1, ecp)

      npair = b%nshell*(b%nshell + 1)/2
      allocate (pr_i(npair), pr_j(npair))
      p = 0
      do i = 1, b%nshell
         do j = 1, i
            p = p + 1
            pr_i(p) = i
            pr_j(p) = j
         end do
      end do

      allocate (block_buffer(SLOT_SIZE, npair, max(natm_ecp, 1)), conv_fail(npair, max(natm_ecp, 1)))
      block_buffer = 0.0_dp
      conv_fail = 0

      !$acc enter data copyin(ws, pr_i, pr_j, b, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%sh_ao) &
      !$acc            copyin(grp_l, atm_ecp, atm_grp0, atm_grp1, b%at_r) &
      !$acc            create(block_buffer, conv_fail)
      if (natm_ecp > 0) then
         do concurrent(p=1:npair, a=1:natm_ecp)
            call ecp_pair_atom_item(pr_i(p), pr_j(p), atm_ecp(a), atm_grp0(a), atm_grp1(a), &
                                    ngroup, grp_l, &
                                    b%nshell, b%sh_l, b%sh_np, b%maxnp, b%sh_e, b%sh_c, b%sh_r, &
                                    b%natm, b%at_r, &
                                    rs, ws, ur_cache, &
                                    block_buffer(:, p, a), conv_fail(p, a))
         end do
      end if
      !$acc update self(block_buffer, conv_fail)
      !$acc exit data delete(ws, pr_i, pr_j, b%sh_l, b%sh_np, b%sh_e, b%sh_c, b%sh_r, b%sh_ao, b%at_r, b, &
      !$acc                  grp_l, atm_ecp, atm_grp0, atm_grp1, rs, ur_cache, block_buffer, conv_fail)

      nfail = sum(conv_fail)
      if (nfail > 0) then
         call error%set(ERROR_VALIDATION, "trc_ecp_matrix: the ECP radial quadrature did not "// &
                        "converge within this port's ECP_LEVEL_MAX for at least one shell "// &
                        "pair/channel; raise trc_ecp_quadrature's ECP_LEVEL_MAX")
         return
      end if

      ! K3: reduce the per-atom slots and write the Cartesian block (and
      ! its symmetric mirror -- U is symmetric, scalar ECPs are real, and
      ! only i>=j pairs were built).
      do p = 1, npair
         call deposit_block(b, mcart, pr_i(p), pr_j(p), block_buffer(:, p, :), natm_ecp)
      end do

      deallocate (ur_cache, grp_atom, grp_l, grp_row0, grp_row1)
      deallocate (atm_ecp, atm_grp0, atm_grp1, pr_i, pr_j, block_buffer, conv_fail)
   end subroutine ecp_matrix_cart

   subroutine deposit_block(b, mcart, i, j, slots, natm_ecp)
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(inout) :: mcart(b%nao, b%nao)
      integer, intent(in) :: i, j, natm_ecp
      real(dp), intent(in) :: slots(SLOT_SIZE, max(natm_ecp, 1))
      integer :: li, lj, nfi, nfj, ai, aj, mi, mj, a
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
            mcart(ai + mi - 1, aj + mj - 1) = mcart(ai + mi - 1, aj + mj - 1) + acc
            if (i /= j) mcart(aj + mj - 1, ai + mi - 1) = mcart(aj + mj - 1, ai + mi - 1) + acc
         end do
      end do
   end subroutine deposit_block

   !
   ! Group the ECP rows into (atom, l) channels: consecutive rows sharing
   ! both -- see trc_ecp_types.F90 for why the reader always emits a
   ! channel's rows consecutively.
   !
   subroutine build_groups(ecp, grp_atom, grp_l, grp_row0, grp_row1, ngroup)
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

   !
   ! K2's item: one (orbital shell pair, ECP atom). Sums every one of that
   ! atom's channels (local and projected) into its OWN output slot,
   ! `slot` -- SLOT_SIZE = NFMAX*NFMAX reals, column-major (mi, mj) with mi the
   ! FASTER index, zero-padded past the pair's actual (nfi, nfj). No
   ! atomics: nothing else ever writes this slot.
   !
   pure subroutine ecp_pair_atom_item(ish, jsh, atm_id, g0, g1, &
                                      ngroup, grp_l, &
                                      nshell, sh_l, sh_np, maxnp, sh_e, sh_c, sh_r, &
                                      natm, at_r, &
                                      rs, ws, ur_cache, &
                                      slot, conv_fail)
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
                                ri, rj, rc, rs, ws, ur_cache(:, g), local_fail)
         else
            call ecp_type2_item(gctr, li, lj, npi, npj, &
                                sh_e(1:npi, ish), sh_c(1:npi, ish), sh_e(1:npj, jsh), sh_c(1:npj, jsh), &
                                ri, rj, rc, grp_l(g), rs, ws, ur_cache(:, g), local_fail)
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

   !
   ! The type-1 (local) term for one (shell pair, ECP channel). Per-primitive-
   ! pair convergence, exactly as libfint's `ecp_type1_cart`: each (ip,jp)'s
   ! own radial integral is refined level by level and the angular
   ! contraction runs only once, after every pair has converged (or the
   ! level ceiling is hit, which sets `conv_fail`).
   !
   pure subroutine ecp_type1_item(gctr, li, lj, npi, npj, ei, ci, ej, cj, ri, rj, rc, &
                                  rs, ws, ur, conv_fail)
      !$acc routine seq
      integer, intent(in) :: li, lj, npi, npj
      real(dp), intent(in) :: ei(npi), ci(npi), ej(npj), cj(npj)
      real(dp), intent(in) :: ri(3), rj(3), rc(3)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(inout) :: gctr(NFMAX, NFMAX)
      integer, intent(out) :: conv_fail

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: rca(3), rcb(3), rij(3)
      real(dp) :: cei(NPMAX), cej(NPMAX)
      real(dp) :: rad_all(0:NPMAX*NPMAX*D2MAX - 1), plast(0:D2MAX - 1)
      real(dp) :: ur_level(0:ECP_NRS - 1), rs_level(0:ECP_NRS - 1)
      real(dp) :: rad_ang(0:D3MAX - 1), rad_ang_all(0:D3MAX - 1)
      real(dp) :: ifac(0:NFMAX*D1MAX*D1MAX*D1MAX - 1), jfac(0:NFMAX*D1MAX*D1MAX*D1MAX - 1)
      integer :: lilj1, d1, d2, d3, di1, dj1, pr
      integer :: ip, jp, mi, mj, n, i1, i2, i3, j1, j2, j3, ix, iy, iz, jx, jy, jz
      integer :: level, start, step, nrs0, nrs
      real(dp) :: wtscale, r2ca, r2cb, fac
      integer :: converged(NPMAX, NPMAX)
      logical :: pair_close, all_conv

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

      converged = 0
      rad_all(0:npi*npj*d2 - 1) = 0.0_dp
      all_conv = .true.
      nrs0 = 2**ECP_LEVEL0 - 1
      step = 2**(ECP_LEVEL_MAX - ECP_LEVEL0)
      start = step - 1
      wtscale = real(step, dp)

      do level = ECP_LEVEL0, ECP_LEVEL_MAX
         nrs = nrs0
         do n = 0, nrs - 1
            ur_level(n) = ur(start + n*step)*ws(start + n*step)*wtscale
            ! rs is densified the SAME way ur is: ecp_radial_type1_rad_part reads
            ! both at plain index n, not at this level's own (start, step)
            ! -- passing (start, step) with the UNdensified rs here (as an
            ! earlier version of this routine did) paired ur_level(n) (this
            ! level's point) against rs(start + n*step), which is a
            ! DIFFERENT, uninvolved grid point once step > 1: the two must
            ! read the SAME densified index.
            rs_level(n) = rs(start + n*step)
         end do
         all_conv = .true.
         do ip = 1, npi
            do jp = 1, npj
               ! A primitive pair whose centre-shift prefactor has already
               ! decayed past this floor contributes negligibly to gctr
               ! (fac*rad_ang below) however rad_all's own magnitude looks:
               ! `ecp_radial_type1_rad_part`'s exponent-shift (complete the
               ! square around kaij, rather than evaluate exp(-aij*rn^2)
               ! directly) trades a huge-but-finite intermediate for one
               ! that would otherwise underflow outright, and that
               ! intermediate's own growth tracks the SAME exponents this
               ! prefactor decays by, so it never overtakes a floor this
               ! far below any value that could matter (measured: a
               ! prefactor of 8.4e-251 next to a "val" of 8.6e+187, net
               ! contribution ~1e-63) -- so this never needs refining, and
               ! never blocks all_conv. A tight core primitive far from the
               ! ECP centre (ecp_ints_ih_far's 25 bohr I-H, ecp_ints_csh's
               ! Cs 1s next to the same-atom local channel) reaches this
               ! routinely.
               if (converged(ip, jp) >= 2 .or. abs(cei(ip)*cej(jp)) < 1.0e-150_dp) cycle
               pr = ((ip - 1)*npj + (jp - 1))*d2
               do n = 0, d2 - 1
                  plast(n) = rad_all(pr + n)
                  rad_all(pr + n) = rad_all(pr + n)*0.5_dp
               end do
               rij = ei(ip)*rca + ej(jp)*rcb
               call ecp_radial_type1_rad_part(rad_all(pr:pr + d2 - 1), li + lj, &
                                      sqrt(sum(rij**2))*2.0_dp, ei(ip) + ej(jp), &
                                      ur_level, rs_level, 0, 1, nrs)
               pair_close = .true.
               do n = 0, d2 - 1
                  if (.not. close_enough(plast(n), rad_all(pr + n))) then
                     pair_close = .false.
                     exit
                  end if
               end do
               if (pair_close) then
                  converged(ip, jp) = converged(ip, jp) + 1
                  if (converged(ip, jp) < 2) all_conv = .false.
               else
                  converged(ip, jp) = 0
                  all_conv = .false.
               end if
            end do
         end do
         if (all_conv) exit
         nrs0 = 2**level - 1
         step = 2**(ECP_LEVEL_MAX - level)
         start = (start - 1)/2
         wtscale = wtscale*0.5_dp
      end do
      conv_fail = 0
      if (.not. all_conv) conv_fail = 1
      rad_ang_all = 0.0_dp
      do ip = 1, npi
         do jp = 1, npj
            pr = ((ip - 1)*npj + (jp - 1))*d2
            rij = ei(ip)*rca + ej(jp)*rcb
            call ecp_type1_rad_ang(rad_ang, li + lj, rij, rad_all(pr:pr + d2 - 1))
            fac = cei(ip)*cej(jp)
            do n = 0, d3 - 1
               rad_ang_all(n) = rad_ang_all(n) + fac*rad_ang(n)
            end do
         end do
      end do

      call ecp_radial_type1_static_facs(ifac, li, rca)
      call ecp_radial_type1_static_facs(jfac, lj, rcb)

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

   !
   ! The type-2 (semi-local) term for one (shell pair, ECP channel of
   ! angular momentum lc). No nc/ic loop (trc_basis_t never holds a general
   ! contraction), so the convergence state is indexed by `lab` alone, and
   ! rad_all carries no (ic,jc) offset.
   !
   pure subroutine ecp_type2_item(gctr, li, lj, npi, npj, ei, ci, ej, cj, ri, rj, rc, lc, &
                                  rs, ws, ur, conv_fail)
      !$acc routine seq
      integer, intent(in) :: li, lj, npi, npj, lc
      real(dp), intent(in) :: ei(npi), ci(npi), ej(npj), cj(npj)
      real(dp), intent(in) :: ri(3), rj(3), rc(3)
      real(dp), intent(in) :: rs(0:ECP_NRS - 1), ws(0:ECP_NRS - 1), ur(0:ECP_NRS - 1)
      real(dp), intent(inout) :: gctr(NFMAX, NFMAX)
      integer, intent(out) :: conv_fail

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: rca(3), rcb(3)
      real(dp) :: dca, dcb, common_fac
      real(dp) :: radi(0:ECP_NRS*LC1MAX - 1), radj(0:ECP_NRS*LC1MAX - 1)
      real(dp) :: ur_level(0:ECP_NRS - 1)
      real(dp) :: rad_all(0:D3MAX2 - 1), plast(0:D2MAX2 - 1)
      real(dp) :: angi(0:(ECP_ORB_LMAX + 1)*NFMAX*(2*ECP_LMAX + 1)*LC1MAX - 1)
      real(dp) :: angj(0:(ECP_ORB_LMAX + 1)*NFMAX*(2*ECP_LMAX + 1)*LC1MAX - 1)
      real(dp) :: buffer(0:NFMAX*(2*ECP_LMAX + 1) - 1)
      real(dp) :: rpow(0:ECP_NRS - 1)
      integer :: lilj1, lilc1, ljlc1, d2, im, mq, dlc
      integer :: level, nrs0, start, step, nrs, lab
      integer :: i, j, n, mi, mj
      real(dp) :: wtscale, s, acc
      logical :: all_conv, pair_close
      integer :: converged(0:D1MAX - 1)

      lilj1 = li + lj + 1
      lilc1 = li + lc + 1
      ljlc1 = lj + lc + 1
      d2 = lilc1*ljlc1
      dlc = 2*lc + 1
      im = ncart_l(li)*dlc
      mq = dlc*ljlc1
      common_fac = 16.0_dp*PI*PI

      rca = rc - ri
      rcb = rc - rj
      dca = sqrt(sum(rca**2))
      dcb = sqrt(sum(rcb**2))

      converged = 0
      rad_all(0:lilj1*d2 - 1) = 0.0_dp
      all_conv = .true.
      nrs0 = 2**ECP_LEVEL0 - 1
      step = 2**(ECP_LEVEL_MAX - ECP_LEVEL0)
      start = step - 1
      wtscale = real(step, dp)

      do level = ECP_LEVEL0, ECP_LEVEL_MAX
         nrs = nrs0
         call ecp_radial_type2_facs_rad(radi, li, npi, ei, ci, lc, dca, rs, start, step, nrs)
         call ecp_radial_type2_facs_rad(radj, lj, npj, ej, cj, lc, dcb, rs, start, step, nrs)
         do n = 0, nrs - 1
            ur_level(n) = ur(start + n*step)*ws(start + n*step)*wtscale
         end do

         all_conv = .true.
         do lab = 0, li + lj
            if (converged(lab) >= 2) cycle
            ! r^lab at this level's own grid points, computed fresh per lab
            ! (rather than an incremental ladder) because an already-
            ! converged lab's `cycle` above would otherwise skip the
            ! ladder's advance and desynchronise every lab after it.
            do n = 0, nrs - 1
               rpow(n) = rs(start + n*step)**lab
            end do
            do n = 0, d2 - 1
               plast(n) = rad_all(lab*d2 + n)
            end do
            do i = 0, lilc1 - 1
               do j = 0, ljlc1 - 1
                  ! Halved, then the new (finer) points' contribution added
                  ! -- the same nested-quadrature refinement type 1's own
                  ! rad_all halving is, and for the same reason: the
                  ! previous level's points reappear in this one at twice
                  ! the weight.
                  s = rad_all(lab*d2 + i*ljlc1 + j)*0.5_dp
                  do n = 0, nrs - 1
                     s = s + ur_level(n)*rpow(n)*radi(n*lilc1 + i)*radj(n*ljlc1 + j)
                  end do
                  rad_all(lab*d2 + i*ljlc1 + j) = s
               end do
            end do
            ! An absolutely negligible lab-block (every entry already below
            ! ECP_NEGLIGIBLE) is accepted without the relative test:
            ! common_fac (16*pi^2) and the angular contraction below can
            ! only scale it further, never promote it into significance,
            ! and a term already this small can fail the 1e-12 RELATIVE
            ! test indefinitely (ecp_ints_ih_far's 25 bohr I-H) without
            ! ever being wrong to within this port's own tolerance.
            if (maxval(abs(rad_all(lab*d2:lab*d2 + d2 - 1))) < ECP_NEGLIGIBLE) then
               pair_close = .true.
            else
               pair_close = .true.
               do n = 0, d2 - 1
                  if (.not. close_enough(plast(n), rad_all(lab*d2 + n))) then
                     pair_close = .false.
                     exit
                  end if
               end do
            end if
            if (pair_close) then
               converged(lab) = converged(lab) + 1
               if (converged(lab) < 2) all_conv = .false.
            else
               converged(lab) = 0
               all_conv = .false.
            end if
         end do

         if (all_conv) exit
         nrs0 = 2**level - 1
         step = 2**(ECP_LEVEL_MAX - level)
         start = (start - 1)/2
         wtscale = wtscale*0.5_dp
      end do
      conv_fail = 0
      if (.not. all_conv) conv_fail = 1
      call ecp_type2_facs_ang(angi, li, lc, rca)
      call ecp_type2_facs_ang(angj, lj, lc, rcb)

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
