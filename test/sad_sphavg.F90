!
! Invariants of trc_sad.F90's spherical averaging (`average_atom_density`),
! checked directly rather than through the energy it leaves unchanged.
!
! PART 1 -- THE AVERAGING ALGEBRA, s/p/d/f, on a free O/cc-pVTZ atom (both
! its spherical AO space and its forced-Cartesian one). Reads only basis
! bookkeeping and trc_sph_tables' SPH_T (l = 0..4 regardless of
! TERCO_LMAX/ONE_E_LMAX), so it needs no integral kernel for any angular
! momentum and runs at every TERCO_LMAX. For every shell pair it builds TWO
! densities that share each block's trace but differ in its off-diagonal
! content -- standing in for the two different orthonormal bases a
! degenerate eigensolver could have returned for the identical physical
! atom -- and checks that `average_atom_density` sends both to the exact
! same result (determinism), that a same-l block becomes exactly
! (trace/(2l+1)) times the identity and a different-l one exactly zero
! (block proportionality), and that every block's trace survives (the
! electron-count argument in trc_sad.F90's own header, checked here as
! plain arithmetic rather than assumed).
!
! PART 2 -- THE FULL PIPELINE: `trc_sad_build` on the free O atom alone,
! Tr(D S) against Z = 8 with the REAL overlap from `trc_1e`. Gated on this
! build's one- and two-electron ceilings actually reaching f (TRC_MAXL,
! trc_tables' LMAX) -- the shipped TERCO_LMAX=2 configuration does not, and
! Part 1 above is what covers f there.
!
program sad_sphavg
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_1e, trc_bind_device, TRC_MAXL
   use trc_tables, only: TRC_ERI_LMAX => LMAX
   use trc_sph_tables, only: SPH_T
   use trc_sph, only: sph_cart_to_sph
   use trc_basis_json, only: trc_basis_from_json
   use trc_sad, only: trc_sad_build, average_atom_density
   use trc_error, only: error_t
   use trc_test_basis, only: read_xyz
   implicit none

   integer :: natm, dev
   integer, allocatable :: zint(:)
   real(dp), allocatable :: at_r(:, :)
   type(trc_basis_t) :: ab_sph, ab_cart
   type(error_t) :: err
   logical :: ok

   dev = trc_bind_device(0)
   call read_xyz('o_atom.xyz', natm, zint, at_r)
   ok = .true.

   call trc_basis_from_json('basis_sets/cc-pvtz.json', natm, zint, at_r, ab_sph, err)
   if (err%has_error()) call die(err)
   call trc_basis_from_json('basis_sets/cc-pvtz.json', natm, zint, at_r, ab_cart, err, force_cart=.true.)
   if (err%has_error()) call die(err)

   ! l_cap = 3 (s/p/d/f, this basis' own ceiling) for the spherical case:
   ! there is no contaminant there at any l, so `average_atom_density`'s
   ! simple replacement is exact end to end. l_cap = 1 (s/p only) for the
   ! Cartesian case: see `check_algebra`'s own header for why its
   ! determinism/block-proportionality fixture is deliberately NOT pushed
   ! past p there -- SAD's own open-shell atoms (O, C, N, F, ...) all have
   ! their Hund frontier at p, which is exactly where this holds.
   call check_algebra(ab_sph, .true., 3, ok)
   call check_algebra(ab_cart, .false., 1, ok)

   if (TRC_MAXL >= 3 .and. TRC_ERI_LMAX >= 3) then
      call check_pipeline(ok)
   else
      print '(a)', "sad_sphavg: this build's one-/two-electron ceiling is below f "// &
         "(TRC_MAXL, trc_tables' LMAX) -- skipping the full trc_sad_build(O/cc-pVTZ) check"
   end if

   call ab_sph%release(); call ab_cart%release()
   print '(a)', merge("sad_sphavg: PASS", "sad_sphavg: FAIL", ok)
   if (.not. ok) stop 1

contains

   subroutine die(e)
      type(error_t), intent(inout) :: e
      print '(a)', "sad_sphavg: "//e%get_message()
      stop 1
   end subroutine die

   !
   ! See this program's own header, PART 1. `l_cap` caps which shell pairs
   ! the synthetic fixture (and the expected-answer check) touches at all:
   ! `average_atom_density`'s simple same-trace replacement
   ! (`cart_to_sph_atom`/average/`sph_to_cart_atom`) is exact in the
   ! SPHERICAL representation for any l (no contaminant exists there), but
   ! in the CARTESIAN one it is exact only where T is the identity --
   ! l <= 1, where Cartesian AND spherical are the SAME representation (see
   ! `average_atom_density`'s header: T^T T is not the identity for l >= 2,
   ! so the averaging step itself, not just the electron count, runs in a
   ! metric-distorted space there). SAD's own open-shell atoms (O, C, N, F,
   ! ...) all have their Hund frontier at p, so capping the Cartesian
   ! check at l_cap = 1 tests exactly the regime that determinism bug
   ! lived in, without overclaiming determinism for a Cartesian d or f
   ! shell's own internal contaminant -- which `average_atom_density`
   ! still leaves untouched rather than discarding (see its header), so
   ! electron count survives there regardless; PART 2's real pipeline
   ! check below covers that block's Tr(D S) end to end instead.
   !
   subroutine check_algebra(ab, sph, l_cap, ok)
      type(trc_basis_t), intent(in) :: ab
      logical, intent(in) :: sph
      integer, intent(in) :: l_cap
      logical, intent(inout) :: ok
      real(dp), allocatable :: ds1(:, :), ds2(:, :), d1(:, :), d2(:, :), back1(:, :)
      integer :: n, ish, jsh, li, lj, bi, bj, nsi, nsj
      real(dp) :: maxdiff

      n = ab%nao_sph
      allocate (ds1(n, n), ds2(n, n))
      ds1 = 0.0_dp; ds2 = 0.0_dp
      do jsh = 1, ab%nshell
         do ish = 1, ab%nshell
            if (jsh < ish) cycle
            li = ab%sh_l(ish); lj = ab%sh_l(jsh)
            if (li > l_cap .or. lj > l_cap) cycle
            bi = ab%sh_ao_sph(ish); bj = ab%sh_ao_sph(jsh)
            nsi = 2*li + 1; nsj = 2*lj + 1
            if (li == lj) then
               call fill_block(nsi, base_of(ish, jsh), 0.37_dp, ds1(bi:bi + nsi - 1, bj:bj + nsi - 1))
               call fill_block(nsi, base_of(ish, jsh), -1.91_dp, ds2(bi:bi + nsi - 1, bj:bj + nsi - 1))
            else
               ! A different-l block: no physical density looks like this,
               ! only the zeroing is being exercised, so a plain nonzero
               ! constant (different in the two) is enough.
               ds1(bi:bi + nsi - 1, bj:bj + nsj - 1) = 0.731_dp
               ds2(bi:bi + nsi - 1, bj:bj + nsj - 1) = -0.218_dp
            end if
            if (jsh /= ish) then
               ds1(bj:bj + nsj - 1, bi:bi + nsi - 1) = transpose(ds1(bi:bi + nsi - 1, bj:bj + nsj - 1))
               ds2(bj:bj + nsj - 1, bi:bi + nsi - 1) = transpose(ds2(bi:bi + nsi - 1, bj:bj + nsj - 1))
            end if
         end do
      end do

      ! Both ds1 and ds2 have NO contaminant at all (every block came
      ! straight out of `fill_block`/the cross-l junk above): for the
      ! Cartesian case, d1 := local_sph_to_cart(ds1) puts them into the
      ! Cartesian representation EXACTLY at T's range, and nowhere else.
      if (sph) then
         allocate (d1(n, n), d2(n, n))
         d1 = ds1; d2 = ds2
      else
         allocate (d1(ab%nao, ab%nao), d2(ab%nao, ab%nao))
         call local_sph_to_cart(ab, ds1, d1)
         call local_sph_to_cart(ab, ds2, d2)
      end if

      call average_atom_density(ab, sph, d1)
      call average_atom_density(ab, sph, d2)

      maxdiff = maxval(abs(d1 - d2))
      print '(a,l1,a,es10.2)', "sad_sphavg: determinism (sph=", sph, ")  max|D1-D2| = ", maxdiff
      if (maxdiff > 1.0e-12_dp) ok = .false.

      ! Block proportionality: compare d1 against the KNOWN averaged
      ! answer (every same-l block exactly base_of(ish,jsh) times the
      ! identity, every different-l block exactly zero), folded to
      ! Cartesian by the SAME forward transform used to build the fixture
      ! -- not a round trip through `local_cart_to_sph`, which (like
      ! `cart_to_sph_atom` -- see `average_atom_density`'s header) is not
      ! the inverse of `local_sph_to_cart`/`sph_to_cart_atom` for l >= 2.
      allocate (back1(n, n))
      call build_expected(ab, l_cap, back1)
      if (sph) then
         maxdiff = maxval(abs(d1 - back1))
      else
         block
            real(dp), allocatable :: expect_cart(:, :)
            allocate (expect_cart(ab%nao, ab%nao))
            call local_sph_to_cart(ab, back1, expect_cart)
            maxdiff = maxval(abs(d1 - expect_cart))
            deallocate (expect_cart)
         end block
      end if
      print '(a,es10.2)', "sad_sphavg: block proportionality  max|D - expected| = ", maxdiff
      if (maxdiff > 1.0e-12_dp) ok = .false.
      deallocate (back1)
      deallocate (ds1, ds2, d1, d2)
   end subroutine check_algebra

   !
   ! The known answer: base_of(ish,jsh) times the identity for a same-l
   ! block, zero for a different-l one -- what `average_sph_blocks` must
   ! produce from the fixture `check_algebra` built.
   !
   subroutine build_expected(ab, l_cap, d)
      type(trc_basis_t), intent(in) :: ab
      integer, intent(in) :: l_cap
      real(dp), intent(out) :: d(ab%nao_sph, ab%nao_sph)
      integer :: ish, jsh, li, lj, bi, bj, nsi, r
      d = 0.0_dp
      do jsh = 1, ab%nshell
         do ish = 1, ab%nshell
            li = ab%sh_l(ish); lj = ab%sh_l(jsh)
            if (li /= lj .or. li > l_cap) cycle
            bi = ab%sh_ao_sph(ish); bj = ab%sh_ao_sph(jsh)
            nsi = 2*li + 1
            do r = 0, nsi - 1
               d(bi + r, bj + r) = base_of(ish, jsh)
            end do
         end do
      end do
   end subroutine build_expected

   !
   ! The trace each diagonal/off-diagonal block is filled with before
   ! averaging -- shared by the fixture above and the check below, so the
   ! two cannot silently drift apart.
   !
   pure real(dp) function base_of(ish, jsh)
      !! Symmetric in (ish, jsh): the (jsh, ish) block is the TRANSPOSE of
      !! (ish, jsh) (same trace), so both must look up the same base.
      integer, intent(in) :: ish, jsh
      base_of = real(min(ish, jsh), dp) + real(max(ish, jsh), dp)*0.1_dp
   end function base_of

   subroutine fill_block(ns, base, pert, blk)
      !! A symmetric ns x ns block with trace ns*base and an off-diagonal
      !! perturbation that depends only on `pert` -- so two calls with the
      !! same `base` and different `pert` share a trace but not a shape.
      !! `i + j` (not `i - j`) is the off-diagonal pattern deliberately:
      !! it is symmetric in (i, j) by construction, so this needs no
      !! separate symmetrisation step that could cancel it back out.
      integer, intent(in) :: ns
      real(dp), intent(in) :: base, pert
      real(dp), intent(out) :: blk(ns, ns)
      integer :: i, j
      do j = 1, ns
         do i = 1, ns
            if (i == j) then
               blk(i, j) = base
            else
               blk(i, j) = pert*real(i + j, dp)/real(ns, dp)
            end if
         end do
      end do
   end subroutine fill_block

   !
   ! Plain SPH_T contraction, written independently of trc_sad.F90's own
   ! (private) copy, so this test is not checking the implementation
   ! against itself -- see trc_sph.F90's header for the same formula's
   ! canonical, device-resident form. Only the forward (sph -> cart)
   ! direction is needed here: see `average_atom_density`'s header for why
   ! the reverse direction is not this transform's inverse for l >= 2, and
   ! so is never used as a check.
   !
   subroutine local_sph_to_cart(ab, ds, dc)
      type(trc_basis_t), intent(in) :: ab
      real(dp), intent(in)  :: ds(ab%nao_sph, ab%nao_sph)
      real(dp), intent(out) :: dc(ab%nao, ab%nao)
      integer :: ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj, p, q, r, s
      real(dp) :: acc, inner
      dc = 0.0_dp
      do jsh = 1, ab%nshell
         do ish = 1, ab%nshell
            li = ab%sh_l(ish); lj = ab%sh_l(jsh)
            nci = (li + 1)*(li + 2)/2; ncj = (lj + 1)*(lj + 2)/2
            nsi = 2*li + 1; nsj = 2*lj + 1
            ai = ab%sh_ao(ish); aj = ab%sh_ao(jsh)
            bi = ab%sh_ao_sph(ish); bj = ab%sh_ao_sph(jsh)
            do q = 1, ncj
               do p = 1, nci
                  acc = 0.0_dp
                  do s = 1, nsj
                     inner = 0.0_dp
                     do r = 1, nsi
                        inner = inner + SPH_T(p, r, li)*ds(bi + r - 1, bj + s - 1)
                     end do
                     acc = acc + inner*SPH_T(q, s, lj)
                  end do
                  dc(ai + p - 1, aj + q - 1) = acc
               end do
            end do
         end do
      end do
   end subroutine local_sph_to_cart

   !
   ! See this program's own header, PART 2.
   !
   subroutine check_pipeline(ok)
      logical, intent(inout) :: ok
      integer :: natm2, iv
      integer, allocatable :: z2(:)
      real(dp), allocatable :: r2(:, :), dguess(:, :), smat(:, :), tmat(:, :), vmat(:, :), ssph(:, :)
      type(trc_basis_t) :: b
      type(trc_pairlist_t) :: pl
      type(error_t) :: e
      real(dp) :: trace

      call read_xyz('o_atom.xyz', natm2, z2, r2)
      do iv = 1, 2
         if (iv == 1) then
            call trc_basis_from_json('basis_sets/cc-pvtz.json', natm2, z2, r2, b, e)
         else
            call trc_basis_from_json('basis_sets/cc-pvtz.json', natm2, z2, r2, b, e, force_cart=.true.)
         end if
         if (e%has_error()) call die(e)
         call trc_sad_build(b, dguess, e, verbose=.true.)
         if (e%has_error()) call die(e)
         call b%to_device()
         call pl%build(b, 1.0e-12_dp); call pl%to_device()
         allocate (smat(b%nao, b%nao), tmat(b%nao, b%nao), vmat(b%nao, b%nao))
         call trc_1e(b, pl, smat, tmat, vmat)
         if (b%spherical) then
            allocate (ssph(b%nao_sph, b%nao_sph))
            call sph_cart_to_sph(b, smat, ssph)
            trace = sum(dguess*ssph)
            deallocate (ssph)
         else
            trace = sum(dguess*smat)
         end if
         print '(a,l1,a,f18.12)', "sad_sphavg: free O/cc-pVTZ  spherical=", b%spherical, "  Tr(D S) = ", trace
         if (abs(trace - 8.0_dp) > 1.0e-10_dp) ok = .false.
         deallocate (smat, tmat, vmat, dguess)
         call pl%release(); call b%release()
      end do
   end subroutine check_pipeline

end program sad_sphavg
