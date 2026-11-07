!
! SYMMETRY-UNIQUE PAIR LIST: design prototype (docs/SYMMETRY_UNIQUE_PAIRS.md).
! Registered DISABLED; not part of the suite. CPU, plain Fortran.
!
! For one molecule (cc-pVDZ or any basis JSON, Cartesian caller space) it
! takes terco's own symmetry-aware engine (`eri%build(..., sym=geom)`: bins
! closed under the group, representatives first) and answers, with numbers,
! the questions of the design study:
!
!  1. COUNTS (integral free). Pairs, orbits (representatives), primitive
!     pairs, and the number of quartet orbits the skeleton has to evaluate,
!     against the number a "representatives on both sides" enumeration
!     reaches. The orbit count is exact: Burnside on the pair-of-pairs
!     (G acts on P x P: N = (1/|G|) sum_g fix(g)^2; unordered {p,q} with the
!     bra<->ket swap: N = (1/2|G|) sum_g (fix(g)^2 + npair)).
!  2. MEMORY MODEL. Bytes of the pair data held by the merged scheme (2a),
!     taken from the sizes of the production arrays, and by the unique-pair
!     variants: U1 (primitive pairs and Schwarz bounds for representatives
!     only, bin tables expanded from them) and U2 (everything on the device
!     representatives only, the image read in the kernel).
!  3. RECONSTRUCTION. For every pair p of the bins that is not a
!     representative, the primitive-pair data is rebuilt from its
!     representative's alone: zeta, coefficient and exponent unchanged,
!     the product centre P and the shell centres A, B sign-flipped by the
!     operation's mask. Compared with the directly built data (max
!     deviation, primitive counts, ordering).
!  4. EXACTNESS. The engine's pair data of every non-representative bin pair
!     is replaced by the reconstruction, the bin-interleaved tables are
!     rebuilt from it, and the symmetrized skeleton G is compared with the
!     full G (gate 5e-11*max(1,|G|)): J+K and K only, a random invariant
!     density.
!
! usage: spike_unique <xyz> <basis.json> <count|verify> [thresh]
!
! The counts and the reconstruction check read the FULL pair data the merged
! build forms: run with TRC_UNIQUE_PAIRS=0 (the engine's default is now the
! unique-pair build, which forms the representatives' only).
!
module spike_unique_mod
   use pic_types, only: dp, int64
   use trc_eri_data, only: trc_shell_pairs_t
   implicit none
contains

   pure integer function flip_sign(mask, d) result(s)
      integer, intent(in) :: mask, d
      s = merge(-1, 1, btest(mask, d - 1))
   end function flip_sign

   subroutine copy_pairs(src, dst)
      type(trc_shell_pairs_t), intent(in) :: src
      type(trc_shell_pairs_t), intent(out) :: dst
      dst = src
   end subroutine copy_pairs

   subroutine symmetrize_ident(n, x)
      integer, intent(in) :: n
      real(dp), intent(inout) :: x(n, n)
      x = 0.5_dp*(x + transpose(x))
   end subroutine symmetrize_ident

end module spike_unique_mod

program spike_unique
   use pic_types, only: dp, int64
   use spike_unique_mod
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_bins, only: bins_build_pair_quartets
   use trc_basis_json, only: trc_basis_from_json
   use trc_bin_tiles, only: bin_tiles_set_enabled
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_detect
   use trc_sym_pairs, only: trc_sym_pairs_t, sym_pairs_build, sym_pairs_tri
   use trc_sym_scf, only: trc_sym_scf_t
   use trc_test_basis, only: read_xyz
   implicit none

   character(len=256) :: xyz, bfile, mode, arg
   real(dp) :: thresh
   integer :: nat, i, j, p, g, nop, np, nsh, dev, nrep
   integer, allocatable :: elem(:)
   real(dp), allocatable :: r(:, :)
   type(trc_basis_t) :: bas
   type(trc_sym_geom_t) :: geom
   type(trc_sym_ao_t) :: ao
   type(trc_eri_t) :: eri
   type(trc_sym_pairs_t) :: tab
   type(error_t) :: err
   logical :: verify
   integer, allocatable :: nfix(:)
   real(dp) :: t0, t1

   dev = trc_bind_device(0)
   call get_command_argument(1, xyz)
   call get_command_argument(2, bfile)
   call get_command_argument(3, mode)
   thresh = 1.0e-10_dp
   call get_command_argument(4, arg)
   if (len_trim(arg) > 0) read (arg, *) thresh
   verify = index(mode, 'verify') > 0

   call read_xyz(trim(xyz), nat, elem, r)
   geom = sym_detect(r, elem, err, basis_id=elem)
   if (err%has_error()) stop 'sym_detect failed'
   call trc_basis_from_json(trim(bfile), nat, elem, geom%r_std, bas, err, force_cart=.true.)
   if (err%has_error()) stop 'basis failed'
   bas%spherical = .false.
   call bas%to_device()
   print '(a,a,a,i0,a,a)', 'molecule ', trim(xyz), '  natm ', nat, '  group ', trim(geom%group%name)

   call cpu_time(t0)
   call eri%build(bas, thresh, sym=geom)
   call cpu_time(t1)
   print '(a,f8.3,a,5f8.3)', 'build (sym) cpu s ', t1 - t0, '  stages ', eri%t_build
   if (.not. eri%orb%built) stop 'no orbit tables (C1?)'
   nop = eri%orb%nop
   np = eri%bins%npair
   nsh = eri%nbas

   tab = sym_pairs_build(eri%orb%shmap, eri%bins%sp_i(1:np), eri%bins%sp_j(1:np), eri%bins%sp_q(1:np), err)
   if (err%has_error()) stop 'sym_pairs_build failed'
   nrep = tab%nrep
   print '(a,l1,a,i0,a,i0,a,i0,a,i0)', 'engine: unique-pair build ', eri%upair, ';  shell pairs formed ', eri%upair_nshell_pairs, &
      ' of ', nsh*nsh, ';  primitive pairs stored ', eri%upair_npp_stored, ' for ', eri%upair_npp_bins
   if (eri%dec%active) then
      print '(a,i0,a,f10.2,a,f10.2)', 'hgp_build_pairs keep scratch (max primitives ', size(eri%pb%sh_e, 1), '): MB all pairs ', &
         real(size(eri%pb%sh_e, 1), dp)**2*nsh*nsh/1.0e6_dp, '  wanted only ', &
         real(size(eri%pb%sh_e, 1), dp)**2*eri%upair_nshell_pairs/1.0e6_dp
   else
      print '(a,i0,a,f10.2,a,f10.2)', 'hgp_build_pairs keep scratch (max primitives ', size(bas%sh_e, 1), '): MB all pairs ', &
         real(size(bas%sh_e, 1), dp)**2*nsh*nsh/1.0e6_dp, '  wanted only ', &
         real(size(bas%sh_e, 1), dp)**2*eri%upair_nshell_pairs/1.0e6_dp
   end if
   call counts()
   if (verify) call reconstruct_and_verify()

contains

   integer function pkey(a, b) result(k)
      integer, intent(in) :: a, b
      k = (a - 1)*nsh + b
   end function pkey

   subroutine counts()
      integer(int64) :: npp_all, npp_rep, npp_rep_or, npp_stored
      integer :: cnt(max(nop, 1)), k, h
      integer(int64) :: sumfix2, nq_ord, nq_unord, nrr
      integer(int64) :: b_pairs_full, b_pairs_u1, b_q_full, b_q_u2, b_bins_full, b_bins_u2, nqp_rep
      logical, allocatable :: need(:)
      integer :: pu, ndist
      real(dp), parameter :: MB = 1.0e6_dp

      cnt = tab%size_counts()
      print '(a,i0,a,i0,a,i0,a,i0)', 'nop ', nop, '  shells ', nsh, '  pairs (bins) ', np, '  orbits ', nrep
      print '(a,8i9)', 'orbits by size (1..nop): ', cnt(1:nop)
      ndist = 0
      do g = 1, nop
         if (.not. any([(all(eri%orb%shmap(:, g) == eri%orb%shmap(:, h)), h=1, g - 1)])) ndist = ndist + 1
      end do
      print '(a,i0,a,i0)', 'effective order (distinct shell permutations) ', ndist, ' of ', nop

      ! primitive pairs
      npp_all = 0
      npp_rep = 0
      do p = 1, np
         npp_all = npp_all + eri%pairs%n(pkey(eri%bins%sp_i(p), eri%bins%sp_j(p)))
         if (tab%rep(p) == p) npp_rep = npp_rep + eri%pairs%n(pkey(eri%bins%sp_i(p), eri%bins%sp_j(p)))
      end do
      ! orientations the images need: the representative's own, and its reverse
      ! when an image is stored reversed
      allocate (need(nsh*nsh))
      need = .false.
      do p = 1, np
         if (tab%rep(p) /= p) cycle
         need(pkey(eri%bins%sp_i(p), eri%bins%sp_j(p))) = .true.
      end do
      do p = 1, np
         if (tab%rep(p) == p) cycle
         call source_of(p, pu, h)
         need(pu) = .true.
      end do
      npp_rep_or = 0
      do k = 1, nsh*nsh
         if (need(k)) npp_rep_or = npp_rep_or + eri%pairs%n(k)
      end do
      npp_stored = sum(int(eri%pairs%n, int64))
      print '(a,i0,a,i0,a,i0,a,i0)', 'primitive pairs: production pairs object (both orientations) ', npp_stored, &
         ';  bins, all pairs ', npp_all, ';  bins, representatives ', npp_rep, ';  unique list incl. reversed orientations ', npp_rep_or

      ! quartet orbits (Burnside on the pair-of-pairs)
      allocate (nfix(nop))
      nfix = 0
      do g = 1, nop
         do p = 1, np
            if (tab%img(p, g) == p) nfix(g) = nfix(g) + 1
         end do
      end do
      sumfix2 = 0
      do g = 1, nop
         sumfix2 = sumfix2 + int(nfix(g), int64)**2
      end do
      nq_ord = sumfix2/nop
      nq_unord = (sumfix2 + int(nop, int64)*np)/(2*int(nop, int64))
      nrr = int(nrep, int64)*(nrep + 1)/2
      print '(a,i0)', 'quartets, full unordered pair-of-pairs (npair(npair+1)/2): ', int(np, int64)*(np + 1)/2
      print '(a,i0)', 'quartet orbits, unordered {ab,cd} (Burnside):            ', nq_unord
      print '(a,i0,a,f8.5)', 'quartets with both pairs representatives:                ', nrr, &
         '   fraction of the orbits ', real(nrr, dp)/real(nq_unord, dp)
      print '(a,f8.3,a,f8.3)', 'full/orbits (ideal reduction of a per-quartet scheme) ', &
         real(int(np, int64)*(np + 1)/2, dp)/real(nq_unord, dp), ';  full/(rep x rep) ', real(int(np, int64)*(np + 1)/2, dp)/real(nrr, dp)

      ! memory model, bytes (production sizes)
      nqp_rep = 0
      do p = 1, np
         if (tab%rep(p) == p) nqp_rep = nqp_rep + eri%pairs%n(pkey(eri%bins%sp_i(p), eri%bins%sp_j(p)))
      end do
      b_pairs_full = 8_int64*size(eri%pairs%p) + 8_int64*size(eri%pairs%r) + 8_int64*size(eri%pairs%ra) &
                     + 8_int64*size(eri%pairs%rb) + 8_int64*size(eri%pairs%c) + 8_int64*size(eri%pairs%a) &
                     + 4_int64*size(eri%pairs%off) + 4_int64*size(eri%pairs%n)
      b_pairs_u1 = 88_int64*npp_rep_or + 8_int64*nsh*nsh
      b_q_full = 8_int64*(size(eri%quartets%z) + size(eri%quartets%c) + size(eri%quartets%a) + size(eri%quartets%px) &
                          + size(eri%quartets%py) + size(eri%quartets%pz)) &
                 + 8_int64*(size(eri%quartets%ax) + size(eri%quartets%ay) + size(eri%quartets%az) + size(eri%quartets%abx) &
                            + size(eri%quartets%aby) + size(eri%quartets%abz)) + 4_int64*size(eri%quartets%n)
      b_q_u2 = 48_int64*nqp_rep + 52_int64*nrep
      b_bins_full = 20_int64*np
      b_bins_u2 = 20_int64*nrep
      print '(a)', 'pair data, MB (production sizes; the n_shell^2 tables are listed apart):'
      print '(a,f10.2)', '  pairs object (hgp prims p,r,ra,rb,c,a; off,n)  merged 2a          ', real(b_pairs_full, dp)/MB
      print '(a,f10.2)', '  pairs object, unique list (U1/U2)                                 ', real(b_pairs_u1, dp)/MB
      print '(a,f10.2)', '  bin-interleaved tables (z,c,a,P; A,AB,n)      merged 2a, U1       ', real(b_q_full, dp)/MB
      print '(a,f10.2)', '  bin-interleaved tables, representatives only   U2                 ', real(b_q_u2, dp)/MB
      print '(a,f10.2)', '  bins (sp_i,sp_j,sp_q)                          merged 2a, U1       ', real(b_bins_full, dp)/MB
      print '(a,f10.2)', '  bins, representatives only                     U2                 ', real(b_bins_u2, dp)/MB
      print '(a,f10.2)', '  dsh (nshell^2 doubles), unchanged by all schemes                    ', 8.0_dp*nsh*nsh/MB
      print '(a,f10.2)', '  Q triangle (nshell^2/2 doubles), host transient                     ', 4.0_dp*nsh*nsh/MB
      print '(a,f10.2)', '  sym_pairs_build transient: img(npair,nop)+pidx(nshell^2/2) int32    ', &
         (4.0_dp*np*nop + 2.0_dp*nsh*nsh)/MB
      print '(a,f10.2,a,f10.2,a,f10.2)', '  TOTAL pair data   merged 2a ', real(b_pairs_full + b_q_full + b_bins_full, dp)/MB, &
         '   U1 ', real(b_pairs_u1 + b_q_full + b_bins_full, dp)/MB, '   U2 ', real(b_pairs_u1 + b_q_u2 + b_bins_u2, dp)/MB
   end subroutine counts

   subroutine source_of(p, key_src, gsel)
      !! The pairs-object key of the representative's data that pair p is the
      !! image of, and the operation (a mask index into the group) that maps it.
      integer, intent(in) :: p
      integer, intent(out) :: key_src, gsel
      integer :: rp, ri, rj, x, y, u, v, gg
      rp = tab%rep(p)
      ri = eri%bins%sp_i(rp)
      rj = eri%bins%sp_j(rp)
      x = eri%bins%sp_i(p)
      y = eri%bins%sp_j(p)
      key_src = 0
      gsel = 0
      do gg = 1, nop
         u = eri%orb%shmap(x, gg)
         v = eri%orb%shmap(y, gg)
         ! operations are involutions: g sends x to u, so it sends u back to x
         if (u == ri .and. v == rj) then
            key_src = pkey(ri, rj)
            gsel = gg
            return
         else if (u == rj .and. v == ri) then
            key_src = pkey(rj, ri)
            gsel = gg
            return
         end if
      end do
      stop 'source_of: no operation maps the pair onto its representative'
   end subroutine source_of

   subroutine reconstruct_and_verify()
      type(trc_shell_pairs_t) :: pu
      integer :: kk, key_dst, key_src, gsel, mask, d, nn, o_src, o_dst, ndiff_n, nmis, ncomp
      real(dp) :: dmax_r, dmax_a, dmax_c, v, dpos
      real(dp), allocatable :: dc(:, :), gf(:, :), gs(:, :), gk(:, :), gkf(:, :), w(:, :)
      real(dp) :: gmax, e_jk, e_k, rawerr
      type(trc_sym_ao_t) :: aos
      type(trc_sym_scf_t) :: sc
      real(dp), allocatable :: ident(:, :)
      integer :: n

      call copy_pairs(eri%pairs, pu)
      dmax_r = 0.0_dp
      dmax_a = 0.0_dp
      dmax_c = 0.0_dp
      ndiff_n = 0
      nmis = 0
      ncomp = 0
      do p = 1, np
         if (tab%rep(p) == p) cycle
         key_dst = pkey(eri%bins%sp_i(p), eri%bins%sp_j(p))
         call source_of(p, key_src, gsel)
         mask = geom%group%ops(gsel)
         nn = eri%pairs%n(key_src)
         if (eri%pairs%n(key_dst) /= nn) then
            ndiff_n = ndiff_n + 1
            cycle
         end if
         o_src = eri%pairs%off(key_src)
         o_dst = eri%pairs%off(key_dst)
         do kk = 1, nn
            ncomp = ncomp + 1
            ! unchanged: zeta, combined coefficient, i-side exponent
            if (eri%pairs%p(o_dst + kk) /= eri%pairs%p(o_src + kk) .or. &
                eri%pairs%c(o_dst + kk) /= eri%pairs%c(o_src + kk) .or. &
                eri%pairs%a(o_dst + kk) /= eri%pairs%a(o_src + kk)) then
               dmax_c = max(dmax_c, abs(eri%pairs%p(o_dst + kk) - eri%pairs%p(o_src + kk)), &
                            abs(eri%pairs%c(o_dst + kk) - eri%pairs%c(o_src + kk)) &
                            /max(abs(eri%pairs%c(o_src + kk)), 1.0e-300_dp), &
                            abs(eri%pairs%a(o_dst + kk) - eri%pairs%a(o_src + kk)))
               nmis = nmis + 1
            end if
            do d = 1, 3
               v = flip_sign(mask, d)*eri%pairs%r(o_src + kk, d)
               dmax_r = max(dmax_r, abs(v - eri%pairs%r(o_dst + kk, d)))
               v = flip_sign(mask, d)*eri%pairs%ra(o_src + kk, d)
               dmax_a = max(dmax_a, abs(v - eri%pairs%ra(o_dst + kk, d)))
               v = flip_sign(mask, d)*eri%pairs%rb(o_src + kk, d)
               dmax_a = max(dmax_a, abs(v - eri%pairs%rb(o_dst + kk, d)))
               ! the reconstruction itself
               pu%r(o_dst + kk, d) = flip_sign(mask, d)*eri%pairs%r(o_src + kk, d)
               pu%ra(o_dst + kk, d) = flip_sign(mask, d)*eri%pairs%ra(o_src + kk, d)
               pu%rb(o_dst + kk, d) = flip_sign(mask, d)*eri%pairs%rb(o_src + kk, d)
            end do
            pu%p(o_dst + kk) = eri%pairs%p(o_src + kk)
            pu%c(o_dst + kk) = eri%pairs%c(o_src + kk)
            pu%a(o_dst + kk) = eri%pairs%a(o_src + kk)
         end do
      end do
      print '(a,i0,a,i0,a,i0)', 'RECON: primitive pairs rebuilt ', ncomp, ';  pairs whose primitive count differs ', ndiff_n, &
         ';  primitives with zeta/c/a not bit-identical ', nmis
      print '(a,es10.3,a,es10.3,a,es10.3)', 'RECON: max|P - flip(P_rep)| ', dmax_r, '   max|A,B - flip| ', dmax_a, &
         '   max zeta/c/a deviation ', dmax_c

      ! swap the engine's pair data for the reconstruction and rebuild its tables
      n = bas%nao
      allocate (dc(n, n), gf(n, n), gs(n, n), gk(n, n), gkf(n, n), w(n, n), ident(n, n))
      aos = sym_ao_build(bas, geom, err)
      sc%ao = aos
      ident = 0.0_dp
      do i = 1, n
         ident(i, i) = 1.0_dp
      end do
      call sc%prepare(ident, 1, 0.0_dp, err)
      call random_inv(n, aos, sc, dc)

      call eri%fock(bas, dc, gf)
      gmax = maxval(abs(gf))
      call eri%fock(bas, dc, gk, k_scale=1.0_dp, j_scale=0.0_dp)
      gkf = gk

      eri%pairs = pu
      call bins_build_pair_quartets(eri%bins, eri%pairs, eri%quartets)
      call eri%fock(bas, dc, gs, skeleton=.true.)
      rawerr = maxval(abs(gs - gf))
      call sc%symmetrize(n, gs, w)
      e_jk = maxval(abs(gs - gf))
      call eri%fock(bas, dc, gk, k_scale=1.0_dp, j_scale=0.0_dp, skeleton=.true.)
      call sc%symmetrize(n, gk, w)
      e_k = maxval(abs(gk - gkf))
      print '(a,es10.3,a,es10.3,a,es10.3,a,es10.3)', 'EXACT: max|G| ', gmax, '  raw skeleton err ', rawerr, &
         '  symmetrized err (J+K) ', e_jk, '  (K only) ', e_k
      print '(a,es10.3,a,l1)', 'EXACT: gate 5e-11*max(1,|G|) = ', 5.0e-11_dp*max(1.0_dp, gmax), &
         '   passed ', (e_jk <= 5.0e-11_dp*max(1.0_dp, gmax)) .and. (e_k <= 5.0e-11_dp*max(1.0_dp, gmax))
      dpos = 0.0_dp
   end subroutine reconstruct_and_verify

   subroutine random_inv(n, aos, sc, d)
      integer, intent(in) :: n
      type(trc_sym_ao_t), intent(in) :: aos
      type(trc_sym_scf_t), intent(inout) :: sc
      real(dp), intent(out) :: d(n, n)
      real(dp), allocatable :: w(:, :)
      integer(int64) :: st
      integer :: a, b
      st = 12345_int64
      do b = 1, n
         do a = 1, b
            st = mod(st*6364136223846793005_int64 + 1442695040888963407_int64, 2_int64**62)
            d(a, b) = 0.4_dp*(real(mod(st/1024_int64, 20001_int64), dp)/10000.0_dp - 1.0_dp)
            d(b, a) = d(a, b)
         end do
      end do
      allocate (w(n, n))
      call sc%symmetrize(n, d, w)
      if (aos%nao < 0) stop
   end subroutine random_inv

end program spike_unique
