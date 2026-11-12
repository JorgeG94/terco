!
! Shell-pair orbit tables of the skeleton Fock build (trc_sym_pairs), checked
! without a single integral. For water (C2v), B2H6 (D2h, nonplanar),
! naphthalene (D2h, planar) and benzene (D6h detected as D2h, via S12), with
! cc-pVDZ (decontracted by the four-centre build) and, for the two small
! molecules, 6-31G* (not decontracted), Cartesian and spherical caller AO
! spaces, on the bins of the production ERI object:
!
!   - the shell map is a permutation, each operation an involution, the
!     identity first; images agree with the shell map;
!   - every orbit is closed under the group and complete (the pairs that share
!     a representative are exactly the images of one pair), the
!     representatives partition the pair set, orbit sizes divide |G|;
!   - Q is invariant under every operation and never below the raw bound;
!     dsh is invariant, idempotent and never below the raw block maximum;
!   - bins closed on the representative (a pair joins the bin of its
!     representative) keep the size bucket and the angular momenta, and the
!     weighted rep-quartet count equals the full count, for every admitted
!     bin pair and in total (the S13 prototype's integral-free identity,
!     with the 1/2 of ab = cd and of the diagonal bin pair);
!   - the full admitted count reproduces the S13 prototype / production
!     `eri%nwork` (14627, 360415, 20671728, 4741001 for cc-pVDZ at 1e-10), and
!     the orbit counts (pairs, representatives, orbits per size) its output;
!   - a pair list that is not closed is refused, a Q symmetrized before the
!     selection is not, and malformed input is refused.
!
! Registered twice, as check_sym_ci is: `check_sym_pairs` (mode "main") and
! `check_sym_pairs_f` (mode "f", the uncontracted l = 0..3 set, which prints
! SKIP when this build's LMAX is below 3).
!
! NO INTERNAL PROCEDURES: the helpers are module procedures.
!
module check_sym_pairs_suite
   use pic_types, only: dp, int64
   use pic_strings, only: to_string
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t
   use trc_basis_json, only: trc_basis_from_json
   use trc_eri, only: trc_eri_t
   use trc_error, only: error_t
   use trc_sym_ao, only: trc_sym_ao_t, sym_ao_build
   use trc_sym_geom, only: trc_sym_geom_t, sym_geom_detect
   use trc_sym_pairs, only: trc_sym_pairs_t, sym_pairs_build, sym_pairs_shell_map, &
                            sym_pairs_symmetrize_q, sym_pairs_symmetrize_dsh, sym_pairs_tri
   use trc_test_basis, only: read_xyz, build_aux
   implicit none
   private
   public :: collect_sym_pairs, g_mode

   character(len=8) :: g_mode = "main"
      !! "main" or "f", set from argv before run_testsuite.

   real(dp), parameter :: THRESH = 1.0e-10_dp  !! the production threshold the prototype counted at
   character(len=32), parameter :: CCPVDZ = "basis_sets/cc-pvdz.json"
   character(len=32), parameter :: B631GS = "basis_sets/6-31gs.json"

   type :: expect_t
      !! What the S13 prototype printed for cc-pVDZ at 1e-10 (-1: not recorded).
      integer(int64) :: nfull = -1  !! admitted quartets of the full enumeration (= eri%nwork)
      integer :: npair = -1  !! significant shell pairs of the working basis
      integer :: nrep = -1  !! orbit representatives
      integer(int64) :: nraw = -1  !! quartets of the skeleton plan before the diagonal fold (reps x all, cheaper side)
   end type expect_t

contains

   subroutine collect_sym_pairs(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      if (g_mode == "f") then
         testsuite = [new_unittest("uncontracted l = 0..3: orbits, Q, dsh, sum of weights", test_lmax3)]
      else
         testsuite = [ &
                     new_unittest("water C2v cc-pVDZ, 6-31G*: orbits and weights", test_water), &
                     new_unittest("B2H6 D2h (nonplanar) cc-pVDZ, 6-31G*: orbits and weights", test_b2h6), &
                     new_unittest("benzene D2h (via S12) cc-pVDZ: orbits and weights", test_benzene), &
                     new_unittest("naphthalene D2h cc-pVDZ: orbits and weights", test_naphthalene), &
                     new_unittest("a pair set that is not closed is refused; Q symmetrized first is not", test_closure), &
                     new_unittest("malformed input is refused", test_refusals) &
                     ]
      end if
   end subroutine collect_sym_pairs

   subroutine setup_geom(name, elem, geom, expect_group, error)
      character(len=*), intent(in) :: name, expect_group
      integer, allocatable, intent(out) :: elem(:)
      type(trc_sym_geom_t), intent(out) :: geom
      type(error_type), allocatable, intent(out) :: error
      integer :: nat
      real(dp), allocatable :: r(:, :)
      type(error_t) :: err
      call read_xyz(trim(name)//".xyz", nat, elem, r)
      geom = sym_geom_detect(r, elem, err, basis_id=elem)
      call check(error, .not. err%has_error(), "sym_geom_detect: "//err%get_message())
      if (allocated(error)) return
      call check(error, geom%group%name, expect_group, "group of "//trim(name))
   end subroutine setup_geom

   subroutine run_json(name, expect_group, basis_json, ex, error)
      !! One molecule and basis file, Cartesian and spherical.
      character(len=*), intent(in) :: name, expect_group, basis_json
      type(expect_t), intent(in) :: ex
      type(error_type), allocatable, intent(out) :: error
      type(trc_sym_geom_t) :: geom
      type(trc_basis_t) :: bas
      type(error_t) :: err
      integer, allocatable :: elem(:)
      integer :: ik
      call setup_geom(name, elem, geom, expect_group, error)
      if (allocated(error)) return
      do ik = 1, 2
         call trc_basis_from_json(basis_json, size(elem), elem, geom%r_std, bas, err, force_cart=(ik == 1))
         call check(error, .not. err%has_error(), "basis: "//err%get_message())
         if (allocated(error)) return
         bas%spherical = (ik == 2)
         call analyse(bas, geom, trim(name)//" "//trim(basis_json)//merge(" sph ", " cart", ik == 2), ex, error)
         call bas%release()
         if (allocated(error)) return
      end do
   end subroutine run_json

   subroutine test_water(error)
      type(error_type), allocatable, intent(out) :: error
      call run_json("water", "C2v", CCPVDZ, expect_t(14627_int64, 171, 126, 9968_int64), error)
      if (allocated(error)) return
      call run_json("water", "C2v", B631GS, expect_t(), error)
   end subroutine test_water

   subroutine test_b2h6(error)
      type(error_type), allocatable, intent(out) :: error
      call run_json("b2h6_d2h", "D2h", CCPVDZ, expect_t(360415_int64, 866, 288, 107763_int64), error)
      if (allocated(error)) return
      call run_json("b2h6_d2h", "D2h", B631GS, expect_t(), error)
   end subroutine test_b2h6

   subroutine test_benzene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_json("tops_benzene", "D2h", CCPVDZ, expect_t(4741001_int64, 3315, 944, 1342393_int64), error)
   end subroutine test_benzene

   subroutine test_naphthalene(error)
      type(error_type), allocatable, intent(out) :: error
      call run_json("naphthalene_d2h", "D2h", CCPVDZ, expect_t(20671728_int64, 7279, 1911, 5451483_int64), error)
   end subroutine test_naphthalene

   subroutine test_lmax3(error)
      !! Uncontracted l = 0..3 shells on every atom, Cartesian and spherical.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(trc_sym_geom_t) :: geom
      integer, allocatable :: elem(:), sh_l(:), sh_np(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :)
      character(len=12), parameter :: MOL(2) = [character(len=12) :: "water", "b2h6_d2h"]
      character(len=4), parameter :: GRP(2) = [character(len=4) :: "C2v", "D2h"]
      integer :: im, ik, nsh, maxnp
      do im = 1, 2
         call setup_geom(trim(MOL(im)), elem, geom, trim(GRP(im)), error)
         if (allocated(error)) return
         do ik = 1, 2
            call build_aux(size(elem), elem, geom%r_std, 3, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
            call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, size(elem), real(elem, dp), geom%r_std, maxnp)
            bas%spherical = (ik == 2)
            call analyse(bas, geom, trim(MOL(im))//" l<=3 "//merge("sph ", "cart", ik == 2), expect_t(), error)
            call bas%release()
            if (allocated(error)) return
         end do
      end do
   end subroutine test_lmax3

   subroutine analyse(bas, geom, label, ex, error)
      !! Build the production bins of `bas`, the orbit tables of their pairs,
      !! and run every identity of the header.
      type(trc_basis_t), intent(inout), target :: bas
      type(trc_sym_geom_t), intent(in) :: geom
      character(len=*), intent(in) :: label
      type(expect_t), intent(in) :: ex
      type(error_type), allocatable, intent(out) :: error
      type(trc_eri_t), target :: eri
      type(trc_basis_t), pointer :: wb
      type(trc_basis_t) :: wbc
      type(trc_sym_ao_t) :: aow
      type(trc_sym_pairs_t) :: tab
      type(error_t) :: err
      integer :: np, ncnt(8), n

      call bas%to_device()
      call eri%build(bas, THRESH)
      if (eri%dec%active) then
         wb => eri%pb
      else
         wb => bas
      end if
      wbc = wb
      wbc%spherical = .false.  ! the four-centre build works in Cartesian functions
      aow = sym_ao_build(wbc, geom, err)
      call check(error, .not. err%has_error(), label//": sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      np = eri%bins%npair
      block
         integer, allocatable :: shmap(:, :)
         call sym_pairs_shell_map(wbc, aow, geom, shmap, err)
         call check(error, .not. err%has_error(), label//": shell map: "//err%get_message())
         if (allocated(error)) return
         tab = sym_pairs_build(shmap, eri%bins%sp_i(1:np), eri%bins%sp_j(1:np), eri%bins%sp_q(1:np), err)
         call check(error, .not. err%has_error(), label//": sym_pairs_build: "//err%get_message())
         if (allocated(error)) return
      end block

      call check_shell_map(tab, label, error)
      if (allocated(error)) return
      call check_orbits(tab, eri%bins%sp_i(1:np), eri%bins%sp_j(1:np), label, error)
      if (allocated(error)) return
      call check_screening(tab, eri%bins%sp_q(1:np), label, error)
      if (allocated(error)) return
      call check_weights(tab, eri, label, ex%nfull, ex%nraw, error)
      if (allocated(error)) return

      if (ex%npair >= 0) call check(error, np, ex%npair, label//": significant pairs")
      if (allocated(error)) return
      if (ex%nrep >= 0) call check(error, tab%nrep, ex%nrep, label//": representatives")
      if (allocated(error)) return
      ncnt = 0
      n = size(tab%size_counts())
      ncnt(1:n) = tab%size_counts()
      print "(a,a,a,i0,a,i0,a,i0,a,8i9)", "  ", label, ": nshell ", wbc%nshell, " pairs ", np, " orbits ", tab%nrep, &
         "  orbits of size 1..nop ", ncnt(1:n)
      call eri%release()
   end subroutine analyse

   subroutine check_shell_map(tab, label, error)
      type(trc_sym_pairs_t), intent(in) :: tab
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      integer :: s, g, nbad
      nbad = count(tab%shmap(:, 1) /= [(s, s=1, tab%nshell)])
      call check(error, nbad == 0, label//": operation 1 is not the identity on shells")
      if (allocated(error)) return
      nbad = 0
      do g = 1, tab%nop
         do s = 1, tab%nshell
            if (tab%shmap(tab%shmap(s, g), g) /= s) nbad = nbad + 1
         end do
      end do
      call check(error, nbad == 0, label//": "//to_string(nbad)//" shells are not sent back by applying an operation twice")
   end subroutine check_shell_map

   subroutine check_orbits(tab, sp_i, sp_j, label, error)
      !! Closure, completeness and partition.
      type(trc_sym_pairs_t), intent(in) :: tab
      integer, intent(in) :: sp_i(:), sp_j(:)
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: members(:)
      integer :: p, g, r, nbad(5), nd
      integer(int64) :: k1, k2

      nbad = 0
      allocate (members(tab%npair))
      members = 0
      do p = 1, tab%npair
         r = tab%rep(p)
         if (r < 1 .or. r > tab%npair) then
            nbad(1) = nbad(1) + 1
            cycle
         end if
         members(r) = members(r) + 1
         if (tab%rep(r) /= r) nbad(1) = nbad(1) + 1       ! the representative represents itself
         if (modulo(tab%nop, tab%osz(p)) /= 0 .or. tab%osz(p) < 1) nbad(2) = nbad(2) + 1
         nd = 0
         do g = 1, tab%nop
            ! the image is the pair the shell map gives
            k1 = sym_pairs_tri(tab%shmap(sp_i(p), g), tab%shmap(sp_j(p), g))
            k2 = sym_pairs_tri(sp_i(tab%img(p, g)), sp_j(tab%img(p, g)))
            if (k1 /= k2) nbad(3) = nbad(3) + 1
            ! closed: same representative and size as the image
            if (tab%rep(tab%img(p, g)) /= r .or. tab%osz(tab%img(p, g)) /= tab%osz(p)) nbad(4) = nbad(4) + 1
            if (.not. any(tab%img(p, 1:g - 1) == tab%img(p, g))) nd = nd + 1
         end do
         if (nd /= tab%osz(p)) nbad(2) = nbad(2) + 1
      end do
      ! complete: the pairs that share a representative are exactly its osz images
      do p = 1, tab%npair
         if (tab%rep(p) == p .and. members(p) /= tab%osz(p)) nbad(5) = nbad(5) + 1
      end do
      call check(error, all(nbad == 0), label//": orbit defects (rep/partition, size, image, closure, completeness) "// &
                 to_string(nbad(1))//" "//to_string(nbad(2))//" "//to_string(nbad(3))//" "//to_string(nbad(4))//" "// &
                 to_string(nbad(5)))
      if (allocated(error)) return
      call check(error, sum(members) == tab%npair, label//": the representatives do not partition the pairs")
      if (allocated(error)) return
      call check(error, sum(tab%osz(tab%rep_list)) == tab%npair, label//": sum of orbit sizes over representatives /= pairs")
      if (allocated(error)) return
      call check(error, count(tab%rep == [(p, p=1, tab%npair)]) == tab%nrep, label//": nrep")
      if (allocated(error)) return
      call check(error, size(tab%rep_list) == tab%nrep, label//": rep_list size")
      if (allocated(error)) return
      do p = 2, tab%nrep
         if (sym_pairs_tri(sp_i(tab%rep_list(p)), sp_j(tab%rep_list(p))) <= &
             sym_pairs_tri(sp_i(tab%rep_list(p - 1)), sp_j(tab%rep_list(p - 1)))) then
            call check(error, .false., label//": rep_list is not in ascending canonical key")
            return
         end if
      end do
      ! the representative has the smallest key of its orbit
      nbad = 0
      do p = 1, tab%npair
         k1 = sym_pairs_tri(sp_i(p), sp_j(p))
         k2 = sym_pairs_tri(sp_i(tab%rep(p)), sp_j(tab%rep(p)))
         if (k2 > k1) nbad(1) = nbad(1) + 1
      end do
      call check(error, nbad(1) == 0, label//": a representative is not the smallest key of its orbit")
   end subroutine check_orbits

   subroutine check_screening(tab, sp_q, label, error)
      !! Q and dsh are invariant under every operation.
      type(trc_sym_pairs_t), intent(in) :: tab
      real(dp), intent(in) :: sp_q(:)
      character(len=*), intent(in) :: label
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: d(:, :), d0(:, :)
      integer :: p, g, i, j, nbad
      nbad = 0
      do p = 1, tab%npair
         if (tab%q(p) < sp_q(p)) nbad = nbad + 1
         do g = 1, tab%nop
            if (tab%q(tab%img(p, g)) /= tab%q(p)) nbad = nbad + 1
         end do
      end do
      call check(error, nbad == 0, label//": "//to_string(nbad)//" Q values below the raw bound or not orbit-invariant")
      if (allocated(error)) return

      allocate (d(tab%nshell, tab%nshell), d0(tab%nshell, tab%nshell))
      do j = 1, tab%nshell
         do i = 1, tab%nshell
            d(i, j) = abs(cos(real(2*i + 3*j, dp)) + cos(real(2*j + 3*i, dp)))
         end do
      end do
      d0 = d
      call sym_pairs_symmetrize_dsh(tab%shmap, d)
      nbad = 0
      do g = 1, tab%nop
         do j = 1, tab%nshell
            do i = 1, tab%nshell
               if (d(tab%shmap(i, g), tab%shmap(j, g)) /= d(i, j)) nbad = nbad + 1
            end do
         end do
      end do
      call check(error, nbad == 0, label//": dsh is not invariant after symmetrization ("//to_string(nbad)//")")
      if (allocated(error)) return
      call check(error, all(d >= d0), label//": symmetrized dsh below the raw maximum")
      if (allocated(error)) return
      d0 = d
      call sym_pairs_symmetrize_dsh(tab%shmap, d)
      call check(error, all(d == d0), label//": dsh symmetrization is not idempotent")
   end subroutine check_screening

   subroutine check_weights(tab, eri, label, nfull_ref, nraw_ref, error)
      !! Bins closed on the representative; the weighted rep-quartet count
      !! equals the full count for every admitted bin pair and in total.
      type(trc_sym_pairs_t), intent(in) :: tab
      type(trc_eri_t), intent(in) :: eri
      character(len=*), intent(in) :: label
      integer(int64), intent(in) :: nfull_ref, nraw_ref
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: binorig(:), cb(:), live(:)
      integer(int64), allocatable :: nall(:), nrp(:), wsum(:)
      integer :: nb, k, p, a, b, ba, bb, nlive, smax_keep, nbad_type, nbad_pair
      integer(int64) :: full2, sk2a, sk2b, tot_full2, tot_sk2, raw, tot_raw, tot_nfull

      nb = eri%bins%nbin
      allocate (binorig(tab%npair), cb(tab%npair), nall(nb), nrp(nb), wsum(nb))
      do k = 1, nb
         do p = eri%bins%bin_off(k) + 1, eri%bins%bin_off(k) + eri%bins%bin_cnt(k)
            binorig(p) = k
         end do
      end do
      nall = 0
      nrp = 0
      wsum = 0
      nbad_type = 0
      do p = 1, tab%npair
         ba = binorig(p)
         bb = binorig(tab%rep(p))
         cb(p) = bb
         ! an orbit shares the size bucket, the angular momenta and the general side
         if (eri%bins%bin_s(ba) /= eri%bins%bin_s(bb) .or. eri%bins%bin_la(ba) /= eri%bins%bin_la(bb) .or. &
             eri%bins%bin_lb(ba) /= eri%bins%bin_lb(bb) .or. eri%bins%bin_g(ba) /= eri%bins%bin_g(bb)) then
           nbad_type = nbad_type + 1
         end if
         nall(bb) = nall(bb) + 1
         if (tab%rep(p) == p) then
            nrp(bb) = nrp(bb) + 1
            wsum(bb) = wsum(bb) + tab%osz(p)
         end if
      end do
      call check(error, nbad_type == 0, label//": "//to_string(nbad_type)//" pairs differ from their representative in bin type")
      if (allocated(error)) return

      nlive = count(nall > 0)
      allocate (live(nlive))
      live = pack([(k, k=1, nb)], nall > 0)
      smax_keep = int(-log10(THRESH))*eri%bins%sres
      tot_full2 = 0
      tot_sk2 = 0
      tot_raw = 0
      tot_nfull = 0
      nbad_pair = 0
      do a = 1, nlive
         ba = live(a)
         ! the orbits of a closed bin add up to the bin
         if (wsum(ba) /= nall(ba)) nbad_pair = nbad_pair + 1
         do b = 1, a
            bb = live(b)
            if (eri%bins%bin_s(ba) + eri%bins%bin_s(bb) > smax_keep) cycle
            if (a == b) then
               ! reps x all, scale 1/2; the triangle with 1/2 on ab = cd is nall^2 / 2
               full2 = nall(ba)*nall(ba)
               sk2a = wsum(ba)*nall(ba)
               raw = nrp(ba)*nall(ba)
               tot_nfull = tot_nfull + nall(ba)*(nall(ba) + 1)/2
            else
               full2 = 2*nall(ba)*nall(bb)
               sk2a = 2*wsum(ba)*nall(bb)
               sk2b = 2*wsum(bb)*nall(ba)
               if (sk2a /= sk2b) nbad_pair = nbad_pair + 1
               raw = min(nrp(ba)*nall(bb), nrp(bb)*nall(ba))
               tot_nfull = tot_nfull + nall(ba)*nall(bb)
            end if
            if (sk2a /= full2) nbad_pair = nbad_pair + 1
            tot_full2 = tot_full2 + full2
            tot_sk2 = tot_sk2 + sk2a
            tot_raw = tot_raw + raw
         end do
      end do
      call check(error, nbad_pair == 0, label//": "//to_string(nbad_pair)// &
                 " admitted bin pairs where the weighted count /= the full one")
      if (allocated(error)) return
      call check(error, tot_sk2 == tot_full2, label//": total weighted count /= full count")
      if (allocated(error)) return
      if (nfull_ref >= 0) then
         call check(error, tot_nfull == nfull_ref, label//": admitted quartets "//to_string(tot_nfull)// &
                    " against the prototype's "//to_string(nfull_ref))
         if (allocated(error)) return
      end if
      if (nraw_ref >= 0) then
         call check(error, tot_raw == nraw_ref, label//": rep quartets "//to_string(tot_raw)// &
                    " against the prototype's "//to_string(nraw_ref))
         if (allocated(error)) return
      end if
      ! the skeleton must actually save quartets (|G| > 1 and some orbit larger than 1)
      call check(error, tot_raw < tot_nfull, label//": the rep quartets are not fewer than the full ones")
      print "(a,a,a,i0,a,i0,a,f6.2)", "  ", label, ": admitted ", tot_nfull, &
         "  rep quartets (before the diagonal fold) ", tot_raw, &
         "  ratio ", real(tot_nfull, dp)/real(max(tot_raw, 1_int64), dp)
   end subroutine check_weights

   subroutine fill_q(nsh, q)
      !! A deterministic triangular Q array in (0, 1] that does not know the group.
      integer, intent(in) :: nsh
      real(dp), intent(out) :: q(:)
      integer :: i, j
      do i = 1, nsh
         do j = 1, i
            q(int(sym_pairs_tri(i, j))) = 0.5_dp + 0.5_dp*cos(real(7*i + 11*j, dp)*0.37_dp)
         end do
      end do
   end subroutine fill_q

   subroutine select_pairs(q, nsh, cut, sp_i, sp_j, sp_q)
      real(dp), intent(in) :: q(:), cut
      integer, intent(in) :: nsh
      integer, allocatable, intent(out) :: sp_i(:), sp_j(:)
      real(dp), allocatable, intent(out) :: sp_q(:)
      integer :: i, j, n
      real(dp) :: qmax
      qmax = maxval(q(1:nsh*(nsh + 1)/2))
      allocate (sp_i(nsh*(nsh + 1)/2), sp_j(nsh*(nsh + 1)/2), sp_q(nsh*(nsh + 1)/2))
      n = 0
      do i = 1, nsh
         do j = 1, i
            if (q(int(sym_pairs_tri(i, j)))*qmax > cut) then
               n = n + 1
               sp_i(n) = i
               sp_j(n) = j
               sp_q(n) = q(int(sym_pairs_tri(i, j)))
            end if
         end do
      end do
      sp_i = sp_i(1:n)
      sp_j = sp_j(1:n)
      sp_q = sp_q(1:n)
   end subroutine select_pairs

   subroutine small_water(bas, geom, wbc, aow, shmap, error)
      type(trc_basis_t), intent(out) :: bas
      type(trc_sym_geom_t), intent(out) :: geom
      type(trc_basis_t), intent(out) :: wbc
      type(trc_sym_ao_t), intent(out) :: aow
      integer, allocatable, intent(out) :: shmap(:, :)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: elem(:)
      type(error_t) :: err
      call setup_geom("water", elem, geom, "C2v", error)
      if (allocated(error)) return
      call trc_basis_from_json(B631GS, size(elem), elem, geom%r_std, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), "basis: "//err%get_message())
      if (allocated(error)) return
      bas%spherical = .false.
      wbc = bas
      aow = sym_ao_build(wbc, geom, err)
      call check(error, .not. err%has_error(), "sym_ao_build: "//err%get_message())
      if (allocated(error)) return
      call sym_pairs_shell_map(wbc, aow, geom, shmap, err)
      call check(error, .not. err%has_error(), "shell map: "//err%get_message())
   end subroutine small_water

   subroutine test_closure(error)
      !! A threshold cut on a Q that is not orbit-invariant leaves an open pair
      !! set, which is refused; the cut on the symmetrized Q is closed.
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas, wbc
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: aow
      type(trc_sym_pairs_t) :: tab
      type(error_t) :: err
      integer, allocatable :: shmap(:, :), sp_i(:), sp_j(:)
      real(dp), allocatable :: q(:), qs(:), sp_q(:)
      integer :: nsh, nq
      real(dp), parameter :: CUT = 0.6_dp

      call small_water(bas, geom, wbc, aow, shmap, error)
      if (allocated(error)) return
      nsh = wbc%nshell
      nq = nsh*(nsh + 1)/2
      allocate (q(nq), qs(nq))
      call fill_q(nsh, q)
      call select_pairs(q, nsh, CUT, sp_i, sp_j, sp_q)
      tab = sym_pairs_build(shmap, sp_i, sp_j, sp_q, err)
      call check(error, err%has_error(), "an open pair set was accepted")
      if (allocated(error)) return

      qs = q
      call sym_pairs_symmetrize_q(shmap, qs)
      call check(error, all(qs >= q), "symmetrized Q below the raw Q")
      if (allocated(error)) return
      call select_pairs(qs, nsh, CUT, sp_i, sp_j, sp_q)
      call check(error, size(sp_i) > 0 .and. size(sp_i) < nq, "the synthetic cut selects all or no pairs")
      if (allocated(error)) return
      err = error_t()
      tab = sym_pairs_build(shmap, sp_i, sp_j, sp_q, err)
      call check(error, .not. err%has_error(), "closed pair set refused: "//err%get_message())
      if (allocated(error)) return
      call check(error, all(tab%q == sp_q), "Q of a symmetrized list changed")
      if (allocated(error)) return
      q = qs
      call sym_pairs_symmetrize_q(shmap, qs)
      call check(error, all(qs == q), "Q symmetrization is not idempotent")
      call bas%release()
   end subroutine test_closure

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas, wbc
      type(trc_sym_geom_t) :: geom
      type(trc_sym_ao_t) :: aow
      type(trc_sym_pairs_t) :: tab
      type(error_t) :: err
      integer, allocatable :: shmap(:, :)
      integer :: nsh

      call small_water(bas, geom, wbc, aow, shmap, error)
      if (allocated(error)) return
      nsh = wbc%nshell
      ! a pair listed twice
      tab = sym_pairs_build(shmap, [1, 1], [1, 1], [1.0_dp, 1.0_dp], err)
      call check(error, err%has_error(), "a duplicated pair was accepted")
      if (allocated(error)) return
      ! a shell that does not exist
      err = error_t()
      tab = sym_pairs_build(shmap, [nsh + 1], [1], [1.0_dp], err)
      call check(error, err%has_error(), "a nonexistent shell was accepted")
      if (allocated(error)) return
      ! sizes that disagree
      err = error_t()
      tab = sym_pairs_build(shmap, [1, 2], [1, 1], [1.0_dp], err)
      call check(error, err%has_error(), "mismatched array sizes were accepted")
      if (allocated(error)) return
      ! an AO map built on another basis
      err = error_t()
      wbc%spherical = .true.
      call sym_pairs_shell_map(wbc, aow, geom, shmap, err)
      call check(error, err%has_error(), "an AO map of another AO space was accepted")
      call bas%release()
   end subroutine test_refusals

end module check_sym_pairs_suite

program check_sym_pairs
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use trc_tables, only: LMAX
   use check_sym_pairs_suite, only: collect_sym_pairs, g_mode
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, g_mode)
   if (trim(g_mode) == "f" .and. LMAX < 3) then
      print "(a,i0,a)", "  SKIP: the l = 0..3 set needs LMAX = 3 and this build has LMAX = ", LMAX, "."
      stop 0
   end if
   stat = 0
   call run_testsuite(collect_sym_pairs, error_unit, stat, parallel=.false.)
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_sym_pairs
