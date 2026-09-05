!
! Water cluster in 6-31G -- a real basis on a real geometry, which is where the
! synthetic lattice stops being informative.
!
! 6-31G water is 9 shells and 13 basis functions per molecule (O: three s shells
! and two p shells, 9 bf; H: two s shells, 2 bf).  So 75 waters is 675 shells
! and 975 basis functions, and 150 waters is 1350/1950 -- the size Jorge
! remembered, and the one the LibERI paper's (H2O)128/6-31G row sits next to.
!
! Note 6-31G has NO d functions, so this exercises the s,p regime exclusively:
! the deepest end of the memory-bound region the Asadchev-Valeev intensity
! argument describes, and the case HGP should be best at.
!
! Geometry: oxygens on a cubic lattice at the density of liquid water (~3.1 A
! nearest-neighbour), each molecule given a different deterministic orientation
! so the cluster is not pathologically symmetric -- symmetry would make the
! screening unrepresentatively effective.
!
program bench_binwater
   use trc_boys, only: dp, boys_init
   use trc_tables, only: LTOT, NHERM_MAX, tables_init
   use trc_cart, only: cart_init
   use trc_batch, only: build_pairs, eri_batch, ncart, common_fac_sp
   use trc_hgp, only: build_pairs_hgp, hgp_batch

   use trc_screen, only: schwarz_bounds
   use trc_bins, only: pair_bins_t, build_binned_pairs
   use trc_binkernel, only: fock_bins, trc_plan_t
   implicit none

   integer, parameter :: NREP = 3
   real(dp), parameter :: ANG = 1.8897261254578281_dp   !! angstrom -> bohr
   real(dp), parameter :: PI = 3.14159265358979323846_dp

   integer  :: NW = 75, ENGINE = 1
   character(len=256) :: XYZ = ''   !! optional geometry file (Angstrom)
   integer :: nat_file
   real(dp), allocatable :: at_r(:,:)
   integer,  allocatable :: at_z(:)
   real(dp) :: THRESH = 1.0e-10_dp
   character(len=256) :: DMFILE = ''
   character(len=256) :: GFILE = ''
   !! 6-31G(d): one uncontracted d per heavy atom. Selected by the
   !! environment rather than a positional argument so existing command lines
   !! keep working. Needs an LMAX=2 build.
   logical :: POLAR = .false.

   integer :: nat, nbas, nao, npp, nhpp, nq, i, j, k, l, iq, maxnp
   integer, allocatable :: sh_l(:), sh_np(:), ao_off(:), sh_at(:)
   real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), cfac(:)
   integer, allocatable :: pp_off(:), pp_n(:), hp_off(:), hp_n(:)
   real(dp), allocatable :: pp_p(:), pp_r(:, :), pp_c(:), pp_e(:, :)
   real(dp), allocatable :: hp_p(:), hp_r(:, :), hp_ra(:, :), hp_rb(:, :), hp_c(:)
   integer, allocatable :: q_i(:), q_j(:), q_k(:), q_l(:), q_off(:)
   real(dp), allocatable :: out(:), rscr(:, :), qs(:), dmax(:, :)
   integer(kind=8) :: nrect
   integer :: nout, r, nchunk, ic, c, l0, h0, nlaunch
   integer(kind=8) :: nwork, nkept
   type(pair_bins_t) :: bins
   real(dp), allocatable :: dmat(:,:), jmat(:,:), kmat(:,:), dsh(:,:)
   logical :: DSCREEN = .true.
   integer, parameter :: CHUNK = 400000
   real(dp) :: t0, t1, t, best, worst

   call read_args()
   call build_water_631g()

   call tables_init(); call cart_init(); call boys_init()

   call build_pairs(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, THRESH, &
                    pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, npp)
   call build_pairs_hgp(nbas, sh_l, sh_np, sh_e, sh_c, sh_r, cfac, &
                        hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, nhpp)

   allocate (qs(nbas*(nbas + 1)/2), dmax(nbas, nbas))
   dmax = 1.0_dp
   call tick(t0)
   call schwarz_bounds(nbas, npp, sh_l, pp_off, pp_n, pp_p, pp_r, pp_c, pp_e, qs)
   call tick(t1)
   print '(a,f8.3,a)', '  schwarz bounds  : ', t1 - t0, ' s'

   call tick(t0)
   call build_binned_pairs(nbas, sh_l, sh_np, sh_r, qs, THRESH, bins)
   call tick(t1)
   print '(a,f8.3,a)', '  bin the pairs   : ', t1 - t0, ' s'

   allocate (dmat(nao, nao), jmat(nao, nao), kmat(nao, nao))
   if (len_trim(DMFILE) > 0) then
      call read_density()
   else
      ! Model density: flat, |D| in [0.3, 0.5] everywhere, no spatial decay.
      ! Fine for validating arithmetic, useless for density screening -- see
      ! scripts/dump_density.py.
      do j = 1, nao
         do i = 1, nao
            dmat(i, j) = 0.4_dp + 0.1_dp*cos(real(3*i + 5*j, dp))
         end do
      end do
      dmat = 0.5_dp*(dmat + transpose(dmat))
   end if
   jmat = 0.0_dp; kmat = 0.0_dp

   ! Shell-block density bound for the per-quartet density screen.  With the
   ! flat model density above every entry is ~0.5 and this screens nothing,
   ! which is the honest answer -- density screening only pays on a real,
   ! spatially decaying density.  DSCREEN=0 disables it outright.
   allocate (dsh(nbas, nbas))
   if (DSCREEN) then
      do j = 1, nbas
         do i = 1, nbas
            dsh(i, j) = maxval(abs(dmat(ao_off(i):ao_off(i) + ncart(sh_l(i)) - 1, &
                                        ao_off(j):ao_off(j) + ncart(sh_l(j)) - 1)))
         end do
      end do
   else
      dsh = huge(1.0_dp)*1.0e-30_dp
   end if

   print '(a)', ''
   print '(a,i0,a,i0,a,i0,a,i0)', '  waters ', NW, '   atoms ', nat, &
      '   shells ', nbas, '   nao ', nao
   print '(a,i0,a,i0,a,i0)', '  primitive pairs ', nhpp, '   shell pairs ', bins%npair, &
      '   live bins ', bins%nlive
   block
      integer :: kk, mx
      integer(kind=8) :: seg
      mx = 0
      do kk = 1, bins%nbin
         mx = max(mx, bins%bin_cnt(kk))
      end do
      seg = int(mx, 8)*int(mx + 1, 8)/2
      print '(a,i0,a,i0,a,l1)', '  largest bin ', mx, '  -> max segment ', seg, &
         '   fits int32: ', seg < 2147483647_8
   end block
   print '(a,es9.2)', '  schwarz thresh  ', THRESH
   print '(a)', ''


   !$acc enter data copyin(sh_l, ao_off, dmat, dsh, bins, bins%sp_i, bins%sp_j, bins%sp_q, &
   !$acc                   hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c) &
   !$acc            copyin(jmat, kmat)

   call fock_bins(bins, nbas, nhpp, nao, sh_l, ao_off, THRESH, .false., &
                     1.0_dp, 1.0_dp, .false., dsh, &
                  hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                  1, dmat, jmat, kmat, 0, 1, nlaunch, nwork, nkept)
   !$acc wait
   jmat = 0.0_dp; kmat = 0.0_dp
   !$acc update device(jmat, kmat)

   ! How much of a Fock build is the launch plan? It is rebuilt on every
   ! call today, so if it is a material fraction the plan object saves it on
   ! every SCF iteration; if it is noise, the object is worth having for the
   ! API's sake and not for speed. Measure rather than assume.
   block
      type(trc_plan_t) :: pln
      real(dp) :: tp0, tp1
      call tick(tp0)
      call pln%build(bins, THRESH, .false.)
      call tick(tp1)
      print '(a,f9.5,a,i0,a,i0)', '  plan build      ', tp1 - tp0, &
         ' s   segments ', pln%nseg, '   work ', pln%nwork
      call pln%release()
   end block

   best = huge(1.0_dp); worst = 0.0_dp
   do r = 0, NREP
      ! Zero per repetition: without this the accumulator holds NREP+1 copies
      ! and any checksum taken afterwards depends on the repetition count.
      jmat = 0.0_dp; kmat = 0.0_dp
      !$acc update device(jmat, kmat)
      call tick(t0)
      call fock_bins(bins, nbas, nhpp, nao, sh_l, ao_off, THRESH, .false., &
                     1.0_dp, 1.0_dp, .false., dsh, &
                     hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, &
                     1, dmat, jmat, kmat, 0, 1, nlaunch, nwork)
      !$acc wait
      call tick(t1)
      t = t1 - t0
      if (r > 0) then
         best = min(best, t); worst = max(worst, t)
      end if
   end do

   !$acc update self(jmat, kmat)
   !$acc exit data delete(sh_l, ao_off, dmat, bins%sp_i, bins%sp_j, bins, &
   !$acc                  hp_off, hp_n, hp_p, hp_r, hp_ra, hp_rb, hp_c, jmat, kmat)

   print '(a,i0,a,i0)', '  kernel launches ', nlaunch, '   work items ', nwork
   print '(a,i0,a,f6.2,a)', '  survive screen  ', nkept, '  (', &
      100.0_dp*real(nkept, dp)/real(max(nwork, 1_8), dp), '% of work items)'
   print '(a,f9.4,a,f5.1,a)', '  FOCK BUILD      ', best, ' s   [+', &
      100.0_dp*(worst/best - 1.0_dp), '%]'
   print '(a,es10.3)', '  quartets/s      ', real(nwork, dp)/best
   !
   ! Checksums of the accumulated matrix.  The correctness tests run on 21
   ! shell pairs; this is 44,846, and the only way to check the fast path at
   ! that scale is against our own slower path on the same input.
   !
   ! The folded kernel accumulates only one triangle's worth of each
   ! permutation and at twice scale, so the physical G = 2J - K is
   ! (jmat + jmat^T)/4.  Verified elementwise against gpu4pyscf's 2*vj - vk at
   ! Gly30 (1273 functions): agreement 6.4e-7 relative, the residual being
   ! screening policy rather than arithmetic.
   !
   ! G = 2J - K = raw + raw^T, for ONE Fock build.
   !
   ! This used to read 0.25*(raw + raw^T), which was right only because the
   ! timing loop below runs NREP+1 = 4 builds and never zeroed `jmat` between
   ! them: the accumulator held four copies and the 0.25 cancelled them. The
   ! timings were unaffected -- each repetition is a full build -- but the
   ! printed checksum silently depended on NREP, and comparing it against
   ! another code only worked for NREP = 3.
   !
   ! The loop now zeroes `jmat` per repetition, so this factor is the physical
   ! one and the checksum means the same thing whatever NREP is.
   !
   jmat = jmat + transpose(jmat)
   print '(a,es24.16)', '  sum(G)          ', sum(jmat)
   print '(a,es24.16)', '  sum|G|          ', sum(abs(jmat))

   ! Dump G so it can be compared elementwise against another code, which is a
   ! far stronger check at this size than any summary norm.
   if (len_trim(GFILE) > 0) then
      block
         integer :: ug
         open (newunit=ug, file=trim(GFILE), status='replace', &
               access='stream', form='unformatted')
         write (ug) int(nao, kind=8)
         write (ug) jmat
         close (ug)
         print '(a,a)', '  wrote G to      ', trim(GFILE)
      end block
   end if
   print '(a,es24.16)', '  max|G|          ', maxval(abs(jmat))

contains

   subroutine read_args()
      character(len=32) :: a
      if (command_argument_count() >= 1) then
         call get_command_argument(1, a); read (a, *) NW
      end if
      if (command_argument_count() >= 2) then
         call get_command_argument(2, a); read (a, *) THRESH
      end if
      if (command_argument_count() >= 3) then
         call get_command_argument(3, a); read (a, *) ENGINE
      end if
      if (command_argument_count() >= 4) call get_command_argument(4, XYZ)
      if (command_argument_count() >= 5) call get_command_argument(5, DMFILE)
      if (command_argument_count() >= 6) call get_command_argument(6, GFILE)
      block
         character(len=32) :: bs
         integer :: ls
         call get_environment_variable('TRC_BASIS', bs, ls)
         POLAR = (ls > 0 .and. (index(bs, '*') > 0 .or. index(bs, 'd') > 0))
         if (POLAR) print '(a)', '  basis           : 6-31G(d)'
      end block
   end subroutine read_args

   !
   ! 6-31G, standard exponents and coefficients.  Contraction coefficients are
   ! scaled by the libcint primitive normalisation here, matching what env
   ! carries; build_pairs then applies common_fac_sp.
   !
   subroutine build_water_631g()
      real(dp) :: oS6e(6), oS6c(6), oS3e(3), oS3c(3), oP3c(3), oS1e(1)
      real(dp) :: hS3e(3), hS3c(3), hS1e(1)
      integer  :: iw, nx, ny, nz, ish, ia
      real(dp) :: o(3), h1(3), h2(3), th, ph, rr, aa, ca, sa

      oS6e = [5484.6717000_dp, 825.2349500_dp, 188.0469600_dp, &
              52.9645000_dp, 16.8975700_dp, 5.7996353_dp]
      oS6c = [0.0018311_dp, 0.0139501_dp, 0.0684451_dp, &
              0.2327143_dp, 0.4701930_dp, 0.3585209_dp]
      oS3e = [15.5396160_dp, 3.5999336_dp, 1.0137618_dp]
      oS3c = [-0.1107775_dp, -0.1480263_dp, 1.1307670_dp]
      oP3c = [0.0708743_dp, 0.3397528_dp, 0.7271586_dp]
      oS1e = [0.2700058_dp]
      hS3e = [18.7311370_dp, 2.8253937_dp, 0.6401217_dp]
      hS3c = [0.0334946_dp, 0.2347269_dp, 0.8137573_dp]
      hS1e = [0.1612778_dp]

      if (len_trim(XYZ) > 0) then
         call read_xyz()
         nat = nat_file
         nbas = 0
         do iw = 1, nat
            if (at_z(iw) > 10) then
               nbas = nbas + 7      ! second row: four s shells, three p
               if (POLAR) nbas = nbas + 1
            else if (at_z(iw) > 2) then
               nbas = nbas + 5      ! first row: three s shells, two p
               if (POLAR) nbas = nbas + 1
            else
               nbas = nbas + 2      ! H: two s shells
            end if
         end do
         NW = count(at_z == 8)
      else
         nat = 3*NW
         nbas = 9*NW
      end if
      allocate (sh_l(nbas), sh_np(nbas), sh_r(3, nbas), sh_at(nbas), cfac(nbas), ao_off(nbas))
      maxnp = 6
      allocate (sh_e(maxnp, nbas), sh_c(maxnp, nbas))
      sh_e = 0.0_dp; sh_c = 0.0_dp

      ish = 0; nao = 0
      if (len_trim(XYZ) > 0) then
         ! Real geometry: place shells by element rather than assuming O,H,H.
         do iw = 1, nat
            o = at_r(:, iw)
            call put_atom(ish, at_z(iw), o)
         end do
         return
      end if
      do iw = 0, NW - 1
         ! oxygens on a cubic lattice at liquid-water density
         nx = mod(iw, 5); ny = mod(iw/5, 5); nz = iw/25
         o = [3.1_dp*real(nx, dp), 3.1_dp*real(ny, dp), 3.1_dp*real(nz, dp)]*ANG

         ! a different deterministic orientation per molecule
         th = 0.7_dp*real(iw, dp); ph = 0.37_dp*real(iw, dp)
         rr = 0.957_dp*ANG
         aa = 104.5_dp*PI/180.0_dp
         ca = cos(th); sa = sin(th)
         h1 = o + rr*[ca*cos(ph), sa*cos(ph), sin(ph)]
         h2 = o + rr*[cos(th + aa)*cos(ph), sin(th + aa)*cos(ph), sin(ph + 0.3_dp)]

         call put(ish, 0, 6, oS6e, oS6c, o)
         call put(ish, 0, 3, oS3e, oS3c, o)
         call put(ish, 1, 3, oS3e, oP3c, o)
         call put(ish, 0, 1, oS1e, [1.0_dp], o)
         call put(ish, 1, 1, oS1e, [1.0_dp], o)
         call put(ish, 0, 3, hS3e, hS3c, h1)
         call put(ish, 0, 1, hS1e, [1.0_dp], h1)
         call put(ish, 0, 3, hS3e, hS3c, h2)
         call put(ish, 0, 1, hS1e, [1.0_dp], h2)
      end do
   end subroutine build_water_631g

   !
   ! Read a converged density written by scripts/dump_density.py: one int64
   ! holding nao, then nao*nao float64 in terco's own AO order.
   !
   subroutine read_density()
      integer :: u, ios, ii, jj
      integer(kind=8) :: nao_file
      real(dp) :: amax, amean
      open (newunit=u, file=trim(DMFILE), status='old', action='read', &
            access='stream', form='unformatted', iostat=ios)
      if (ios /= 0) then
         print '(a,a)', '  cannot open density file: ', trim(DMFILE)
         stop 1
      end if
      read (u) nao_file
      if (int(nao_file) /= nao) then
         print '(a,i0,a,i0)', '  density nao mismatch: file has ', nao_file, &
            ', basis has ', nao
         stop 1
      end if
      read (u) dmat
      close (u)
      amax = 0.0_dp; amean = 0.0_dp
      do jj = 1, nao
         do ii = 1, nao
            amax = max(amax, abs(dmat(ii, jj)))
            amean = amean + abs(dmat(ii, jj))
         end do
      end do
      amean = amean/real(nao, dp)**2
      print '(a,a)', '  density         : ', trim(DMFILE)
      print '(a,es10.3,a,es10.3)', '    max|D| ', amax, '   mean|D| ', amean
   end subroutine read_density

   !> Minimal xyz reader: count, comment, then `sym x y z` in Angstrom.
   subroutine read_xyz()
      integer :: u, i, ios
      character(len=8) :: sym
      character(len=256) :: line
      real(dp) :: x, y, z
      open (newunit=u, file=trim(XYZ), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print '(a,a)', '  cannot open geometry file: ', trim(XYZ)
         stop 1
      end if
      read (u, *) nat_file
      read (u, '(a)') line
      allocate (at_r(3, nat_file), at_z(nat_file))
      do i = 1, nat_file
         read (u, *) sym, x, y, z
         at_r(:, i) = [x, y, z]*ANG
         select case (trim(adjustl(sym)))
         case ('H', 'h'); at_z(i) = 1
         case ('C', 'c'); at_z(i) = 6
         case ('N', 'n'); at_z(i) = 7
         case ('O', 'o'); at_z(i) = 8
         case ('P', 'p'); at_z(i) = 15
         case default
            print '(a,a)', '  unsupported element in geometry: ', trim(sym)
            stop 1
         end select
      end do
      close (u)
      print '(a,a,a,i0,a)', '  geometry        : ', trim(XYZ), '  (', nat_file, ' atoms)'
   end subroutine read_xyz

   !
   ! 6-31G for H, C, N, O.  Standard Basis Set Exchange values; the SP shells
   ! are split into separate s and p sharing exponents, which is what libcint's
   ! layout wants and what gen_reference already assumes.
   !
   subroutine put_atom(ish, z, r)
      integer,  intent(inout) :: ish
      integer,  intent(in)    :: z
      real(dp), intent(in)    :: r(3)
      real(dp) :: e6(6), c6(6), e3(3), cs3(3), cp3(3), e1(1)

      select case (z)
      case (1)
         call put(ish, 0, 3, [18.7311370_dp, 2.8253937_dp, 0.6401217_dp], &
                  [0.0334946_dp, 0.2347269_dp, 0.8137573_dp], r)
         call put(ish, 0, 1, [0.1612778_dp], [1.0_dp], r)
         return
      case (6)
         e6 = [3047.5249000_dp, 457.36951000_dp, 103.94869000_dp, &
               29.210155000_dp, 9.2866630000_dp, 3.1639270000_dp]
         c6 = [0.0018347_dp, 0.0140373_dp, 0.0688426_dp, &
               0.2321844_dp, 0.4679413_dp, 0.3623120_dp]
         e3 = [7.8682724_dp, 1.8812885_dp, 0.5442493_dp]
         cs3 = [-0.1193324_dp, -0.1608542_dp, 1.1434564_dp]
         cp3 = [0.0689991_dp, 0.3164240_dp, 0.7443083_dp]
         e1 = [0.1687144_dp]
      case (7)
         e6 = [4173.5110000_dp, 627.45790000_dp, 142.90210000_dp, &
               40.234330000_dp, 12.820210000_dp, 4.3904370000_dp]
         c6 = [0.0018348_dp, 0.0139950_dp, 0.0685870_dp, &
               0.2322410_dp, 0.4690700_dp, 0.3604550_dp]
         e3 = [11.626358_dp, 2.7162800_dp, 0.7722180_dp]
         cs3 = [-0.1149610_dp, -0.1691180_dp, 1.1458520_dp]
         cp3 = [0.0675797_dp, 0.3239073_dp, 0.7408951_dp]
         e1 = [0.2120313_dp]
      case (8)
         e6 = [5484.6717000_dp, 825.23495000_dp, 188.04696000_dp, &
               52.964500000_dp, 16.897570000_dp, 5.7996353000_dp]
         c6 = [0.0018311_dp, 0.0139501_dp, 0.0684451_dp, &
               0.2327143_dp, 0.4701930_dp, 0.3585209_dp]
         e3 = [15.5396160_dp, 3.5999336_dp, 1.0137618_dp]
         cs3 = [-0.1107775_dp, -0.1480263_dp, 1.1307670_dp]
         cp3 = [0.0708743_dp, 0.3397528_dp, 0.7271586_dp]
         e1 = [0.2700058_dp]
      case (15)
         ! P.  Second row: 1s(6), 2sp(6), 3sp(3), 3sp(1) -> seven shells, 13 bf.
         call put(ish, 0, 6, [19413.3000000_dp, 2909.4200000_dp, 661.3640000_dp, &
                              185.75900000_dp, 59.194300000_dp, 20.031000000_dp], &
                  [0.0018516_dp, 0.0142062_dp, 0.0699995_dp, &
                   0.2400789_dp, 0.4847617_dp, 0.3351998_dp], r)
         e6 = [339.47800000_dp, 81.010100000_dp, 25.878000000_dp, &
               9.4522100000_dp, 3.6656600000_dp, 1.4674600000_dp]
         call put(ish, 0, 6, e6, [-0.0027822_dp, -0.0360499_dp, -0.1166310_dp, &
                                  0.0968328_dp, 0.6144180_dp, 0.4037980_dp], r)
         call put(ish, 1, 6, e6, [0.0045646_dp, 0.0336936_dp, 0.1397549_dp, &
                                  0.3393617_dp, 0.4509206_dp, 0.2385858_dp], r)
         e3 = [2.1562300_dp, 0.7489970_dp, 0.2831450_dp]
         call put(ish, 0, 3, e3, [-0.2529241_dp, 0.0328518_dp, 1.0812548_dp], r)
         call put(ish, 1, 3, e3, [-0.0177653_dp, 0.2740582_dp, 0.7854216_dp], r)
         call put(ish, 0, 1, [0.0998317_dp], [1.0_dp], r)
         call put(ish, 1, 1, [0.0998317_dp], [1.0_dp], r)
         ! 6-31G(d) polarisation: P's d exponent is 0.55, not the 0.8 the
         ! first row uses.
         if (POLAR) call put(ish, 2, 1, [0.55_dp], [1.0_dp], r)
         return
      case default
         print '(a,i0)', '  no 6-31G data for Z = ', z
         stop 1
      end select

      call put(ish, 0, 6, e6, c6, r)
      call put(ish, 0, 3, e3, cs3, r)
      call put(ish, 1, 3, e3, cp3, r)
      call put(ish, 0, 1, e1, [1.0_dp], r)
      call put(ish, 1, 1, e1, [1.0_dp], r)
      ! 6-31G(d): a single uncontracted d, exponent 0.8 across the first row.
      ! Cartesian, so six components -- which is what the published N counts
      ! (Gly30: 121*15 + 92*2 = 1999).
      if (POLAR) call put(ish, 2, 1, [0.8_dp], [1.0_dp], r)
   end subroutine put_atom

   subroutine put(ish, l, np, e, c, r)
      integer,  intent(inout) :: ish
      integer,  intent(in)    :: l, np
      real(dp), intent(in)    :: e(:), c(:), r(3)
      integer :: k
      ish = ish + 1
      sh_l(ish) = l; sh_np(ish) = np; sh_r(:, ish) = r
      cfac(ish) = common_fac_sp(l)
      do k = 1, np
         sh_e(k, ish) = e(k)
         sh_c(k, ish) = c(k)*gto_norm(l, e(k))
      end do
      ao_off(ish) = nao + 1
      nao = nao + ncart(l)
   end subroutine put

   !
   ! libcint's CINTgto_norm: 1/sqrt(gaussian_int(2l+2, 2a)),
   ! gaussian_int(n,a) = 0.5 a^(-(n+1)/2) Gamma((n+1)/2)
   !
   pure real(dp) function gto_norm(l, a)
      integer,  intent(in) :: l
      real(dp), intent(in) :: a
      real(dp) :: nn, gi
      nn = real(2*l + 2, dp)
      gi = 0.5_dp*(2.0_dp*a)**(-(nn + 1.0_dp)/2.0_dp)*gamma((nn + 1.0_dp)/2.0_dp)
      gto_norm = 1.0_dp/sqrt(gi)
   end function gto_norm

   subroutine tick(tt)
      real(dp), intent(out) :: tt
      integer(kind=8) :: c, rate
      call system_clock(c, rate)
      tt = real(c, dp)/real(rate, dp)
   end subroutine tick

end program bench_binwater
