!
! The McMurchie-Davidson J engine: the Coulomb matrix without a quartet.
!
! J_ab = sum_cd (ab|cd) D_cd is formed in Hermite space, in three stages
! (Palethorpe and Barca, JCTC 21, 9388 (2025); the algebra and the kernels
! are scripts/gen_md_j.py's, checked in numpy by scripts/md/reference.py):
!
!   ket      per primitive pair Q:  Dt[Q] = (-1)^|tau| c_Q E^cd D_cd
!   Hermite  per (P, Q):            Jt[P] += [P|Q] Dt[Q]
!   bra      per shell pair:        J_ab = sum_P c_P E^ab Jt[P]
!
! Only the Hermite stage is quadratic in pairs, and it depends on the bra
! and ket only through l_a + l_b and l_c + l_d. No Cartesian quartet is
! ever formed, which is what every quartet algorithm here (HGP, Rys,
! rotated axis) ran out of registers and local memory on.
!
! SCREENING. Each (P, Q) is skipped when I_P ||rho_Q|| <= cut, with I_P the
! largest Coulomb self-norm of P's component products and ||rho_Q|| the
! Coulomb norm of Q's share of the density. Coulomb is an inner product, so
! by Cauchy-Schwarz that product bounds what the pair adds to any J
! element: the screen is rigorous and density-weighted in one number, and
! it is formed by the ket stage itself every build.
!
! DATA. A primitive pair list in (l_a + l_b, l_a) type order, so each lp
! block is contiguous; Dt and Jt share one flat array layout, an lp block
! at a time, each Hermite row of a block contiguous over the block's pairs.
! Everything lives on the device from `build` to `release`.
!
! P <-> Q SYMMETRY (the default; CMake TERCO_MDJ_SYM_LMAX). The Hermite stage evaluates each
! [P|Q] twice: once in the [lp|lq] launch for Jt[P], once in [lq|lp] for
! Jt[Q]. gen_md_j.py can instead emit symmetric kernels (trc_mdj_s<lp><lq>,
! lp >= lq, CMake TERCO_MDJ_SYM_LMAX >= 0) that evaluate each unordered pair
! once and add [P|Q] Dt[Q] to Jt[P] (in registers) and [Q|P] Dt[P] to Jt[Q]
! (an FP64 atomic per Hermite slot of Q per pair, on the smaller block);
! [Q|P]_{t+tau} = (-1)^|t+tau| [P|Q]_{t+tau}. Each direction keeps its own
! screen (I_P ||rho_Q|| for Jt[P], I_Q ||rho_P|| for Jt[Q]); R is formed if
! either passes, so J equals the gather form's to round-off (measured
! max|dJ| ~1e-12 of max|J| 58.6, gly20 cc-pVDZ). Two layout details mattered:
! each thread starts its ket chunk at its own offset (all of a warp's
! atomics on one address made the stage 2.4x slower), and a diagonal block
! pairs P with Q = P + s (mod n), s = 0..n/2, so no thread idles on a triangle.
! Tiling 2-4 bra pairs per thread to share the scatter's atomics doubled the
! registers and lost.
!
! Why it is not 2x: the screen is density-weighted per direction, so the two
! launches skip different pairs, and only pairs kept by both are shared. On
! water20 cc-pVDZ with the SCF density the union of the two kept sets is
! 0.553 of what the gather form evaluates (0.54-0.56 in every class): that is
! the ceiling. The kernels are FP64-pipe bound (ncu, [1|0]: math-pipe
! throttle dominant, atomics ~7% of L2), and the symmetric one costs more per
! evaluation: 1.6x the instructions, 20.3 against 24.4 active lanes per warp,
! 28.6% against 40% occupancy, and 200-255 registers with spills in the
! heavy classes.
!
! Measured (V100; GPU partly shared with another job, so +-10%): J, SCF
! builds 2-4 (nsys), gly20 cc-pVDZ 0.426 s gather against 0.384 s symmetric.
! In-process interleaved A/B on a synthetic density, all-symmetric over
! all-gather (median): gly20 cc-pVDZ 0.85, water20 cc-pVDZ 0.74, gly50
! cc-pVDZ 0.91 (0.67 with symmetry up to lp + lq = 5-7, noisy), gly20
! cc-pVTZ ~1.00 at every cut, water40 cc-pVTZ 1.15 (worse at every cut).
! Not a clear win where J matters (cc-pVTZ), so the default was gather.
!
! It is the default now. Once boys_eval stopped calling erf and exp for every
! distant pair, the gather form's second evaluation of each kept pair was
! the cost left to save, and the symmetric kernels win on a dedicated V100
! (RHF, 1e-11, full builds, Fock s/build, gather -> symmetric): gly50
! cc-pVDZ 4.49 -> 4.23, def2-SVP 2.98 -> 2.75, gly30 def2-TZVP 20.7 -> 20.6,
! water40 cc-pVTZ 11.8 -> 11.8, gly20 cc-pVTZ 16.5 -> 16.5; energies equal to
! 1e-12. Compile time, ninja: LMAX=2 71 s -> 15 s of MD J modules; LMAX=3
! 327 s -> 338 s, the slowest module 34 s -> 61 s.
! TERCO_MDJ_GATHER_ALL keeps both sets so TRC_MDJ_SYM_LMAX races them at run
! time. Left to try: the symmetric off-diagonal kernels reading Dt[P] from
! memory (built, fewer registers, e.g. [4|1] 192 -> 150, not timed on a
! quiet GPU); a per-class cut chosen by measurement; a screen that is
! symmetric by construction (a direction-free bound), which would make
! every pair shared and lift the 0.55 ceiling at a looser screen.
!
! REFERENCES (inspiration; the formulation and code here are terco's own)
!   L. E. McMurchie, E. R. Davidson, J. Comput. Phys. 26, 218 (1978).
!   E. Palethorpe, G. M. J. Barca, J. Chem. Theory Comput. 21, 9388 (2025).
!   I. S. Ufimtsev, T. J. Martinez, J. Chem. Theory Comput. 4, 222 (2008);
!     5, 1004 (2009).
!   A. Asadchev, E. F. Valeev, J. Phys. Chem. A 127, 10889 (2023);
!     J. Chem. Phys. 160, 244109 (2024).
!
module trc_mdj
   use pic_types, only: dp, int64
   use trc_basis_types, only: trc_basis_t
   use trc_bins, only: pair_bins_t
   use trc_mdj_kernels, only: mdj_ket, mdj_bra, mdj_bound, mdj_hermite, mdj_hermite_sym, &
                              MDJ_LPMAX, MDJ_SYM_LMAX, MDJ_GATHER_ALL
   use pic_profiler, only: profiler_start, profiler_stop
   implicit none
   private
   public :: trc_mdj_t

   integer, parameter :: MDJ_THREADS_TARGET = 262144
   !! Threads a Hermite launch aims for: past a V100's 80 SMs x 2048 resident
   !! threads, and an A100's 108 x 2048. One thread per bra primitive pair
   !! gave a heavy class 128 threads, each walking every ket pair: 27 ms
   !! for (dd|ss) on w20 cc-pVDZ, against 2 ms for all of L = 6 in HGP.
   integer, parameter :: MDJ_QLEN_MIN = 32
   !! Fewest kets per chunk: each chunk ends in one atomic per Hermite slot.

   real(dp), parameter :: MDJ_MARGIN_DEFAULT = 1.0e-3_dp
   !! The (P, Q) cut sits this far below the screening threshold: a
   !! contracted quartet is a sum over primitive quartets, and the same
   !! margin as the HGP kernels' primitive bound (TRC_PRIM_MARGIN).

   type :: trc_mdj_t
      logical :: built = .false.
      logical :: on_device = .false.
      integer :: nhp = 0
      !! primitive pairs
      integer :: nsp = 0
      !! shell pairs with at least one kept primitive pair
      integer :: ndt = 0
      !! length of the Hermite arrays dt and jt
      integer :: ntype = 0
      real(dp) :: cut = 0.0_dp
      integer :: sym_lmax = MDJ_SYM_LMAX
      !! classes with lp + lq up to this use the symmetric Hermite kernels
      real(dp), allocatable :: mp_z(:), mp_px(:), mp_py(:), mp_pz(:)
      !! exponent p and centre P of each primitive pair
      real(dp), allocatable :: mp_ax(:), mp_ay(:), mp_az(:), mp_bx(:), mp_by(:), mp_bz(:)
      !! P - A and P - B
      real(dp), allocatable :: mp_c(:)
      !! c_a c_b exp(-ab/p |AB|^2)
      real(dp), allocatable :: mp_ip(:), mp_nq(:)
      !! bra bound I_P (geometry) and ket norm ||rho_Q|| (each build)
      integer,  allocatable :: mp_sp(:)
      !! shell pair of each primitive pair
      integer,  allocatable :: sp_ao1(:), sp_ao2(:), sp_k0(:), sp_nk(:)
      !! first AO of each shell, first primitive pair less one, count
      integer,  allocatable :: sp_sh1(:), sp_sh2(:)
      !! the two shells, for the density bound of an exchange quartet
      integer,  allocatable :: sp_as1(:), sp_as2(:)
      !! first spherical AO of each shell (trc_mdg's spherical exchange)
      integer :: nao_sph = 0
      real(dp), allocatable :: sp_q(:)
      !! the bins' Schwarz bound of the pair
      integer,  allocatable :: ty_la(:), ty_lb(:), ty_s0(:), ty_ns(:), ty_k0(:), ty_nk(:)
      !! per (la, lb) type: shell pairs and primitive pairs, each first-less-one and count
      integer :: cl_k0(0:MDJ_LPMAX) = 0, cl_nk(0:MDJ_LPMAX) = 0, cl_hb(0:MDJ_LPMAX) = 0
      !! per lp block: first primitive pair less one, count, first Hermite slot less one
      real(dp), allocatable :: dt(:), jt(:)
      !
      ! The skeleton Fock build's bra (docs/SYMMETRY_FOCK.md 4.5). Built only
      ! when `build` is given the pair orbit weights. The orbit
      ! representatives' primitive pairs are COPIED after the nhp ordinary ones
      ! (so the exchange engine's view of the lists is unchanged), the copy's
      ! c already multiplied by the orbit size. The Hermite stage then runs P
      ! over the copies only and Q over all pairs, and the bra stage reads
      ! the copies' blocks of jt: the same kernels, other offsets.
      !
      logical :: has_reps = .false.
      !! the representative copies exist
      integer :: nhp_all = 0
      !! primitive pairs with the copies: the length of the mp_ arrays
      integer :: nrp = 0
      !! representative shell pairs (the copies)
      integer, allocatable :: sp_w(:)
      !! per shell pair: its orbit size when it is an orbit representative, else 0
      integer,  allocatable :: rsp_ao1(:), rsp_ao2(:), rsp_k0(:), rsp_nk(:), rsp_w(:)
      !! first AO of each shell, first primitive pair less one, count, of the copies
      integer,  allocatable :: rty_s0(:), rty_ns(:), rty_k0(:), rty_nk(:)
      !! per type: the copies' shell pairs and primitive pairs, as ty_s0.. are for the full lists
      integer :: cr_k0(0:MDJ_LPMAX) = 0, cr_nk(0:MDJ_LPMAX) = 0, cr_hb(0:MDJ_LPMAX) = 0
      !! per lp block of copies: first primitive pair less one, count, first Hermite slot less one
      integer(int64) :: npq = 0
      !! (P, Q) primitive pairs the last `coulomb` handed the Hermite stage: the
      !! work counter (a skeleton call's is a fraction of a full one's)
   contains
      procedure :: build => mdj_build
      procedure :: coulomb => mdj_coulomb
      procedure :: release => mdj_release
   end type trc_mdj_t

contains

   subroutine mdj_build(this, b, bins, thresh, pair_w)
      !! The primitive pair list over `bins`' shell pairs of `b`, its bounds,
      !! and the Hermite work arrays, all on the device.
      class(trc_mdj_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      !! the basis the bins were built over (the split one on a general contraction)
      type(pair_bins_t), intent(in) :: bins
      real(dp), intent(in) :: thresh
      integer, intent(in), optional :: pair_w(:)
      !! per bins pair: its orbit size if it is an orbit representative, else
      !! 0 (trc_bins' bins_rep_weights). Given, the skeleton bra exists and
      !! `coulomb`'s skeleton option can be used.

      integer :: t, lp

      call this%release()
      call order_pairs(this, b, bins)
      if (present(pair_w)) then
         call fill_prims(this, b, bins, thresh, pair_w)
      else
         call fill_prims(this, b, bins, thresh)
      end if

      this%ndt = 0
      do lp = 0, MDJ_LPMAX
         this%cl_hb(lp) = this%ndt
         this%ndt = this%ndt + this%cl_nk(lp)*nherm(lp)
      end do
      do lp = 0, MDJ_LPMAX
         this%cr_hb(lp) = this%ndt
         this%ndt = this%ndt + this%cr_nk(lp)*nherm(lp)
      end do
      allocate (this%dt(max(this%ndt, 1)), this%jt(max(this%ndt, 1)))
      allocate (this%mp_ip(max(this%nhp_all, 1)), this%mp_nq(max(this%nhp_all, 1)))
      this%cut = thresh*margin()
      this%sym_lmax = sym_lmax_choice()

      !$acc enter data copyin(this%mp_z, this%mp_px, this%mp_py, this%mp_pz, &
      !$acc                   this%mp_ax, this%mp_ay, this%mp_az, &
      !$acc                   this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, this%mp_sp, &
      !$acc                   this%sp_ao1, this%sp_ao2, this%sp_k0, this%sp_nk, &
      !$acc                   this%sp_sh1, this%sp_sh2, this%sp_q, this%sp_as1, this%sp_as2, &
      !$acc                   this%rsp_ao1, this%rsp_ao2, this%rsp_k0, this%rsp_nk) &
      !$acc            create(this%mp_ip, this%mp_nq, this%dt, this%jt)
      this%on_device = .true.

      do t = 1, this%ntype
         if (this%ty_nk(t) == 0) cycle
         call mdj_bound(this%ty_la(t), this%ty_lb(t), this%ty_nk(t), this%ty_k0(t), this%nhp_all, &
                        this%mp_z, this%mp_ax, this%mp_ay, this%mp_az, &
                        this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, this%mp_ip)
      end do
      if (this%has_reps) then
         ! The copies' bound I_P from the plain c, THEN the orbit size into c:
         ! the weight is a factor of the bra stage, not of the screen.
         do t = 1, this%ntype
            if (this%rty_nk(t) == 0) cycle
            call mdj_bound(this%ty_la(t), this%ty_lb(t), this%rty_nk(t), this%rty_k0(t), this%nhp_all, &
                           this%mp_z, this%mp_ax, this%mp_ay, this%mp_az, &
                           this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, this%mp_ip)
         end do
         call weight_copies(this)
      end if
      this%built = .true.
   end subroutine mdj_build

   subroutine mdj_coulomb(this, nao, dmat, jmat, rank, nranks, skeleton)
      !! jmat = J[dmat], both (nao, nao), Cartesian, already on the device.
      !! Every element outside the kept shell pairs is set to zero.
      !!
      !! With rank and nranks, this rank's share only: every nranks-th ket
      !! chunk of each Hermite launch. The bra stage is linear in Jt, so the
      !! shares' jmat sum over ranks to J. The ket stage is replicated (it is
      !! linear in pairs, the Hermite stage quadratic).
      class(trc_mdj_t), intent(inout) :: this
      integer,  intent(in)  :: nao
      real(dp), intent(in)  :: dmat(nao, nao)
      real(dp), intent(inout) :: jmat(nao, nao)
      integer,  intent(in), optional :: rank, nranks
      logical,  intent(in), optional :: skeleton
      !! Only the bra of the orbit representatives, each block times its orbit
      !! size: jmat is then the SKELETON J, which the caller's symmetrization
      !! turns into J (docs/SYMMETRY_FOCK.md 4.5). The ket side stays whole.
      !! Needs `build` to have been given the orbit weights. Default false.

      integer :: t, lp, lq, lpt, qfix, c0, cstep
      logical :: skel, sym

      skel = .false.
      if (present(skeleton)) skel = skeleton
      if (skel .and. .not. this%has_reps) then
        error stop "terco: the skeleton MD J needs mdj%build to be given the orbit weights"
      end if
      this%npq = 0
      qfix = qlen_override()
      c0 = 0
      cstep = 1
      if (present(rank) .and. present(nranks)) then
         c0 = rank
         cstep = max(1, nranks)
      end if
      ! NVTX ranges (TERCO_NVTX builds; no-ops otherwise), each closed after a
      ! device sync so its length is the stage's GPU time.
      call profiler_start("mdj: ket", nvtx_only=.true.)
      call zero_device(this%ndt, this%jt)
      call zero_device(nao*nao, jmat)
      do t = 1, this%ntype
         if (this%ty_nk(t) == 0) cycle
         lpt = this%ty_la(t) + this%ty_lb(t)
         call mdj_ket(this%ty_la(t), this%ty_lb(t), this%ty_nk(t), this%ty_k0(t), &
                      this%ty_k0(t) - this%cl_k0(lpt), this%cl_nk(lpt), this%cl_hb(lpt), this%nhp_all, &
                      this%mp_z, this%mp_ax, this%mp_ay, this%mp_az, &
                      this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, this%mp_sp, &
                      this%nsp, this%sp_ao1, this%sp_ao2, nao, dmat, this%ndt, this%dt, this%mp_nq)
      end do
      !$acc wait
      call profiler_stop("mdj: ket")
      call profiler_start("mdj: hermite", nvtx_only=.true.)
      do lp = 0, MDJ_LPMAX
         if (this%cl_nk(lp) == 0) cycle
         do lq = 0, MDJ_LPMAX
            if (this%cl_nk(lq) == 0) cycle
            if (skel) then
               ! The skeleton bra is the orbit representatives' weighted
               ! copies, not the ket, so there is no P <-> Q symmetry to use:
               ! one gather launch per (lp, lq).
               if (this%cr_nk(lp) == 0) cycle
               call mdj_hermite(lp, lq, this%cr_nk(lp), this%cr_k0(lp), this%cr_hb(lp), &
                                this%cl_nk(lq), this%cl_k0(lq), this%cl_hb(lq), &
                                chunk_length(this%cr_nk(lp), this%cl_nk(lq), .false., qfix, cstep), c0, cstep, &
                                this%nhp_all, &
                                this%mp_z, this%mp_px, this%mp_py, this%mp_pz, &
                                this%mp_ip, this%mp_nq, this%ndt, this%dt, this%jt, this%cut)
               this%npq = this%npq + int(this%cr_nk(lp), int64)*this%cl_nk(lq)
            else
               ! A symmetric class is one launch for lq <= lp, filling both
               ! blocks' Jt; the others are a gather launch for each order.
               sym = lp + lq <= this%sym_lmax
               if (sym .and. lq > lp) cycle
               if (sym) then
                  call mdj_hermite_sym(lp, lq, this%cl_nk(lp), this%cl_k0(lp), this%cl_hb(lp), &
                                       this%cl_nk(lq), this%cl_k0(lq), this%cl_hb(lq), &
                                       chunk_length(this%cl_nk(lp), this%cl_nk(lq), lp == lq, qfix, cstep), &
                                       c0, cstep, this%nhp_all, &
                                       this%mp_z, this%mp_px, this%mp_py, this%mp_pz, &
                                       this%mp_ip, this%mp_nq, this%ndt, this%dt, this%jt, this%cut)
                  ! `npq` counts the (P, Q) pairs whose Jt is filled, so a
                  ! symmetric launch counts both blocks, as the two gather
                  ! launches it replaces would: the full build's count stays
                  ! comparable with the skeleton's (check_sym_fock).
                  if (lp == lq) then
                     this%npq = this%npq + int(this%cl_nk(lp), int64)*this%cl_nk(lq)
                  else
                     this%npq = this%npq + 2*int(this%cl_nk(lp), int64)*this%cl_nk(lq)
                  end if
               else
                  call mdj_hermite(lp, lq, this%cl_nk(lp), this%cl_k0(lp), this%cl_hb(lp), &
                                   this%cl_nk(lq), this%cl_k0(lq), this%cl_hb(lq), &
                                   chunk_length(this%cl_nk(lp), this%cl_nk(lq), .false., qfix, cstep), &
                                   c0, cstep, this%nhp_all, &
                                   this%mp_z, this%mp_px, this%mp_py, this%mp_pz, &
                                   this%mp_ip, this%mp_nq, this%ndt, this%dt, this%jt, this%cut)
                  this%npq = this%npq + int(this%cl_nk(lp), int64)*this%cl_nk(lq)
               end if
            end if
         end do
      end do
      !$acc wait
      call profiler_stop("mdj: hermite")
      call profiler_start("mdj: bra", nvtx_only=.true.)
      do t = 1, this%ntype
         lpt = this%ty_la(t) + this%ty_lb(t)
         if (skel) then
            if (this%rty_ns(t) == 0) cycle
            call mdj_bra(this%ty_la(t), this%ty_lb(t), this%rty_ns(t), this%rty_s0(t), &
                         this%cr_k0(lpt), this%cr_nk(lpt), this%cr_hb(lpt), this%nhp_all, &
                         this%mp_z, this%mp_ax, this%mp_ay, this%mp_az, &
                         this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, &
                         this%nrp, this%rsp_ao1, this%rsp_ao2, this%rsp_k0, this%rsp_nk, &
                         this%ndt, this%jt, nao, jmat)
         else
            if (this%ty_ns(t) == 0) cycle
            call mdj_bra(this%ty_la(t), this%ty_lb(t), this%ty_ns(t), this%ty_s0(t), &
                         this%cl_k0(lpt), this%cl_nk(lpt), this%cl_hb(lpt), this%nhp_all, &
                         this%mp_z, this%mp_ax, this%mp_ay, this%mp_az, &
                         this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, &
                         this%nsp, this%sp_ao1, this%sp_ao2, this%sp_k0, this%sp_nk, &
                         this%ndt, this%jt, nao, jmat)
         end if
      end do
      !$acc wait
      call profiler_stop("mdj: bra")
   end subroutine mdj_coulomb

   subroutine order_pairs(this, b, bins)
      !! Group the bins' shell pairs into (la, lb) types, ordered by
      !! (la + lb, la) so each lp block of primitive pairs is contiguous.
      type(trc_mdj_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(pair_bins_t), intent(in) :: bins

      integer :: la, lb, lp, t, s, n
      integer, allocatable :: tyof(:, :)

      allocate (tyof(0:MDJ_LPMAX, 0:MDJ_LPMAX))
      tyof = 0
      this%ntype = 0
      do lp = 0, MDJ_LPMAX
         do la = lp, 0, -1
            lb = lp - la
            if (lb > la .or. la > MDJ_LPMAX/2) cycle
            this%ntype = this%ntype + 1
            tyof(la, lb) = this%ntype
         end do
      end do
      allocate (this%ty_la(this%ntype), this%ty_lb(this%ntype), this%ty_s0(this%ntype), &
                this%ty_ns(this%ntype), this%ty_k0(this%ntype), this%ty_nk(this%ntype))
      do la = 0, MDJ_LPMAX/2
         do lb = 0, la
            t = tyof(la, lb)
            this%ty_la(t) = la
            this%ty_lb(t) = lb
         end do
      end do
      this%ty_ns = 0
      do s = 1, bins%npair
         t = tyof(b%sh_l(bins%sp_i(s)), b%sh_l(bins%sp_j(s)))
         this%ty_ns(t) = this%ty_ns(t) + 1
      end do
      n = 0
      do t = 1, this%ntype
         this%ty_s0(t) = n
         n = n + this%ty_ns(t)
      end do
      deallocate (tyof)
   end subroutine order_pairs

   subroutine fill_prims(this, b, bins, thresh, pair_w)
      !! Shell pairs into type order and their primitive pairs after them.
      !!
      !! A primitive pair is kept on the HGP kernels' own rule
      !! (hgp_build_pairs): its s-type sqrt((ab|ab)) above thresh 10^-3, so
      !! the two engines see the same pairs. A shell pair with none left is
      !! dropped; its J block is zero to that threshold.
      type(trc_mdj_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(pair_bins_t), intent(in) :: bins
      real(dp), intent(in) :: thresh
      integer, intent(in), optional :: pair_w(:)
      !! see mdj_build

      real(dp), parameter :: SQRT2_PI54 = 1.4142135623730951_dp*4.1827004988466890_dp
      integer :: s, t, i, j, ki, kj, k, n, lp, ia, ib, nrep_hp
      integer, allocatable :: pos(:), sord(:), nk(:)
      real(dp) :: ea, eb, p, mu, ab2, cut
      real(dp) :: pc(3)

      cut = thresh*1.0e-3_dp
      ! Type order of the bins' pairs, keeping the bins' order within a type.
      allocate (pos(this%ntype), sord(bins%npair), nk(bins%npair))
      pos = this%ty_s0
      do s = 1, bins%npair
         t = type_of(this, b%sh_l(bins%sp_i(s)), b%sh_l(bins%sp_j(s)))
         pos(t) = pos(t) + 1
         sord(pos(t)) = s
      end do
      ! Count the kept primitive pairs of each.
      this%nhp = 0
      do n = 1, bins%npair
         s = sord(n)
         i = bins%sp_i(s)
         j = bins%sp_j(s)
         ab2 = sum((b%sh_r(:, i) - b%sh_r(:, j))**2)
         nk(n) = 0
         do ki = 1, b%sh_np(i)
            do kj = 1, b%sh_np(j)
               if (keep(b%sh_e(ki, i), b%sh_e(kj, j), b%sh_c(ki, i), b%sh_c(kj, j), ab2, cut)) then
                 nk(n) = nk(n) + 1
               end if
            end do
         end do
         this%nhp = this%nhp + nk(n)
      end do

      ! Drop the shell pairs with nothing kept, renumbering within each type.
      this%nsp = count(nk > 0)
      this%has_reps = present(pair_w)
      nrep_hp = 0
      if (this%has_reps) then
         do n = 1, bins%npair
            if (pair_w(sord(n)) > 0) nrep_hp = nrep_hp + nk(n)
         end do
      end if
      this%nhp_all = this%nhp + nrep_hp
      allocate (this%sp_ao1(max(this%nsp, 1)), this%sp_ao2(max(this%nsp, 1)), &
                this%sp_k0(max(this%nsp, 1)), this%sp_nk(max(this%nsp, 1)), &
                this%sp_sh1(max(this%nsp, 1)), this%sp_sh2(max(this%nsp, 1)), this%sp_q(max(this%nsp, 1)), &
                this%sp_as1(max(this%nsp, 1)), this%sp_as2(max(this%nsp, 1)))
      this%nao_sph = b%nao_sph
      allocate (this%mp_z(max(this%nhp_all, 1)), this%mp_px(max(this%nhp_all, 1)), &
                this%mp_py(max(this%nhp_all, 1)), this%mp_pz(max(this%nhp_all, 1)), &
                this%mp_ax(max(this%nhp_all, 1)), this%mp_ay(max(this%nhp_all, 1)), &
                this%mp_az(max(this%nhp_all, 1)), this%mp_bx(max(this%nhp_all, 1)), &
                this%mp_by(max(this%nhp_all, 1)), this%mp_bz(max(this%nhp_all, 1)), &
                this%mp_c(max(this%nhp_all, 1)), this%mp_sp(max(this%nhp_all, 1)))
      allocate (this%sp_w(max(this%nsp, 1)))
      this%sp_w = 0
      ia = 0
      k = 0
      do t = 1, this%ntype
         ib = ia
         this%ty_k0(t) = k
         do n = this%ty_s0(t) + 1, this%ty_s0(t) + this%ty_ns(t)
            if (nk(n) == 0) cycle
            s = sord(n)
            i = bins%sp_i(s)
            j = bins%sp_j(s)
            ia = ia + 1
            this%sp_ao1(ia) = b%sh_ao(i)
            this%sp_ao2(ia) = b%sh_ao(j)
            this%sp_k0(ia) = k
            this%sp_nk(ia) = nk(n)
            this%sp_sh1(ia) = i
            this%sp_sh2(ia) = j
            this%sp_as1(ia) = b%sh_ao_sph(i)
            this%sp_as2(ia) = b%sh_ao_sph(j)
            this%sp_q(ia) = bins%sp_q(s)
            if (present(pair_w)) this%sp_w(ia) = pair_w(s)
            ab2 = sum((b%sh_r(:, i) - b%sh_r(:, j))**2)
            do ki = 1, b%sh_np(i)
               ea = b%sh_e(ki, i)
               do kj = 1, b%sh_np(j)
                  eb = b%sh_e(kj, j)
                  if (.not. keep(ea, eb, b%sh_c(ki, i), b%sh_c(kj, j), ab2, cut)) cycle
                  k = k + 1
                  p = ea + eb
                  mu = ea*eb/p
                  pc = (ea*b%sh_r(:, i) + eb*b%sh_r(:, j))/p
                  this%mp_z(k) = p
                  this%mp_px(k) = pc(1)
                  this%mp_py(k) = pc(2)
                  this%mp_pz(k) = pc(3)
                  this%mp_ax(k) = pc(1) - b%sh_r(1, i)
                  this%mp_ay(k) = pc(2) - b%sh_r(2, i)
                  this%mp_az(k) = pc(3) - b%sh_r(3, i)
                  this%mp_bx(k) = pc(1) - b%sh_r(1, j)
                  this%mp_by(k) = pc(2) - b%sh_r(2, j)
                  this%mp_bz(k) = pc(3) - b%sh_r(3, j)
                  this%mp_c(k) = b%sh_c(ki, i)*b%sh_c(kj, j)*exp(-mu*ab2)
                  this%mp_sp(k) = ia
               end do
            end do
         end do
         this%ty_s0(t) = ib
         this%ty_ns(t) = ia - ib
         this%ty_nk(t) = k - this%ty_k0(t)
      end do

      this%cl_k0 = 0
      this%cl_nk = 0
      do t = 1, this%ntype
         lp = this%ty_la(t) + this%ty_lb(t)
         if (this%cl_nk(lp) == 0) this%cl_k0(lp) = this%ty_k0(t)
         this%cl_nk(lp) = this%cl_nk(lp) + this%ty_nk(t)
      end do
      call copy_reps(this)
      deallocate (pos, sord, nk)

   contains

      pure logical function keep(a, bb, ca, cb, r2, c) result(res)
         real(dp), intent(in) :: a, bb, ca, cb, r2, c
         real(dp) :: pp
         pp = a + bb
         res = abs(ca*cb)*exp(-(a*bb/pp)*r2)*SQRT2_PI54/(pp*sqrt(sqrt(2.0_dp*pp))) > c
      end function keep
   end subroutine fill_prims

   subroutine copy_reps(this)
      !! The orbit representatives' primitive pairs, again, after the nhp
      !! ordinary ones: type by type, so each lp block of copies is
      !! contiguous, with their own shell-pair arrays. A copy's c is the
      !! plain one until weight_copies (the bound needs the plain one).
      type(trc_mdj_t), intent(inout) :: this
      integer :: t, lp, ia, kk, k, kn, ir
      allocate (this%rty_s0(this%ntype), this%rty_ns(this%ntype), this%rty_k0(this%ntype), this%rty_nk(this%ntype))
      this%rty_s0 = 0
      this%rty_ns = 0
      this%rty_k0 = this%nhp
      this%rty_nk = 0
      this%nrp = 0
      if (this%has_reps) this%nrp = count(this%sp_w(1:this%nsp) > 0)
      allocate (this%rsp_ao1(max(this%nrp, 1)), this%rsp_ao2(max(this%nrp, 1)), &
                this%rsp_k0(max(this%nrp, 1)), this%rsp_nk(max(this%nrp, 1)), this%rsp_w(max(this%nrp, 1)))
      this%cr_k0 = 0
      this%cr_nk = 0
      if (.not. this%has_reps) return
      ir = 0
      kn = this%nhp
      do t = 1, this%ntype
         this%rty_s0(t) = ir
         this%rty_k0(t) = kn
         do ia = this%ty_s0(t) + 1, this%ty_s0(t) + this%ty_ns(t)
            if (this%sp_w(ia) == 0) cycle
            ir = ir + 1
            this%rsp_ao1(ir) = this%sp_ao1(ia)
            this%rsp_ao2(ir) = this%sp_ao2(ia)
            this%rsp_nk(ir) = this%sp_nk(ia)
            this%rsp_w(ir) = this%sp_w(ia)
            this%rsp_k0(ir) = kn
            do kk = 1, this%sp_nk(ia)
               k = this%sp_k0(ia) + kk
               kn = kn + 1
               this%mp_z(kn) = this%mp_z(k)
               this%mp_px(kn) = this%mp_px(k)
               this%mp_py(kn) = this%mp_py(k)
               this%mp_pz(kn) = this%mp_pz(k)
               this%mp_ax(kn) = this%mp_ax(k)
               this%mp_ay(kn) = this%mp_ay(k)
               this%mp_az(kn) = this%mp_az(k)
               this%mp_bx(kn) = this%mp_bx(k)
               this%mp_by(kn) = this%mp_by(k)
               this%mp_bz(kn) = this%mp_bz(k)
               this%mp_c(kn) = this%mp_c(k)
               this%mp_sp(kn) = ir
            end do
         end do
         this%rty_ns(t) = ir - this%rty_s0(t)
         this%rty_nk(t) = kn - this%rty_k0(t)
      end do
      do t = 1, this%ntype
         lp = this%ty_la(t) + this%ty_lb(t)
         if (this%cr_nk(lp) == 0) this%cr_k0(lp) = this%rty_k0(t)
         this%cr_nk(lp) = this%cr_nk(lp) + this%rty_nk(t)
      end do
   end subroutine copy_reps

   subroutine weight_copies(this)
      !! The orbit size into each copy's c, host and device, after the bound
      !! was formed from the plain c.
      type(trc_mdj_t), intent(inout) :: this
      integer :: ir, kk
      do ir = 1, this%nrp
         do kk = 1, this%rsp_nk(ir)
            this%mp_c(this%rsp_k0(ir) + kk) = this%mp_c(this%rsp_k0(ir) + kk)*real(this%rsp_w(ir), dp)
         end do
      end do
      !$acc update device(this%mp_c(this%nhp+1:this%nhp_all))
   end subroutine weight_copies

   pure integer function chunk_length(np, nq, half, qfix, nshare) result(qlen)
      !! Kets per thread for an np x nq Hermite launch: enough chunks to reach
      !! MDJ_THREADS_TARGET threads on each of nshare ranks, never under
      !! MDJ_QLEN_MIN kets each. A symmetric launch of one block against
      !! itself (half) gives each P only np/2 + 1 kets.
      integer, intent(in) :: np, nq, qfix, nshare
      logical, intent(in) :: half
      integer :: nchunk, nk
      if (qfix > 0) then
         qlen = qfix
         return
      end if
      nk = nq
      if (half) nk = nq/2 + 1
      nchunk = max(1, (MDJ_THREADS_TARGET*nshare + np - 1)/np)
      qlen = max(MDJ_QLEN_MIN, (nk + nchunk - 1)/nchunk)
   end function chunk_length

   integer function qlen_override() result(res)
      !! TRC_MDJ_QLEN fixes the kets per chunk for every launch, for sweeping
      !! it; unset or not positive leaves chunk_length to choose.
      character(len=32) :: e
      integer :: ios
      e = " "
      res = 0
      call get_environment_variable("TRC_MDJ_QLEN", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) res
         if (ios /= 0) res = 0
      end if
   end function qlen_override

   pure integer function type_of(this, la, lb) result(t)
      type(trc_mdj_t), intent(in) :: this
      integer, intent(in) :: la, lb
      do t = 1, this%ntype
         if (this%ty_la(t) == la .and. this%ty_lb(t) == lb) return
      end do
      t = 0
   end function type_of

   pure integer function nherm(l) result(n)
      integer, intent(in) :: l
      n = (l + 1)*(l + 2)*(l + 3)/6
   end function nherm

   real(dp) function margin() result(res)
      !! MDJ_MARGIN_DEFAULT unless TRC_MDJ_MARGIN says otherwise; 0 turns
      !! the (P, Q) screen off, for measuring what it costs in accuracy.
      character(len=32) :: e
      integer :: ios
      e = " "
      call get_environment_variable("TRC_MDJ_MARGIN", e)
      res = MDJ_MARGIN_DEFAULT
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) res
         if (ios /= 0) res = MDJ_MARGIN_DEFAULT
      end if
   end function margin

   integer function sym_lmax_choice() result(res)
      !! MDJ_SYM_LMAX, or less with TRC_MDJ_SYM_LMAX in a build whose every
      !! class also has its gather kernels (gen_md_j.py --gather-all), for
      !! racing the two forms; -1 is all gather.
      character(len=32) :: e
      integer :: ios, v
      res = MDJ_SYM_LMAX
      if (.not. MDJ_GATHER_ALL) return
      e = " "
      call get_environment_variable("TRC_MDJ_SYM_LMAX", e)
      if (len_trim(e) > 0) then
         read (e, *, iostat=ios) v
         if (ios == 0) res = max(-1, min(v, MDJ_SYM_LMAX))
      end if
   end function sym_lmax_choice

   subroutine zero_device(n, a)
      integer,  intent(in)    :: n
      real(dp), intent(inout) :: a(n)
      integer :: i
      do concurrent(i = 1:n)
         a(i) = 0.0_dp
      end do
   end subroutine zero_device

   subroutine mdj_release(this)
      class(trc_mdj_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%mp_z, this%mp_px, this%mp_py, this%mp_pz, &
         !$acc                  this%mp_ax, this%mp_ay, this%mp_az, &
         !$acc                  this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, this%mp_sp, &
         !$acc                  this%sp_ao1, this%sp_ao2, this%sp_k0, this%sp_nk, &
         !$acc                  this%sp_sh1, this%sp_sh2, this%sp_q, this%sp_as1, this%sp_as2, &
         !$acc                  this%rsp_ao1, this%rsp_ao2, this%rsp_k0, this%rsp_nk, &
         !$acc                  this%mp_ip, this%mp_nq, this%dt, this%jt)
         this%on_device = .false.
      end if
      if (allocated(this%mp_z)) deallocate (this%mp_z, this%mp_px, this%mp_py, this%mp_pz, &
                                            this%mp_ax, this%mp_ay, this%mp_az, &
                                            this%mp_bx, this%mp_by, this%mp_bz, this%mp_c, this%mp_sp)
      if (allocated(this%sp_ao1)) deallocate (this%sp_ao1, this%sp_ao2, this%sp_k0, this%sp_nk, &
                                              this%sp_sh1, this%sp_sh2, this%sp_q, this%sp_as1, this%sp_as2)
      if (allocated(this%sp_w)) deallocate (this%sp_w)
      if (allocated(this%rsp_ao1)) deallocate (this%rsp_ao1, this%rsp_ao2, this%rsp_k0, this%rsp_nk, this%rsp_w)
      if (allocated(this%rty_s0)) deallocate (this%rty_s0, this%rty_ns, this%rty_k0, this%rty_nk)
      if (allocated(this%mp_ip)) deallocate (this%mp_ip, this%mp_nq)
      if (allocated(this%dt)) deallocate (this%dt, this%jt)
      if (allocated(this%ty_la)) deallocate (this%ty_la, this%ty_lb, this%ty_s0, this%ty_ns, &
                                             this%ty_k0, this%ty_nk)
      this%nhp = 0
      this%nhp_all = 0
      this%nsp = 0
      this%nrp = 0
      this%ndt = 0
      this%ntype = 0
      this%has_reps = .false.
      this%built = .false.
   end subroutine mdj_release

end module trc_mdj
