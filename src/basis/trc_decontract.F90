!
! Splitting a general contraction into primitive shells, and carrying the
! contraction as a matrix transform instead.
!
module trc_decontract
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use trc_basis_types, only: trc_basis_t
   implicit none
   private

   public :: decon_t, decontract_basis, decontract_expand, decontract_fold, decontract_release
   public :: decontract_expand_host, decontract_fold_host
   public :: trc_decontract_options_t, decontract_whole_screen
   public :: decontract_options_for_exchange, decontract_options_equal
   public :: DECON_CUTOFF_DEFAULT, DECON_DIFFUSE_DEFAULT, DECON_SHELL_COST_DEFAULT

   real(dp), parameter :: DECON_CUTOFF_DEFAULT = 0.3_dp
   !! Exponents at or below this are diffuse. GPU4PySCF's `diffuse_cutoff`.

   logical, parameter :: DECON_DIFFUSE_DEFAULT = .false.
   !! SPLIT DIFFUSE PRIMITIVES OUT OF SEGMENTED SHELLS -- OFF BY DEFAULT.
   !!
   !! GPU4PySCF moves every primitive at or below the cutoff out of a
   !! contracted shell into a shell of its own, so the tight remainder
   !! screens on its own bound. In pyscf's own data for H, C, N, O and Si
   !! (what trc_basis_json's pruning turns the BSE files into) the only
   !! contracted shells holding a primitive at or below 0.3 are Si's 6-31G*
   !! sp; at 0.5 H's cc-pVDZ and def2-SVPD s join them; at 0.1 nothing
   !! does. It loses wherever it acts (V100, RHF, 1e-11, full builds, mean
   !! of builds 2-10): msn/6-31G* 3.76 -> 4.32 s at 0.3; gly20/cc-pVDZ
   !! 1.31 -> 1.41, gly20/def2-SVPD 8.82 -> 9.48, water40/aug-cc-pVDZ
   !! 25.8 -> 29.7 s at 0.5. Every class with an s shell grows by 8-11%
   !! (nsys, gly20/cc-pVDZ): more shells, each split one carrying a
   !! screening amplitude of 2 (`build_maps`), so its pairs are judged four
   !! times tighter, and the transform of D and G every build.

   real(dp), parameter :: DECON_WHOLE_UNEVEN_MAX = 128.0_dp
   !! A general contraction is kept whole only if no column's coefficient
   !! falls more than this far below the largest the group gives that
   !! primitive (group_unevenness): the kernels' primitive cut is then
   !! tightened by at most 128**4 = 2.7e8 (`decontract_whole_screen`). In
   !! cc-pVDZ that keeps O's s (36 in pyscf's data, 71 in the pruned BSE
   !! file), Si p (27), P p (32), S p (55) and the Fe/Zn d (6.3) whole, and
   !! splits every second-row s (700-1000), Cl p (319) and the
   !! transition-metal s and p.

   real(dp), parameter :: DECON_SHELL_COST_DEFAULT = 2.0_dp
   !! THE PRICE OF ONE MORE SHELL, IN PRIMITIVES, FOR THE GENERAL-CONTRACTION
   !! SPLIT. A general contraction is split only when
   !!
   !!   sum over the shells it becomes of (nprim + shell_cost)
   !!
   !! is smaller than the same sum over the contracted shells: the split
   !! removes the primitives the columns repeat and adds shells, and every
   !! shell costs pair and quartet bookkeeping and a Schwarz bound carrying
   !! the transform's amplitude (`build_maps`, up to 12 on cc-pVDZ oxygen),
   !! whatever its primitive count. At 2 the 8-primitive, 2-column s of C,
   !! N and O (16 + 2c against 8 + 8c) stays whole. A group too uneven to
   !! keep (DECON_WHOLE_UNEVEN_MAX) splits whatever the price, and
   !! `md_coulomb` splits every group. 0 splits every general contraction.
   !! That is not quite what this module did before the price existed: a
   !! group with exactly one multi-primitive column then kept that whole
   !! column, shared primitives included, where PLAN_SEGMENT keeps only the
   !! primitives no other column uses.
   !!
   !! WHAT A SHELL KEPT WHOLE GUARANTEES (`screen_whole`, on by default).
   !! Three screens act. Two of them act on primitives, and for those the
   !! guarantee is by construction (`decontract_whole_screen`): the pair
   !! prune and the kernels' primitive-quartet cut drop nothing the split
   !! would have kept. The third is the shell-pair Schwarz test. There the
   !! bound carries the amplitude the split's primitives would have
   !! (`whole_amplitude`). That matches the two screens' strictness but is
   !! not a proof, so the claim there is measured, not derived. Converged
   !! RHF at 1e-11, against the same system at 1e-14, kept whole is as
   !! accurate as all split: see the commit that set DECON_WHOLE_UNEVEN_MAX
   !! for the table (water, glycine, silica, P/S/Cl and Fe/Zn systems), and
   !! check_decontract_policy for the loose-threshold case that fails
   !! without `screen_whole`.
   !! Measured against 0, which is main before the price (V100, RHF,
   !! 1e-11, full builds, mean of builds 2-10, s per build): gly20/cc-pVDZ
   !! 1.945 -> 1.484, aug-cc-pVDZ 27.67 -> 19.24; water40/cc-pVDZ 0.990 ->
   !! 0.849, aug-cc-pVDZ 36.38 -> 27.08; NA-040/cc-pVDZ 2.101 -> 1.764,
   !! aug-cc-pVDZ 48.48 -> 34.63; w172/cc-pVDZ 16.83 -> 11.37; the silica
   !! fragment and slice, cc-pVDZ, 9.66 -> 7.67 and 6.85 -> 5.56. Without
   !! the primitive-level screen (`decontract_whole_screen`) these were
   !! 1.370, 18.71, 0.788, 26.58, 1.586, 33.49, 10.78, 7.47 and 5.64, at
   !! accuracy that depended on the policy.

   type :: trc_decontract_options_t
      !! How the working basis is split. TRC_DECON_CUTOFF, TRC_DECON_DIFFUSE
      !! (0 or 1), TRC_DECON_SHELL_COST and TRC_DECON_SCREEN_WHOLE (0 or 1)
      !! override the fields at build time (`trc_eri`, rank 0's environment
      !! on every rank), so a policy can be swept without a rebuild. The
      !! integrals over the split basis are exact either way; what changes
      !! is the cost and what the screens drop (DECON_SHELL_COST_DEFAULT).
      real(dp) :: cutoff = DECON_CUTOFF_DEFAULT
      !! Exponents at or below this are diffuse.
      logical :: diffuse = DECON_DIFFUSE_DEFAULT
      !! Move diffuse primitives out of segmented shells into their own.
      real(dp) :: shell_cost = DECON_SHELL_COST_DEFAULT
      !! See DECON_SHELL_COST_DEFAULT.
      logical :: screen_whole = .true.
      !! A general contraction kept whole is screened as its split would be:
      !! its primitives carry the split's coefficients in the pair prune and
      !! the kernels' primitive cut is tightened to match
      !! (`decontract_whole_screen`), and its shell pairs carry the split's
      !! amplitude (`whole_amplitude`). See DECON_SHELL_COST_DEFAULT for
      !! what that guarantees. TRC_DECON_SCREEN_WHOLE=0 turns it off, which
      !! is the contracted shell's own screen and is measurably less
      !! accurate at the same threshold.
      logical :: md_coulomb = .false.
      !! The Coulomb matrix will come from the McMurchie-Davidson J engine
      !! (trc_mdj: a basis with f or higher, or a functional without exact
      !! exchange). Its cost goes with primitive pairs, not shells, so every
      !! general contraction is split, whatever `shell_cost` says. Measured
      !! on gly20/cc-pVTZ (V100, RHF, 1e-11, nsys over a 3-build run):
      !! keeping C, N and O's s whole takes the MD J kernels from 6.5 to
      !! 11.4 s and the exchange bins from 28.8 to 29.4 s, 17.1 -> 19.3 s
      !! per build; water20, water40 and NA-040 in cc-pVTZ lose 1-5% the
      !! same way.
   end type trc_decontract_options_t

   integer, parameter :: PLAN_WHOLE = 0
   !! Keep every member of the group as it is.
   integer, parameter :: PLAN_SEGMENT = 1
   !! One column keeps the primitives only it uses as a segmented shell;
   !! every other primitive becomes a shell of one, shared by the columns.
   integer, parameter :: PLAN_PRIMITIVES = 2
   !! Every primitive becomes a shell of one.

   type :: decon_t
      logical :: active = .false.
      integer :: nao_c = 0, nao_p = 0, nnz = 0
      integer,  allocatable :: f_off(:)   !! (nao_p + 1)
      integer,  allocatable :: f_c(:)     !! (nnz) contracted AO
      real(dp), allocatable :: f_a(:)     !! (nnz) coefficient
      integer,  allocatable :: t_off(:)   !! (nao_c + 1)
      integer,  allocatable :: t_p(:)     !! (nnz) primitive AO
      real(dp), allocatable :: t_a(:)     !! (nnz)
      !
      real(dp), allocatable :: amp(:)     !! (pb%nshell)
      !! Per primitive SHELL, the factor its Schwarz bound must be multiplied
      !! by for screening in this basis to bound the error in the CONTRACTED
      !! matrix. Derived, with the reasoning, at the end of `build_maps`.
      !! Host only -- screening is built once, on the host.
      logical :: screen_whole = .false.
      !! The working basis's general contractions kept whole screen as the
      !! split would (`decontract_whole_screen`). False for a transform that
      !! never went through `decontract_basis`, which is the old screen.
      logical :: on_device = .false.
   end type decon_t

   type :: decon_build_t
    !! The primitive basis and the transform `emit_segment`/`emit_primitive`
    !! fill, one shell/triple at a time, bundled so those two internal
    !! helpers stay under the dummy-argument limit -- `decontract_basis`
    !! owns the allocation (sized once, an upper bound on every primitive
    !! of every shell becoming its own) and reads the result back out
    !! through `nsh_p`/`ntr`, exactly as it did through the separate
    !! arrays before.
      integer, allocatable :: p_l(:), p_np(:)
      real(dp), allocatable :: p_e(:, :), p_c(:, :), p_r(:, :)
      logical, allocatable :: p_single(:)
      !! A shell of one primitive that the columns of its group share, as
      !! opposed to a contracted shell that happens to have one primitive:
      !! only the former may be reused by a later column.
      real(dp), allocatable :: p_screen(:)
      !! Extra screening amplitude per emitted shell: `whole_amplitude` for a
      !! member of a general contraction kept whole, 1 otherwise.
      integer :: nsh_p = 0
      ! one triple per (contracted shell, emitted shell) pair
      integer, allocatable :: tr_mu(:), tr_s(:)
      real(dp), allocatable :: tr_a(:)
      integer :: ntr = 0
   end type decon_build_t

contains

   function decontract_options_for_exchange(base, exx, exx_lr) result(opts)
      !! The policy an engine for a method with these exact-exchange
      !! fractions must be built with: `base`, with `md_coulomb` on when there
      !! is no exact exchange (no bins at all; J is the MD engine's). The one
      !! rule every builder -- the SCF's, a gradient's fresh one, a response
      !! solve's -- applies, so an engine reused from the SCF and one built
      !! fresh are the same engine.
      type(trc_decontract_options_t), intent(in) :: base
      real(dp), intent(in) :: exx, exx_lr
      type(trc_decontract_options_t) :: opts
      opts = base
      if (exx == 0.0_dp .and. exx_lr == 0.0_dp) opts%md_coulomb = .true.
   end function decontract_options_for_exchange

   pure logical function decontract_options_equal(a, b) result(res)
      !! The two policies build the same working basis.
      type(trc_decontract_options_t), intent(in) :: a, b
      res = a%cutoff == b%cutoff .and. (a%diffuse .eqv. b%diffuse) .and. &
            a%shell_cost == b%shell_cost .and. (a%screen_whole .eqv. b%screen_whole) .and. &
            (a%md_coulomb .eqv. b%md_coulomb)
   end function decontract_options_equal

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! Build the primitive basis and the transform.
   !
   ! `active` comes back false for a basis with nothing to split, and then
   ! `pb` and the maps are untouched: the caller must go on using `b`.
   !
   subroutine decontract_basis(b, pb, m, options)
      type(trc_basis_t), intent(in)    :: b
      type(trc_basis_t), intent(inout) :: pb
      type(decon_t),     intent(inout) :: m
      type(trc_decontract_options_t), intent(in), optional :: options
      !! The split policy; the defaults if absent.

      integer, allocatable :: grp(:)          ! group id per shell
      type(decon_build_t) :: bld
      type(trc_decontract_options_t) :: opts
      integer :: t
      integer :: maxnp_p
      integer :: ng, g, i, n
      integer, parameter :: LTOP = 8
      real(dp) :: fl(0:LTOP)

      if (present(options)) opts = options

      call decontract_release(m)
      m%screen_whole = opts%screen_whole
      call group_shells(b, grp, ng)

      ! Nothing shares an exponent with anything, or nothing that does is
      ! worth the split: a segmented basis as far as this module cares. Leave.
      if (.not. any_group_splits(b, grp, ng, opts)) then
         m%active = .false.
         ! The caller keeps its basis, but a general contraction kept whole
         ! still screens as its split would have: `amp` then lives over the
         ! caller's shells, and only when some shell needs it.
         if (opts%screen_whole) then
            allocate (m%amp(b%nshell))
            m%amp = 1.0_dp
            do g = 1, ng
               call group_whole_amplitudes(b, grp, g, m%amp)
            end do
            if (all(m%amp == 1.0_dp)) deallocate (m%amp)
         end if
         deallocate (grp)
         return
      end if

      ! Upper bounds: every primitive of every shell could become its own
      ! shell, and each contributes one triple.
      n = 0
      do i = 1, b%nshell
         n = n + b%sh_np(i)
      end do
      allocate (bld%p_l(n + b%nshell), bld%p_np(n + b%nshell), bld%p_r(3, n + b%nshell))
      allocate (bld%p_e(b%maxnp, n + b%nshell), bld%p_c(b%maxnp, n + b%nshell))
      allocate (bld%tr_mu(n + b%nshell), bld%tr_s(n + b%nshell), bld%tr_a(n + b%nshell))
      allocate (bld%p_single(n + b%nshell), bld%p_screen(n + b%nshell))
      bld%p_screen = 1.0_dp
      bld%p_e = 0.0_dp
      bld%p_c = 0.0_dp
      bld%nsh_p = 0
      bld%ntr = 0

      do g = 1, ng
         call split_group(b, grp, g, opts, bld)
      end do

      maxnp_p = 1
      do i = 1, bld%nsh_p
         maxnp_p = max(maxnp_p, bld%p_np(i))
      end do
      call pb%build(bld%nsh_p, bld%p_l(1:bld%nsh_p), bld%p_np(1:bld%nsh_p), bld%p_e(1:maxnp_p, 1:bld%nsh_p), &
                    bld%p_c(1:maxnp_p, 1:bld%nsh_p), bld%p_r(:, 1:bld%nsh_p), &
                    b%natm, b%at_z, b%at_r, maxnp_p, at_elem=b%at_elem)

      !
      ! BUILD IT TWICE, TO CANCEL THE CONSTRUCTOR'S OWN SCALING.
      !
      ! `basis_build` multiplies every coefficient it is given by
      ! common_fac_sp(l), so a shell handed back its own coefficients comes
      ! out that factor too large and alpha has to be that factor too small.
      ! Correct, but 1/common_fac_sp is 3.5 for s, and a quartet of four such
      ! shells scales the screening error by 3.5^4 = 158: water/cc-pVDZ came
      ! out 1.2 mHa from the reference with every integral individually right.
      !
      ! So learn the factor from the basis just built -- rather than hard
      ! coding a constant that lives in another module and could change --
      ! divide it out, and build again. Now the split shells carry exactly
      ! the coefficients the contracted ones did, every alpha is at most one,
      ! and screening in the primitive basis means what it says.
      !
      fl = 1.0_dp
      do i = 1, bld%nsh_p
         if (bld%p_c(1, i) /= 0.0_dp) fl(bld%p_l(i)) = pb%sh_c(1, i)/bld%p_c(1, i)
      end do
      do i = 1, bld%nsh_p
         if (fl(bld%p_l(i)) /= 0.0_dp) bld%p_c(:, i) = bld%p_c(:, i)/fl(bld%p_l(i))
      end do
      call pb%build(bld%nsh_p, bld%p_l(1:bld%nsh_p), bld%p_np(1:bld%nsh_p), bld%p_e(1:maxnp_p, 1:bld%nsh_p), &
                    bld%p_c(1:maxnp_p, 1:bld%nsh_p), bld%p_r(:, 1:bld%nsh_p), &
                    b%natm, b%at_z, b%at_r, maxnp_p, at_elem=b%at_elem)

      ! The coefficient a split shell ENDS UP with is not the one handed to
      ! the constructor: `basis_build` folds common_fac_sp into it. Deriving
      ! alpha from the built basis rather than from the input is exact
      ! whatever the constructor does to the numbers on the way through --
      ! computing it beforehand left every s and p shell short by that
      ! factor, which is a 4e-9 error in G and a silent one.
      do t = 1, bld%ntr
         bld%tr_a(t) = bld%tr_a(t)/pb%sh_c(1, bld%tr_s(t))
      end do

      call build_maps(b, pb, bld%tr_mu, bld%tr_s, bld%tr_a, bld%ntr, m)
      if (opts%screen_whole) m%amp = m%amp*bld%p_screen(1:bld%nsh_p)
      m%active = .true.
      block
         character(len=8) :: dbge
         integer :: c, k
         real(dp) :: amx, asum, smx
         dbge = " "
         call get_environment_variable("TRC_BUILD_TIMING", dbge)
         if (len_trim(dbge) > 0) then
            amx = 0.0_dp
            smx = 0.0_dp
            do c = 1, m%nao_c
               asum = 0.0_dp
               do k = m%t_off(c), m%t_off(c + 1) - 1
                  amx = max(amx, abs(m%t_a(k)))
                  asum = asum + abs(m%t_a(k))
               end do
               smx = max(smx, asum)
            end do
            block
               character(len=1024) :: logline
               write (logline, "(a,i0,a,i0,a,i0,a,i0,a,es10.3,a,es10.3,a,f6.2)") &
                  "  [decon] shells ", b%nshell, " -> ", pb%nshell, &
                  "   nao ", m%nao_c, " -> ", m%nao_p, &
                  "   max|C| ", amx, "   max col sum|C| ", smx, &
                  "   max amp ", maxval(m%amp)
               call global_logger%info(trim(logline))
            end block
         end if
      end block
      deallocate (bld%p_l, bld%p_np, bld%p_e, bld%p_c, bld%p_r, bld%tr_mu, bld%tr_s, bld%tr_a, bld%p_single, bld%p_screen, grp)
   end subroutine decontract_basis

   !
   ! Shells at one centre and one angular momentum that share an exponent,
   ! transitively, are one group. Sharing is what a general contraction is;
   ! a segmented basis puts every shell in a group of its own.
   !
   subroutine group_shells(b, grp, ng)
      type(trc_basis_t), intent(in) :: b
      integer, allocatable, intent(out) :: grp(:)
      integer, intent(out) :: ng

      integer :: i, j, gi, gj
      integer, allocatable :: root(:)

      allocate (root(b%nshell), grp(b%nshell))
      do i = 1, b%nshell
         root(i) = i
      end do
      do i = 1, b%nshell
         do j = i + 1, b%nshell
            if (b%sh_l(i) /= b%sh_l(j)) cycle
            if (maxval(abs(b%sh_r(:, i) - b%sh_r(:, j))) > 1.0e-10_dp) cycle
            if (.not. shares_exponent(b, i, j)) cycle
            gi = find(root, i)
            gj = find(root, j)
            if (gi /= gj) root(gj) = gi
         end do
      end do
      ! Compact the roots into 1..ng.
      grp = 0
      ng = 0
      do i = 1, b%nshell
         gi = find(root, i)
         if (grp(gi) == 0) then
            ng = ng + 1
            grp(gi) = ng
         end if
      end do
      do i = 1, b%nshell
         grp(i) = grp(find(root, i))
      end do
      deallocate (root)
   end subroutine group_shells

   recursive integer function find(root, i) result(r)
      integer, intent(inout) :: root(:)
      integer, intent(in) :: i
      if (root(i) == i) then
         r = i
      else
         r = find(root, root(i))
         root(i) = r
      end if
   end function find

   logical function shares_exponent(b, i, j) result(res)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: i, j
      integer :: ki, kj
      res = .false.
      do ki = 1, b%sh_np(i)
         do kj = 1, b%sh_np(j)
            if (same_exp(b%sh_e(ki, i), b%sh_e(kj, j))) then
               res = .true.
               return
            end if
         end do
      end do
   end function shares_exponent

   pure logical function same_exp(a, c) result(res)
      real(dp), intent(in) :: a, c
      res = abs(a - c) <= 1.0e-10_dp*max(1.0_dp, abs(a))
   end function same_exp

   !
   ! Would any group actually be split? If not the caller keeps the basis it
   ! has, and a segmented run is bit-for-bit what it was before this module
   ! existed.
   !
   logical function any_group_splits(b, grp, ng, opts) result(res)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ng
      integer, intent(in) :: grp(:)
      type(trc_decontract_options_t), intent(in) :: opts
      integer :: g
      res = .false.
      do g = 1, ng
         if (group_plan(b, grp, g, opts) /= PLAN_WHOLE) then
            res = .true.
            return
         end if
      end do
   end function any_group_splits

   !
   ! What one group becomes: PLAN_WHOLE, PLAN_SEGMENT or PLAN_PRIMITIVES.
   !
   ! A group of one shell is a segmented shell. It is split only when
   ! `opts%diffuse` asks for it and it has a diffuse primitive among several:
   ! the compact primitives stay one segment (PLAN_SEGMENT), the diffuse
   ! ones go out on their own, GPU4PySCF's rule for `nctr == 1`.
   !
   ! A general contraction -- more than one shell sharing exponents -- takes
   ! GPU4PySCF's `_optimize_contraction` rule: if exactly one column holds
   ! more than one compact primitive, that column keeps the primitives no
   ! other column uses as a segment and the rest become shared shells of one
   ! (PLAN_SEGMENT); otherwise every primitive does (PLAN_PRIMITIVES). That
   ! rule only says how to split. Whether to is the price in
   ! DECON_SHELL_COST_DEFAULT: the split happens only if it is cheaper than
   ! the contracted shells by that measure.
   !
   integer function group_plan(b, grp, g, opts) result(plan)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: g
      integer, intent(in) :: grp(:)
      type(trc_decontract_options_t), intent(in) :: opts

      integer :: mem(b%nshell)
      logical :: in_segment(b%maxnp)
      integer :: nmem, i, mu, ncomp, nmulti, i0, nunique, nseg
      real(dp) :: cost_whole, cost_split, price

      call group_members(b, grp, g, mem, nmem)
      plan = PLAN_WHOLE

      if (nmem == 1) then
         mu = mem(1)
         if (.not. opts%diffuse .or. b%sh_np(mu) < 2) return
         ncomp = count(b%sh_e(1:b%sh_np(mu), mu) > opts%cutoff)
         if (ncomp < b%sh_np(mu)) plan = PLAN_SEGMENT
         return
      end if

      price = opts%shell_cost
      if (opts%md_coulomb) price = 0.0_dp
      ! A column whose primitives the group weights very unevenly is split
      ! whatever the price: kept whole it would need the kernels' primitive
      ! cut tightened by the fourth power of that unevenness for the whole
      ! group's screen to stay inside the split's (`decontract_whole_screen`).
      if (group_unevenness(b, mem, nmem) > DECON_WHOLE_UNEVEN_MAX) price = 0.0_dp
      nmulti = 0
      i0 = 0
      cost_whole = 0.0_dp
      do i = 1, nmem
         mu = mem(i)
         cost_whole = cost_whole + real(b%sh_np(mu), dp) + price
         if (count(b%sh_e(1:b%sh_np(mu), mu) > opts%cutoff) > 1) then
            nmulti = nmulti + 1
            i0 = mu
         end if
      end do
      nunique = group_unique_exponents(b, mem, nmem)

      nseg = 0
      if (nmulti == 1) then
         call segment_mask(b, mem, nmem, i0, opts, in_segment)
         nseg = count(in_segment(1:b%sh_np(i0)))
      end if
      if (nseg >= 2) then
         ! One segment of nseg primitives, every other primitive on its own.
         cost_split = real(nseg, dp) + price + real(nunique - nseg, dp)*(1.0_dp + price)
         if (cost_split < cost_whole) plan = PLAN_SEGMENT
      else
         cost_split = real(nunique, dp)*(1.0_dp + price)
         if (cost_split < cost_whole) plan = PLAN_PRIMITIVES
      end if
   end function group_plan

   real(dp) function group_unevenness(b, mem, nmem) result(f)
      !! The largest, over the members and their primitives, of
      !! max_mu |c_k,mu| / |c_k,nu|: how much smaller a column's coefficient
      !! is than the one the split would give that primitive. 1 for a group
      !! of one.
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem
      integer, intent(in) :: mem(:)
      integer :: i, k
      f = 1.0_dp
      do i = 1, nmem
         do k = 1, b%sh_np(mem(i))
            if (b%sh_c(k, mem(i)) == 0.0_dp) cycle
            f = max(f, primitive_cbig(b, mem, nmem, b%sh_e(k, mem(i)))/abs(b%sh_c(k, mem(i))))
         end do
      end do
   end function group_unevenness

   real(dp) function primitive_cbig(b, mem, nmem, e) result(cbig)
      !! The largest |coefficient| any member gives the exponent `e`: what the
      !! split scales that primitive's shell by.
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem
      integer, intent(in) :: mem(:)
      real(dp), intent(in) :: e
      integer :: i, k
      cbig = 0.0_dp
      do i = 1, nmem
         do k = 1, b%sh_np(mem(i))
            if (same_exp(b%sh_e(k, mem(i)), e)) cbig = max(cbig, abs(b%sh_c(k, mem(i))))
         end do
      end do
   end function primitive_cbig

   !
   ! THE PRIMITIVE SCREEN OF A GENERAL CONTRACTION KEPT WHOLE.
   !
   ! Two screens act on primitives, not shells, and neither reads `amp`:
   ! the pair prune in hgp_build_pairs, which judges a primitive pair by
   ! `camp`, and the kernels' primitive-quartet cut, which judges it by the
   ! pair's own coefficient product against `pcut`. Split, a primitive's
   ! shell carries cbig, the largest coefficient any column gives it. Kept
   ! whole, it carries its column's own coefficient, which can be far
   ! smaller, so both screens drop more.
   !
   ! Over the working basis `wb` (what the integrals are evaluated in), for
   ! every group of more than one shell -- which after decontract_basis is
   ! exactly a general contraction kept whole, the split ones sharing no
   ! exponent any more:
   !
   ! - `camp(k, mu)` is cbig_k, so the prune keeps exactly the primitive
   !   pairs the split would;
   ! - `uneven` is the largest cbig_k/|c_k,mu| (group_unevenness). The
   !   kernels' cut tests a product of four coefficients, so the caller
   !   multiplies `pcut` by uneven**(-4): every primitive quartet the whole
   !   shells drop then has a cbig product at or below the split's own cut,
   !   which the split would have dropped too.
   !
   ! Everything else gets camp = sh_c and uneven = 1, the old screen.
   !
   subroutine decontract_whole_screen(wb, camp, uneven)
      type(trc_basis_t), intent(in) :: wb
      real(dp), intent(out) :: camp(:, :)
      !! (wb%maxnp, wb%nshell)
      real(dp), intent(out) :: uneven
      integer, allocatable :: grp(:)
      integer :: mem(wb%nshell)
      integer :: ng, g, nmem, i, k
      camp = wb%sh_c
      uneven = 1.0_dp
      call group_shells(wb, grp, ng)
      do g = 1, ng
         call group_members(wb, grp, g, mem, nmem)
         if (nmem < 2) cycle
         uneven = max(uneven, group_unevenness(wb, mem, nmem))
         do i = 1, nmem
            do k = 1, wb%sh_np(mem(i))
               camp(k, mem(i)) = primitive_cbig(wb, mem, nmem, wb%sh_e(k, mem(i)))
            end do
         end do
      end do
      deallocate (grp)
   end subroutine decontract_whole_screen

   subroutine group_members(b, grp, g, mem, nmem)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: g
      integer, intent(in) :: grp(:)
      integer, intent(out) :: mem(:)
      integer, intent(out) :: nmem
      integer :: i
      nmem = 0
      do i = 1, b%nshell
         if (grp(i) == g) then
            nmem = nmem + 1
            mem(nmem) = i
         end if
      end do
   end subroutine group_members

   integer function group_unique_exponents(b, mem, nmem) result(n)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem
      integer, intent(in) :: mem(:)
      integer :: i, k, ii, kk
      logical :: seen
      n = 0
      do i = 1, nmem
         do k = 1, b%sh_np(mem(i))
            seen = .false.
            do ii = 1, i
               do kk = 1, merge(k - 1, b%sh_np(mem(ii)), ii == i)
                  if (same_exp(b%sh_e(kk, mem(ii)), b%sh_e(k, mem(i)))) seen = .true.
               end do
            end do
            if (.not. seen) n = n + 1
         end do
      end do
   end function group_unique_exponents

   !
   ! Which primitives of shell `mu` stay in its segment: the ones no other
   ! member of the group uses, and, when `opts%diffuse` is on, only the
   ! compact ones among those. Everything else becomes a shared shell of one.
   !
   subroutine segment_mask(b, mem, nmem, mu, opts, in_segment)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem, mu
      integer, intent(in) :: mem(:)
      type(trc_decontract_options_t), intent(in) :: opts
      logical, intent(out) :: in_segment(:)
      integer :: k, i, kk
      in_segment = .false.
      do k = 1, b%sh_np(mu)
         in_segment(k) = .true.
         if (opts%diffuse .and. b%sh_e(k, mu) <= opts%cutoff) in_segment(k) = .false.
         do i = 1, nmem
            if (mem(i) == mu) cycle
            do kk = 1, b%sh_np(mem(i))
               if (same_exp(b%sh_e(kk, mem(i)), b%sh_e(k, mu))) in_segment(k) = .false.
            end do
         end do
      end do
   end subroutine segment_mask

   !
   ! THE SCREEN OF A GENERAL CONTRACTION KEPT WHOLE.
   !
   ! Split, every primitive p of the group becomes a shell of one scaled by
   ! the largest coefficient any column gives it, and its Schwarz bound is
   ! multiplied by amp_p, the sum over the columns nu it feeds of
   ! a_nu = sum_q |c_q,nu| / max_mu |c_q,mu| (`build_maps`). That is what
   ! makes `thresh` bound the error of the contracted matrix there, and it
   ! screens far tighter than the contracted shells do by their own bound:
   ! amp is 7 to 12 on cc-pVDZ oxygen's s, and water10/cc-pVDZ at 1e-11
   ! lands 4.1e-10 Eh from the 1e-14 energy split against 1.0e-09 whole.
   !
   ! So that the split policy changes the cost and not the accuracy, a
   ! column kept whole carries the largest amplitude any of its primitives
   ! would have had split, A_mu = max over p in mu of amp_p, and its pairs
   ! are judged by Q_mu Q_nu A_mu A_nu. That is not a proof that it drops
   ! no more than the split: the contracted bound and the primitive ones
   ! bound different sums. It matches the two screens' strictness, and it
   ! only ever tightens the screen. The two primitive-level screens are
   ! matched by construction instead (`decontract_whole_screen`). Measured
   ! (V100, RHF, 1e-11, converged) after hgp_build_pairs' sort was fixed:
   ! the silica slice with Si's p kept whole is 1.2e-08 Eh from the 1e-14
   ! energy with none of the three, 3.9e-09 with them, 4.4e-09 all split.
   !
   real(dp) function whole_amplitude(b, mem, nmem, mu) result(a)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem, mu
      integer, intent(in) :: mem(:)
      integer :: k, i, kk
      real(dp) :: amp_p
      a = 1.0_dp
      do k = 1, b%sh_np(mu)
         amp_p = 0.0_dp
         do i = 1, nmem
            do kk = 1, b%sh_np(mem(i))
               if (same_exp(b%sh_e(kk, mem(i)), b%sh_e(k, mu))) then
                  amp_p = amp_p + column_sum(b, mem, nmem, mem(i))
               end if
            end do
         end do
         a = max(a, amp_p)
      end do
   end function whole_amplitude

   real(dp) function column_sum(b, mem, nmem, nu) result(a)
      !! a_nu: the column sum of |C| contracted shell `nu` would have split,
      !! each primitive scaled by the largest coefficient the group gives it.
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem, nu
      integer, intent(in) :: mem(:)
      integer :: k, i, kk
      real(dp) :: cbig
      a = 0.0_dp
      do k = 1, b%sh_np(nu)
         cbig = 0.0_dp
         do i = 1, nmem
            do kk = 1, b%sh_np(mem(i))
               if (same_exp(b%sh_e(kk, mem(i)), b%sh_e(k, nu))) cbig = max(cbig, abs(b%sh_c(kk, mem(i))))
            end do
         end do
         if (cbig > 0.0_dp) a = a + abs(b%sh_c(k, nu))/cbig
      end do
   end function column_sum

   subroutine group_whole_amplitudes(b, grp, g, amp)
      !! `whole_amplitude` for every member of group `g` when it has more than
      !! one, into `amp` over the caller's shells.
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: g
      integer, intent(in) :: grp(:)
      real(dp), intent(inout) :: amp(:)
      integer :: mem(b%nshell)
      integer :: nmem, i
      call group_members(b, grp, g, mem, nmem)
      if (nmem < 2) return
      do i = 1, nmem
         amp(mem(i)) = whole_amplitude(b, mem, nmem, mem(i))
      end do
   end subroutine group_whole_amplitudes

   !
   ! Emit the shells one group becomes, and the triples that map the group's
   ! contracted shells onto them.
   !
   subroutine split_group(b, grp, g, opts, bld)
      type(trc_basis_t), intent(in) :: b
      integer,  intent(in) :: g
      integer,  intent(in) :: grp(:)
      type(trc_decontract_options_t), intent(in) :: opts
      type(decon_build_t), intent(inout) :: bld

      integer :: nmem, i, mu, k, plan, i0, first
      integer :: mem(b%nshell)
      logical :: in_segment(b%maxnp)

      call group_members(b, grp, g, mem, nmem)
      plan = group_plan(b, grp, g, opts)

      ! A segmented shell, or a general contraction whose split costs more
      ! than it saves: keep every member whole, mapped one to one, exactly as
      ! the unsplit path would have them.
      if (plan == PLAN_WHOLE) then
         in_segment = .true.
         do i = 1, nmem
            call emit_segment(b, mem(i), in_segment(1:b%sh_np(mem(i))), bld)
            if (nmem > 1 .and. opts%screen_whole) then
               bld%p_screen(bld%nsh_p) = whole_amplitude(b, mem, nmem, mem(i))
            end if
         end do
         return
      end if

      ! The column that keeps a segment: the only one with more than one
      ! compact primitive, or the group's one shell.
      i0 = 0
      if (plan == PLAN_SEGMENT) then
         do i = 1, nmem
            if (nmem == 1 .or. count(b%sh_e(1:b%sh_np(mem(i)), mem(i)) > opts%cutoff) > 1) i0 = mem(i)
         end do
         if (nmem == 1) then
            ! A segmented shell: the compact primitives stay together.
            in_segment = .false.
            in_segment(1:b%sh_np(i0)) = b%sh_e(1:b%sh_np(i0), i0) > opts%cutoff
         else
            call segment_mask(b, mem, nmem, i0, opts, in_segment)
         end if
         if (count(in_segment(1:b%sh_np(i0))) > 0) then
            call emit_segment(b, i0, in_segment(1:b%sh_np(i0)), bld)
         end if
      end if

      ! Everything else becomes one shell per primitive. An exponent already
      ! emitted for this group is reused, so the columns of a general
      ! contraction share their primitive shells -- which is the entire point.
      first = bld%nsh_p + 1
      do i = 1, nmem
         mu = mem(i)
         do k = 1, b%sh_np(mu)
            if (mu == i0) then
               if (in_segment(k)) cycle
            end if
            call emit_primitive(b, mem, nmem, mu, k, first, bld)
         end do
      end do
   end subroutine split_group

   !
   ! Primitive `k` of contracted shell `mu`, as a shared shell of one: found
   ! among the ones this group has emitted since `first`, or emitted now.
   !
   subroutine emit_primitive(b, mem, nmem, mu, k, first, bld)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nmem, mu, k, first
      integer, intent(in) :: mem(:)
      type(decon_build_t), intent(inout) :: bld
      integer :: s, kk, kc
      real(dp) :: e, cbig

      e = b%sh_e(k, mu)
      s = 0
      do kk = first, bld%nsh_p
         if (.not. bld%p_single(kk)) cycle
         if (same_exp(bld%p_e(1, kk), e)) then
            s = kk
            exit
         end if
      end do
      if (s == 0) then
         ! Scale the primitive shell by the largest coefficient any column
         ! gives it, so every alpha lands in [-1, 1] and the transform stays
         ! well conditioned whatever the basis does with tight exponents.
         cbig = 0.0_dp
         do kk = 1, nmem
            do kc = 1, b%sh_np(mem(kk))
               if (same_exp(b%sh_e(kc, mem(kk)), e)) then
                  if (abs(b%sh_c(kc, mem(kk))) > abs(cbig)) cbig = b%sh_c(kc, mem(kk))
               end if
            end do
         end do
         bld%nsh_p = bld%nsh_p + 1
         bld%p_l(bld%nsh_p) = b%sh_l(mu)
         bld%p_np(bld%nsh_p) = 1
         bld%p_e(1, bld%nsh_p) = e
         bld%p_c(1, bld%nsh_p) = cbig
         bld%p_r(:, bld%nsh_p) = b%sh_r(:, mu)
         bld%p_single(bld%nsh_p) = .true.
         s = bld%nsh_p
      end if
      ! Numerator only; the denominator is the coefficient the emitted shell
      ! actually ends up with, filled in by the caller once the basis is
      ! built.
      bld%ntr = bld%ntr + 1
      bld%tr_mu(bld%ntr) = mu
      bld%tr_s(bld%ntr) = s
      bld%tr_a(bld%ntr) = b%sh_c(k, mu)
   end subroutine emit_primitive

   !
   ! Contracted shell `mu` restricted to the primitives `keep` marks, with
   ! their own coefficients: the whole shell when every one is kept.
   !
   subroutine emit_segment(b, mu, keep, bld)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: mu
      logical, intent(in) :: keep(:)
      type(decon_build_t), intent(inout) :: bld
      integer :: k, n, kfirst
      bld%nsh_p = bld%nsh_p + 1
      bld%p_l(bld%nsh_p) = b%sh_l(mu)
      bld%p_r(:, bld%nsh_p) = b%sh_r(:, mu)
      bld%p_single(bld%nsh_p) = .false.
      n = 0
      kfirst = 0
      do k = 1, b%sh_np(mu)
         if (.not. keep(k)) cycle
         if (kfirst == 0) kfirst = k
         n = n + 1
         bld%p_e(n, bld%nsh_p) = b%sh_e(k, mu)
         bld%p_c(n, bld%nsh_p) = b%sh_c(k, mu)
      end do
      bld%p_np(bld%nsh_p) = n
      ! NOT alpha = 1. The constructor scales every coefficient it is given
      ! by common_fac_sp, so a shell handed back its own coefficients comes
      ! out that factor too large -- uniformly, across all its primitives.
      ! Recording the first kept coefficient and dividing by the one the
      ! built shell ends up with recovers the ratio without this module
      ! having to know the constant. Assuming alpha = 1 here put every s and
      ! p shell out by 0.28 and 0.49 and was invisible to check_gc, which
      ! compares terco's two paths to each other and not to anything
      ! external.
      bld%ntr = bld%ntr + 1
      bld%tr_mu(bld%ntr) = mu
      bld%tr_s(bld%ntr) = bld%nsh_p
      bld%tr_a(bld%ntr) = b%sh_c(kfirst, mu)
   end subroutine emit_segment

   !
   ! Expand the per-shell triples into per-AO CSR, both directions.
   !
   subroutine build_maps(b, pb, tr_mu, tr_s, tr_a, ntr, m)
      type(trc_basis_t), intent(in) :: b, pb
      integer,  intent(in) :: ntr
      integer,  intent(in) :: tr_mu(:), tr_s(:)
      real(dp), intent(in) :: tr_a(:)
      type(decon_t), intent(inout) :: m

      integer :: t, comp, nc, nnz, p, c, i
      integer, allocatable :: cnt(:), pos(:)

      m%nao_c = b%nao
      m%nao_p = pb%nao
      nnz = 0
      do t = 1, ntr
         nnz = nnz + ncart(b%sh_l(tr_mu(t)))
      end do
      m%nnz = nnz
      allocate (m%f_off(m%nao_p + 1), m%f_c(nnz), m%f_a(nnz))
      allocate (m%t_off(m%nao_c + 1), m%t_p(nnz), m%t_a(nnz))

      ! forward: rows are primitive AOs
      allocate (cnt(max(m%nao_p, m%nao_c) + 1), pos(max(m%nao_p, m%nao_c) + 1))
      cnt = 0
      do t = 1, ntr
         nc = ncart(b%sh_l(tr_mu(t)))
         do comp = 0, nc - 1
            p = pb%sh_ao(tr_s(t)) + comp
            cnt(p) = cnt(p) + 1
         end do
      end do
      m%f_off(1) = 1
      do i = 1, m%nao_p
         m%f_off(i + 1) = m%f_off(i) + cnt(i)
      end do
      pos(1:m%nao_p) = m%f_off(1:m%nao_p)
      do t = 1, ntr
         nc = ncart(b%sh_l(tr_mu(t)))
         do comp = 0, nc - 1
            p = pb%sh_ao(tr_s(t)) + comp
            c = b%sh_ao(tr_mu(t)) + comp
            m%f_c(pos(p)) = c
            m%f_a(pos(p)) = tr_a(t)
            pos(p) = pos(p) + 1
         end do
      end do

      ! transpose: rows are contracted AOs
      cnt = 0
      do t = 1, ntr
         nc = ncart(b%sh_l(tr_mu(t)))
         do comp = 0, nc - 1
            c = b%sh_ao(tr_mu(t)) + comp
            cnt(c) = cnt(c) + 1
         end do
      end do
      m%t_off(1) = 1
      do i = 1, m%nao_c
         m%t_off(i + 1) = m%t_off(i) + cnt(i)
      end do
      pos(1:m%nao_c) = m%t_off(1:m%nao_c)
      do t = 1, ntr
         nc = ncart(b%sh_l(tr_mu(t)))
         do comp = 0, nc - 1
            p = pb%sh_ao(tr_s(t)) + comp
            c = b%sh_ao(tr_mu(t)) + comp
            m%t_p(pos(c)) = p
            m%t_a(pos(c)) = tr_a(t)
            pos(c) = pos(c) + 1
         end do
      end do
      !
      ! THE SCREENING AMPLITUDE.
      !
      ! A primitive quartet's contribution to the matrix the CALLER sees is
      ! not the primitive integral, it is that integral folded back:
      !
      !   dG_uv = c_pu c_qv (pq|rs) Dp_rs
      !
      ! and one contracted quartet (uv|ls) is the sum of every primitive
      ! quartet its four shells can form. If each of those is dropped at
      ! `thresh`, the contracted quartet is wrong by
      !
      !   a_u a_v a_l a_s thresh,   a_u = sum_p |c_pu|
      !
      ! -- the COLUMN SUMS of |C|, not its largest entry. max|C| is 1 here,
      ! which is why the per-quartet bound looks safe and is not: on
      ! cc-pVDZ water max a_u is 7.03, so 7.03^4 = 2400 dropped primitive
      ! quartets' worth of error lands in one contracted quartet. Measured:
      ! E_RHF 1.04e-05 off the pyscf reference, against 4.2e-10 for the
      ! unsplit basis, and the whole of it in the kernel's exact Schwarz
      ! test. Tightening that test by 1e-4 by hand recovered 7e-13.
      !
      ! So the bound a pair is judged by carries the amplitude of both its
      ! shells, and `thresh` means in this basis what it means in the
      ! caller's. Per shell rather than one global factor, because the
      ! shells that need it are the general-contraction s and p -- every d
      ! in cc-pVDZ is a single primitive with a_u = 1 and must not be
      ! screened any harder than before.
      !
      ! SUMMED over the contracted shells a primitive feeds, not maxed.
      ! The max bounds the error in any ONE element of G; a general
      ! contraction points the same primitive at two or three contracted
      ! functions, one dropped pair damages all of them, and an energy
      ! adds those errors up. The max left water at 3.1e-09 -- the right
      ! shape, an order short; the sum lands where the unsplit basis does.
      !
      allocate (m%amp(pb%nshell))
      m%amp = 0.0_dp
      block
         real(dp), allocatable :: acol(:)
         real(dp) :: amu
         integer :: k, mu, sp, na
         allocate (acol(m%nao_c))
         do c = 1, m%nao_c
            acol(c) = 0.0_dp
            do k = m%t_off(c), m%t_off(c + 1) - 1
               acol(c) = acol(c) + abs(m%t_a(k))
            end do
         end do
         do t = 1, ntr
            mu = tr_mu(t)
            sp = tr_s(t)
            na = ncart(b%sh_l(mu))
            amu = maxval(acol(b%sh_ao(mu):b%sh_ao(mu) + na - 1))
            m%amp(sp) = m%amp(sp) + amu
         end do
         deallocate (acol)
         ! A shell no transform entry names cannot happen -- every emitted
         ! shell has one -- but a zero here would silently screen it away
         ! entirely, so it does not get the chance.
         do sp = 1, pb%nshell
            if (m%amp(sp) <= 0.0_dp) m%amp(sp) = 1.0_dp
         end do
      end block

      deallocate (cnt, pos)
   end subroutine build_maps

   !
   ! D_p = C D C^T, on the device, as a gather over the sparse rows of C.
   !
   ! Both matrices must already be present. `nao_p^2 x (nnz/nao_p)^2` work,
   ! about twenty million operations on the silica slice, against a Fock build
   ! of teraflops.
   !
   subroutine decontract_expand(m, d, dp_out)
      type(decon_t), intent(in) :: m
      real(dp), intent(in)  :: d(m%nao_c, m%nao_c)
      real(dp), intent(out) :: dp_out(m%nao_p, m%nao_p)
      ! Component arrays, not `m%...`, inside the loop: naming a derived type
      ! there makes nvfortran map the whole thing and fail at "partially
      ! present". Same reason as fold_dsh_arrays in trc_bins.
      call expand_kernel(m%nao_c, m%nao_p, m%nnz, m%f_off, m%f_c, m%f_a, d, dp_out)
   end subroutine decontract_expand

   subroutine expand_kernel(nao_c, nao_p, nnz, f_off, f_c, f_a, d, dp_out)
      integer,  intent(in) :: nao_c, nao_p, nnz
      integer,  intent(in) :: f_off(nao_p + 1), f_c(nnz)
      real(dp), intent(in) :: f_a(nnz), d(nao_c, nao_c)
      real(dp), intent(out) :: dp_out(nao_p, nao_p)
      integer :: ij, i, j, ii, jj
      real(dp) :: acc
      ! One flat index, not `do concurrent(j=..., i=...)`. With the
      ! two-dimensional header nvfortran did not give every iteration its own
      ! accumulator despite the `local` clause, and the transform came back
      ! wrong by order one -- while the identical host loop was exact, which
      ! is what made it look like a mapping problem rather than a codegen one.
      do concurrent(ij = 0:nao_p*nao_p - 1) local(acc, i, j, ii, jj)
         j = ij/nao_p + 1
         i = ij - (j - 1)*nao_p + 1
         acc = 0.0_dp
         do ii = f_off(i), f_off(i + 1) - 1
            do jj = f_off(j), f_off(j + 1) - 1
               acc = acc + f_a(ii)*f_a(jj)*d(f_c(ii), f_c(jj))
            end do
         end do
         dp_out(i, j) = acc
      end do
   end subroutine expand_kernel

   !
   ! G = C^T G_p C, the fold back. A gather again, over the transpose map, so
   ! no atomics and no ordering.
   !
   subroutine decontract_fold(m, gp, g)
      type(decon_t), intent(in) :: m
      real(dp), intent(in)  :: gp(m%nao_p, m%nao_p)
      real(dp), intent(out) :: g(m%nao_c, m%nao_c)
      call fold_kernel(m%nao_c, m%nao_p, m%nnz, m%t_off, m%t_p, m%t_a, gp, g)
   end subroutine decontract_fold

   subroutine fold_kernel(nao_c, nao_p, nnz, t_off, t_p, t_a, gp, g)
      integer,  intent(in) :: nao_c, nao_p, nnz
      integer,  intent(in) :: t_off(nao_c + 1), t_p(nnz)
      real(dp), intent(in) :: t_a(nnz), gp(nao_p, nao_p)
      real(dp), intent(out) :: g(nao_c, nao_c)
      integer :: ij, i, j, ii, jj
      real(dp) :: acc
      ! Flattened for the same reason as expand_kernel above.
      do concurrent(ij = 0:nao_c*nao_c - 1) local(acc, i, j, ii, jj)
         j = ij/nao_c + 1
         i = ij - (j - 1)*nao_c + 1
         acc = 0.0_dp
         do ii = t_off(i), t_off(i + 1) - 1
            do jj = t_off(j), t_off(j + 1) - 1
               acc = acc + t_a(ii)*t_a(jj)*gp(t_p(ii), t_p(jj))
            end do
         end do
         g(i, j) = acc
      end do
   end subroutine fold_kernel

   !
   ! The same two transforms with the loops left on the host.
   !
   ! `eri_fock`, `eri_fock_nosym` and `eri_fock_many` take a host density and
   ! return a host matrix; running `do concurrent` over arrays that were never
   ! mapped would mean an implicit copy of both in and out for twenty
   ! megaflops of work. `eri_fock_resident` uses the device pair above,
   ! because there everything already lives on the GPU and pulling it down
   ! would defeat the whole routine.
   !
   subroutine decontract_expand_host(m, d, dp_out)
      type(decon_t), intent(in) :: m
      real(dp), intent(in)  :: d(m%nao_c, m%nao_c)
      real(dp), intent(out) :: dp_out(m%nao_p, m%nao_p)
      integer :: i, j, ii, jj
      real(dp) :: acc
      do j = 1, m%nao_p
         do i = 1, m%nao_p
            acc = 0.0_dp
            do ii = m%f_off(i), m%f_off(i + 1) - 1
               do jj = m%f_off(j), m%f_off(j + 1) - 1
                  acc = acc + m%f_a(ii)*m%f_a(jj)*d(m%f_c(ii), m%f_c(jj))
               end do
            end do
            dp_out(i, j) = acc
         end do
      end do
   end subroutine decontract_expand_host

   subroutine decontract_fold_host(m, gp, g)
      type(decon_t), intent(in) :: m
      real(dp), intent(in)  :: gp(m%nao_p, m%nao_p)
      real(dp), intent(out) :: g(m%nao_c, m%nao_c)
      integer :: i, j, ii, jj
      real(dp) :: acc
      do j = 1, m%nao_c
         do i = 1, m%nao_c
            acc = 0.0_dp
            do ii = m%t_off(i), m%t_off(i + 1) - 1
               do jj = m%t_off(j), m%t_off(j + 1) - 1
                  acc = acc + m%t_a(ii)*m%t_a(jj)*gp(m%t_p(ii), m%t_p(jj))
               end do
            end do
            g(i, j) = acc
         end do
      end do
   end subroutine decontract_fold_host

   subroutine decontract_release(m)
      type(decon_t), intent(inout) :: m
      if (m%on_device) then
         !$acc exit data delete(m%f_off, m%f_c, m%f_a, m%t_off, m%t_p, m%t_a)
         !$acc exit data delete(m)
         m%on_device = .false.
      end if
      if (allocated(m%f_off)) deallocate (m%f_off, m%f_c, m%f_a)
      if (allocated(m%t_off)) deallocate (m%t_off, m%t_p, m%t_a)
      if (allocated(m%amp)) deallocate (m%amp)
      m%active = .false.
      m%nao_c = 0
      m%nao_p = 0
      m%nnz = 0
   end subroutine decontract_release

end module trc_decontract
