!
! Splitting a general contraction into primitive shells, and carrying the
! contraction as a matrix transform instead.
!
module trc_decontract
   use pic_types, only: dp
   use pic_logger, only: global_logger
   use trc_api, only: trc_basis_t
   implicit none
   private

   public :: decon_t, decontract_basis, decon_expand, decon_fold, decon_release
   public :: decon_expand_host, decon_fold_host
   public :: DECON_CUTOFF_DEFAULT, DECON_MINPRIM_DEFAULT

   real(dp), parameter :: DECON_CUTOFF_DEFAULT = 0.3_dp
   !! Exponents at or below this are diffuse. GPU4PySCF's `diffuse_cutoff`.

   !
   integer, parameter :: DECON_MINPRIM_DEFAULT = 0
   !! SPLIT ONLY THE GROUPS THAT PAY FOR IT -- OFF BY DEFAULT.

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
      logical :: on_device = .false.
   end type decon_t

contains

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
   subroutine decontract_basis(b, pb, m, cutoff, minprim)
      type(trc_basis_t), intent(in)    :: b
      type(trc_basis_t), intent(inout) :: pb
      type(decon_t),     intent(inout) :: m
      real(dp), intent(in), optional :: cutoff
      !! Diffuse threshold; DECON_CUTOFF_DEFAULT if absent.
      integer, intent(in), optional :: minprim
      !! Smallest compact primitive count that makes a split worth its
      !! screen; DECON_MINPRIM_DEFAULT if absent. Zero splits everything.

      integer, allocatable :: grp(:)          ! group id per shell
      integer, allocatable :: p_l(:), p_np(:)
      real(dp), allocatable :: p_e(:, :), p_c(:, :), p_r(:, :)
      ! one triple per (contracted shell, emitted shell) pair
      integer, allocatable :: tr_mu(:), tr_s(:)
      real(dp), allocatable :: tr_a(:)
      logical, allocatable :: tr_whole(:)
      integer :: t
      integer :: ntr, nsh_p, maxnp_p
      integer :: ng, g, i, n
      real(dp) :: cut
      integer, parameter :: LTOP = 8
      real(dp) :: fl(0:LTOP)
      integer :: mnp

      cut = DECON_CUTOFF_DEFAULT
      if (present(cutoff)) cut = cutoff
      mnp = DECON_MINPRIM_DEFAULT
      if (present(minprim)) mnp = minprim

      call decon_release(m)
      call group_shells(b, grp, ng)

      ! Nothing shares an exponent with anything, or nothing that does is
      ! worth the split: a segmented basis as far as this module cares. Leave.
      if (.not. any_group_splits(b, grp, ng, cut, mnp)) then
         m%active = .false.
         return
      end if

      ! Upper bounds: every primitive of every shell could become its own
      ! shell, and each contributes one triple.
      n = 0
      do i = 1, b%nshell
         n = n + b%sh_np(i)
      end do
      allocate (p_l(n + b%nshell), p_np(n + b%nshell), p_r(3, n + b%nshell))
      allocate (p_e(b%maxnp, n + b%nshell), p_c(b%maxnp, n + b%nshell))
      allocate (tr_mu(n + b%nshell), tr_s(n + b%nshell), tr_a(n + b%nshell))
      allocate (tr_whole(n + b%nshell))
      p_e = 0.0_dp
      p_c = 0.0_dp
      nsh_p = 0
      ntr = 0

      do g = 1, ng
         call split_group(b, grp, g, cut, mnp, p_l, p_np, p_e, p_c, p_r, nsh_p, &
                          tr_mu, tr_s, tr_a, tr_whole, ntr)
      end do

      maxnp_p = 1
      do i = 1, nsh_p
         maxnp_p = max(maxnp_p, p_np(i))
      end do
      call pb%build(nsh_p, p_l(1:nsh_p), p_np(1:nsh_p), p_e(1:maxnp_p, 1:nsh_p), &
                    p_c(1:maxnp_p, 1:nsh_p), p_r(:, 1:nsh_p), &
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
      do i = 1, nsh_p
         if (p_c(1, i) /= 0.0_dp) fl(p_l(i)) = pb%sh_c(1, i)/p_c(1, i)
      end do
      do i = 1, nsh_p
         if (fl(p_l(i)) /= 0.0_dp) p_c(:, i) = p_c(:, i)/fl(p_l(i))
      end do
      call pb%build(nsh_p, p_l(1:nsh_p), p_np(1:nsh_p), p_e(1:maxnp_p, 1:nsh_p), &
                    p_c(1:maxnp_p, 1:nsh_p), p_r(:, 1:nsh_p), &
                    b%natm, b%at_z, b%at_r, maxnp_p, at_elem=b%at_elem)

      ! The coefficient a split shell ENDS UP with is not the one handed to
      ! the constructor: `basis_build` folds common_fac_sp into it. Deriving
      ! alpha from the built basis rather than from the input is exact
      ! whatever the constructor does to the numbers on the way through --
      ! computing it beforehand left every s and p shell short by that
      ! factor, which is a 4e-9 error in G and a silent one.
      do t = 1, ntr
         tr_a(t) = tr_a(t)/pb%sh_c(1, tr_s(t))
      end do

      call build_maps(b, pb, tr_mu, tr_s, tr_a, ntr, m)
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
      deallocate (p_l, p_np, p_e, p_c, p_r, tr_mu, tr_s, tr_a, tr_whole, grp)
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
   logical function any_group_splits(b, grp, ng, cut, minprim) result(res)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ng, minprim
      integer, intent(in) :: grp(:)
      real(dp), intent(in) :: cut
      integer :: g
      res = .false.
      do g = 1, ng
         if (group_pays(b, grp, g, cut, minprim)) then
            res = .true.
            return
         end if
      end do
   end function any_group_splits

   !
   ! Is this group worth splitting?
   !
   ! It has to BE a general contraction -- more than one shell sharing the
   ! exponents -- and its work has to be worth the screen the split tightens.
   ! The proxy for the work is the compact primitive count of the group's
   ! biggest contracted shell, which is what the unsplit path pays twice over:
   ! once in the primitive pairs of every column, and once in `build_pairs`
   ! and `schwarz_bounds`, where msn's setup went 2.0 s -> 9.7 s with the
   ! split off. See DECON_MINPRIM_DEFAULT for the measurement.
   !
   logical function group_pays(b, grp, g, cut, minprim) result(res)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: g, minprim
      integer, intent(in) :: grp(:)
      real(dp), intent(in) :: cut
      integer :: i, k, nmem, ncomp, best
      nmem = 0
      best = 0
      do i = 1, b%nshell
         if (grp(i) /= g) cycle
         nmem = nmem + 1
         ncomp = 0
         do k = 1, b%sh_np(i)
            if (b%sh_e(k, i) > cut) ncomp = ncomp + 1
         end do
         best = max(best, ncomp)
      end do
      res = (nmem > 1) .and. (best >= minprim)
   end function group_pays

   !
   ! Emit the shells one group becomes, and the triples that map the group's
   ! contracted shells onto them.
   !
   subroutine split_group(b, grp, g, cut, minprim, p_l, p_np, p_e, p_c, p_r, nsh_p, &
                          tr_mu, tr_s, tr_a, tr_whole, ntr)
      type(trc_basis_t), intent(in) :: b
      integer,  intent(in) :: g, minprim
      integer,  intent(in) :: grp(:)
      real(dp), intent(in) :: cut
      integer,  intent(inout) :: nsh_p, ntr
      integer,  intent(inout) :: p_l(:), p_np(:), tr_mu(:), tr_s(:)
      real(dp), intent(inout) :: p_e(:, :), p_c(:, :), p_r(:, :), tr_a(:)
      logical,  intent(inout) :: tr_whole(:)

      integer :: nmem, i, mu, k, kk, kc, s, ncomp, nmulti, i0
      integer :: mem(b%nshell)
      real(dp) :: e, cbig
      logical :: seen

      nmem = 0
      do i = 1, b%nshell
         if (grp(i) == g) then
            nmem = nmem + 1
            mem(nmem) = i
         end if
      end do

      ! A segmented shell, or a general contraction whose columns are too
      ! cheap to be worth the tightened screen: keep every member whole,
      ! mapped one to one, exactly as the unsplit path would have them.
      if (.not. group_pays(b, grp, g, cut, minprim)) then
         do i = 1, nmem
            call emit_whole(b, mem(i), p_l, p_np, p_e, p_c, p_r, nsh_p, &
                            tr_mu, tr_s, tr_a, tr_whole, ntr)
         end do
         return
      end if

      ! Does exactly one column hold more than one compact primitive?
      nmulti = 0
      i0 = 0
      do i = 1, nmem
         mu = mem(i)
         ncomp = 0
         do k = 1, b%sh_np(mu)
            if (b%sh_e(k, mu) > cut) ncomp = ncomp + 1
         end do
         if (ncomp > 1) then
            nmulti = nmulti + 1
            i0 = mu
         end if
      end do

      do i = 1, nmem
         mu = mem(i)
         if (nmulti == 1 .and. mu == i0) then
            ! The one real contraction in the group survives intact.
            call emit_whole(b, mu, p_l, p_np, p_e, p_c, p_r, nsh_p, &
                            tr_mu, tr_s, tr_a, tr_whole, ntr)
            cycle
         end if
         ! Everything else becomes one shell per primitive. An exponent
         ! already emitted for this group is reused, so the columns of a
         ! general contraction share their primitive shells -- which is the
         ! entire point.
         do k = 1, b%sh_np(mu)
            e = b%sh_e(k, mu)
            s = 0
            do kk = nsh_p, 1, -1
               if (p_np(kk) /= 1) cycle
               if (p_l(kk) /= b%sh_l(mu)) cycle
               if (maxval(abs(p_r(:, kk) - b%sh_r(:, mu))) > 1.0e-10_dp) cycle
               if (same_exp(p_e(1, kk), e)) then
                  s = kk
                  exit
               end if
            end do
            if (s == 0) then
               ! Scale the primitive shell by the largest coefficient any
               ! column gives it, so every alpha lands in [-1, 1] and the
               ! transform stays well conditioned whatever the basis does
               ! with tight exponents.
               cbig = 0.0_dp
               do kk = 1, nmem
                  do kc = 1, b%sh_np(mem(kk))
                     if (same_exp(b%sh_e(kc, mem(kk)), e)) then
                        if (abs(b%sh_c(kc, mem(kk))) > abs(cbig)) cbig = b%sh_c(kc, mem(kk))
                     end if
                  end do
               end do
               nsh_p = nsh_p + 1
               p_l(nsh_p) = b%sh_l(mu)
               p_np(nsh_p) = 1
               p_e(1, nsh_p) = e
               p_c(1, nsh_p) = cbig
               p_r(:, nsh_p) = b%sh_r(:, mu)
               s = nsh_p
            end if
            ! Numerator only; the denominator is the coefficient the
            ! emitted shell actually ends up with, filled in by the caller
            ! once the basis is built.
            ntr = ntr + 1
            tr_mu(ntr) = mu
            tr_s(ntr) = s
            tr_a(ntr) = b%sh_c(k, mu)
            tr_whole(ntr) = .false.
         end do
      end do
   end subroutine split_group

   subroutine emit_whole(b, mu, p_l, p_np, p_e, p_c, p_r, nsh_p, &
                         tr_mu, tr_s, tr_a, tr_whole, ntr)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: mu
      integer,  intent(inout) :: nsh_p, ntr
      integer,  intent(inout) :: p_l(:), p_np(:), tr_mu(:), tr_s(:)
      real(dp), intent(inout) :: p_e(:, :), p_c(:, :), p_r(:, :), tr_a(:)
      logical,  intent(inout) :: tr_whole(:)
      integer :: k
      nsh_p = nsh_p + 1
      p_l(nsh_p) = b%sh_l(mu)
      p_np(nsh_p) = b%sh_np(mu)
      p_r(:, nsh_p) = b%sh_r(:, mu)
      do k = 1, b%sh_np(mu)
         p_e(k, nsh_p) = b%sh_e(k, mu)
         p_c(k, nsh_p) = b%sh_c(k, mu)
      end do
      ! NOT alpha = 1. The constructor scales every coefficient it is given
      ! by common_fac_sp, so a shell handed back its own coefficients comes
      ! out that factor too large -- uniformly, across all its primitives.
      ! Recording the first coefficient and dividing by the one the built
      ! shell ends up with recovers the ratio without this module having to
      ! know the constant. Assuming alpha = 1 here put every s and p shell
      ! out by 0.28 and 0.49 and was invisible to check_gc, which compares
      ! terco's two paths to each other and not to anything external.
      ntr = ntr + 1
      tr_mu(ntr) = mu
      tr_s(ntr) = nsh_p
      tr_a(ntr) = b%sh_c(1, mu)
      tr_whole(ntr) = .true.
   end subroutine emit_whole

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
   subroutine decon_expand(m, d, dp_out)
      type(decon_t), intent(in) :: m
      real(dp), intent(in)  :: d(m%nao_c, m%nao_c)
      real(dp), intent(out) :: dp_out(m%nao_p, m%nao_p)
      ! Component arrays, not `m%...`, inside the loop: naming a derived type
      ! there makes nvfortran map the whole thing and fail at "partially
      ! present". Same reason as fold_dsh_arrays in trc_bins.
      call expand_kernel(m%nao_c, m%nao_p, m%nnz, m%f_off, m%f_c, m%f_a, d, dp_out)
   end subroutine decon_expand

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
   subroutine decon_fold(m, gp, g)
      type(decon_t), intent(in) :: m
      real(dp), intent(in)  :: gp(m%nao_p, m%nao_p)
      real(dp), intent(out) :: g(m%nao_c, m%nao_c)
      call fold_kernel(m%nao_c, m%nao_p, m%nnz, m%t_off, m%t_p, m%t_a, gp, g)
   end subroutine decon_fold

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
   subroutine decon_expand_host(m, d, dp_out)
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
   end subroutine decon_expand_host

   subroutine decon_fold_host(m, gp, g)
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
   end subroutine decon_fold_host

   subroutine decon_release(m)
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
   end subroutine decon_release

end module trc_decontract
